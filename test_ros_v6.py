"""測試 ROS 參數 v6 - 恢復時間策略

核心洞察：
問題不在於 ROS 累積本身，而在於**恢復不足**：
- L6D6: 每天 6h UVA + 18h 恢復 → 充分恢復
- H12D3: 每天 12h UVA + 12h 恢復 → 恢復不足

策略：讓**持續照射時間**影響 ROS 清除效率
- 長時間連續照射會耗竭清除酵素
- 但這個效應在恢復期間可以恢復
"""
import numpy as np
from scipy.integrate import solve_ivp
import sys
sys.path.insert(0, '/home/kasm-user/projects/uva-simulation')

from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment
from lettuce_uva_carbon_complete_model import SunParams as BaseSunParams
from lettuce_uva_carbon_complete_model import sun_derivatives_final


class UVAParamsV6(BaseSunParams):
    def __init__(self, params):
        super().__init__()
        for key, val in params.items():
            setattr(self, key, val)


def calculate_dynamic_dw_fw_ratio(Stress, p):
    stress_effect = p.ldmc_stress_sensitivity * Stress / (p.K_ldmc + Stress + 1e-9)
    ratio = p.dw_fw_ratio_base * (1.0 + stress_effect)
    return min(ratio, p.dw_fw_ratio_max)


def uva_derivatives_v6(t, state, p, env):
    """ROS v6 - 加入清除能力耗竭狀態變量"""
    X_d, C_buf, LAI, Anth, Stress, ROS_mM, Depletion = state

    # 邊界
    X_d = max(X_d, 1e-9)
    C_buf = max(C_buf, 0.0)
    LAI = max(LAI, 0.1)
    Anth = max(Anth, 0.0)
    Stress = max(Stress, 0.0)
    ROS_mM = max(ROS_mM, 0.0)
    Depletion = max(Depletion, 0.0)  # 清除能力耗竭程度 (0=完全健康, 高=耗竭)

    # 時間
    hour = (t / 3600.0) % 24.0
    day_from_sowing = t / 86400.0
    light_on = env['light_on_hour']
    light_off = env['light_off_hour']
    is_day = light_on <= hour < light_off if light_on <= light_off else (hour >= light_on or hour < light_off)
    I_base = env['I_day'] if is_day else 0.0
    Tc = env['T_day'] if is_day else env['T_night']

    # UVA
    uva_on = env.get('uva_on', False)
    I_UVA = 0.0
    if uva_on:
        uva_start_day = env.get('uva_start_day', 29)
        uva_end_day = env.get('uva_end_day', 35)
        uva_hour_on = env.get('uva_hour_on', 10)
        uva_hour_off = env.get('uva_hour_off', 16)
        uva_intensity = env.get('uva_intensity', 22.0)
        if uva_hour_on <= uva_hour_off:
            if uva_start_day <= day_from_sowing <= uva_end_day:
                if uva_hour_on <= hour < uva_hour_off:
                    I_UVA = uva_intensity
        else:
            if hour >= uva_hour_on:
                if uva_start_day <= day_from_sowing <= uva_end_day:
                    I_UVA = uva_intensity
            elif hour < uva_hour_off:
                day_session_started = day_from_sowing - 1
                if uva_start_day <= day_session_started <= uva_end_day:
                    I_UVA = uva_intensity

    is_night_uva = (I_UVA > 0) and (not is_day)

    # 合併光照
    I_effective = I_base + I_UVA * p.par_conversion_factor
    env_mod = env.copy()
    env_mod['I_override'] = I_effective
    env_mod['T_override'] = Tc
    env_mod['is_day_override'] = is_day

    # Sun 基礎
    base_state = [X_d, C_buf, LAI]
    dXd_dt_base, dCbuf_dt, dLAI_dt_base = sun_derivatives_final(t, base_state, p, env_mod)

    # Vulnerability
    base_vuln = (p.LAI_ref_vuln / LAI) ** p.n_vuln
    vulnerability = p.cap_vuln * base_vuln / (p.cap_vuln + base_vuln)

    # Stress 非線性
    nonlinear_factor = 1.0 + p.stress_nonlinear_coeff * Stress / (p.K_nonlinear + Stress + 1e-9)

    # ROS 放大
    ros_amplify = 1.0 + p.ros_damage_gain * ROS_mM / (p.K_ros_damage + ROS_mM + 1e-9)

    # 夜間
    circadian_penalty = p.circadian_disruption_factor if is_night_uva else 1.0

    # 損傷
    damage_rate = p.stress_damage_coeff * I_UVA * vulnerability * nonlinear_factor * ros_amplify * circadian_penalty

    # 修復
    repair_capacity = p.base_repair_capacity + p.carbon_repair_bonus * C_buf / (p.K_carbon + C_buf + 1e-9)
    repair_rate = p.stress_repair_coeff * Stress * repair_capacity
    dStress_dt = damage_rate - repair_rate

    # ===== ROS v6: 清除能力耗竭 =====
    # 清除能力受耗竭程度影響
    clearance_efficiency = 1.0 / (1.0 + Depletion / p.K_depletion)

    # ROS 生成
    ros_prod = p.ros_prod_coeff * I_UVA

    # ROS 清除 (受耗竭影響)
    ROS_pos = max(ROS_mM, 0.0)
    ros_apx = p.ros_vmax_apx * ROS_pos / (p.K_ros_apx + ROS_pos + 1e-9) * clearance_efficiency
    ros_cat = p.ros_vmax_cat * ROS_pos / (p.K_ros_cat + ROS_pos + 1e-9) * clearance_efficiency
    ros_decay = p.ros_passive_decay * ROS_pos

    dROS_dt = ros_prod - (ros_apx + ros_cat + ros_decay)
    if ROS_mM <= 0 and dROS_dt < 0:
        dROS_dt = 0.0

    # 耗竭動態
    # - UVA 照射時：耗竭增加
    # - 不照射時：耗竭恢復
    if I_UVA > 0:
        depletion_rate = p.depletion_rate_uva * I_UVA / 22.0  # 正規化到標準強度
    else:
        depletion_rate = 0.0
    depletion_recovery = p.depletion_recovery_rate * Depletion
    dDepletion_dt = depletion_rate - depletion_recovery

    # Stress 抑制
    stress_inhibition = Stress / (p.K_stress + Stress + 1e-9)
    xd_reduction = p.stress_photosynthesis_inhibition * stress_inhibition
    lai_reduction = p.stress_lai_inhibition * stress_inhibition
    dXd_dt = dXd_dt_base * (1.0 - xd_reduction) if dXd_dt_base > 0 else dXd_dt_base
    dLAI_dt = dLAI_dt_base * (1.0 - lai_reduction) if dLAI_dt_base > 0 else dLAI_dt_base

    # 花青素
    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress, p)
    FW_kg_m2 = X_d / dw_fw_ratio
    base_synth = (p.base_anth_rate_light if is_day else p.base_anth_rate_dark)
    stress_induced = p.V_max_anth * Stress / (p.K_stress_anth + Stress + 1e-12)
    synthesis_rate = FW_kg_m2 * (base_synth + stress_induced)
    degradation = p.k_deg * Anth
    dAnth_dt = synthesis_rate - degradation

    # 碳消耗
    repair_carbon_consumption = repair_rate * p.repair_carbon_cost
    anth_carbon_consumption = synthesis_rate * p.anth_carbon_cost
    dCbuf_dt = dCbuf_dt - repair_carbon_consumption - anth_carbon_consumption

    return np.array([dXd_dt, dCbuf_dt, dLAI_dt, dAnth_dt, dStress_dt, dROS_dt, dDepletion_dt])


def test_params(params, label):
    p = UVAParamsV6(params)

    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    C_buf_init = Xd_init * 0.1
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6
    # 7 狀態: [X_d, C_buf, LAI, Anth, Stress, ROS, Depletion]
    initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0, 0.0]

    transplant_day = SIMULATION['transplant_offset']
    simulation_days = SIMULATION['days']
    t_start = transplant_day * 86400
    t_end = (transplant_day + simulation_days) * 86400

    print(f"\n{label}")
    print("-" * 130)
    print(f"{'Treatment':<10} {'FW_obs':>7} {'FW_sim':>7} {'FW_Err':>7}   {'Anth_obs':>8} {'Anth_sim':>8} {'Anth_Err':>8}   {'Stress':>7} {'ROS':>6} {'Depl':>6}")
    print("-" * 130)

    fw_errs = []
    anth_errs = []

    for treatment in ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']:
        env = get_env_for_treatment(treatment)
        target = TARGETS.get(treatment, {'FW': 0, 'Anth': 0})

        sol = solve_ivp(
            uva_derivatives_v6,
            (t_start, t_end),
            initial_state,
            args=(p, env),
            method='RK45',
            max_step=60,
            t_eval=np.array([t_end])
        )

        if sol.success:
            Xd_f, Cbuf_f, LAI_f, Anth_f, Stress_f, ROS_f, Depl_f = sol.y[:, -1]
            dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress_f, p)
            X_total = Xd_f / ENV_BASE['plant_density']
            FW_sim = X_total / dw_fw_ratio * 1000
            FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
            Anth_sim = (Anth_f / FW_total_kg) * 1e6

            FW_obs = target['FW']
            Anth_obs = target.get('Anth', 0)

            fw_err = (FW_sim - FW_obs) / FW_obs * 100
            anth_err = (Anth_sim - Anth_obs) / Anth_obs * 100 if Anth_obs > 0 else 0

            fw_errs.append(abs(fw_err))
            anth_errs.append(abs(anth_err))

            status_fw = "✓" if abs(fw_err) < 5 else "✗"
            status_anth = "✓" if abs(anth_err) < 10 else "✗"
            print(f"{treatment:<10} {FW_obs:>6.1f}g {FW_sim:>6.1f}g {fw_err:>+6.1f}%{status_fw}  "
                  f"{Anth_obs:>7.1f} {Anth_sim:>7.1f} {anth_err:>+7.1f}%{status_anth}  {Stress_f:>7.2f} {ROS_f:>6.2f} {Depl_f:>6.2f}")
        else:
            print(f"{treatment:<10} 失敗: {sol.message}")
            fw_errs.append(100)
            anth_errs.append(100)

    print("-" * 130)
    fw_ok = sum(1 for e in fw_errs if e < 5)
    anth_ok = sum(1 for e in anth_errs if e < 10)
    print(f"FW <5%: {fw_ok}/6, Anth <10%: {anth_ok}/6, Total: {fw_ok + anth_ok}/12")
    return fw_ok + anth_ok


# 基礎參數
BASE_PARAMS = {
    'c_alpha': 0.555,
    'par_conversion_factor': 1.0,
    'stress_damage_coeff': 0.66e-6,
    'stress_repair_coeff': 1.0e-5,
    'stress_nonlinear_coeff': 8.0,
    'K_nonlinear': 0.8,
    'LAI_ref_vuln': 6.5,
    'n_vuln': 8,
    'cap_vuln': 100.0,
    'circadian_disruption_factor': 3.8,
    'stress_photosynthesis_inhibition': 0.66,
    'stress_lai_inhibition': 0.66,
    'K_stress': 1.9,
    'base_repair_capacity': 0.5,
    'carbon_repair_bonus': 0.5,
    'K_carbon': 0.001,
    'repair_carbon_cost': 1.0e-6,
    'base_anth_rate_light': 2.0e-10,
    'base_anth_rate_dark': 1.0e-10,
    'V_max_anth': 2.35e-11,
    'K_stress_anth': 0.30,
    'k_deg': 3.02e-6,
    'anth_carbon_cost': 0.0,
    'dw_fw_ratio_base': 0.05,
    'ldmc_stress_sensitivity': 1.0,
    'K_ldmc': 50.0,
    'dw_fw_ratio_max': 0.12,
    # ROS 基礎參數
    'ros_prod_coeff': 1e-5,
    'ros_vmax_apx': 0.001,
    'ros_vmax_cat': 0.002,
    'K_ros_apx': 0.074,
    'K_ros_cat': 50.0,
    'ros_passive_decay': 1e-5,
    'ros_damage_gain': 50.0,
    'K_ros_damage': 2.0,
    # v6 新參數: 清除能力耗竭
    'depletion_rate_uva': 0.001,      # UVA 造成耗竭的速率 [1/s]
    'depletion_recovery_rate': 5e-5,   # 耗竭恢復速率 [1/s]
    'K_depletion': 1.0,               # 耗竭半飽和常數
}

print("=" * 130)
print("ROS v6 - 清除能力耗竭策略")
print("說明: 長時間 UVA 照射會耗竭清除酵素，但恢復期間可以恢復")
print("=" * 130)

# 測試耗竭速率
for rate in [0.0005, 0.001, 0.002]:
    params = BASE_PARAMS.copy()
    params['depletion_rate_uva'] = rate
    test_params(params, f"耗竭速率 = {rate}")

# 調整恢復速率
print("\n" + "=" * 130)
print("調整恢復速率")
print("=" * 130)

for recovery in [1e-5, 3e-5, 5e-5]:
    params = BASE_PARAMS.copy()
    params['depletion_rate_uva'] = 0.001
    params['depletion_recovery_rate'] = recovery
    test_params(params, f"耗竭=0.001, 恢復={recovery}")
