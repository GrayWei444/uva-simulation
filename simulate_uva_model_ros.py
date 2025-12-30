"""
萵苣生長與UVA效應整合模型 - ROS 版本 v1.0
================================================
版本: v1.0 (已校準)
日期: 2025-12-29

目的:
  - 結合 ROS 動態與日內能量耗竭機制
  - 新增 ROS 狀態變量 (mM) + APX/CAT Michaelis-Menten 飽和動力學
  - 保留 v6.9 的日內能量耗竭機制 (depletion_factor)
  - 達成 12/12 精度目標

狀態變量 (6 個):
  [X_d, C_buf, LAI, Anth, Stress, ROS_mM]

文獻支持:
  - APX Km = 0.074 mM (Corpas et al. 2020, Ishikawa et al. 1998)
  - CAT Km = 50 mM (Corpas et al. 2020, Mhamdi et al. 2010)

達成精度:
  - FW: 6/6 <5%
  - Anth: 6/6 <10%
  - 總達標: 12/12 (100%)
"""

import numpy as np
from scipy.integrate import solve_ivp

from lettuce_uva_carbon_complete_model import SunParams as BaseSunParams
from lettuce_uva_carbon_complete_model import sun_derivatives_final
from model_config import ENV_BASE, TREATMENT_CONFIGS, TARGETS, ODE_SETTINGS, SIMULATION, get_env_for_treatment


def softplus(x, k=1.0):
    """Smooth max(0, x) 用於避免硬閾值。"""
    return np.log(1 + np.exp(np.clip(k * x, -500, 500))) / k


# ==============================================================================
# 參數定義 (內嵌)
# ==============================================================================
ALL_PARAMS = {
    # 基礎光合
    'c_alpha': 0.555,

    # UVA-PAR
    'par_conversion_factor': 1.0,

    # Stress 動態
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

    # 碳修復
    'base_repair_capacity': 0.5,
    'carbon_repair_bonus': 0.5,
    'K_carbon': 0.001,
    'repair_carbon_cost': 1.0e-6,

    # 花青素
    'base_anth_rate_light': 2.0e-10,
    'base_anth_rate_dark': 1.0e-10,
    'V_max_anth': 2.35e-11,
    'K_stress_anth': 0.30,
    'k_deg': 3.02e-6,
    'anth_carbon_cost': 0.0,

    # LDMC
    'dw_fw_ratio_base': 0.05,
    'ldmc_stress_sensitivity': 1.0,
    'K_ldmc': 50.0,
    'dw_fw_ratio_max': 0.12,

    # ROS 生成/清除
    # 參考: Corpas et al. 2020 (Frontiers in Plant Science)
    #       Ishikawa et al. 1998 (Plant Cell Physiol)
    #       Mhamdi et al. 2010 (J Exp Bot)
    'ros_prod_coeff': 1.0e-5,    # ROS 生成係數 [mM/(W·s)]
    'ros_passive_decay': 1.0e-4, # 被動衰減 [1/s]
    'ros_vmax_apx': 0.02,        # APX Vmax [mM/s]
    'ros_vmax_cat': 0.05,        # CAT Vmax [mM/s]
    'K_ros_apx': 0.074,          # APX Km [mM] (文獻值)
    'K_ros_cat': 50.0,           # CAT Km [mM] (文獻值)
    'ros_damage_gain': 1.0,      # ROS 對損傷的放大係數 (溫和)
    'K_ros_damage': 10.0,        # ROS 放大半飽和 [mM]

    # 日內能量耗竭 (來自 v6.9)
    'E_ref': 475.2,              # 參考能量 (kJ/m²) ≈ 6h @ 22 W/m²
    'E_scale': 237.6,            # 能量尺度 (kJ/m²)
    'k_depletion': 54.0,         # 耗竭放大係數
    'depletion_power': 2.0,      # 非線性指數
    'softplus_sharpness': 3.0,   # Softplus 銳度
}


class UVAParams(BaseSunParams):
    """UVA 參數類別 (ROS 版本 v1.0)"""

    def __init__(self, params=None):
        super().__init__()
        if params is None:
            params = ALL_PARAMS

        # 從 params 字典設置所有參數
        for key, val in params.items():
            setattr(self, key, val)


def calculate_dynamic_dw_fw_ratio(Stress, p):
    """根據 Stress 計算動態 DW:FW 比例（LDMC 效應）"""
    stress_effect = p.ldmc_stress_sensitivity * Stress / (p.K_ldmc + Stress + 1e-9)
    ratio = p.dw_fw_ratio_base * (1.0 + stress_effect)
    return min(ratio, p.dw_fw_ratio_max)


# ==============================================================================
# 微分方程
# ==============================================================================
def uva_sun_derivatives_ros(t, state, p, env):
    """
    UVA 整合微分方程（ROS 版本 v1.0）

    結合:
    1. ROS 動態 (Michaelis-Menten 清除)
    2. 日內能量耗竭 (depletion_factor)

    狀態: [X_d, C_buf, LAI, Anth, Stress, ROS_mM]
    """
    X_d, C_buf, LAI, Anth, Stress, ROS_mM = state

    # 邊界保護
    X_d = max(X_d, 1e-9)
    C_buf = max(C_buf, 0.0)
    LAI = max(LAI, 0.1)
    Anth = max(Anth, 0.0)
    Stress = max(Stress, 0.0)
    ROS_mM = max(ROS_mM, 0.0)

    # 時間與日夜
    hour = (t / 3600.0) % 24.0
    day_from_sowing = t / 86400.0
    light_on = env['light_on_hour']
    light_off = env['light_off_hour']
    is_day = light_on <= hour < light_off if light_on <= light_off else (hour >= light_on or hour < light_off)
    I_base = env['I_day'] if is_day else 0.0
    Tc = env['T_day'] if is_day else env['T_night']

    # UVA 排程
    uva_on = env.get('uva_on', False)
    uva_hour_on = env.get('uva_hour_on', 10)
    uva_hour_off = env.get('uva_hour_off', 16)
    uva_intensity = env.get('uva_intensity', 22.0)
    I_UVA = 0.0
    in_uva_window = False

    if uva_on:
        uva_start_day = env.get('uva_start_day', 29)
        uva_end_day = env.get('uva_end_day', 35)
        if uva_hour_on <= uva_hour_off:
            if uva_start_day <= day_from_sowing <= uva_end_day:
                if uva_hour_on <= hour < uva_hour_off:
                    I_UVA = uva_intensity
                    in_uva_window = True
        else:
            if hour >= uva_hour_on:
                if uva_start_day <= day_from_sowing <= uva_end_day:
                    I_UVA = uva_intensity
                    in_uva_window = True
            elif hour < uva_hour_off:
                day_session_started = day_from_sowing - 1
                if uva_start_day <= day_session_started <= uva_end_day:
                    I_UVA = uva_intensity
                    in_uva_window = True

    is_night_uva = (I_UVA > 0) and (not is_day)

    # ===== 日內累積能量 (kJ/m²) =====
    if in_uva_window:
        if uva_hour_on <= uva_hour_off:
            hours_exposed_today = hour - uva_hour_on
        else:
            if hour >= uva_hour_on:
                hours_exposed_today = hour - uva_hour_on
            else:
                hours_exposed_today = (24 - uva_hour_on) + hour
    else:
        hours_exposed_today = 0.0

    E_daily = hours_exposed_today * 3600 * uva_intensity / 1000  # kJ/m²

    # ===== 日內能量耗竭因子 (來自 v6.9) =====
    depletion_raw = softplus((E_daily - p.E_ref) / p.E_scale, p.softplus_sharpness)
    depletion_factor = 1.0 + p.k_depletion * (depletion_raw ** p.depletion_power)

    # UVA-PAR 合併到基礎光照
    I_effective = I_base + I_UVA * p.par_conversion_factor
    env_mod = env.copy()
    env_mod['I_override'] = I_effective
    env_mod['T_override'] = Tc
    env_mod['is_day_override'] = is_day

    # Sun 基礎模型
    base_state = [X_d, C_buf, LAI]
    dXd_dt_base, dCbuf_dt, dLAI_dt_base = sun_derivatives_final(t, base_state, p, env_mod)

    # Vulnerability
    base_vuln = (p.LAI_ref_vuln / LAI) ** p.n_vuln
    vulnerability = p.cap_vuln * base_vuln / (p.cap_vuln + base_vuln)

    # Stress 非線性
    nonlinear_factor = 1.0 + p.stress_nonlinear_coeff * Stress / (p.K_nonlinear + Stress + 1e-9)

    # ROS 放大 (溫和)
    ros_amplify = 1.0 + p.ros_damage_gain * ROS_mM / (p.K_ros_damage + ROS_mM + 1e-9)

    # 夜間加成
    circadian_penalty = p.circadian_disruption_factor if is_night_uva else 1.0

    # 損傷率 (結合 ROS 放大和日內耗竭)
    damage_rate = (p.stress_damage_coeff * I_UVA * vulnerability *
                   nonlinear_factor * ros_amplify * circadian_penalty * depletion_factor)

    # 碳依賴修復
    repair_capacity = p.base_repair_capacity + p.carbon_repair_bonus * C_buf / (p.K_carbon + C_buf + 1e-9)
    repair_rate = p.stress_repair_coeff * Stress * repair_capacity
    dStress_dt = damage_rate - repair_rate

    # ===== ROS 動態 =====
    ros_prod = p.ros_prod_coeff * I_UVA
    ROS_pos = max(ROS_mM, 0.0)
    ros_apx = p.ros_vmax_apx * ROS_pos / (p.K_ros_apx + ROS_pos + 1e-9)
    ros_cat = p.ros_vmax_cat * ROS_pos / (p.K_ros_cat + ROS_pos + 1e-9)
    ros_decay = p.ros_passive_decay * ROS_pos
    dROS_dt = ros_prod - (ros_apx + ros_cat + ros_decay)

    # 防止 ROS 變負
    if ROS_mM <= 0 and dROS_dt < 0:
        dROS_dt = 0.0

    # Stress 對生長抑制
    stress_inhibition = Stress / (p.K_stress + Stress + 1e-9)
    xd_reduction = p.stress_photosynthesis_inhibition * stress_inhibition
    lai_reduction = p.stress_lai_inhibition * stress_inhibition
    dXd_dt = dXd_dt_base * (1.0 - xd_reduction) if dXd_dt_base > 0 else dXd_dt_base
    dLAI_dt = dLAI_dt_base * (1.0 - lai_reduction) if dLAI_dt_base > 0 else dLAI_dt_base

    # 花青素 (FW-based)
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

    return np.array([dXd_dt, dCbuf_dt, dLAI_dt, dAnth_dt, dStress_dt, dROS_dt])


# ==============================================================================
# 驗證函數
# ==============================================================================
def validate_all_treatments():
    """驗證所有處理組的結果"""
    p = UVAParams()

    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    C_buf_init = Xd_init * 0.1
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6
    initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]

    transplant_day = SIMULATION['transplant_offset']
    simulation_days = SIMULATION['days']
    t_start = transplant_day * 86400
    t_end = (transplant_day + simulation_days) * 86400

    print("=" * 100)
    print("ROS 版本 v1.0 - 完整驗證結果")
    print("=" * 100)
    print(f"{'Treatment':<10} {'FW_obs':>7} {'FW_sim':>7} {'FW_Err':>7}   "
          f"{'Anth_obs':>8} {'Anth_sim':>8} {'Anth_Err':>8}   {'Stress':>7} {'ROS':>6}")
    print("-" * 100)

    fw_errs = []
    anth_errs = []

    for treatment in ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']:
        env = get_env_for_treatment(treatment)
        target = TARGETS.get(treatment, {'FW': 0, 'Anth': 0})

        sol = solve_ivp(
            uva_sun_derivatives_ros,
            (t_start, t_end),
            initial_state,
            args=(p, env),
            method='RK45',
            max_step=120,
            t_eval=np.array([t_end])
        )

        if sol.success:
            Xd_f, Cbuf_f, LAI_f, Anth_f, Stress_f, ROS_f = sol.y[:, -1]
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

            s1 = "ok" if abs(fw_err) < 5 else "NG"
            s2 = "ok" if abs(anth_err) < 10 else "NG"
            print(f"{treatment:<10} {FW_obs:>6.1f}g {FW_sim:>6.1f}g {fw_err:>+6.1f}%{s1} "
                  f"{Anth_obs:>7.1f} {Anth_sim:>7.1f} {anth_err:>+7.1f}%{s2} {Stress_f:>7.2f} {ROS_f:>6.3f}")
        else:
            print(f"{treatment:<10} Failed: {sol.message}")
            fw_errs.append(100)
            anth_errs.append(100)

    print("=" * 100)
    fw_ok = sum(1 for e in fw_errs if e < 5)
    anth_ok = sum(1 for e in anth_errs if e < 10)
    print(f"FW <5%: {fw_ok}/6, Anth <10%: {anth_ok}/6, Total: {fw_ok + anth_ok}/12")

    return fw_ok + anth_ok


# ==============================================================================
# 主程式
# ==============================================================================
if __name__ == "__main__":
    validate_all_treatments()
