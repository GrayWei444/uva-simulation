"""測試調整後的 ROS 參數"""
import numpy as np
from scipy.integrate import solve_ivp
import sys
sys.path.insert(0, '/home/kasm-user/projects/uva-simulation')

from simulate_uva_model_ros import ALL_PARAMS, UVAParams, uva_sun_derivatives_ros, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TARGETS, ODE_SETTINGS, SIMULATION, get_env_for_treatment

def test_params(params_override, label):
    params = ALL_PARAMS.copy()
    params.update(params_override)
    p = UVAParams(params)

    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    C_buf_init = Xd_init * 0.1
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6
    ROS_init_mM = 0.01
    initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0, ROS_init_mM]

    transplant_day = SIMULATION['transplant_offset']
    simulation_days = SIMULATION['days']
    t_start = transplant_day * 86400
    t_end = (transplant_day + simulation_days) * 86400

    print(f"\n{label}")
    print("-" * 110)
    print(f"{'Treatment':<10} {'FW_obs':>7} {'FW_sim':>7} {'FW_Err':>7}   {'Anth_obs':>8} {'Anth_sim':>8} {'Anth_Err':>8}   {'Stress':>7} {'ROS':>8}")
    print("-" * 110)

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
            max_step=60,
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

            print(f"{treatment:<10} {FW_obs:>6.1f}g {FW_sim:>6.1f}g {fw_err:>+6.1f}%   "
                  f"{Anth_obs:>7.1f} {Anth_sim:>7.1f} {anth_err:>+7.1f}%   {Stress_f:>7.2f} {ROS_f:>8.4f}")
        else:
            print(f"{treatment:<10} 模擬失敗")
            fw_errs.append(100)
            anth_errs.append(100)

    print("-" * 110)
    fw_ok = sum(1 for e in fw_errs if e < 5)
    anth_ok = sum(1 for e in anth_errs if e < 10)
    print(f"FW <5%: {fw_ok}/6, Anth <10%: {anth_ok}/6, Total: {fw_ok + anth_ok}/12")
    return fw_ok + anth_ok

print("=" * 110)
print("ROS 參數調整測試")
print("=" * 110)

# 測試 1: 提高 ros_prod_coeff
test_params({
    'ros_prod_coeff': 5e-4,  # 從 1e-5 提高 50 倍
}, "測試 1: ros_prod_coeff = 5e-4")

# 測試 2: 進一步提高
test_params({
    'ros_prod_coeff': 1e-3,
}, "測試 2: ros_prod_coeff = 1e-3")

# 測試 3: 同時調整 damage_gain
test_params({
    'ros_prod_coeff': 5e-4,
    'ros_damage_gain': 50.0,
}, "測試 3: ros_prod=5e-4, gain=50")

# 測試 4: 降低 APX Vmax (讓 ROS 更容易累積)
test_params({
    'ros_prod_coeff': 1e-4,
    'ros_vmax_apx': 0.005,  # 從 0.02 降到 0.005
    'ros_damage_gain': 30.0,
}, "測試 4: 降低 APX Vmax")

# 測試 5: 綜合調整
test_params({
    'ros_prod_coeff': 2e-4,
    'ros_vmax_apx': 0.008,
    'ros_vmax_cat': 0.02,
    'ros_damage_gain': 40.0,
    'K_ros_damage': 0.5,
}, "測試 5: 綜合調整")
