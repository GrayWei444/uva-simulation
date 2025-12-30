"""快速 ROS 測試 - 比較不同策略"""
import numpy as np
from scipy.integrate import solve_ivp
import sys
sys.path.insert(0, '/home/kasm-user/projects/uva-simulation')

from simulate_uva_model_ros import ALL_PARAMS, UVAParams, uva_sun_derivatives_ros, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment

# 只測試關鍵處理組
TREATMENTS = ['CK', 'L6D6', 'H12D3']

def test_single(params_override, label):
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
    initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]

    transplant_day = SIMULATION['transplant_offset']
    simulation_days = SIMULATION['days']
    t_start = transplant_day * 86400
    t_end = (transplant_day + simulation_days) * 86400

    print(f"\n{label}")
    print(f"{'Treat':<8} {'FW_obs':>6} {'FW_sim':>6} {'Err':>6} | {'S':>5} {'ROS':>6}")

    for treatment in TREATMENTS:
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
            FW_obs = target['FW']
            fw_err = (FW_sim - FW_obs) / FW_obs * 100
            print(f"{treatment:<8} {FW_obs:>5.1f}g {FW_sim:>5.1f}g {fw_err:>+5.1f}% | {Stress_f:>5.1f} {ROS_f:>6.3f}")
        else:
            print(f"{treatment:<8} 失敗")

print("=" * 60)
print("快速比較: L6D6 vs H12D3")
print("目標: L6D6 ≈ 91g (+0%), H12D3 ≈ 61g (-30%)")
print("=" * 60)

# 核心問題：ROS 如何讓 H12D3 有更多損傷？
# 嘗試不同 K_ros_damage

for K in [0.5, 1.0, 2.0, 5.0]:
    test_single({
        'ros_prod_coeff': 1e-5,
        'ros_vmax_apx': 0.0005,
        'ros_vmax_cat': 0.001,
        'ros_passive_decay': 5e-6,
        'ros_damage_gain': 80.0,
        'K_ros_damage': K,
    }, f"K_ros_damage = {K}")
