"""測試 ROS 參數 v3 - 調整策略"""
import numpy as np
from scipy.integrate import solve_ivp
import sys
sys.path.insert(0, '/home/kasm-user/projects/uva-simulation')

from simulate_uva_model_ros import ALL_PARAMS, UVAParams, uva_sun_derivatives_ros, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment

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
    ROS_init_mM = 0.0
    initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0, ROS_init_mM]

    transplant_day = SIMULATION['transplant_offset']
    simulation_days = SIMULATION['days']
    t_start = transplant_day * 86400
    t_end = (transplant_day + simulation_days) * 86400

    print(f"\n{label}")
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

            print(f"{treatment:<10} FW:{FW_sim:>5.1f}({fw_err:>+5.1f}%) Anth:{Anth_sim:>5.1f}({anth_err:>+5.1f}%) S:{Stress_f:>5.1f} ROS:{ROS_f:>5.3f}")
        else:
            print(f"{treatment:<10} 失敗")
            fw_errs.append(100)
            anth_errs.append(100)

    fw_ok = sum(1 for e in fw_errs if e < 5)
    anth_ok = sum(1 for e in anth_errs if e < 10)
    print(f"達標: FW {fw_ok}/6, Anth {anth_ok}/6, Total {fw_ok + anth_ok}/12")
    return fw_ok + anth_ok

print("=" * 110)
print("ROS v3 - 新策略: 高 K_ros_damage 讓低 ROS 時損傷小")
print("=" * 110)

# 關鍵洞察:
# - L6D6 (6h/天) 的 ROS 穩態應該低於 K_ros_damage → 損傷放大 ~1
# - H12D3 (12h/天) 的 ROS 穩態應該高於 K_ros_damage → 損傷放大顯著

# 測試系列: 找到正確的 K_ros_damage
for K in [2.0, 3.0, 5.0, 8.0]:
    test_params({
        'ros_prod_coeff': 2e-5,
        'ros_vmax_apx': 0.0002,
        'ros_vmax_cat': 0.0005,
        'ros_passive_decay': 2e-6,
        'ros_damage_gain': 60.0,
        'K_ros_damage': K,
    }, f"K_ros_damage = {K}")

# 測試: 進一步調整生成和清除
print("\n" + "=" * 110)
print("調整生成/清除平衡")
print("=" * 110)

test_params({
    'ros_prod_coeff': 1e-5,   # 降低生成
    'ros_vmax_apx': 0.0001,   # 降低清除
    'ros_vmax_cat': 0.0002,
    'ros_passive_decay': 1e-6,
    'ros_damage_gain': 80.0,
    'K_ros_damage': 3.0,
}, "低生成低清除")

test_params({
    'ros_prod_coeff': 3e-5,   # 提高生成
    'ros_vmax_apx': 0.0003,
    'ros_vmax_cat': 0.0008,
    'ros_passive_decay': 3e-6,
    'ros_damage_gain': 50.0,
    'K_ros_damage': 4.0,
}, "高生成高清除")
