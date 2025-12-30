"""測試 ROS 參數 v4 - 累積效應策略

核心洞察：
1. L6D6 (6h/天 × 6天 = 36h 總照射) 應該有最小損傷
2. H12D3 (12h/天 × 3天 = 36h 總照射) 應該有最大損傷
3. 兩者總照射時間相同，但 H12D3 每天照射更長

策略：讓 ROS 在長時間照射時累積更高，不是靠總時間而是靠「每天持續時間」
"""
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
    print("-" * 120)
    print(f"{'Treatment':<10} {'FW_obs':>7} {'FW_sim':>7} {'FW_Err':>7}   {'Anth_obs':>8} {'Anth_sim':>8} {'Anth_Err':>8}   {'Stress':>7} {'ROS':>8}")
    print("-" * 120)

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

            status_fw = "✓" if abs(fw_err) < 5 else "✗"
            status_anth = "✓" if abs(anth_err) < 10 else "✗"
            print(f"{treatment:<10} {FW_obs:>6.1f}g {FW_sim:>6.1f}g {fw_err:>+6.1f}%{status_fw}  "
                  f"{Anth_obs:>7.1f} {Anth_sim:>7.1f} {anth_err:>+7.1f}%{status_anth}  {Stress_f:>7.2f} {ROS_f:>8.4f}")
        else:
            print(f"{treatment:<10} 模擬失敗")
            fw_errs.append(100)
            anth_errs.append(100)

    print("-" * 120)
    fw_ok = sum(1 for e in fw_errs if e < 5)
    anth_ok = sum(1 for e in anth_errs if e < 10)
    print(f"FW <5%: {fw_ok}/6, Anth <10%: {anth_ok}/6, Total: {fw_ok + anth_ok}/12")
    return fw_ok + anth_ok

print("=" * 120)
print("ROS v4 - 測試極低清除率讓 ROS 累積")
print("=" * 120)

# 策略 1: 極低 APX Vmax，讓 ROS 在 12h 照射時累積很高
# 但 6h 照射後有 18h 恢復時間，ROS 會降回來
test_params({
    'ros_prod_coeff': 1e-5,      # 標準生成
    'ros_vmax_apx': 0.0001,      # 極低清除 (原本 0.02)
    'ros_vmax_cat': 0.0002,
    'ros_passive_decay': 1e-6,   # 極低被動衰減
    'ros_damage_gain': 100.0,    # 高放大
    'K_ros_damage': 5.0,         # 高半飽和讓低 ROS 時損傷小
}, "策略 1: 極低清除 + 高 K_ros_damage")

# 策略 2: 更激進
test_params({
    'ros_prod_coeff': 2e-5,
    'ros_vmax_apx': 5e-5,
    'ros_vmax_cat': 1e-4,
    'ros_passive_decay': 5e-7,
    'ros_damage_gain': 150.0,
    'K_ros_damage': 8.0,
}, "策略 2: 更低清除 + 更高 K_ros_damage")

# 策略 3: 嘗試平衡
test_params({
    'ros_prod_coeff': 5e-6,      # 較低生成
    'ros_vmax_apx': 2e-5,
    'ros_vmax_cat': 5e-5,
    'ros_passive_decay': 2e-7,
    'ros_damage_gain': 200.0,
    'K_ros_damage': 10.0,
}, "策略 3: 低生成低清除 + 極高 K_ros_damage")

# 策略 4: 測試 ROS 是否真的會累積
# 先看看 12h 照射下 ROS 會達到多高
print("\n" + "=" * 120)
print("策略 4: 極端測試 - 幾乎無清除")
print("=" * 120)
test_params({
    'ros_prod_coeff': 1e-5,
    'ros_vmax_apx': 1e-6,        # 幾乎無清除
    'ros_vmax_cat': 1e-6,
    'ros_passive_decay': 1e-8,
    'ros_damage_gain': 50.0,
    'K_ros_damage': 50.0,        # 非常高，讓低 ROS 時幾乎無損傷
}, "策略 4: 幾乎無清除 + 極高 K")
