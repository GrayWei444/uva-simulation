"""
測試 ROS 動態對模型的實際貢獻

問題：simulate_uva_model_ros.py 達到 12/12，但 ROS 狀態接近 0
這意味著 ROS 動態可能沒有實際貢獻

測試：
1. 完整 ROS 模型
2. 關閉 ROS 放大 (ros_damage_gain = 0)
3. 比較結果
"""

import numpy as np
from scipy.integrate import solve_ivp
from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment
from simulate_uva_model_ros import UVAParams, uva_sun_derivatives_ros, calculate_dynamic_dw_fw_ratio, ALL_PARAMS


def simulate_treatment(treatment, params):
    """模擬單一處理組"""
    p = UVAParams(params)
    env = get_env_for_treatment(treatment)

    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    # 6 state variables for ROS model
    initial_state = [Xd_init, Xd_init * 0.1, (dw_init_g / 0.01) * 0.04,
                     5.0 * fw_init_g * ENV_BASE['plant_density'] / 1000 / 1e6, 0.0, 0.0]
    t_start = SIMULATION['transplant_offset'] * 86400
    t_end = (SIMULATION['transplant_offset'] + SIMULATION['days']) * 86400

    sol = solve_ivp(uva_sun_derivatives_ros, (t_start, t_end), initial_state,
                    args=(p, env), method='RK45', max_step=3600, t_eval=[t_end])

    if sol.success:
        Xd_f, _, _, Anth_f, Stress_f, ROS_f = sol.y[:, -1]
        dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress_f, p)
        FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
        FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
        Anth_sim = Anth_f / FW_total_kg * 1e6
        return FW_sim, Anth_sim, Stress_f, ROS_f
    return None, None, None, None


print("=" * 80)
print("測試 ROS 動態對模型的實際貢獻")
print("=" * 80)

# 測試 1：完整 ROS 模型
print("\n1. 完整 ROS 模型 (ros_damage_gain = 1.0):")
print("-" * 80)
params_full = ALL_PARAMS.copy()
for t in ['CK', 'L6D6', 'H12D3']:
    FW, Anth, Stress, ROS = simulate_treatment(t, params_full)
    target = TARGETS[t]
    fw_err = (FW - target['FW']) / target['FW'] * 100
    anth_err = (Anth - target['Anth']) / target['Anth'] * 100
    print(f"{t:<8} FW:{FW:>5.1f}g({fw_err:>+5.1f}%) Anth:{Anth:>5.1f}({anth_err:>+5.1f}%) "
          f"Stress:{Stress:>5.1f} ROS:{ROS:.4f}")

# 測試 2：關閉 ROS 放大
print("\n2. 關閉 ROS 放大 (ros_damage_gain = 0.0):")
print("-" * 80)
params_no_ros = ALL_PARAMS.copy()
params_no_ros['ros_damage_gain'] = 0.0
for t in ['CK', 'L6D6', 'H12D3']:
    FW, Anth, Stress, ROS = simulate_treatment(t, params_no_ros)
    target = TARGETS[t]
    fw_err = (FW - target['FW']) / target['FW'] * 100
    anth_err = (Anth - target['Anth']) / target['Anth'] * 100
    print(f"{t:<8} FW:{FW:>5.1f}g({fw_err:>+5.1f}%) Anth:{Anth:>5.1f}({anth_err:>+5.1f}%) "
          f"Stress:{Stress:>5.1f} ROS:{ROS:.4f}")

# 測試 3：完全關閉 ROS 生成
print("\n3. 完全關閉 ROS 生成 (ros_prod_coeff = 0):")
print("-" * 80)
params_no_ros_prod = ALL_PARAMS.copy()
params_no_ros_prod['ros_prod_coeff'] = 0.0
params_no_ros_prod['ros_damage_gain'] = 0.0
for t in ['CK', 'L6D6', 'H12D3']:
    FW, Anth, Stress, ROS = simulate_treatment(t, params_no_ros_prod)
    target = TARGETS[t]
    fw_err = (FW - target['FW']) / target['FW'] * 100
    anth_err = (Anth - target['Anth']) / target['Anth'] * 100
    print(f"{t:<8} FW:{FW:>5.1f}g({fw_err:>+5.1f}%) Anth:{Anth:>5.1f}({anth_err:>+5.1f}%) "
          f"Stress:{Stress:>5.1f} ROS:{ROS:.4f}")

print("\n" + "=" * 80)
print("結論")
print("=" * 80)
print("""
如果三種測試結果幾乎相同，說明：
  → ROS 動態在 simulate_uva_model_ros.py 中沒有實際貢獻
  → 真正起作用的是 depletion_factor (來自 v6.9)
  → ROS 只是「裝飾性」的狀態變量

這解釋了為什麼 ROS 欄位總是接近 0。
""")

# 測試 4：增加 ROS 的貢獻
print("\n4. 增加 ROS 貢獻 (ros_damage_gain = 10.0, ros_prod_coeff = 1e-4):")
print("-" * 80)
params_more_ros = ALL_PARAMS.copy()
params_more_ros['ros_damage_gain'] = 10.0
params_more_ros['ros_prod_coeff'] = 1e-4  # 10x
for t in ['CK', 'L6D6', 'H12D3']:
    FW, Anth, Stress, ROS = simulate_treatment(t, params_more_ros)
    if FW is not None:
        target = TARGETS[t]
        fw_err = (FW - target['FW']) / target['FW'] * 100
        anth_err = (Anth - target['Anth']) / target['Anth'] * 100
        ok = "✓" if abs(fw_err) < 5 and abs(anth_err) < 10 else "✗"
        print(f"{t:<8} FW:{FW:>5.1f}g({fw_err:>+5.1f}%) Anth:{Anth:>5.1f}({anth_err:>+5.1f}%) "
              f"Stress:{Stress:>5.1f} ROS:{ROS:.4f} {ok}")
    else:
        print(f"{t:<8} Failed")
