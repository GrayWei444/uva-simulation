"""ROS v7 - 最終測試版本

回顧問題：
1. ROS Michaelis-Menten 動力學無法區分 L6D6 vs H12D3
2. 因為兩者總照射時間相同 (36h)
3. 自催化和耗竭策略也無法解決

新思路：也許 ROS 版本需要不同的損傷機制
- 不是用 ROS 來放大損傷
- 而是直接讓損傷累積依賴於連續照射時間

或者：接受 v6.9 的成功，ROS 版本暫時擱置
"""
import numpy as np
from scipy.integrate import solve_ivp
import sys
sys.path.insert(0, '/home/kasm-user/projects/uva-simulation')

from simulate_uva_model_ros import ALL_PARAMS, UVAParams, uva_sun_derivatives_ros, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment

def test_all(params_override, label):
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
            max_step=120,  # 增加 max_step 加速
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

            s1 = "✓" if abs(fw_err) < 5 else "✗"
            s2 = "✓" if abs(anth_err) < 10 else "✗"
            print(f"{treatment:<10} FW:{FW_sim:>5.1f}g({fw_err:>+5.1f}%{s1}) Anth:{Anth_sim:>5.1f}({anth_err:>+5.1f}%{s2}) S:{Stress_f:>5.2f} R:{ROS_f:>6.3f}")
        else:
            print(f"{treatment:<10} 失敗")
            fw_errs.append(100)
            anth_errs.append(100)

    fw_ok = sum(1 for e in fw_errs if e < 5)
    anth_ok = sum(1 for e in anth_errs if e < 10)
    print(f"達標: FW {fw_ok}/6, Anth {anth_ok}/6, Total {fw_ok + anth_ok}/12")
    return fw_ok + anth_ok

print("=" * 100)
print("ROS v7 - 最終嘗試")
print("=" * 100)

# 回到原始 ROS 參數看基線
test_all({}, "原始 ROS 參數")

# 嘗試極高 K_ros_damage 讓 ROS 幾乎不影響損傷
# 這樣就回到了 v6.9 的行為（沒有 ROS 影響）
test_all({
    'ros_damage_gain': 1.0,
    'K_ros_damage': 100.0,
}, "禁用 ROS 放大")

# 如果禁用 ROS 後 L6D6 正常，問題就是 ROS 放大太強
# 嘗試非常溫和的 ROS 放大
test_all({
    'ros_prod_coeff': 1e-5,
    'ros_damage_gain': 5.0,
    'K_ros_damage': 10.0,
}, "溫和 ROS: gain=5, K=10")
