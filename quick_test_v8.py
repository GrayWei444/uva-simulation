"""快速測試 v8.0 參數"""
import numpy as np
from scipy.integrate import solve_ivp
from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment
from simulate_uva_model_v8_ros import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio, ALL_PARAMS

def test_params(params_override):
    params = ALL_PARAMS.copy()
    params.update(params_override)
    p = UVAParams(params)

    results = []
    for t in ['CK', 'L6D6', 'H12D3']:  # 只測試關鍵組
        env = get_env_for_treatment(t)
        target = TARGETS[t]

        fw_init_g = SIMULATION['initial_fw_g']
        dw_init_g = fw_init_g * p.dw_fw_ratio_base
        Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
        initial_state = [Xd_init, Xd_init * 0.1, (dw_init_g / 0.01) * 0.04,
                         5.0 * fw_init_g * ENV_BASE['plant_density'] / 1000 / 1e6, 0.0]
        t_start = SIMULATION['transplant_offset'] * 86400
        t_end = (SIMULATION['transplant_offset'] + SIMULATION['days']) * 86400

        sol = solve_ivp(uva_sun_derivatives, (t_start, t_end), initial_state,
                        args=(p, env), method='RK45', max_step=3600, t_eval=[t_end])

        if sol.success:
            Xd_f, _, _, Anth_f, Stress_f = sol.y[:, -1]
            dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress_f, p)
            FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
            fw_err = (FW_sim - target['FW']) / target['FW'] * 100
            results.append((t, FW_sim, fw_err, Stress_f))

    return results

# 測試不同 stress_damage_coeff
print("測試 stress_damage_coeff:")
print("-" * 60)
for sdc in [1e-8, 5e-8, 1e-7, 2e-7, 5e-7]:
    results = test_params({'stress_damage_coeff': sdc})
    print(f"sdc={sdc:.0e}:", end=" ")
    for t, fw, err, stress in results:
        print(f"{t}:{fw:.1f}g({err:+.0f}%,S={stress:.1f})", end=" ")
    print()

print("\n" + "=" * 60)
print("問題分析:")
print("""
v8.0 的問題：
- day_factor 在 6h 時已經約 4，導致 L6D6 累積大量 Stress
- 而 v7.0 的 day_factor = 3.8 在整體模型中效果較溫和

根本原因：
- v7.0 的 day_factor 乘在 stress_damage_coeff 上
- 但 v7.0 已經針對其他參數做了校準
- v8.0 直接用相同的 day_factor 形狀會產生不同結果

解決方案：
1. 大幅降低 stress_damage_coeff
2. 或者改變 v8.0 的設計，使其在小 factor 範圍工作
""")

# 嘗試極低的 stress_damage_coeff
print("\n極低 stress_damage_coeff 測試:")
print("-" * 60)
for sdc in [1e-9, 5e-9, 1e-8, 2e-8]:
    results = test_params({'stress_damage_coeff': sdc})
    print(f"sdc={sdc:.0e}:", end=" ")
    for t, fw, err, stress in results:
        ok = "✓" if abs(err) < 5 else "✗"
        print(f"{t}:{fw:.1f}g({err:+.0f}%{ok})", end=" ")
    print()
