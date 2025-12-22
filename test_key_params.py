"""
測試關鍵參數組合
=================
快速測試幾個精選的參數組合
"""

import numpy as np
from scipy.integrate import solve_ivp

from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TARGETS, ODE_SETTINGS, SIMULATION, get_env_for_treatment


def run_all_treatments(p):
    """運行所有處理組"""
    treatments = ['CK', 'L6D6', 'L6D6-N', 'H12D3', 'VL3D12', 'L6D12']
    results = {}

    for treatment in treatments:
        env = get_env_for_treatment(treatment)
        fw_init_g = SIMULATION['initial_fw_g']
        dw_init_g = fw_init_g * p.dw_fw_ratio_base
        Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
        LAI_init = (dw_init_g / 0.01) * 0.04
        fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
        Anth_init = 5.0 * fw_total_init / 1e6
        initial_state = [Xd_init, 0.0, LAI_init, Anth_init, 0.0, 0.0]

        sol = solve_ivp(
            fun=uva_sun_derivatives,
            t_span=(0, SIMULATION['days'] * 86400),
            y0=initial_state,
            args=(p, env),
            method=ODE_SETTINGS['method'],
            max_step=ODE_SETTINGS['max_step']
        )

        if sol.success:
            sim_stress = sol.y[4, -1]
            sim_E_stress = sol.y[5, -1]
            sim_dw_fw_ratio = calculate_dynamic_dw_fw_ratio(sim_stress, p)
            sim_dw_per_plant = sol.y[0, -1] / ENV_BASE['plant_density'] * 1000
            sim_fw_per_plant = sim_dw_per_plant / sim_dw_fw_ratio
            sim_anth_kg_m2 = sol.y[3, -1]
            sim_fw_kg_m2 = sim_fw_per_plant * ENV_BASE['plant_density'] / 1000
            sim_anth_ppm = sim_anth_kg_m2 * 1e6 / (sim_fw_kg_m2 + 1e-9)

            results[treatment] = {
                'FW_sim': sim_fw_per_plant,
                'Anth_sim': sim_anth_ppm,
                'Stress': sim_stress,
                'E_stress': sim_E_stress
            }
        else:
            results[treatment] = None

    return results


def print_results(params_name, V_max, K_E, n_hill, results):
    """顯示結果"""
    print(f"\n{'=' * 100}")
    print(f"參數組合: {params_name}")
    print(f"V_max_anth={V_max:.2e}, K_E_anth={K_E:.1f}, n_hill_anth={n_hill:.1f}")
    print(f"{'=' * 100}")

    print(f"\n{'Treatment':<10} {'FW_obs':<8} {'FW_sim':<8} {'FW_Err':<9} "
          f"{'Anth_obs':<10} {'Anth_sim':<10} {'Anth_Err':<11} "
          f"{'Stress':<8} {'E_stress':<10}")
    print("-" * 100)

    fw_errors = []
    anth_errors = []

    for treatment in ['CK', 'L6D6', 'L6D6-N', 'H12D3', 'VL3D12', 'L6D12']:
        target = TARGETS[treatment]
        r = results[treatment]

        if r:
            fw_err = (r['FW_sim'] - target['FW']) / target['FW'] * 100
            anth_err = (r['Anth_sim'] - target['Anth']) / target['Anth'] * 100

            fw_errors.append(abs(fw_err))
            anth_errors.append(abs(anth_err))

            print(f"{treatment:<10} {target['FW']:<8.1f} {r['FW_sim']:<8.1f} {fw_err:<+9.1f} "
                  f"{target['Anth']:<10.1f} {r['Anth_sim']:<10.1f} {anth_err:<+11.1f} "
                  f"{r['Stress']:<8.2f} {r['E_stress']:<10.2f}")
        else:
            print(f"{treatment:<10} 模擬失敗")

    if fw_errors:
        print("-" * 100)
        print(f"{'Mean |Err|':<10} {'':<8} {'':<8} {np.mean(fw_errors):<9.1f} "
              f"{'':<10} {'':<10} {np.mean(anth_errors):<11.1f}")
        print(f"{'Max |Err|':<10} {'':<8} {'':<8} {np.max(fw_errors):<9.1f} "
              f"{'':<10} {'':<10} {np.max(anth_errors):<11.1f}")

        # 重點關注
        vl_anth_err = abs((results['VL3D12']['Anth_sim'] - TARGETS['VL3D12']['Anth']) / TARGETS['VL3D12']['Anth'] * 100)
        l6_anth_err = abs((results['L6D12']['Anth_sim'] - TARGETS['L6D12']['Anth']) / TARGETS['L6D12']['Anth'] * 100)

        print("-" * 100)
        print(f"VL3D12 Anth 誤差: {vl_anth_err:.1f}% (v5.7: 34.1%) {'✓✓ 大幅改善' if vl_anth_err < 20 else '✓ 有改善' if vl_anth_err < 34.1 else '❌ 未改善'}")
        print(f"L6D12 Anth 誤差:  {l6_anth_err:.1f}% (v5.7: 40.0%) {'✓✓ 大幅改善' if l6_anth_err < 25 else '✓ 有改善' if l6_anth_err < 40.0 else '❌ 未改善'}")

        return np.mean(anth_errors), np.max(anth_errors), vl_anth_err, l6_anth_err
    else:
        return None, None, None, None


# 測試關鍵參數組合
param_sets = [
    ("當前參數 (v5.9 預設)", 3.0e-10, 50.0, 2.0),
    ("降低 K_E (更敏感)", 3.0e-10, 30.0, 2.0),
    ("提高 K_E (更不敏感)", 3.0e-10, 70.0, 2.0),
    ("提高 n (更陡峭)", 3.0e-10, 50.0, 3.0),
    ("降低 V_max (降低誘導)", 2.0e-10, 50.0, 2.0),
    ("提高 V_max (增強誘導)", 4.0e-10, 50.0, 2.0),
    ("組合A (低K高n)", 3.0e-10, 35.0, 2.5),
    ("組合B (高K低n)", 3.0e-10, 60.0, 1.8),
]

print("=" * 100)
print("測試關鍵參數組合 - v5.9 快速評估")
print("=" * 100)
print(f"\n總共 {len(param_sets)} 組參數需要測試")

best_score = float('inf')
best_params = None
all_scores = []

for i, (name, V_max, K_E, n_hill) in enumerate(param_sets):
    print(f"\n[{i+1}/{len(param_sets)}] 測試: {name}...")

    p = UVAParams()
    p.V_max_anth = V_max
    p.K_E_anth = K_E
    p.n_hill_anth = n_hill

    results = run_all_treatments(p)
    mean_anth, max_anth, vl_err, l6_err = print_results(name, V_max, K_E, n_hill, results)

    if mean_anth is not None:
        # 綜合評分: VL3D12 和 L6D12 各佔 40%, 平均誤差佔 20%
        score = 0.4 * vl_err + 0.4 * l6_err + 0.2 * mean_anth
        all_scores.append((score, name, V_max, K_E, n_hill, mean_anth, max_anth, vl_err, l6_err))

        print(f"綜合評分: {score:.2f} (越低越好)")

        if score < best_score:
            best_score = score
            best_params = (name, V_max, K_E, n_hill)
            print("  ⭐ 目前最佳!")

# 總結
print("\n" + "=" * 100)
print("測試完成 - 結果總結")
print("=" * 100)

if all_scores:
    # 排序
    all_scores.sort(key=lambda x: x[0])

    print("\n排名 (依綜合評分):")
    print(f"\n{'排名':<5} {'參數組合':<25} {'評分':<10} {'VL3D12':<10} {'L6D12':<10} {'Anth Mean':<12} {'Anth Max':<12}")
    print("-" * 100)

    for rank, (score, name, V_max, K_E, n_hill, mean_anth, max_anth, vl_err, l6_err) in enumerate(all_scores[:5]):
        print(f"{rank+1:<5} {name:<25} {score:<10.2f} {vl_err:<10.1f} {l6_err:<10.1f} {mean_anth:<12.1f} {max_anth:<12.1f}")

    # 最佳參數
    _, best_name, best_V, best_K, best_n, _, _, best_vl, best_l6 = all_scores[0]

    print("\n" + "=" * 100)
    print(f"✅ 最佳參數組合: {best_name}")
    print("=" * 100)
    print(f"\nV_max_anth = {best_V:.3e} kg/m²/s")
    print(f"K_E_anth = {best_K:.2f} Stress·day")
    print(f"n_hill_anth = {best_n:.2f}")
    print(f"\n改善效果:")
    print(f"  VL3D12: 34.1% → {best_vl:.1f}% ({34.1 - best_vl:+.1f}%)")
    print(f"  L6D12:  40.0% → {best_l6:.1f}% ({40.0 - best_l6:+.1f}%)")

    print("\n" + "=" * 100)
    print("請將最佳參數更新到 simulate_uva_model.py (第 166-171 行)")
    print("=" * 100)
