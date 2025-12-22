"""
快速校準腳本 - 使用網格搜尋
================================
快速測試不同參數組合，找到較佳的花青素參數
"""

import numpy as np
from scipy.integrate import solve_ivp
import sys

from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TARGETS, ODE_SETTINGS, SIMULATION, get_env_for_treatment


def run_simulation(treatment, p):
    """執行單一處理組的模擬"""
    env = get_env_for_treatment(treatment)

    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6
    initial_state = [Xd_init, 0.0, LAI_init, Anth_init, 0.0, 0.0]

    try:
        sol = solve_ivp(
            fun=uva_sun_derivatives,
            t_span=(0, SIMULATION['days'] * 86400),
            y0=initial_state,
            args=(p, env),
            method=ODE_SETTINGS['method'],
            max_step=ODE_SETTINGS['max_step']
        )

        if not sol.success:
            return None, None, None, None

        sim_stress = sol.y[4, -1]
        sim_E_stress = sol.y[5, -1]
        sim_dw_fw_ratio = calculate_dynamic_dw_fw_ratio(sim_stress, p)
        sim_dw_per_plant = sol.y[0, -1] / ENV_BASE['plant_density'] * 1000
        sim_fw_per_plant = sim_dw_per_plant / sim_dw_fw_ratio
        sim_anth_kg_m2 = sol.y[3, -1]
        sim_fw_kg_m2 = sim_fw_per_plant * ENV_BASE['plant_density'] / 1000
        sim_anth_ppm = sim_anth_kg_m2 * 1e6 / (sim_fw_kg_m2 + 1e-9)

        return sim_fw_per_plant, sim_anth_ppm, sim_stress, sim_E_stress
    except Exception as e:
        print(f"  警告: {treatment} 模擬失敗 - {e}")
        return None, None, None, None


def evaluate_params(V_max, K_E, n_hill):
    """評估一組參數"""
    p = UVAParams()
    p.V_max_anth = V_max
    p.K_E_anth = K_E
    p.n_hill_anth = n_hill

    treatments = ['CK', 'L6D6', 'L6D6-N', 'H12D3', 'VL3D12', 'L6D12']
    weights = {'CK': 1.0, 'L6D6': 1.0, 'L6D6-N': 1.0, 'H12D3': 1.0, 'VL3D12': 2.0, 'L6D12': 2.0}

    total_error = 0.0
    total_weight = 0.0
    results = {}

    for treatment in treatments:
        target = TARGETS[treatment]
        sim_fw, sim_anth, sim_stress, sim_E_stress = run_simulation(treatment, p)

        if sim_fw is None:
            return 1e6, None

        results[treatment] = {
            'FW_sim': sim_fw, 'Anth_sim': sim_anth,
            'Stress': sim_stress, 'E_stress': sim_E_stress
        }

        # 花青素誤差 (主要目標)
        anth_error = ((sim_anth - target['Anth']) / target['Anth']) ** 2
        # 鮮重誤差 (次要目標，確保不破壞 FW 校準)
        fw_error = ((sim_fw - target['FW']) / target['FW']) ** 2

        # 加權: 花青素 80%, 鮮重 20%
        combined_error = 0.8 * anth_error + 0.2 * fw_error
        total_error += weights[treatment] * combined_error
        total_weight += weights[treatment]

    return total_error / total_weight, results


def grid_search():
    """網格搜尋最佳參數"""

    print("=" * 80)
    print("快速網格搜尋校準 v5.9")
    print("=" * 80)

    # 定義搜尋網格 (粗網格)
    V_max_values = [2.0e-10, 3.0e-10, 4.0e-10, 5.0e-10]
    K_E_values = [30.0, 40.0, 50.0, 60.0, 80.0]
    n_hill_values = [1.5, 2.0, 2.5, 3.0]

    print(f"\n搜尋範圍:")
    print(f"  V_max_anth: {len(V_max_values)} 個值")
    print(f"  K_E_anth: {len(K_E_values)} 個值")
    print(f"  n_hill_anth: {len(n_hill_values)} 個值")
    print(f"  總組合數: {len(V_max_values) * len(K_E_values) * len(n_hill_values)}")

    best_error = float('inf')
    best_params = None
    best_results = None

    total_combinations = len(V_max_values) * len(K_E_values) * len(n_hill_values)
    current = 0

    print("\n開始搜尋...")
    print("-" * 80)

    for V_max in V_max_values:
        for K_E in K_E_values:
            for n_hill in n_hill_values:
                current += 1
                error, results = evaluate_params(V_max, K_E, n_hill)

                if error < best_error:
                    best_error = error
                    best_params = (V_max, K_E, n_hill)
                    best_results = results
                    print(f"  [{current}/{total_combinations}] 新最佳! V_max={V_max:.2e}, K_E={K_E:.1f}, n={n_hill:.1f}, Error={error:.6f}")
                elif current % 10 == 0:
                    print(f"  [{current}/{total_combinations}] 搜尋中...")

    print("-" * 80)
    print("\n" + "=" * 80)
    print("搜尋完成！")
    print("=" * 80)

    if best_params:
        V_max_best, K_E_best, n_hill_best = best_params

        print(f"\n最佳參數:")
        print(f"  V_max_anth = {V_max_best:.3e} kg/m²/s")
        print(f"  K_E_anth = {K_E_best:.2f} Stress·day")
        print(f"  n_hill_anth = {n_hill_best:.2f}")
        print(f"\n目標函數值: {best_error:.6f}")

        # 顯示詳細結果
        print("\n" + "=" * 80)
        print("所有處理組驗證結果")
        print("=" * 80)

        print(f"\n{'Treatment':<10} {'FW_obs':<8} {'FW_sim':<8} {'FW_Err':<9} "
              f"{'Anth_obs':<10} {'Anth_sim':<10} {'Anth_Err':<11} "
              f"{'Stress':<8} {'E_stress':<10}")
        print("-" * 100)

        fw_errors = []
        anth_errors = []

        for treatment in ['CK', 'L6D6', 'L6D6-N', 'H12D3', 'VL3D12', 'L6D12']:
            target = TARGETS[treatment]
            r = best_results[treatment]

            fw_err = (r['FW_sim'] - target['FW']) / target['FW'] * 100
            anth_err = (r['Anth_sim'] - target['Anth']) / target['Anth'] * 100

            fw_errors.append(abs(fw_err))
            anth_errors.append(abs(anth_err))

            print(f"{treatment:<10} {target['FW']:<8.1f} {r['FW_sim']:<8.1f} {fw_err:<+9.1f} "
                  f"{target['Anth']:<10.1f} {r['Anth_sim']:<10.1f} {anth_err:<+11.1f} "
                  f"{r['Stress']:<8.2f} {r['E_stress']:<10.2f}")

        print("-" * 100)
        print(f"{'Mean |Err|':<10} {'':<8} {'':<8} {np.mean(fw_errors):<9.1f} "
              f"{'':<10} {'':<10} {np.mean(anth_errors):<11.1f}")
        print(f"{'Max |Err|':<10} {'':<8} {'':<8} {np.max(fw_errors):<9.1f} "
              f"{'':<10} {'':<10} {np.max(anth_errors):<11.1f}")

        print("=" * 80)

        # 比較 v5.7
        print("\nv5.7 vs v5.9 比較:")
        print("-" * 40)
        print(f"{'指標':<20} {'v5.7':<10} {'v5.9':<10} {'改善':<10}")
        print("-" * 40)
        print(f"{'FW Mean Error':<20} {'2.3%':<10} {f'{np.mean(fw_errors):.1f}%':<10} {'-' if np.mean(fw_errors) > 2.3 else '✓':<10}")
        print(f"{'FW Max Error':<20} {'5.0%':<10} {f'{np.max(fw_errors):.1f}%':<10} {'-' if np.max(fw_errors) > 5.0 else '✓':<10}")
        print(f"{'Anth Mean Error':<20} {'15.6%':<10} {f'{np.mean(anth_errors):.1f}%':<10} {'✓' if np.mean(anth_errors) < 15.6 else '-':<10}")
        print(f"{'Anth Max Error':<20} {'40.0%':<10} {f'{np.max(anth_errors):.1f}%':<10} {'✓' if np.max(anth_errors) < 40.0 else '-':<10}")

        # VL3D12 和 L6D12 專項檢查
        vl_anth_err = abs((best_results['VL3D12']['Anth_sim'] - TARGETS['VL3D12']['Anth']) / TARGETS['VL3D12']['Anth'] * 100)
        l6_anth_err = abs((best_results['L6D12']['Anth_sim'] - TARGETS['L6D12']['Anth']) / TARGETS['L6D12']['Anth'] * 100)

        print(f"{'VL3D12 Anth Error':<20} {'34.1%':<10} {f'{vl_anth_err:.1f}%':<10} {'✓✓' if vl_anth_err < 20.0 else '✓' if vl_anth_err < 34.1 else '-':<10}")
        print(f"{'L6D12 Anth Error':<20} {'40.0%':<10} {f'{l6_anth_err:.1f}%':<10} {'✓✓' if l6_anth_err < 25.0 else '✓' if l6_anth_err < 40.0 else '-':<10}")
        print("=" * 80)

        return best_params, best_results
    else:
        print("\n未找到有效參數")
        return None, None


if __name__ == "__main__":
    best_params, best_results = grid_search()

    if best_params:
        V_max, K_E, n_hill = best_params
        print("\n請將以下參數更新到 simulate_uva_model.py (第 166-171 行):")
        print("-" * 60)
        print(f"self.V_max_anth = {V_max:.3e}  # kg/m²/s")
        print(f"self.K_E_anth = {K_E:.2f}  # Stress·day")
        print(f"self.n_hill_anth = {n_hill:.2f}  # Hill 係數")
        print("-" * 60)
