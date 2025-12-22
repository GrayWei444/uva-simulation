#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
v6.0 快速參數校準 - 網格搜索
"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import get_env_for_treatment, SIMULATION, ENV_BASE, TARGETS

def simulate_all_treatments(params_dict):
    """模擬所有處理組"""
    p = UVAParams()
    p.V_max_anth = params_dict['V_max']
    p.K_stress_anth = params_dict['K_stress']
    p.n_stress_anth = params_dict['n_stress']
    p.k_deg = params_dict['k_deg']

    treatments = ['CK', 'L6D6', 'L6D6-N', 'H12D3', 'VL3D12', 'L6D12']
    results = []

    for treatment in treatments:
        env = get_env_for_treatment(treatment)

        fw_init_g = SIMULATION['initial_fw_g']
        dw_init_g = fw_init_g * p.dw_fw_ratio_base
        Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
        LAI_init = (dw_init_g / 0.01) * 0.04
        fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
        Anth_init = 5.0 * fw_total_init / 1e6
        y0 = [Xd_init, 0.0, LAI_init, Anth_init, 0.0, 0.0]

        sol = solve_ivp(
            lambda t, y: uva_sun_derivatives(t, y, p, env),
            (0, SIMULATION['days'] * 86400),
            y0,
            method='RK45',
            max_step=60
        )

        if not sol.success:
            return None

        X_d_final = sol.y[0, -1]
        Anth_final = sol.y[3, -1]
        Stress_final = sol.y[4, -1]

        dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress_final, p)
        FW_sim = (X_d_final / ENV_BASE['plant_density'] * 1000) / dw_fw_ratio
        fw_kg_m2 = FW_sim * ENV_BASE['plant_density'] / 1000
        Anth_sim = Anth_final * 1e6 / (fw_kg_m2 + 1e-9)

        target = TARGETS[treatment]
        anth_err = abs((Anth_sim - target['Anth']) / target['Anth'] * 100)

        results.append({
            'treatment': treatment,
            'Anth_sim': Anth_sim,
            'Anth_exp': target['Anth'],
            'error': anth_err
        })

    return results


print("=" * 80)
print("v6.0 快速參數校準 - 網格搜索")
print("=" * 80)

# 粗網格搜索
V_max_list = [1.5e-10, 2.0e-10, 2.5e-10, 3.0e-10]
K_list = [10.0, 15.0, 20.0, 25.0]
n_list = [1.5, 2.0, 2.5, 3.0]
k_deg_list = [2.0e-6, 2.5e-6, 3.0e-6]

best_error = 1000
best_params = None
count = 0
total = len(V_max_list) * len(K_list) * len(n_list) * len(k_deg_list)

print(f"\n總共測試 {total} 種組合...\n")

for V_max in V_max_list:
    for K in K_list:
        for n in n_list:
            for k_deg in k_deg_list:
                count += 1
                params = {
                    'V_max': V_max,
                    'K_stress': K,
                    'n_stress': n,
                    'k_deg': k_deg
                }

                results = simulate_all_treatments(params)
                if results is None:
                    continue

                mean_error = np.mean([r['error'] for r in results])

                if mean_error < best_error:
                    best_error = mean_error
                    best_params = params.copy()
                    print(f"[{count}/{total}] 新最佳! Mean={mean_error:.2f}%")
                    print(f"  V_max={V_max:.2e}, K={K:.1f}, n={n:.1f}, k_deg={k_deg:.2e}")

                if count % 20 == 0:
                    print(f"[{count}/{total}] 當前最佳 Mean={best_error:.2f}%")

print("\n" + "=" * 80)
print("校準完成！")
print("=" * 80)

print("\n最佳參數:")
print(f"  V_max_anth = {best_params['V_max']:.2e}")
print(f"  K_stress_anth = {best_params['K_stress']:.2f}")
print(f"  n_stress_anth = {best_params['n_stress']:.2f}")
print(f"  k_deg = {best_params['k_deg']:.2e}")

# 最終驗證
print("\n" + "=" * 80)
print("最終驗證:")
print("=" * 80)
results = simulate_all_treatments(best_params)

print(f"{'Treatment':<10} {'Anth_sim':>9} {'Anth_exp':>9} {'Error':>8}")
print("-" * 40)

for r in results:
    marker = '✓' if r['error'] < 5.0 else '✗'
    err_sign = '+' if r['Anth_sim'] > r['Anth_exp'] else ''
    print(f"{r['treatment']:<10} {r['Anth_sim']:>9.1f} {r['Anth_exp']:>9.1f} {err_sign}{r['error']:>6.1f}% {marker}")

print("-" * 40)
errors = [r['error'] for r in results]
print(f"Mean |Error|: {np.mean(errors):.2f}%")
print(f"Max |Error|:  {np.max(errors):.2f}%")
print(f"誤差 < 5% 的處理組: {sum(1 for e in errors if e < 5.0)}/6")
