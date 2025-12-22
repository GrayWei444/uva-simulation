#!/usr/bin/env python3
"""
v6.1 Stress 參數校準
目標: 在保持 L6D6 和 CK 優秀表現的同時，改善高劑量組預測
"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import (
    ENV_BASE, TREATMENT_CONFIGS, TARGETS, ODE_SETTINGS, SIMULATION,
    get_env_for_treatment
)
from params_config import ALL_PARAMS

def test_params(param_changes):
    """測試參數組合"""

    # 創建修改後的參數
    params = ALL_PARAMS.copy()
    params.update(param_changes)

    p = UVAParams(params)

    results = []
    treatments = ['CK', 'L6D6', 'L6D6-N', 'H12D3', 'VL3D12', 'L6D12']

    for treatment_name in treatments:
        env = get_env_for_treatment(treatment_name)
        target = TARGETS.get(treatment_name, {'FW': 0, 'Anth': 0})

        # 初始條件
        fw_init_g = SIMULATION['initial_fw_g']
        dw_init_g = fw_init_g * p.dw_fw_ratio_base
        X_d_init = dw_init_g / 1000.0 * ENV_BASE['plant_density']
        C_buf_init = X_d_init * 0.1
        LAI_init = (dw_init_g / 0.01) * 0.04
        fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
        Anth_init = 5.0 * fw_total_init / 1e6

        y0 = [X_d_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]

        # 模擬
        simulation_days = SIMULATION['days']
        transplant_day = SIMULATION['transplant_offset']
        t_start = transplant_day * 86400
        t_end = (transplant_day + simulation_days) * 86400

        sol = solve_ivp(
            fun=lambda t, y: uva_sun_derivatives(t, y, p, env),
            t_span=(t_start, t_end),
            y0=y0,
            method=ODE_SETTINGS['method'],
            max_step=ODE_SETTINGS['max_step'],
            dense_output=True
        )

        if not sol.success:
            continue

        # 結果
        y_final = sol.y[:, -1]
        X_d, C_buf, LAI, Anth, Stress, E_stress = y_final

        dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress, p)
        dw_per_plant = X_d / ENV_BASE['plant_density'] * 1000
        fw_per_plant = dw_per_plant / dw_fw_ratio

        fw_kg_m2 = fw_per_plant * ENV_BASE['plant_density'] / 1000
        anth_ppm = Anth * 1e6 / (fw_kg_m2 + 1e-9)

        fw_err = ((fw_per_plant - target['FW']) / target['FW'] * 100) if target['FW'] > 0 else 0
        anth_err = ((anth_ppm - target['Anth']) / target['Anth'] * 100) if target['Anth'] > 0 else 0

        results.append({
            'treatment': treatment_name,
            'fw_pred': fw_per_plant,
            'fw_target': target['FW'],
            'fw_err': fw_err,
            'anth_pred': anth_ppm,
            'anth_target': target['Anth'],
            'anth_err': anth_err,
            'stress': Stress,
        })

    # 計算總誤差
    fw_errors = [abs(r['fw_err']) for r in results]
    avg_fw_err = np.mean(fw_errors)
    max_fw_err = np.max(fw_errors)

    return results, avg_fw_err, max_fw_err

# 測試不同參數組合
print('=' * 80)
print('v6.1 Stress 參數校準')
print('=' * 80)
print()

param_sets = [
    {
        'name': 'v6.1 當前',
        'changes': {},
    },
    {
        'name': 'v6.1.1 - 提高損傷係數',
        'changes': {
            'stress_damage_coeff': 1.2e-6,  # 從 0.8e-6 提高 50%
        },
    },
    {
        'name': 'v6.1.2 - 增強非線性',
        'changes': {
            'stress_nonlinear_coeff': 3.0,  # 從 2.0 提高 50%
        },
    },
    {
        'name': 'v6.1.3 - 降低 E_50',
        'changes': {
            'E_50': 600.0,  # 從 800.0 降低 25%
        },
    },
    {
        'name': 'v6.1.4 - 組合調整',
        'changes': {
            'stress_damage_coeff': 1.0e-6,  # +25%
            'stress_nonlinear_coeff': 2.5,   # +25%
            'E_50': 700.0,                   # -12.5%
        },
    },
    {
        'name': 'v6.1.5 - 激進調整',
        'changes': {
            'stress_damage_coeff': 1.5e-6,  # +87.5%
            'stress_nonlinear_coeff': 3.5,   # +75%
            'E_50': 500.0,                   # -37.5%
        },
    },
]

best_set = None
best_avg_err = float('inf')

for param_set in param_sets:
    print(f'測試: {param_set["name"]}')

    # 顯示參數變化
    if param_set['changes']:
        for key, value in param_set['changes'].items():
            print(f'  {key} = {value}')

    results, avg_err, max_err = test_params(param_set['changes'])

    print(f'  平均誤差: {avg_err:.1f}%')
    print(f'  最大誤差: {max_err:.1f}%')

    # 顯示關鍵組別
    ck_result = next(r for r in results if r['treatment'] == 'CK')
    l6d6_result = next(r for r in results if r['treatment'] == 'L6D6')
    h12d3_result = next(r for r in results if r['treatment'] == 'H12D3')

    print(f'  CK: {ck_result["fw_pred"]:.2f}g ({ck_result["fw_err"]:+.1f}%)')
    print(f'  L6D6: {l6d6_result["fw_pred"]:.2f}g ({l6d6_result["fw_err"]:+.1f}%)')
    print(f'  H12D3: {h12d3_result["fw_pred"]:.2f}g ({h12d3_result["fw_err"]:+.1f}%)')

    # 記錄最佳
    if avg_err < best_avg_err:
        best_avg_err = avg_err
        best_set = {
            'name': param_set['name'],
            'changes': param_set['changes'],
            'results': results,
            'avg_err': avg_err,
            'max_err': max_err,
        }

    print()

print('=' * 80)
print('最佳參數組合')
print('=' * 80)
print()

print(f'名稱: {best_set["name"]}')
print(f'平均誤差: {best_set["avg_err"]:.1f}%')
print(f'最大誤差: {best_set["max_err"]:.1f}%')
print()

print('參數變化:')
if best_set['changes']:
    for key, value in best_set['changes'].items():
        print(f'  {key}: {ALL_PARAMS[key]} → {value}')
else:
    print('  (無變化)')

print()
print('詳細結果:')
print(f'{"處理組":<10s} {"預測 FW":>10s} {"目標 FW":>10s} {"誤差":>8s}')
print('-' * 80)

for r in best_set['results']:
    marker = '✅' if abs(r['fw_err']) < 5 else ('✓' if abs(r['fw_err']) < 10 else '⚠️')
    print(f'{r["treatment"]:<10s} {r["fw_pred"]:10.2f} {r["fw_target"]:10.1f} {r["fw_err"]:+7.1f}% {marker}')

print()
print('=' * 80)
