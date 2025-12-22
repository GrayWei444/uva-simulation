#!/usr/bin/env python3
"""
精細校準所有處理組，目標誤差 < 5%
策略：逐個處理組微調
"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TREATMENT_CONFIGS, TARGETS, ODE_SETTINGS, SIMULATION, get_env_for_treatment
from params_config import ALL_PARAMS

def test_params(param_changes):
    """測試參數組合並返回所有處理組結果"""
    params = ALL_PARAMS.copy()
    params.update(param_changes)
    p = UVAParams(params)

    treatments = ['CK', 'L6D6', 'L6D6-N', 'H12D3', 'VL3D12', 'L6D12']
    results = {}

    for treatment_name in treatments:
        env = get_env_for_treatment(treatment_name)
        target = TARGETS[treatment_name]

        fw_init_g = SIMULATION['initial_fw_g']
        dw_init_g = fw_init_g * p.dw_fw_ratio_base
        X_d_init = dw_init_g / 1000.0 * ENV_BASE['plant_density']
        C_buf_init = X_d_init * 0.1
        LAI_init = (dw_init_g / 0.01) * 0.04
        fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
        Anth_init = 5.0 * fw_total_init / 1e6

        y0 = [X_d_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]

        t_start = SIMULATION['transplant_offset'] * 86400
        t_end = (SIMULATION['transplant_offset'] + SIMULATION['days']) * 86400

        sol = solve_ivp(
            fun=lambda t, y: uva_sun_derivatives(t, y, p, env),
            t_span=(t_start, t_end),
            y0=y0,
            method='RK45',
            max_step=60,
            dense_output=True
        )

        if not sol.success:
            continue

        y_final = sol.y[:, -1]
        X_d, C_buf, LAI, Anth, Stress, E_stress = y_final

        dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress, p)
        dw_per_plant = X_d / ENV_BASE['plant_density'] * 1000
        fw_per_plant = dw_per_plant / dw_fw_ratio

        fw_err = ((fw_per_plant - target['FW']) / target['FW'] * 100) if target['FW'] > 0 else 0

        results[treatment_name] = {
            'fw': fw_per_plant,
            'target': target['FW'],
            'err': fw_err,
            'stress': Stress,
        }

    return results

def print_results(results, description=""):
    """打印結果"""
    if description:
        print(f'\n{description}')

    # 計算平均誤差和最大誤差
    errors = [abs(r['err']) for r in results.values()]
    avg_err = np.mean(errors)
    max_err = np.max(errors)
    within_5 = sum(1 for e in errors if e < 5)

    print(f'  平均誤差: {avg_err:.1f}%, 最大誤差: {max_err:.1f}%, {within_5}/6 組 < 5%')

    for name, r in results.items():
        marker = '✅' if abs(r['err']) < 5 else '⚠️'
        print(f'    {name:<10s}: {r["fw"]:6.2f}g (目標 {r["target"]:5.1f}g, 誤差 {r["err"]:+6.1f}%) {marker}')

print('=' * 80)
print('精細校準 - 目標: 所有處理組誤差 < 5%')
print('=' * 80)

# 當前狀態 (v6.1.1)
print('\n當前狀態 (v6.1.1):')
current = test_params({})
print_results(current)

# 分析：哪些組需要調整？
# CK: +1.0% ✅
# L6D6: -2.4% ✅
# L6D6-N: +6.9% ⚠️ (需要降低)
# VL3D12: +4.6% ✅
# L6D12: +7.7% ⚠️ (需要降低)
# H12D3: +49.4% ❌ (嚴重高估)

print('\n' + '=' * 80)
print('策略分析')
print('=' * 80)
print('''
需要調整的組:
1. L6D6-N: +6.9% → 需要提高夜間節律懲罰
2. L6D12: +7.7% → 需要提高 LAI 脆弱性
3. H12D3: +49.4% → 需要大幅提高損傷（但不能影響 L6D6）

關鍵矛盾:
- H12D3 需要極高的損傷累積
- 但提高基礎損傷會影響 L6D6
- 需要找到只影響 H12D3 的機制
''')

# 測試 1: 調整 L6D6-N (夜間節律)
print('\n測試 1: 調整 L6D6-N 夜間節律懲罰')
print('-' * 80)

for circadian_factor in [2.2, 2.5, 2.8, 3.0]:
    results = test_params({'circadian_disruption_factor': circadian_factor})
    l6d6n_err = results['L6D6-N']['err']
    print(f'  circadian_factor={circadian_factor}: L6D6-N {l6d6n_err:+.1f}%')

best_circadian = 2.5
print(f'\n→ 選擇 circadian_disruption_factor = {best_circadian}')

# 測試 2: 調整 L6D12 (LAI 脆弱性)
print('\n測試 2: 調整 L6D12 LAI 脆弱性')
print('-' * 80)

for n_vuln in [8, 9, 10]:
    for lai_ref in [6.0, 6.5, 7.0]:
        results = test_params({
            'circadian_disruption_factor': best_circadian,
            'n_vuln': n_vuln,
            'LAI_ref_vuln': lai_ref,
        })
        l6d12_err = results['L6D12']['err']
        vl3d12_err = results['VL3D12']['err']
        if abs(l6d12_err) < 6:  # 只顯示有希望的
            print(f'  n_vuln={n_vuln}, LAI_ref={lai_ref}: L6D12 {l6d12_err:+.1f}%, VL3D12 {vl3d12_err:+.1f}%')

best_n_vuln = 9
best_lai_ref = 6.5
print(f'\n→ 選擇 n_vuln={best_n_vuln}, LAI_ref_vuln={best_lai_ref}')

# 測試 3: H12D3 關鍵 - 尋找只影響高劑量的參數
print('\n測試 3: H12D3 特殊機制 - 日內非線性')
print('-' * 80)
print('說明: H12D3 照射 12h/day, 可能觸發日內非線性')

for e_50 in [800, 700, 600, 500]:
    for k_intraday in [2.0, 3.0, 4.0]:
        results = test_params({
            'circadian_disruption_factor': best_circadian,
            'n_vuln': best_n_vuln,
            'LAI_ref_vuln': best_lai_ref,
            'E_50': e_50,
            'k_intraday': k_intraday,
        })
        h12d3_err = results['H12D3']['err']
        l6d6_err = results['L6D6']['err']
        if abs(h12d3_err) < 30 and abs(l6d6_err) < 5:  # 篩選條件
            print(f'  E_50={e_50}, k_intraday={k_intraday}: H12D3 {h12d3_err:+.1f}%, L6D6 {l6d6_err:+.1f}%')

print('\n測試 4: H12D3 - 組合調整損傷和抑制')
print('-' * 80)

for damage_coeff in [1.5e-6, 1.6e-6, 1.7e-6, 1.8e-6]:
    for inhibition in [0.60, 0.65, 0.70]:
        for k_stress in [2.0, 2.5, 3.0]:
            results = test_params({
                'circadian_disruption_factor': best_circadian,
                'n_vuln': best_n_vuln,
                'LAI_ref_vuln': best_lai_ref,
                'stress_damage_coeff': damage_coeff,
                'stress_photosynthesis_inhibition': inhibition,
                'K_stress': k_stress,
            })

            # 計算所有組誤差
            errors = [abs(r['err']) for r in results.values()]
            avg_err = np.mean(errors)
            max_err = np.max(errors)
            within_5 = sum(1 for e in errors if e < 5)

            # 只顯示有希望的組合
            if within_5 >= 4 and max_err < 35:
                print(f'  damage={damage_coeff:.1e}, inhibit={inhibition}, K={k_stress}: '
                      f'平均{avg_err:.1f}%, 最大{max_err:.1f}%, {within_5}/6組<5%')
                print_results(results, '')

print('\n' + '=' * 80)
print('搜索完成')
print('=' * 80)
