#!/usr/bin/env python3
"""
按不同機制分別校準各處理組
策略：
1. L6D6: 降低基礎損傷，讓它接近目標 (當前 -8.4%)
2. L6D6-N: 用 circadian_disruption_factor 調整（夜間節律）
3. VL3D12, L6D12: 用 LAI 脆弱性調整 (n_vuln, LAI_ref_vuln)
4. H12D3: 用非線性累積調整 (stress_nonlinear_coeff, K_nonlinear)
"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TREATMENT_CONFIGS, TARGETS, ODE_SETTINGS, SIMULATION, get_env_for_treatment
from params_config import ALL_PARAMS

def test_params(param_changes, description=""):
    """測試參數組合"""
    params = ALL_PARAMS.copy()
    params.update(param_changes)
    p = UVAParams(params)

    treatments = ['CK', 'L6D6', 'L6D6-N', 'H12D3', 'VL3D12', 'L6D12']
    results = []

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

        results.append({
            'treatment': treatment_name,
            'fw': fw_per_plant,
            'target': target['FW'],
            'err': fw_err,
            'stress': Stress,
            'LAI': LAI,
        })

    return results

print('=' * 80)
print('分階段校準各處理組')
print('=' * 80)
print()

# 階段 1: 調整 L6D6 (降低基礎損傷)
print('階段 1: 調整 L6D6 (當前 -8.4% → 目標 -2%)')
print('-' * 80)

for damage_coeff in [1.6e-6, 1.4e-6, 1.2e-6]:
    results = test_params({
        'stress_damage_coeff': damage_coeff,
    })

    l6d6 = next(r for r in results if r['treatment'] == 'L6D6')
    ck = next(r for r in results if r['treatment'] == 'CK')

    print(f'  damage_coeff={damage_coeff:.1e}: CK {ck["err"]:+.1f}%, L6D6 {l6d6["err"]:+.1f}%')

print()

# 選擇最佳值
best_damage_coeff = 1.4e-6
print(f'→ 選擇 stress_damage_coeff = {best_damage_coeff:.1e}')
print()

# 階段 2: 調整 L6D6-N (夜間節律懲罰)
print('階段 2: 調整 L6D6-N 夜間節律 (當前 -3.5%，已經不錯)')
print('-' * 80)

for circadian_factor in [1.5, 2.0, 2.5]:
    results = test_params({
        'stress_damage_coeff': best_damage_coeff,
        'circadian_disruption_factor': circadian_factor,
    })

    l6d6n = next(r for r in results if r['treatment'] == 'L6D6-N')

    print(f'  circadian_factor={circadian_factor}: L6D6-N {l6d6n["err"]:+.1f}%')

print()
print('→ 保持 circadian_disruption_factor = 2.0 (L6D6-N 已經很好)')
print()

# 階段 3: 調整 VL3D12, L6D12 (LAI 脆弱性)
print('階段 3: 調整 VL3D12, L6D12 (LAI 脆弱性)')
print('-' * 80)
print('說明: 長期照射 → LAI 增長 → 脆弱性降低 → 需要提高脆弱性敏感度')
print()

for n_vuln in [7, 8, 9, 10]:
    results = test_params({
        'stress_damage_coeff': best_damage_coeff,
        'n_vuln': n_vuln,
    })

    vl3d12 = next(r for r in results if r['treatment'] == 'VL3D12')
    l6d12 = next(r for r in results if r['treatment'] == 'L6D12')
    l6d6 = next(r for r in results if r['treatment'] == 'L6D6')

    print(f'  n_vuln={n_vuln}: L6D6 {l6d6["err"]:+.1f}%, VL3D12 {vl3d12["err"]:+.1f}%, L6D12 {l6d12["err"]:+.1f}%')

print()

# 測試降低 LAI_ref_vuln (提高脆弱性)
print('測試降低 LAI_ref_vuln (讓高 LAI 也脆弱):')
for lai_ref in [7.5, 6.0, 5.0]:
    results = test_params({
        'stress_damage_coeff': best_damage_coeff,
        'LAI_ref_vuln': lai_ref,
    })

    vl3d12 = next(r for r in results if r['treatment'] == 'VL3D12')
    l6d12 = next(r for r in results if r['treatment'] == 'L6D12')
    l6d6 = next(r for r in results if r['treatment'] == 'L6D6')

    print(f'  LAI_ref={lai_ref}: L6D6 {l6d6["err"]:+.1f}%, VL3D12 {vl3d12["err"]:+.1f}% (LAI={vl3d12["LAI"]:.2f}), L6D12 {l6d12["err"]:+.1f}% (LAI={l6d12["LAI"]:.2f})')

print()
best_lai_ref = 6.0
print(f'→ 選擇 LAI_ref_vuln = {best_lai_ref}')
print()

# 階段 4: 調整 H12D3 (非線性累積)
print('階段 4: 調整 H12D3 (非線性累積傷害)')
print('-' * 80)
print('說明: 短期高強度 → 提高非線性係數讓 Stress 快速累積')
print()

for nonlinear_coeff in [2.0, 3.0, 4.0, 5.0, 6.0]:
    for k_nonlinear in [2.0, 1.5, 1.0]:
        results = test_params({
            'stress_damage_coeff': best_damage_coeff,
            'LAI_ref_vuln': best_lai_ref,
            'stress_nonlinear_coeff': nonlinear_coeff,
            'K_nonlinear': k_nonlinear,
        })

        h12d3 = next(r for r in results if r['treatment'] == 'H12D3')
        l6d6 = next(r for r in results if r['treatment'] == 'L6D6')

        if abs(h12d3['err']) < 20:  # 只顯示有希望的組合
            print(f'  nonlinear={nonlinear_coeff}, K={k_nonlinear}: '
                  f'L6D6 {l6d6["err"]:+.1f}%, H12D3 {h12d3["err"]:+.1f}% (Stress={h12d3["stress"]:.2f})')

print()
print('=' * 80)
