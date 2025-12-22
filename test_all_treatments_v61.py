#!/usr/bin/env python3
"""
測試所有處理組 (v6.1 - Bug Fix: day_from_sowing 計算錯誤)
"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import (
    ENV_BASE, TREATMENT_CONFIGS, TARGETS, ODE_SETTINGS, SIMULATION,
    get_env_for_treatment
)

def simulate_treatment(treatment_name):
    """模擬單一處理組"""
    p = UVAParams()
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

    # 模擬時間
    simulation_days = SIMULATION['days']
    transplant_day = SIMULATION['transplant_offset']
    t_start = transplant_day * 86400
    t_end = (transplant_day + simulation_days) * 86400

    # ODE 求解
    sol = solve_ivp(
        fun=lambda t, y: uva_sun_derivatives(t, y, p, env),
        t_span=(t_start, t_end),
        y0=y0,
        method=ODE_SETTINGS['method'],
        max_step=ODE_SETTINGS['max_step'],
        dense_output=True
    )

    if not sol.success:
        return None

    # 最終結果
    y_final = sol.y[:, -1]
    X_d, C_buf, LAI, Anth, Stress, E_stress = y_final

    # 鮮重
    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress, p)
    dw_per_plant = X_d / ENV_BASE['plant_density'] * 1000
    fw_per_plant = dw_per_plant / dw_fw_ratio

    # 花青素
    fw_kg_m2 = fw_per_plant * ENV_BASE['plant_density'] / 1000
    anth_ppm = Anth * 1e6 / (fw_kg_m2 + 1e-9)

    # 誤差
    fw_err = ((fw_per_plant - target['FW']) / target['FW'] * 100) if target['FW'] > 0 else 0
    anth_err = ((anth_ppm - target['Anth']) / target['Anth'] * 100) if target['Anth'] > 0 else 0

    return {
        'treatment': treatment_name,
        'fw_pred': fw_per_plant,
        'fw_target': target['FW'],
        'fw_err': fw_err,
        'anth_pred': anth_ppm,
        'anth_target': target['Anth'],
        'anth_err': anth_err,
        'stress': Stress,
        'e_stress': E_stress,
    }

# 測試所有處理組
treatments = ['CK', 'L6D6', 'L6D6-N', 'H12D3', 'VL3D12', 'L6D12']

print('=' * 80)
print('v6.1 所有處理組預測結果')
print('=' * 80)
print()

results = []
for treatment in treatments:
    result = simulate_treatment(treatment)
    if result:
        results.append(result)

# 輸出表格
print(f'{"處理組":<10s} {"預測 FW":>10s} {"目標 FW":>10s} {"誤差":>8s} {"預測 Anth":>12s} {"目標 Anth":>12s} {"誤差":>8s}')
print('-' * 80)

for r in results:
    print(f'{r["treatment"]:<10s} {r["fw_pred"]:10.2f} {r["fw_target"]:10.1f} {r["fw_err"]:+7.1f}% '
          f'{r["anth_pred"]:12.1f} {r["anth_target"]:12.1f} {r["anth_err"]:+7.1f}%')

print()
print('=' * 80)
print('統計摘要')
print('=' * 80)

fw_errors = [abs(r['fw_err']) for r in results]
anth_errors = [abs(r['anth_err']) for r in results]

print(f'  鮮重平均誤差: {np.mean(fw_errors):.1f}%')
print(f'  鮮重最大誤差: {np.max(fw_errors):.1f}%')
print(f'  花青素平均誤差: {np.mean(anth_errors):.1f}%')
print(f'  花青素最大誤差: {np.max(anth_errors):.1f}%')

print()
print('=' * 80)
print('與 v6.0 對比')
print('=' * 80)

# v6.0 結果 (from V60_FINAL_CALIBRATION_REPORT.md)
v60_fw_errors = {
    'CK': 1.0,
    'L6D6': 18.3,
    'L6D6-N': 9.2,
    'H12D3': 11.8,
    'VL3D12': 11.5,
    'L6D12': 21.1,
}

print()
print('鮮重誤差改善:')
for r in results:
    t = r['treatment']
    v60_err = v60_fw_errors.get(t, 0)
    v61_err = abs(r['fw_err'])
    improvement = v60_err - v61_err
    marker = '✅' if improvement > 5 else ('✓' if improvement > 0 else '❌')
    print(f'  {t:<10s}: v6.0 {v60_err:5.1f}% → v6.1 {v61_err:5.1f}% (改善 {improvement:+5.1f}%) {marker}')

print()
print('=' * 80)
