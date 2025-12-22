#!/usr/bin/env python3
"""
分析 v6.1 的 Stress 累積情況
找出為什麼高劑量組被高估
"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import (
    ENV_BASE, TREATMENT_CONFIGS, TARGETS, ODE_SETTINGS, SIMULATION,
    get_env_for_treatment
)

def analyze_treatment(treatment_name):
    """分析單一處理組的 Stress 動態"""
    p = UVAParams()
    env = get_env_for_treatment(treatment_name)

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

    # 分析 UVA 期間的 Stress
    uva_start_day = env.get('uva_start_day', 35)
    uva_end_day = env.get('uva_end_day', 36)

    # 採樣點
    sample_days = np.linspace(uva_start_day, uva_end_day, 50)
    stress_samples = []

    for day in sample_days:
        t = day * 86400
        if t <= t_end:
            y = sol.sol(t)
            stress_samples.append(y[4])

    avg_stress = np.mean(stress_samples) if stress_samples else 0
    max_stress = np.max(stress_samples) if stress_samples else 0

    # 最終結果
    y_final = sol.y[:, -1]
    X_d, C_buf, LAI, Anth, Stress_final, E_stress = y_final

    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress_final, p)
    dw_per_plant = X_d / ENV_BASE['plant_density'] * 1000
    fw_per_plant = dw_per_plant / dw_fw_ratio

    return {
        'treatment': treatment_name,
        'avg_stress': avg_stress,
        'max_stress': max_stress,
        'final_stress': Stress_final,
        'E_stress': E_stress,
        'fw': fw_per_plant,
        'uva_intensity': env.get('uva_intensity', 0),
        'uva_hours': env.get('uva_hour_off', 0) - env.get('uva_hour_on', 0),
        'uva_days': uva_end_day - uva_start_day,
    }

# 分析所有 UVA 處理組
treatments = ['L6D6', 'L6D6-N', 'H12D3', 'VL3D12', 'L6D12']

print('=' * 100)
print('v6.1 Stress 累積分析')
print('=' * 100)
print()

results = []
for treatment in treatments:
    result = analyze_treatment(treatment)
    results.append(result)

# 表格
print(f'{"處理組":<10s} {"強度":>8s} {"時數/天":>8s} {"天數":>6s} {"總劑量":>10s} '
      f'{"平均Stress":>12s} {"最大Stress":>12s} {"最終FW":>10s} {"目標FW":>10s} {"誤差":>8s}')
print('-' * 100)

for r in results:
    total_dose = r['uva_intensity'] * r['uva_hours'] * r['uva_days']
    target_fw = TARGETS[r['treatment']]['FW']
    error = (r['fw'] - target_fw) / target_fw * 100

    print(f'{r["treatment"]:<10s} {r["uva_intensity"]:8.1f} {r["uva_hours"]:8.0f} {r["uva_days"]:6.0f} {total_dose:10.0f} '
          f'{r["avg_stress"]:12.6f} {r["max_stress"]:12.6f} {r["fw"]:10.2f} {target_fw:10.1f} {error:+7.1f}%')

print()
print('=' * 100)
print('診斷')
print('=' * 100)
print()

# 找出問題
print('觀察:')
print()

# 按總劑量排序
results_sorted = sorted(results, key=lambda x: x['uva_intensity'] * x['uva_hours'] * x['uva_days'])

print('1. 總劑量 vs 平均 Stress:')
for r in results_sorted:
    total_dose = r['uva_intensity'] * r['uva_hours'] * r['uva_days']
    print(f'   {r["treatment"]:<10s}: 總劑量 {total_dose:6.0f} → 平均 Stress {r["avg_stress"]:.6f}')

print()
print('2. Stress vs 鮮重抑制:')
for r in results_sorted:
    target_fw = TARGETS[r['treatment']]['FW']
    error = (r['fw'] - target_fw) / target_fw * 100
    marker = '✓' if abs(error) < 10 else '⚠️' if abs(error) < 20 else '❌'
    print(f'   {r["treatment"]:<10s}: Stress {r["final_stress"]:6.4f} → FW誤差 {error:+6.1f}% {marker}')

print()
print('問題:')
print('  - H12D3 (總劑量 2904) 的平均 Stress 太低')
print('  - 導致鮮重被高估 (+48.7%)')
print('  - 需要提高損傷累積，特別是高劑量/長時照射的情況')

print()
print('建議調整方向:')
print('  1. 提高 stress_damage_coeff (基礎損傷率)')
print('  2. 增強 stress_nonlinear_coeff (高 Stress 時加速累積)')
print('  3. 降低 E_50 (讓長時照射觸發日內非線性)')
print('  4. 提高 k_intraday (日內非線性放大係數)')

print()
print('=' * 100)
