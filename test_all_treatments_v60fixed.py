#!/usr/bin/env python3
"""
測試所有處理組 - v6.0 Bug 修復後
"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TREATMENT_CONFIGS, TARGETS, ODE_SETTINGS, SIMULATION

def simulate_treatment(treatment_name):
    """模擬單一處理組"""

    p = UVAParams()
    env = ENV_BASE.copy()
    env.update(TREATMENT_CONFIGS[treatment_name])

    # 初始條件 (v6.0: 修復 C_buf 初始化)
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    C_buf_init = Xd_init * 0.1  # BUG FIX!
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6

    y0 = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]

    # 模擬時間 (v6.0: 從移植日開始)
    simulation_days = SIMULATION['days']
    transplant_day = SIMULATION['transplant_offset']
    t_start = transplant_day * 86400
    t_end = (transplant_day + simulation_days) * 86400

    # ODE 求解
    sol = solve_ivp(
        fun=uva_sun_derivatives,
        t_span=(t_start, t_end),
        y0=y0,
        args=(p, env),
        method=ODE_SETTINGS['method'],
        max_step=ODE_SETTINGS['max_step']
    )

    if sol.success:
        sim_stress = sol.y[4, -1]
        sim_dw_fw_ratio = calculate_dynamic_dw_fw_ratio(sim_stress, p)
        sim_dw_per_plant = sol.y[0, -1] / ENV_BASE['plant_density'] * 1000
        sim_fw_per_plant = sim_dw_per_plant / sim_dw_fw_ratio

        sim_anth_kg_m2 = sol.y[3, -1]
        sim_fw_kg_m2 = sim_fw_per_plant * ENV_BASE['plant_density'] / 1000
        sim_anth_ppm = sim_anth_kg_m2 * 1e6 / (sim_fw_kg_m2 + 1e-9)

        return sim_fw_per_plant, sim_anth_ppm, sim_stress
    else:
        return 0, 0, 0

print('=' * 80)
print('v6.0 Bug 修復後 - 所有處理組測試')
print('=' * 80)
print()
print('修復內容:')
print('  1. C_buf 初始化: 0.0 → X_d * 0.1')
print('  2. 模擬起始時間: Day 0 → Day 14 (移植日)')
print()
print('-' * 80)
print(f'{"處理組":<10s} {"預測FW":>10s} {"目標FW":>10s} {"誤差":>10s} {"預測Anth":>12s} {"目標Anth":>12s} {"誤差":>10s}')
print('-' * 80)

# 測試所有處理組
treatments = ['CK', 'L6D6', 'L6D6-N', 'H12D3', 'VL3D12', 'L6D12']

results = []
for treatment in treatments:
    target = TARGETS[treatment]
    sim_fw, sim_anth, sim_stress = simulate_treatment(treatment)

    fw_err = ((sim_fw - target['FW']) / target['FW'] * 100) if target['FW'] > 0 else 0
    anth_err = ((sim_anth - target['Anth']) / target['Anth'] * 100) if target['Anth'] > 0 else 0

    print(f'{treatment:<10s} {sim_fw:10.1f} {target["FW"]:10.1f} {fw_err:+9.1f}% {sim_anth:12.1f} {target["Anth"]:12.1f} {anth_err:+9.1f}%')

    results.append({
        'treatment': treatment,
        'fw_err': fw_err,
        'anth_err': anth_err
    })

print('-' * 80)
print()

# 統計
fw_errors = [abs(r['fw_err']) for r in results]
anth_errors = [abs(r['anth_err']) for r in results]

print('統計:')
print('-' * 80)
print(f'  鮮重誤差: 平均 {np.mean(fw_errors):.1f}%, 最大 {np.max(fw_errors):.1f}%')
print(f'  花青素誤差: 平均 {np.mean(anth_errors):.1f}%, 最大 {np.max(anth_errors):.1f}%')
print()

# 分析
if np.mean(fw_errors) < 10:
    print('✓ 鮮重預測良好 (平均誤差 < 10%)')
else:
    print('✗ 鮮重預測仍需改進')

if np.mean(anth_errors) < 20:
    print('✓ 花青素預測可接受 (平均誤差 < 20%)')
else:
    print('⚠️  花青素預測需要校準')

print()
print('=' * 80)
