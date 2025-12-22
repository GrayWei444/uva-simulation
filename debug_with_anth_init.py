#!/usr/bin/env python3
"""
Debug with Anth_init (與 trace_full_period.py 相同的初始條件)
"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TREATMENT_CONFIGS, ODE_SETTINGS, SIMULATION

# L6D6 配置
p = UVAParams()
env = ENV_BASE.copy()
env.update(TREATMENT_CONFIGS['L6D6'])

# 初始條件 (與 trace_full_period.py 完全相同)
fw_init_g = SIMULATION['initial_fw_g']
dw_init_g = fw_init_g * p.dw_fw_ratio_base
X_d_init = dw_init_g / 1000.0 * env['plant_density']
C_buf_init = X_d_init * 0.1
LAI_init = (dw_init_g / 0.01) * 0.04
fw_total_init = fw_init_g * env['plant_density'] / 1000
Anth_init = 5.0 * fw_total_init / 1e6  # ← 關鍵差異

y0 = [X_d_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]

print('=' * 80)
print('Debug with Anth_init')
print('=' * 80)
print()
print('初始條件:')
print(f'  X_d_init = {X_d_init}')
print(f'  C_buf_init = {C_buf_init}')
print(f'  LAI_init = {LAI_init}')
print(f'  Anth_init = {Anth_init} ← 與之前的 0.0 不同')
print()

# 模擬到 Day 21
t_start = 14 * 86400
t_end = 21 * 86400

sol = solve_ivp(
    fun=lambda t, y: uva_sun_derivatives(t, y, p, env),
    t_span=(t_start, t_end),
    y0=y0,
    method='RK45',
    max_step=60,
    dense_output=True
)

# 檢查 Day 15, 18, 21
test_days = [14, 15, 18, 21]

print(f'{"Day":<6s} {"X_d":>12s} {"C_buf":>14s} {"LAI":>8s} {"Anth":>12s} {"Stress":>10s} {"FW (g)":>10s}')
print('-' * 80)

for day in test_days:
    t = day * 86400
    y = sol.sol(t)
    X_d, C_buf, LAI, Anth, Stress, E_stress = y

    dw_fw = calculate_dynamic_dw_fw_ratio(Stress, p)
    dw_per_plant = X_d / env['plant_density'] * 1000
    fw = dw_per_plant / dw_fw

    print(f'{day:<6d} {X_d:12.8f} {C_buf:14.10f} {LAI:8.4f} {Anth:12.8e} {Stress:10.6f} {fw:10.2f}')

print()
print('=' * 80)
print('診斷:')
print('=' * 80)

y_day21 = sol.sol(21 * 86400)
Stress_day21 = y_day21[4]

if Stress_day21 > 1.0:
    print(f'✗✗✗ Bug 確認！')
    print(f'    Day 21 的 Stress = {Stress_day21:.6f}')
    print(f'    這與 trace_full_period.py 的結果一致')
    print()
    print('下一步: 檢查為何 Anth_init 會影響 Stress 累積')
else:
    print(f'✓ Stress = {Stress_day21:.6f} (正常)')
    print('  與之前 Anth_init = 0 的測試一致')

print()
print('=' * 80)
