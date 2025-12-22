#!/usr/bin/env python3
"""
Debug Day 15 結束時的 Stress 累積
為什麼從 Day 14.5 的 dStress/dt = 0，變成 Day 15 的 Stress = 30.77？
"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives
from model_config import ENV_BASE, TREATMENT_CONFIGS, ODE_SETTINGS, SIMULATION

# L6D6 配置
p = UVAParams()
env = ENV_BASE.copy()
env.update(TREATMENT_CONFIGS['L6D6'])

# 初始條件
fw_init_g = SIMULATION['initial_fw_g']
dw_init_g = fw_init_g * p.dw_fw_ratio_base
X_d_init = dw_init_g / 1000.0 * env['plant_density']
C_buf_init = X_d_init * 0.1
LAI_init = (dw_init_g / 0.01) * 0.04

y0 = [X_d_init, C_buf_init, LAI_init, 0.0, 0.0, 0.0]

# 模擬到 Day 15 結束
t_start = 14 * 86400
t_end = 15 * 86400

sol = solve_ivp(
    fun=lambda t, y: uva_sun_derivatives(t, y, p, env),
    t_span=(t_start, t_end),
    y0=y0,
    method='RK45',
    max_step=60,
    dense_output=True
)

print('=' * 80)
print('Debug Day 15 結束時的 Stress')
print('=' * 80)
print()

# 檢查多個時間點
test_hours = [0, 6, 10, 11, 12, 13, 14, 16, 18, 24]

print(f'{"Hour":>6s} {"days":>8s} {"I_UVA":>8s} {"Stress":>12s} {"dStress/dt":>15s}')
print('-' * 80)

for h in test_hours:
    t = t_start + h * 3600
    if t > t_end:
        continue

    y = sol.sol(t)
    X_d, C_buf, LAI, Anth, Stress, E_stress = y

    # 計算導數
    derivs = uva_sun_derivatives(t, y, p, env)
    dStress_dt = derivs[4]

    # 計算 I_UVA
    days_from_transplant = t / 86400.0
    hour = (t / 3600.0) % 24

    uva_on = env.get('uva_on', False)
    uva_start_day = env.get('uva_start_day', 29)
    uva_end_day = env.get('uva_end_day', 35)
    uva_hour_on = env.get('uva_hour_on', 10)
    uva_hour_off = env.get('uva_hour_off', 16)
    uva_intensity = env.get('uva_intensity', 22.0)

    I_UVA = 0.0
    if uva_on and uva_start_day <= days_from_transplant < uva_end_day:
        if uva_hour_on <= uva_hour_off:
            in_uva_window = uva_hour_on <= hour < uva_hour_off
        else:
            in_uva_window = hour >= uva_hour_on or hour < uva_hour_off

        if in_uva_window:
            I_UVA = uva_intensity

    print(f'{h:6.1f} {days_from_transplant:8.2f} {I_UVA:8.1f} {Stress:12.6f} {dStress_dt:15.8e}')

print()
print('=' * 80)
print('Day 15 結束時的狀態:')
print('=' * 80)

y_final = sol.sol(t_end)
X_d, C_buf, LAI, Anth, Stress, E_stress = y_final

print(f'  X_d = {X_d:.6f}')
print(f'  C_buf = {C_buf:.8f}')
print(f'  LAI = {LAI:.4f}')
print(f'  Stress = {Stress:.6f}')
print(f'  E_stress = {E_stress:.4f}')
print()

if Stress > 1.0:
    print(f'✗✗✗ Bug 確認！')
    print(f'    Day 15 的 Stress = {Stress:.6f} > 0')
    print(f'    但 I_UVA 全天應該都是 0 (uva_start_day = 29)')
    print()
    print('需要檢查:')
    print('  1. uva_sun_derivatives 的 I_UVA 計算')
    print('  2. Stress 累積機制是否有其他來源')
    print('  3. E_stress 累積是否正確')
else:
    print(f'✓ Stress = {Stress:.6f} (正常)')

print()
print('=' * 80)
