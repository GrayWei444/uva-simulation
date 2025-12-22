#!/usr/bin/env python3
"""
Debug Day 15 的異常 Stress 累積
找出為什麼 Day 15 (UVA 還沒開始) 就有巨大的 Stress
"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams
from model_config import ENV_BASE, TREATMENT_CONFIGS, ODE_SETTINGS, SIMULATION

# 先模擬到 Day 15
p = UVAParams()
env = ENV_BASE.copy()
env.update(TREATMENT_CONFIGS['L6D6'])

fw_init_g = SIMULATION['initial_fw_g']
dw_init_g = fw_init_g * p.dw_fw_ratio_base
X_d_init = dw_init_g / 1000.0 * env['plant_density']
C_buf_init = X_d_init * 0.1
LAI_init = (dw_init_g / 0.01) * 0.04

y0 = [X_d_init, C_buf_init, LAI_init, 0.0, 0.0, 0.0]

t_start = 14 * 86400
t_end = 16 * 86400

sol = solve_ivp(
    fun=lambda t, y: [0]*6,  # 暫時用零導數，先不模擬
    t_span=(t_start, t_end),
    y0=y0,
    method='RK45',
    max_step=60,
    dense_output=True
)

# Day 14.5 (中午) 的狀態
t_test = 14.5 * 86400
y_test = y0  # 使用初始值

print('=' * 80)
print('Debug Day 14.5 (開始後第一個中午)')
print('=' * 80)
print()

print('狀態:')
X_d, C_buf, LAI, Anth, Stress, E_stress = y_test
print(f'  X_d = {X_d}')
print(f'  C_buf = {C_buf}')
print(f'  LAI = {LAI}')
print(f'  Stress = {Stress}')
print()

# 手動計算導數，追蹤 Stress 來源
days_from_transplant = t_test / 86400.0
hour = (t_test / 3600.0) % 24

print(f'時間:')
print(f'  days_from_transplant = {days_from_transplant}')
print(f'  hour = {hour}')
print()

# UVA 判斷
uva_on = env.get('uva_on', False)
uva_start_day = env.get('uva_start_day', 29)
uva_end_day = env.get('uva_end_day', 35)
uva_hour_on = env.get('uva_hour_on', 10)
uva_hour_off = env.get('uva_hour_off', 16)
uva_intensity = env.get('uva_intensity', 22.0)

I_UVA = 0.0

if uva_on and uva_start_day <= days_from_transplant < uva_end_day:
    print(f'✗ Bug 發現！在 Day {days_from_transplant:.2f} 觸發了 UVA')
    print(f'  uva_start_day = {uva_start_day}')
    print(f'  uva_end_day = {uva_end_day}')

    if uva_hour_on <= uva_hour_off:
        in_uva_window = uva_hour_on <= hour < uva_hour_off
    else:
        in_uva_window = hour >= uva_hour_on or hour < uva_hour_off

    if in_uva_window:
        I_UVA = uva_intensity
        print(f'  in_uva_window = True')
        print(f'  I_UVA = {I_UVA} W/m²')
else:
    print(f'✓ UVA 判斷正確: I_UVA = 0')

print()

# 損傷率
print(f'損傷率計算:')
print(f'  stress_damage_coeff = {p.stress_damage_coeff}')
print(f'  I_UVA = {I_UVA}')

if I_UVA == 0:
    print(f'  → 損傷率應該是 0')
else:
    print(f'  → 損傷率 = {p.stress_damage_coeff} × {I_UVA} × ... > 0')

print()
print('=' * 80)
print('測試: 直接調用 uva_sun_derivatives')
print('=' * 80)
print()

# 導入函數
from simulate_uva_model import uva_sun_derivatives

derivs = uva_sun_derivatives(t_test, y_test, p, env)

print(f'導數:')
print(f'  dX_d/dt = {derivs[0]:.10e}')
print(f'  dC_buf/dt = {derivs[1]:.10e}')
print(f'  dLAI/dt = {derivs[2]:.10e}')
print(f'  dAnth/dt = {derivs[3]:.10e}')
print(f'  dStress/dt = {derivs[4]:.10e}  ← 關鍵！')
print(f'  dE_stress/dt = {derivs[5]:.10e}')
print()

if derivs[4] > 1e-10:
    print(f'✗✗✗ Bug 確認！')
    print(f'    Day 14.5 的 dStress/dt = {derivs[4]:.10e} > 0')
    print(f'    但此時不應該有 UVA 照射 (uva_start_day = {uva_start_day})')
    print(f'    這表示即使 I_UVA = 0，模型仍在累積 Stress')
else:
    print(f'✓ dStress/dt ≈ 0 或 < 0 (修復)')

print()
print('=' * 80)
