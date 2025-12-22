#!/usr/bin/env python3
"""
Debug Day 16-17 的 Stress 躍升
從 debug_with_anth_init.py 我們知道:
  Day 15: Stress = 0.000000
  Day 18: Stress = 122.397664

現在找出 Stress 在哪一天開始累積
"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TREATMENT_CONFIGS, ODE_SETTINGS, SIMULATION

# L6D6 配置
p = UVAParams()
env = ENV_BASE.copy()
env.update(TREATMENT_CONFIGS['L6D6'])

# 初始條件 (with Anth_init)
fw_init_g = SIMULATION['initial_fw_g']
dw_init_g = fw_init_g * p.dw_fw_ratio_base
X_d_init = dw_init_g / 1000.0 * env['plant_density']
C_buf_init = X_d_init * 0.1
LAI_init = (dw_init_g / 0.01) * 0.04
fw_total_init = fw_init_g * env['plant_density'] / 1000
Anth_init = 5.0 * fw_total_init / 1e6

y0 = [X_d_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]

# 模擬到 Day 18
t_start = 14 * 86400
t_end = 18 * 86400

sol = solve_ivp(
    fun=lambda t, y: uva_sun_derivatives(t, y, p, env),
    t_span=(t_start, t_end),
    y0=y0,
    method='RK45',
    max_step=60,
    dense_output=True
)

print('=' * 80)
print('Day 14-18 每半天的 Stress 變化')
print('=' * 80)
print()

# 每半天檢查一次
test_times = []
for day in range(14, 19):
    for hour in [0, 12]:
        test_times.append((day, hour))

print(f'{"Day":<6s} {"Hour":>6s} {"Stress":>12s} {"I_UVA":>8s} {"dStress/dt":>15s}')
print('-' * 80)

for day, hour in test_times:
    t = (day * 86400) + (hour * 3600)
    if t > t_end:
        break

    y = sol.sol(t)
    X_d, C_buf, LAI, Anth, Stress, E_stress = y

    # 計算導數
    derivs = uva_sun_derivatives(t, y, p, env)
    dStress_dt = derivs[4]

    # 計算 I_UVA
    days_from_transplant = t / 86400.0
    h = (t / 3600.0) % 24

    uva_on = env.get('uva_on', False)
    uva_start_day = env.get('uva_start_day', 29)
    uva_end_day = env.get('uva_end_day', 35)
    uva_hour_on = env.get('uva_hour_on', 10)
    uva_hour_off = env.get('uva_hour_off', 16)
    uva_intensity = env.get('uva_intensity', 22.0)

    I_UVA = 0.0
    if uva_on and uva_start_day <= days_from_transplant < uva_end_day:
        if uva_hour_on <= uva_hour_off:
            in_uva_window = uva_hour_on <= h < uva_hour_off
        else:
            in_uva_window = h >= uva_hour_on or h < uva_hour_off

        if in_uva_window:
            I_UVA = uva_intensity

    print(f'{day:<6d} {hour:6d} {Stress:12.6f} {I_UVA:8.1f} {dStress_dt:15.8e}')

print()
print('=' * 80)
print('診斷:')
print('=' * 80)

# 找出 Stress 首次 > 0.01 的時間點
for day in range(14, 19):
    for h in range(0, 24, 3):  # 每3小時檢查一次
        t = (day * 86400) + (h * 3600)
        if t > t_end:
            break

        y = sol.sol(t)
        Stress = y[4]

        if Stress > 0.01:
            print(f'✗ Stress 首次 > 0.01 在 Day {day} Hour {h}')
            print(f'  Stress = {Stress:.6f}')
            print()

            # 檢查前一個時間點
            t_prev = t - (3 * 3600)
            if t_prev >= t_start:
                y_prev = sol.sol(t_prev)
                Stress_prev = y_prev[4]
                print(f'  前3小時 (Day {day} Hour {h-3}): Stress = {Stress_prev:.6f}')
                print(f'  增長: {Stress - Stress_prev:.6f}')
            break
    else:
        continue
    break

print()
print('=' * 80)
