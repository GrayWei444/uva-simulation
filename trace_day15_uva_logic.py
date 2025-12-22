#!/usr/bin/env python3
"""
深入追蹤 Day 15 的 UVA 邏輯
檢查 in_uva_window, hours_elapsed, E_elapsed 的計算
"""

import numpy as np
from simulate_uva_model import UVAParams, uva_sun_derivatives
from model_config import ENV_BASE, TREATMENT_CONFIGS, SIMULATION

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

print('=' * 80)
print('Day 15 UVA 邏輯追蹤')
print('=' * 80)
print()

print('L6D6 配置:')
print(f'  uva_on = {env.get("uva_on")}')
print(f'  uva_intensity = {env.get("uva_intensity")} W/m²')
print(f'  uva_start_day = {env.get("uva_start_day")}')
print(f'  uva_end_day = {env.get("uva_end_day")}')
print(f'  uva_hour_on = {env.get("uva_hour_on")}')
print(f'  uva_hour_off = {env.get("uva_hour_off")}')
print()

# 檢查一整天的每個小時
print(f'{"Hour":>6s} {"day":>8s} {"I_UVA":>8s} {"in_window":>10s} {"hours_elapsed":>15s} {"E_elapsed":>12s}')
print('-' * 80)

for h in range(24):
    t = (14 * 86400) + (h * 3600)  # Day 14, hour h

    # 手動重現邏輯
    day_from_sowing = t / 86400.0
    hour = (t / 3600.0) % 24

    uva_on = env.get('uva_on', False)
    uva_start_day = env.get('uva_start_day', 29)
    uva_end_day = env.get('uva_end_day', 35)
    uva_hour_on = env.get('uva_hour_on', 10)
    uva_hour_off = env.get('uva_hour_off', 16)
    uva_intensity = env.get('uva_intensity', 22.0)

    I_UVA = 0.0
    in_uva_window = False

    if uva_on and uva_start_day <= day_from_sowing < uva_end_day:
        if uva_hour_on <= uva_hour_off:
            in_uva_window = uva_hour_on <= hour < uva_hour_off
        else:
            in_uva_window = hour >= uva_hour_on or hour < uva_hour_off

        if in_uva_window:
            I_UVA = uva_intensity

    # 計算 hours_elapsed (simulate_uva_model.py line 294-303)
    if in_uva_window:
        if uva_hour_on <= uva_hour_off:
            hours_elapsed = hour - uva_hour_on
        else:
            if hour >= uva_hour_on:
                hours_elapsed = hour - uva_hour_on
            else:
                hours_elapsed = (24 - uva_hour_on) + hour
    else:
        hours_elapsed = 0.0

    # BUG LINE: I_UVA_config
    I_UVA_config = env.get('I_UVA', 11.0)  # ← simulate_uva_model.py line 308
    E_elapsed = I_UVA_config * hours_elapsed * 3.6

    print(f'{h:6.0f} {day_from_sowing:8.2f} {I_UVA:8.1f} {str(in_uva_window):>10s} {hours_elapsed:15.1f} {E_elapsed:12.1f}')

print()
print('=' * 80)
print('診斷:')
print('=' * 80)
print()

print(f'env.get("I_UVA", 11.0) = {env.get("I_UVA", 11.0)}  ← BUG!')
print(f'應該使用: env.get("uva_intensity", 22.0) = {env.get("uva_intensity", 22.0)}')
print()

print('但是 Day 14 (UVA 還沒開始):')
print('  - in_uva_window = False (因為 day < uva_start_day)')
print('  - hours_elapsed = 0')
print('  - E_elapsed = I_UVA_config * 0 * 3.6 = 0')
print('  → 即使 I_UVA_config 錯誤，E_elapsed 仍是 0')
print()

print('所以 I_UVA_config bug 不是 Day 15-21 Stress 累積的原因...')
print('Bug 一定在其他地方！')

print()
print('=' * 80)
