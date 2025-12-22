#!/usr/bin/env python3
"""
Debug L6D6 模擬，檢查光照強度是否正確傳遞
"""

import numpy as np
from simulate_uva_model import UVAParams, uva_sun_derivatives
from model_config import ENV_BASE, TREATMENT_CONFIGS

# 設置參數
p = UVAParams()
env = ENV_BASE.copy()
env.update(TREATMENT_CONFIGS['L6D6'])

# 測試時間點: Day 30 (UVA 照射期間), 12:00 (中午)
t_test = (30 * 86400) + (12 * 3600)

# 初始狀態 (隨意設定，只是為了測試)
y_test = [0.1, 0.01, 5.0, 0.0, 0.0, 0.0]

print('=' * 70)
print('Debug L6D6 - 檢查光照強度傳遞')
print('=' * 70)
print()
print(f'測試時間: Day 30, 12:00 (應該在 UVA 照射窗口內)')
print(f'env["I_day"] = {env["I_day"]} W/m²')
print(f'env["uva_intensity"] = {env.get("uva_intensity", 0)} W/m²')
print(f'par_conversion_factor = {p.par_conversion_factor}')
print()

# 在 uva_sun_derivatives 中加入 debug 輸出
# 修改函數暫時輸出中間值

# 手動計算應該得到的值
day_from_sowing = t_test / 86400
hour = (t_test / 3600) % 24

print(f'計算結果:')
print(f'  day_from_sowing = {day_from_sowing:.1f}')
print(f'  hour = {hour:.1f}')
print(f'  在 UVA 窗口內? {env["uva_start_day"]} <= {day_from_sowing:.1f} < {env["uva_end_day"]}')
print(f'  在 UVA 時段內? {env["uva_hour_on"]} <= {hour:.1f} < {env["uva_hour_off"]}')
print()

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

I_base = env['I_day']
I_gain_par = I_UVA * p.par_conversion_factor
I_effective = I_base + I_gain_par

print(f'預期結果:')
print(f'  I_UVA = {I_UVA} W/m²')
print(f'  I_base = {I_base} W/m²')
print(f'  I_gain_par = {I_gain_par} W/m²')
print(f'  I_effective = {I_effective} W/m²')
print()

if I_effective == 79:
    print('✓ 光照強度正確！應為 57 + 22 = 79 W/m²')
else:
    print(f'✗ 問題！I_effective 應為 79 W/m²，實際為 {I_effective} W/m²')
