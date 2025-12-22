#!/usr/bin/env python3
"""
檢查 CK 組的 Stress 累積
CK 組沒有 UVA，Stress 應該是 0 才對！
"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives
from model_config import ENV_BASE, TREATMENT_CONFIGS, SIMULATION

def simulate_treatment(treatment_name):
    """模擬處理組"""

    p = UVAParams()
    env = ENV_BASE.copy()
    env.update(TREATMENT_CONFIGS[treatment_name])

    plant_density = env['plant_density']

    # 初始條件
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * 0.05
    X_d_init = dw_init_g / 1000.0 * plant_density
    C_buf_init = X_d_init * 0.1
    LAI_init = 0.5

    y0 = [X_d_init, C_buf_init, LAI_init, 0.0, 0.0, 0.0]

    # 模擬時間
    t_start = 14 * 86400
    t_end = 35 * 86400

    # ODE 求解
    sol = solve_ivp(
        fun=lambda t, y: uva_sun_derivatives(t, y, p, env),
        t_span=(t_start, t_end),
        y0=y0,
        method='RK45',
        max_step=60,
        dense_output=True
    )

    return sol

print('=' * 80)
print('檢查 CK 組的 Stress 累積')
print('=' * 80)
print()

print('CK 組配置:')
ck_env = ENV_BASE.copy()
ck_env.update(TREATMENT_CONFIGS['CK'])
print(f'  I_day = {ck_env["I_day"]} W/m²')
print(f'  uva_on = {ck_env.get("uva_on", False)}')
print(f'  uva_intensity = {ck_env.get("uva_intensity", 0)} W/m²')
print()

print('模擬 CK...')
sol_ck = simulate_treatment('CK')
print('完成')
print()

# 關鍵時間點
key_times = {
    'Day 14 (開始)': 14 * 86400,
    'Day 21': 21 * 86400,
    'Day 28': 28 * 86400,
    'Day 35 (結束)': 35 * 86400,
}

print('CK 組關鍵時間點的狀態:')
print('-' * 80)
print(f'{"時間點":<20s} {"X_d":>12s} {"C_buf":>12s} {"LAI":>8s} {"Stress":>12s}')
print('-' * 80)

for time_name, t in key_times.items():
    y = sol_ck.sol(t)
    X_d, C_buf, LAI, Anth, Stress, E_stress = y
    print(f'{time_name:<20s} {X_d:12.6f} {C_buf:12.8f} {LAI:8.4f} {Stress:12.8f}')

print()

# 檢查 Day 30 的導數
t_test = (30 * 86400) + (12 * 3600)
y_test = sol_ck.sol(t_test)

X_d, C_buf, LAI, Anth, Stress, E_stress = y_test

print('Day 30 12:00 的狀態:')
print('-' * 80)
print(f'  X_d     = {X_d:.8f} kg/m²')
print(f'  C_buf   = {C_buf:.10f} kg C/m²')
print(f'  LAI     = {LAI:.6f}')
print(f'  Stress  = {Stress:.8f}')
print()

# 計算 UVA 強度
p = UVAParams()
env = ENV_BASE.copy()
env.update(TREATMENT_CONFIGS['CK'])

days_from_transplant = t_test / 86400.0
hour = (t_test / 3600.0) % 24

uva_on = env.get('uva_on', False)
uva_start_day = env.get('uva_start_day', 29)
uva_end_day = env.get('uva_end_day', 35)
uva_hour_on = env.get('uva_hour_on', 10)
uva_hour_off = env.get('uva_hour_off', 16)
uva_intensity = env.get('uva_intensity', 0)

I_UVA = 0.0
if uva_on and uva_start_day <= days_from_transplant < uva_end_day:
    if uva_hour_on <= hour < uva_hour_off:
        I_UVA = uva_intensity

print('UVA 照射檢查:')
print('-' * 80)
print(f'  uva_on = {uva_on}')
print(f'  days_from_transplant = {days_from_transplant:.2f}')
print(f'  hour = {hour:.2f}')
print(f'  I_UVA = {I_UVA} W/m²')
print()

if Stress > 0.01:
    print(f'✗✗✗ 問題發現！CK 組的 Stress = {Stress:.4f}')
    print(f'    CK 組完全沒有 UVA (I_UVA = {I_UVA})，Stress 應該是 0！')
    print()
    print(f'可能原因:')
    print(f'  1. Stress 初始值不是 0')
    print(f'  2. 損傷率計算有問題 (即使 I_UVA=0 也累積 Stress)')
    print(f'  3. 修復率太慢，無法抵消微小的數值誤差')
else:
    print(f'✓ CK 組 Stress = {Stress:.8f} ≈ 0')

print()
print('=' * 80)
