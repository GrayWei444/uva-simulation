#!/usr/bin/env python3
"""
追蹤整個生長期間的狀態變化
特別關注 C_buf, Stress, 和生長速率的動態
"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TREATMENT_CONFIGS, ODE_SETTINGS, SIMULATION

def simulate_and_track(treatment_name):
    """模擬並記錄關鍵狀態"""

    p = UVAParams()
    env = ENV_BASE.copy()
    env.update(TREATMENT_CONFIGS[treatment_name])

    # 初始條件
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    X_d_init = dw_init_g / 1000.0 * env['plant_density']
    C_buf_init = X_d_init * 0.1
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * env['plant_density'] / 1000
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

    return sol, p, env

print('=' * 80)
print('完整生長期間狀態追蹤')
print('=' * 80)
print()

# 模擬 CK 和 L6D6
print('模擬 CK 和 L6D6...')
sol_ck, p_ck, env_ck = simulate_and_track('CK')
sol_l6d6, p_l6d6, env_l6d6 = simulate_and_track('L6D6')
print('完成')
print()

# 關鍵時間點
key_days = [14, 21, 28, 29, 30, 31, 32, 35]

print('關鍵時間點比較 (CK vs L6D6):')
print('=' * 80)
print(f'{"Day":<6s} {"組別":<6s} {"X_d":>12s} {"C_buf":>14s} {"LAI":>8s} {"Stress":>10s} {"FW (g)":>10s}')
print('-' * 80)

for day in key_days:
    t = day * 86400

    # CK
    y_ck = sol_ck.sol(t)
    X_d_ck, C_buf_ck, LAI_ck, Anth_ck, Stress_ck, E_stress_ck = y_ck
    dw_fw_ck = calculate_dynamic_dw_fw_ratio(Stress_ck, p_ck)
    dw_per_plant_ck = X_d_ck / env_ck['plant_density'] * 1000
    fw_ck = dw_per_plant_ck / dw_fw_ck

    # L6D6
    y_l6d6 = sol_l6d6.sol(t)
    X_d_l6d6, C_buf_l6d6, LAI_l6d6, Anth_l6d6, Stress_l6d6, E_stress_l6d6 = y_l6d6
    dw_fw_l6d6 = calculate_dynamic_dw_fw_ratio(Stress_l6d6, p_l6d6)
    dw_per_plant_l6d6 = X_d_l6d6 / env_l6d6['plant_density'] * 1000
    fw_l6d6 = dw_per_plant_l6d6 / dw_fw_l6d6

    # 輸出
    marker = ' ← UVA' if day >= 29 else ''
    print(f'{day:<6d} {"CK":<6s} {X_d_ck:12.8f} {C_buf_ck:14.10f} {LAI_ck:8.4f} {Stress_ck:10.6f} {fw_ck:10.2f}')
    print(f'{"":6s} {"L6D6":<6s} {X_d_l6d6:12.8f} {C_buf_l6d6:14.10f} {LAI_l6d6:8.4f} {Stress_l6d6:10.6f} {fw_l6d6:10.2f}{marker}')
    print()

# 分析最終結果
print('=' * 80)
print('最終結果分析')
print('=' * 80)
print()

t_final = 35 * 86400
y_ck_final = sol_ck.sol(t_final)
y_l6d6_final = sol_l6d6.sol(t_final)

X_d_ck_final = y_ck_final[0]
C_buf_ck_final = y_ck_final[1]
Stress_ck_final = y_ck_final[4]

X_d_l6d6_final = y_l6d6_final[0]
C_buf_l6d6_final = y_l6d6_final[1]
Stress_l6d6_final = y_l6d6_final[4]

dw_fw_ck_final = calculate_dynamic_dw_fw_ratio(Stress_ck_final, p_ck)
dw_fw_l6d6_final = calculate_dynamic_dw_fw_ratio(Stress_l6d6_final, p_l6d6)

fw_ck_final = (X_d_ck_final / env_ck['plant_density'] * 1000) / dw_fw_ck_final
fw_l6d6_final = (X_d_l6d6_final / env_l6d6['plant_density'] * 1000) / dw_fw_l6d6_final

print('CK (Day 35):')
print(f'  X_d = {X_d_ck_final:.8f} kg/m²')
print(f'  C_buf = {C_buf_ck_final:.10f} kg C/m²')
print(f'  Stress = {Stress_ck_final:.6f}')
print(f'  DW:FW = {dw_fw_ck_final:.6f}')
print(f'  鮮重 = {fw_ck_final:.2f} g/plant')
print()

print('L6D6 (Day 35):')
print(f'  X_d = {X_d_l6d6_final:.8f} kg/m²')
print(f'  C_buf = {C_buf_l6d6_final:.10f} kg C/m²')
print(f'  Stress = {Stress_l6d6_final:.6f}')
print(f'  DW:FW = {dw_fw_l6d6_final:.6f}')
print(f'  鮮重 = {fw_l6d6_final:.2f} g/plant')
print()

print('對比:')
print(f'  X_d: L6D6/CK = {X_d_l6d6_final/X_d_ck_final:.4f}x ({(X_d_l6d6_final/X_d_ck_final-1)*100:+.1f}%)')
print(f'  鮮重: L6D6/CK = {fw_l6d6_final/fw_ck_final:.4f}x ({(fw_l6d6_final/fw_ck_final-1)*100:+.1f}%)')
print()

# 診斷
print('=' * 80)
print('關鍵發現')
print('=' * 80)
print()

# 檢查 C_buf
if C_buf_ck_final < 0:
    print(f'⚠️  CK C_buf 變負: {C_buf_ck_final:.10f}')

if C_buf_l6d6_final < 0:
    print(f'⚠️  L6D6 C_buf 變負: {C_buf_l6d6_final:.10f}')

# 檢查 C_buf 何時變負
print()
print('追蹤 C_buf 何時變負:')
for day in range(14, 36):
    t = day * 86400
    C_buf_ck_day = sol_ck.sol(t)[1]
    C_buf_l6d6_day = sol_l6d6.sol(t)[1]

    if C_buf_ck_day < 0 and sol_ck.sol((day-1)*86400)[1] >= 0:
        print(f'  CK: C_buf 在 Day {day} 變負')

    if C_buf_l6d6_day < 0 and sol_l6d6.sol((day-1)*86400)[1] >= 0:
        print(f'  L6D6: C_buf 在 Day {day} 變負')

print()
print('=' * 80)
