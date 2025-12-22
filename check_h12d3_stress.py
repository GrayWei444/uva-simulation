#!/usr/bin/env python3
"""檢查 H12D3 的 Stress 累積"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, get_env_for_treatment, TARGETS, SIMULATION

p = UVAParams()
env = get_env_for_treatment('H12D3')

# 初始條件
fw_init_g = SIMULATION['initial_fw_g']
dw_init_g = fw_init_g * p.dw_fw_ratio_base
X_d_init = dw_init_g / 1000.0 * ENV_BASE['plant_density']
C_buf_init = X_d_init * 0.1
LAI_init = (dw_init_g / 0.01) * 0.04
fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
Anth_init = 5.0 * fw_total_init / 1e6

y0 = [X_d_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]

# 模擬
t_start = SIMULATION['transplant_offset'] * 86400
t_end = (SIMULATION['transplant_offset'] + SIMULATION['days']) * 86400

sol = solve_ivp(
    fun=lambda t, y: uva_sun_derivatives(t, y, p, env),
    t_span=(t_start, t_end),
    y0=y0,
    method='RK45',
    max_step=60,
    dense_output=True
)

print('H12D3 配置:')
print(f'  UVA: 22 W/m², 12h/day, Day 32-35 (3 days)')
print(f'  目標鮮重: 60.6g (比 CK 的 87g 低 30.3%)')
print()

# 關鍵時間點
days = [28, 31, 32, 33, 34, 35]

print(f'{"Day":<6s} {"Stress":>10s} {"FW":>10s} {"說明":<30s}')
print('-' * 60)

for day in days:
    t = day * 86400
    y = sol.sol(t)
    X_d, C_buf, LAI, Anth, Stress, E_stress = y

    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress, p)
    dw_per_plant = X_d / ENV_BASE['plant_density'] * 1000
    fw = dw_per_plant / dw_fw_ratio

    note = ''
    if day == 28:
        note = 'UVA 前'
    elif day == 32:
        note = 'UVA 開始'
    elif day == 35:
        note = 'UVA 結束/模擬結束'

    print(f'{day:<6d} {Stress:10.4f} {fw:10.2f} {note:<30s}')

y_final = sol.y[:, -1]
Stress_final = y_final[4]

print()
print('分析:')
print(f'  最終 Stress = {Stress_final:.4f}')
print(f'  預測鮮重 = {fw:.2f}g')
print(f'  目標鮮重 = 60.6g')
print(f'  誤差 = {(fw - 60.6) / 60.6 * 100:+.1f}%')
print()

# 估算需要多少 Stress 才能達到目標
print('反推: 要達到 60.6g 需要多少 Stress?')
print()

# CK 的預測鮮重是 87.8g
# H12D3 應該是 60.6g，相當於 -31% 抑制
# 使用公式反推 Stress

ck_fw = 87.8
target_fw = 60.6
inhibition_ratio = (ck_fw - target_fw) / ck_fw  # 0.31

print(f'  CK 鮮重: {ck_fw}g')
print(f'  目標鮮重: {target_fw}g')
print(f'  需要抑制: {inhibition_ratio * 100:.1f}%')
print()

# 抑制公式: inhibition = stress_photosynthesis_inhibition * Stress / (K_stress + Stress)
# 假設 inhibition = 0.31, stress_photosynthesis_inhibition = 0.50, K_stress = 3.0
# 0.31 = 0.50 * Stress / (3.0 + Stress)
# 0.31 * (3.0 + Stress) = 0.50 * Stress
# 0.93 + 0.31 * Stress = 0.50 * Stress
# 0.93 = 0.19 * Stress
# Stress = 0.93 / 0.19 = 4.89

inhibition_coeff = p.stress_photosynthesis_inhibition
K_stress = p.K_stress

required_stress = (inhibition_ratio * K_stress) / (inhibition_coeff - inhibition_ratio)

print(f'  當前參數:')
print(f'    stress_photosynthesis_inhibition = {inhibition_coeff}')
print(f'    K_stress = {K_stress}')
print(f'  需要的 Stress ≈ {required_stress:.2f}')
print(f'  實際累積的 Stress = {Stress_final:.4f}')
print(f'  差距: 需要 {required_stress / Stress_final:.1f}x 的 Stress!')
print()

print('結論:')
print('  H12D3 只照 3 天，Stress 累積遠遠不足')
print('  即使大幅提高 stress_damage_coeff，也難以在 3 天內累積足夠 Stress')
print('  需要考慮其他機制，或者接受這個誤差')
