#!/usr/bin/env python3
"""
追蹤 Day 15 Hour 12 的 damage_rate 計算
為什麼 I_UVA = 0 但 dStress/dt > 0？
"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives
from model_config import ENV_BASE, TREATMENT_CONFIGS, SIMULATION

# Softplus helper
def softplus(x, sharpness):
    return np.log(1.0 + np.exp(sharpness * x)) / sharpness

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

# 模擬到 Day 15 Hour 12
t_start = 14 * 86400
t_target = (15 * 86400) + (12 * 3600)

sol = solve_ivp(
    fun=lambda t, y: uva_sun_derivatives(t, y, p, env),
    t_span=(t_start, t_target),
    y0=y0,
    method='RK45',
    max_step=60,
    dense_output=True
)

# Day 15 Hour 12 的狀態
y = sol.sol(t_target)
X_d, C_buf, LAI, Anth, Stress, E_stress = y

print('=' * 80)
print('Day 15 Hour 12 損傷率計算詳細追蹤')
print('=' * 80)
print()

print('狀態:')
print(f'  X_d = {X_d:.8f}')
print(f'  C_buf = {C_buf:.10f}')
print(f'  LAI = {LAI:.4f}')
print(f'  Anth = {Anth:.8e}')
print(f'  Stress = {Stress:.6f}')
print(f'  E_stress = {E_stress:.8f}')
print()

# 手動重現損傷率計算 (simulate_uva_model.py lines 272-346)
t = t_target
hour = (t / 3600.0) % 24
day_from_sowing = t / 86400.0 + p.transplant_day

print(f'時間:')
print(f'  hour = {hour}')
print(f'  day_from_sowing = {day_from_sowing}')
print()

# UVA 判斷
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

print(f'UVA 判斷:')
print(f'  uva_on = {uva_on}')
print(f'  uva_start_day = {uva_start_day}')
print(f'  day_from_sowing = {day_from_sowing}')
print(f'  uva_start_day <= day_from_sowing < uva_end_day = {uva_start_day <= day_from_sowing < uva_end_day}')
print(f'  → in_uva_window = {in_uva_window}')
print(f'  → I_UVA = {I_UVA} W/m²')
print()

# LAI vulnerability
vulnerability = p.cap_vuln * (p.LAI_ref_vuln / LAI) ** p.n_vuln / (p.cap_vuln + (p.LAI_ref_vuln / LAI) ** p.n_vuln)

print(f'LAI Vulnerability:')
print(f'  LAI = {LAI:.4f}')
print(f'  LAI_ref_vuln = {p.LAI_ref_vuln}')
print(f'  vulnerability = {vulnerability:.6f}')
print()

# Intraday factor
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

I_UVA_config = env.get('I_UVA', 11.0)  # ← BUG in simulate_uva_model.py
E_elapsed = I_UVA_config * hours_elapsed * 3.6

normalized_E = (E_elapsed - p.E_50) / p.E_scale
excess_normalized = softplus(normalized_E, p.sharpness_intraday)
intraday_factor = 1.0 + p.k_intraday * (excess_normalized ** p.m_intraday)

print(f'Intraday Factor:')
print(f'  in_uva_window = {in_uva_window}')
print(f'  hours_elapsed = {hours_elapsed}')
print(f'  I_UVA_config = {I_UVA_config} ← env.get("I_UVA", 11.0)')
print(f'  E_elapsed = {E_elapsed:.1f} kJ/m²')
print(f'  normalized_E = {normalized_E:.6f}')
print(f'  intraday_factor = {intraday_factor:.6f}')
print()

# Nonlinear factor
nonlinear_factor = 1.0 + p.stress_nonlinear_coeff * Stress / (p.K_nonlinear + Stress + 1e-9)

print(f'Nonlinear Factor:')
print(f'  Stress = {Stress:.6f}')
print(f'  nonlinear_factor = {nonlinear_factor:.6f}')
print()

# Circadian penalty
is_day = (6 <= hour < 18)  # 簡化判斷
is_night_uva = (I_UVA > 0) and (not is_day)
circadian_penalty = p.circadian_disruption_factor if is_night_uva else 1.0

print(f'Circadian:')
print(f'  is_night_uva = {is_night_uva}')
print(f'  circadian_penalty = {circadian_penalty}')
print()

# Damage rate (simulate_uva_model.py line 332-333)
damage_rate = p.stress_damage_coeff * I_UVA * vulnerability * \
              intraday_factor * nonlinear_factor * circadian_penalty

print(f'=' * 80)
print(f'損傷率計算:')
print(f'=' * 80)
print(f'  damage_rate = stress_damage_coeff × I_UVA × vulnerability × intraday_factor × nonlinear_factor × circadian_penalty')
print(f'              = {p.stress_damage_coeff} × {I_UVA} × {vulnerability:.6f} × {intraday_factor:.6f} × {nonlinear_factor:.6f} × {circadian_penalty}')
print(f'              = {damage_rate:.10e}')
print()

if I_UVA == 0 and damage_rate > 1e-10:
    print(f'✗✗✗ BUG 發現！')
    print(f'    I_UVA = 0 但 damage_rate = {damage_rate:.10e} > 0')
    print(f'    這是因為 I_UVA = 0 時，damage_rate 應該強制為 0')
elif I_UVA == 0:
    print(f'✓ I_UVA = 0, damage_rate = 0 (正確)')
else:
    print(f'I_UVA > 0, damage_rate = {damage_rate:.10e}')

print()
print('=' * 80)
