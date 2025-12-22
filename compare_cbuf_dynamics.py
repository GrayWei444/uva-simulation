#!/usr/bin/env python3
"""
對比 Sun 模型和 UVA 模型的碳緩衝池動態
找出為什麼 UVA 模型的 dC_buf/dt 異常
"""

import numpy as np
from lettuce_uva_carbon_complete_model import SunParams, sun_derivatives_final
from simulate_uva_model import UVAParams, uva_sun_derivatives
from model_config import ENV_BASE, TREATMENT_CONFIGS

# 測試點: Day 30, 12:00 (UVA 照射期間，日間)
t_test = (30 * 86400) + (12 * 3600)

# 相同的狀態
X_d = 0.1
C_buf = 0.01
LAI = 5.0

print('=' * 80)
print('對比 Sun 模型 vs UVA 模型的碳緩衝池動態')
print('=' * 80)
print()
print(f'測試時間: Day 30, 12:00')
print(f'狀態: X_d={X_d}, C_buf={C_buf}, LAI={LAI}')
print()

# Sun 模型
p_sun = SunParams()
p_sun.c_alpha = 0.57
env_sun = ENV_BASE.copy()
env_sun['I_day'] = 79  # 57 + 22

y_sun = [X_d, C_buf, LAI]
derivs_sun = sun_derivatives_final(t_test, y_sun, p_sun, env_sun)

print('Sun 模型輸出:')
print('-' * 80)
print(f'  dX_d/dt  = {derivs_sun[0]:.8e} kg/m²/s')
print(f'  dC_buf/dt = {derivs_sun[1]:.8e} kg C/m²/s')
print(f'  dLAI/dt  = {derivs_sun[2]:.8e} /s')
print()

# UVA 模型
p_uva = UVAParams()
env_uva = ENV_BASE.copy()
env_uva.update(TREATMENT_CONFIGS['L6D6'])

y_uva = [X_d, C_buf, LAI, 0.0, 0.0, 0.0]
derivs_uva = uva_sun_derivatives(t_test, y_uva, p_uva, env_uva)

print('UVA 模型輸出:')
print('-' * 80)
print(f'  dX_d/dt  = {derivs_uva[0]:.8e} kg/m²/s')
print(f'  dC_buf/dt = {derivs_uva[1]:.8e} kg C/m²/s')
print(f'  dLAI/dt  = {derivs_uva[2]:.8e} /s')
print(f'  dAnth/dt = {derivs_uva[3]:.8e} kg/m²/s')
print(f'  dStress/dt = {derivs_uva[4]:.8e} /s')
print()

# 計算差異
print('差異分析:')
print('-' * 80)
dx_diff = (derivs_uva[0] - derivs_sun[0]) / abs(derivs_sun[0]) * 100 if derivs_sun[0] != 0 else 0
dc_diff = (derivs_uva[1] - derivs_sun[1]) / abs(derivs_sun[1]) * 100 if derivs_sun[1] != 0 else 0

print(f'  dX_d/dt 差異:  {dx_diff:+.1f}%')
print(f'  dC_buf/dt 差異: {dc_diff:+.1f}%')
print()

# 分析 dC_buf/dt 的組成
print('UVA 模型的碳消耗分析:')
print('-' * 80)

# 花青素碳消耗
synthesis_rate = p_uva.base_anth_rate_light  # 日間
anth_carbon_consumption = synthesis_rate * p_uva.anth_carbon_cost
print(f'  花青素合成率: {synthesis_rate:.3e} kg/m²/s')
print(f'  花青素碳成本係數: {p_uva.anth_carbon_cost}')
print(f'  花青素碳消耗: {anth_carbon_consumption:.3e} kg C/m²/s')
print()

# 修復碳消耗
stress = 0.0  # L6D6 的 Stress 很小
repair_rate = p_uva.stress_repair_coeff * stress
repair_carbon_consumption = repair_rate * p_uva.repair_carbon_cost
print(f'  Stress: {stress}')
print(f'  修復率: {repair_rate:.3e} /s')
print(f'  修復碳消耗: {repair_carbon_consumption:.3e} kg C/m²/s')
print()

total_carbon_consumption = anth_carbon_consumption + repair_carbon_consumption
print(f'  總碳消耗: {total_carbon_consumption:.3e} kg C/m²/s')
print(f'  Sun 模型 dC_buf/dt: {derivs_sun[1]:.3e} kg C/m²/s')
print()

if anth_carbon_consumption > abs(derivs_sun[1]):
    print('⚠️  花青素碳消耗 > Sun 模型的碳生產！')
    print('    這就是為什麼 UVA 模型的 dC_buf/dt 是負的！')
else:
    print('✓ 碳消耗在合理範圍內')
