#!/usr/bin/env python3
"""
直接對比 sun_derivatives_final 和 uva_sun_derivatives
在完全相同的狀態和環境下的輸出
"""

import numpy as np
from lettuce_uva_carbon_complete_model import SunParams, sun_derivatives_final
from simulate_uva_model import UVAParams, uva_sun_derivatives
from model_config import ENV_BASE

# 測試時間點: Day 30, 12:00
t_test = (30 * 86400) + (12 * 3600)

# 相同的狀態
X_d = 0.1
C_buf = 0.01
LAI = 5.0

print('=' * 80)
print('直接對比 sun_derivatives_final vs uva_sun_derivatives (CK 無 UVA 條件)')
print('=' * 80)
print()
print(f'測試時間: Day 30, 12:00')
print(f'狀態: X_d={X_d}, C_buf={C_buf}, LAI={LAI}')
print()

# ==============================================================================
# Sun 原始模型
# ==============================================================================
print('Sun 原始模型:')
print('-' * 80)

p_sun = SunParams()
p_sun.c_alpha = 0.57

env_sun = ENV_BASE.copy()
env_sun['I_day'] = 57  # CK 的光強度

y_sun = [X_d, C_buf, LAI]
derivs_sun = sun_derivatives_final(t_test, y_sun, p_sun, env_sun)

print(f'  dX_d/dt  = {derivs_sun[0]:.10e} kg/m²/s')
print(f'  dC_buf/dt = {derivs_sun[1]:.10e} kg C/m²/s')
print(f'  dLAI/dt  = {derivs_sun[2]:.10e} /s')
print()

# ==============================================================================
# UVA 模型 (CK 設定: uva_on=False)
# ==============================================================================
print('UVA 模型 (CK 設定, uva_on=False):')
print('-' * 80)

p_uva = UVAParams()
p_uva.c_alpha = 0.57

env_uva = ENV_BASE.copy()
env_uva['I_day'] = 57
env_uva['uva_on'] = False
env_uva['uva_intensity'] = 0

y_uva = [X_d, C_buf, LAI, 0.0, 0.0, 0.0]  # [X_d, C_buf, LAI, Anth, Stress, E_stress]
derivs_uva = uva_sun_derivatives(t_test, y_uva, p_uva, env_uva)

print(f'  dX_d/dt  = {derivs_uva[0]:.10e} kg/m²/s')
print(f'  dC_buf/dt = {derivs_uva[1]:.10e} kg C/m²/s')
print(f'  dLAI/dt  = {derivs_uva[2]:.10e} /s')
print(f'  dAnth/dt = {derivs_uva[3]:.10e} kg/m²/s')
print(f'  dStress/dt = {derivs_uva[4]:.10e} /s')
print()

# ==============================================================================
# 差異分析
# ==============================================================================
print('差異分析:')
print('-' * 80)

dx_diff = derivs_uva[0] - derivs_sun[0]
dc_diff = derivs_uva[1] - derivs_sun[1]
dl_diff = derivs_uva[2] - derivs_sun[2]

dx_pct = (dx_diff / abs(derivs_sun[0]) * 100) if derivs_sun[0] != 0 else 0
dc_pct = (dc_diff / abs(derivs_sun[1]) * 100) if derivs_sun[1] != 0 else 0
dl_pct = (dl_diff / abs(derivs_sun[2]) * 100) if derivs_sun[2] != 0 else 0

print(f'  dX_d/dt 差異:  {dx_diff:+.3e} ({dx_pct:+.1f}%)')
print(f'  dC_buf/dt 差異: {dc_diff:+.3e} ({dc_pct:+.1f}%)')
print(f'  dLAI/dt 差異:  {dl_diff:+.3e} ({dl_pct:+.1f}%)')
print()

if abs(dx_pct) > 1 or abs(dc_pct) > 1 or abs(dl_pct) > 1:
    print('✗✗✗ 問題發現！')
    print('    UVA 模型在 CK 設定 (無 UVA) 下與 Sun 模型不一致！')
    print('    這就是為什麼 CK 預測錯誤的根本原因！')
    print()
    print('可能的原因:')
    print('  1. uva_sun_derivatives 對 Sun 返回值做了不正確的修改')
    print('  2. Stress 相關的計算即使 I_UVA=0 也影響了結果')
    print('  3. 花青素合成消耗碳，即使沒有 UVA 也在合成')
else:
    print('✓ UVA 模型與 Sun 模型一致')

print()
print('=' * 80)
