#!/usr/bin/env python3
"""
深度追蹤 dC_buf/dt 的計算過程
逐步輸出每個環節的值，找出為什麼 UVA 模型的 dC_buf/dt 變成負值
"""

import numpy as np
from lettuce_uva_carbon_complete_model import SunParams, sun_derivatives_final
from simulate_uva_model import UVAParams
from model_config import ENV_BASE, TREATMENT_CONFIGS

# 測試點: Day 30, 12:00 (UVA 照射期間，日間)
t_test = (30 * 86400) + (12 * 3600)

# 相同的狀態
X_d = 0.1
C_buf = 0.01
LAI = 5.0
Anth = 0.0
Stress = 0.0
E_stress = 0.0

print('=' * 80)
print('深度追蹤 dC_buf/dt 計算過程')
print('=' * 80)
print()
print(f'測試時間: Day 30, 12:00')
print(f'狀態: X_d={X_d}, C_buf={C_buf}, LAI={LAI}, Stress={Stress}')
print()

# ==============================================================================
# 第一部分: Sun 模型的 dC_buf/dt
# ==============================================================================
print('=' * 80)
print('第一部分: Sun 原始模型')
print('=' * 80)
print()

p_sun = SunParams()
p_sun.c_alpha = 0.57

env_sun = ENV_BASE.copy()
env_sun['I_day'] = 79  # 57 + 22

y_sun = [X_d, C_buf, LAI]
derivs_sun = sun_derivatives_final(t_test, y_sun, p_sun, env_sun)

print(f'Sun 模型輸出:')
print(f'  dX_d/dt  = {derivs_sun[0]:.10e} kg/m²/s')
print(f'  dC_buf/dt = {derivs_sun[1]:.10e} kg C/m²/s  ← 這個值是正的')
print(f'  dLAI/dt  = {derivs_sun[2]:.10e} /s')
print()

# ==============================================================================
# 第二部分: UVA 模型的 dC_buf/dt (手動逐步計算)
# ==============================================================================
print('=' * 80)
print('第二部分: UVA 模型 (逐步追蹤)')
print('=' * 80)
print()

p_uva = UVAParams()
env_uva = ENV_BASE.copy()
env_uva.update(TREATMENT_CONFIGS['L6D6'])

# --- 步驟1: 時間相關參數 ---
days_from_transplant = t_test / 86400.0
hour = (t_test / 3600.0) % 24
is_day = (hour >= 6) and (hour < 18)
is_night = not is_day

print('步驟1: 時間參數')
print('-' * 80)
print(f'  days_from_transplant = {days_from_transplant:.2f}')
print(f'  hour = {hour:.2f}')
print(f'  is_day = {is_day}')
print()

# --- 步驟2: UVA 照射判斷 ---
uva_on = env_uva.get('uva_on', False)
uva_start_day = env_uva.get('uva_start_day', 29)
uva_end_day = env_uva.get('uva_end_day', 35)
uva_hour_on = env_uva.get('uva_hour_on', 10)
uva_hour_off = env_uva.get('uva_hour_off', 16)
uva_intensity = env_uva.get('uva_intensity', 22.0)

I_UVA = 0.0
in_uva_window = False

if uva_on and uva_start_day <= days_from_transplant < uva_end_day:
    if uva_hour_on <= uva_hour_off:
        in_uva_window = uva_hour_on <= hour < uva_hour_off
    else:
        in_uva_window = hour >= uva_hour_on or hour < uva_hour_off

    if in_uva_window:
        I_UVA = uva_intensity

print('步驟2: UVA 照射判斷')
print('-' * 80)
print(f'  in_uva_window = {in_uva_window}')
print(f'  I_UVA = {I_UVA} W/m²')
print()

# --- 步驟3: 有效光強度 ---
I_base = env_uva['I_day']
I_gain_par = I_UVA * p_uva.par_conversion_factor
I_effective = I_base + I_gain_par

print('步驟3: 有效光強度')
print('-' * 80)
print(f'  I_base = {I_base} W/m²')
print(f'  I_UVA = {I_UVA} W/m²')
print(f'  par_conversion_factor = {p_uva.par_conversion_factor}')
print(f'  I_gain_par = {I_gain_par} W/m²')
print(f'  I_effective = {I_effective} W/m²')
print()

# --- 步驟4: 調用 Sun 模型 ---
print('步驟4: 調用 Sun 模型')
print('-' * 80)

# 建立 Sun 模型的環境
env_for_sun = env_uva.copy()
env_for_sun['I_day'] = I_effective

# UVAParams 繼承自 BaseSunParams，所以直接可以傳給 sun_derivatives_final
y_for_sun = [X_d, C_buf, LAI]
derivs_from_sun = sun_derivatives_final(t_test, y_for_sun, p_uva, env_for_sun)

dXd_dt_base = derivs_from_sun[0]
dCbuf_dt_from_sun = derivs_from_sun[1]
dLAI_dt_base = derivs_from_sun[2]

print(f'  Sun 模型返回的值:')
print(f'    dXd_dt_base = {dXd_dt_base:.10e}')
print(f'    dCbuf_dt_from_sun = {dCbuf_dt_from_sun:.10e}  ← 這裡應該是正值')
print(f'    dLAI_dt_base = {dLAI_dt_base:.10e}')
print()

dCbuf_dt = dCbuf_dt_from_sun
print(f'  初始化: dCbuf_dt = {dCbuf_dt:.10e}')
print()

# --- 步驟5: Stress 動態 (簡化，只計算修復率) ---
print('步驟5: Stress 動態與修復')
print('-' * 80)

C_buf_positive = max(C_buf, 0)
repair_capacity = p_uva.base_repair_capacity + \
                  p_uva.carbon_repair_bonus * C_buf_positive / (p_uva.K_carbon + C_buf_positive + 1e-9)

repair_rate = p_uva.stress_repair_coeff * Stress * repair_capacity

print(f'  Stress = {Stress}')
print(f'  C_buf_positive = {C_buf_positive}')
print(f'  repair_capacity = {repair_capacity}')
print(f'  repair_rate = {repair_rate:.10e} /s')
print()

# --- 步驟6: 花青素合成 ---
print('步驟6: 花青素合成')
print('-' * 80)

day_weight = 1.0 if is_day else 0.0
base_synthesis = day_weight * p_uva.base_anth_rate_light + (1 - day_weight) * p_uva.base_anth_rate_dark

Stress_power_n = Stress ** p_uva.n_stress_anth
K_power_n = p_uva.K_stress_anth ** p_uva.n_stress_anth
uva_induced = p_uva.V_max_anth * Stress_power_n / (K_power_n + Stress_power_n + 1e-12)

synthesis_rate = base_synthesis + uva_induced

print(f'  is_day = {is_day}')
print(f'  base_synthesis = {base_synthesis:.10e} kg/m²/s')
print(f'  uva_induced = {uva_induced:.10e} kg/m²/s')
print(f'  synthesis_rate = {synthesis_rate:.10e} kg/m²/s')
print()

# --- 步驟7: 碳消耗計算 ---
print('步驟7: 碳消耗計算')
print('-' * 80)

repair_carbon_consumption = repair_rate * p_uva.repair_carbon_cost
anth_carbon_consumption = synthesis_rate * p_uva.anth_carbon_cost

print(f'  repair_carbon_cost = {p_uva.repair_carbon_cost}')
print(f'  anth_carbon_cost = {p_uva.anth_carbon_cost}')
print(f'  repair_carbon_consumption = {repair_carbon_consumption:.10e} kg C/m²/s')
print(f'  anth_carbon_consumption = {anth_carbon_consumption:.10e} kg C/m²/s')
print()

total_carbon_consumption = repair_carbon_consumption + anth_carbon_consumption
print(f'  總碳消耗 = {total_carbon_consumption:.10e} kg C/m²/s')
print()

# --- 步驟8: 更新 dCbuf_dt ---
print('步驟8: 更新 dCbuf_dt')
print('-' * 80)

print(f'  更新前: dCbuf_dt = {dCbuf_dt:.10e}')
dCbuf_dt_new = dCbuf_dt - repair_carbon_consumption - anth_carbon_consumption
print(f'  更新後: dCbuf_dt = {dCbuf_dt_new:.10e}')
print()

# --- 步驟9: Stress 對生長的抑制 ---
print('步驟9: Stress 對生長的抑制 (不影響 dCbuf_dt)')
print('-' * 80)

stress_inhibition = Stress / (p_uva.K_stress + Stress + 1e-9)
xd_reduction = p_uva.stress_photosynthesis_inhibition * stress_inhibition
lai_reduction = p_uva.stress_lai_inhibition * stress_inhibition

print(f'  stress_inhibition = {stress_inhibition:.6f}')
print(f'  xd_reduction = {xd_reduction:.6f}')
print(f'  lai_reduction = {lai_reduction:.6f}')
print()

if dXd_dt_base > 0:
    dXd_dt = dXd_dt_base * (1.0 - xd_reduction)
else:
    dXd_dt = dXd_dt_base

if dLAI_dt_base > 0:
    dLAI_dt = dLAI_dt_base * (1.0 - lai_reduction)
else:
    dLAI_dt = dLAI_dt_base

print(f'  最終 dXd_dt = {dXd_dt:.10e}')
print(f'  最終 dLAI_dt = {dLAI_dt:.10e}')
print()

# ==============================================================================
# 第三部分: 對比分析
# ==============================================================================
print('=' * 80)
print('第三部分: 對比分析')
print('=' * 80)
print()

print('dCbuf_dt 變化追蹤:')
print('-' * 80)
print(f'  Sun 模型返回值:       {dCbuf_dt_from_sun:+.10e}')
print(f'  減去修復碳消耗後:     {dCbuf_dt_from_sun - repair_carbon_consumption:+.10e}')
print(f'  減去花青素碳消耗後:   {dCbuf_dt_new:+.10e}')
print()

# 與 Sun 原始模型對比
print('與 Sun 原始模型 (I=79 W/m²) 對比:')
print('-' * 80)
print(f'  Sun 模型:  dCbuf_dt = {derivs_sun[1]:+.10e}')
print(f'  UVA 模型:  dCbuf_dt = {dCbuf_dt_new:+.10e}')
print()

diff_pct = (dCbuf_dt_new - derivs_sun[1]) / abs(derivs_sun[1]) * 100 if derivs_sun[1] != 0 else 0
print(f'  差異: {diff_pct:+.1f}%')
print()

# 診斷
print('診斷:')
print('-' * 80)

if abs(dCbuf_dt_from_sun - derivs_sun[1]) > 1e-12:
    print('  ✗ 問題1: Sun 模型調用有問題！')
    print(f'    UVA 中調用 Sun 返回: {dCbuf_dt_from_sun:.10e}')
    print(f'    直接調用 Sun 返回:   {derivs_sun[1]:.10e}')
    print(f'    差異: {abs(dCbuf_dt_from_sun - derivs_sun[1]):.10e}')
else:
    print('  ✓ Sun 模型調用正確')

print()

if total_carbon_consumption > abs(dCbuf_dt_from_sun) and dCbuf_dt_from_sun > 0:
    print('  ✗ 問題2: 碳消耗過大！')
    print(f'    Sun 模型碳生產: {dCbuf_dt_from_sun:.10e}')
    print(f'    總碳消耗:       {total_carbon_consumption:.10e}')
    print(f'    消耗/生產比:    {total_carbon_consumption / dCbuf_dt_from_sun:.1f}x')
else:
    print('  ✓ 碳消耗在合理範圍內')

print()

if abs(diff_pct) > 10:
    print('  ✗ 問題3: UVA 模型與 Sun 模型差異過大！')
    print(f'    可能原因:')
    if repair_carbon_consumption > 1e-12:
        print(f'      - 修復碳消耗: {repair_carbon_consumption:.10e}')
    if anth_carbon_consumption > 1e-12:
        print(f'      - 花青素碳消耗: {anth_carbon_consumption:.10e}')
else:
    print('  ✓ UVA 模型與 Sun 模型一致')

print()
print('=' * 80)
