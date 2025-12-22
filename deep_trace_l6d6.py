#!/usr/bin/env python3
"""
深度追蹤 L6D6 模擬過程
逐步檢查每個計算環節，找出為什麼預測這麼低
"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives
from lettuce_uva_carbon_complete_model import sun_derivatives_final, SunParams
from model_config import ENV_BASE, TREATMENT_CONFIGS, ODE_SETTINGS, SIMULATION

print('=' * 80)
print('L6D6 深度追蹤診斷')
print('=' * 80)
print()

# 設定參數
p = UVAParams()
env = ENV_BASE.copy()
env.update(TREATMENT_CONFIGS['L6D6'])

# 初始條件
fw_init_g = SIMULATION['initial_fw_g']
dw_init_g = fw_init_g * p.dw_fw_ratio_base
X_d_init = dw_init_g / 1000.0 * env['plant_density']
C_buf_init = X_d_init * 0.1
LAI_init = (dw_init_g / 0.01) * 0.04
fw_total_init = fw_init_g * env['plant_density'] / 1000
Anth_init = 5.0 * fw_total_init / 1e6

y0 = [X_d_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]

print('初始狀態:')
print(f'  X_d = {X_d_init:.8f} kg/m²')
print(f'  C_buf = {C_buf_init:.8f} kg C/m²')
print(f'  LAI = {LAI_init:.4f}')
print()

# 模擬時間
simulation_days = SIMULATION['days']
transplant_day = SIMULATION['transplant_offset']
t_start = transplant_day * 86400
t_end = (transplant_day + simulation_days) * 86400

# 關鍵測試時間點: Day 30, 12:00 (UVA 照射期間)
t_test = (30 * 86400) + (12 * 3600)

print('測試時間點: Day 30, 12:00 (UVA 照射期間)')
print()

# 先進行完整模擬，取得該時間點的狀態
sol = solve_ivp(
    fun=lambda t, y: uva_sun_derivatives(t, y, p, env),
    t_span=(t_start, t_end),
    y0=y0,
    method=ODE_SETTINGS['method'],
    max_step=ODE_SETTINGS['max_step'],
    dense_output=True
)

y_test = sol.sol(t_test)
X_d, C_buf, LAI, Anth, Stress, E_stress = y_test

print('Day 30 12:00 的實際狀態:')
print(f'  X_d = {X_d:.8f} kg/m²')
print(f'  C_buf = {C_buf:.10f} kg C/m²')
print(f'  LAI = {LAI:.6f}')
print(f'  Anth = {Anth:.8e} kg/m²')
print(f'  Stress = {Stress:.8f}')
print(f'  E_stress = {E_stress:.4f}')
print()

# ============================================================================
# 手動計算 uva_sun_derivatives 在這個時間點的輸出
# ============================================================================
print('=' * 80)
print('手動追蹤 uva_sun_derivatives 計算過程')
print('=' * 80)
print()

# 時間參數
days_from_transplant = t_test / 86400.0
hour = (t_test / 3600.0) % 24
is_day = (hour >= 6) and (hour < 18)

print(f'步驟 1: 時間參數')
print(f'  days_from_transplant = {days_from_transplant:.2f}')
print(f'  hour = {hour:.2f}')
print(f'  is_day = {is_day}')
print()

# UVA 照射判斷
uva_on = env.get('uva_on', False)
uva_start_day = env.get('uva_start_day', 29)
uva_end_day = env.get('uva_end_day', 35)
uva_hour_on = env.get('uva_hour_on', 10)
uva_hour_off = env.get('uva_hour_off', 16)
uva_intensity = env.get('uva_intensity', 22.0)

I_UVA = 0.0
in_uva_window = False

if uva_on and uva_start_day <= days_from_transplant < uva_end_day:
    if uva_hour_on <= uva_hour_off:
        in_uva_window = uva_hour_on <= hour < uva_hour_off
    else:
        in_uva_window = hour >= uva_hour_on or hour < uva_hour_off

    if in_uva_window:
        I_UVA = uva_intensity

print(f'步驟 2: UVA 照射判斷')
print(f'  in_uva_window = {in_uva_window}')
print(f'  I_UVA = {I_UVA} W/m²')
print()

# 有效光強度
I_base = env['I_day']
I_gain_par = I_UVA * p.par_conversion_factor
I_effective = I_base + I_gain_par

print(f'步驟 3: 有效光強度計算')
print(f'  I_base (LED) = {I_base} W/m²')
print(f'  I_UVA = {I_UVA} W/m²')
print(f'  par_conversion_factor = {p.par_conversion_factor}')
print(f'  I_gain_par = {I_gain_par} W/m²')
print(f'  I_effective = {I_effective} W/m²  ← 這應該等於 79 W/m²')
print()

# 調用 Sun 模型
print(f'步驟 4: 調用 Sun 原始模型')
print(f'  輸入狀態: X_d={X_d:.6f}, C_buf={C_buf:.8f}, LAI={LAI:.4f}')
print(f'  有效光強: I_effective={I_effective} W/m²')

env_for_sun = env.copy()
env_for_sun['I_day'] = I_effective

y_for_sun = [X_d, C_buf, LAI]
derivs_sun = sun_derivatives_final(t_test, y_for_sun, p, env_for_sun)

dXd_dt_base = derivs_sun[0]
dCbuf_dt_base = derivs_sun[1]
dLAI_dt_base = derivs_sun[2]

print(f'  Sun 模型輸出:')
print(f'    dX_d/dt  = {dXd_dt_base:.10e} kg/m²/s')
print(f'    dC_buf/dt = {dCbuf_dt_base:.10e} kg C/m²/s')
print(f'    dLAI/dt  = {dLAI_dt_base:.10e} /s')
print()

# 對比純 CK (無 UVA) 的 Sun 模型
print(f'步驟 4b: 對比 CK (I=57 W/m²) 的 Sun 模型')
env_ck = env.copy()
env_ck['I_day'] = 57

derivs_ck = sun_derivatives_final(t_test, y_for_sun, p, env_ck)
print(f'  CK Sun 模型輸出:')
print(f'    dX_d/dt  = {derivs_ck[0]:.10e} kg/m²/s')
print(f'    dC_buf/dt = {derivs_ck[1]:.10e} kg C/m²/s')
print()
print(f'  對比:')
print(f'    L6D6 vs CK dX_d/dt 比例: {dXd_dt_base / derivs_ck[0]:.4f}x')
print(f'    L6D6 vs CK dC_buf/dt 比例: {dCbuf_dt_base / derivs_ck[1]:.4f}x')
print()

# Stress 相關計算
print(f'步驟 5: Stress 損傷與修復')
print(f'  當前 Stress = {Stress:.8f}')
print(f'  I_UVA = {I_UVA} W/m²')
print()

# 損傷率計算
# LAI 脆弱性
LAI_ratio = LAI / p.LAI_ref_vuln
vulnerability = min(LAI_ratio ** p.n_vuln, p.cap_vuln)

# 日累積能量
uva_start_hour = 10
if hour >= uva_start_hour and in_uva_window:
    hours_elapsed = hour - uva_start_hour
else:
    hours_elapsed = 0

E_elapsed = I_UVA * hours_elapsed * 3.6  # kJ/m²
normalized_E = (E_elapsed - p.E_50) / p.E_scale

# Softplus
def softplus(x, sharpness):
    return np.log(1.0 + np.exp(sharpness * x)) / sharpness

excess_normalized = softplus(normalized_E, p.sharpness_intraday)
intraday_factor = 1.0 + p.k_intraday * (excess_normalized ** p.m_intraday)

# Stress 非線性
nonlinear_factor = 1.0 + p.stress_nonlinear_coeff * Stress / (p.K_nonlinear + Stress + 1e-9)

# 節律損傷
is_night_uva = (not is_day) and in_uva_window
circadian_penalty = p.circadian_disruption_factor if is_night_uva else 1.0

# 損傷率
damage_rate = p.stress_damage_coeff * I_UVA * vulnerability * intraday_factor * nonlinear_factor * circadian_penalty

print(f'  損傷率計算:')
print(f'    vulnerability = {vulnerability:.6f}')
print(f'    hours_elapsed = {hours_elapsed:.2f} h')
print(f'    E_elapsed = {E_elapsed:.1f} kJ/m²')
print(f'    normalized_E = {normalized_E:.4f}')
print(f'    intraday_factor = {intraday_factor:.6f}')
print(f'    nonlinear_factor = {nonlinear_factor:.6f}')
print(f'    circadian_penalty = {circadian_penalty:.2f}')
print(f'    damage_rate = {damage_rate:.10e} /s')
print()

# 修復率
C_buf_positive = max(C_buf, 0)
repair_capacity = p.base_repair_capacity + p.carbon_repair_bonus * C_buf_positive / (p.K_carbon + C_buf_positive + 1e-9)
repair_rate = p.stress_repair_coeff * Stress * repair_capacity

print(f'  修復率計算:')
print(f'    C_buf_positive = {C_buf_positive:.10f}')
print(f'    repair_capacity = {repair_capacity:.6f}')
print(f'    repair_rate = {repair_rate:.10e} /s')
print()

dStress_dt = damage_rate - repair_rate
print(f'  dStress/dt = {damage_rate:.10e} - {repair_rate:.10e} = {dStress_dt:.10e} /s')
print()

# Stress 抑制
stress_inhibition = Stress / (p.K_stress + Stress + 1e-9)
xd_reduction = p.stress_photosynthesis_inhibition * stress_inhibition
lai_reduction = p.stress_lai_inhibition * stress_inhibition

print(f'步驟 6: Stress 對生長的抑制')
print(f'  stress_inhibition = {stress_inhibition:.8f}')
print(f'  xd_reduction = {xd_reduction:.8f} ({xd_reduction*100:.2f}%)')
print(f'  lai_reduction = {lai_reduction:.8f} ({lai_reduction*100:.2f}%)')
print()

# 應用抑制
if dXd_dt_base > 0:
    dXd_dt = dXd_dt_base * (1.0 - xd_reduction)
else:
    dXd_dt = dXd_dt_base

if dLAI_dt_base > 0:
    dLAI_dt = dLAI_dt_base * (1.0 - lai_reduction)
else:
    dLAI_dt = dLAI_dt_base

print(f'  最終導數:')
print(f'    dX_d/dt  = {dXd_dt_base:.10e} × (1 - {xd_reduction:.4f}) = {dXd_dt:.10e}')
print(f'    dLAI/dt  = {dLAI_dt_base:.10e} × (1 - {lai_reduction:.4f}) = {dLAI_dt:.10e}')
print()

# ============================================================================
# 診斷
# ============================================================================
print('=' * 80)
print('診斷')
print('=' * 80)
print()

# 檢查 1: 有效光強度
if abs(I_effective - 79) < 0.1:
    print('✓ 有效光強度正確: I_effective = 79 W/m²')
else:
    print(f'✗ 有效光強度異常: I_effective = {I_effective} W/m² (應為 79 W/m²)')

# 檢查 2: Sun 模型光合促進效應
sun_ratio = dXd_dt_base / derivs_ck[0]
if sun_ratio > 1.0:
    print(f'✓ Sun 模型顯示促進效應: dX_d/dt 提升 {(sun_ratio-1)*100:.1f}%')
else:
    print(f'✗ Sun 模型無促進效應: dX_d/dt 比例 {sun_ratio:.4f}x')

# 檢查 3: Stress 抑制程度
if xd_reduction > 0.1:
    print(f'⚠️  Stress 抑制顯著: {xd_reduction*100:.1f}% 抑制')
    print(f'    抑制後 dX_d/dt = {dXd_dt:.10e}')
    print(f'    抑制前 dX_d/dt = {dXd_dt_base:.10e}')
    print(f'    損失比例 = {xd_reduction*100:.1f}%')
else:
    print(f'✓ Stress 抑制較小: {xd_reduction*100:.2f}%')

# 檢查 4: 淨效應
final_vs_ck_ratio = dXd_dt / derivs_ck[0]
print()
print(f'淨效應分析:')
print(f'  CK (57 W/m²) dX_d/dt = {derivs_ck[0]:.10e}')
print(f'  L6D6 最終 dX_d/dt    = {dXd_dt:.10e}')
print(f'  比例 = {final_vs_ck_ratio:.4f}x ({(final_vs_ck_ratio-1)*100:+.1f}%)')
print()

if final_vs_ck_ratio < 1.0:
    print(f'✗✗✗ 問題確認！')
    print(f'    L6D6 的最終生長速率比 CK **更低**')
    print(f'    即使有額外 22 W/m² 輻射促進光合')
    print(f'    Stress 抑制效應 ({xd_reduction*100:.1f}%) 完全抵消了促進效應')
    print()
    print(f'    輻射促進: +{(sun_ratio-1)*100:.1f}%')
    print(f'    Stress 抑制: -{xd_reduction*100:.1f}%')
    print(f'    淨效應: {(final_vs_ck_ratio-1)*100:+.1f}%')
else:
    print(f'✓ L6D6 生長速率高於 CK (+{(final_vs_ck_ratio-1)*100:.1f}%)')

print()
print('=' * 80)
