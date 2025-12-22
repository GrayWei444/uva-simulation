#!/usr/bin/env python3
"""
分析 L6D6 的 Stress 累積問題
為什麼 L6D6 (低劑量短期) 的 Stress 會比預期高？
"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives
from model_config import ENV_BASE, TREATMENT_CONFIGS, ODE_SETTINGS, SIMULATION

def simulate_and_track_stress(treatment_name):
    """模擬並追蹤 Stress 動態"""

    p = UVAParams()
    env = ENV_BASE.copy()
    env.update(TREATMENT_CONFIGS[treatment_name])

    plant_density = env['plant_density']

    # 初始條件
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    X_d_init = dw_init_g / 1000.0 * plant_density
    C_buf_init = X_d_init * 0.1
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * plant_density / 1000
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
print('L6D6 Stress 累積分析')
print('=' * 80)
print()

# L6D6 配置
l6d6_config = TREATMENT_CONFIGS['L6D6']
print('L6D6 配置:')
print(f'  UVA 強度: {l6d6_config["uva_intensity"]} W/m²')
print(f'  照射期間: Day {l6d6_config["uva_start_day"]} - {l6d6_config["uva_end_day"]} (共 {l6d6_config["uva_end_day"]-l6d6_config["uva_start_day"]} 天)')
print(f'  照射時間: {l6d6_config["uva_hour_on"]}:00 - {l6d6_config["uva_hour_off"]}:00 (日間 {l6d6_config["uva_hour_off"]-l6d6_config["uva_hour_on"]} 小時)')
print()

# 模擬 L6D6
sol, p, env = simulate_and_track_stress('L6D6')

# 關鍵時間點
key_times = {
    'Day 14 (開始)': 14 * 86400,
    'Day 28 (UVA 前)': 28 * 86400,
    'Day 29 12:00 (UVA 第1天中午)': (29 * 86400) + (12 * 3600),
    'Day 30 12:00 (UVA 第2天中午)': (30 * 86400) + (12 * 3600),
    'Day 32 12:00 (UVA 第4天中午)': (32 * 86400) + (12 * 3600),
    'Day 35 (結束)': 35 * 86400,
}

print('L6D6 關鍵時間點的狀態:')
print('-' * 80)
print(f'{"時間點":<30s} {"X_d":>12s} {"C_buf":>12s} {"LAI":>8s} {"Stress":>12s} {"E_stress":>12s}')
print('-' * 80)

for time_name, t in key_times.items():
    y = sol.sol(t)
    X_d, C_buf, LAI, Anth, Stress, E_stress = y
    print(f'{time_name:<30s} {X_d:12.6f} {C_buf:12.8f} {LAI:8.4f} {Stress:12.8f} {E_stress:12.4f}')

print()

# 計算 UVA 期間的平均 Stress
uva_start_t = 29 * 86400
uva_end_t = 35 * 86400
t_samples = np.linspace(uva_start_t, uva_end_t, 100)
stress_samples = [sol.sol(t)[4] for t in t_samples]
avg_stress = np.mean(stress_samples)
max_stress = np.max(stress_samples)

print('UVA 照射期間 Stress 統計:')
print('-' * 80)
print(f'  平均 Stress: {avg_stress:.6f}')
print(f'  最大 Stress: {max_stress:.6f}')
print(f'  最終 Stress: {sol.sol(uva_end_t)[4]:.6f}')
print()

# 分析日累積能量
print('日累積能量分析 (Day 30):')
print('-' * 80)

# Day 30 的 UVA 照射
uva_start_hour = l6d6_config['uva_hour_on']   # 10:00
uva_end_hour = l6d6_config['uva_hour_off']    # 16:00
uva_hours = uva_end_hour - uva_start_hour     # 6 小時
I_UVA = l6d6_config['uva_intensity']          # 22 W/m²

# 累積能量 [kJ/m²] = I [W/m²] × hours [h] × 3.6 [kJ/Wh]
E_elapsed = I_UVA * uva_hours * 3.6

print(f'  UVA 強度: {I_UVA} W/m²')
print(f'  照射時數: {uva_hours} 小時/天')
print(f'  日累積能量: {E_elapsed:.1f} kJ/m²')
print()

# 對比參數
print('日內非線性參數:')
print('-' * 80)
print(f'  E_50 (半飽和能量): {p.E_50} kJ/m²')
print(f'  E_scale (能量尺度): {p.E_scale} kJ/m²')
print(f'  k_intraday (放大係數): {p.k_intraday}')
print(f'  m_intraday (指數): {p.m_intraday}')
print()

# 計算正規化能量
normalized_E = (E_elapsed - p.E_50) / p.E_scale
print(f'正規化能量: (E - E_50) / E_scale = ({E_elapsed} - {p.E_50}) / {p.E_scale} = {normalized_E:.3f}')
print()

# 計算 softplus
def softplus(x, sharpness):
    return np.log(1.0 + np.exp(sharpness * x)) / sharpness

sp = softplus(normalized_E, p.sharpness_intraday)
intraday_factor = 1.0 + p.k_intraday * (sp ** p.m_intraday)

print(f'Softplus 輸出: {sp:.4f}')
print(f'日內非線性因子: 1 + {p.k_intraday} × {sp:.4f}^{p.m_intraday} = {intraday_factor:.4f}')
print()

# 診斷
print('=' * 80)
print('診斷:')
print('=' * 80)

if E_elapsed < p.E_50:
    print(f'✓ 日累積能量 ({E_elapsed:.1f} kJ/m²) < E_50 ({p.E_50} kJ/m²)')
    print('  → 理論上應該不觸發非線性放大')
    print('  → L6D6 的日內累積逆境應該很低')
else:
    print(f'✗ 日累積能量 ({E_elapsed:.1f} kJ/m²) > E_50 ({p.E_50} kJ/m²)')
    print('  → 觸發了非線性放大')
    print('  → 這可能不合理，因為 L6D6 只照射 6 小時/天')

print()

if avg_stress > 0.1:
    print(f'⚠️  平均 Stress = {avg_stress:.4f} 過高')
    print('   L6D6 (低劑量短期) 不應該累積這麼多 Stress')
    print()
    print('可能原因:')
    print('  1. stress_damage_coeff 仍然太高')
    print('  2. 日內非線性參數設定不當')
    print('  3. 修復能力不足')
    print('  4. LAI 脆弱性機制過於敏感')
else:
    print(f'✓ 平均 Stress = {avg_stress:.4f} 較低')

print()
print('=' * 80)
