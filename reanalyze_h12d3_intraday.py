#!/usr/bin/env python3
"""
重新分析 H12D3 - 當日內動態
正確理解：非線性累積是當日內 Stress 的正反饋，不是跨天累積
"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives
from model_config import ENV_BASE, TREATMENT_CONFIGS, TARGETS, SIMULATION, get_env_for_treatment
from params_config import ALL_PARAMS

def simulate_single_day(treatment_name, day_number, params):
    """模擬單一天的 Stress 動態"""
    p = UVAParams(params)
    env = get_env_for_treatment(treatment_name)

    # 從該天開始的初始狀態（簡化：假設之前已經模擬過）
    # 這裡用簡化的初始值
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    X_d_init = dw_init_g / 1000.0 * ENV_BASE['plant_density'] * 2.5  # Day 32 時已經長大
    C_buf_init = X_d_init * 0.1
    LAI_init = 4.0  # Day 32 時 LAI 較高
    Anth_init = 1e-6
    Stress_init = 0.0  # 每天開始時 Stress = 0

    y0 = [X_d_init, C_buf_init, LAI_init, Anth_init, Stress_init, 0.0]

    # 模擬一天 (24小時)
    t_start = day_number * 86400
    t_end = (day_number + 1) * 86400
    t_eval = np.linspace(t_start, t_end, 100)  # 每天 100 個採樣點

    sol = solve_ivp(
        fun=lambda t, y: uva_sun_derivatives(t, y, p, env),
        t_span=(t_start, t_end),
        y0=y0,
        method='RK45',
        max_step=60,
        dense_output=True,
        t_eval=t_eval
    )

    return sol

print('=' * 80)
print('H12D3 當日內 Stress 動態分析')
print('=' * 80)
print()

print('配置: 22 W/m², 12h/day (06:00-18:00), Day 32-35')
print()

# 測試不同的非線性參數對當日 Stress 峰值的影響
test_configs = [
    {'name': '當前 v6.2', 'stress_nonlinear_coeff': 5.5, 'K_nonlinear': 1.0},
    {'name': '更強非線性', 'stress_nonlinear_coeff': 10.0, 'K_nonlinear': 0.8},
    {'name': '極強非線性', 'stress_nonlinear_coeff': 15.0, 'K_nonlinear': 0.5},
]

for config in test_configs:
    print('=' * 80)
    print(f'{config["name"]}')
    print(f'  stress_nonlinear_coeff = {config["stress_nonlinear_coeff"]}')
    print(f'  K_nonlinear = {config["K_nonlinear"]}')
    print('=' * 80)

    params = ALL_PARAMS.copy()
    params.update({
        'stress_nonlinear_coeff': config['stress_nonlinear_coeff'],
        'K_nonlinear': config['K_nonlinear'],
    })

    # 模擬 Day 32 (UVA 第一天)
    sol = simulate_single_day('H12D3', 32, params)

    # 分析當日動態
    stress_profile = sol.y[4, :]  # Stress
    max_stress = np.max(stress_profile)
    final_stress = stress_profile[-1]  # 一天結束時

    # 找出何時達到最大 Stress
    max_idx = np.argmax(stress_profile)
    max_time = sol.t[max_idx]
    max_hour = (max_time / 3600) % 24

    print(f'\nDay 32 當日動態:')
    print(f'  Stress 峰值: {max_stress:.4f} (在 {max_hour:.1f}:00)')
    print(f'  Stress 終值: {final_stress:.4f} (24:00)')
    print()

    # 顯示關鍵時間點
    print(f'{"Hour":>6s} {"Stress":>10s} {"變化率":>12s}')
    print('-' * 40)

    for i in [0, 20, 40, 60, 80, 99]:  # 選擇代表性時間點
        hour = (sol.t[i] / 3600) % 24
        stress = sol.y[4, i]
        print(f'{hour:6.1f} {stress:10.6f}')

    print()

print('=' * 80)
print('關鍵發現')
print('=' * 80)
print()

print('''
正確理解當日動態：

1. **每天的 Stress 動態**:
   - 06:00 開始照射，Stress 開始累積
   - 18:00 照射結束，Stress 開始修復
   - 24:00 一天結束，Stress 可能部分修復

2. **非線性累積的作用**:
   - 在照射期間 (06:00-18:00)，當 Stress 累積到一定程度
   - nonlinear_factor = 1 + coeff * Stress / (K + Stress)
   - 放大當前的損傷率，加速累積
   - 這是**當日內的正反饋**

3. **關鍵問題**:
   - H12D3 照射 12h/day，比 L6D6 的 6h/day 長
   - 理論上應該累積更多 Stress
   - 但為什麼還是不夠？

4. **可能的原因**:
   - 修復能力太強？
   - 損傷率基礎值太低？
   - 日內非線性 (E_50, k_intraday) 的影響？

需要檢查:
---------
1. 日內能量非線性 (intraday_factor)
2. 修復率 vs 損傷率的平衡
3. LAI 脆弱性 (H12D3 時 LAI 較高，脆弱性降低？)
''')

print()
print('=' * 80)
print('檢查日內能量非線性 (E_50, k_intraday)')
print('=' * 80)
print()

print('H12D3 配置:')
print('  I_UVA = 22 W/m²')
print('  照射時數 = 12 h')
print('  日累積能量 = 22 × 12 × 3.6 = 950.4 kJ/m²')
print()

print('當前參數:')
p = UVAParams()
print(f'  E_50 = {p.E_50} kJ/m²')
print(f'  E_scale = {p.E_scale} kJ/m²')
print(f'  k_intraday = {p.k_intraday}')
print()

# 計算 H12D3 的日內非線性因子
E_elapsed = 22 * 12 * 3.6  # 950.4 kJ/m²
normalized_E = (E_elapsed - p.E_50) / p.E_scale
print(f'正規化能量 = (E - E_50) / E_scale = ({E_elapsed} - {p.E_50}) / {p.E_scale} = {normalized_E:.3f}')

# Softplus
def softplus(x, sharpness):
    if x > 10:  # 避免數值溢出
        return x
    return np.log(1.0 + np.exp(sharpness * x)) / sharpness

sp = softplus(normalized_E, p.sharpness_intraday)
intraday_factor = 1.0 + p.k_intraday * (sp ** p.m_intraday)

print(f'softplus({normalized_E:.3f}, {p.sharpness_intraday}) = {sp:.4f}')
print(f'intraday_factor = 1 + {p.k_intraday} × {sp:.4f}^{p.m_intraday} = {intraday_factor:.4f}')
print()

if normalized_E > 0:
    print('✓ H12D3 觸發了日內非線性！')
    print(f'  損傷率放大 {intraday_factor:.2f}x')
else:
    print('✗ H12D3 沒有觸發日內非線性')
    print('  E_elapsed < E_50，沒有放大效應')

print()
print('=' * 80)
print('對比 L6D6:')
print('=' * 80)

E_l6d6 = 22 * 6 * 3.6  # 475.2 kJ/m²
normalized_E_l6d6 = (E_l6d6 - p.E_50) / p.E_scale
sp_l6d6 = softplus(normalized_E_l6d6, p.sharpness_intraday)
intraday_l6d6 = 1.0 + p.k_intraday * (sp_l6d6 ** p.m_intraday)

print(f'L6D6 日累積能量 = {E_l6d6} kJ/m²')
print(f'L6D6 正規化能量 = {normalized_E_l6d6:.3f}')
print(f'L6D6 intraday_factor = {intraday_l6d6:.4f}')
print()

print(f'H12D3 vs L6D6:')
print(f'  累積能量: {E_elapsed / E_l6d6:.2f}x')
print(f'  日內放大: {intraday_factor / intraday_l6d6:.2f}x')

print()
print('=' * 80)
