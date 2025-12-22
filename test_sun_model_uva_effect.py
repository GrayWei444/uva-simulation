#!/usr/bin/env python3
"""
測試 Sun 原始模型對 UVA 額外輻射的響應
L6D6: Day 29-35 (共6天), 每天 10:00-16:00 (6小時) 增加 22 W/m²
"""

import numpy as np
from scipy.integrate import solve_ivp
from lettuce_uva_carbon_complete_model import SunParams, sun_derivatives_final
from model_config import ENV_BASE, SIMULATION

def sun_derivatives_with_dynamic_I(t, y, p, env, uva_schedule):
    """
    Sun 模型 + 動態輻射調整（模擬 UVA 照射期間的額外輻射）
    """
    # 計算當前時間
    days_since_transplant = t / 86400
    hour = (t / 3600) % 24

    # 判斷是否在 UVA 照射期間
    uva_on = False
    if 'uva_start_day' in uva_schedule and 'uva_end_day' in uva_schedule:
        day_in_range = uva_schedule['uva_start_day'] <= days_since_transplant < uva_schedule['uva_end_day']
        hour_in_range = uva_schedule['uva_hour_on'] <= hour < uva_schedule['uva_hour_off']
        uva_on = day_in_range and hour_in_range

    # 修改環境參數
    env_modified = env.copy()
    if uva_on:
        # UVA 照射期間，增加 22 W/m²
        env_modified['I_day'] = env['I_day'] + 22

    return sun_derivatives_final(t, y, p, env_modified)

def simulate_treatment(p, treatment_name, uva_schedule):
    """模擬單一處理組"""
    env = ENV_BASE.copy()
    plant_density = env['plant_density']

    # 初始條件
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * 0.05
    X_d_init = dw_init_g / 1000.0 * plant_density
    C_buf_init = X_d_init * 0.1
    LAI_init = 0.5
    y0 = [X_d_init, C_buf_init, LAI_init]

    # 模擬時間
    t_start = 14 * 86400
    t_end = 35 * 86400

    # ODE 求解
    sol = solve_ivp(
        fun=lambda t, y: sun_derivatives_with_dynamic_I(t, y, p, env, uva_schedule),
        t_span=(t_start, t_end),
        y0=y0,
        method='RK45',
        max_step=60
    )

    # 計算鮮重
    X_d_final = sol.y[0, -1]
    DW_final_g = X_d_final / plant_density * 1000
    FW_final_g = DW_final_g / 0.05

    return FW_final_g

def main():
    p = SunParams()
    p.c_alpha = 0.57

    print('=' * 70)
    print('Sun 原始模型測試 - 精確模擬 L6D6 的 UVA 額外輻射效應')
    print('=' * 70)
    print()
    print(f'基礎 LED: {ENV_BASE["I_day"]} W/m²')
    print(f'UVA 額外: 22 W/m²')
    print(f'UVA 照射: Day 29-35 (6天), 10:00-16:00 (6小時/天)')
    print()
    print('-' * 70)

    # CK (無 UVA)
    ck_schedule = {}
    fw_ck = simulate_treatment(p, 'CK', ck_schedule)

    # L6D6 (Day 29-35, 10:00-16:00)
    l6d6_schedule = {
        'uva_start_day': 29,
        'uva_end_day': 35,
        'uva_hour_on': 10,
        'uva_hour_off': 16
    }
    fw_l6d6 = simulate_treatment(p, 'L6D6', l6d6_schedule)

    # 計算增益
    fw_increase_g = fw_l6d6 - fw_ck
    fw_increase_pct = (fw_l6d6 - fw_ck) / fw_ck * 100

    print(f'CK 模擬鮮重:   {fw_ck:.1f} g (實驗值: 87.0 g)')
    print(f'L6D6 模擬鮮重: {fw_l6d6:.1f} g (實驗值: 91.4 g)')
    print()
    print(f'UVA 導致的鮮重增加:')
    print(f'  模擬: {fw_increase_g:+.1f} g ({fw_increase_pct:+.1f}%)')
    print(f'  實驗: {91.4 - 87.0:+.1f} g ({(91.4 - 87.0) / 87.0 * 100:+.1f}%)')
    print()

    # 分析
    exp_increase_pct = (91.4 - 87.0) / 87.0 * 100
    ratio = fw_increase_pct / exp_increase_pct if exp_increase_pct > 0 else 0

    print('-' * 70)
    print('分析:')
    print(f'  實驗觀察到的增益: {exp_increase_pct:.1f}%')
    print(f'  Sun 模型預測增益: {fw_increase_pct:.1f}%')
    print(f'  預測/實驗比例: {ratio:.2f}x')
    print()

    if ratio > 0.8 and ratio < 1.2:
        print('  ✓ Sun 模型可以解釋 UVA 的增益效應')
        print('  → par_conversion_factor = 1.0 是合理的')
    elif ratio > 1.2:
        print('  ⚠️ Sun 模型高估了 UVA 的增益效應')
        print(f'  → 建議 par_conversion_factor = {1.0 / ratio:.2f}')
    else:
        print('  ⚠️ Sun 模型低估了 UVA 的增益效應')
        print(f'  → 需要 par_conversion_factor = {1.0 / ratio:.2f} 來匹配實驗')

if __name__ == "__main__":
    main()
