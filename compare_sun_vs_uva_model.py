#!/usr/bin/env python3
"""
對比 Sun 原始模型 vs UVA 模型對 L6D6 的預測
找出為什麼 UVA 模型預測鮮重這麼低
"""

import numpy as np
from scipy.integrate import solve_ivp
from lettuce_uva_carbon_complete_model import SunParams, sun_derivatives_final
from simulate_uva_model import UVAParams, uva_sun_derivatives
from model_config import ENV_BASE, TREATMENT_CONFIGS, SIMULATION

def simulate_sun_model(c_alpha_value):
    """使用純 Sun 模型模擬 L6D6 (動態增加 UVA 輻射)"""
    p = SunParams()
    p.c_alpha = c_alpha_value

    env = ENV_BASE.copy()
    plant_density = env['plant_density']

    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * 0.05
    X_d_init = dw_init_g / 1000.0 * plant_density
    C_buf_init = X_d_init * 0.1
    LAI_init = 0.5
    y0 = [X_d_init, C_buf_init, LAI_init]

    t_start = 14 * 86400
    t_end = 35 * 86400

    def sun_with_uva(t, y):
        # Day 29-35, 10:00-16:00 增加 22 W/m²
        days = t / 86400
        hour = (t / 3600) % 24

        env_mod = env.copy()
        if 29 <= days < 35 and 10 <= hour < 16:
            env_mod['I_day'] = env['I_day'] + 22

        return sun_derivatives_final(t, y, p, env_mod)

    sol = solve_ivp(sun_with_uva, (t_start, t_end), y0, method='RK45', max_step=60)
    X_d_final = sol.y[0, -1]
    FW_final = (X_d_final / plant_density * 1000) / 0.05
    return FW_final

def simulate_uva_model(c_alpha_value):
    """使用 UVA 模型模擬 L6D6"""
    p = UVAParams()
    p.c_alpha = c_alpha_value

    env = ENV_BASE.copy()
    env.update(TREATMENT_CONFIGS['L6D6'])

    plant_density = env['plant_density']
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * 0.05
    X_d_init = dw_init_g / 1000.0 * plant_density
    C_buf_init = X_d_init * 0.1
    LAI_init = 0.5
    y0 = [X_d_init, C_buf_init, LAI_init, 0.0, 0.0, 0.0]  # 6 個狀態變量

    t_start = 14 * 86400
    t_end = 35 * 86400

    sol = solve_ivp(
        lambda t, y: uva_sun_derivatives(t, y, p, env),
        (t_start, t_end), y0, method='RK45', max_step=60
    )

    X_d_final = sol.y[0, -1]
    Stress_final = sol.y[4, -1]
    FW_final = (X_d_final / plant_density * 1000) / 0.05
    return FW_final, Stress_final, sol.y[0, :]

def main():
    c_alpha = 0.57

    print('=' * 80)
    print('對比 Sun 原始模型 vs UVA 模型 - L6D6 預測')
    print('=' * 80)
    print()
    print(f'使用相同的 c_alpha = {c_alpha}')
    print(f'L6D6: Day 29-35, 10:00-16:00 增加 22 W/m²')
    print()
    print('-' * 80)

    # Sun 模型
    fw_sun = simulate_sun_model(c_alpha)

    # UVA 模型
    fw_uva, stress_uva, xd_trajectory = simulate_uva_model(c_alpha)

    print(f'Sun 原始模型: {fw_sun:.1f} g')
    print(f'UVA 模型:     {fw_uva:.1f} g (Stress={stress_uva:.2f})')
    print(f'實驗值:       91.4 g')
    print()
    print(f'差異: {fw_sun - fw_uva:.1f} g ({(fw_sun - fw_uva) / fw_sun * 100:.1f}%)')
    print()

    if abs(fw_sun - fw_uva) < 1:
        print('✓ 兩個模型預測一致')
    else:
        print('✗ 兩個模型預測不一致！需要找出原因')
        print()
        print('可能的原因:')
        print('  1. UVA 模型的 Stress 機制抑制了生長')
        print('  2. UVA 模型的花青素合成消耗了太多碳')
        print('  3. UVA 模型的修復機制消耗了太多碳')
        print('  4. UVA 模型的 LDMC 效應改變了 DW/FW 比例')
        print()
        print(f'當前 Stress = {stress_uva:.2f}')
        if stress_uva > 1:
            print('  → Stress 可能是主要原因')

        # 檢查 X_d 軌跡
        xd_ratio = xd_trajectory[-1] / xd_trajectory[0]
        print(f'X_d 增長倍數 (UVA模型): {xd_ratio:.2f}x')

if __name__ == "__main__":
    main()
