#!/usr/bin/env python3
"""
v6.0 模型驗證腳本
測試移除 par_conversion_factor 放大效應後的預測結果
"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives
from model_config import ENV_BASE, TREATMENT_CONFIGS, TARGETS, ODE_SETTINGS, SIMULATION

def simulate_treatment(p, treatment_name, treatment_config):
    """模擬單一處理組"""
    # 建立環境
    env = ENV_BASE.copy()
    env.update(treatment_config)

    # 初始條件
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    plant_density = ENV_BASE['plant_density']

    X_d_init = dw_init_g / 1000.0 * plant_density
    C_buf_init = X_d_init * 0.1
    LAI_init = 0.5  # 移植時初始 LAI
    Anth_init = 0.0
    Stress_init = 0.0
    E_stress_init = 0.0

    y0 = [X_d_init, C_buf_init, LAI_init, Anth_init, Stress_init, E_stress_init]

    # 模擬時間（從第14天到第35天，共21天）
    t_start = 14 * 86400
    t_end = 35 * 86400
    t_span = (t_start, t_end)

    # ODE 求解
    sol = solve_ivp(
        fun=lambda t, y: uva_sun_derivatives(t, y, p, env),
        t_span=t_span,
        y0=y0,
        method=ODE_SETTINGS['method'],
        max_step=ODE_SETTINGS['max_step']
    )

    # 提取最終結果
    X_d_final = sol.y[0, -1]
    LAI_final = sol.y[2, -1]
    Anth_final = sol.y[3, -1]
    Stress_final = sol.y[4, -1]
    E_stress_final = sol.y[5, -1]

    # 計算 DW/FW 比例
    ldmc_increase = p.ldmc_stress_sensitivity * Stress_final / (p.K_ldmc + Stress_final + 1e-9)
    dw_fw_ratio = min(p.dw_fw_ratio_base * (1 + ldmc_increase), p.dw_fw_ratio_max)

    # 計算鮮重
    DW_final_g = X_d_final / plant_density * 1000
    FW_final_g = DW_final_g / dw_fw_ratio

    # 計算花青素濃度
    Anth_final_ppm = (Anth_final / X_d_final) * 1e6 if X_d_final > 0 else 0

    return {
        'FW_final_g': FW_final_g,
        'Anth_final_ppm': Anth_final_ppm,
        'Stress_final': Stress_final,
        'E_stress_final': E_stress_final,
        'LAI_final': LAI_final
    }

def main():
    p = UVAParams()

    print('=' * 80)
    print('v6.0 模型驗證 (par_conversion_factor = 1.0)')
    print('=' * 80)
    print()
    print(f'par_conversion_factor = {p.par_conversion_factor}')
    print()
    print('Treatment   FW_sim   FW_exp  FW_Err   Anth_sim  Anth_exp  Anth_Err  Stress')
    print('-' * 80)

    fw_errors = []
    anth_errors = []

    for treatment in ['CK', 'L6D6', 'L6D6-N', 'H12D3', 'VL3D12', 'L6D12']:
        config = TREATMENT_CONFIGS[treatment]
        target = TARGETS[treatment]

        result = simulate_treatment(p, treatment, config)

        fw_sim = result['FW_final_g']
        fw_exp = target['FW']
        fw_err = (fw_sim - fw_exp) / fw_exp * 100

        anth_sim = result['Anth_final_ppm']
        anth_exp = target['Anth']
        anth_err = (anth_sim - anth_exp) / anth_exp * 100

        stress = result['Stress_final']

        fw_errors.append(abs(fw_err))
        anth_errors.append(abs(anth_err))

        print(f'{treatment:8s}    {fw_sim:5.1f}g   {fw_exp:5.1f}g  {fw_err:+5.1f}%    {anth_sim:5.1f}     {anth_exp:5.1f}     {anth_err:+5.1f}%   {stress:5.1f}')

    print()
    print('統計結果:')
    print(f'  FW  - Mean Error: {np.mean(fw_errors):.1f}%, Max Error: {np.max(fw_errors):.1f}%')
    print(f'  Anth - Mean Error: {np.mean(anth_errors):.1f}%, Max Error: {np.max(anth_errors):.1f}%')
    print()

if __name__ == "__main__":
    main()
