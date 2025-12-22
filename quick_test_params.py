#!/usr/bin/env python3
"""快速測試單一參數組合"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TREATMENT_CONFIGS, TARGETS, ODE_SETTINGS, SIMULATION, get_env_for_treatment
from params_config import ALL_PARAMS

# 測試參數 - 只增強非線性，不提高基礎損傷
test_params = ALL_PARAMS.copy()
test_params.update({
    'stress_nonlinear_coeff': 5.0,  # 從 2.0 大幅提高到 5.0 (高 Stress 時加速累積)
    'K_nonlinear': 1.0,  # 從 2.0 降至 1.0 (更早觸發非線性)
    'stress_photosynthesis_inhibition': 0.65,  # 從 0.50 提高到 0.65
    'K_stress': 2.0,  # 從 3.0 降至 2.0 (抑制更敏感)
})

p = UVAParams(test_params)

print('測試參數:')
print(f'  stress_nonlinear_coeff = {p.stress_nonlinear_coeff}')
print(f'  K_nonlinear = {p.K_nonlinear}')
print(f'  stress_photosynthesis_inhibition = {p.stress_photosynthesis_inhibition}')
print(f'  K_stress = {p.K_stress}')
print()

# 只測試關鍵處理組
treatments = ['CK', 'L6D6', 'H12D3']

for treatment_name in treatments:
    env = get_env_for_treatment(treatment_name)
    target = TARGETS[treatment_name]

    # 初始條件
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    X_d_init = dw_init_g / 1000.0 * ENV_BASE['plant_density']
    C_buf_init = X_d_init * 0.1
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6

    y0 = [X_d_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]

    # 模擬
    t_start = SIMULATION['transplant_offset'] * 86400
    t_end = (SIMULATION['transplant_offset'] + SIMULATION['days']) * 86400

    sol = solve_ivp(
        fun=lambda t, y: uva_sun_derivatives(t, y, p, env),
        t_span=(t_start, t_end),
        y0=y0,
        method='RK45',
        max_step=60,
        dense_output=True
    )

    y_final = sol.y[:, -1]
    X_d, C_buf, LAI, Anth, Stress, E_stress = y_final

    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress, p)
    dw_per_plant = X_d / ENV_BASE['plant_density'] * 1000
    fw_per_plant = dw_per_plant / dw_fw_ratio

    fw_err = (fw_per_plant - target['FW']) / target['FW'] * 100

    print(f'{treatment_name}: {fw_per_plant:.2f}g (目標 {target["FW"]}g, 誤差 {fw_err:+.1f}%)')

print('\n完成')
