# -*- coding: utf-8 -*-
"""
敏感度分析腳本 (Sensitivity Analysis)
=====================================
版本: v1.0
日期: 2025-12-17
"""

import sys
sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment, ODE_SETTINGS

PARAMS_TO_ANALYZE = {
    'E_50': {'base': 237.6, 'range': 0.30, 'description': '半效能量'},
    'E_scale': {'base': 39.6, 'range': 0.30, 'description': '能量尺度'},
    'stress_damage_coeff': {'base': 3.5e-6, 'range': 0.30, 'description': '損傷係數'},
    'stress_repair_coeff': {'base': 1.0e-5, 'range': 0.30, 'description': '修復係數'},
    'LAI_ref_vuln': {'base': 7.5, 'range': 0.30, 'description': 'LAI參考值'},
    'n_vuln': {'base': 7, 'range': 0.30, 'description': '脆弱性指數'},
    'k_intraday': {'base': 1.5, 'range': 0.30, 'description': '日內非線性係數'},
    'stress_nonlinear_coeff': {'base': 1.5, 'range': 0.30, 'description': 'Stress非線性係數'},
    'circadian_disruption_factor': {'base': 2.0, 'range': 0.50, 'description': '節律損傷因子'},
    'stress_photosynthesis_inhibition': {'base': 0.70, 'range': 0.30, 'description': '光合抑制係數'},
    'K_stress': {'base': 5.0, 'range': 0.30, 'description': '抑制半飽和常數'},
    'par_conversion_factor': {'base': 1.0, 'range': 0.10, 'description': 'UVA-PAR轉換係數'},
}

TEST_TREATMENTS = ['CK', 'L6D6', 'L6D6-N', 'H12D3']

def run_simulation(p, treatment):
    env = get_env_for_treatment(treatment)
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6
    initial_state = [Xd_init, 0.0, LAI_init, Anth_init, 0.0]
    simulation_days = SIMULATION['days']
    try:
        sol = solve_ivp(
            fun=uva_sun_derivatives,
            t_span=(0, simulation_days * 86400),
            y0=initial_state,
            args=(p, env),
            method=ODE_SETTINGS['method'],
            max_step=ODE_SETTINGS['max_step'],
        )
        if sol.success:
            sim_stress = sol.y[4, -1]
            sim_dw_fw_ratio = calculate_dynamic_dw_fw_ratio(sim_stress, p)
            sim_dw_per_plant = sol.y[0, -1] / ENV_BASE['plant_density'] * 1000
            sim_fw_per_plant = sim_dw_per_plant / sim_dw_fw_ratio
            return {'FW': sim_fw_per_plant, 'Stress': sim_stress, 'success': True}
    except:
        pass
    return {'FW': np.nan, 'Stress': np.nan, 'success': False}

def calculate_sensitivity(param_name, param_info, treatments):
    base_value = param_info['base']
    variation = param_info['range']
    low_value = base_value * (1 - variation)
    high_value = base_value * (1 + variation)
    results = {'param': param_name, 'treatments': {}}
    for treatment in treatments:
        p_base = UVAParams()
        base_result = run_simulation(p_base, treatment)
        p_low = UVAParams()
        setattr(p_low, param_name, low_value)
        low_result = run_simulation(p_low, treatment)
        p_high = UVAParams()
        setattr(p_high, param_name, high_value)
        high_result = run_simulation(p_high, treatment)
        if all([base_result['success'], low_result['success'], high_result['success']]):
            fw_base = base_result['FW']
            fw_low = low_result['FW']
            fw_high = high_result['FW']
            delta_x_ratio = 2 * variation
            delta_fw_ratio = (fw_high - fw_low) / fw_base if fw_base > 0 else 0
            sensitivity_fw = delta_fw_ratio / delta_x_ratio if delta_x_ratio > 0 else 0
            stress_base = base_result['Stress']
            stress_low = low_result['Stress']
            stress_high = high_result['Stress']
            delta_stress_ratio = (stress_high - stress_low) / (stress_base + 0.01)
            sensitivity_stress = delta_stress_ratio / delta_x_ratio if delta_x_ratio > 0 else 0
            results['treatments'][treatment] = {
                'FW_base': fw_base, 'FW_low': fw_low, 'FW_high': fw_high, 'FW_sensitivity': sensitivity_fw,
                'Stress_base': stress_base, 'Stress_low': stress_low, 'Stress_high': stress_high, 'Stress_sensitivity': sensitivity_stress,
            }
        else:
            results['treatments'][treatment] = None
    return results

print("="*70)
print("敏感度分析 - 萵苣UVA模型 v5.7")
print("="*70)
print()
print(f"分析參數數量: {len(PARAMS_TO_ANALYZE)}")
print(f"測試處理組: {TEST_TREATMENTS}")
print()

all_results = {}
for param_name, param_info in PARAMS_TO_ANALYZE.items():
    print(f"分析參數: {param_name} ({param_info['description']})...")
    results = calculate_sensitivity(param_name, param_info, TEST_TREATMENTS)
    all_results[param_name] = results

print()
print("="*90)
print("敏感度分析結果 (S = 彈性係數)")
print("="*90)

for treatment in TEST_TREATMENTS:
    print(f"
[{treatment}]")
    print("-"*80)
    print(f'{"參數":<32} {"FW敏感度":>10} {"Stress敏感度":>12} {"FW範圍":>18}')
    print("-"*80)
    sensitivities = []
    for param_name, results in all_results.items():
        if treatment in results['treatments'] and results['treatments'][treatment]:
            data = results['treatments'][treatment]
            fw_sens = data['FW_sensitivity']
            stress_sens = data['Stress_sensitivity']
            fw_range = f"{data['FW_low']:.1f}-{data['FW_high']:.1f}g"
            print(f'{param_name:<32} {fw_sens:>10.3f} {stress_sens:>12.3f} {fw_range:>18}')
            sensitivities.append({'param': param_name, 'fw_sens': abs(fw_sens), 'stress_sens': abs(stress_sens)})
    print()
    sorted_by_fw = sorted(sensitivities, key=lambda x: x['fw_sens'], reverse=True)[:3]
    print(f"  FW最敏感: " + ", ".join([f"{s['param']}({s['fw_sens']:.3f})" for s in sorted_by_fw]))

print()
print("="*70)
print("敏感度分析完成！")
print("="*70)
