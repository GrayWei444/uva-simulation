#!/usr/bin/env python3
"""
Generate Fig 19: Sensitivity Analysis
- Clearer labels and larger fonts
- Better layout
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from simulate_uva_model_v10 import (
    UVAParams, uva_sun_derivatives, nonlinear_damage_factor,
    calculate_dynamic_dw_fw_ratio
)

OUTPUT_DIR = 'paper_figures'

ENV_BASE = {
    'light_on_hour': 6,
    'light_off_hour': 22,
    'I_day': 57,
    'T_day': 25,
    'T_night': 18,
    'CO2_day': 1200,
    'CO2_night': 1200,
    'RH_day': 0.70,
    'RH_night': 0.85,
    'plant_density': 36,
}


def run_simulation_with_params(p, env):
    """Run simulation and return FW, Anth, Stress"""
    fw_init_g = 10
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * env['plant_density']
    C_buf_init = Xd_init * 0.1
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * env['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6

    initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]

    transplant_day = 14
    simulation_days = 21
    t_start = transplant_day * 86400
    harvest_hour = 6
    t_end = (transplant_day + simulation_days) * 86400 + harvest_hour * 3600

    t_eval_points = np.linspace(t_start, t_end, 500)
    sol = solve_ivp(
        uva_sun_derivatives,
        (t_start, t_end),
        initial_state,
        args=(p, env),
        method='RK45',
        max_step=300,
        t_eval=t_eval_points
    )

    if not sol.success:
        return None, None, None

    Xd_f, Cbuf_f, LAI_f, Anth_f, Stress_f, ROS_f = sol.y[:, -1]

    uva_start = env.get('uva_start_day', 35) * 86400
    stress_sum, stress_count = 0.0, 0
    for i in range(len(sol.t)):
        if sol.t[i] >= uva_start:
            stress_sum += sol.y[4, i]
            stress_count += 1
    avg_stress = stress_sum / max(1, stress_count)

    uva_hour_on = env.get('uva_hour_on', 0)
    uva_hour_off = env.get('uva_hour_off', 0)
    hours_daily = uva_hour_off - uva_hour_on if uva_hour_on < uva_hour_off else 24 - uva_hour_on + uva_hour_off
    if not env.get('uva_on', False):
        hours_daily = 0
    nonlin_factor = nonlinear_damage_factor(hours_daily, p)
    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(avg_stress, p, nonlin_factor)

    FW_sim = Xd_f / env['plant_density'] / dw_fw_ratio * 1000
    FW_total_kg = FW_sim / 1000 * env['plant_density']
    Anth_sim = Anth_f / FW_total_kg * 1e6

    return FW_sim, Anth_sim, avg_stress


def generate_fig19():
    """Generate sensitivity analysis figure"""
    print("Generating Fig 19: Sensitivity Analysis...")

    # Baseline: optimal treatment (9h × 5d)
    env_opt = dict(ENV_BASE)
    env_opt['uva_on'] = True
    env_opt['uva_intensity'] = 11.0
    env_opt['uva_start_day'] = 30  # 35 - 5 = 30
    env_opt['uva_end_day'] = 35
    env_opt['uva_hour_on'] = 6
    env_opt['uva_hour_off'] = 15  # 6 + 9 = 15

    p_base = UVAParams()
    fw_base, anth_base, stress_base = run_simulation_with_params(p_base, env_opt)

    # Parameters to test
    param_info = {
        'A_vulnerability': ('LAI Vulnerability (A)', p_base.A_vulnerability),
        'gompertz_threshold': ('Gompertz Threshold (h)', p_base.gompertz_threshold),
        'V_max_anth': ('Anth Synthesis Rate', p_base.V_max_anth),
        'k_deg': ('Anth Degradation Rate', p_base.k_deg),
        'k_stress_decay': ('Stress Decay Rate', p_base.k_stress_decay),
    }

    # Multipliers for sensitivity curve
    multipliers = np.linspace(0.5, 2.0, 20)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # (a) Parameter sensitivity curves - FW
    ax1 = axes[0]
    colors = ['#1f77b4', '#2ca02c', '#d62728', '#9467bd', '#17becf']

    for idx, (param_name, (label, base_val)) in enumerate(param_info.items()):
        fw_changes = []
        for mult in multipliers:
            p_test = UVAParams()
            setattr(p_test, param_name, base_val * mult)
            fw, anth, stress = run_simulation_with_params(p_test, env_opt)
            if fw is not None:
                fw_changes.append((fw - fw_base) / fw_base * 100)
            else:
                fw_changes.append(0)

        ax1.plot(multipliers, fw_changes, '-o', color=colors[idx], linewidth=2,
                 markersize=4, label=label)

    ax1.axhline(y=0, color='gray', linestyle='--', linewidth=1)
    ax1.axvline(x=1.0, color='gray', linestyle='--', linewidth=1)
    ax1.set_xlabel('Parameter Multiplier', fontsize=12)
    ax1.set_ylabel('Fresh Weight Change (%)', fontsize=12)
    ax1.set_title('(a) Parameter Sensitivity on Fresh Weight', fontsize=13, fontweight='bold')
    ax1.legend(loc='best', fontsize=10)
    ax1.set_xlim(0.5, 2.0)
    ax1.grid(True, alpha=0.3)

    # (b) Elasticity heatmap
    ax2 = axes[1]

    # Calculate elasticity (% change in output / % change in input)
    delta = 0.1  # 10% change
    outputs = ['FW', 'Stress', 'Anth']
    elasticity = np.zeros((len(param_info), len(outputs)))

    for i, (param_name, (label, base_val)) in enumerate(param_info.items()):
        # Perturb up
        p_up = UVAParams()
        setattr(p_up, param_name, base_val * (1 + delta))
        fw_up, anth_up, stress_up = run_simulation_with_params(p_up, env_opt)

        # Perturb down
        p_down = UVAParams()
        setattr(p_down, param_name, base_val * (1 - delta))
        fw_down, anth_down, stress_down = run_simulation_with_params(p_down, env_opt)

        if fw_up and fw_down:
            # Elasticity = (Δy/y) / (Δx/x)
            elasticity[i, 0] = ((fw_up - fw_down) / fw_base) / (2 * delta)
            elasticity[i, 1] = ((stress_up - stress_down) / max(stress_base, 0.01)) / (2 * delta)
            elasticity[i, 2] = ((anth_up - anth_down) / anth_base) / (2 * delta)

    im = ax2.imshow(elasticity, cmap='RdBu_r', aspect='auto', vmin=-2, vmax=2)
    ax2.set_xticks(range(len(outputs)))
    ax2.set_xticklabels(outputs, fontsize=11)
    ax2.set_yticks(range(len(param_info)))
    ax2.set_yticklabels([info[0] for info in param_info.values()], fontsize=10)
    ax2.set_title('(b) Elasticity Heatmap', fontsize=13, fontweight='bold')

    # Add text annotations
    for i in range(len(param_info)):
        for j in range(len(outputs)):
            color = 'white' if abs(elasticity[i, j]) > 1 else 'black'
            ax2.text(j, i, f'{elasticity[i, j]:.2f}', ha='center', va='center',
                     fontsize=11, color=color, fontweight='bold')

    cbar = plt.colorbar(im, ax=ax2)
    cbar.set_label('Elasticity', fontsize=11)

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Fig19_sensitivity_analysis.png', dpi=300)
    plt.savefig(f'{OUTPUT_DIR}/Fig19_sensitivity_analysis.pdf')
    plt.close()
    print(f"  Saved: {OUTPUT_DIR}/Fig19_sensitivity_analysis.png/pdf")


if __name__ == '__main__':
    generate_fig19()
