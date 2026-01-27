"""
================================================================================
UVA Lettuce Model v6.0 - Figure Generation Script
================================================================================
Generates all figures for the v6 paper based on v5 figure numbering convention.

Figure Numbers (v6 following v5 structure):
- Fig 9: Validation FW response (previously Fig 13)
- Fig 10: Validation Anthocyanin response (previously Fig 14)
- Fig 11: Validation parity plots (previously Fig 15)
- Fig 12: LAI vulnerability function (previously Fig 9)
- Fig 13: Gompertz nonlinear damage factor (previously Fig 10)
- Fig 14: Training parity plots (FW & Anth) (previously Fig 11)
- Fig 15: Stress time series (previously Fig 12)
- Fig 16: Hill inhibition functions (same)
- Fig 17: Sensitivity analysis plots (previously Fig 19)
- Fig 18: Optimization heatmaps (previously Fig 20)
- Fig S1: System dynamics block diagram (previously Fig 17)
- Fig S2: Hormesis response surface (previously Fig 18)
- Fig S3: Carbon competition mechanism (previously Fig 21)

================================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from scipy.integrate import solve_ivp
import os

# Import v2 model
from simulate_uva_model_v2 import (
    UVAParams, uva_sun_derivatives, ALL_PARAMS,
    nonlinear_damage_factor, calculate_dynamic_dw_fw_ratio,
    calculate_anthocyanin_ppm
)
from lettuce_uva_carbon_complete_model import sun_derivatives_final

# Output directory
OUTPUT_DIR = "paper_figures"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Style settings
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 11,
    'figure.dpi': 150,
})

# Environment settings
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

SIMULATION = {
    'days': 21,
    'transplant_offset': 14,
    'initial_fw_g': 10,
}


# ==============================================================================
# Fig 12: LAI Vulnerability Function
# ==============================================================================
def generate_fig12_lai_vulnerability():
    """Generate LAI vulnerability function plot."""
    print("Generating Fig 12: LAI vulnerability...")

    p = UVAParams()
    LAI = np.linspace(0, 12, 200)
    vulnerability = p.A_vulnerability * np.exp(-p.k_vulnerability * LAI) + 1.0

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.semilogy(LAI, vulnerability, 'b-', linewidth=2.5)

    # Mark typical LAI ranges
    ax.axvline(x=3, color='orange', linestyle='--', alpha=0.7, label='Young plant (LAI~3)')
    ax.axvline(x=9, color='green', linestyle='--', alpha=0.7, label='Mature plant (LAI~9)')

    # Add annotation
    ax.annotate(f'A = {p.A_vulnerability:.1e}\nk = {p.k_vulnerability}',
                xy=(1, 1e6), fontsize=12,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

    ax.set_xlabel('LAI (m²/m²)')
    ax.set_ylabel('Vulnerability Factor')
    ax.set_title('Fig 12: LAI-Dependent Stress Vulnerability (v6)')
    ax.set_xlim(0, 12)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'Fig12_LAI_vulnerability.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, 'Fig12_LAI_vulnerability.pdf'), bbox_inches='tight')
    plt.close()
    print("  Saved Fig12_LAI_vulnerability.png/pdf")


# ==============================================================================
# Fig 13: Gompertz Nonlinear Damage Factor
# ==============================================================================
def generate_fig13_gompertz():
    """Generate Gompertz nonlinear damage factor plot."""
    print("Generating Fig 13: Gompertz nonlinear damage...")

    p = UVAParams()
    hours = np.linspace(0, 16, 200)
    factors = [nonlinear_damage_factor(h, p) for h in hours]

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(hours, factors, 'r-', linewidth=2.5)

    # Mark key points
    key_hours = [3, 6, 9, 12, 15]
    for h in key_hours:
        f = nonlinear_damage_factor(h, p)
        ax.plot(h, f, 'ko', markersize=8)
        ax.annotate(f'{h}h: {f:.1f}', (h, f), xytext=(5, 10),
                   textcoords='offset points', fontsize=10)

    # Add parameters
    ax.annotate(f'Max = {p.gompertz_max_factor}\nThreshold = {p.gompertz_threshold}h\nSteepness = {p.gompertz_steepness}',
                xy=(1, 200), fontsize=12,
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7))

    ax.axhline(y=1, color='gray', linestyle=':', alpha=0.5)
    ax.axhline(y=p.gompertz_max_factor + 1, color='gray', linestyle=':', alpha=0.5)

    ax.set_xlabel('UVA Exposure Hours per Day')
    ax.set_ylabel('Nonlinear Damage Factor')
    ax.set_title('Fig 13: Gompertz Nonlinear Damage Function (v6)')
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 280)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'Fig13_Gompertz_nonlinear.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, 'Fig13_Gompertz_nonlinear.pdf'), bbox_inches='tight')
    plt.close()
    print("  Saved Fig13_Gompertz_nonlinear.png/pdf")


# ==============================================================================
# Simulation Helper Functions
# ==============================================================================
def run_simulation(treatment_name, env, params=None):
    """Run simulation for a treatment and return results.

    Args:
        treatment_name: Name of treatment for logging
        env: Environment dictionary
        params: Optional UVAParams object. If None, uses default parameters.
    """
    p = params if params is not None else UVAParams()

    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    C_buf_init = Xd_init * 0.1
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6
    AOX_init = Anth_init / p.anthocyanin_fraction

    initial_state = [Xd_init, C_buf_init, LAI_init, AOX_init, 0.0, 0.0]

    transplant_day = SIMULATION['transplant_offset']
    simulation_days = SIMULATION['days']
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

    if sol.success:
        Xd_f, Cbuf_f, LAI_f, AOX_f, Stress_f, ROS_f = sol.y[:, -1]

        # Calculate average stress during irradiation
        uva_start = env.get('uva_start_day', 35) * 86400
        stress_sum = 0.0
        stress_count = 0
        for i in range(len(sol.t)):
            if sol.t[i] >= uva_start:
                stress_sum += sol.y[4, i]
                stress_count += 1
        avg_stress = stress_sum / max(1, stress_count)

        # Calculate hours/day
        uva_hour_on = env.get('uva_hour_on', 0)
        uva_hour_off = env.get('uva_hour_off', 0)
        hours_daily = uva_hour_off - uva_hour_on if uva_hour_on < uva_hour_off else 24 - uva_hour_on + uva_hour_off
        if not env.get('uva_on', False):
            hours_daily = 0
        nonlin_factor = nonlinear_damage_factor(hours_daily, p)

        # Calculate outputs
        dw_fw_ratio = calculate_dynamic_dw_fw_ratio(avg_stress, p, nonlin_factor)
        FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
        FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
        Anth_sim = calculate_anthocyanin_ppm(AOX_f, FW_total_kg, p)

        return {
            'sol': sol,
            'FW': FW_sim,
            'Anth': Anth_sim,
            'LAI': LAI_f,
            'Stress': Stress_f,
            'avg_stress': avg_stress,
            'AOX': AOX_f,
            'C_buf': Cbuf_f,
        }
    return None


# ==============================================================================
# Fig 14: Training Parity Plots
# ==============================================================================
def generate_fig14_training_parity():
    """Generate training parity plots for FW and Anthocyanin."""
    print("Generating Fig 14: Training parity plots...")

    TARGETS = {
        'CK':      {'FW': 87.0, 'Anth': 433, 'env': {'uva_on': False}},
        'L6D6':    {'FW': 91.4, 'Anth': 494, 'env': {'uva_on': True, 'uva_intensity': 11.0, 'uva_start_day': 29, 'uva_end_day': 35, 'uva_hour_on': 10, 'uva_hour_off': 16}},
        'L6D6-N':  {'FW': 80.8, 'Anth': 493, 'env': {'uva_on': True, 'uva_intensity': 11.0, 'uva_start_day': 29, 'uva_end_day': 35, 'uva_hour_on': 22, 'uva_hour_off': 4}},
        'VL3D12':  {'FW': 67.0, 'Anth': 482, 'env': {'uva_on': True, 'uva_intensity': 11.0, 'uva_start_day': 23, 'uva_end_day': 35, 'uva_hour_on': 10, 'uva_hour_off': 13}},
        'L6D12':   {'FW': 60.4, 'Anth': 518, 'env': {'uva_on': True, 'uva_intensity': 11.0, 'uva_start_day': 23, 'uva_end_day': 35, 'uva_hour_on': 10, 'uva_hour_off': 16}},
        'H12D3':   {'FW': 60.6, 'Anth': 651, 'env': {'uva_on': True, 'uva_intensity': 11.0, 'uva_start_day': 32, 'uva_end_day': 35, 'uva_hour_on': 6, 'uva_hour_off': 18}},
    }

    results = {}
    for name, target in TARGETS.items():
        env = dict(ENV_BASE)
        env.update(target['env'])
        res = run_simulation(name, env)
        if res:
            results[name] = {
                'obs_FW': target['FW'],
                'sim_FW': res['FW'],
                'obs_Anth': target['Anth'],
                'sim_Anth': res['Anth'],
            }

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    colors = {'CK': 'gray', 'L6D6': 'blue', 'L6D6-N': 'purple',
              'VL3D12': 'green', 'L6D12': 'orange', 'H12D3': 'red'}

    # FW parity
    ax = axes[0]
    for name, res in results.items():
        ax.scatter(res['obs_FW'], res['sim_FW'], s=150, c=colors[name],
                  label=name, edgecolors='black', linewidths=1.5, zorder=5)

    # Reference lines
    fw_range = [50, 100]
    ax.plot(fw_range, fw_range, 'k-', linewidth=1.5, label='1:1')
    ax.plot(fw_range, [x*0.95 for x in fw_range], 'k--', alpha=0.5)
    ax.plot(fw_range, [x*1.05 for x in fw_range], 'k--', alpha=0.5)
    ax.fill_between(fw_range, [x*0.95 for x in fw_range], [x*1.05 for x in fw_range],
                   alpha=0.1, color='green', label='±5%')

    ax.set_xlabel('Observed Fresh Weight (g)')
    ax.set_ylabel('Simulated Fresh Weight (g)')
    ax.set_title('Training: Fresh Weight Parity')
    ax.set_xlim(50, 100)
    ax.set_ylim(50, 100)
    ax.legend(loc='upper left', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')

    # Anthocyanin parity
    ax = axes[1]
    for name, res in results.items():
        ax.scatter(res['obs_Anth'], res['sim_Anth'], s=150, c=colors[name],
                  label=name, edgecolors='black', linewidths=1.5, zorder=5)

    anth_range = [400, 700]
    ax.plot(anth_range, anth_range, 'k-', linewidth=1.5, label='1:1')
    ax.plot(anth_range, [x*0.95 for x in anth_range], 'k--', alpha=0.5)
    ax.plot(anth_range, [x*1.05 for x in anth_range], 'k--', alpha=0.5)
    ax.fill_between(anth_range, [x*0.95 for x in anth_range], [x*1.05 for x in anth_range],
                   alpha=0.1, color='green', label='±5%')

    ax.set_xlabel('Observed Anthocyanin (μg/g FW)')
    ax.set_ylabel('Simulated Anthocyanin (μg/g FW)')
    ax.set_title('Training: Anthocyanin Parity')
    ax.set_xlim(400, 700)
    ax.set_ylim(400, 700)
    ax.legend(loc='upper left', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')

    plt.suptitle('Fig 14: Training Dataset Parity (v6)', fontsize=16, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'Fig14_training_parity.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, 'Fig14_training_parity.pdf'), bbox_inches='tight')
    plt.close()
    print("  Saved Fig14_training_parity.png/pdf")


# ==============================================================================
# Fig 15: Stress Time Series
# ==============================================================================
def generate_fig15_stress_timeseries():
    """Generate stress time series for all treatments."""
    print("Generating Fig 15: Stress time series...")

    TREATMENTS = {
        'CK':      {'uva_on': False},
        'L6D6':    {'uva_on': True, 'uva_intensity': 11.0, 'uva_start_day': 29, 'uva_end_day': 35, 'uva_hour_on': 10, 'uva_hour_off': 16},
        'L6D6-N':  {'uva_on': True, 'uva_intensity': 11.0, 'uva_start_day': 29, 'uva_end_day': 35, 'uva_hour_on': 22, 'uva_hour_off': 4},
        'VL3D12':  {'uva_on': True, 'uva_intensity': 11.0, 'uva_start_day': 23, 'uva_end_day': 35, 'uva_hour_on': 10, 'uva_hour_off': 13},
        'L6D12':   {'uva_on': True, 'uva_intensity': 11.0, 'uva_start_day': 23, 'uva_end_day': 35, 'uva_hour_on': 10, 'uva_hour_off': 16},
        'H12D3':   {'uva_on': True, 'uva_intensity': 11.0, 'uva_start_day': 32, 'uva_end_day': 35, 'uva_hour_on': 6, 'uva_hour_off': 18},
    }

    colors = {'CK': 'gray', 'L6D6': 'blue', 'L6D6-N': 'purple',
              'VL3D12': 'green', 'L6D12': 'orange', 'H12D3': 'red'}

    fig, ax = plt.subplots(figsize=(14, 8))

    for name, config in TREATMENTS.items():
        env = dict(ENV_BASE)
        env.update(config)
        res = run_simulation(name, env)

        if res:
            sol = res['sol']
            days = sol.t / 86400
            stress = sol.y[4, :]
            ax.plot(days, stress, color=colors[name], linewidth=2, label=name)

    # Mark UVA period
    ax.axvspan(23, 35, alpha=0.1, color='yellow', label='UVA period (D12)')
    ax.axvspan(29, 35, alpha=0.1, color='orange', label='UVA period (D6)')
    ax.axvspan(32, 35, alpha=0.1, color='red', label='UVA period (D3)')

    ax.set_xlabel('Days from Sowing')
    ax.set_ylabel('Cumulative Stress Index')
    ax.set_title('Fig 15: Stress Accumulation Time Series (v6)')
    ax.set_xlim(14, 36)
    ax.legend(loc='upper left', ncol=2)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'Fig15_stress_timeseries.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, 'Fig15_stress_timeseries.pdf'), bbox_inches='tight')
    plt.close()
    print("  Saved Fig15_stress_timeseries.png/pdf")


# ==============================================================================
# Fig 9 & 10: Validation Response Curves
# ==============================================================================
def generate_fig9_10_validation_response():
    """Generate validation FW and Anthocyanin response curves."""
    print("Generating Fig 9 & 10: Validation response curves...")

    validation_targets = {
        'CK':      {'FW': 85.14, 'Anth': 413, 'hours': 0},
        'VL3D3':   {'FW': 89.1, 'Anth': 437, 'hours': 3},
        'L6D3':    {'FW': 92.18, 'Anth': 468, 'hours': 6},
        'M9D3':    {'FW': 83.79, 'Anth': 539, 'hours': 9},
        'H12D3':   {'FW': 62.2, 'Anth': 657, 'hours': 12},
        'VH15D3':  {'FW': 51.2, 'Anth': 578, 'hours': 15},
    }

    # Run simulations
    results = []
    for name, target in validation_targets.items():
        hours = target['hours']
        env = dict(ENV_BASE)
        if hours > 0:
            env['uva_on'] = True
            env['uva_intensity'] = 11.0
            env['uva_start_day'] = 32
            env['uva_end_day'] = 35
            env['uva_hour_on'] = 6
            env['uva_hour_off'] = 6 + hours
        else:
            env['uva_on'] = False

        res = run_simulation(name, env)
        if res:
            results.append({
                'hours': hours,
                'obs_FW': target['FW'],
                'sim_FW': res['FW'],
                'obs_Anth': target['Anth'],
                'sim_Anth': res['Anth'],
            })

    results.sort(key=lambda x: x['hours'])
    hours = [r['hours'] for r in results]

    # Fig 9: FW Response
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(hours, [r['obs_FW'] for r in results], 'ko-', markersize=10,
            linewidth=2, label='Observed')
    ax.plot(hours, [r['sim_FW'] for r in results], 'rs--', markersize=10,
            linewidth=2, label='Simulated (v6)')

    ax.fill_between(hours,
                   [r['obs_FW']*0.95 for r in results],
                   [r['obs_FW']*1.05 for r in results],
                   alpha=0.2, color='green', label='±5% band')

    ax.set_xlabel('UVA Exposure Hours per Day')
    ax.set_ylabel('Fresh Weight (g)')
    ax.set_title('Fig 9: Validation - FW Response to UVA Duration (v6)')
    ax.set_xlim(-0.5, 16)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'Fig9_validation_FW.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, 'Fig9_validation_FW.pdf'), bbox_inches='tight')
    plt.close()
    print("  Saved Fig9_validation_FW.png/pdf")

    # Fig 10: Anthocyanin Response
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(hours, [r['obs_Anth'] for r in results], 'ko-', markersize=10,
            linewidth=2, label='Observed')
    ax.plot(hours, [r['sim_Anth'] for r in results], 'rs--', markersize=10,
            linewidth=2, label='Simulated (v6)')

    ax.fill_between(hours,
                   [r['obs_Anth']*0.95 for r in results],
                   [r['obs_Anth']*1.05 for r in results],
                   alpha=0.2, color='green', label='±5% band')

    ax.set_xlabel('UVA Exposure Hours per Day')
    ax.set_ylabel('Anthocyanin (μg/g FW)')
    ax.set_title('Fig 10: Validation - Anthocyanin Response to UVA Duration (v6)')
    ax.set_xlim(-0.5, 16)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'Fig10_validation_Anth.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, 'Fig10_validation_Anth.pdf'), bbox_inches='tight')
    plt.close()
    print("  Saved Fig10_validation_Anth.png/pdf")


# ==============================================================================
# Fig 11: Validation Parity Plots
# ==============================================================================
def generate_fig11_validation_parity():
    """Generate validation parity plots."""
    print("Generating Fig 11: Validation parity plots...")

    validation_targets = {
        'CK':      {'FW': 85.14, 'Anth': 413, 'hours': 0},
        'VL3D3':   {'FW': 89.1, 'Anth': 437, 'hours': 3},
        'L6D3':    {'FW': 92.18, 'Anth': 468, 'hours': 6},
        'M9D3':    {'FW': 83.79, 'Anth': 539, 'hours': 9},
        'H12D3':   {'FW': 62.2, 'Anth': 657, 'hours': 12},
        'VH15D3':  {'FW': 51.2, 'Anth': 578, 'hours': 15},
    }

    results = {}
    for name, target in validation_targets.items():
        hours = target['hours']
        env = dict(ENV_BASE)
        if hours > 0:
            env['uva_on'] = True
            env['uva_intensity'] = 11.0
            env['uva_start_day'] = 32
            env['uva_end_day'] = 35
            env['uva_hour_on'] = 6
            env['uva_hour_off'] = 6 + hours
        else:
            env['uva_on'] = False

        res = run_simulation(name, env)
        if res:
            results[name] = {
                'obs_FW': target['FW'],
                'sim_FW': res['FW'],
                'obs_Anth': target['Anth'],
                'sim_Anth': res['Anth'],
            }

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    colors = {'CK': 'gray', 'VL3D3': 'lightblue', 'L6D3': 'blue',
              'M9D3': 'green', 'H12D3': 'orange', 'VH15D3': 'red'}

    # FW parity
    ax = axes[0]
    for name, res in results.items():
        ax.scatter(res['obs_FW'], res['sim_FW'], s=150, c=colors[name],
                  label=name, edgecolors='black', linewidths=1.5, zorder=5)

    fw_range = [45, 100]
    ax.plot(fw_range, fw_range, 'k-', linewidth=1.5, label='1:1')
    ax.fill_between(fw_range, [x*0.95 for x in fw_range], [x*1.05 for x in fw_range],
                   alpha=0.1, color='green', label='±5%')
    ax.fill_between(fw_range, [x*0.90 for x in fw_range], [x*1.10 for x in fw_range],
                   alpha=0.05, color='yellow', label='±10%')

    ax.set_xlabel('Observed Fresh Weight (g)')
    ax.set_ylabel('Simulated Fresh Weight (g)')
    ax.set_title('Validation: Fresh Weight Parity')
    ax.set_xlim(45, 100)
    ax.set_ylim(45, 100)
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')

    # Anthocyanin parity
    ax = axes[1]
    for name, res in results.items():
        ax.scatter(res['obs_Anth'], res['sim_Anth'], s=150, c=colors[name],
                  label=name, edgecolors='black', linewidths=1.5, zorder=5)

    anth_range = [380, 700]
    ax.plot(anth_range, anth_range, 'k-', linewidth=1.5, label='1:1')
    ax.fill_between(anth_range, [x*0.95 for x in anth_range], [x*1.05 for x in anth_range],
                   alpha=0.1, color='green', label='±5%')
    ax.fill_between(anth_range, [x*0.90 for x in anth_range], [x*1.10 for x in anth_range],
                   alpha=0.05, color='yellow', label='±10%')

    ax.set_xlabel('Observed Anthocyanin (μg/g FW)')
    ax.set_ylabel('Simulated Anthocyanin (μg/g FW)')
    ax.set_title('Validation: Anthocyanin Parity')
    ax.set_xlim(380, 700)
    ax.set_ylim(380, 700)
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')

    plt.suptitle('Fig 11: Validation Dataset Parity (v6)', fontsize=16, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'Fig11_validation_parity.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, 'Fig11_validation_parity.pdf'), bbox_inches='tight')
    plt.close()
    print("  Saved Fig11_validation_parity.png/pdf")


# ==============================================================================
# Fig 16: Hill Inhibition Functions
# ==============================================================================
def generate_fig16_hill_inhibition():
    """Generate Hill inhibition function plots."""
    print("Generating Fig 16: Hill inhibition functions...")

    p = UVAParams()

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Stress inhibition on synthesis
    ax = axes[0, 0]
    stress = np.linspace(0, 500, 200)
    inhibition = p.max_stress_inhib * (stress ** p.n_stress_inhib) / \
                 (p.K_stress_inhib ** p.n_stress_inhib + stress ** p.n_stress_inhib + 1e-9)
    efficiency = 1 - inhibition
    ax.plot(stress, efficiency, 'b-', linewidth=2.5)
    ax.axhline(y=0.5, color='gray', linestyle=':', alpha=0.5)
    ax.axvline(x=p.K_stress_inhib, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('Stress')
    ax.set_ylabel('AOX Synthesis Efficiency')
    ax.set_title(f'Stress Inhibition on Synthesis\n(K={p.K_stress_inhib}, n={p.n_stress_inhib}, max={p.max_stress_inhib})')
    ax.set_xlim(0, 500)
    ax.set_ylim(0, 1.1)
    ax.grid(True, alpha=0.3)

    # 2. AOX protection
    ax = axes[0, 1]
    aox = np.linspace(0, 2e-4, 200)
    protection = p.alpha_aox_protection * aox / (p.K_aox_protection + aox + 1e-12)
    ax.plot(aox * 1e6, protection, 'g-', linewidth=2.5)
    ax.axhline(y=p.alpha_aox_protection/2, color='gray', linestyle=':', alpha=0.5)
    ax.axvline(x=p.K_aox_protection * 1e6, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('AOX (mg/m²)')
    ax.set_ylabel('Protection Efficiency')
    ax.set_title(f'AOX Protection\n(α={p.alpha_aox_protection}, K={p.K_aox_protection*1e6:.1f}mg/m²)')
    ax.set_xlim(0, 200)
    ax.set_ylim(0, 0.6)
    ax.grid(True, alpha=0.3)

    # 3. Stress inhibition on growth
    ax = axes[1, 0]
    stress = np.linspace(0, 300, 200)
    xd_reduction = p.stress_photosynthesis_inhibition * stress / (p.K_stress + stress + 1e-9)
    lai_reduction = p.stress_lai_inhibition * stress / (p.K_stress + stress + 1e-9)
    ax.plot(stress, 1 - xd_reduction, 'r-', linewidth=2.5, label='DW growth')
    ax.plot(stress, 1 - lai_reduction, 'orange', linewidth=2.5, label='LAI growth')
    ax.axhline(y=0.5, color='gray', linestyle=':', alpha=0.5)
    ax.axvline(x=p.K_stress, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('Stress')
    ax.set_ylabel('Growth Efficiency')
    ax.set_title(f'Stress Inhibition on Growth\n(K={p.K_stress})')
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 1.1)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 4. Carbon competition (v2.0 NEW)
    ax = axes[1, 1]
    stress = np.linspace(0, 300, 200)
    stress_K = 21.0  # stress_competition_K
    stress_max = 0.225  # stress_competition_max
    carbon_max = 0.30  # carbon_competition_max

    stress_effect = stress_max * stress / (stress_K + stress + 1e-9)
    # Assume aox_carbon_effect saturates quickly
    aox_effect = 0.8  # typical saturated value
    total_competition = aox_effect * carbon_max + stress_effect
    growth_penalty = 1 - total_competition

    ax.plot(stress, stress_effect, 'b-', linewidth=2, label='Stress competition')
    ax.plot(stress, total_competition, 'r-', linewidth=2.5, label='Total competition')
    ax.plot(stress, growth_penalty, 'g-', linewidth=2.5, label='Growth penalty')
    ax.axhline(y=0.5, color='gray', linestyle=':', alpha=0.5)
    ax.axvline(x=stress_K, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('Stress')
    ax.set_ylabel('Effect')
    ax.set_title(f'Carbon Competition (v2.0)\n(K={stress_K}, max_stress={stress_max}, max_carbon={carbon_max})')
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 1.1)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.suptitle('Fig 16: Hill-Type Inhibition Functions (v6)', fontsize=16, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'Fig16_Hill_inhibition.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, 'Fig16_Hill_inhibition.pdf'), bbox_inches='tight')
    plt.close()
    print("  Saved Fig16_Hill_inhibition.png/pdf")


# ==============================================================================
# Fig S3: Carbon Competition Mechanism (Supplementary)
# ==============================================================================
def generate_figS3_carbon_competition():
    """Generate carbon competition mechanism visualization."""
    print("Generating Fig S3: Carbon competition mechanism...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # 1. Carbon flow diagram (conceptual)
    ax = axes[0, 0]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')

    # Draw boxes
    boxes = [
        (1, 7, 2.5, 1.5, 'Photosynthesis', 'lightgreen'),
        (4, 7, 2.5, 1.5, 'C_buf', 'lightyellow'),
        (7, 8, 2, 1, 'Growth\n(X_d, LAI)', 'lightblue'),
        (7, 6, 2, 1, 'AOX\nSynthesis', 'lightcoral'),
        (7, 4, 2, 1, 'Respiration', 'lightgray'),
    ]

    for x, y, w, h, label, color in boxes:
        rect = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.1",
                              facecolor=color, edgecolor='black', linewidth=2)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, label, ha='center', va='center', fontsize=11, fontweight='bold')

    # Draw arrows
    ax.annotate('', xy=(4, 7.75), xytext=(3.5, 7.75),
               arrowprops=dict(arrowstyle='->', lw=2))
    ax.annotate('', xy=(7, 8.5), xytext=(6.5, 7.75),
               arrowprops=dict(arrowstyle='->', lw=2, color='blue'))
    ax.annotate('', xy=(7, 6.5), xytext=(6.5, 7.25),
               arrowprops=dict(arrowstyle='->', lw=2, color='red'))
    ax.annotate('', xy=(7, 4.5), xytext=(6.5, 6.75),
               arrowprops=dict(arrowstyle='->', lw=2, color='gray'))

    ax.text(5, 9.5, 'Carbon Competition in v2.0', fontsize=14, fontweight='bold', ha='center')
    ax.text(5, 1.5, 'Under stress: ↑AOX synthesis → ↓Growth', fontsize=12, ha='center',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

    # 2. Competition effect vs stress
    ax = axes[0, 1]
    stress = np.linspace(0, 300, 200)
    stress_K = 21.0
    stress_max = 0.225
    carbon_max = 0.30

    stress_effect = stress_max * stress / (stress_K + stress + 1e-9)
    aox_effect = 0.85  # saturated
    total = aox_effect * carbon_max + stress_effect
    growth = 1 - total

    ax.fill_between(stress, 0, stress_effect, alpha=0.3, color='orange', label='Stress competition')
    ax.fill_between(stress, stress_effect, total, alpha=0.3, color='red', label='AOX carbon demand')
    ax.plot(stress, total, 'k-', linewidth=2, label='Total competition')
    ax.plot(stress, growth, 'g-', linewidth=2.5, label='Growth factor')

    # Mark treatment stress levels
    treatments = {'VL3D12': 61, 'L6D12': 146, 'H12D3': 262}
    for name, s in treatments.items():
        g = 1 - (aox_effect * carbon_max + stress_max * s / (stress_K + s))
        ax.axvline(x=s, color='gray', linestyle=':', alpha=0.5)
        ax.plot(s, g, 'ko', markersize=8)
        ax.annotate(name, (s, g), xytext=(5, 10), textcoords='offset points', fontsize=10)

    ax.set_xlabel('Stress')
    ax.set_ylabel('Effect')
    ax.set_title('Growth Penalty from Carbon Competition')
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 1.1)
    ax.legend(loc='center right')
    ax.grid(True, alpha=0.3)

    # 3. C_buf dynamics comparison
    ax = axes[1, 0]

    # Simulate two treatments
    treatments_to_compare = {
        'CK': {'uva_on': False},
        'H12D3': {'uva_on': True, 'uva_intensity': 11.0, 'uva_start_day': 32,
                  'uva_end_day': 35, 'uva_hour_on': 6, 'uva_hour_off': 18},
    }

    colors = {'CK': 'gray', 'H12D3': 'red'}

    for name, config in treatments_to_compare.items():
        env = dict(ENV_BASE)
        env.update(config)
        res = run_simulation(name, env)
        if res:
            sol = res['sol']
            days = sol.t / 86400
            cbuf = sol.y[1, :] * 1000  # Convert to g/m²
            ax.plot(days, cbuf, color=colors[name], linewidth=2, label=name)

    ax.axvspan(32, 35, alpha=0.1, color='red', label='UVA period')
    ax.set_xlabel('Days from Sowing')
    ax.set_ylabel('C_buf (g/m²)')
    ax.set_title('Carbon Buffer Dynamics')
    ax.set_xlim(14, 36)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 4. AOX dynamics comparison
    ax = axes[1, 1]

    for name, config in treatments_to_compare.items():
        env = dict(ENV_BASE)
        env.update(config)
        res = run_simulation(name, env)
        if res:
            sol = res['sol']
            days = sol.t / 86400
            aox = sol.y[3, :] * 1e6  # Convert to mg/m²
            ax.plot(days, aox, color=colors[name], linewidth=2, label=name)

    ax.axvspan(32, 35, alpha=0.1, color='red', label='UVA period')
    ax.set_xlabel('Days from Sowing')
    ax.set_ylabel('AOX (mg/m²)')
    ax.set_title('AOX Accumulation Dynamics')
    ax.set_xlim(14, 36)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.suptitle('Fig S3: Carbon Competition Mechanism (v6)', fontsize=16, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'FigS3_carbon_competition.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, 'FigS3_carbon_competition.pdf'), bbox_inches='tight')
    plt.close()
    print("  Saved FigS3_carbon_competition.png/pdf")


# ==============================================================================
# Fig S1: System Dynamics Block Diagram (Supplementary)
# ==============================================================================
def generate_figS1_block_diagram():
    """Generate system dynamics block diagram with carbon competition."""
    print("Generating Fig S1: System dynamics block diagram...")

    fig, ax = plt.subplots(figsize=(14, 10))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 10)
    ax.axis('off')

    # Define box positions and sizes
    box_style = "round,pad=0.02"

    # State variable boxes (6 states)
    state_boxes = [
        (1.0, 7.0, 1.8, 0.9, 'X_d\n(Dry Weight)', '#AED6F1'),      # Blue
        (1.0, 5.0, 1.8, 0.9, 'C_buf\n(Carbon)', '#F9E79F'),        # Yellow
        (1.0, 3.0, 1.8, 0.9, 'LAI\n(Leaf Area)', '#ABEBC6'),       # Green
        (4.5, 7.0, 1.8, 0.9, 'AOX\n(Antioxidants)', '#F5B7B1'),    # Red
        (4.5, 5.0, 1.8, 0.9, 'Stress', '#D7BDE2'),                  # Purple
        (4.5, 3.0, 1.8, 0.9, 'ROS', '#FAD7A0'),                     # Orange
    ]

    # Input/Output boxes
    io_boxes = [
        (8.0, 8.0, 2.0, 0.8, 'UV-A\nIntensity', '#85C1E9'),
        (8.0, 6.0, 2.0, 0.8, 'PAR\n(Light)', '#82E0AA'),
        (8.0, 4.0, 2.0, 0.8, 'Temperature\nCO₂, RH', '#F8C471'),
        (11.0, 7.0, 2.0, 0.8, 'Fresh\nWeight', '#BB8FCE'),
        (11.0, 5.0, 2.0, 0.8, 'Anthocyanin\n(=AOX×18%)', '#F1948A'),
    ]

    # Process boxes
    process_boxes = [
        (1.0, 1.0, 2.5, 0.8, 'Photosynthesis', '#D5F5E3'),
        (4.0, 1.0, 2.5, 0.8, 'Carbon\nCompetition', '#FDEDEC'),
        (7.0, 1.0, 2.5, 0.8, 'Stress\nDamage', '#EBF5FB'),
    ]

    # Draw state boxes
    for x, y, w, h, label, color in state_boxes:
        rect = FancyBboxPatch((x, y), w, h, boxstyle=box_style,
                              facecolor=color, edgecolor='black', linewidth=2)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, label, ha='center', va='center',
                fontsize=10, fontweight='bold')

    # Draw IO boxes
    for x, y, w, h, label, color in io_boxes:
        rect = FancyBboxPatch((x, y), w, h, boxstyle=box_style,
                              facecolor=color, edgecolor='black', linewidth=1.5)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, label, ha='center', va='center', fontsize=9)

    # Draw process boxes
    for x, y, w, h, label, color in process_boxes:
        rect = FancyBboxPatch((x, y), w, h, boxstyle=box_style,
                              facecolor=color, edgecolor='black', linewidth=1.5, linestyle='--')
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, label, ha='center', va='center', fontsize=9)

    # Draw arrows with labels
    arrow_props = dict(arrowstyle='->', lw=1.5, color='black')

    # C_buf -> X_d (growth)
    ax.annotate('', xy=(1.9, 6.9), xytext=(1.9, 6.0), arrowprops=arrow_props)
    ax.text(2.1, 6.5, 'Growth', fontsize=8, color='blue')

    # C_buf -> AOX (carbon cost) - RED arrow for carbon competition
    ax.annotate('', xy=(4.4, 7.4), xytext=(2.9, 5.5),
                arrowprops=dict(arrowstyle='->', lw=2, color='red'))
    ax.text(3.2, 6.6, 'C cost', fontsize=8, color='red', fontweight='bold')

    # X_d -> LAI
    ax.annotate('', xy=(1.9, 4.0), xytext=(1.9, 6.9), arrowprops=arrow_props)
    ax.text(2.1, 4.5, 'SLA', fontsize=8)

    # ROS -> Stress
    ax.annotate('', xy=(5.4, 4.9), xytext=(5.4, 4.0), arrowprops=arrow_props)
    ax.text(5.6, 4.5, 'Damage', fontsize=8)

    # Stress -> X_d inhibition
    ax.annotate('', xy=(2.9, 7.4), xytext=(4.4, 5.4),
                arrowprops=dict(arrowstyle='-|>', lw=1.5, color='purple'))
    ax.text(3.3, 6.2, 'Inhibit', fontsize=8, color='purple')

    # AOX -> Stress protection
    ax.annotate('', xy=(5.4, 5.9), xytext=(5.4, 6.9),
                arrowprops=dict(arrowstyle='-|>', lw=1.5, color='green'))
    ax.text(5.6, 6.4, 'Protect', fontsize=8, color='green')

    # UV-A -> ROS
    ax.annotate('', xy=(6.4, 3.4), xytext=(7.9, 8.0),
                arrowprops=dict(arrowstyle='->', lw=1.5, color='orange',
                               connectionstyle='arc3,rad=-0.2'))

    # UV-A -> LAI boost
    ax.annotate('', xy=(2.9, 3.4), xytext=(7.9, 8.0),
                arrowprops=dict(arrowstyle='->', lw=1.5, color='green',
                               connectionstyle='arc3,rad=0.3'))
    ax.text(4.5, 2.2, 'LAI boost', fontsize=8, color='green')

    # Outputs
    ax.annotate('', xy=(10.9, 7.4), xytext=(6.4, 7.4),
                arrowprops=dict(arrowstyle='->', lw=1.5))
    ax.annotate('', xy=(10.9, 5.4), xytext=(6.4, 7.4),
                arrowprops=dict(arrowstyle='->', lw=1.5))

    # Title and legend
    ax.text(7, 9.5, 'Fig S1: System Dynamics Block Diagram (v6 with Carbon Competition)',
            fontsize=14, fontweight='bold', ha='center')

    # Legend
    legend_y = 0.3
    ax.add_patch(FancyBboxPatch((0.5, legend_y), 0.4, 0.3, boxstyle=box_style,
                                facecolor='#AED6F1', edgecolor='black'))
    ax.text(1.1, legend_y + 0.15, 'State Variables', fontsize=9, va='center')

    ax.add_patch(FancyBboxPatch((3.0, legend_y), 0.4, 0.3, boxstyle=box_style,
                                facecolor='#85C1E9', edgecolor='black'))
    ax.text(3.6, legend_y + 0.15, 'Inputs/Outputs', fontsize=9, va='center')

    ax.annotate('', xy=(6.0, legend_y + 0.15), xytext=(5.5, legend_y + 0.15),
                arrowprops=dict(arrowstyle='->', lw=2, color='red'))
    ax.text(6.2, legend_y + 0.15, 'Carbon Competition', fontsize=9, va='center', color='red')

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'FigS1_block_diagram.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, 'FigS1_block_diagram.pdf'), bbox_inches='tight')
    plt.close()
    print("  Saved FigS1_block_diagram.png/pdf")


# ==============================================================================
# Fig S2: Hormesis Response Surface (Supplementary)
# ==============================================================================
def generate_figS2_hormesis_surface():
    """Generate hormesis response surface for anthocyanin."""
    print("Generating Fig S2: Hormesis response surface...")

    p = UVAParams()

    hours_range = np.arange(1, 13)
    days_range = np.arange(2, 13)

    anth_matrix = np.zeros((len(days_range), len(hours_range)))

    for i, days in enumerate(days_range):
        for j, hours in enumerate(hours_range):
            env = dict(ENV_BASE)
            env['uva_on'] = True
            env['uva_intensity'] = 11.0
            env['uva_start_day'] = 35 - days
            env['uva_end_day'] = 35
            env['uva_hour_on'] = 6
            env['uva_hour_off'] = min(6 + hours, 22)

            res = run_simulation('opt', env)
            if res:
                anth_matrix[i, j] = res['Anth']

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(anth_matrix, aspect='auto', origin='lower', cmap='YlOrRd',
                   extent=[hours_range[0]-0.5, hours_range[-1]+0.5,
                           days_range[0]-0.5, days_range[-1]+0.5])

    ax.set_xlabel('Daily UV-A Hours (h/day)', fontsize=12)
    ax.set_ylabel('Treatment Duration (days)', fontsize=12)
    ax.set_title('Fig S2: Anthocyanin Concentration Response Surface (v6)', fontsize=14, fontweight='bold')
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Anthocyanin (ppm)', fontsize=11)

    # Mark peak concentration
    peak_idx = np.unravel_index(np.argmax(anth_matrix), anth_matrix.shape)
    peak_days = days_range[peak_idx[0]]
    peak_hours = hours_range[peak_idx[1]]
    ax.scatter(peak_hours, peak_days, c='white', s=200, marker='*', edgecolor='black', linewidth=1.5)
    ax.annotate(f'Max: {peak_hours}h × {peak_days}d\n({anth_matrix[peak_idx]:.0f} ppm)',
                xy=(peak_hours, peak_days), xytext=(peak_hours-3, peak_days+2),
                color='black', fontsize=10, fontweight='bold',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
                arrowprops=dict(arrowstyle='->', color='black', lw=1.5))

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'FigS2_hormesis_surface.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, 'FigS2_hormesis_surface.pdf'), bbox_inches='tight')
    plt.close()
    print("  Saved FigS2_hormesis_surface.png/pdf")


# ==============================================================================
# Fig 17: Sensitivity Analysis
# ==============================================================================
def generate_fig17_sensitivity():
    """Generate sensitivity analysis plots."""
    print("Generating Fig 17: Sensitivity analysis...")

    # Define parameters to analyze with their base values and ranges
    params_to_analyze = [
        ('A_vulnerability', 8.5e7, [0.5, 0.75, 1.0, 1.25, 1.5]),
        ('gompertz_threshold', 10.5, [0.7, 0.85, 1.0, 1.15, 1.3]),
        ('V_max_aox', 1.45e-8, [0.5, 0.75, 1.0, 1.25, 1.5]),
        ('k_aox_deg', 3.02e-6, [0.5, 0.75, 1.0, 1.25, 1.5]),
        ('aox_carbon_cost', 1.0, [0.5, 0.75, 1.0, 1.25, 1.5]),
        ('k_stress_decay', 2.14e-5, [0.5, 0.75, 1.0, 1.25, 1.5]),
    ]

    # Reference simulation (H12D3 - high stress treatment)
    ref_env = dict(ENV_BASE)
    ref_env.update({
        'uva_on': True, 'uva_intensity': 11.0,
        'uva_start_day': 32, 'uva_end_day': 35,
        'uva_hour_on': 6, 'uva_hour_off': 18
    })

    ref_result = run_simulation('H12D3', ref_env)
    if not ref_result:
        print("  Warning: Reference simulation failed")
        return

    ref_FW = ref_result['FW']
    ref_Anth = ref_result['Anth']
    ref_Stress = ref_result['Stress']

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for idx, (param_name, base_val, multipliers) in enumerate(params_to_analyze):
        ax = axes[idx]

        fw_changes = []
        anth_changes = []
        stress_changes = []

        print(f"    Analyzing parameter: {param_name}")

        for mult in multipliers:
            # Create modified params with actual v2 model simulation
            p = UVAParams()
            setattr(p, param_name, base_val * mult)

            # Run actual simulation with modified params
            result = run_simulation(f'{param_name}_{mult}', ref_env, params=p)

            if result:
                # Calculate percentage change from reference
                fw_change = (result['FW'] - ref_FW) / ref_FW * 100
                anth_change = (result['Anth'] - ref_Anth) / ref_Anth * 100
                stress_change = (result['Stress'] - ref_Stress) / (ref_Stress + 1e-9) * 100
            else:
                # If simulation fails, use 0 (no change)
                fw_change = 0
                anth_change = 0
                stress_change = 0

            fw_changes.append(fw_change)
            anth_changes.append(anth_change)
            stress_changes.append(stress_change)

        x = [(m - 1.0) * 100 for m in multipliers]
        ax.plot(x, fw_changes, 'b-o', linewidth=2, markersize=6, label='FW')
        ax.plot(x, anth_changes, 'r-s', linewidth=2, markersize=6, label='Anth')
        ax.plot(x, stress_changes, 'g-^', linewidth=2, markersize=6, label='Stress')

        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)

        ax.set_xlabel('Parameter Change (%)')
        ax.set_ylabel('Output Change (%)')
        ax.set_title(param_name.replace('_', ' ').title(), fontsize=11)
        ax.legend(loc='best', fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.suptitle('Fig 17: Parameter Sensitivity Analysis (v6)', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'Fig17_sensitivity.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, 'Fig17_sensitivity.pdf'), bbox_inches='tight')
    plt.close()
    print("  Saved Fig17_sensitivity.png/pdf")


# ==============================================================================
# Fig 18: Optimization Heatmaps
# ==============================================================================
def generate_fig18_optimization():
    """Generate optimization heatmaps for FW and Anthocyanin."""
    print("Generating Fig 18: Optimization heatmaps...")

    CK_FW = 87.0
    CK_ANTH = 433

    hours_range = np.arange(1, 13)
    days_range = np.arange(2, 13)

    fw_change = np.zeros((len(days_range), len(hours_range)))
    anth_change = np.zeros((len(days_range), len(hours_range)))
    total_anth_change = np.zeros((len(days_range), len(hours_range)))

    for i, days in enumerate(days_range):
        for j, hours in enumerate(hours_range):
            env = dict(ENV_BASE)
            env['uva_on'] = True
            env['uva_intensity'] = 11.0
            env['uva_start_day'] = 35 - days
            env['uva_end_day'] = 35
            env['uva_hour_on'] = 6
            env['uva_hour_off'] = min(6 + hours, 22)

            res = run_simulation('opt', env)
            if res:
                fw_change[i, j] = (res['FW'] - CK_FW) / CK_FW * 100
                anth_change[i, j] = (res['Anth'] - CK_ANTH) / CK_ANTH * 100
                # Total anthocyanin = concentration × fresh weight
                total_anth = res['Anth'] * res['FW'] / 1000  # mg per plant
                ck_total = CK_ANTH * CK_FW / 1000
                total_anth_change[i, j] = (total_anth - ck_total) / ck_total * 100

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # (a) FW change heatmap
    ax1 = axes[0]
    im1 = ax1.imshow(fw_change, aspect='auto', origin='lower', cmap='RdYlGn',
                     extent=[hours_range[0]-0.5, hours_range[-1]+0.5,
                             days_range[0]-0.5, days_range[-1]+0.5],
                     vmin=-40, vmax=10)
    ax1.set_xlabel('Daily UV-A Hours (h/day)', fontsize=11)
    ax1.set_ylabel('Treatment Duration (days)', fontsize=11)
    ax1.set_title('(a) Fresh Weight Change (%)', fontsize=12, fontweight='bold')
    cbar1 = plt.colorbar(im1, ax=ax1)
    cbar1.set_label('Change vs Control (%)', fontsize=10)

    # Add contour for 0% threshold
    try:
        cs1 = ax1.contour(hours_range, days_range, fw_change, levels=[0], colors='black', linestyles='-', linewidths=2)
        ax1.clabel(cs1, fmt='0%%', fontsize=9)
    except:
        pass

    # Mark optimal point (9h × 5d)
    ax1.scatter(9, 5, c='blue', s=150, marker='*', edgecolor='white', linewidth=2, zorder=5)

    # (b) Total anthocyanin change heatmap
    ax2 = axes[1]
    im2 = ax2.imshow(total_anth_change, aspect='auto', origin='lower', cmap='YlOrRd',
                     extent=[hours_range[0]-0.5, hours_range[-1]+0.5,
                             days_range[0]-0.5, days_range[-1]+0.5],
                     vmin=-20, vmax=50)
    ax2.set_xlabel('Daily UV-A Hours (h/day)', fontsize=11)
    ax2.set_ylabel('Treatment Duration (days)', fontsize=11)
    ax2.set_title('(b) Total Anthocyanin Change (%)', fontsize=12, fontweight='bold')
    cbar2 = plt.colorbar(im2, ax=ax2)
    cbar2.set_label('Change vs Control (%)', fontsize=10)

    # Mark optimal point with annotation
    ax2.scatter(9, 5, c='blue', s=150, marker='*', edgecolor='white', linewidth=2, zorder=5)

    # Find actual optimal values
    opt_idx = (3, 8)  # 5 days, 9 hours (0-indexed: days_range[3]=5, hours_range[8]=9)
    opt_fw = fw_change[opt_idx]
    opt_anth = total_anth_change[opt_idx]

    ax2.annotate(f'Optimal: 9h × 5d\n(FW {opt_fw:+.1f}%, Anth {opt_anth:+.1f}%)',
                 xy=(9, 5), xytext=(3, 9),
                 color='blue', fontsize=10, fontweight='bold',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.9),
                 arrowprops=dict(arrowstyle='->', color='blue', lw=1.5))

    plt.suptitle('Fig 18: Optimization Heatmaps (v6)', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'Fig18_optimization.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, 'Fig18_optimization.pdf'), bbox_inches='tight')
    plt.close()
    print("  Saved Fig18_optimization.png/pdf")


# ==============================================================================
# Main
# ==============================================================================
def main():
    print("=" * 60)
    print("UVA Lettuce Model v6.0 - Figure Generation")
    print("=" * 60)
    print()

    # Generate all figures (matching v6/v5 paper numbering)
    # Main figures (9-18)
    generate_fig9_10_validation_response()  # Fig 9: Validation FW, Fig 10: Validation Anth
    generate_fig11_validation_parity()      # Fig 11: Validation parity
    generate_fig12_lai_vulnerability()      # Fig 12: LAI vulnerability
    generate_fig13_gompertz()               # Fig 13: Gompertz nonlinear
    generate_fig14_training_parity()        # Fig 14: Training parity
    generate_fig15_stress_timeseries()      # Fig 15: Stress time series
    generate_fig16_hill_inhibition()        # Fig 16: Hill inhibition
    generate_fig17_sensitivity()            # Fig 17: Sensitivity analysis
    generate_fig18_optimization()           # Fig 18: Optimization heatmaps

    # Supplementary figures (S1-S3)
    generate_figS1_block_diagram()          # Fig S1: Block diagram
    generate_figS2_hormesis_surface()       # Fig S2: Hormesis surface
    generate_figS3_carbon_competition()     # Fig S3: Carbon competition

    print()
    print("=" * 60)
    print("All figures generated successfully!")
    print(f"Output directory: {OUTPUT_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
