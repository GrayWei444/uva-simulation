#!/usr/bin/env python3
"""
Generate Paper Figures for UVA-Lettuce Model
=============================================
Figures matching the manuscript requirements.

Output directory: paper_figures/ (not tracked by git)

Figures (matching paper numbering):
- Fig 9: LAI vulnerability function
- Fig 10: Nonlinear damage amplification (Gompertz)
- Fig 11: Training set parity plots (FW and Anth)
- Fig 12: Stress time series
- Fig 13: Validation FW response curve
- Fig 14: Validation Anth response curve (hormesis)
- Fig 15: Validation parity plots
- Fig 16: Hill-type inhibition function
- Fig 18: Hormesis response surface
- Fig 19: Sensitivity analysis
- Fig 20: Optimization heatmaps

Usage:
    python generate_paper_figures.py
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.integrate import solve_ivp

# Import model components
from simulate_uva_model_v10 import (
    UVAParams, uva_sun_derivatives, nonlinear_damage_factor,
    calculate_dynamic_dw_fw_ratio
)

# Create output directory
OUTPUT_DIR = 'paper_figures'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Set matplotlib style for publication
plt.rcParams.update({
    'font.size': 10,
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Liberation Sans', 'Arial', 'Helvetica'],
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 1.0,
})

# ==============================================================================
# Data Configuration
# ==============================================================================

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

# Training set (Table 5)
TRAINING_DATA = {
    'CK':      {'FW_obs': 87.0, 'Anth_obs': 433},
    'L6D6':    {'FW_obs': 91.4, 'Anth_obs': 494},
    'L6D6-N':  {'FW_obs': 80.8, 'Anth_obs': 493},
    'VL3D12':  {'FW_obs': 67.0, 'Anth_obs': 482},
    'L6D12':   {'FW_obs': 60.4, 'Anth_obs': 518},
    'H12D3':   {'FW_obs': 60.6, 'Anth_obs': 651},
}

# Validation set (Table 5a)
VALIDATION_DATA = {
    'CK':      {'FW_obs': 85.2, 'Anth_obs': 413, 'hours': 0},
    'VL3D3':   {'FW_obs': 89.0, 'Anth_obs': 437, 'hours': 3},
    'L6D3':    {'FW_obs': 92.2, 'Anth_obs': 468, 'hours': 6},
    'M9D3':    {'FW_obs': 83.8, 'Anth_obs': 539, 'hours': 9},
    'H12D3':   {'FW_obs': 62.2, 'Anth_obs': 657, 'hours': 12},
    'VH15D3':  {'FW_obs': 51.3, 'Anth_obs': 578, 'hours': 15},
}

TREATMENT_CONFIGS = {
    'CK':      {'uva_on': False},
    'L6D6':    {'uva_on': True, 'uva_intensity': 11.0, 'uva_start_day': 29, 'uva_end_day': 35, 'uva_hour_on': 10, 'uva_hour_off': 16},
    'L6D6-N':  {'uva_on': True, 'uva_intensity': 11.0, 'uva_start_day': 29, 'uva_end_day': 35, 'uva_hour_on': 22, 'uva_hour_off': 4},
    'H12D3':   {'uva_on': True, 'uva_intensity': 11.0, 'uva_start_day': 32, 'uva_end_day': 35, 'uva_hour_on': 6, 'uva_hour_off': 18},
    'VL3D12':  {'uva_on': True, 'uva_intensity': 11.0, 'uva_start_day': 23, 'uva_end_day': 35, 'uva_hour_on': 10, 'uva_hour_off': 13},
    'L6D12':   {'uva_on': True, 'uva_intensity': 11.0, 'uva_start_day': 23, 'uva_end_day': 35, 'uva_hour_on': 10, 'uva_hour_off': 16},
}


def get_env_for_treatment(treatment):
    env = ENV_BASE.copy()
    if treatment in TREATMENT_CONFIGS:
        env.update(TREATMENT_CONFIGS[treatment])
    return env


def run_simulation(treatment, p, env):
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    C_buf_init = Xd_init * 0.1
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6

    initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]

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

    return sol


def calculate_fw_anth(sol, p, env):
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

    FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
    FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
    Anth_sim = Anth_f / FW_total_kg * 1e6

    return FW_sim, Anth_sim, avg_stress


# ==============================================================================
# Figure 9: LAI Vulnerability Function
# ==============================================================================
def generate_fig9_lai_vulnerability():
    print("Generating Fig 9: LAI Vulnerability Function...")

    p = UVAParams()

    # v(LAI) = A * exp(-k * LAI) + 1
    LAI_vals = np.linspace(0, 12, 100)
    vulnerability = p.A_vulnerability * np.exp(-p.k_vulnerability * LAI_vals) + 1

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(LAI_vals, vulnerability, 'b-', linewidth=2)

    # Mark treatment onset LAI values - use legend to avoid overlap
    treatment_lai = [
        ('VL3D12, L6D12 (Day 23)', 3.5, '#2ca02c'),
        ('L6D6, L6D6-N (Day 29)', 7.0, '#1f77b4'),
        ('H12D3 (Day 32)', 8.5, '#d62728'),
    ]

    for label, lai, color in treatment_lai:
        vuln = p.A_vulnerability * np.exp(-p.k_vulnerability * lai) + 1
        ax.scatter(lai, vuln, s=120, color=color, zorder=5, edgecolor='black', linewidth=1.5, label=label)

    ax.set_xlabel('LAI at Treatment Onset (m²/m²)')
    ax.set_ylabel('Vulnerability Factor v(LAI)')
    ax.set_title('LAI Vulnerability Function: v(LAI) = A·exp(−k·LAI) + 1')
    ax.set_xlim(0, 12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc='upper right', fontsize=10, framealpha=0.9)

    formula = f'A = {p.A_vulnerability:.1e}, k = {p.k_vulnerability}'
    ax.text(0.6, 0.85, formula, transform=ax.transAxes, fontsize=10,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Fig9_LAI_vulnerability.png')
    plt.savefig(f'{OUTPUT_DIR}/Fig9_LAI_vulnerability.pdf')
    plt.close()
    print(f"  Saved: {OUTPUT_DIR}/Fig9_LAI_vulnerability.png/pdf")


# ==============================================================================
# Figure 10: Nonlinear Damage Amplification (Gompertz)
# ==============================================================================
def generate_fig10_gompertz():
    print("Generating Fig 10: Gompertz Nonlinear Damage...")

    p = UVAParams()
    hours = np.linspace(0, 18, 200)
    factors = [nonlinear_damage_factor(h, p) for h in hours]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(hours, factors, 'b-', linewidth=2)

    # Mark key points
    key_hours = [3, 6, 9, 12, 15]
    key_factors = [nonlinear_damage_factor(h, p) for h in key_hours]
    ax.scatter(key_hours, key_factors, c='red', s=80, zorder=5, edgecolor='black')

    for h, f in zip(key_hours, key_factors):
        ax.annotate(f'{h}h: {f:.1f}', xy=(h, f), xytext=(h+0.5, f+15),
                    fontsize=9, ha='left')

    # Mark threshold
    ax.axvline(x=p.gompertz_threshold, color='gray', linestyle='--', linewidth=1)
    ax.text(p.gompertz_threshold + 0.2, 10, f'Threshold = {p.gompertz_threshold}h',
            fontsize=9, color='gray')

    ax.set_xlabel('Daily UVA Exposure (hours)')
    ax.set_ylabel('Nonlinear Damage Factor')
    ax.set_title('Gompertz Nonlinear Damage Amplification')
    ax.set_xlim(0, 18)
    ax.set_ylim(0, 260)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    formula = r'$f = 1 + F_{max} \cdot \exp(-\exp(-k_{steep} \cdot (h - H_{threshold})))$'
    ax.text(0.35, 0.15, formula, transform=ax.transAxes, fontsize=10,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Fig10_Gompertz_nonlinear.png')
    plt.savefig(f'{OUTPUT_DIR}/Fig10_Gompertz_nonlinear.pdf')
    plt.close()
    print(f"  Saved: {OUTPUT_DIR}/Fig10_Gompertz_nonlinear.png/pdf")


# ==============================================================================
# Figure 11: Training Set Parity Plots
# ==============================================================================
def generate_fig11_training_parity():
    print("Generating Fig 11: Training Set Parity Plots...")

    p = UVAParams()

    fw_obs, fw_pred = [], []
    anth_obs, anth_pred = [], []
    labels = []

    for treatment in ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']:
        env = get_env_for_treatment(treatment)
        sol = run_simulation(treatment, p, env)
        if sol.success:
            fw, anth, _ = calculate_fw_anth(sol, p, env)
            fw_pred.append(fw)
            anth_pred.append(anth)
            fw_obs.append(TRAINING_DATA[treatment]['FW_obs'])
            anth_obs.append(TRAINING_DATA[treatment]['Anth_obs'])
            labels.append(treatment)

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # FW parity
    ax1 = axes[0]
    ax1.scatter(fw_obs, fw_pred, c='#4472C4', s=80, edgecolor='black', linewidth=0.5)
    for i, label in enumerate(labels):
        ax1.annotate(label, xy=(fw_obs[i], fw_pred[i]), xytext=(5, 5),
                     textcoords='offset points', fontsize=8)

    lims = [50, 100]
    ax1.plot(lims, lims, 'k-', linewidth=1, label='1:1 line')
    ax1.fill_between(lims, [l*0.95 for l in lims], [l*1.05 for l in lims],
                     alpha=0.2, color='green', label='±5% band')
    ax1.set_xlim(lims)
    ax1.set_ylim(lims)
    ax1.set_xlabel('Observed FW (g/plant)')
    ax1.set_ylabel('Predicted FW (g/plant)')
    ax1.set_title('(a) Fresh Weight')
    ax1.legend(loc='lower right')
    ax1.set_aspect('equal')

    # Anth parity
    ax2 = axes[1]
    ax2.scatter(anth_obs, anth_pred, c='#ED7D31', s=80, edgecolor='black', linewidth=0.5)
    for i, label in enumerate(labels):
        ax2.annotate(label, xy=(anth_obs[i], anth_pred[i]), xytext=(5, 5),
                     textcoords='offset points', fontsize=8)

    lims = [400, 700]
    ax2.plot(lims, lims, 'k-', linewidth=1, label='1:1 line')
    ax2.fill_between(lims, [l*0.95 for l in lims], [l*1.05 for l in lims],
                     alpha=0.2, color='green', label='±5% band')
    ax2.set_xlim(lims)
    ax2.set_ylim(lims)
    ax2.set_xlabel('Observed Anthocyanin (ppm)')
    ax2.set_ylabel('Predicted Anthocyanin (ppm)')
    ax2.set_title('(b) Anthocyanin')
    ax2.legend(loc='lower right')
    ax2.set_aspect('equal')

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Fig11_training_parity.png')
    plt.savefig(f'{OUTPUT_DIR}/Fig11_training_parity.pdf')
    plt.close()
    print(f"  Saved: {OUTPUT_DIR}/Fig11_training_parity.png/pdf")


# ==============================================================================
# Figure 12: Stress Time Series
# ==============================================================================
def generate_fig12_stress_timeseries():
    print("Generating Fig 12: Stress Time Series...")

    p = UVAParams()

    fig, ax = plt.subplots(figsize=(10, 6))

    treatments = ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']
    colors = {'CK': '#2E7D32', 'L6D6': '#1976D2', 'L6D6-N': '#7B1FA2',
              'VL3D12': '#F57C00', 'L6D12': '#C62828', 'H12D3': '#00ACC1'}
    linestyles = {'CK': '-', 'L6D6': '-', 'L6D6-N': '--',
                  'VL3D12': '-', 'L6D12': '--', 'H12D3': '-'}

    # Mark treatment onset days
    onset_days = {'VL3D12': 23, 'L6D12': 23, 'L6D6': 29, 'L6D6-N': 29, 'H12D3': 32}

    for treatment in treatments:
        env = get_env_for_treatment(treatment)
        sol = run_simulation(treatment, p, env)

        if sol.success:
            days = sol.t / 86400  # Convert to days from sowing
            ax.plot(days, sol.y[4, :], label=treatment,
                    color=colors[treatment], linewidth=1.5,
                    linestyle=linestyles[treatment])

    # Mark onset times
    for treatment, day in onset_days.items():
        ax.axvline(x=day, color=colors[treatment], linestyle=':', alpha=0.5, linewidth=1)

    ax.set_xlabel('Days After Sowing')
    ax.set_ylabel('Stress')
    ax.set_title('Stress Accumulation Time Series')
    ax.legend(loc='upper left', ncol=2)
    ax.set_xlim(14, 35)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Fig12_stress_timeseries.png')
    plt.savefig(f'{OUTPUT_DIR}/Fig12_stress_timeseries.pdf')
    plt.close()
    print(f"  Saved: {OUTPUT_DIR}/Fig12_stress_timeseries.png/pdf")


# ==============================================================================
# Figure 13: Validation FW Response Curve
# ==============================================================================
def generate_fig13_validation_fw():
    print("Generating Fig 13: Validation FW Response...")

    p = UVAParams()

    hours_list = [0, 3, 6, 9, 12, 15]
    fw_obs = [VALIDATION_DATA[t]['FW_obs'] for t in ['CK', 'VL3D3', 'L6D3', 'M9D3', 'H12D3', 'VH15D3']]
    fw_pred = []

    for hours in hours_list:
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

        sol = run_simulation('test', p, env)
        if sol.success:
            fw, _, _ = calculate_fw_anth(sol, p, env)
            fw_pred.append(fw)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(hours_list, fw_obs, 'o-', color='#4472C4', linewidth=2, markersize=10, label='Observed')
    ax.plot(hours_list, fw_pred, 's--', color='#ED7D31', linewidth=2, markersize=10, label='Predicted')

    ax.set_xlabel('Daily UVA Exposure (hours)')
    ax.set_ylabel('Fresh Weight (g/plant)')
    ax.set_title('Validation: Fresh Weight Response to Daily UV-A Dose')
    ax.legend()
    ax.set_xlim(-0.5, 16)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Fig13_validation_FW.png')
    plt.savefig(f'{OUTPUT_DIR}/Fig13_validation_FW.pdf')
    plt.close()
    print(f"  Saved: {OUTPUT_DIR}/Fig13_validation_FW.png/pdf")


# ==============================================================================
# Figure 14: Validation Anth Response Curve (Hormesis)
# ==============================================================================
def generate_fig14_validation_anth():
    print("Generating Fig 14: Validation Anth Response (Hormesis)...")

    p = UVAParams()

    hours_list = [0, 3, 6, 9, 12, 15]
    anth_obs = [VALIDATION_DATA[t]['Anth_obs'] for t in ['CK', 'VL3D3', 'L6D3', 'M9D3', 'H12D3', 'VH15D3']]
    anth_pred = []

    for hours in hours_list:
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

        sol = run_simulation('test', p, env)
        if sol.success:
            _, anth, _ = calculate_fw_anth(sol, p, env)
            anth_pred.append(anth)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(hours_list, anth_obs, 'o-', color='#4472C4', linewidth=2, markersize=10, label='Observed')
    ax.plot(hours_list, anth_pred, 's--', color='#ED7D31', linewidth=2, markersize=10, label='Predicted')

    # Mark hormesis peak
    peak_idx = anth_obs.index(max(anth_obs))
    ax.annotate('Hormesis peak', xy=(hours_list[peak_idx], anth_obs[peak_idx]),
                xytext=(hours_list[peak_idx]-2, anth_obs[peak_idx]+30),
                arrowprops=dict(arrowstyle='->', color='gray'), fontsize=9)

    ax.set_xlabel('Daily UVA Exposure (hours)')
    ax.set_ylabel('Anthocyanin (ppm)')
    ax.set_title('Validation: Anthocyanin Response Showing Hormesis')
    ax.legend()
    ax.set_xlim(-0.5, 16)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Fig14_validation_Anth.png')
    plt.savefig(f'{OUTPUT_DIR}/Fig14_validation_Anth.pdf')
    plt.close()
    print(f"  Saved: {OUTPUT_DIR}/Fig14_validation_Anth.png/pdf")


# ==============================================================================
# Figure 15: Validation Parity Plots
# ==============================================================================
def generate_fig15_validation_parity():
    print("Generating Fig 15: Validation Parity Plots...")

    p = UVAParams()

    hours_list = [0, 3, 6, 9, 12, 15]
    treatments = ['CK', 'VL3D3', 'L6D3', 'M9D3', 'H12D3', 'VH15D3']

    fw_obs = [VALIDATION_DATA[t]['FW_obs'] for t in treatments]
    anth_obs = [VALIDATION_DATA[t]['Anth_obs'] for t in treatments]
    fw_pred, anth_pred = [], []

    for hours in hours_list:
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

        sol = run_simulation('test', p, env)
        if sol.success:
            fw, anth, _ = calculate_fw_anth(sol, p, env)
            fw_pred.append(fw)
            anth_pred.append(anth)

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # FW parity
    ax1 = axes[0]
    ax1.scatter(fw_obs, fw_pred, c='#4472C4', s=80, edgecolor='black', linewidth=0.5)
    for i, label in enumerate(treatments):
        ax1.annotate(label, xy=(fw_obs[i], fw_pred[i]), xytext=(5, 5),
                     textcoords='offset points', fontsize=8)

    lims = [45, 100]
    ax1.plot(lims, lims, 'k-', linewidth=1, label='1:1 line')
    ax1.fill_between(lims, [l*0.95 for l in lims], [l*1.05 for l in lims],
                     alpha=0.15, color='green', label='±5% band')
    ax1.fill_between(lims, [l*0.90 for l in lims], [l*1.10 for l in lims],
                     alpha=0.1, color='orange', label='±10% band')
    ax1.set_xlim(lims)
    ax1.set_ylim(lims)
    ax1.set_xlabel('Observed FW (g/plant)')
    ax1.set_ylabel('Predicted FW (g/plant)')
    ax1.set_title('(a) Fresh Weight')
    ax1.legend(loc='lower right', fontsize=8)
    ax1.set_aspect('equal')

    # Anth parity
    ax2 = axes[1]
    ax2.scatter(anth_obs, anth_pred, c='#ED7D31', s=80, edgecolor='black', linewidth=0.5)
    for i, label in enumerate(treatments):
        ax2.annotate(label, xy=(anth_obs[i], anth_pred[i]), xytext=(5, 5),
                     textcoords='offset points', fontsize=8)

    lims = [380, 700]
    ax2.plot(lims, lims, 'k-', linewidth=1, label='1:1 line')
    ax2.fill_between(lims, [l*0.95 for l in lims], [l*1.05 for l in lims],
                     alpha=0.15, color='green', label='±5% band')
    ax2.fill_between(lims, [l*0.90 for l in lims], [l*1.10 for l in lims],
                     alpha=0.1, color='orange', label='±10% band')
    ax2.set_xlim(lims)
    ax2.set_ylim(lims)
    ax2.set_xlabel('Observed Anthocyanin (ppm)')
    ax2.set_ylabel('Predicted Anthocyanin (ppm)')
    ax2.set_title('(b) Anthocyanin')
    ax2.legend(loc='lower right', fontsize=8)
    ax2.set_aspect('equal')

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Fig15_validation_parity.png')
    plt.savefig(f'{OUTPUT_DIR}/Fig15_validation_parity.pdf')
    plt.close()
    print(f"  Saved: {OUTPUT_DIR}/Fig15_validation_parity.png/pdf")


# ==============================================================================
# Figure 16: Hill-Type Inhibition Function
# ==============================================================================
def generate_fig16_hill_inhibition():
    print("Generating Fig 16: Hill-Type Inhibition...")

    p = UVAParams()

    nonlin_factors = np.linspace(1, 300, 200)
    K = 800
    n = 1.5
    efficiencies = 1.0 / (1.0 + (nonlin_factors / K) ** n)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(nonlin_factors, efficiencies * 100, 'b-', linewidth=2)

    # Mark key points
    key_hours = [3, 6, 9, 12, 15]
    key_nonlin = [nonlinear_damage_factor(h, p) for h in key_hours]
    key_eff = [1.0 / (1.0 + (f / K) ** n) * 100 for f in key_nonlin]

    ax.scatter(key_nonlin, key_eff, c='red', s=80, zorder=5, edgecolor='black')

    for h, f, e in zip(key_hours, key_nonlin, key_eff):
        ax.annotate(f'{h}h\n({e:.1f}%)', xy=(f, e), xytext=(f+10, e+2),
                    fontsize=8, ha='left')

    ax.set_xlabel('Nonlinear Damage Factor')
    ax.set_ylabel('Anthocyanin Synthesis Efficiency (%)')
    ax.set_title('Hill-Type Inhibition: Efficiency vs Damage Factor')
    ax.set_xlim(0, 300)
    ax.set_ylim(50, 105)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    formula = r'$\eta = \frac{1}{1 + (f/K)^n}$, K=800, n=1.5'
    ax.text(0.55, 0.85, formula, transform=ax.transAxes, fontsize=11,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Fig16_Hill_inhibition.png')
    plt.savefig(f'{OUTPUT_DIR}/Fig16_Hill_inhibition.pdf')
    plt.close()
    print(f"  Saved: {OUTPUT_DIR}/Fig16_Hill_inhibition.png/pdf")


# ==============================================================================
# Figure 18: Hormesis Response Surface - Anthocyanin Concentration
# ==============================================================================
def generate_fig18_hormesis_surface():
    print("Generating Fig 18: Hormesis Response Surface...")

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

            sol = run_simulation('opt', p, env)
            if sol.success:
                _, anth, _ = calculate_fw_anth(sol, p, env)
                anth_matrix[i, j] = anth

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(anth_matrix, aspect='auto', origin='lower', cmap='YlOrRd',
                   extent=[hours_range[0]-0.5, hours_range[-1]+0.5,
                           days_range[0]-0.5, days_range[-1]+0.5])

    ax.set_xlabel('Daily UV-A Hours (h/day)', fontsize=11)
    ax.set_ylabel('Treatment Duration (days)', fontsize=11)
    ax.set_title('Anthocyanin Concentration Response Surface', fontsize=12, fontweight='bold')
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
    plt.savefig(f'{OUTPUT_DIR}/Fig18_hormesis_surface.png', dpi=300)
    plt.savefig(f'{OUTPUT_DIR}/Fig18_hormesis_surface.pdf')
    plt.close()
    print(f"  Saved: {OUTPUT_DIR}/Fig18_hormesis_surface.png/pdf")


# ==============================================================================
# Figure 20: Optimization Heatmaps - FW Change and Anth Concentration Change
# ==============================================================================
def generate_fig20_optimization_heatmap():
    print("Generating Fig 20: Optimization Heatmaps...")

    p = UVAParams()
    CK_FW = 87.0
    CK_ANTH = 433

    hours_range = np.arange(1, 13)
    days_range = np.arange(2, 13)

    fw_change = np.zeros((len(days_range), len(hours_range)))
    anth_change = np.zeros((len(days_range), len(hours_range)))

    for i, days in enumerate(days_range):
        for j, hours in enumerate(hours_range):
            start_day = 35 - days
            env = dict(ENV_BASE)
            env['uva_on'] = True
            env['uva_intensity'] = 11.0
            env['uva_start_day'] = start_day
            env['uva_end_day'] = 35
            env['uva_hour_on'] = 6
            env['uva_hour_off'] = min(6 + hours, 22)

            sol = run_simulation('opt', p, env)
            if sol.success:
                fw, anth, _ = calculate_fw_anth(sol, p, env)
                fw_change[i, j] = (fw - CK_FW) / CK_FW * 100
                anth_change[i, j] = (anth - CK_ANTH) / CK_ANTH * 100

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

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

    # Add contour for 0% (no yield loss) threshold
    cs1 = ax1.contour(hours_range, days_range, fw_change, levels=[0], colors='black', linestyles='-', linewidths=2)
    ax1.clabel(cs1, fmt='0%%', fontsize=9)

    # Mark optimal point on FW plot too
    ax1.scatter(9, 5, c='blue', s=150, marker='*', edgecolor='white', linewidth=2, zorder=5)

    # (b) Anthocyanin concentration change heatmap
    ax2 = axes[1]
    im2 = ax2.imshow(anth_change, aspect='auto', origin='lower', cmap='YlOrRd',
                     extent=[hours_range[0]-0.5, hours_range[-1]+0.5,
                             days_range[0]-0.5, days_range[-1]+0.5],
                     vmin=0, vmax=60)
    ax2.set_xlabel('Daily UV-A Hours (h/day)', fontsize=11)
    ax2.set_ylabel('Treatment Duration (days)', fontsize=11)
    ax2.set_title('(b) Anthocyanin Change (%)', fontsize=12, fontweight='bold')
    cbar2 = plt.colorbar(im2, ax=ax2)
    cbar2.set_label('Change vs Control (%)', fontsize=10)

    # Mark optimal (9h × 5d) with annotation
    ax2.scatter(9, 5, c='blue', s=150, marker='*', edgecolor='white', linewidth=2, zorder=5)
    ax2.annotate('Optimal: 9h × 5d\n(FW +1.2%, Anth +39.7%)',
                 xy=(9, 5), xytext=(3, 9),
                 color='blue', fontsize=10, fontweight='bold',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.9),
                 arrowprops=dict(arrowstyle='->', color='blue', lw=1.5))

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Fig20_optimization_heatmap.png', dpi=300)
    plt.savefig(f'{OUTPUT_DIR}/Fig20_optimization_heatmap.pdf')
    plt.close()
    print(f"  Saved: {OUTPUT_DIR}/Fig20_optimization_heatmap.png/pdf")


# ==============================================================================
# Main
# ==============================================================================
if __name__ == '__main__':
    print("=" * 60)
    print("Generating Paper Figures for UVA-Lettuce Model")
    print("=" * 60)
    print(f"Output directory: {OUTPUT_DIR}/")
    print()

    generate_fig9_lai_vulnerability()
    generate_fig10_gompertz()
    generate_fig11_training_parity()
    generate_fig12_stress_timeseries()
    generate_fig13_validation_fw()
    generate_fig14_validation_anth()
    generate_fig15_validation_parity()
    generate_fig16_hill_inhibition()
    generate_fig18_hormesis_surface()
    generate_fig20_optimization_heatmap()

    print()
    print("=" * 60)
    print("All figures generated successfully!")
    print(f"Output files saved to: {OUTPUT_DIR}/")
    print("=" * 60)
