"""
Paper Figure Generation Script (v5.7)
======================================
Generates publication-ready figures for the UVA-lettuce model paper.
All labels in English for international publication.

Figures:
- Fig 14: LAI Vulnerability Function
- Fig 15: Intraday Energy Factor (energy-based formula)
- Fig 16: Model Prediction vs Observation (1:1 plot)
- Fig 17: Stress Time Series for all treatments
- Fig 18: Stress-Growth Inhibition Curve
- Supplementary: Bar chart comparison

Author: Gray (with Claude assistance)
Version: v5.7
Date: 2025-12-16
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# Create output directory
OUTPUT_DIR = "paper_figures"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =============================================================================
# Publication-quality settings
# =============================================================================
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 1.0,
})

# =============================================================================
# Color scheme for treatments (consistent across all figures)
# =============================================================================
COLORS = {
    'CK': '#2E7D32',       # Dark green
    'L6D6': '#1976D2',     # Blue
    'L6D6-N': '#7B1FA2',   # Purple
    'H12D3': '#D32F2F',    # Red
    'VL3D12': '#F57C00',   # Orange
    'L6D12': '#795548',    # Brown
}

MARKERS = {
    'CK': 'o',
    'L6D6': 's',
    'L6D6-N': '^',
    'H12D3': 'D',
    'VL3D12': 'v',
    'L6D12': 'p',
}

# =============================================================================
# Model v5.7 calibrated parameters
# =============================================================================
PARAMS = {
    # LAI vulnerability
    'LAI_ref_vuln': 7.5,
    'n_vuln': 7,
    'cap_vuln': 100,

    # Intraday factor (v5.7 - energy-based)
    'k_intraday': 1.5,
    'm_intraday': 2.0,
    'sharpness_intraday': 3.0,
    'E_50': 237.6,       # kJ/m² (≈6h @ 11 W/m²)
    'E_scale': 39.6,     # kJ/m² (≈1h @ 11 W/m²)
    'I_UVA_ref': 11.0,   # W/m² (reference UVA intensity)

    # Stress-growth inhibition
    'stress_photosynthesis_inhibition': 0.70,
    'K_stress': 5.0,

    # Circadian
    'circadian_disruption_factor': 2.0,

    # PAR conversion
    'par_conversion_factor': 1.0,  # v6.0: 移除放大效應
}

# =============================================================================
# Experimental data (v5.6 calibration results)
# =============================================================================
DATA = {
    'CK':     {'FW_obs': 87.0, 'FW_sim': 87.6, 'Anth_obs': 43.3, 'Anth_sim': 45.7, 'Stress': 0.0},
    'L6D6':   {'FW_obs': 91.4, 'FW_sim': 88.3, 'Anth_obs': 49.4, 'Anth_sim': 45.3, 'Stress': 2.0},
    'L6D6-N': {'FW_obs': 80.8, 'FW_sim': 81.4, 'Anth_obs': 49.3, 'Anth_sim': 49.2, 'Stress': 4.8},
    'H12D3':  {'FW_obs': 60.6, 'FW_sim': 62.5, 'Anth_obs': 65.1, 'Anth_sim': 68.7, 'Stress': 36.1},
    'VL3D12': {'FW_obs': 67.0, 'FW_sim': 63.6, 'Anth_obs': 48.2, 'Anth_sim': 64.6, 'Stress': 6.5},
    'L6D12':  {'FW_obs': 60.4, 'FW_sim': 61.0, 'Anth_obs': 51.8, 'Anth_sim': 72.5, 'Stress': 11.2},
}


# =============================================================================
# Figure 14: LAI Vulnerability Function
# =============================================================================
def figure_14_lai_vulnerability():
    """Generate Figure 14: LAI Vulnerability Function"""
    print("Generating Figure 14: LAI Vulnerability Function...")

    LAI_ref = PARAMS['LAI_ref_vuln']
    n = PARAMS['n_vuln']
    cap = PARAMS['cap_vuln']

    # LAI range
    LAI = np.linspace(0.5, 10, 200)

    # Vulnerability function: cap * (LAI_ref/LAI)^n / (cap + (LAI_ref/LAI)^n)
    ratio = (LAI_ref / LAI) ** n
    vulnerability = cap * ratio / (cap + ratio)

    fig, ax = plt.subplots(figsize=(6, 4.5))

    ax.plot(LAI, vulnerability, 'b-', linewidth=2, label='Vulnerability function')
    ax.axvline(LAI_ref, color='gray', linestyle='--', alpha=0.7, label=f'LAI_ref = {LAI_ref}')
    ax.axhline(50, color='gray', linestyle=':', alpha=0.5)

    # Mark key point at LAI_ref
    vuln_at_ref = cap * 1 / (cap + 1)
    ax.plot(LAI_ref, vuln_at_ref, 'ro', markersize=8, label=f'At LAI_ref: {vuln_at_ref:.1f}')

    ax.set_xlabel('Leaf Area Index (LAI, m$^2$/m$^2$)')
    ax.set_ylabel('Vulnerability Factor (dimensionless)')
    ax.set_title('LAI-Dependent Vulnerability to UVA Stress')
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 105)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # Add annotation
    ax.annotate('Young plants\n(high vulnerability)',
                xy=(2, 95), fontsize=9, ha='center',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax.annotate('Mature plants\n(low vulnerability)',
                xy=(9, 20), fontsize=9, ha='center',
                bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Fig14_LAI_vulnerability.png')
    plt.savefig(f'{OUTPUT_DIR}/Fig14_LAI_vulnerability.pdf')
    plt.close()
    print("  Saved: Fig14_LAI_vulnerability.png/pdf")


# =============================================================================
# Figure 15: Intraday Energy Factor (v5.7 - energy-based formula)
# =============================================================================
def figure_15_intraday_factor():
    """Generate Figure 15: Intraday Energy Factor (energy-based)"""
    print("Generating Figure 15: Intraday Energy Factor...")

    k = PARAMS['k_intraday']
    m = PARAMS['m_intraday']
    sharpness = PARAMS['sharpness_intraday']
    E_50 = PARAMS['E_50']
    E_scale = PARAMS['E_scale']
    I_UVA = PARAMS['I_UVA_ref']

    # Energy range (0 to 14 hours equivalent at 11 W/m²)
    E = np.linspace(0, 14 * I_UVA * 3.6, 200)  # kJ/m²
    hours = E / (I_UVA * 3.6)  # Convert back to hours for x-axis

    # Softplus function with sharpness: log(1 + exp(s*x)) / s
    def softplus(x, s=1.0):
        return np.log1p(np.exp(np.clip(s * x, -500, 500))) / s

    # Energy-based intraday factor: 1 + k * softplus((E - E_50) / E_scale, sharpness)^m
    normalized_E = (E - E_50) / E_scale
    excess = softplus(normalized_E, sharpness)
    intraday_factor = 1 + k * (excess ** m)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # Panel (a): Factor vs Hours (at reference intensity)
    ax1.plot(hours, intraday_factor, 'b-', linewidth=2)
    threshold_hours = E_50 / (I_UVA * 3.6)
    ax1.axvline(threshold_hours, color='r', linestyle='--', alpha=0.7,
               label=f'E$_{{50}}$ = {E_50:.1f} kJ/m² ({threshold_hours:.1f}h @ {I_UVA} W/m²)')
    ax1.axhline(1, color='gray', linestyle=':', alpha=0.5)

    # Safe and risk zones
    ax1.axvspan(0, threshold_hours, alpha=0.1, color='green')
    ax1.axvspan(threshold_hours, 14, alpha=0.1, color='red')

    # Mark specific points
    for h in [3, 6, 12]:
        E_h = h * I_UVA * 3.6
        norm_E = (E_h - E_50) / E_scale
        f = 1 + k * (softplus(norm_E, sharpness) ** m)
        ax1.plot(h, f, 'ko', markersize=6)
        ax1.annotate(f'{h}h: {f:.2f}', xy=(h, f), xytext=(h+0.5, f+0.3),
                   fontsize=8)

    ax1.set_xlabel(f'Daily UVA Irradiation Duration (hours @ {I_UVA} W/m²)')
    ax1.set_ylabel('Intraday Damage Factor (dimensionless)')
    ax1.set_title('(a) Factor vs Irradiation Time')
    ax1.set_xlim(0, 14)
    ax1.set_ylim(0.8, max(intraday_factor) * 1.1)
    ax1.legend(loc='upper left', fontsize=8)
    ax1.grid(True, alpha=0.3)

    # Panel (b): Factor vs Energy (showing universality)
    ax2.plot(E, intraday_factor, 'b-', linewidth=2, label=f'{I_UVA} W/m²')

    # Add curves for different intensities to show universality
    for I_alt, color, ls in [(5.5, 'green', '--'), (22, 'orange', ':')]:
        hours_alt = E / (I_alt * 3.6)
        ax2.plot(E, intraday_factor, color=color, linestyle=ls, linewidth=1.5,
                label=f'{I_alt} W/m² (same curve)', alpha=0.7)

    ax2.axvline(E_50, color='r', linestyle='--', alpha=0.7,
               label=f'E$_{{50}}$ = {E_50:.1f} kJ/m²')
    ax2.axhline(1, color='gray', linestyle=':', alpha=0.5)

    # Mark key energies
    for E_mark, label in [(118.8, '3h'), (237.6, '6h'), (475.2, '12h')]:
        idx = np.argmin(np.abs(E - E_mark))
        ax2.plot(E_mark, intraday_factor[idx], 'ko', markersize=6)
        ax2.annotate(f'{label}\n({E_mark:.0f} kJ/m²)',
                    xy=(E_mark, intraday_factor[idx]),
                    xytext=(E_mark+30, intraday_factor[idx]+0.5),
                    fontsize=7)

    ax2.set_xlabel('Cumulative UVA Energy (kJ/m²)')
    ax2.set_ylabel('Intraday Damage Factor (dimensionless)')
    ax2.set_title('(b) Factor vs Energy (Universal)')
    ax2.set_xlim(0, 600)
    ax2.set_ylim(0.8, max(intraday_factor) * 1.1)
    ax2.legend(loc='upper left', fontsize=8)
    ax2.grid(True, alpha=0.3)

    # Add annotation about universality
    ax2.annotate('Energy-based formula:\nsame curve for all intensities',
                xy=(400, 2), fontsize=8, style='italic',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Fig15_intraday_energy_factor.png')
    plt.savefig(f'{OUTPUT_DIR}/Fig15_intraday_energy_factor.pdf')
    plt.close()
    print("  Saved: Fig15_intraday_energy_factor.png/pdf")


# =============================================================================
# Figure 16: Model Prediction vs Observation (1:1 plot)
# =============================================================================
def figure_16_model_vs_observation():
    """Generate Figure 16: 1:1 comparison plot"""
    print("Generating Figure 16: Model vs Observation...")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    # Fresh Weight comparison
    for treatment, data in DATA.items():
        ax1.scatter(data['FW_obs'], data['FW_sim'],
                   c=COLORS[treatment], marker=MARKERS[treatment],
                   s=100, label=treatment, edgecolors='black', linewidth=0.5)

    # 1:1 line
    fw_range = [55, 95]
    ax1.plot(fw_range, fw_range, 'k--', linewidth=1, label='1:1 line')

    # +/-5% lines
    ax1.fill_between(fw_range, [x*0.95 for x in fw_range], [x*1.05 for x in fw_range],
                     alpha=0.2, color='gray', label='5% range')

    ax1.set_xlabel('Observed Fresh Weight (g/plant)')
    ax1.set_ylabel('Simulated Fresh Weight (g/plant)')
    ax1.set_title('(a) Fresh Weight Validation')
    ax1.set_xlim(fw_range)
    ax1.set_ylim(fw_range)
    ax1.set_aspect('equal')
    ax1.legend(loc='lower right', fontsize=8)
    ax1.grid(True, alpha=0.3)

    # Calculate R2 and RMSE for FW
    fw_obs = [d['FW_obs'] for d in DATA.values()]
    fw_sim = [d['FW_sim'] for d in DATA.values()]
    ss_res = sum((o - s)**2 for o, s in zip(fw_obs, fw_sim))
    ss_tot = sum((o - np.mean(fw_obs))**2 for o in fw_obs)
    r2_fw = 1 - ss_res / ss_tot
    rmse_fw = np.sqrt(ss_res / len(fw_obs))
    ax1.text(58, 90, f'R$^2$ = {r2_fw:.3f}\nRMSE = {rmse_fw:.2f} g', fontsize=9,
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Anthocyanin comparison
    for treatment, data in DATA.items():
        ax2.scatter(data['Anth_obs'], data['Anth_sim'],
                   c=COLORS[treatment], marker=MARKERS[treatment],
                   s=100, label=treatment, edgecolors='black', linewidth=0.5)

    # 1:1 line
    anth_range = [40, 80]
    ax2.plot(anth_range, anth_range, 'k--', linewidth=1, label='1:1 line')

    # +/-20% lines
    ax2.fill_between(anth_range, [x*0.80 for x in anth_range], [x*1.20 for x in anth_range],
                     alpha=0.2, color='gray', label='20% range')

    ax2.set_xlabel('Observed Anthocyanin (ppm)')
    ax2.set_ylabel('Simulated Anthocyanin (ppm)')
    ax2.set_title('(b) Anthocyanin Content Validation')
    ax2.set_xlim(anth_range)
    ax2.set_ylim(anth_range)
    ax2.set_aspect('equal')
    ax2.legend(loc='lower right', fontsize=8)
    ax2.grid(True, alpha=0.3)

    # Calculate R2 for Anthocyanin
    anth_obs = [d['Anth_obs'] for d in DATA.values()]
    anth_sim = [d['Anth_sim'] for d in DATA.values()]
    ss_res = sum((o - s)**2 for o, s in zip(anth_obs, anth_sim))
    ss_tot = sum((o - np.mean(anth_obs))**2 for o in anth_obs)
    r2_anth = 1 - ss_res / ss_tot
    rmse_anth = np.sqrt(ss_res / len(anth_obs))
    ax2.text(42, 75, f'R$^2$ = {r2_anth:.3f}\nRMSE = {rmse_anth:.2f} ppm', fontsize=9,
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Fig16_model_validation.png')
    plt.savefig(f'{OUTPUT_DIR}/Fig16_model_validation.pdf')
    plt.close()
    print("  Saved: Fig16_model_validation.png/pdf")


# =============================================================================
# Figure 17: Stress Time Series (conceptual)
# =============================================================================
def figure_17_stress_dynamics():
    """Generate Figure 17: Stress dynamics over time (conceptual)"""
    print("Generating Figure 17: Stress Dynamics...")

    # Treatment timeline information
    treatments_info = {
        'CK':     {'start': 0, 'end': 0, 'final_stress': 0.0},
        'L6D6':   {'start': 6, 'end': 11, 'final_stress': 2.0},
        'L6D6-N': {'start': 6, 'end': 11, 'final_stress': 4.8},
        'H12D3':  {'start': 9, 'end': 11, 'final_stress': 36.1},
        'VL3D12': {'start': 0, 'end': 11, 'final_stress': 6.5},
        'L6D12':  {'start': 0, 'end': 11, 'final_stress': 11.2},
    }

    fig, ax = plt.subplots(figsize=(8, 5))

    days = np.linspace(0, 12, 200)

    for treatment, info in treatments_info.items():
        if treatment == 'CK':
            stress = np.zeros_like(days)
        else:
            stress = np.zeros_like(days)
            start_day = info['start']
            final_stress = info['final_stress']
            treatment_days = info['end'] - info['start'] + 1

            for i, d in enumerate(days):
                if d >= start_day:
                    progress = min((d - start_day) / treatment_days, 1.0)
                    stress[i] = final_stress * progress ** 0.7

        ax.plot(days, stress, color=COLORS[treatment], linewidth=2,
                label=f'{treatment} (final: {info["final_stress"]:.1f})')

    ax.set_xlabel('Days After Transplanting')
    ax.set_ylabel('Accumulated Stress (dimensionless)')
    ax.set_title('Stress Accumulation Dynamics Across Treatments')
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 40)
    ax.legend(loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.3)

    # Add treatment period annotations
    ax.axvspan(0, 6, alpha=0.05, color='blue')
    ax.axvspan(6, 12, alpha=0.05, color='orange')
    ax.axvline(6, color='gray', linestyle=':', alpha=0.5)
    ax.text(3, 38, 'Early phase', ha='center', fontsize=9, style='italic')
    ax.text(9, 38, 'UVA treatment phase', ha='center', fontsize=9, style='italic')

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Fig17_stress_dynamics.png')
    plt.savefig(f'{OUTPUT_DIR}/Fig17_stress_dynamics.pdf')
    plt.close()
    print("  Saved: Fig17_stress_dynamics.png/pdf")


# =============================================================================
# Figure 18: Stress-Growth Inhibition Curve
# =============================================================================
def figure_18_stress_inhibition():
    """Generate Figure 18: Stress vs Growth Inhibition relationship"""
    print("Generating Figure 18: Stress-Growth Inhibition Curve...")

    max_inhibition = PARAMS['stress_photosynthesis_inhibition']
    K_stress = PARAMS['K_stress']

    # Stress range
    stress = np.linspace(0, 50, 200)

    # Michaelis-Menten type inhibition
    inhibition = max_inhibition * stress / (K_stress + stress)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    # Panel (a): Inhibition curve
    ax1.plot(stress, inhibition * 100, 'b-', linewidth=2)
    ax1.axhline(max_inhibition * 100, color='r', linestyle='--', alpha=0.7,
                label=f'Max inhibition = {max_inhibition*100:.0f}%')
    ax1.axvline(K_stress, color='gray', linestyle=':', alpha=0.7,
                label=f'K_stress = {K_stress}')

    # Mark treatment points
    for treatment, data in DATA.items():
        s = data['Stress']
        inh = max_inhibition * s / (K_stress + s) * 100
        ax1.scatter(s, inh, c=COLORS[treatment], marker=MARKERS[treatment],
                   s=80, edgecolors='black', linewidth=0.5, zorder=5)

    ax1.set_xlabel('Accumulated Stress (dimensionless)')
    ax1.set_ylabel('Photosynthesis Inhibition (%)')
    ax1.set_title('(a) Stress-Inhibition Function')
    ax1.set_xlim(0, 50)
    ax1.set_ylim(0, 80)
    ax1.legend(loc='lower right', fontsize=8)
    ax1.grid(True, alpha=0.3)

    # Panel (b): Stress vs Fresh Weight
    for treatment, data in DATA.items():
        ax2.scatter(data['Stress'], data['FW_obs'], c=COLORS[treatment],
                   marker=MARKERS[treatment], s=100, edgecolors='black',
                   linewidth=0.5, label=f'{treatment}')

    # Theoretical curve based on CK baseline
    ck_fw = DATA['CK']['FW_sim']
    stress_theory = np.linspace(0, 40, 100)
    fw_theory = ck_fw * (1 - max_inhibition * stress_theory / (K_stress + stress_theory))
    ax2.plot(stress_theory, fw_theory, 'k--', linewidth=1.5, alpha=0.7,
             label='Theoretical curve')

    ax2.set_xlabel('Accumulated Stress (dimensionless)')
    ax2.set_ylabel('Fresh Weight (g/plant)')
    ax2.set_title('(b) Stress vs Fresh Weight')
    ax2.set_xlim(-1, 42)
    ax2.set_ylim(50, 100)
    ax2.legend(loc='upper right', fontsize=7)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/Fig18_stress_inhibition.png')
    plt.savefig(f'{OUTPUT_DIR}/Fig18_stress_inhibition.pdf')
    plt.close()
    print("  Saved: Fig18_stress_inhibition.png/pdf")


# =============================================================================
# Supplementary Figure: Bar chart comparison
# =============================================================================
def figure_sup_bar_comparison():
    """Generate supplementary bar chart comparing all treatments"""
    print("Generating Supplementary Figure: Bar Comparison...")

    treatments = list(DATA.keys())
    x = np.arange(len(treatments))
    width = 0.35

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # Fresh Weight comparison
    fw_obs = [DATA[t]['FW_obs'] for t in treatments]
    fw_sim = [DATA[t]['FW_sim'] for t in treatments]

    ax1.bar(x - width/2, fw_obs, width, label='Observed', color='steelblue', edgecolor='black')
    ax1.bar(x + width/2, fw_sim, width, label='Simulated', color='lightsteelblue', edgecolor='black')

    ax1.set_xlabel('Treatment')
    ax1.set_ylabel('Fresh Weight (g/plant)')
    ax1.set_title('(a) Fresh Weight: Observed vs Simulated')
    ax1.set_xticks(x)
    ax1.set_xticklabels(treatments)
    ax1.legend()
    ax1.set_ylim(0, 100)
    ax1.grid(True, alpha=0.3, axis='y')

    # Add error labels
    for i, (obs, sim) in enumerate(zip(fw_obs, fw_sim)):
        err = (sim - obs) / obs * 100
        color = 'red' if abs(err) > 5 else 'green'
        ax1.annotate(f'{err:+.1f}%', xy=(x[i], max(obs, sim) + 2),
                    ha='center', fontsize=8, color=color)

    # Anthocyanin comparison
    anth_obs = [DATA[t]['Anth_obs'] for t in treatments]
    anth_sim = [DATA[t]['Anth_sim'] for t in treatments]

    ax2.bar(x - width/2, anth_obs, width, label='Observed', color='purple', edgecolor='black')
    ax2.bar(x + width/2, anth_sim, width, label='Simulated', color='plum', edgecolor='black')

    ax2.set_xlabel('Treatment')
    ax2.set_ylabel('Anthocyanin Content (ppm)')
    ax2.set_title('(b) Anthocyanin: Observed vs Simulated')
    ax2.set_xticks(x)
    ax2.set_xticklabels(treatments)
    ax2.legend()
    ax2.set_ylim(0, 85)
    ax2.grid(True, alpha=0.3, axis='y')

    # Add error labels
    for i, (obs, sim) in enumerate(zip(anth_obs, anth_sim)):
        err = (sim - obs) / obs * 100
        color = 'red' if abs(err) > 20 else 'green'
        ax2.annotate(f'{err:+.1f}%', xy=(x[i], max(obs, sim) + 2),
                    ha='center', fontsize=8, color=color)

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/FigS1_bar_comparison.png')
    plt.savefig(f'{OUTPUT_DIR}/FigS1_bar_comparison.pdf')
    plt.close()
    print("  Saved: FigS1_bar_comparison.png/pdf")


# =============================================================================
# Summary table figure
# =============================================================================
def figure_summary_table():
    """Generate a summary table as a figure"""
    print("Generating Summary Table Figure...")

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.axis('off')

    # Table data
    columns = ['Treatment', 'FW_obs (g)', 'FW_sim (g)', 'FW Err (%)',
               'Anth_obs', 'Anth_sim', 'Anth Err (%)', 'Stress']

    table_data = []
    for t, d in DATA.items():
        fw_err = (d['FW_sim'] - d['FW_obs']) / d['FW_obs'] * 100
        anth_err = (d['Anth_sim'] - d['Anth_obs']) / d['Anth_obs'] * 100
        table_data.append([
            t, f"{d['FW_obs']:.1f}", f"{d['FW_sim']:.1f}", f"{fw_err:+.1f}",
            f"{d['Anth_obs']:.1f}", f"{d['Anth_sim']:.1f}", f"{anth_err:+.1f}",
            f"{d['Stress']:.1f}"
        ])

    # Add summary row
    mean_fw_err = np.mean([abs((d['FW_sim'] - d['FW_obs']) / d['FW_obs'] * 100) for d in DATA.values()])
    mean_anth_err = np.mean([abs((d['Anth_sim'] - d['Anth_obs']) / d['Anth_obs'] * 100) for d in DATA.values()])
    table_data.append(['Mean |Err|', '-', '-', f'{mean_fw_err:.1f}', '-', '-', f'{mean_anth_err:.1f}', '-'])

    table = ax.table(cellText=table_data, colLabels=columns, loc='center',
                     cellLoc='center', colColours=['lightgray']*8)
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.2, 1.5)

    # Color code errors
    for i in range(len(table_data)):
        if i < len(DATA):
            # FW error column (index 3)
            fw_err = float(table_data[i][3])
            color = 'lightgreen' if abs(fw_err) <= 5 else 'lightyellow' if abs(fw_err) <= 10 else 'lightcoral'
            table[(i+1, 3)].set_facecolor(color)

            # Anth error column (index 6)
            anth_err = float(table_data[i][6])
            color = 'lightgreen' if abs(anth_err) <= 10 else 'lightyellow' if abs(anth_err) <= 20 else 'lightcoral'
            table[(i+1, 6)].set_facecolor(color)

    ax.set_title('Model v5.7 Calibration Results Summary', fontsize=12, fontweight='bold', pad=20)

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/FigS2_summary_table.png')
    plt.savefig(f'{OUTPUT_DIR}/FigS2_summary_table.pdf')
    plt.close()
    print("  Saved: FigS2_summary_table.png/pdf")


# =============================================================================
# Main execution
# =============================================================================
def generate_all_figures():
    """Generate all paper figures"""
    print("=" * 60)
    print("Generating Paper Figures for UVA-Lettuce Model v5.7")
    print("=" * 60)
    print(f"Output directory: {OUTPUT_DIR}/")
    print()

    # Main figures
    figure_14_lai_vulnerability()
    figure_15_intraday_factor()
    figure_16_model_vs_observation()
    figure_17_stress_dynamics()
    figure_18_stress_inhibition()

    # Supplementary figures
    figure_sup_bar_comparison()
    figure_summary_table()

    print()
    print("=" * 60)
    print("All figures generated successfully!")
    print(f"Files saved to: {OUTPUT_DIR}/")
    print("=" * 60)


if __name__ == "__main__":
    generate_all_figures()
