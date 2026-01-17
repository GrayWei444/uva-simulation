"""
Quick generation of Fig20 optimization heatmap (v10.39)
Uses fewer sampling points to speed up calculation
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial Unicode MS', 'sans-serif']
rcParams['axes.unicode_minus'] = False

from simulate_uva_model_v10 import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio, nonlinear_damage_factor
from model_config import ENV_BASE, SIMULATION

def simulate_treatment(hours_per_day, num_days):
    """Simulate specified UVA treatment protocol"""
    p = UVAParams()

    harvest_day = 35
    start_day = harvest_day - num_days
    end_day = start_day + num_days - 1

    uva_hour_on = 10
    uva_hour_off = 10 + hours_per_day

    env = ENV_BASE.copy()
    env['uva_on'] = True
    env['uva_start_day'] = start_day
    env['uva_end_day'] = end_day
    env['uva_hour_on'] = uva_hour_on
    env['uva_hour_off'] = uva_hour_off
    env['uva_intensity'] = 11.0  # W/m²

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
    t_end = (transplant_day + simulation_days) * 86400

    sol = solve_ivp(
        uva_sun_derivatives,
        (t_start, t_end),
        initial_state,
        args=(p, env),
        method='RK45',
        max_step=600  # Increase step size to speed up
    )

    if not sol.success:
        return None, None, False

    Xd_f, Cbuf_f, LAI_f, Anth_f, Stress_f, ROS_f = sol.y[:, -1]

    # Calculate average Stress
    uva_start_time = start_day * 86400
    stress_sum = 0
    stress_count = 0
    for i in range(len(sol.t)):
        if sol.t[i] >= uva_start_time:
            stress_sum += sol.y[4, i]
            stress_count += 1
    avg_stress = stress_sum / max(1, stress_count)

    nonlin_factor = nonlinear_damage_factor(hours_per_day, p)
    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(avg_stress, p, nonlin_factor)

    FW = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
    FW_total_kg = FW / 1000 * ENV_BASE['plant_density']
    Anth = Anth_f / FW_total_kg * 1e6

    return FW, Anth, True


if __name__ == "__main__":
    print("=" * 60)
    print("Generating Fig20: UVA Optimization Heatmap (v10.39)")
    print("=" * 60)

    # Control group
    print("\n[1] Simulating control group...")
    p = UVAParams()
    env_ck = ENV_BASE.copy()
    env_ck['uva_on'] = False

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
    t_end = (transplant_day + simulation_days) * 86400

    sol = solve_ivp(
        uva_sun_derivatives,
        (t_start, t_end),
        initial_state,
        args=(p, env_ck),
        method='RK45',
        max_step=600
    )

    Xd_f = sol.y[0, -1]
    Anth_f = sol.y[3, -1]
    Stress_f = sol.y[4, -1]
    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress_f, p)
    FW_ck = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
    FW_total_kg = FW_ck / 1000 * ENV_BASE['plant_density']
    Anth_ck = Anth_f / FW_total_kg * 1e6

    print(f"   Control: FW = {FW_ck:.1f} g, Anth = {Anth_ck:.1f} ug/g")

    # Search range
    hours_range = [3, 6, 9, 12]  # Key hours
    days_range = [3, 4, 5, 6, 7, 8, 9, 10, 12]  # Key days

    print(f"\n[2] Search range: {hours_range} hours/day, {days_range} days")

    results = []
    total = len(hours_range) * len(days_range)
    count = 0

    print("\n[3] Executing optimization search...")
    for hours in hours_range:
        for days in days_range:
            count += 1
            print(f"   {count}/{total}: {hours}h x {days}d", end="")

            FW, Anth, success = simulate_treatment(hours, days)

            if success:
                anth_total = Anth * FW
                anth_total_ck = Anth_ck * FW_ck

                fw_change = (FW - FW_ck) / FW_ck * 100
                anth_change = (Anth - Anth_ck) / Anth_ck * 100
                anth_total_change = (anth_total - anth_total_ck) / anth_total_ck * 100

                results.append({
                    'hours': hours,
                    'days': days,
                    'FW': FW,
                    'Anth': Anth,
                    'fw_change': fw_change,
                    'anth_change': anth_change,
                    'anth_total_change': anth_total_change,
                })
                print(f" -> FW={FW:.1f}g ({fw_change:+.1f}%), Anth={Anth:.0f} ({anth_change:+.1f}%)")
            else:
                print(" -> Failed")

    # Display best results
    print("\n" + "=" * 60)
    print("Best strategies (FW loss >= -5%, highest total anthocyanin)")
    print("=" * 60)

    safe_results = [r for r in results if r['fw_change'] >= -5]
    safe_results.sort(key=lambda x: x['anth_total_change'], reverse=True)

    print(f"\n{'Rank':<4} {'Hours/Day':<8} {'Days':<6} {'FW(g)':<10} {'Anth':<10} {'FW Change':<10} {'Anth Change':<10} {'Total Change':<10}")
    print("-" * 80)
    for i, r in enumerate(safe_results[:5]):
        print(f"{i+1:<4} {r['hours']:<8} {r['days']:<6} {r['FW']:<10.1f} {r['Anth']:<10.0f} {r['fw_change']:>+8.1f}% {r['anth_change']:>+8.1f}% {r['anth_total_change']:>+8.1f}%")

    # Generate heatmap
    print("\n[4] Generating heatmap...")

    hours_list = sorted(set(r['hours'] for r in results))
    days_list = sorted(set(r['days'] for r in results))

    fw_matrix = np.full((len(days_list), len(hours_list)), np.nan)
    anth_matrix = np.full((len(days_list), len(hours_list)), np.nan)

    for r in results:
        i = days_list.index(r['days'])
        j = hours_list.index(r['hours'])
        fw_matrix[i, j] = r['fw_change']
        anth_matrix[i, j] = r['anth_total_change']

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # FW change heatmap
    im1 = axes[0].imshow(fw_matrix, cmap='RdYlGn', aspect='auto',
                          vmin=-30, vmax=10)
    axes[0].set_xticks(range(len(hours_list)))
    axes[0].set_xticklabels(hours_list)
    axes[0].set_yticks(range(len(days_list)))
    axes[0].set_yticklabels(days_list)
    axes[0].set_xlabel('Hours per day', fontsize=12)
    axes[0].set_ylabel('Days of treatment', fontsize=12)
    axes[0].set_title('Fresh Weight Change (%)', fontsize=14)
    plt.colorbar(im1, ax=axes[0])

    # Label values
    for i in range(len(days_list)):
        for j in range(len(hours_list)):
            val = fw_matrix[i, j]
            if not np.isnan(val):
                color = 'white' if abs(val) > 15 else 'black'
                axes[0].text(j, i, f'{val:.0f}', ha='center', va='center',
                            fontsize=9, color=color)

    # Anthocyanin total change heatmap
    im2 = axes[1].imshow(anth_matrix, cmap='YlOrRd', aspect='auto',
                          vmin=0, vmax=50)
    axes[1].set_xticks(range(len(hours_list)))
    axes[1].set_xticklabels(hours_list)
    axes[1].set_yticks(range(len(days_list)))
    axes[1].set_yticklabels(days_list)
    axes[1].set_xlabel('Hours per day', fontsize=12)
    axes[1].set_ylabel('Days of treatment', fontsize=12)
    axes[1].set_title('Anthocyanin Total Change (%)', fontsize=14)
    plt.colorbar(im2, ax=axes[1])

    # Label values
    for i in range(len(days_list)):
        for j in range(len(hours_list)):
            val = anth_matrix[i, j]
            if not np.isnan(val):
                color = 'white' if val > 30 else 'black'
                axes[1].text(j, i, f'{val:.0f}', ha='center', va='center',
                            fontsize=9, color=color)

    plt.suptitle('Fig. 20. UVA Optimization Heatmap (v10.39, 11 W/m²)', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig('fig20_optimization_heatmap.png', dpi=150, bbox_inches='tight')
    print("   Generated: fig20_optimization_heatmap.png")
    plt.close()

    print("\nDone!")
