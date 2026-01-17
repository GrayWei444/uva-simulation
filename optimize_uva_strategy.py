"""
================================================================================
UVA Irradiation Strategy Optimization Analysis
================================================================================
Objective: Find the optimal UVA irradiation time (hours per day) and duration combination
Evaluation metrics:
  1. Maximize anthocyanin content
  2. Maintain or increase fresh weight (no reduction)
  3. Composite score = anthocyanin increase ratio × (1 + fresh weight change ratio)

================================================================================
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Set font configuration
rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial Unicode MS', 'sans-serif']
rcParams['axes.unicode_minus'] = False

# Import model
from simulate_uva_model_v10 import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio, nonlinear_damage_factor
from model_config import ENV_BASE, SIMULATION

def simulate_treatment(hours_per_day, num_days, start_day=None, night_mode=False):
    """
    Simulate a specified UVA treatment protocol

    Parameters:
    - hours_per_day: Hours of irradiation per day (1-16)
    - num_days: Total days of irradiation (1-15)
    - start_day: Day to start irradiation (from seeding), default is 35 - num_days
    - night_mode: Whether to use nighttime irradiation mode

    Returns:
    - FW: Fresh weight (g)
    - Anth: Anthocyanin concentration (ug/g FW)
    - success: Whether simulation succeeded
    """
    p = UVAParams()

    # Calculate start day (ensure harvest at Day 35)
    harvest_day = 35
    if start_day is None:
        start_day = harvest_day - num_days

    end_day = start_day + num_days - 1

    # Set irradiation time period
    if night_mode:
        # Nighttime irradiation: starts at 22:00
        uva_hour_on = 22
        uva_hour_off = (22 + hours_per_day) % 24
    else:
        # Daytime irradiation: starts at 10:00
        uva_hour_on = 10
        uva_hour_off = 10 + hours_per_day

    # Build environment parameters
    env = ENV_BASE.copy()
    env['uva_on'] = True
    env['uva_start_day'] = start_day
    env['uva_end_day'] = end_day
    env['uva_hour_on'] = uva_hour_on
    env['uva_hour_off'] = uva_hour_off
    env['uva_intensity'] = 11.0  # W/m² (consistent with paper)

    # Initial conditions
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    C_buf_init = Xd_init * 0.1
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6

    initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]

    # Simulation time
    transplant_day = SIMULATION['transplant_offset']
    simulation_days = SIMULATION['days']
    t_start = transplant_day * 86400
    t_end = (transplant_day + simulation_days) * 86400

    # Solve ODE
    sol = solve_ivp(
        uva_sun_derivatives,
        (t_start, t_end),
        initial_state,
        args=(p, env),
        method='RK45',
        max_step=300
    )

    if not sol.success:
        return None, None, False

    Xd_f, Cbuf_f, LAI_f, Anth_f, Stress_f, ROS_f = sol.y[:, -1]

    # Calculate average Stress (consistent with main model - only during UVA irradiation period)
    uva_start = start_day * 86400
    stress_sum = 0
    stress_count = 0
    for i in range(len(sol.t)):
        if sol.t[i] >= uva_start:
            stress_sum += sol.y[4, i]
            stress_count += 1
    avg_stress = stress_sum / max(1, stress_count)

    # Calculate nonlinear_factor
    nonlin_factor = nonlinear_damage_factor(hours_per_day, p)

    # Use avg_stress and nonlin_factor to calculate dw_fw_ratio (consistent with main model)
    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(avg_stress, p, nonlin_factor)

    FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
    FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
    Anth_sim = Anth_f / FW_total_kg * 1e6

    return FW_sim, Anth_sim, True


def run_optimization():
    """Execute optimization search"""

    print("=" * 80)
    print("UVA Irradiation Strategy Optimization Analysis")
    print("=" * 80)

    # First simulate control group
    print("\n[1] Simulating control group (no UVA)...")

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
        max_step=300
    )

    Xd_f, Cbuf_f, LAI_f, Anth_f, Stress_f, ROS_f = sol.y[:, -1]
    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress_f, p)
    FW_ck = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
    FW_total_kg = FW_ck / 1000 * ENV_BASE['plant_density']
    Anth_ck = Anth_f / FW_total_kg * 1e6

    print(f"   Control: FW = {FW_ck:.1f} g, Anth = {Anth_ck:.1f} ug/g FW")

    # Search parameter range
    hours_range = range(1, 13)  # 1-12 hours/day
    days_range = range(2, 13)   # 2-12 days

    print(f"\n[2] Search range: {min(hours_range)}-{max(hours_range)} hours/day, {min(days_range)}-{max(days_range)} days")
    print("    (Daytime irradiation mode, start day calculated backwards from Day 35)")

    # Store results
    results = []

    print("\n[3] Executing optimization search...")
    total = len(hours_range) * len(days_range)
    count = 0

    for hours in hours_range:
        for days in days_range:
            count += 1
            if count % 20 == 0:
                print(f"    Progress: {count}/{total} ({100*count/total:.0f}%)")

            FW, Anth, success = simulate_treatment(hours, days)

            if success:
                # Calculate total anthocyanin (ug/plant)
                anth_total = Anth * FW  # concentration × fresh weight
                anth_total_ck = Anth_ck * FW_ck

                # Calculate relative changes
                fw_change = (FW - FW_ck) / FW_ck * 100
                anth_change = (Anth - Anth_ck) / Anth_ck * 100  # concentration change
                anth_total_change = (anth_total - anth_total_ck) / anth_total_ck * 100  # total change

                # New score: total anthocyanin change (considers both concentration and fresh weight)
                score = anth_total_change

                # Safe score: total anthocyanin change when fresh weight doesn't decrease more than 5%
                score_safe = anth_total_change if fw_change >= -5 else -999

                results.append({
                    'hours': hours,
                    'days': days,
                    'FW': FW,
                    'Anth': Anth,
                    'anth_total': anth_total,
                    'fw_change': fw_change,
                    'anth_change': anth_change,
                    'anth_total_change': anth_total_change,
                    'score': score,
                    'score_safe': score_safe,
                    'total_hours': hours * days
                })

    print(f"\n[4] Search completed, {len(results)} valid combinations found")

    # Sort results
    results_sorted_score = sorted(results, key=lambda x: x['score'], reverse=True)
    results_sorted_safe = sorted(results, key=lambda x: x['score_safe'], reverse=True)
    results_sorted_anth = sorted(results, key=lambda x: x['anth_change'], reverse=True)

    # Display results
    print("\n" + "=" * 80)
    print("Optimization Results")
    print("=" * 80)

    print("\n[Strategy 1] Maximum Total Anthocyanin (Concentration x Fresh Weight):")
    print("-" * 100)
    print(f"{'Rank':<4} {'Hours/Day':<10} {'Days':<6} {'Total Hrs':<10} {'FW(g)':<10} {'Anth Conc':<12} {'Anth Total':<12} {'FW Change':<12} {'Total Chg':<10}")
    print("-" * 100)
    for i, r in enumerate(results_sorted_score[:10]):
        print(f"{i+1:<4} {r['hours']:<10} {r['days']:<6} {r['total_hours']:<10} "
              f"{r['FW']:<10.1f} {r['Anth']:<12.1f} {r['anth_total']:<12.0f} {r['fw_change']:>+10.1f}% {r['anth_total_change']:>+8.1f}%")

    print("\n[Strategy 2] Maximum Total Anthocyanin with FW >= -5%:")
    print("-" * 100)
    print(f"{'Rank':<4} {'Hours/Day':<10} {'Days':<6} {'Total Hrs':<10} {'FW(g)':<10} {'Anth Conc':<12} {'Anth Total':<12} {'FW Change':<12} {'Total Chg':<10}")
    print("-" * 100)
    safe_results = [r for r in results_sorted_safe if r['score_safe'] > -999]
    for i, r in enumerate(safe_results[:10]):
        print(f"{i+1:<4} {r['hours']:<10} {r['days']:<6} {r['total_hours']:<10} "
              f"{r['FW']:<10.1f} {r['Anth']:<12.1f} {r['anth_total']:<12.0f} {r['fw_change']:>+10.1f}% {r['anth_total_change']:>+8.1f}%")

    print("\n[Strategy 3] Maximum Anthocyanin Concentration (ignoring fresh weight):")
    print("-" * 100)
    print(f"{'Rank':<4} {'Hours/Day':<10} {'Days':<6} {'Total Hrs':<10} {'FW(g)':<10} {'Anth Conc':<12} {'Anth Total':<12} {'FW Change':<12} {'Conc Chg':<10}")
    print("-" * 100)
    for i, r in enumerate(results_sorted_anth[:10]):
        print(f"{i+1:<4} {r['hours']:<10} {r['days']:<6} {r['total_hours']:<10} "
              f"{r['FW']:<10.1f} {r['Anth']:<12.1f} {r['anth_total']:<12.0f} {r['fw_change']:>+10.1f}% {r['anth_change']:>+8.1f}%")

    # Create heatmap
    print("\n[5] Generating heatmap...")

    # Prepare heatmap data
    hours_list = sorted(set(r['hours'] for r in results))
    days_list = sorted(set(r['days'] for r in results))

    fw_matrix = np.zeros((len(days_list), len(hours_list)))
    anth_total_matrix = np.zeros((len(days_list), len(hours_list)))

    for r in results:
        i = days_list.index(r['days'])
        j = hours_list.index(r['hours'])
        fw_matrix[i, j] = r['fw_change']
        anth_total_matrix[i, j] = r['anth_total_change']

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Fresh weight change heatmap
    im1 = axes[0].imshow(fw_matrix, cmap='RdYlGn', aspect='auto',
                         extent=[min(hours_list)-0.5, max(hours_list)+0.5,
                                max(days_list)+0.5, min(days_list)-0.5])
    axes[0].set_xlabel('Hours per Day')
    axes[0].set_ylabel('Number of Days')
    axes[0].set_title('Fresh Weight Change (%)')
    axes[0].set_xticks(hours_list)
    axes[0].set_yticks(days_list)
    plt.colorbar(im1, ax=axes[0])

    # Total anthocyanin change heatmap
    im2 = axes[1].imshow(anth_total_matrix, cmap='RdYlGn', aspect='auto',
                         extent=[min(hours_list)-0.5, max(hours_list)+0.5,
                                max(days_list)+0.5, min(days_list)-0.5])
    axes[1].set_xlabel('Hours per Day')
    axes[1].set_ylabel('Number of Days')
    axes[1].set_title('Total Anthocyanin Change (%)\n(Concentration x Fresh Weight)')
    axes[1].set_xticks(hours_list)
    axes[1].set_yticks(days_list)
    plt.colorbar(im2, ax=axes[1])

    # Mark best point (highest total anthocyanin)
    best = results_sorted_score[0]
    for ax in axes:
        ax.plot(best['hours'], best['days'], 'w*', markersize=15, markeredgecolor='black')

    # Mark safe best point (highest total anthocyanin with FW >= -5%)
    best_safe = safe_results[0] if safe_results else None
    if best_safe:
        for ax in axes:
            ax.plot(best_safe['hours'], best_safe['days'], 'y*', markersize=15, markeredgecolor='black')

    plt.tight_layout()
    plt.savefig('optimization_heatmap.png', dpi=150, bbox_inches='tight')
    print("    Heatmap saved: optimization_heatmap.png")

    # Final recommendations
    print("\n" + "=" * 80)
    print("Final Recommendations")
    print("=" * 80)
    print(f"\nControl (CK): FW = {FW_ck:.1f} g, Anth = {Anth_ck:.1f} ug/g, Total = {Anth_ck * FW_ck:.0f} ug/plant")

    best_overall = results_sorted_score[0]
    best_safe = safe_results[0] if safe_results else None

    print(f"\n* Best Total Anthocyanin Strategy (white star):")
    print(f"  - Daily irradiation: {best_overall['hours']} hours")
    print(f"  - Duration: {best_overall['days']} days")
    print(f"  - Start day: Day {35 - best_overall['days']}")
    print(f"  - Expected fresh weight: {best_overall['FW']:.1f} g ({best_overall['fw_change']:+.1f}%)")
    print(f"  - Expected anthocyanin concentration: {best_overall['Anth']:.1f} ug/g ({best_overall['anth_change']:+.1f}%)")
    print(f"  - Expected total anthocyanin: {best_overall['anth_total']:.0f} ug/plant ({best_overall['anth_total_change']:+.1f}%)")

    if best_safe:
        print(f"\n** Best Safe Strategy (yellow star, FW reduction <= 5%):")
        print(f"  - Daily irradiation: {best_safe['hours']} hours")
        print(f"  - Duration: {best_safe['days']} days")
        print(f"  - Start day: Day {35 - best_safe['days']}")
        print(f"  - Expected fresh weight: {best_safe['FW']:.1f} g ({best_safe['fw_change']:+.1f}%)")
        print(f"  - Expected anthocyanin concentration: {best_safe['Anth']:.1f} ug/g ({best_safe['anth_change']:+.1f}%)")
        print(f"  - Expected total anthocyanin: {best_safe['anth_total']:.0f} ug/plant ({best_safe['anth_total_change']:+.1f}%)")

    return results


if __name__ == "__main__":
    results = run_optimization()
