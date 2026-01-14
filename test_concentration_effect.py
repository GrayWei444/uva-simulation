#!/usr/bin/env python3
"""
分析花青素濃度 = 花青素總量 / 鮮重

關鍵洞察:
- 花青素濃度 (ppm) = 花青素總量 (μg) / 鮮重 (g)
- 當 FW 下降很多時，即使總量下降，濃度也可能上升 (濃縮效應)

H12D3 的高濃度可能來自:
1. 花青素合成本身並沒有那麼高
2. 但因為 FW 大幅下降 (水分損失)，造成「濃縮」
"""

import numpy as np
from scipy.integrate import solve_ivp

from simulate_uva_model_v10 import (
    UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio,
    nonlinear_damage_factor
)
from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment


def run_simulation_detailed(treatment_name):
    """執行模擬並返回詳細資訊"""
    p = UVAParams()
    env = get_env_for_treatment(treatment_name)

    # 初始條件
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    C_buf_init = Xd_init * 0.1
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6

    initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]

    # 時間設定
    transplant_day = SIMULATION['transplant_offset']
    simulation_days = SIMULATION['days']
    t_start = transplant_day * 86400
    t_end = (transplant_day + simulation_days) * 86400

    # 求解 ODE
    sol = solve_ivp(
        uva_sun_derivatives,
        (t_start, t_end),
        initial_state,
        args=(p, env),
        method='RK45',
        max_step=300
    )

    if not sol.success:
        return None

    # 提取結果
    Xd_f, Cbuf_f, LAI_f, Anth_f, Stress_f, ROS_f = sol.y[:, -1]

    # 計算平均 Stress
    uva_start_day = env.get('uva_start_day', 29)
    uva_start = uva_start_day * 86400
    stress_sum = 0
    stress_count = 0
    for i in range(len(sol.t)):
        if sol.t[i] >= uva_start:
            stress_sum += sol.y[4, i]
            stress_count += 1
    avg_stress = stress_sum / max(1, stress_count)

    # 計算 nonlinear_factor
    uva_hour_on = env.get('uva_hour_on', 0)
    uva_hour_off = env.get('uva_hour_off', 0)
    if env.get('uva_on', False):
        daily_hours = uva_hour_off - uva_hour_on if uva_hour_on < uva_hour_off else 24 - uva_hour_on + uva_hour_off
    else:
        daily_hours = 0
    nonlin_factor = nonlinear_damage_factor(daily_hours, p)

    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(avg_stress, p, nonlin_factor)

    FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000  # g/plant
    FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']  # kg/m²
    Anth_conc = Anth_f / FW_total_kg * 1e6  # ppm = μg/g

    # 花青素總量 (μg/plant)
    Anth_total_per_plant = Anth_conc * FW_sim  # μg/plant

    return {
        'FW_g': FW_sim,
        'Anth_conc': Anth_conc,  # ppm (μg/g)
        'Anth_total': Anth_total_per_plant,  # μg/plant
        'Anth_raw': Anth_f,  # ODE 輸出的原始值 (kg/m²)
        'dw_fw_ratio': dw_fw_ratio,
        'avg_stress': avg_stress,
        'nonlin_factor': nonlin_factor,
        'daily_hours': daily_hours
    }


def analyze_concentration_effect():
    """分析濃縮效應"""

    print("=" * 130)
    print("花青素濃度分析: 濃縮效應 vs 合成效率")
    print("=" * 130)

    print("\n【公式說明】")
    print("  花青素濃度 (ppm) = 花青素總量 (μg/plant) / 鮮重 (g/plant)")
    print("  當 FW ↓↓ 而 Anth_total 不變時，濃度會 ↑ (濃縮效應)")
    print()

    # 測試組別
    test_treatments = [
        ('CK_val', 0),
        ('VL3D3', 3),
        ('L6D3', 6),
        ('M9D3', 9),
        ('H12D3_val', 12),
        ('VH15D3', 15),
    ]

    print("-" * 130)
    header = f"{'處理組':>10} | {'時數':>4} | {'觀測FW':>8} | {'模擬FW':>8} | {'FW變化':>8} | {'觀測Anth':>8} | {'模擬Anth':>8} | {'誤差':>8} | {'Anth總量':>10} | {'DW/FW':>8}"
    print(header)
    print("-" * 130)

    results = {}
    for treatment, hours in test_treatments:
        obs_fw = TARGETS[treatment]['FW']
        obs_anth = TARGETS[treatment]['Anth']

        result = run_simulation_detailed(treatment)
        results[treatment] = result

        if treatment == 'CK_val':
            fw_change = 0
        else:
            fw_change = (result['FW_g'] - results['CK_val']['FW_g']) / results['CK_val']['FW_g'] * 100

        error = (result['Anth_conc'] - obs_anth) / obs_anth * 100

        print(f"{treatment:>10} | {hours:>4}h | {obs_fw:>8.1f} | {result['FW_g']:>8.1f} | {fw_change:>+7.1f}% | {obs_anth:>8} | {result['Anth_conc']:>8.1f} | {error:>+7.1f}% | {result['Anth_total']:>10.0f} | {result['dw_fw_ratio']:>8.4f}")

    print("-" * 130)

    # 計算濃縮效應
    ck_result = results['CK_val']
    ck_fw = ck_result['FW_g']
    ck_anth_total = ck_result['Anth_total']

    print("\n【濃縮效應分解】")
    print("將花青素濃度增益分解為：合成增益 + 濃縮增益")
    print("-" * 100)
    header2 = f"{'處理組':>10} | {'Anth濃度':>10} | {'vs CK':>10} | {'Anth總量增益':>12} | {'濃縮增益':>12} | {'合成增益':>12}"
    print(header2)
    print("-" * 100)

    for treatment, hours in test_treatments:
        if treatment == 'CK_val':
            continue

        result = results[treatment]

        # 實際濃度增益
        conc_gain = result['Anth_conc'] - ck_result['Anth_conc']

        # 花青素總量增益
        total_gain_pct = (result['Anth_total'] - ck_anth_total) / ck_anth_total * 100

        # 如果 FW 維持 CK 水平，濃度會是多少？
        # 這代表純粹的「合成增益」
        anth_if_ck_fw = result['Anth_total'] / ck_fw
        synthesis_gain = anth_if_ck_fw - ck_result['Anth_conc']

        # 濃縮增益 = 總增益 - 合成增益
        concentration_gain = conc_gain - synthesis_gain

        print(f"{treatment:>10} | {result['Anth_conc']:>10.1f} | {conc_gain:>+9.1f} | {total_gain_pct:>+11.1f}% | {concentration_gain:>+11.1f} | {synthesis_gain:>+11.1f}")

    print("-" * 100)

    print("\n【關鍵發現】")
    print("1. 合成增益: 實際花青素產量的增加 (μg/plant)")
    print("2. 濃縮增益: 因 FW 下降造成的「假性」濃度增加")
    print("3. H12D3/VH15D3 的高濃度很大程度來自濃縮效應 (FW 下降 27-44%)")
    print()
    print("【生物學意義】")
    print("長時間 UVA 照射 → 水分損失 → FW 下降 → 花青素「濃縮」")
    print("這解釋了為什麼 H12D3 濃度 > M9D3，即使合成效率可能已經下降")


if __name__ == "__main__":
    analyze_concentration_effect()
