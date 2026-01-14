#!/usr/bin/env python3
"""
分析轉折點問題

目標: 找出為什麼模型在 M9D3 沒有正確轉折

觀測數據顯示:
- M9D3 (9h) 是花青素總量的峰值
- 之後 H12D3, VH15D3 應該下降

但移除非對稱高斯後，M9D3 預測偏高 30%
說明模型沒有正確捕捉 9h 的轉折機制
"""

import numpy as np
from scipy.integrate import solve_ivp

from simulate_uva_model_v10 import (
    UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio,
    nonlinear_damage_factor, calculate_nonlin_anth_efficiency
)
from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment


def run_detailed_simulation(treatment_name):
    """執行詳細模擬，追蹤所有關鍵變數"""
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

    # 計算每日照射時數
    uva_hour_on = env.get('uva_hour_on', 0)
    uva_hour_off = env.get('uva_hour_off', 0)
    if env.get('uva_on', False):
        daily_hours = uva_hour_off - uva_hour_on if uva_hour_on < uva_hour_off else 24 - uva_hour_on + uva_hour_off
    else:
        daily_hours = 0

    # 計算各種因子
    nonlin_factor = nonlinear_damage_factor(daily_hours, p)
    nonlin_eff = calculate_nonlin_anth_efficiency(nonlin_factor, p)

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

    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(avg_stress, p, nonlin_factor)

    FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
    FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
    Anth_conc = Anth_f / FW_total_kg * 1e6
    Anth_total = Anth_conc * FW_sim

    return {
        'daily_hours': daily_hours,
        'nonlin_factor': nonlin_factor,
        'nonlin_eff': nonlin_eff,
        'avg_stress': avg_stress,
        'final_stress': Stress_f,
        'dw_fw_ratio': dw_fw_ratio,
        'FW_g': FW_sim,
        'Anth_conc': Anth_conc,
        'Anth_total': Anth_total,
        'Anth_raw': Anth_f,
        'LAI': LAI_f,
        'ROS': ROS_f
    }


def analyze_transition():
    """分析轉折點"""

    print("=" * 140)
    print("轉折點分析: 為什麼 M9D3 是峰值？")
    print("=" * 140)

    treatments = [
        ('CK_val', 0),
        ('VL3D3', 3),
        ('L6D3', 6),
        ('M9D3', 9),
        ('H12D3_val', 12),
        ('VH15D3', 15),
    ]

    print("\n【觀測數據 vs 模擬數據】")
    print("-" * 140)
    header = f"{'處理組':>10} | {'時數':>4} | {'觀測FW':>7} | {'模擬FW':>7} | {'觀測Anth':>8} | {'模擬Anth':>8} | {'NonlinF':>8} | {'NonlinEff':>9} | {'AvgStress':>10} | {'DW/FW':>7} | {'Anth總量':>10}"
    print(header)
    print("-" * 140)

    results = {}
    for treatment, hours in treatments:
        obs_fw = TARGETS[treatment]['FW']
        obs_anth = TARGETS[treatment]['Anth']
        obs_total = obs_fw * obs_anth  # 觀測花青素總量

        result = run_detailed_simulation(treatment)
        results[treatment] = result

        print(f"{treatment:>10} | {hours:>4}h | {obs_fw:>7.1f} | {result['FW_g']:>7.1f} | {obs_anth:>8} | {result['Anth_conc']:>8.1f} | {result['nonlin_factor']:>8.1f} | {result['nonlin_eff']*100:>8.1f}% | {result['avg_stress']:>10.2f} | {result['dw_fw_ratio']:>7.4f} | {result['Anth_total']:>10.0f}")

    print("-" * 140)

    # 分析觀測數據的花青素總量
    print("\n【觀測數據的花青素總量分析】")
    print("-" * 80)
    print(f"{'處理組':>10} | {'時數':>4} | {'觀測FW':>8} | {'觀測Anth':>10} | {'觀測總量':>12} | {'vs CK':>10}")
    print("-" * 80)

    ck_obs_total = TARGETS['CK_val']['FW'] * TARGETS['CK_val']['Anth']
    for treatment, hours in treatments:
        obs_fw = TARGETS[treatment]['FW']
        obs_anth = TARGETS[treatment]['Anth']
        obs_total = obs_fw * obs_anth
        vs_ck = (obs_total - ck_obs_total) / ck_obs_total * 100
        print(f"{treatment:>10} | {hours:>4}h | {obs_fw:>8.1f} | {obs_anth:>10} | {obs_total:>12.0f} | {vs_ck:>+9.1f}%")

    print("-" * 80)

    print("\n【關鍵洞察】")
    print("1. 觀測數據顯示 M9D3 (9h) 花青素總量最高")
    print("2. H12D3 總量下降，但濃度因 FW 下降而「濃縮」")
    print("3. 模型需要在 9h 後開始抑制花青素合成")
    print()
    print("【問題診斷】")
    print("目前的非對稱高斯只針對 M9D3 做 50% 抑制，但這不是正確的機制")
    print("正確的機制應該是: 9h 後，花青素合成效率開始單調下降")


if __name__ == "__main__":
    analyze_transition()
