#!/usr/bin/env python3
"""
調整 Stress 抑制參數

目標:
- M9D3 (avgS=69): 需要約 23% 抑制，使 Anth 從 703 降到 539
- L6D3 (avgS=6): 幾乎無抑制
- H12D3 (avgS=405): 高抑制 (濃縮效應補償)
- VH15D3 (avgS=938): 最高抑制

使用 softplus 連續過渡:
effective_S = softplus((S - threshold) / scale) × scale
inhibition = max_inhib × (effective_S^n) / (K^n + effective_S^n)
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

from simulate_uva_model_v10 import (
    UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio,
    nonlinear_damage_factor
)
from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment


def softplus(x):
    return np.log(1.0 + np.exp(np.clip(x, -50, 50)))


def stress_inhibition_v2(stress, threshold, scale, K, n, max_inhib):
    """
    改進的 Stress 抑制函數 (使用 softplus 連續過渡)
    """
    effective_stress = softplus((stress - threshold) / scale) * scale
    if effective_stress < 1e-9:
        return 0.0
    inhibition = max_inhib * (effective_stress ** n) / (K ** n + effective_stress ** n)
    return inhibition


def run_simulation_with_stress_inhib(treatment_name, threshold, scale, K_inhib, n_inhib, max_inhib):
    """執行模擬，使用修改的 Stress 抑制參數"""
    p = UVAParams()
    env = get_env_for_treatment(treatment_name)

    # 計算每日照射時數
    uva_hour_on = env.get('uva_hour_on', 0)
    uva_hour_off = env.get('uva_hour_off', 0)
    if env.get('uva_on', False):
        daily_hours = uva_hour_off - uva_hour_on if uva_hour_on < uva_hour_off else 24 - uva_hour_on + uva_hour_off
    else:
        daily_hours = 0

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

    # Monkey patch: 修改 stress 抑制和非線性抑制
    import simulate_uva_model_v10 as model

    # 保存原函數
    original_nonlin_func = model.calculate_nonlin_anth_efficiency

    # 禁用非對稱高斯 (讓 Stress 抑制來做)
    model.calculate_nonlin_anth_efficiency = lambda nf, p: 1.0

    # 修改 Stress 抑制參數
    original_K = p.K_stress_inhib
    original_n = p.n_stress_inhib
    original_max = p.max_stress_inhib

    # 注意: 我們不能直接修改 ODE 內部的公式
    # 所以這種方法不會生效...我們需要另一種方式

    # 求解 ODE
    sol = solve_ivp(
        uva_sun_derivatives,
        (t_start, t_end),
        initial_state,
        args=(p, env),
        method='RK45',
        max_step=300
    )

    # 恢復原函數
    model.calculate_nonlin_anth_efficiency = original_nonlin_func

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

    nonlin_factor = nonlinear_damage_factor(daily_hours, p)
    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(avg_stress, p, nonlin_factor)

    FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
    FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
    Anth_sim = Anth_f / FW_total_kg * 1e6

    return {
        'FW_g': FW_sim,
        'Anth': Anth_sim,
        'avg_stress': avg_stress,
        'daily_hours': daily_hours
    }


def analyze_required_inhibition():
    """分析需要多少 Stress 抑制才能達到觀測值"""
    print("=" * 100)
    print("分析所需的 Stress 抑制程度")
    print("=" * 100)

    # 目前模型 (無非對稱高斯) 的預測
    from test_remove_asymmetric import run_simulation

    treatments = [
        ('VL3D3', 3, 3.03),
        ('L6D3', 6, 6.06),
        ('M9D3', 9, 68.9),
        ('H12D3_val', 12, 405.5),
        ('VH15D3', 15, 937.6),
    ]

    print("\n【所需抑制程度計算】")
    print("-" * 100)
    print(f"{'處理組':>10} | {'時數':>4} | {'avgStress':>10} | {'觀測Anth':>10} | {'無抑制預測':>10} | {'所需抑制':>10}")
    print("-" * 100)

    for treatment, hours, avg_stress in treatments:
        obs_anth = TARGETS[treatment]['Anth']

        # 執行無非對稱高斯的模擬
        result = run_simulation(treatment, disable_asymmetric=True)
        pred_anth = result['Anth']

        # 計算所需抑制
        # 如果 pred > obs，需要抑制 (obs - pred) / obs
        # 簡化計算: 假設花青素合成與效率成正比
        if pred_anth > obs_anth:
            required_eff = obs_anth / pred_anth
            required_inhib = 1 - required_eff
        else:
            required_inhib = 0

        print(f"{treatment:>10} | {hours:>4}h | {avg_stress:>10.1f} | {obs_anth:>10} | {pred_anth:>10.1f} | {required_inhib*100:>9.1f}%")

    print("-" * 100)

    print("\n【關鍵發現】")
    print("M9D3 需要約 23% 抑制 (avgStress=69)")
    print("這意味著需要調整 Stress 抑制公式，使其在 S=69 時產生足夠抑制")
    print()
    print("【建議的參數調整】")

    # 計算不同參數組合
    print("\n嘗試不同的 K_stress_inhib 值:")
    for K in [20, 30, 40, 50, 60, 70]:
        n = 2.0
        max_inhib = 0.80
        stress = 68.9
        inhibition = max_inhib * (stress ** n) / (K ** n + stress ** n)
        print(f"  K={K}: M9D3 抑制 = {inhibition*100:.1f}%")


if __name__ == "__main__":
    analyze_required_inhibition()
