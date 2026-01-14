#!/usr/bin/env python3
"""
調整單調遞減效率函數參數

目標:
- M9D3 (9h): 轉折點，效率開始下降
- H12D3 (12h): 效率繼續下降
- VH15D3 (15h): 效率最低

機制: efficiency = 1 - max_inhib × sigmoid((hours - threshold) × steepness)

需要找到合適的 threshold, steepness, max_inhib 使得:
1. 3h, 6h: 效率 ≈ 100%
2. 9h: 效率開始下降，使 M9D3 Anth 接近觀測值 539
3. 12h, 15h: 效率繼續單調下降
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

from simulate_uva_model_v10 import (
    UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio,
    nonlinear_damage_factor
)
from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment


def monotonic_efficiency(hours, threshold, steepness, max_inhib):
    """
    單調遞減效率函數

    efficiency = 1 - max_inhib × sigmoid((hours - threshold) × steepness)

    當 hours < threshold: efficiency ≈ 1
    當 hours > threshold: efficiency 單調下降
    """
    x = (hours - threshold) * steepness
    sigmoid = 1 / (1 + np.exp(-x))
    efficiency = 1.0 - max_inhib * sigmoid
    return efficiency


def run_simulation_with_efficiency(treatment_name, eff_func):
    """執行模擬，使用自定義效率函數"""
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

    # Monkey patch 效率函數
    import simulate_uva_model_v10 as model
    original_func = model.calculate_nonlin_anth_efficiency

    # 使用自定義效率
    eff = eff_func(daily_hours)
    model.calculate_nonlin_anth_efficiency = lambda nf, p: eff

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
    model.calculate_nonlin_anth_efficiency = original_func

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
        'efficiency': eff,
        'daily_hours': daily_hours
    }


def objective(params):
    """優化目標函數"""
    threshold, steepness, max_inhib = params

    # 限制參數範圍
    if threshold < 6 or threshold > 12:
        return 1e10
    if steepness < 0.1 or steepness > 2.0:
        return 1e10
    if max_inhib < 0.1 or max_inhib > 0.8:
        return 1e10

    def eff_func(hours):
        return monotonic_efficiency(hours, threshold, steepness, max_inhib)

    # 計算各組的誤差
    treatments = ['M9D3', 'H12D3_val', 'VH15D3']
    total_error = 0

    for treatment in treatments:
        obs_anth = TARGETS[treatment]['Anth']
        result = run_simulation_with_efficiency(treatment, eff_func)
        if result:
            error = ((result['Anth'] - obs_anth) / obs_anth) ** 2
            total_error += error
        else:
            return 1e10

    return total_error


def grid_search():
    """網格搜索找最佳參數"""
    print("=" * 100)
    print("網格搜索: 單調遞減效率函數參數")
    print("=" * 100)

    best_error = float('inf')
    best_params = None

    # 參數範圍
    thresholds = np.arange(7.0, 10.0, 0.5)
    steepnesses = np.arange(0.3, 1.5, 0.2)
    max_inhibs = np.arange(0.2, 0.6, 0.1)

    print(f"\n搜索範圍:")
    print(f"  threshold: {thresholds}")
    print(f"  steepness: {steepnesses}")
    print(f"  max_inhib: {max_inhibs}")
    print()

    for threshold in thresholds:
        for steepness in steepnesses:
            for max_inhib in max_inhibs:
                def eff_func(hours):
                    return monotonic_efficiency(hours, threshold, steepness, max_inhib)

                # 計算誤差
                treatments = ['M9D3', 'H12D3_val', 'VH15D3']
                errors = []
                valid = True

                for treatment in treatments:
                    obs_anth = TARGETS[treatment]['Anth']
                    result = run_simulation_with_efficiency(treatment, eff_func)
                    if result:
                        error = (result['Anth'] - obs_anth) / obs_anth * 100
                        errors.append(error)
                    else:
                        valid = False
                        break

                if valid:
                    total_error = sum(e**2 for e in errors)
                    if total_error < best_error:
                        best_error = total_error
                        best_params = (threshold, steepness, max_inhib)
                        print(f"  新最佳: threshold={threshold:.1f}, steepness={steepness:.1f}, max_inhib={max_inhib:.1f}")
                        print(f"         誤差: M9D3={errors[0]:+.1f}%, H12D3={errors[1]:+.1f}%, VH15D3={errors[2]:+.1f}%")

    print("\n" + "=" * 100)
    print(f"最佳參數: threshold={best_params[0]:.1f}, steepness={best_params[1]:.1f}, max_inhib={best_params[2]:.1f}")
    print("=" * 100)

    return best_params


def test_best_params(threshold, steepness, max_inhib):
    """測試最佳參數"""
    print("\n【最佳參數測試結果】")
    print("-" * 100)

    def eff_func(hours):
        return monotonic_efficiency(hours, threshold, steepness, max_inhib)

    # 顯示效率曲線
    print("\n效率曲線:")
    print(f"{'時數':>6} | {'效率':>10}")
    print("-" * 20)
    for hours in [0, 3, 6, 9, 12, 15]:
        eff = eff_func(hours)
        print(f"{hours:>5}h | {eff*100:>9.1f}%")
    print("-" * 20)

    # 測試所有組別
    treatments = [
        ('CK_val', 0),
        ('VL3D3', 3),
        ('L6D3', 6),
        ('M9D3', 9),
        ('H12D3_val', 12),
        ('VH15D3', 15),
    ]

    print("\n模擬結果:")
    print("-" * 100)
    header = f"{'處理組':>10} | {'時數':>4} | {'觀測Anth':>8} | {'模擬Anth':>10} | {'誤差':>10} | {'效率':>10}"
    print(header)
    print("-" * 100)

    for treatment, hours in treatments:
        obs_anth = TARGETS[treatment]['Anth']
        result = run_simulation_with_efficiency(treatment, eff_func)
        if result:
            error = (result['Anth'] - obs_anth) / obs_anth * 100
            print(f"{treatment:>10} | {hours:>4}h | {obs_anth:>8} | {result['Anth']:>10.1f} | {error:>+9.1f}% | {result['efficiency']*100:>9.1f}%")

    print("-" * 100)


if __name__ == "__main__":
    # 網格搜索
    best_params = grid_search()

    # 測試最佳參數
    test_best_params(*best_params)
