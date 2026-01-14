#!/usr/bin/env python3
"""
測試移除非對稱高斯效率函數的影響

比較:
1. 當前模型 (有非對稱高斯)
2. 移除非對稱高斯 (效率固定為 1.0)
"""

import numpy as np
from scipy.integrate import solve_ivp

# 導入模型
from simulate_uva_model_v10 import (
    UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio,
    nonlinear_damage_factor, calculate_nonlin_anth_efficiency
)
from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment

# 照射時數 -> nonlinear_factor (from Gompertz)
HOURS_TO_NONLIN = {0: 1.0, 3: 1.0, 6: 1.8, 9: 70.2, 12: 188.7, 15: 235.5}


def run_simulation(treatment_name, disable_asymmetric=False):
    """執行單次模擬並返回結果"""
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

    # 如果要禁用非對稱高斯，暫時 monkey patch
    if disable_asymmetric:
        import simulate_uva_model_v10 as model
        original_func = model.calculate_nonlin_anth_efficiency
        model.calculate_nonlin_anth_efficiency = lambda nf, p: 1.0

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
    if disable_asymmetric:
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

    # 計算 nonlinear_factor
    uva_hour_on = env.get('uva_hour_on', 0)
    uva_hour_off = env.get('uva_hour_off', 0)
    if env.get('uva_on', False):
        daily_hours = uva_hour_off - uva_hour_on if uva_hour_on < uva_hour_off else 24 - uva_hour_on + uva_hour_off
    else:
        daily_hours = 0
    nonlin_factor = nonlinear_damage_factor(daily_hours, p)

    # 如果禁用非對稱高斯，確認 efficiency 是 1.0
    if disable_asymmetric:
        nonlin_eff = 1.0
    else:
        nonlin_eff = calculate_nonlin_anth_efficiency(nonlin_factor, p)

    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(avg_stress, p, nonlin_factor)

    FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
    FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
    Anth_sim = Anth_f / FW_total_kg * 1e6

    return {
        'FW_g': FW_sim,
        'Anth': Anth_sim,
        'LAI': LAI_f,
        'Stress': Stress_f,
        'avg_stress': avg_stress,
        'nonlin_factor': nonlin_factor,
        'nonlin_eff': nonlin_eff,
        'daily_hours': daily_hours
    }


def run_comparison():
    """比較有/無非對稱高斯的模擬結果"""

    print("=" * 100)
    print("測試移除非對稱高斯效率函數的影響")
    print("=" * 100)

    # 先顯示非對稱高斯對各小時的效率影響
    print("\n【非對稱高斯效率函數效果】")
    print("-" * 60)
    print(f"{'照射時數':>8} | {'NonlinFactor':>12} | {'當前效率':>10} | {'移除後效率':>10}")
    print("-" * 60)

    p = UVAParams()
    for hours in [0, 3, 6, 9, 12, 15]:
        nonlin = HOURS_TO_NONLIN[hours]
        current_eff = calculate_nonlin_anth_efficiency(nonlin, p)
        removed_eff = 1.0  # 移除後效率固定為 1.0
        print(f"{hours:>7}h | {nonlin:>12.1f} | {current_eff*100:>9.1f}% | {removed_eff*100:>9.1f}%")

    print("-" * 60)

    # 關鍵組別測試 (受非對稱高斯影響的)
    test_treatments = [
        ('M9D3', 9),
        ('H12D3_val', 12),
        ('VH15D3', 15),
    ]

    print("\n【高照射時數組別 (主要受影響)】")
    print("-" * 100)
    header = f"{'處理組':>10} | {'時數':>4} | {'觀測Anth':>8} | {'當前模擬':>10} | {'當前誤差':>10} | {'移除後模擬':>10} | {'移除後誤差':>10}"
    print(header)
    print("-" * 100)

    for treatment, hours in test_treatments:
        obs_anth = TARGETS[treatment]['Anth']

        # 當前模型
        result_current = run_simulation(treatment, disable_asymmetric=False)
        anth_current = result_current['Anth']
        error_current = (anth_current - obs_anth) / obs_anth * 100

        # 移除非對稱高斯
        result_removed = run_simulation(treatment, disable_asymmetric=True)
        anth_removed = result_removed['Anth']
        error_removed = (anth_removed - obs_anth) / obs_anth * 100

        print(f"{treatment:>10} | {hours:>4}h | {obs_anth:>8} | {anth_current:>10.1f} | {error_current:>+9.1f}% | {anth_removed:>10.1f} | {error_removed:>+9.1f}%")

    print("-" * 100)

    # 也測試低時數組別確認不受影響
    low_hour_treatments = [
        ('VL3D3', 3),
        ('L6D3', 6),
        ('L6D6', 6),
    ]

    print("\n【低照射時數組別 (應不受影響)】")
    print("-" * 100)
    print(header)
    print("-" * 100)

    for treatment, hours in low_hour_treatments:
        obs_anth = TARGETS[treatment]['Anth']

        # 當前模型
        result_current = run_simulation(treatment, disable_asymmetric=False)
        anth_current = result_current['Anth']
        error_current = (anth_current - obs_anth) / obs_anth * 100

        # 移除非對稱高斯
        result_removed = run_simulation(treatment, disable_asymmetric=True)
        anth_removed = result_removed['Anth']
        error_removed = (anth_removed - obs_anth) / obs_anth * 100

        print(f"{treatment:>10} | {hours:>4}h | {obs_anth:>8} | {anth_current:>10.1f} | {error_current:>+9.1f}% | {anth_removed:>10.1f} | {error_removed:>+9.1f}%")

    print("-" * 100)

    print("\n【結論】")
    print("1. 如果移除後誤差仍在 ±15% 內，建議移除此機制")
    print("2. 非對稱高斯的生物學解釋不合理 (9h 是最低點，但 12h/15h 恢復)")
    print("3. 這是一個純數學補丁，缺乏生物學基礎")


if __name__ == "__main__":
    run_comparison()
