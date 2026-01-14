#!/usr/bin/env python3
"""
測試調整 Gompertz 參數對 M9D3 的影響

目標:
- 讓 9h 的 nonlinear_factor 更高
- 這樣 M9D3 的 Stress 更高
- Stress 抑制機制會降低花青素
- 實現「轉折點」效果
"""

import numpy as np
from scipy.integrate import solve_ivp

from simulate_uva_model_v10 import (
    UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
)
from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment


def gompertz_factor(hours, threshold, steepness, max_factor):
    """計算 Gompertz 非線性因子"""
    exponent = -steepness * (hours - threshold)
    factor = 1.0 + max_factor * np.exp(-np.exp(exponent))
    return factor


def run_simulation_with_gompertz(treatment_name, threshold, steepness, max_factor):
    """使用修改的 Gompertz 參數執行模擬"""
    import simulate_uva_model_v10 as model

    # 保存原參數
    original_threshold = model.ALL_PARAMS['gompertz_threshold']
    original_steepness = model.ALL_PARAMS['gompertz_steepness']
    original_max = model.ALL_PARAMS['gompertz_max_factor']

    # 修改參數
    model.ALL_PARAMS['gompertz_threshold'] = threshold
    model.ALL_PARAMS['gompertz_steepness'] = steepness
    model.ALL_PARAMS['gompertz_max_factor'] = max_factor

    # 同時禁用非對稱高斯
    original_nonlin_func = model.calculate_nonlin_anth_efficiency
    model.calculate_nonlin_anth_efficiency = lambda nf, p: 1.0

    try:
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

        # 重新導入以獲取更新後的 nonlinear_damage_factor
        from simulate_uva_model_v10 import nonlinear_damage_factor
        nonlin_factor = nonlinear_damage_factor(daily_hours, p)

        dw_fw_ratio = calculate_dynamic_dw_fw_ratio(avg_stress, p, nonlin_factor)

        FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
        FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
        Anth_sim = Anth_f / FW_total_kg * 1e6

        return {
            'FW_g': FW_sim,
            'Anth': Anth_sim,
            'avg_stress': avg_stress,
            'nonlin_factor': nonlin_factor,
            'daily_hours': daily_hours
        }

    finally:
        # 恢復原參數
        model.ALL_PARAMS['gompertz_threshold'] = original_threshold
        model.ALL_PARAMS['gompertz_steepness'] = original_steepness
        model.ALL_PARAMS['gompertz_max_factor'] = original_max
        model.calculate_nonlin_anth_efficiency = original_nonlin_func


def main():
    print("=" * 120)
    print("測試調整 Gompertz 參數")
    print("=" * 120)

    # 原始參數
    orig_threshold = 9.5
    orig_steepness = 0.5
    orig_max = 250

    print("\n【原始 Gompertz 參數】")
    print(f"threshold={orig_threshold}, steepness={orig_steepness}, max_factor={orig_max}")
    print("-" * 60)
    for hours in [3, 6, 9, 12, 15]:
        factor = gompertz_factor(hours, orig_threshold, orig_steepness, orig_max)
        print(f"  {hours}h: nonlin_factor = {factor:.1f}")

    # 測試原始參數 (無非對稱高斯)
    treatments = [
        ('VL3D3', 3),
        ('L6D3', 6),
        ('M9D3', 9),
        ('H12D3_val', 12),
        ('VH15D3', 15),
    ]

    print("\n【原始參數模擬結果 (無非對稱高斯)】")
    print("-" * 120)
    header = f"{'處理組':>10} | {'時數':>4} | {'觀測Anth':>8} | {'模擬Anth':>10} | {'誤差':>10} | {'avgStress':>10} | {'NonlinF':>10}"
    print(header)
    print("-" * 120)

    for treatment, hours in treatments:
        obs_anth = TARGETS[treatment]['Anth']
        result = run_simulation_with_gompertz(treatment, orig_threshold, orig_steepness, orig_max)
        if result:
            error = (result['Anth'] - obs_anth) / obs_anth * 100
            print(f"{treatment:>10} | {hours:>4}h | {obs_anth:>8} | {result['Anth']:>10.1f} | {error:>+9.1f}% | {result['avg_stress']:>10.1f} | {result['nonlin_factor']:>10.1f}")

    # 測試新參數: threshold=8.0, steepness=0.6
    new_threshold = 8.0
    new_steepness = 0.6
    new_max = 250

    print("\n\n【新 Gompertz 參數】")
    print(f"threshold={new_threshold}, steepness={new_steepness}, max_factor={new_max}")
    print("-" * 60)
    for hours in [3, 6, 9, 12, 15]:
        factor = gompertz_factor(hours, new_threshold, new_steepness, new_max)
        print(f"  {hours}h: nonlin_factor = {factor:.1f}")

    print("\n【新參數模擬結果】")
    print("-" * 120)
    print(header)
    print("-" * 120)

    for treatment, hours in treatments:
        obs_anth = TARGETS[treatment]['Anth']
        result = run_simulation_with_gompertz(treatment, new_threshold, new_steepness, new_max)
        if result:
            error = (result['Anth'] - obs_anth) / obs_anth * 100
            print(f"{treatment:>10} | {hours:>4}h | {obs_anth:>8} | {result['Anth']:>10.1f} | {error:>+9.1f}% | {result['avg_stress']:>10.1f} | {result['nonlin_factor']:>10.1f}")

    print("-" * 120)


if __name__ == "__main__":
    main()
