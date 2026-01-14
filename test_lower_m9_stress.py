#!/usr/bin/env python3
"""
測試降低 M9D3 的 Stress

策略: 讓 Gompertz 轉折點往後移
- 原始: threshold=9.5 → 9h 的 nonlin_factor=70.2
- 新: threshold=10.5 或更高 → 9h 的 nonlin_factor 更低

這樣 M9D3 的 Stress 累積更少，花青素合成自然降低
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


def run_simulation_modified(treatment_name, threshold, steepness, max_factor):
    """使用修改的參數執行模擬"""
    import simulate_uva_model_v10 as model

    # 保存原參數
    orig = {
        'threshold': model.ALL_PARAMS['gompertz_threshold'],
        'steepness': model.ALL_PARAMS['gompertz_steepness'],
        'max_factor': model.ALL_PARAMS['gompertz_max_factor'],
    }

    # 修改參數
    model.ALL_PARAMS['gompertz_threshold'] = threshold
    model.ALL_PARAMS['gompertz_steepness'] = steepness
    model.ALL_PARAMS['gompertz_max_factor'] = max_factor

    # 禁用非對稱高斯
    orig_nonlin_func = model.calculate_nonlin_anth_efficiency
    model.calculate_nonlin_anth_efficiency = lambda nf, p: 1.0

    try:
        p = UVAParams()
        env = get_env_for_treatment(treatment_name)

        uva_hour_on = env.get('uva_hour_on', 0)
        uva_hour_off = env.get('uva_hour_off', 0)
        if env.get('uva_on', False):
            daily_hours = uva_hour_off - uva_hour_on if uva_hour_on < uva_hour_off else 24 - uva_hour_on + uva_hour_off
        else:
            daily_hours = 0

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
            max_step=300
        )

        if not sol.success:
            return None

        Xd_f, Cbuf_f, LAI_f, Anth_f, Stress_f, ROS_f = sol.y[:, -1]

        uva_start_day = env.get('uva_start_day', 29)
        uva_start = uva_start_day * 86400
        stress_sum = 0
        stress_count = 0
        for i in range(len(sol.t)):
            if sol.t[i] >= uva_start:
                stress_sum += sol.y[4, i]
                stress_count += 1
        avg_stress = stress_sum / max(1, stress_count)

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
        }

    finally:
        model.ALL_PARAMS['gompertz_threshold'] = orig['threshold']
        model.ALL_PARAMS['gompertz_steepness'] = orig['steepness']
        model.ALL_PARAMS['gompertz_max_factor'] = orig['max_factor']
        model.calculate_nonlin_anth_efficiency = orig_nonlin_func


def main():
    print("=" * 120)
    print("測試降低 M9D3 Stress (讓 Gompertz 轉折點往後移)")
    print("=" * 120)

    treatments = [
        ('L6D3', 6),
        ('M9D3', 9),
        ('H12D3_val', 12),
        ('VH15D3', 15),
    ]

    # 測試不同 threshold
    test_params = [
        (9.5, 0.5, 250, "原始"),
        (10.0, 0.5, 250, "threshold +0.5"),
        (10.5, 0.5, 250, "threshold +1.0"),
        (11.0, 0.5, 250, "threshold +1.5"),
        (10.5, 0.6, 250, "threshold +1.0, steepness +0.1"),
    ]

    for threshold, steepness, max_factor, label in test_params:
        print(f"\n【{label}: threshold={threshold}, steepness={steepness}】")
        print("-" * 60)
        print("NonlinFactor: ", end="")
        for hours in [6, 9, 12, 15]:
            nf = gompertz_factor(hours, threshold, steepness, max_factor)
            print(f"{hours}h={nf:.1f}  ", end="")
        print()
        print("-" * 100)
        header = f"{'處理組':>10} | {'觀測Anth':>8} | {'模擬Anth':>10} | {'誤差':>10} | {'avgStress':>10}"
        print(header)
        print("-" * 100)

        for treatment, hours in treatments:
            obs_anth = TARGETS[treatment]['Anth']
            result = run_simulation_modified(treatment, threshold, steepness, max_factor)
            if result:
                error = (result['Anth'] - obs_anth) / obs_anth * 100
                print(f"{treatment:>10} | {obs_anth:>8} | {result['Anth']:>10.1f} | {error:>+9.1f}% | {result['avg_stress']:>10.1f}")

    print("\n" + "=" * 120)


if __name__ == "__main__":
    main()
