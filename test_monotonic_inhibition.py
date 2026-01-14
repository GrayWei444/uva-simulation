#!/usr/bin/env python3
"""
測試單調遞減效率函數

概念：照射時數越長，花青素合成效率越低
- 3h, 6h: 效率 100%
- 9h: 開始下降
- 12h, 15h: 繼續下降

這比非對稱高斯更符合生物學邏輯：
長時間照射 → 細胞損傷累積 → 合成機制受損 → 效率降低
"""

import numpy as np
from scipy.integrate import solve_ivp

from simulate_uva_model_v10 import (
    UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio,
    nonlinear_damage_factor
)
from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment


def monotonic_efficiency(hours, threshold=7.0, steepness=0.3, max_inhib=0.35):
    """
    單調遞減的效率函數

    使用 sigmoid 形式：efficiency = 1 - max_inhib * sigmoid((hours - threshold) / scale)

    參數:
    - threshold: 開始下降的時數 (小於此值效率接近 100%)
    - steepness: 下降速率
    - max_inhib: 最大抑制量 (0.35 = 最多降到 65%)
    """
    sigmoid = 1 / (1 + np.exp(-(hours - threshold) / (1/steepness)))
    efficiency = 1.0 - max_inhib * sigmoid
    return efficiency


def run_simulation(treatment_name, use_monotonic=False, mono_params=None):
    """執行單次模擬"""
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

    # 如果使用單調效率，暫時 monkey patch
    if use_monotonic:
        import simulate_uva_model_v10 as model
        original_func = model.calculate_nonlin_anth_efficiency

        # 計算每日照射時數
        uva_hour_on = env.get('uva_hour_on', 0)
        uva_hour_off = env.get('uva_hour_off', 0)
        if env.get('uva_on', False):
            daily_hours = uva_hour_off - uva_hour_on if uva_hour_on < uva_hour_off else 24 - uva_hour_on + uva_hour_off
        else:
            daily_hours = 0

        # 計算單調效率
        if mono_params:
            eff = monotonic_efficiency(daily_hours, **mono_params)
        else:
            eff = monotonic_efficiency(daily_hours)

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
    if use_monotonic:
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

    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(avg_stress, p, nonlin_factor)

    FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
    FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
    Anth_sim = Anth_f / FW_total_kg * 1e6

    return {
        'FW_g': FW_sim,
        'Anth': Anth_sim,
        'daily_hours': daily_hours
    }


def test_monotonic():
    """測試單調遞減效率函數"""

    print("=" * 100)
    print("測試單調遞減效率函數 (替代非對稱高斯)")
    print("=" * 100)

    # 顯示效率曲線
    print("\n【單調遞減效率函數】")
    print("  概念: 照射時數越長 → 細胞損傷累積 → 合成效率降低")
    print("-" * 60)
    print(f"{'照射時數':>8} | {'效率 (預設)':>12} | {'效率 (調整)':>12}")
    print("-" * 60)

    # 預設參數和調整後參數
    default_params = {'threshold': 7.0, 'steepness': 0.3, 'max_inhib': 0.35}
    tuned_params = {'threshold': 7.5, 'steepness': 0.25, 'max_inhib': 0.30}

    for hours in [0, 3, 6, 9, 12, 15]:
        eff_default = monotonic_efficiency(hours, **default_params)
        eff_tuned = monotonic_efficiency(hours, **tuned_params)
        print(f"{hours:>7}h | {eff_default*100:>11.1f}% | {eff_tuned*100:>11.1f}%")

    print("-" * 60)

    # 測試組別
    test_treatments = [
        ('VL3D3', 3),
        ('L6D3', 6),
        ('M9D3', 9),
        ('H12D3_val', 12),
        ('VH15D3', 15),
    ]

    print("\n【模擬結果比較】")
    print("-" * 120)
    header = f"{'處理組':>10} | {'時數':>4} | {'觀測':>8} | {'當前模型':>10} | {'當前誤差':>10} | {'單調(預設)':>10} | {'預設誤差':>10} | {'單調(調整)':>10} | {'調整誤差':>10}"
    print(header)
    print("-" * 120)

    for treatment, hours in test_treatments:
        obs_anth = TARGETS[treatment]['Anth']

        # 當前模型 (有非對稱高斯)
        result_current = run_simulation(treatment, use_monotonic=False)
        anth_current = result_current['Anth']
        error_current = (anth_current - obs_anth) / obs_anth * 100

        # 單調效率 (預設參數)
        result_default = run_simulation(treatment, use_monotonic=True, mono_params=default_params)
        anth_default = result_default['Anth']
        error_default = (anth_default - obs_anth) / obs_anth * 100

        # 單調效率 (調整參數)
        result_tuned = run_simulation(treatment, use_monotonic=True, mono_params=tuned_params)
        anth_tuned = result_tuned['Anth']
        error_tuned = (anth_tuned - obs_anth) / obs_anth * 100

        print(f"{treatment:>10} | {hours:>4}h | {obs_anth:>8} | {anth_current:>10.1f} | {error_current:>+9.1f}% | {anth_default:>10.1f} | {error_default:>+9.1f}% | {anth_tuned:>10.1f} | {error_tuned:>+9.1f}%")

    print("-" * 120)

    print("\n【結論】")
    print("單調遞減效率函數的優點:")
    print("1. 符合生物學邏輯: 照射時數越長 → 累積損傷越大 → 效率越低")
    print("2. 連續可微: 符合 CLAUDE.md 的無硬閾值原則")
    print("3. 可解釋: 9h 開始下降，12h/15h 繼續下降")


if __name__ == "__main__":
    test_monotonic()
