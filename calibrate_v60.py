#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
v6.0 花青素參數自動校準腳本
目標: 所有處理組花青素誤差 < 5%
"""

import numpy as np
from scipy.optimize import differential_evolution
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import get_env_for_treatment, SIMULATION, ENV_BASE, TARGETS

def simulate_treatment(treatment, params_dict):
    """
    模擬單一處理組

    參數:
        treatment: 處理組名稱
        params_dict: 字典 {'V_max': ..., 'K_stress': ..., 'n_stress': ..., 'k_deg': ...}

    回傳:
        (FW_sim, Anth_sim, Stress_final)
    """
    p = UVAParams()

    # 更新參數
    p.V_max_anth = params_dict['V_max']
    p.K_stress_anth = params_dict['K_stress']
    p.n_stress_anth = params_dict['n_stress']
    p.k_deg = params_dict['k_deg']

    env = get_env_for_treatment(treatment)

    # 初始狀態
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6
    y0 = [Xd_init, 0.0, LAI_init, Anth_init, 0.0, 0.0]

    # 模擬
    sol = solve_ivp(
        lambda t, y: uva_sun_derivatives(t, y, p, env),
        (0, SIMULATION['days'] * 86400),
        y0,
        method='RK45',
        max_step=60
    )

    if not sol.success:
        return None, None, None

    X_d_final = sol.y[0, -1]
    Anth_final = sol.y[3, -1]
    Stress_final = sol.y[4, -1]

    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress_final, p)
    FW_sim = (X_d_final / ENV_BASE['plant_density'] * 1000) / dw_fw_ratio

    fw_kg_m2 = FW_sim * ENV_BASE['plant_density'] / 1000
    Anth_sim = Anth_final * 1e6 / (fw_kg_m2 + 1e-9)

    return FW_sim, Anth_sim, Stress_final


def objective_function(x):
    """
    目標函數: 最小化花青素誤差

    x = [V_max, K_stress, n_stress, k_deg]
    """
    params_dict = {
        'V_max': x[0],
        'K_stress': x[1],
        'n_stress': x[2],
        'k_deg': x[3]
    }

    treatments = ['CK', 'L6D6', 'L6D6-N', 'H12D3', 'VL3D12', 'L6D12']
    errors = []

    for treatment in treatments:
        FW_sim, Anth_sim, _ = simulate_treatment(treatment, params_dict)

        if Anth_sim is None:
            return 1e6  # 懲罰失敗的模擬

        target = TARGETS[treatment]
        anth_err = abs((Anth_sim - target['Anth']) / target['Anth'] * 100)
        errors.append(anth_err)

    # 返回平均絕對誤差
    return np.mean(errors)


def main():
    print("=" * 80)
    print("v6.0 花青素參數自動校準")
    print("=" * 80)
    print("\n搜尋範圍:")
    print("  V_max:    1e-10 ~ 5e-10 kg/m²/s")
    print("  K_stress: 5.0 ~ 30.0")
    print("  n_stress: 1.5 ~ 3.5")
    print("  k_deg:    1e-6 ~ 5e-6 1/s")
    print()

    # 定義參數邊界
    bounds = [
        (1e-10, 5e-10),  # V_max
        (5.0, 30.0),     # K_stress
        (1.5, 3.5),      # n_stress
        (1e-6, 5e-6)     # k_deg
    ]

    print("開始全域優化 (差分進化法)...")
    print("這可能需要幾分鐘...\n")

    result = differential_evolution(
        objective_function,
        bounds,
        strategy='best1bin',
        maxiter=100,
        popsize=15,
        tol=0.01,
        mutation=(0.5, 1.0),
        recombination=0.7,
        seed=42,
        polish=True,
        workers=1,
        updating='deferred',
        disp=True
    )

    if result.success:
        print("\n" + "=" * 80)
        print("校準完成！")
        print("=" * 80)

        best_params = {
            'V_max': result.x[0],
            'K_stress': result.x[1],
            'n_stress': result.x[2],
            'k_deg': result.x[3]
        }

        print("\n最佳參數:")
        print(f"  V_max_anth = {best_params['V_max']:.2e}  # kg/m²/s")
        print(f"  K_stress_anth = {best_params['K_stress']:.2f}")
        print(f"  n_stress_anth = {best_params['n_stress']:.2f}")
        print(f"  k_deg = {best_params['k_deg']:.2e}  # 1/s")
        print(f"\n平均花青素誤差: {result.fun:.2f}%")

        # 驗證所有處理組
        print("\n" + "=" * 80)
        print("驗證結果:")
        print("=" * 80)
        treatments = ['CK', 'L6D6', 'L6D6-N', 'H12D3', 'VL3D12', 'L6D12']

        print(f"{'Treatment':<10} {'Anth_sim':>9} {'Anth_exp':>9} {'Error':>8}")
        print("-" * 40)

        all_errors = []
        for treatment in treatments:
            FW_sim, Anth_sim, _ = simulate_treatment(treatment, best_params)
            target = TARGETS[treatment]
            error = (Anth_sim - target['Anth']) / target['Anth'] * 100
            all_errors.append(abs(error))

            marker = '✓' if abs(error) < 5.0 else '✗'
            print(f"{treatment:<10} {Anth_sim:>9.1f} {target['Anth']:>9.1f} {error:>7.1f}% {marker}")

        print("-" * 40)
        print(f"Mean |Error|: {np.mean(all_errors):.2f}%")
        print(f"Max |Error|:  {np.max(all_errors):.2f}%")
        print(f"誤差 < 5% 的處理組: {sum(1 for e in all_errors if e < 5.0)}/6")

    else:
        print("\n優化失敗:", result.message)
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
