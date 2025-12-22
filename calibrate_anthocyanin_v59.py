"""
花青素參數校準腳本 v5.9
========================
目的: 校準新的 Stress 累積能量驅動的花青素合成機制

變更內容:
- 移除硬閾值 stress_threshold_anth = 15.0
- 新增狀態變量 E_stress (Stress 累積能量)
- 使用 Hill 函數: uva_induced = V_max × E^n / (K^n + E^n)

需要校準的參數:
- V_max_anth: 最大誘導合成率 [kg/m²/s]
- K_E_anth: 半飽和 Stress 累積能量 [Stress·day]
- n_hill_anth: Hill 係數 (控制曲線陡度)

目標:
- 改善 VL3D12 (誤差 +34%) 和 L6D12 (誤差 +40%) 的花青素預測
- 保持其他處理組的準確度
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution
import sys

# 導入模型
from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import (
    ENV_BASE, TREATMENT_CONFIGS, TARGETS, ODE_SETTINGS, SIMULATION,
    get_env_for_treatment
)


def run_simulation(treatment, p):
    """執行單一處理組的模擬"""
    env = get_env_for_treatment(treatment)

    # 初始條件
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6

    # v5.9: 新增 E_stress 初始值
    initial_state = [Xd_init, 0.0, LAI_init, Anth_init, 0.0, 0.0]

    simulation_days = SIMULATION['days']

    sol = solve_ivp(
        fun=uva_sun_derivatives,
        t_span=(0, simulation_days * 86400),
        y0=initial_state,
        args=(p, env),
        method=ODE_SETTINGS['method'],
        max_step=ODE_SETTINGS['max_step']
    )

    if not sol.success:
        return None, None, None, None

    # 計算結果
    sim_stress = sol.y[4, -1]
    sim_E_stress = sol.y[5, -1]  # v5.9: 新增
    sim_dw_fw_ratio = calculate_dynamic_dw_fw_ratio(sim_stress, p)
    sim_dw_per_plant = sol.y[0, -1] / ENV_BASE['plant_density'] * 1000
    sim_fw_per_plant = sim_dw_per_plant / sim_dw_fw_ratio

    # 花青素
    sim_anth_kg_m2 = sol.y[3, -1]
    sim_fw_kg_m2 = sim_fw_per_plant * ENV_BASE['plant_density'] / 1000
    sim_anth_ppm = sim_anth_kg_m2 * 1e6 / (sim_fw_kg_m2 + 1e-9)

    return sim_fw_per_plant, sim_anth_ppm, sim_stress, sim_E_stress


def objective_function(params):
    """
    目標函數: 最小化花青素誤差

    params: [V_max_anth, K_E_anth, n_hill_anth]
    """
    V_max_anth, K_E_anth, n_hill_anth = params

    # 更新參數
    p = UVAParams()
    p.V_max_anth = V_max_anth
    p.K_E_anth = K_E_anth
    p.n_hill_anth = n_hill_anth

    # 定義處理組和權重
    # 給予 VL3D12 和 L6D12 更高權重 (因為這兩組誤差最大)
    treatments_weights = {
        'CK': 1.0,
        'L6D6': 1.0,
        'L6D6-N': 1.0,
        'H12D3': 1.0,
        'VL3D12': 2.0,  # 高權重
        'L6D12': 2.0,   # 高權重
    }

    total_error = 0.0
    total_weight = 0.0

    for treatment, weight in treatments_weights.items():
        target = TARGETS[treatment]
        sim_fw, sim_anth, sim_stress, sim_E_stress = run_simulation(treatment, p)

        if sim_fw is None:
            return 1e6  # 模擬失敗

        # 計算花青素誤差 (平方誤差)
        anth_error = ((sim_anth - target['Anth']) / target['Anth']) ** 2

        # 同時考慮鮮重誤差 (較小權重，確保不破壞 FW 校準)
        fw_error = ((sim_fw - target['FW']) / target['FW']) ** 2

        # 加權誤差: 花青素權重 70%, 鮮重權重 30%
        combined_error = 0.7 * anth_error + 0.3 * fw_error

        total_error += weight * combined_error
        total_weight += weight

    return total_error / total_weight


def calibrate_anthocyanin_params():
    """使用差分進化法校準花青素參數"""

    print("=" * 80)
    print("花青素參數校準 v5.9 (Stress 累積能量驅動)")
    print("=" * 80)
    print("\n目標: 改善 VL3D12 和 L6D12 的花青素預測")
    print("方法: 差分進化法 (全域優化)")
    print("\n校準參數:")
    print("  - V_max_anth: 最大誘導合成率 [kg/m²/s]")
    print("  - K_E_anth: 半飽和 Stress 累積能量 [Stress·day]")
    print("  - n_hill_anth: Hill 係數")

    # 參數範圍
    bounds = [
        (1e-10, 1e-8),   # V_max_anth [kg/m²/s]
        (10.0, 200.0),   # K_E_anth [Stress·day]
        (1.0, 4.0),      # n_hill_anth [-]
    ]

    print("\n參數搜尋範圍:")
    print(f"  - V_max_anth: [{bounds[0][0]:.2e}, {bounds[0][1]:.2e}] kg/m²/s")
    print(f"  - K_E_anth: [{bounds[1][0]}, {bounds[1][1]}] Stress·day")
    print(f"  - n_hill_anth: [{bounds[2][0]}, {bounds[2][1]}]")

    print("\n開始優化...")

    result = differential_evolution(
        objective_function,
        bounds,
        strategy='best1bin',
        maxiter=50,
        popsize=15,
        tol=1e-4,
        seed=42,
        disp=True,
        workers=1
    )

    if result.success:
        print("\n" + "=" * 80)
        print("校準成功!")
        print("=" * 80)

        V_max_opt, K_E_opt, n_hill_opt = result.x

        print(f"\n最佳參數:")
        print(f"  - V_max_anth = {V_max_opt:.3e} kg/m²/s")
        print(f"  - K_E_anth = {K_E_opt:.2f} Stress·day")
        print(f"  - n_hill_anth = {n_hill_opt:.2f}")
        print(f"\n目標函數值 (加權RMSE): {result.fun:.6f}")

        # 驗證結果
        print("\n" + "=" * 80)
        print("驗證結果 (所有處理組)")
        print("=" * 80)

        p_opt = UVAParams()
        p_opt.V_max_anth = V_max_opt
        p_opt.K_E_anth = K_E_opt
        p_opt.n_hill_anth = n_hill_opt

        print(f"\n{'Treatment':<10} {'FW_obs':<8} {'FW_sim':<8} {'FW_Err':<8} "
              f"{'Anth_obs':<10} {'Anth_sim':<10} {'Anth_Err':<10} "
              f"{'Stress':<8} {'E_stress':<10}")
        print("-" * 90)

        fw_errors = []
        anth_errors = []

        for treatment in ['CK', 'L6D6', 'L6D6-N', 'H12D3', 'VL3D12', 'L6D12']:
            target = TARGETS[treatment]
            sim_fw, sim_anth, sim_stress, sim_E_stress = run_simulation(treatment, p_opt)

            if sim_fw is not None:
                fw_err = (sim_fw - target['FW']) / target['FW'] * 100
                anth_err = (sim_anth - target['Anth']) / target['Anth'] * 100

                fw_errors.append(abs(fw_err))
                anth_errors.append(abs(anth_err))

                print(f"{treatment:<10} {target['FW']:<8.1f} {sim_fw:<8.1f} {fw_err:<+8.1f} "
                      f"{target['Anth']:<10.1f} {sim_anth:<10.1f} {anth_err:<+10.1f} "
                      f"{sim_stress:<8.2f} {sim_E_stress:<10.2f}")

        print("-" * 90)
        print(f"{'Mean |Err|':<10} {'':<8} {'':<8} {np.mean(fw_errors):<8.1f} "
              f"{'':<10} {'':<10} {np.mean(anth_errors):<10.1f}")
        print(f"{'Max |Err|':<10} {'':<8} {'':<8} {np.max(fw_errors):<8.1f} "
              f"{'':<10} {'':<10} {np.max(anth_errors):<10.1f}")

        return p_opt
    else:
        print(f"\n優化失敗: {result.message}")
        return None


if __name__ == "__main__":
    calibrate_anthocyanin_params()
