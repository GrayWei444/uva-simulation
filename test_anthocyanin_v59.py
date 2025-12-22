"""
花青素機制測試腳本 v5.9
======================
快速測試新的 Stress 累積能量驅動機制

使用方法:
1. 手動調整 UVAParams 中的花青素參數
2. 運行此腳本查看所有處理組的預測結果
3. 逐步調整參數直到滿意

不使用優化算法，便於理解參數對結果的影響
"""

import numpy as np
from scipy.integrate import solve_ivp

# 導入模型
from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import (
    ENV_BASE, TREATMENT_CONFIGS, TARGETS, ODE_SETTINGS, SIMULATION,
    get_env_for_treatment
)


def run_all_treatments(p):
    """運行所有處理組並顯示結果"""

    treatments = ['CK', 'L6D6', 'L6D6-N', 'H12D3', 'VL3D12', 'L6D12']

    print("=" * 100)
    print("所有處理組模擬結果 (v5.9 - Stress 累積能量驅動花青素)")
    print("=" * 100)

    print(f"\n當前花青素參數:")
    print(f"  - base_anth_rate_light = {p.base_anth_rate_light:.3e} kg/m²/s")
    print(f"  - base_anth_rate_dark = {p.base_anth_rate_dark:.3e} kg/m²/s")
    print(f"  - V_max_anth = {p.V_max_anth:.3e} kg/m²/s")
    print(f"  - K_E_anth = {p.K_E_anth:.2f} Stress·day")
    print(f"  - n_hill_anth = {p.n_hill_anth:.2f}")
    print(f"  - k_deg = {p.k_deg:.3e} 1/s")

    print(f"\n{'Treatment':<10} {'FW_obs':<8} {'FW_sim':<8} {'FW_Err':<9} "
          f"{'Anth_obs':<10} {'Anth_sim':<10} {'Anth_Err':<11} "
          f"{'Stress':<8} {'E_stress':<10}")
    print("-" * 100)

    fw_errors = []
    anth_errors = []
    results = {}

    for treatment in treatments:
        target = TARGETS[treatment]
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

        # 執行模擬
        sol = solve_ivp(
            fun=uva_sun_derivatives,
            t_span=(0, simulation_days * 86400),
            y0=initial_state,
            args=(p, env),
            method=ODE_SETTINGS['method'],
            max_step=ODE_SETTINGS['max_step']
        )

        if sol.success:
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

            # 誤差
            fw_err = (sim_fw_per_plant - target['FW']) / target['FW'] * 100
            anth_err = (sim_anth_ppm - target['Anth']) / target['Anth'] * 100

            fw_errors.append(abs(fw_err))
            anth_errors.append(abs(anth_err))

            # 儲存結果
            results[treatment] = {
                'FW_sim': sim_fw_per_plant,
                'Anth_sim': sim_anth_ppm,
                'Stress': sim_stress,
                'E_stress': sim_E_stress
            }

            print(f"{treatment:<10} {target['FW']:<8.1f} {sim_fw_per_plant:<8.1f} {fw_err:<+9.1f} "
                  f"{target['Anth']:<10.1f} {sim_anth_ppm:<10.1f} {anth_err:<+11.1f} "
                  f"{sim_stress:<8.2f} {sim_E_stress:<10.2f}")
        else:
            print(f"{treatment:<10} 模擬失敗: {sol.message}")

    print("-" * 100)
    print(f"{'Mean |Err|':<10} {'':<8} {'':<8} {np.mean(fw_errors):<9.1f} "
          f"{'':<10} {'':<10} {np.mean(anth_errors):<11.1f}")
    print(f"{'Max |Err|':<10} {'':<8} {'':<8} {np.max(fw_errors):<9.1f} "
          f"{'':<10} {'':<10} {np.max(anth_errors):<11.1f}")
    print("=" * 100)

    # 詳細分析
    print("\n" + "=" * 100)
    print("E_stress 與花青素關係分析")
    print("=" * 100)

    print(f"\n{'Treatment':<10} {'E_stress':<12} {'Anth_obs':<12} {'Anth_sim':<12} "
          f"{'Hill函數值':<12} {'誘導合成率':<15}")
    print("-" * 80)

    for treatment in treatments:
        if treatment in results:
            r = results[treatment]
            target = TARGETS[treatment]
            E = r['E_stress']

            # 計算 Hill 函數值
            E_power_n = E ** p.n_hill_anth
            K_power_n = p.K_E_anth ** p.n_hill_anth
            hill_value = E_power_n / (K_power_n + E_power_n + 1e-12)

            # 誘導合成率
            uva_induced = p.V_max_anth * hill_value

            print(f"{treatment:<10} {E:<12.2f} {target['Anth']:<12.1f} {r['Anth_sim']:<12.1f} "
                  f"{hill_value:<12.3f} {uva_induced:<15.3e}")

    print("=" * 100)

    return results


def test_parameter_sensitivity():
    """測試參數敏感度"""

    print("\n" + "=" * 100)
    print("參數敏感度測試")
    print("=" * 100)

    # 測試 K_E_anth 的影響
    print("\n測試 K_E_anth (半飽和常數) 的影響:")
    print(f"{'K_E_anth':<10} {'VL3D12_Anth':<15} {'L6D12_Anth':<15} {'H12D3_Anth':<15}")
    print("-" * 60)

    for K_E in [20, 30, 40, 50, 60, 80, 100]:
        p = UVAParams()
        p.K_E_anth = K_E
        p.V_max_anth = 3.0e-10
        p.n_hill_anth = 2.0

        results_vl = run_single_treatment('VL3D12', p)
        results_l6 = run_single_treatment('L6D12', p)
        results_h12 = run_single_treatment('H12D3', p)

        if results_vl and results_l6 and results_h12:
            print(f"{K_E:<10.0f} {results_vl['Anth_sim']:<15.1f} "
                  f"{results_l6['Anth_sim']:<15.1f} {results_h12['Anth_sim']:<15.1f}")

    # 測試 n_hill_anth 的影響
    print("\n測試 n_hill_anth (Hill 係數) 的影響:")
    print(f"{'n_hill':<10} {'VL3D12_Anth':<15} {'L6D12_Anth':<15} {'H12D3_Anth':<15}")
    print("-" * 60)

    for n_hill in [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]:
        p = UVAParams()
        p.K_E_anth = 40.0
        p.V_max_anth = 3.0e-10
        p.n_hill_anth = n_hill

        results_vl = run_single_treatment('VL3D12', p)
        results_l6 = run_single_treatment('L6D12', p)
        results_h12 = run_single_treatment('H12D3', p)

        if results_vl and results_l6 and results_h12:
            print(f"{n_hill:<10.1f} {results_vl['Anth_sim']:<15.1f} "
                  f"{results_l6['Anth_sim']:<15.1f} {results_h12['Anth_sim']:<15.1f}")

    print("=" * 100)


def run_single_treatment(treatment, p):
    """運行單一處理組（用於敏感度測試）"""
    env = get_env_for_treatment(treatment)

    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6
    initial_state = [Xd_init, 0.0, LAI_init, Anth_init, 0.0, 0.0]

    sol = solve_ivp(
        fun=uva_sun_derivatives,
        t_span=(0, SIMULATION['days'] * 86400),
        y0=initial_state,
        args=(p, env),
        method=ODE_SETTINGS['method'],
        max_step=ODE_SETTINGS['max_step']
    )

    if sol.success:
        sim_stress = sol.y[4, -1]
        sim_E_stress = sol.y[5, -1]
        sim_dw_fw_ratio = calculate_dynamic_dw_fw_ratio(sim_stress, p)
        sim_dw_per_plant = sol.y[0, -1] / ENV_BASE['plant_density'] * 1000
        sim_fw_per_plant = sim_dw_per_plant / sim_dw_fw_ratio
        sim_anth_kg_m2 = sol.y[3, -1]
        sim_fw_kg_m2 = sim_fw_per_plant * ENV_BASE['plant_density'] / 1000
        sim_anth_ppm = sim_anth_kg_m2 * 1e6 / (sim_fw_kg_m2 + 1e-9)

        return {
            'FW_sim': sim_fw_per_plant,
            'Anth_sim': sim_anth_ppm,
            'Stress': sim_stress,
            'E_stress': sim_E_stress
        }
    return None


if __name__ == "__main__":
    # 測試當前參數
    p = UVAParams()

    print("測試 1: 當前參數設定")
    run_all_treatments(p)

    # 測試建議參數
    print("\n\n測試 2: 建議參數設定 (基於理論分析)")
    p2 = UVAParams()
    p2.V_max_anth = 3.0e-10
    p2.K_E_anth = 40.0
    p2.n_hill_anth = 2.5
    run_all_treatments(p2)

    # 參數敏感度測試
    test_parameter_sensitivity()
