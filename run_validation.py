# -*- coding: utf-8 -*-
"""
================================================================================
處理組驗證腳本 (Treatment Validation Script)
================================================================================
版本: v3.0
日期: 2025-12-07

功能說明:
---------
驗證所有處理組的模擬結果是否符合實驗目標值。
輸出包含:
  1. 每個處理組的鮮重 (FW) 和花青素 (Anth) 模擬值與目標值比較
  2. 訓練集的平均誤差統計

處理組列表:
-----------
  全部: CK, L6D6, L6D6-N, H12D3, VL3D12, L6D12

用法:
-----
    python run_validation.py

依賴模組:
---------
  - model_config: 共用設定 (環境參數、處理組配置、目標值)
  - simulate_uva_model: UVA效應模型 v5.0 (5個狀態變量)
================================================================================
"""

import numpy as np
from scipy.integrate import solve_ivp

# 匯入共用設定模組
from model_config import (
    ENV_BASE, TREATMENT_CONFIGS, TARGETS, ODE_SETTINGS, SIMULATION,
    get_env_for_treatment
)

# 匯入模型 (v5.0 簡化版)
from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio


def calculate_initial_conditions(p):
    """
    計算模擬的初始條件 (v5.0: 5個狀態變量)

    參數:
    -----
    p : UVAParams
        模型參數物件

    回傳:
    -----
    list
        5個狀態變量的初始值
        [X_d, C_buf, LAI, Anth, Stress]
    """
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    LAI_init = (dw_init_g / 0.01) * 0.04
    Anth_init = 2.0 * (fw_init_g * ENV_BASE['plant_density'] / 1000) / 1000

    # v5.0: 5個狀態變量 [X_d, C_buf, LAI, Anth, Stress]
    return [Xd_init, 0.0, LAI_init, Anth_init, 0.0]


def simulate_treatment(treatment, p, initial_state):
    """
    對單一處理組執行模擬

    參數:
    -----
    treatment : str
        處理組代碼 (如 'CK', 'L6D6' 等)
    p : UVAParams
        模型參數物件
    initial_state : list
        初始狀態變量 (5個)

    回傳:
    -----
    dict or None
        模擬結果字典，包含 sim_fw, sim_anth, fw_err, anth_err, final_stress
        模擬失敗時回傳 None
    """
    env = get_env_for_treatment(treatment)

    # 執行 ODE 求解
    sol = solve_ivp(
        fun=uva_sun_derivatives,
        t_span=(0, SIMULATION['days'] * 86400),
        y0=initial_state,
        args=(p, env),
        method=ODE_SETTINGS['method'],
        max_step=ODE_SETTINGS['max_step'],
        t_eval=np.array([SIMULATION['days'] * 86400])
    )

    if not sol.success:
        return None

    # v5.0: Stress 在 index 4
    sim_Stress = sol.y[4, -1]

    # 使用 Stress 計算 DW:FW 比例
    sim_dw_fw_ratio = calculate_dynamic_dw_fw_ratio(sim_Stress, p)
    sim_dw_per_plant = sol.y[0, -1] / ENV_BASE['plant_density'] * 1000
    sim_fw = sim_dw_per_plant / sim_dw_fw_ratio

    # 計算花青素濃度 (ppm)
    sim_fw_kg_per_m2 = sim_fw * ENV_BASE['plant_density'] / 1000
    sim_anth_mg_per_m2 = sol.y[3, -1] * 1000
    sim_anth_ppm = sim_anth_mg_per_m2 / (sim_fw_kg_per_m2 + 1e-9)

    # 計算誤差
    target = TARGETS[treatment]
    fw_err = (sim_fw - target['FW']) / target['FW'] * 100
    anth_err = (sim_anth_ppm - target['Anth']) / target['Anth'] * 100

    return {
        'sim_fw': sim_fw,
        'sim_anth': sim_anth_ppm,
        'fw_err': fw_err,
        'anth_err': anth_err,
        'final_stress': sim_Stress
    }


def run_validation():
    """
    驗證所有處理組

    執行流程:
    ---------
    1. 初始化模型參數
    2. 計算初始條件
    3. 對每個處理組執行模擬
    4. 輸出結果表格
    5. 計算統計摘要

    回傳:
    -----
    dict
        所有處理組的模擬結果
    """
    # 初始化
    p = UVAParams()
    initial_state = calculate_initial_conditions(p)

    # 輸出標頭
    print('=' * 90)
    print('Treatment Validation - Model v5.0 (Simplified)')
    print('=' * 90)
    print(f"Environment: CO2_day={ENV_BASE['CO2_day']}ppm, I_day={ENV_BASE['I_day']} umol/m2/s")
    print(f"ODE: method={ODE_SETTINGS['method']}, max_step={ODE_SETTINGS['max_step']}s")
    print('=' * 90)
    print(f'Treatment  | Target_FW | Sim_FW   | FW_Err   | Target_Anth | Sim_Anth   | Anth_Err | Stress')
    print('-' * 90)

    # 執行模擬
    results = {}
    all_treatments = ['CK', 'L6D6', 'L6D6-N', 'H12D3', 'VL3D12', 'L6D12']

    for treatment in all_treatments:
        result = simulate_treatment(treatment, p, initial_state)

        if result is not None:
            results[treatment] = result
            target = TARGETS[treatment]
            print(f'{treatment:<10} | {target["FW"]:>9.1f} | {result["sim_fw"]:>8.2f} | '
                  f'{result["fw_err"]:>+7.1f}% | {target["Anth"]:>11.1f} | '
                  f'{result["sim_anth"]:>10.2f} | {result["anth_err"]:>+7.1f}% | {result["final_stress"]:>6.1f}')
        else:
            print(f'{treatment:<10} | FAILED')

    print('=' * 90)

    # 統計摘要 (全部 6 組)
    all_fw_errs = [abs(results[t]['fw_err']) for t in all_treatments if t in results]
    all_anth_errs = [abs(results[t]['anth_err']) for t in all_treatments if t in results]

    print("\nSummary Statistics (All 6 treatments):")
    print(f"    Mean |FW_Err|:   {np.mean(all_fw_errs):.1f}%")
    print(f"    Mean |Anth_Err|: {np.mean(all_anth_errs):.1f}%")
    print(f"    Max  |FW_Err|:   {np.max(all_fw_errs):.1f}%")
    print(f"    Max  |Anth_Err|: {np.max(all_anth_errs):.1f}%")
    print('=' * 90)

    return results


if __name__ == '__main__':
    run_validation()
