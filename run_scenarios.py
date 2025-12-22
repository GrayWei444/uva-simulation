# -*- coding: utf-8 -*-
"""
================================================================================
夜間抑制參數場景模擬 (Night Inhibition Scenario Analysis)
================================================================================
版本: v2.0
日期: 2025-12-05

功能說明:
---------
對 L6D6 (日間) 和 L6D6-N (夜間) 處理組進行三種情境分析，
以量化夜間抑制參數的不確定性對模型預測的影響。

場景設定:
---------
  - Low: 低抑制 (快速恢復) - 半衰期約6小時
  - Medium: 中等抑制 (文獻中值) - 半衰期約12小時
  - High: 高抑制 (慢速恢復) - 半衰期約24小時

文獻依據:
---------
  - Bennie et al. (2016) J. Ecol. DOI:10.1111/1365-2745.12551
  - Deng et al. (2025) Biology DOI:10.3390/biology14050571
  - Harmer (2009) Annu. Rev. Plant Biol.
  - Covington et al. (2008) Genome Biol.

輸出用途:
---------
  - 論文補充資料中的不確定性帶
  - 驗證觀測值是否落在模型預測範圍內

用法:
-----
    python run_scenarios.py

依賴模組:
---------
  - model_config: 共用設定 (環境參數、場景參數)
  - simulate_uva_model: UVA效應模型 (參數、微分方程)
================================================================================
"""

import numpy as np
from scipy.integrate import solve_ivp

# 匯入共用設定
from model_config import (
    ENV_BASE, TREATMENT_CONFIGS, TARGETS, ODE_SETTINGS, SIMULATION,
    NIGHT_INHIBITION_SCENARIOS, get_env_for_treatment
)

# 匯入模型
from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio


class ScenarioParams(UVAParams):
    """
    場景模擬專用參數類別

    繼承已校準的 UVAParams，只覆蓋夜間抑制參數以進行敏感度分析。
    其他所有參數保持已校準的值不變。

    參數:
    -----
    scenario : str
        場景名稱 ('Low', 'Medium', 'High')

    覆蓋的參數:
    -----------
    - night_uva_base_inhibition: 夜間UVA對基礎花青素合成的抑制強度
    - circadian_inhibition_decay: 抑制效應的衰減速率
    """

    def __init__(self, scenario='Medium'):
        super().__init__()

        # 只覆蓋夜間抑制參數
        scenario_config = NIGHT_INHIBITION_SCENARIOS.get(
            scenario,
            NIGHT_INHIBITION_SCENARIOS['Medium']
        )
        self.night_uva_base_inhibition = scenario_config['night_uva_base_inhibition']
        self.circadian_inhibition_decay = scenario_config['circadian_inhibition_decay']


def calculate_initial_conditions(p):
    """
    計算模擬的初始條件

    參數:
    -----
    p : UVAParams or ScenarioParams
        模型參數物件

    回傳:
    -----
    list
        8個狀態變量的初始值
    """
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    LAI_init = (dw_init_g / 0.01) * 0.04
    Anth_init = 2.0 * (fw_init_g * ENV_BASE['plant_density'] / 1000) / 1000

    return [Xd_init, 0.0, LAI_init, Anth_init, 0.0, 0.0, 0.0, 0.0]


def run_simulation(treatment, scenario='Medium'):
    """
    運行單一處理組的場景模擬

    使用已校準的 uva_sun_derivatives 微分方程，只覆蓋夜間抑制參數。

    參數:
    -----
    treatment : str
        處理組代碼 ('L6D6' 或 'L6D6-N')
    scenario : str
        場景名稱 ('Low', 'Medium', 'High')

    回傳:
    -----
    tuple (float, float) or (None, None)
        (鮮重 g/plant, 花青素 ppm)
        模擬失敗時回傳 (None, None)
    """
    p = ScenarioParams(scenario)
    env = get_env_for_treatment(treatment)
    initial_state = calculate_initial_conditions(p)

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
        return None, None

    # 計算鮮重
    sim_D_UVA = sol.y[4, -1]
    sim_dw_fw_ratio = calculate_dynamic_dw_fw_ratio(sim_D_UVA, p)
    sim_dw_per_plant = sol.y[0, -1] / ENV_BASE['plant_density'] * 1000
    sim_fw = sim_dw_per_plant / sim_dw_fw_ratio

    # 計算花青素濃度 (ppm)
    sim_fw_kg_per_m2 = sim_fw * ENV_BASE['plant_density'] / 1000
    sim_anth_mg_per_m2 = sol.y[3, -1] * 1000
    sim_anth_ppm = sim_anth_mg_per_m2 / (sim_fw_kg_per_m2 + 1e-9)

    return sim_fw, sim_anth_ppm


def run_all_scenarios():
    """
    對 L6D6 和 L6D6-N 執行所有場景模擬

    執行流程:
    ---------
    1. 輸出場景參數設定
    2. 對每個處理組和場景組合執行模擬
    3. 計算並輸出誤差
    4. 輸出不確定性範圍 (供論文使用)

    回傳:
    -----
    dict
        巢狀字典結構: results[treatment][scenario] = {FW, Anth, FW_err, Anth_err}
    """
    treatments = ['L6D6', 'L6D6-N']
    scenarios = ['Low', 'Medium', 'High']

    # 標頭
    print("=" * 80)
    print("夜間抑制參數場景分析 (Night Inhibition Scenario Analysis)")
    print("=" * 80)

    # 場景參數設定
    print("\n場景參數設定:")
    print("-" * 70)
    print(f"{'場景':<10} {'night_uva_base_inhibition':<25} {'circadian_inhibition_decay':<25}")
    print("-" * 70)
    for scenario, config in NIGHT_INHIBITION_SCENARIOS.items():
        print(f"{scenario:<10} {config['night_uva_base_inhibition']:<25.2f} "
              f"{config['circadian_inhibition_decay']:<25.2e}")
    print("-" * 70)

    # 模擬結果
    results = {}
    print("\n模擬結果:")
    print("=" * 80)
    print(f"{'處理組':<10} {'場景':<10} {'FW (g)':<12} {'FW誤差%':<12} "
          f"{'Anth (ppm)':<12} {'Anth誤差%':<12}")
    print("=" * 80)

    for treatment in treatments:
        target = TARGETS[treatment]
        results[treatment] = {}

        for scenario in scenarios:
            fw, anth = run_simulation(treatment, scenario)

            if fw is None:
                print(f"{treatment:<10} {scenario:<10} {'模擬失敗':<12}")
                continue

            # 計算誤差
            fw_err = (fw - target['FW']) / target['FW'] * 100
            anth_err = (anth - target['Anth']) / target['Anth'] * 100

            results[treatment][scenario] = {
                'FW': fw,
                'Anth': anth,
                'FW_err': fw_err,
                'Anth_err': anth_err
            }

            print(f"{treatment:<10} {scenario:<10} {fw:<12.2f} {fw_err:<+12.1f} "
                  f"{anth:<12.2f} {anth_err:<+12.1f}")

    print("=" * 80)

    # 不確定性範圍 (供論文使用)
    print("\n不確定性範圍 (論文用):")
    print("-" * 70)

    for treatment in treatments:
        if treatment not in results:
            continue

        target = TARGETS[treatment]
        fw_range = [results[treatment][s]['FW'] for s in scenarios if s in results[treatment]]
        anth_range = [results[treatment][s]['Anth'] for s in scenarios if s in results[treatment]]
        fw_err_range = [results[treatment][s]['FW_err'] for s in scenarios if s in results[treatment]]
        anth_err_range = [results[treatment][s]['Anth_err'] for s in scenarios if s in results[treatment]]

        if fw_range and anth_range:
            print(f"\n{treatment} (目標: FW={target['FW']}g, Anth={target['Anth']}ppm):")
            print(f"  FW範圍: {min(fw_range):.2f} - {max(fw_range):.2f} g "
                  f"(誤差: {min(fw_err_range):+.1f}% 至 {max(fw_err_range):+.1f}%)")
            print(f"  Anth範圍: {min(anth_range):.2f} - {max(anth_range):.2f} ppm "
                  f"(誤差: {min(anth_err_range):+.1f}% 至 {max(anth_err_range):+.1f}%)")

    print("\n" + "=" * 80)

    return results


if __name__ == '__main__':
    results = run_all_scenarios()
