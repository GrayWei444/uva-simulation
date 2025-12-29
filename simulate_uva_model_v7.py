"""
萵苣生長與UVA效應整合模型
===========================
版本: v7.1 (純連續機制 - 數據驅動參數擬合)
日期: 2025-12-29

目標: 使用純連續冪函數實現所有機制，參數由數據自動擬合

核心設計哲學：
  「所有機制都是連續函數，參數由實驗數據驅動擬合」
  - 3h、6h、12h 只是連續曲線上的點
  - k_day 和 n_day 由優化器從數據中自動確定，非人工預設
  - 模型具有可預測性

v7.1 核心機制 (純連續 + 數據驅動):
1. 日內逆境: factor = 1 + k_day × hours^n_day
   - 參數由 calibrate_day_factor.py 自動擬合
   - n_day ≈ 6.6 是數據要求的，不是人為設定
   - 物理意義: 高冪次反映 ROS 引發的正反饋級聯效應
     - 蛋白質氧化 → 酶失活 → 清除能力下降
     - 脂質過氧化 → 膜損傷 → 離子洩漏
     - 多層效應累積產生高度非線性響應

2. 夜間節律: circ = 1 + k_night × hours_in_dark^n_night
   - k_night = 0.15, n_night = 2.0
   - 連續函數，取代原來的硬閾值

3. LAI 脆弱性: vuln = cap × (LAI_ref/LAI)^n / (cap + (LAI_ref/LAI)^n)
   - Sigmoid 型連續函數

4. Stress 非線性: nl = 1 + k × Stress / (K + Stress)
   - Michaelis-Menten 型連續函數

達標率: 12/12 (100%)
- FW: 6/6 <5%
- Anth: 6/6 <10%

修改紀錄:
- v6.9: 修復容量耗竭機制
- v7.0: 純連續機制
- v7.1: 數據驅動參數擬合 ⭐⭐⭐
        k_day, n_day 由優化器自動確定
        新增 calibrate_day_factor.py 校準腳本
        n_day = 6.6 (優化器找到的最佳值)

===============================================================================
"""

import numpy as np
from scipy.integrate import solve_ivp

# 導入基礎 Sun 模型
from lettuce_uva_carbon_complete_model import SunParams as BaseSunParams
from lettuce_uva_carbon_complete_model import sun_derivatives_final


# ==============================================================================
# 模型參數定義 (v7.1 - 數據驅動擬合版本)
# ==============================================================================

ALL_PARAMS = {
    # ──────────────────────────────────────────────────────────────────────────
    # 1. 基礎光合參數
    # ──────────────────────────────────────────────────────────────────────────
    'c_alpha': 0.555,                    # 光合效率 [-]

    # ──────────────────────────────────────────────────────────────────────────
    # 2. UVA-PAR 轉換參數
    # ──────────────────────────────────────────────────────────────────────────
    'par_conversion_factor': 1.0,        # UVA 轉 PAR 係數 [-]

    # ──────────────────────────────────────────────────────────────────────────
    # 3. Stress 損傷與修復參數
    # ──────────────────────────────────────────────────────────────────────────
    # 3.1 基礎損傷修復
    'stress_damage_coeff': 0.66e-6,      # 損傷係數 [1/(W/m²·s)]
    'stress_repair_coeff': 1.0e-5,       # 修復係數 [1/s]
    'stress_nonlinear_coeff': 8.0,       # Stress 累積非線性係數 [-]
    'K_nonlinear': 0.8,                  # Stress 非線性半飽和常數 [-]

    # 3.2 LAI 脆弱性（幼苗更易受損）- 連續 Sigmoid
    'LAI_ref_vuln': 6.5,                 # 脆弱性參考 LAI [m²/m²]
    'n_vuln': 8,                         # 脆弱性指數 [-]
    'cap_vuln': 100.0,                   # 脆弱性上限 [-]

    # 3.3 日內逆境 (v7.1: 數據驅動擬合)
    # 公式: day_factor = 1 + k_day × hours^n_day
    # 參數由 calibrate_day_factor.py 優化器自動確定，非人工預設
    # 物理意義: n ≈ 6.6 反映 ROS 引發的正反饋級聯效應
    #   - 蛋白質氧化 → 酶失活 → 清除能力下降
    #   - 脂質過氧化 → 膜損傷 → 離子洩漏
    #   - 多層效應累積產生高度非線性響應
    'k_day': 2.29e-5,                    # 日內逆境係數 [-] (數據驅動擬合)
    'n_day': 6.6,                        # 日內逆境冪次 [-] (數據驅動擬合)

    # 3.4 夜間節律 (v7.0 改為連續函數)
    # 公式: circ = 1 + k_night × hours_in_dark^n_night
    # 生理意義: 進入暗期越久，節律干擾越嚴重
    # 取代原來的硬閾值 `3.8 if is_night else 1.0`
    'k_night': 0.15,                     # 夜間節律係數 [-]
    'n_night': 2.0,                      # 夜間節律冪次 [-]

    # 3.5 Stress 對生長的抑制
    'stress_photosynthesis_inhibition': 0.66,  # 光合抑制係數 [-]
    'stress_lai_inhibition': 0.66,             # LAI 抑制係數 [-]
    'K_stress': 1.9,                           # 抑制半飽和常數 [-]

    # ──────────────────────────────────────────────────────────────────────────
    # 4. 碳修復參數
    # ──────────────────────────────────────────────────────────────────────────
    'base_repair_capacity': 0.5,         # 基礎修復能力 [-]
    'carbon_repair_bonus': 0.5,          # 碳池修復加成 [-]
    'K_carbon': 0.001,                   # 碳修復半飽和常數 [kg C/m²]
    'repair_carbon_cost': 1.0e-6,        # 修復碳成本 [kg C/Stress]

    # ──────────────────────────────────────────────────────────────────────────
    # 5. 花青素參數 (v7.0 優化後 - CK 和 L6D6 優先校準)
    # ──────────────────────────────────────────────────────────────────────────
    'base_anth_rate_light': 1.85e-10,    # 日間基礎合成率 [kg Anth/(kg FW·s)]
    'base_anth_rate_dark': 9.25e-11,     # 夜間基礎合成率 [kg Anth/(kg FW·s)]
    'V_max_anth': 3.50e-11,              # 最大誘導合成率 [kg Anth/(kg FW·s)]
    'K_stress_anth': 0.10,               # Stress 半飽和常數 [-]
    'k_deg': 3.02e-6,                    # 花青素降解率 [1/s]
    'anth_carbon_cost': 0.0,             # 合成碳成本 [kg C/kg Anth]

    # ──────────────────────────────────────────────────────────────────────────
    # 6. LDMC 參數
    # ──────────────────────────────────────────────────────────────────────────
    'dw_fw_ratio_base': 0.05,            # 基礎 DW:FW 比例 [-]
    'ldmc_stress_sensitivity': 1.0,      # LDMC 對 Stress 敏感度 [-]
    'K_ldmc': 50.0,                      # LDMC 半飽和 Stress 值 [-]
    'dw_fw_ratio_max': 0.12,             # 最大 DW:FW 比例 [-]
}


# ==============================================================================
# 模型參數類別
# ==============================================================================

class UVAParams(BaseSunParams):
    """
    UVA 效應模型參數類別 (v7.0 純連續版本)

    所有機制都是連續可微函數，無任何硬閾值
    狀態變量: [X_d, C_buf, LAI, Anth, Stress]
    """

    def __init__(self, params=None):
        super().__init__()

        if params is None:
            params = ALL_PARAMS

        # 1. 基礎光合參數
        self.c_alpha = params['c_alpha']

        # 2. UVA-PAR 轉換參數
        self.par_conversion_factor = params['par_conversion_factor']

        # 3. Stress 損傷與修復參數
        self.stress_damage_coeff = params['stress_damage_coeff']
        self.stress_repair_coeff = params['stress_repair_coeff']
        self.stress_nonlinear_coeff = params['stress_nonlinear_coeff']
        self.K_nonlinear = params['K_nonlinear']
        self.LAI_ref_vuln = params['LAI_ref_vuln']
        self.n_vuln = params['n_vuln']
        self.cap_vuln = params['cap_vuln']

        # 日內逆境 (純冪函數)
        self.k_day = params['k_day']
        self.n_day = params['n_day']

        # 夜間節律 (連續冪函數)
        self.k_night = params['k_night']
        self.n_night = params['n_night']

        self.stress_photosynthesis_inhibition = params['stress_photosynthesis_inhibition']
        self.stress_lai_inhibition = params['stress_lai_inhibition']
        self.K_stress = params['K_stress']

        # 4. 碳修復參數
        self.base_repair_capacity = params['base_repair_capacity']
        self.carbon_repair_bonus = params['carbon_repair_bonus']
        self.K_carbon = params['K_carbon']
        self.repair_carbon_cost = params['repair_carbon_cost']

        # 5. 花青素參數
        self.base_anth_rate_light = params['base_anth_rate_light']
        self.base_anth_rate_dark = params['base_anth_rate_dark']
        self.V_max_anth = params['V_max_anth']
        self.K_stress_anth = params['K_stress_anth']
        self.k_deg = params['k_deg']
        self.anth_carbon_cost = params['anth_carbon_cost']

        # 6. LDMC 參數
        self.dw_fw_ratio_base = params['dw_fw_ratio_base']
        self.ldmc_stress_sensitivity = params['ldmc_stress_sensitivity']
        self.K_ldmc = params['K_ldmc']
        self.dw_fw_ratio_max = params['dw_fw_ratio_max']


# ==============================================================================
# 輔助函數
# ==============================================================================

def calculate_dynamic_dw_fw_ratio(Stress, p):
    """
    根據 Stress 計算動態 DW:FW 比例（LDMC 效應）
    連續 Michaelis-Menten 型函數
    """
    stress_effect = p.ldmc_stress_sensitivity * Stress / (p.K_ldmc + Stress + 1e-9)
    ratio = p.dw_fw_ratio_base * (1.0 + stress_effect)
    return min(ratio, p.dw_fw_ratio_max)


# ==============================================================================
# 核心微分方程 (v7.0 純連續版本)
# ==============================================================================

def uva_sun_derivatives(t, state, p, env):
    """
    UVA效應整合模型的微分方程 (v7.0 純連續版本)

    v7.0 核心改進:
    - 日內逆境: 純冪函數 1 + k × hours^n
    - 夜間節律: 連續函數 1 + k × hours_in_dark^n
    - 無任何硬閾值或參考點

    狀態變量 (5個):
    [X_d, C_buf, LAI, Anth, Stress]
    """

    # =========================================================================
    # 步驟1: 解包狀態變量
    # =========================================================================
    X_d, C_buf, LAI, Anth, Stress = state

    X_d = max(X_d, 1e-9)
    C_buf = max(C_buf, 0)
    LAI = max(LAI, 0.1)
    Anth = max(Anth, 0)
    Stress = max(Stress, 0)

    # =========================================================================
    # 步驟2: 計算時間相關變量
    # =========================================================================
    hour = (t / 3600.0) % 24.0
    day_from_sowing = t / 86400.0

    # =========================================================================
    # 步驟3: 判斷日夜狀態
    # =========================================================================
    light_on = env['light_on_hour']
    light_off = env['light_off_hour']

    if light_on <= light_off:
        is_day = light_on <= hour < light_off
    else:
        is_day = hour >= light_on or hour < light_off

    if is_day:
        I_base = env['I_day']
        Tc = env['T_day']
    else:
        I_base = 0.0
        Tc = env['T_night']

    # =========================================================================
    # 步驟4: 計算 UVA 強度和當天照射時間
    # =========================================================================
    uva_on = env.get('uva_on', False)
    uva_start_day = env.get('uva_start_day', 29)
    uva_end_day = env.get('uva_end_day', 35)
    uva_hour_on = env.get('uva_hour_on', 10)
    uva_hour_off = env.get('uva_hour_off', 16)
    uva_intensity = env.get('uva_intensity', 22.0)

    I_UVA = 0.0
    in_uva_window = False
    hours_today = 0.0

    if uva_on:
        if uva_hour_on <= uva_hour_off:
            if uva_start_day <= day_from_sowing <= uva_end_day:
                if uva_hour_on <= hour < uva_hour_off:
                    I_UVA = uva_intensity
                    in_uva_window = True
                    hours_today = hour - uva_hour_on
        else:
            if hour >= uva_hour_on:
                if uva_start_day <= day_from_sowing <= uva_end_day:
                    I_UVA = uva_intensity
                    in_uva_window = True
                    hours_today = hour - uva_hour_on
            elif hour < uva_hour_off:
                day_session_started = day_from_sowing - 1
                if uva_start_day <= day_session_started <= uva_end_day:
                    I_UVA = uva_intensity
                    in_uva_window = True
                    hours_today = (24 - uva_hour_on) + hour

    # =========================================================================
    # 步驟5: 連續機制 1 - 日內逆境 (純冪函數)
    # =========================================================================
    # factor = 1 + k × hours^n
    # 這是完全連續的函數，沒有任何參考點或閾值
    # hours=0: factor=1, hours=6: factor=3.8, hours=12: factor=359
    day_factor = 1.0 + p.k_day * (hours_today ** p.n_day)

    # =========================================================================
    # 步驟6: 連續機制 2 - 夜間節律 (連續冪函數)
    # =========================================================================
    # 計算進入暗期多久
    if light_on <= light_off:
        if light_on <= hour < light_off:
            hours_in_dark = 0.0
        else:
            if hour >= light_off:
                hours_in_dark = hour - light_off
            else:
                hours_in_dark = hour + (24 - light_off)
    else:
        if hour >= light_on or hour < light_off:
            hours_in_dark = 0.0
        else:
            hours_in_dark = hour - light_off if hour >= light_off else hour + (24 - light_off)

    # 夜間 UVA 時使用連續節律因子
    if I_UVA > 0 and hours_in_dark > 0:
        circadian_factor = 1.0 + p.k_night * (hours_in_dark ** p.n_night)
    else:
        circadian_factor = 1.0

    # =========================================================================
    # 步驟7: 呼叫基礎 Sun 模型
    # =========================================================================
    I_effective = I_base + I_UVA * p.par_conversion_factor

    env_modified = env.copy()
    env_modified['I_override'] = I_effective
    env_modified['T_override'] = Tc
    env_modified['is_day_override'] = is_day

    base_state = [X_d, C_buf, LAI]
    dXd_dt_base, dCbuf_dt, dLAI_dt_base = sun_derivatives_final(t, base_state, p, env_modified)

    # =========================================================================
    # 步驟8: 連續機制 3 - LAI 脆弱性 (Sigmoid 型)
    # =========================================================================
    base_vuln = (p.LAI_ref_vuln / LAI) ** p.n_vuln
    vulnerability = p.cap_vuln * base_vuln / (p.cap_vuln + base_vuln)

    # =========================================================================
    # 步驟9: 連續機制 4 - Stress 非線性 (Michaelis-Menten 型)
    # =========================================================================
    nonlinear_factor = 1.0 + p.stress_nonlinear_coeff * Stress / (p.K_nonlinear + Stress + 1e-9)

    # =========================================================================
    # 步驟10: 損傷率 (所有連續因子相乘)
    # =========================================================================
    damage_rate = (p.stress_damage_coeff * I_UVA * vulnerability *
                   day_factor * nonlinear_factor * circadian_factor)

    # =========================================================================
    # 步驟11: 修復率
    # =========================================================================
    C_buf_positive = max(C_buf, 0)
    repair_capacity = (p.base_repair_capacity +
                       p.carbon_repair_bonus * C_buf_positive / (p.K_carbon + C_buf_positive + 1e-9))
    repair_rate = p.stress_repair_coeff * Stress * repair_capacity

    dStress_dt = damage_rate - repair_rate

    # =========================================================================
    # 步驟12: Stress 對生長的抑制
    # =========================================================================
    stress_inhibition = Stress / (p.K_stress + Stress + 1e-9)
    xd_reduction = p.stress_photosynthesis_inhibition * stress_inhibition
    lai_reduction = p.stress_lai_inhibition * stress_inhibition

    dXd_dt = dXd_dt_base * (1.0 - xd_reduction) if dXd_dt_base > 0 else dXd_dt_base
    dLAI_dt = dLAI_dt_base * (1.0 - lai_reduction) if dLAI_dt_base > 0 else dLAI_dt_base

    # =========================================================================
    # 步驟13: 花青素合成
    # =========================================================================
    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress, p)
    FW_kg_m2 = X_d / dw_fw_ratio

    base_synthesis = p.base_anth_rate_light if is_day else p.base_anth_rate_dark
    stress_induced = p.V_max_anth * Stress / (p.K_stress_anth + Stress + 1e-12)
    synthesis_rate = FW_kg_m2 * (base_synthesis + stress_induced)

    degradation = p.k_deg * Anth
    dAnth_dt = synthesis_rate - degradation

    # =========================================================================
    # 步驟14: 碳消耗
    # =========================================================================
    repair_carbon_consumption = repair_rate * p.repair_carbon_cost
    anth_carbon_consumption = synthesis_rate * p.anth_carbon_cost
    dCbuf_dt = dCbuf_dt - repair_carbon_consumption - anth_carbon_consumption

    return np.array([dXd_dt, dCbuf_dt, dLAI_dt, dAnth_dt, dStress_dt])


# ==============================================================================
# 主程式
# ==============================================================================

if __name__ == "__main__":
    from model_config import (
        ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment
    )

    p = UVAParams()

    # 計算 day_factor 數值
    factor_6h = 1 + p.k_day * (6 ** p.n_day)
    factor_12h = 1 + p.k_day * (12 ** p.n_day)

    print("=" * 80)
    print("萵苣生長與UVA效應整合模型 v7.1 (數據驅動參數擬合)")
    print("=" * 80)
    print("\n核心機制 (全部連續，參數由數據擬合):")
    print(f"  1. 日內逆境: factor = 1 + {p.k_day:.2e} × hours^{p.n_day:.1f}")
    print(f"     - 6h:  factor = {factor_6h:.1f}")
    print(f"     - 12h: factor = {factor_12h:.1f}")
    print(f"     - 12h/6h 比值 = {factor_12h/factor_6h:.1f}")
    print(f"     - n = {p.n_day:.1f} 由優化器從數據自動擬合")
    print("  2. 夜間節律: circ = 1 + 0.15 × hours_in_dark^2")
    print("  3. LAI 脆弱性: sigmoid 型連續函數")
    print("  4. Stress 非線性: Michaelis-Menten 型連續函數")
    print("=" * 80)

    # 測試所有處理組
    fw_errs = []
    anth_errs = []

    for treatment in ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']:
        env = get_env_for_treatment(treatment)
        target = TARGETS.get(treatment, {'FW': 0, 'Anth': 0})

        # 初始條件
        fw_init_g = SIMULATION['initial_fw_g']
        dw_init_g = fw_init_g * p.dw_fw_ratio_base
        Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
        C_buf_init = Xd_init * 0.1
        LAI_init = (dw_init_g / 0.01) * 0.04
        fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
        Anth_init = 5.0 * fw_total_init / 1e6
        initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0]

        # 模擬時間
        transplant_day = SIMULATION['transplant_offset']
        simulation_days = SIMULATION['days']
        t_start = transplant_day * 86400
        t_end = (transplant_day + simulation_days) * 86400

        # 求解
        sol = solve_ivp(
            uva_sun_derivatives,
            (t_start, t_end),
            initial_state,
            args=(p, env),
            method='RK45',
            max_step=300,
            t_eval=np.array([t_end])
        )

        if sol.success:
            Xd_f, Cbuf_f, LAI_f, Anth_f, Stress_f = sol.y[:, -1]
            dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress_f, p)
            FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
            FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
            Anth_sim = Anth_f / FW_total_kg * 1e6

            FW_obs = target['FW']
            Anth_obs = target['Anth']

            fw_err = (FW_sim - FW_obs) / FW_obs * 100
            anth_err = (Anth_sim - Anth_obs) / Anth_obs * 100

            fw_errs.append(abs(fw_err))
            anth_errs.append(abs(anth_err))

            s1 = "✓" if abs(fw_err) < 5 else "✗"
            s2 = "✓" if abs(anth_err) < 10 else "✗"
            print(f"{treatment:<8} FW:{FW_sim:>5.1f}g({fw_err:>+5.1f}%{s1}) "
                  f"Anth:{Anth_sim:>5.1f}({anth_err:>+5.1f}%{s2}) "
                  f"LAI:{LAI_f:>4.1f} S:{Stress_f:>5.1f}")
        else:
            print(f"{treatment:<8} 模擬失敗: {sol.message}")

    print("-" * 80)
    fw_ok = sum(1 for e in fw_errs if e < 5)
    anth_ok = sum(1 for e in anth_errs if e < 10)
    print(f"達標: FW {fw_ok}/6, Anth {anth_ok}/6, Total {fw_ok + anth_ok}/12")
