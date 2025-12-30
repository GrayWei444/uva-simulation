"""
萵苣生長與UVA效應整合模型
===========================
版本: v7.2 (機制修正版)
日期: 2025-12-30

目標: 修正機制設計，使生物學意義更合理

v7.2 核心修正：
1. **修復飽和因子 (repair_saturation)**：移到修復項
   - 機制：當日照射時間越長，修復能力越飽和
   - 公式：repair_factor = 1 / (1 + k_day × hours_today^n_day)
   - 生物學：照射消耗修復資源（酶、抗氧化物），導致修復能力下降

2. **夜間節律損傷 (circadian_damage)**：改為加項
   - 機制：只有夜間照射才有額外損傷，日間照射不受影響
   - 公式：damage = base_damage + circadian_damage（僅夜間 UVA）
   - 生物學：夜間 UVA 打斷生理節律，造成額外傷害

3. **ROS 正反饋 (nonlinear_factor)**：保留在損傷項
   - 機制：已累積的 Stress 會加速新損傷
   - 公式：nonlinear_factor = 1 + α × Stress / (K + Stress)
   - 生物學：ROS 級聯效應——傷害加速新傷害

Stress 方程：
  dStress/dt = (base_damage × nonlinear_factor + circadian_damage)
               - (repair_rate × repair_factor)

修改紀錄:
- v7.1: 數據驅動參數擬合
- v7.2: 機制修正 ⭐⭐⭐
        - intraday_factor 移到修復項（修復飽和）
        - circadian_penalty 改為加項（僅夜間）

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
    # 3. Stress 損傷參數
    # ──────────────────────────────────────────────────────────────────────────
    # 3.1 基礎損傷
    'stress_damage_coeff': 0.66e-6,      # 基礎損傷係數 k_damage [1/(W/m²·s)]

    # 3.2 ROS 正反饋 (nonlinear_factor)
    # 公式: nonlinear_factor = 1 + α × Stress / (K + Stress)
    # 文獻: Foyer & Noctor (2005) - ROS cascade amplification
    'stress_nonlinear_coeff': 8.0,       # ROS 正反饋係數 α [-]
    'K_nonlinear': 0.8,                  # 半飽和常數 K [-]

    # 3.3 LAI 脆弱性 (vulnerability)
    # 公式: vuln = cap × (LAI_ref/LAI)^n / (cap + (LAI_ref/LAI)^n)
    # 生物學: 幼苗防禦機制不成熟，更易受損
    'LAI_ref_vuln': 6.5,                 # 參考 LAI [m²/m²]
    'n_vuln': 8,                         # 脆弱性指數 [-]
    'cap_vuln': 100.0,                   # 脆弱性上限 [-]

    # 3.4 花青素保護 (anth_protection)
    # 公式: protection = α_anth × Anth / (K_anth_prot + Anth)
    # 文獻: Xu & Rothstein (2018) DOI:10.1080/15592324.2018.1451708
    #       "anthocyanins protect plants by scavenging excess ROS"
    # 文獻: Cerqueira et al. (2023) DOI:10.1007/s11103-023-01362-4
    #       "anthocyanins act as efficient antioxidants"
    'alpha_anth_protection': 0.5,        # 最大保護效率 (0~1) [-]
    'K_anth_protection': 5e-6,           # 保護半飽和常數 [kg Anth/m²]

    # 3.5 夜間節律損傷 (circadian_damage) - 加項
    # 公式: circ_damage = k_circ × I_UVA × hours_in_dark^n_circ（僅夜間 UVA）
    # 生物學: 夜間 UVA 打斷生理節律，造成額外損傷
    'k_circadian': 1e-7,                 # 夜間節律損傷係數 [-]
    'n_circadian': 2.0,                  # 夜間節律損傷冪次 [-]

    # ──────────────────────────────────────────────────────────────────────────
    # 4. Stress 修復參數
    # ──────────────────────────────────────────────────────────────────────────
    # 4.1 LAI 依賴修復 (repair ∝ LAI)
    # 公式: repair = k_repair × LAI × repair_saturation_factor
    # 文獻: Puthiyaveetil et al. (2014) DOI:10.1093/pcp/pcu062
    #       "decreased biomass and altered PSII functionality"
    # 文獻: Järvi et al. (2015) DOI:10.1016/j.bbabio.2015.01.006
    #       "a high number of auxiliary proteins assist D1 turnover"
    'k_repair_lai': 1e-6,                # LAI 修復係數 [1/(m²/m²·s)]

    # 4.2 修復飽和因子 (照射時間越長，修復能力越飽和)
    # 公式: sat_factor = 1 / (1 + k_sat × hours_today^n_sat)
    # 生物學: 照射消耗修復資源（酶、抗氧化物）
    'k_repair_sat': 1e-4,                # 修復飽和係數 [-]
    'n_repair_sat': 4.0,                 # 修復飽和冪次 [-]

    # ──────────────────────────────────────────────────────────────────────────
    # 5. 花青素消耗參數
    # ──────────────────────────────────────────────────────────────────────────
    # 公式: consumption = k_consume × Anth × base_damage
    # 文獻: Enaru et al. (2021) DOI:10.3390/antiox10121967
    #       "anthocyanins undergo oxidative degradation during ROS scavenging"
    'k_anth_consumption': 1e-3,          # 花青素消耗係數 [1/damage]

    # ──────────────────────────────────────────────────────────────────────────
    # 6. Stress 對生長的抑制
    # ──────────────────────────────────────────────────────────────────────────
    'stress_photosynthesis_inhibition': 0.66,  # 光合抑制係數 [-]
    'stress_lai_inhibition': 0.66,             # LAI 抑制係數 [-]
    'K_stress': 1.9,                           # 抑制半飽和常數 [-]

    # ──────────────────────────────────────────────────────────────────────────
    # 7. 花青素合成參數
    # ──────────────────────────────────────────────────────────────────────────
    'base_anth_rate_light': 1.85e-10,    # 日間基礎合成率 [kg Anth/(kg FW·s)]
    'base_anth_rate_dark': 9.25e-11,     # 夜間基礎合成率 [kg Anth/(kg FW·s)]
    'V_max_anth': 3.50e-11,              # 最大誘導合成率 [kg Anth/(kg FW·s)]
    'K_stress_anth': 0.10,               # Stress 半飽和常數 [-]
    'k_deg': 3.02e-6,                    # 花青素自然降解率 [1/s]
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

        # 3. Stress 損傷參數
        self.stress_damage_coeff = params['stress_damage_coeff']
        self.stress_nonlinear_coeff = params['stress_nonlinear_coeff']
        self.K_nonlinear = params['K_nonlinear']
        self.LAI_ref_vuln = params['LAI_ref_vuln']
        self.n_vuln = params['n_vuln']
        self.cap_vuln = params['cap_vuln']

        # 花青素保護 (v7.2)
        self.alpha_anth_protection = params['alpha_anth_protection']
        self.K_anth_protection = params['K_anth_protection']

        # 夜間節律損傷 (v7.2: 加項)
        self.k_circadian = params['k_circadian']
        self.n_circadian = params['n_circadian']

        # 4. Stress 修復參數 (v7.2: LAI 依賴)
        self.k_repair_lai = params['k_repair_lai']
        self.k_repair_sat = params['k_repair_sat']
        self.n_repair_sat = params['n_repair_sat']

        # 5. 花青素消耗
        self.k_anth_consumption = params['k_anth_consumption']

        # 6. Stress 對生長的抑制
        self.stress_photosynthesis_inhibition = params['stress_photosynthesis_inhibition']
        self.stress_lai_inhibition = params['stress_lai_inhibition']
        self.K_stress = params['K_stress']

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
    # 步驟5: 修復飽和因子 (v7.2: 移到修復項使用)
    # =========================================================================
    # repair_factor = 1 / (1 + k × hours^n)
    # 照射時間越長，修復能力越飽和（越接近 0）
    # - hours=0: repair_factor = 1.0（正常修復）
    # - hours=6: repair_factor ≈ 0.7（部分飽和）
    # - hours=12: repair_factor ≈ 0.3（嚴重飽和）
    repair_saturation_factor = 1.0 / (1.0 + p.k_repair_sat * (hours_today ** p.n_repair_sat))

    # =========================================================================
    # 步驟6: 夜間節律損傷 (v7.2: 改為加項)
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

    # 夜間 UVA 時才有額外節律損傷（加項，非乘項）
    # - 日間照射: circadian_damage = 0
    # - 夜間照射: circadian_damage = k × I_UVA × hours_in_dark^n
    if I_UVA > 0 and hours_in_dark > 0:
        circadian_damage = p.k_circadian * I_UVA * (hours_in_dark ** p.n_circadian)
    else:
        circadian_damage = 0.0

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
    # 步驟10: 花青素保護因子 (v7.2)
    # =========================================================================
    # 公式: protection = α × Anth / (K + Anth)
    # 文獻: Xu & Rothstein (2018) - anthocyanins scavenge excess ROS
    # 文獻: Cerqueira et al. (2023) - anthocyanins act as efficient antioxidants
    anth_protection = p.alpha_anth_protection * Anth / (p.K_anth_protection + Anth + 1e-12)

    # =========================================================================
    # 步驟11: 損傷率 (v7.2: 花青素減少損傷 + 夜間節律加項)
    # =========================================================================
    # 基礎損傷：k_damage × I_UVA × vulnerability × nonlinear_factor × (1 - protection)
    # - vulnerability: LAI 脆弱性（幼苗更易受損）
    # - nonlinear_factor: ROS 正反饋（Stress 加速新損傷）
    # - (1 - protection): 花青素保護效應（減少損傷）
    base_damage = p.stress_damage_coeff * I_UVA * vulnerability * nonlinear_factor
    protected_damage = base_damage * (1.0 - anth_protection)

    # 總損傷 = 受保護損傷 + 夜間節律損傷（加項）
    damage_rate = protected_damage + circadian_damage

    # =========================================================================
    # 步驟12: 修復率 (v7.2: LAI 依賴 + 修復飽和)
    # =========================================================================
    # 公式: repair = k_repair × LAI × repair_saturation_factor
    # 文獻: Puthiyaveetil et al. (2014) - biomass correlates with repair capacity
    # 文獻: Järvi et al. (2015) - auxiliary proteins required for D1 turnover
    # - LAI: 植株越大，修復資源越多
    # - repair_saturation_factor: 照射越久，修復能力越飽和
    repair_rate = p.k_repair_lai * LAI * repair_saturation_factor

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
    # 步驟14: 花青素動態 (v7.2: 合成 - 自然降解 - ROS消耗)
    # =========================================================================
    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress, p)
    FW_kg_m2 = X_d / dw_fw_ratio

    # 合成：基礎合成 + Stress誘導合成
    base_synthesis = p.base_anth_rate_light if is_day else p.base_anth_rate_dark
    stress_induced = p.V_max_anth * Stress / (p.K_stress_anth + Stress + 1e-12)
    synthesis_rate = FW_kg_m2 * (base_synthesis + stress_induced)

    # 自然降解
    natural_degradation = p.k_deg * Anth

    # ROS 消耗：花青素清除 ROS 時被氧化降解
    # 公式: consumption = k_consume × Anth × base_damage
    # 文獻: Enaru et al. (2021) DOI:10.3390/antiox10121967
    #       "anthocyanins undergo oxidative degradation during ROS scavenging"
    ros_consumption = p.k_anth_consumption * Anth * base_damage

    dAnth_dt = synthesis_rate - natural_degradation - ros_consumption

    # =========================================================================
    # 步驟15: 碳消耗
    # =========================================================================
    anth_carbon_consumption = synthesis_rate * p.anth_carbon_cost
    dCbuf_dt = dCbuf_dt - anth_carbon_consumption

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
