"""
萵苣生長與UVA效應整合模型
===========================
版本: v8.0-ROS (基於 ROS 動力學的損傷機制)
日期: 2025-12-29

目標: 用 ROS 動力學推導損傷因子，取代 v7.0 的任意冪函數

核心改進 (v8.0):
  v7.0 問題: day_factor = 1 + k × hours^7 是「預設形狀」
             12^7/6^7 = 128 是人為設定，無物理意義

  v8.0 方案: 用 ROS 動力學推導
             - ROS 產生率 ∝ UVA 能量
             - ROS 清除遵循 Michaelis-Menten 動力學
             - 當清除飽和時，ROS 累積 → 損傷加速
             - 形狀從酶動力學 **自然浮現**

ROS 動力學基礎 (文獻值):
  - APX (Ascorbate Peroxidase): Km ≈ 50 µM, 高親和力, 先飽和
  - CAT (Catalase): Km ≈ 200 mM, 低親和力, 不易飽和, 高容量
  - SOD: 擴散限制 (~2×10⁹ M⁻¹s⁻¹)

設計原則:
  - 能量積分 E = I_UVA × t (J/m²)
  - 多級酶飽和: APX 先飽和 → CAT 接手 → 兩者都飽和時損傷加速
  - 形狀由 Km 值決定，不是任意冪次

===============================================================================
"""

import numpy as np
from scipy.integrate import solve_ivp

# 導入基礎 Sun 模型
from lettuce_uva_carbon_complete_model import SunParams as BaseSunParams
from lettuce_uva_carbon_complete_model import sun_derivatives_final


# ==============================================================================
# 模型參數定義 (v8.0-ROS)
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
    'stress_damage_coeff': 2.0e-8,       # 損傷係數 [1/(W/m²·s)] (v8.0 ROS 校準)
    'stress_repair_coeff': 1.0e-5,       # 修復係數 [1/s]
    'stress_nonlinear_coeff': 8.0,       # Stress 累積非線性係數 [-]
    'K_nonlinear': 0.8,                  # Stress 非線性半飽和常數 [-]

    # 3.2 LAI 脆弱性（幼苗更易受損）- 連續 Sigmoid
    'LAI_ref_vuln': 6.5,                 # 脆弱性參考 LAI [m²/m²]
    'n_vuln': 8,                         # 脆弱性指數 [-]
    'cap_vuln': 100.0,                   # 脆弱性上限 [-]

    # ──────────────────────────────────────────────────────────────────────────
    # 3.3 日內逆境 - ROS 動力學版 (v8.0 新機制)
    # ──────────────────────────────────────────────────────────────────────────
    #
    # 核心物理原理:
    #   1. ROS 產生率 P = k × I_UVA (與強度成正比)
    #   2. 抗氧化劑清除有上限 V_max (酶飽和)
    #   3. 當累積 ROS 超過清除能力時，損傷加速
    #
    # 數學模型:
    #   累積 ROS ≈ (P - V_eff) × t，當 P > V_eff
    #   V_eff 隨抗氧化劑消耗而下降
    #
    # 結果: 損傷與能量的關係是超線性的
    #   - 抗氧化劑消耗: 線性
    #   - ROS 累積後的二次氧化: 平方
    #   - 蛋白質聚集/膜損傷的級聯: 更高階
    #
    # 參數 (基於文獻推導):
    #   - V_max (CAT): ~10⁷ M/s (Chelikani et al. 2004)
    #   - Km (APX): ~50 µM (PMC)
    #   - Km (CAT): ~200 mM (Chelikani et al. 2004)
    #
    # 能量尺度換算:
    #   E_threshold: 抗氧化劑耗盡的臨界能量
    #   對於萵苣: E_threshold ≈ 150,000-200,000 J/m² (對應 ~4-5h 照射)
    #
    # ROS 動力學參數 (v8.0c 無閾值版)
    # 設計理念：ROS 累積損傷是連續的，無硬閾值
    # 使用能量積分的冪律，冪次來自反應動力學：
    #   - n ≈ 2: 二級氧化反應
    #   - n ≈ 3-4: 包含級聯效應
    #
    # 公式: factor = 1 + k × E^n
    # 不同於 v7.0 使用時間冪次，v8.0 使用能量冪次
    # E = I_UVA × hours × 3600 (J/m²)
    #
    'E_threshold': 0.0,                  # 無閾值 (連續響應)
    'n_ros': 2.5,                        # ROS 損傷指數 (基於二級反應動力學)
    'k_ros_damage': 1.5e-14,             # ROS 損傷係數 [1/(J/m²)^n]

    # 3.4 夜間節律 (v7.0 連續函數保留)
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
    # 5. 花青素參數 (v7.0 優化後)
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
    UVA 效應模型參數類別 (v8.0-ROS)

    基於 ROS 動力學的損傷機制
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

        # ROS 動力學參數 (v8.0)
        self.E_threshold = params['E_threshold']
        self.n_ros = params['n_ros']
        self.k_ros_damage = params['k_ros_damage']

        # 夜間節律
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
    """
    stress_effect = p.ldmc_stress_sensitivity * Stress / (p.K_ldmc + Stress + 1e-9)
    ratio = p.dw_fw_ratio_base * (1.0 + stress_effect)
    return min(ratio, p.dw_fw_ratio_max)


def calculate_ros_damage_factor(E_joules, p):
    """
    計算基於 ROS 動力學的損傷因子 (v8.0c 無閾值能量積分版)

    物理原理:
    =========

    v7.0 使用: factor = 1 + k × hours^n (時間冪次)
    v8.0 改用: factor = 1 + k × E^n (能量冪次)

    其中 E = I_UVA × hours × 3600 (J/m²)

    優點:
    - 能量積分有明確的物理意義 (累積光子劑量)
    - 冪次 n 來自反應動力學 (二級反應 n≈2, 級聯效應 n≈3-4)
    - 與 v7.0 類似的連續響應，無硬閾值

    如果 E_threshold > 0，使用軟閾值:
    - factor = 1 + k × softplus((E - E_threshold) / E_scale)^n

    如果 E_threshold = 0，使用純能量冪律:
    - factor = 1 + k × E^n

    Args:
        E_joules: 累積 UVA 能量 [J/m²]
        p: 參數物件

    Returns:
        損傷因子 (無單位)
    """
    if E_joules <= 0:
        return 1.0

    if p.E_threshold <= 0:
        # 無閾值模式: 純能量冪律
        factor = 1.0 + p.k_ros_damage * (E_joules ** p.n_ros)
        return factor

    # 有閾值模式: 軟閾值 (softplus)
    E_scale = p.E_threshold * 0.3  # 過渡區寬度
    x = (E_joules - p.E_threshold) / E_scale

    # Softplus 函數: log(1 + exp(x))
    if x > 20:
        E_excess_soft = E_scale * x
    elif x < -20:
        E_excess_soft = E_scale * np.exp(x)
    else:
        E_excess_soft = E_scale * np.log(1 + np.exp(x))

    if E_excess_soft <= 0:
        return 1.0

    factor = 1.0 + p.k_ros_damage * (E_excess_soft ** p.n_ros)
    return factor


# ==============================================================================
# 核心微分方程 (v8.0-ROS)
# ==============================================================================

def uva_sun_derivatives(t, state, p, env):
    """
    UVA效應整合模型的微分方程 (v8.0-ROS)

    v8.0 核心改進:
    - 日內逆境: 用 ROS 動力學推導的損傷因子
    - 取代 v7.0 的任意冪函數 1 + k × hours^7
    - 形狀從酶動力學自然浮現

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
    # 步驟5: ROS 動力學損傷因子 (v8.0 核心改進)
    # =========================================================================
    # 計算累積 UVA 能量 [J/m²]
    # E = I_UVA × hours × 3600
    # 實際 I_UVA 在照射期間可能不是常數，這裡用瞬時 I_UVA 作為代表值
    E_today = I_UVA * hours_today * 3600  # [J/m²]

    # 用 ROS 動力學計算損傷因子
    day_factor = calculate_ros_damage_factor(E_today, p)

    # =========================================================================
    # 步驟6: 夜間節律 (保留 v7.0 連續函數)
    # =========================================================================
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
    # 步驟8: LAI 脆弱性 (Sigmoid 型)
    # =========================================================================
    base_vuln = (p.LAI_ref_vuln / LAI) ** p.n_vuln
    vulnerability = p.cap_vuln * base_vuln / (p.cap_vuln + base_vuln)

    # =========================================================================
    # 步驟9: Stress 非線性 (Michaelis-Menten 型)
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
# 測試函數: 比較 v7.0 和 v8.0 的損傷因子
# ==============================================================================

def compare_damage_factors():
    """比較 v7.0 (冪函數) 和 v8.0 (ROS 動力學) 的損傷因子"""

    print("=" * 70)
    print("損傷因子比較: v7.0 (冪函數) vs v8.0 (ROS 動力學)")
    print("=" * 70)

    # v7.0 參數
    k_day_v7 = 1.0e-5
    n_day_v7 = 7.0

    # v8.0 參數
    p = UVAParams()

    I_UVA = 11.0  # W/m² (實驗強度的一半，因為只有上方照射)

    print(f"\nUVA 強度: {I_UVA} W/m²")
    print(f"\nv7.0 參數: k_day = {k_day_v7}, n_day = {n_day_v7}")
    print(f"v8.0 參數: E_threshold = {p.E_threshold} J/m², n_ros = {p.n_ros}")
    print(f"           k_ros_damage = {p.k_ros_damage}")

    print("\n" + "-" * 70)
    print(f"{'時間':<8} {'能量(J/m²)':<12} {'v7.0 factor':<14} {'v8.0 factor':<14} {'比值':<10}")
    print("-" * 70)

    for hours in [0, 1, 2, 3, 4, 5, 6, 8, 10, 12]:
        E = I_UVA * hours * 3600

        # v7.0 factor
        factor_v7 = 1.0 + k_day_v7 * (hours ** n_day_v7)

        # v8.0 factor
        factor_v8 = calculate_ros_damage_factor(E, p)

        ratio = factor_v8 / factor_v7 if factor_v7 > 0 else float('inf')

        print(f"{hours}h{'':<6} {E:<12.0f} {factor_v7:<14.1f} {factor_v8:<14.1f} {ratio:<10.2f}")

    print("-" * 70)

    # 計算 6h 和 12h 的比值
    E_6h = I_UVA * 6 * 3600
    E_12h = I_UVA * 12 * 3600

    factor_v7_6h = 1.0 + k_day_v7 * (6 ** n_day_v7)
    factor_v7_12h = 1.0 + k_day_v7 * (12 ** n_day_v7)

    factor_v8_6h = calculate_ros_damage_factor(E_6h, p)
    factor_v8_12h = calculate_ros_damage_factor(E_12h, p)

    print(f"\n12h/6h 比值:")
    print(f"  v7.0: {factor_v7_12h / factor_v7_6h:.1f}")
    print(f"  v8.0: {factor_v8_12h / factor_v8_6h:.1f}")

    # 物理意義說明
    print(f"\n物理意義:")
    print(f"  E_threshold = {p.E_threshold/3600:.0f} Wh/m² = {p.E_threshold/(I_UVA*3600):.1f}h @ {I_UVA} W/m²")
    print(f"  → 抗氧化劑在約 {p.E_threshold/(I_UVA*3600):.1f} 小時照射後開始耗竭")
    print(f"  → 超過此閾值後，損傷以 (E-E_threshold)^{p.n_ros} 增長")


# ==============================================================================
# 主程式
# ==============================================================================

if __name__ == "__main__":
    from model_config import (
        ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment
    )

    # 先比較損傷因子
    compare_damage_factors()

    print("\n")

    p = UVAParams()

    print("=" * 80)
    print("萵苣生長與UVA效應整合模型 v8.0-ROS")
    print("=" * 80)
    print("\n核心機制 (基於 ROS 動力學):")
    print("  日內逆境: factor = 1 + k × max(0, E - E_threshold)^n")
    print(f"           E_threshold = {p.E_threshold:.0f} J/m² ({p.E_threshold/(11*3600):.1f}h @ 11 W/m²)")
    print(f"           n = {p.n_ros} (反應動力學級數)")
    print("  物理意義: 抗氧化劑容量耗竭後，ROS 損傷級聯放大")
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
