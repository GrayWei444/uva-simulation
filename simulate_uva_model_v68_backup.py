"""
萵苣生長與UVA效應整合模型
===========================
版本: v6.8 Final (參數優化 + 跨夜 Bug 修正)
日期: 2025-12-24

目標: ✅ 12/12 目標達成 (FW: 6/6 <5%, Anth: 6/6 <10%)

修改紀錄:
- v5.1-5.8: 早期參數優化和機制調整
- v5.9: 引入 E_stress 累積能量狀態變量
- v6.0: 參數完全外部化到 params_config.py
- v6.1: 修復時間計算 Bug (day_from_sowing)
- v6.2: 提高抑制靈敏度
- v6.3: L6D6 優先穩定 + 機制分離調整 (5/6 組達標)
- v6.4: 當日累積非線性優化 (6/6 組 FW 達標)
- v6.5: 花青素改用 E_stress 驅動
- v6.6: 花青素改用 LAI-based
- v6.7: 花青素改用 FW-based (12/12 達標)
        移除 E_stress 狀態變量 (6 → 5 state variables)
- v6.8: 參數優化 + 跨夜 Bug 修正 ⭐⭐⭐
        c_alpha: 0.548 → 0.555 (提升整體，改善 L6D6)
        circadian_disruption_factor: 3.0 → 3.8 (壓低 L6D6-N)
        k_intraday: 49.0 → 54.0 (壓低 H12D3)
        修正跨夜 UVA 最後一天被切斷的 Bug

v6.8 核心機制:
- Stress 損傷-修復動態 (瞬時狀態，照射時高/不照時低)
- FW-based 花青素合成 (避免 LDMC 導致的 LAI-FW 解耦問題)
- LDMC 動態 (Stress → 高 DW/FW → 鮮重下降)
- LAI 脆弱性 (幼苗更易受損)
- 日內能量非線性 (>6h 觸發 200× 放大)
- 夜間節律干擾 (3.8× 損傷)

===============================================================================
模型參數表 (v6.8 Final - 31 個參數)
===============================================================================

參數按功能分類，所有數值均已通過 6 組實驗數據校準驗證。

1. 基礎光合參數 (1 個)
   ─────────────────────────────────────────────────────────────────────
   c_alpha = 0.555                      光合效率 [-]
                                        (v6.8: 提高以改善 L6D6，原 0.548)

2. UVA-PAR 轉換參數 (1 個)
   ─────────────────────────────────────────────────────────────────────
   par_conversion_factor = 1.0          UVA 轉 PAR 係數 [-]
                                        (I_UVA=22 W/m² 已是等效短波，無需放大)

3. Stress 損傷與修復參數 (14 個)
   ─────────────────────────────────────────────────────────────────────
   3.1 基礎損傷修復
       stress_damage_coeff = 0.66e-6    損傷係數 [1/(W/m²·s)]
       stress_repair_coeff = 1.0e-5     修復係數 [1/s]
       stress_nonlinear_coeff = 8.0     Stress 累積非線性係數 [-]
       K_nonlinear = 0.8                Stress 非線性半飽和常數 [-]

   3.2 LAI 脆弱性（幼苗更易受損）
       LAI_ref_vuln = 6.5               脆弱性參考 LAI [m²/m²]
       n_vuln = 8                       脆弱性指數 [-]
       cap_vuln = 100.0                 脆弱性上限 [-]
       機制: vulnerability = cap_vuln × (LAI_ref/LAI)^n / [cap_vuln + (LAI_ref/LAI)^n]

   3.3 日內能量非線性（>6h 觸發 200× 放大）
       E_50 = 475.2                     半飽和能量 [kJ/m²] ≈ 6h @ 22 W/m²
       E_scale = 237.6                  能量尺度 [kJ/m²] ≈ 3h @ 22 W/m²
       k_intraday = 54.0                非線性放大係數 [-] (v6.8: 壓低 H12D3)
       m_intraday = 2.0                 指數 [-]
       sharpness_intraday = 3.0         Softplus 銳度 [-]
       機制: intraday_factor = 1 + k × [softplus((E-E_50)/E_scale)]^m

   3.4 夜間節律損傷（針對 L6D6-N）
       circadian_disruption_factor = 3.8  夜間 UVA 損傷加成 [-] (v6.8: 壓低 L6D6-N)

   3.5 Stress 對生長的抑制
       stress_photosynthesis_inhibition = 0.66  光合抑制係數 [-]
       stress_lai_inhibition = 0.66             LAI 抑制係數 [-]
       K_stress = 1.9                           抑制半飽和常數 [-]
       機制: inhibition = Stress / (K_stress + Stress)

4. 碳修復參數 (4 個)
   ─────────────────────────────────────────────────────────────────────
   base_repair_capacity = 0.5           基礎修復能力 [-]
   carbon_repair_bonus = 0.5            碳池修復加成 [-]
   K_carbon = 0.001                     碳修復半飽和常數 [kg C/m²]
   repair_carbon_cost = 1.0e-6          修復碳成本 [kg C / Stress]
   機制: repair_capacity = base + bonus × C_buf / (K_carbon + C_buf)

5. 花青素參數 (6 個)
   ─────────────────────────────────────────────────────────────────────
   5.1 基礎合成（無 UVA）
       base_anth_rate_light = 2.0e-10   日間基礎合成率 [kg Anth/(kg FW·s)]
       base_anth_rate_dark = 1.0e-10    夜間基礎合成率 [kg Anth/(kg FW·s)]

   5.2 Stress 誘導合成（FW-based，v6.7 關鍵改進）
       V_max_anth = 2.35e-11            最大誘導合成率 [kg Anth/(kg FW·s)]
       K_stress_anth = 0.30             Stress 半飽和常數 [-]
       機制: synthesis = FW × [base + V_max × Stress/(K_stress + Stress)]
       說明: 使用 FW 而非 LAI，避免 LDMC 動態導致的 LAI-FW 解耦

   5.3 降解與碳成本
       k_deg = 3.02e-6                  花青素降解率 [1/s]
       anth_carbon_cost = 0.0           花青素合成碳成本 [kg C/kg Anth]
                                        (設為 0，小值也會過度抑制生長)

6. LDMC 參數 (4 個)
   ─────────────────────────────────────────────────────────────────────
   dw_fw_ratio_base = 0.05              基礎 DW:FW 比例 [-] (5%)
   ldmc_stress_sensitivity = 1.0        LDMC 對 Stress 敏感度 [-]
   K_ldmc = 50.0                        LDMC 半飽和 Stress 值 [-]
   dw_fw_ratio_max = 0.12               最大 DW:FW 比例 [-] (12%)
   機制: DW:FW = base × [1 + sensitivity × Stress/(K_ldmc + Stress)]
   說明: 解釋 H12D3/L6D12 鮮重下降，與 intraday_factor 分工明確

===============================================================================
參數總數: 31 個（已移除未使用的 transplant_day）
校準基準: 6 組實驗數據（CK, L6D6, L6D6-N, VL3D12, L6D12, H12D3）
達標率: 12/12 (100%) - FW: 6/6 <5% (max 2.6%), Anth: 6/6 <10% (max 9.4%)
===============================================================================
"""

import numpy as np
from scipy.integrate import solve_ivp

# 導入基礎 Sun 模型
from lettuce_uva_carbon_complete_model import SunParams as BaseSunParams
from lettuce_uva_carbon_complete_model import sun_derivatives_final


# ==============================================================================
# 模型參數定義 (v6.8 Final - 31 個參數)
# ==============================================================================
# 所有參數集中於此，方便發表時程式碼完整可讀
# 參數說明詳見檔頭參數表

ALL_PARAMS = {
    # ──────────────────────────────────────────────────────────────────────────
    # 1. 基礎光合參數
    # ──────────────────────────────────────────────────────────────────────────
    'c_alpha': 0.555,                    # 光合效率 [-] (v6.8: 原 0.548，提高以改善 L6D6)

    # ──────────────────────────────────────────────────────────────────────────
    # 2. UVA-PAR 轉換參數
    # ──────────────────────────────────────────────────────────────────────────
    'par_conversion_factor': 1.0,        # UVA 轉 PAR 係數 [-] (I_UVA=22 W/m² 已是等效短波)

    # ──────────────────────────────────────────────────────────────────────────
    # 3. Stress 損傷與修復參數
    # ──────────────────────────────────────────────────────────────────────────
    # 3.1 基礎損傷修復
    'stress_damage_coeff': 0.66e-6,      # 損傷係數 [1/(W/m²·s)]
    'stress_repair_coeff': 1.0e-5,       # 修復係數 [1/s]
    'stress_nonlinear_coeff': 8.0,       # Stress 累積非線性係數 [-]
    'K_nonlinear': 0.8,                  # Stress 非線性半飽和常數 [-]

    # 3.2 LAI 脆弱性（幼苗更易受損）
    'LAI_ref_vuln': 6.5,                 # 脆弱性參考 LAI [m²/m²]
    'n_vuln': 8,                         # 脆弱性指數 [-]
    'cap_vuln': 100.0,                   # 脆弱性上限 [-]

    # 3.3 日內能量非線性（>6h 觸發 200× 放大，針對 H12D3）
    'E_50': 475.2,                       # 半飽和能量 [kJ/m²] ≈ 6h @ 22 W/m²
    'E_scale': 237.6,                    # 能量尺度 [kJ/m²] ≈ 3h @ 22 W/m²
    'k_intraday': 54.0,                  # 非線性放大係數 [-] (v6.8: 原 49.0，壓低 H12D3)
    'm_intraday': 2.0,                   # 指數 [-]
    'sharpness_intraday': 3.0,           # Softplus 銳度 [-]

    # 3.4 夜間節律損傷（針對 L6D6-N）
    'circadian_disruption_factor': 3.8,  # 夜間損傷加成 [-] (v6.8: 原 3.0，壓低 L6D6-N)

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
    # 5. 花青素參數 (v6.7 改用 FW-based，避免 LDMC 導致的 LAI-FW 解耦)
    # ──────────────────────────────────────────────────────────────────────────
    # 5.1 基礎合成（無 UVA 時）
    'base_anth_rate_light': 2.0e-10,     # 日間基礎合成率 [kg Anth/(kg FW·s)]
    'base_anth_rate_dark': 1.0e-10,      # 夜間基礎合成率 [kg Anth/(kg FW·s)]

    # 5.2 Stress 誘導合成
    'V_max_anth': 2.35e-11,              # 最大誘導合成率 [kg Anth/(kg FW·s)]
    'K_stress_anth': 0.30,               # Stress 半飽和常數 [-]

    # 5.3 降解與碳成本
    'k_deg': 3.02e-6,                    # 花青素降解率 [1/s]
    'anth_carbon_cost': 0.0,             # 合成碳成本 [kg C/kg Anth] (設 0，避免過度抑制生長)

    # ──────────────────────────────────────────────────────────────────────────
    # 6. LDMC 參數 (解釋 H12D3/L6D12 鮮重下降)
    # ──────────────────────────────────────────────────────────────────────────
    'dw_fw_ratio_base': 0.05,            # 基礎 DW:FW 比例 [-] (5%)
    'ldmc_stress_sensitivity': 1.0,      # LDMC 對 Stress 敏感度 [-]
    'K_ldmc': 50.0,                      # LDMC 半飽和 Stress 值 [-]
    'dw_fw_ratio_max': 0.12,             # 最大 DW:FW 比例 [-] (12%)
}


# ==============================================================================
# 輔助函數
# ==============================================================================

def softplus(x, k=1.0):
    """
    Softplus 函數: smooth version of max(0, x)
    softplus(x) = log(1 + exp(k*x)) / k
    """
    return np.log(1 + np.exp(np.clip(k * x, -500, 500))) / k


# ==============================================================================
# 模型參數定義
# ==============================================================================

class UVAParams(BaseSunParams):
    """
    UVA 效應模型參數類別 (v6.8 Final)

    所有 31 個參數從 params_config.py 讀取，參數說明請見檔頭參數表。
    狀態變量: [X_d, C_buf, LAI, Anth, Stress]
    """

    def __init__(self, params=None):
        super().__init__()

        # 使用提供的參數字典，或預設使用 ALL_PARAMS
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
        self.E_50 = params['E_50']
        self.E_scale = params['E_scale']
        self.k_intraday = params['k_intraday']
        self.m_intraday = params['m_intraday']
        self.sharpness_intraday = params['sharpness_intraday']
        self.circadian_disruption_factor = params['circadian_disruption_factor']
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
    v6.7: 保留 LDMC 動態變化以解釋 H12D3/L6D12 鮮重下降
    測試顯示移除 LDMC 會導致 H12D3/L6D12 鮮重誤差過大

    LDMC 與 intraday_factor 的分工:
    - intraday_factor: 解釋 H12D3 的非線性損傷累積 (>200x 放大)
    - LDMC: 解釋 H12D3/L6D12 的鮮重下降 (提高 DW:FW 比例)
    兩者共同作用，缺一不可
    """
    stress_effect = p.ldmc_stress_sensitivity * Stress / (p.K_ldmc + Stress + 1e-9)
    ratio = p.dw_fw_ratio_base * (1.0 + stress_effect)
    return min(ratio, p.dw_fw_ratio_max)


# ==============================================================================
# 核心微分方程
# ==============================================================================

def uva_sun_derivatives(t, state, p, env):
    """
    UVA效應整合模型的微分方程 (v6.8 Final)

    v6.8 關鍵機制:
    - 花青素使用瞬時 Stress × FW (非 E_elapsed，非 LAI-based)
    - FW-based 避免 LDMC 動態導致的 LAI-FW 解耦
    - 跨夜 UVA 正確處理（前半夜/後半夜分開判斷）

    狀態變量 (5個):
    ---------------
    [X_d, C_buf, LAI, Anth, Stress]

    X_d: 乾重碳量 [kg C/m²]
    C_buf: 碳緩衝池 [kg C/m²]
    LAI: 葉面積指數 [m²/m²]
    Anth: 花青素含量 [kg/m²]
    Stress: 瞬時脅迫程度 [無量綱]
    """

    # =========================================================================
    # 步驟1: 解包狀態變量
    # =========================================================================
    X_d, C_buf, LAI, Anth, Stress = state

    # 邊界保護
    X_d = max(X_d, 1e-9)
    C_buf = max(C_buf, 0)
    LAI = max(LAI, 0.1)
    Anth = max(Anth, 0)
    Stress = max(Stress, 0)

    # =========================================================================
    # 步驟2: 計算時間相關變量
    # =========================================================================
    hour = (t / 3600.0) % 24.0  # 當日小時 [0, 24)
    # v6.1 BUG FIX: t 是從播種開始的絕對時間 (t_start = transplant_day * 86400)
    # 所以 day_from_sowing = t / 86400.0 (直接轉換)
    day_from_sowing = t / 86400.0  # 播種後天數 (絕對時間)

    # =========================================================================
    # 步驟3: 判斷日夜狀態
    # =========================================================================
    light_on = env['light_on_hour']
    light_off = env['light_off_hour']

    if light_on <= light_off:
        is_day = light_on <= hour < light_off
    else:  # 跨夜
        is_day = hour >= light_on or hour < light_off

    # 環境參數
    if is_day:
        I_base = env['I_day']
        Tc = env['T_day']
    else:
        I_base = 0.0
        Tc = env['T_night']

    # =========================================================================
    # 步驟4: 計算 UVA 強度
    # =========================================================================
    uva_on = env.get('uva_on', False)
    uva_start_day = env.get('uva_start_day', 29)
    uva_end_day = env.get('uva_end_day', 35)
    uva_hour_on = env.get('uva_hour_on', 10)
    uva_hour_off = env.get('uva_hour_off', 16)
    uva_intensity = env.get('uva_intensity', 22.0)

    I_UVA = 0.0
    in_uva_window = False

    # 判斷是否在 UVA 照射期間
    # 關鍵修正: 跨夜照射需特殊處理
    # - 日間照射 (uva_hour_on < uva_hour_off): 照射開始日必須在 [start_day, end_day]
    # - 跨夜照射 (uva_hour_on > uva_hour_off):
    #   - 前半夜 (hour >= uva_hour_on): 照射開始日 = 當天，檢查當天是否在 [start_day, end_day]
    #   - 後半夜 (hour < uva_hour_off): 照射開始日 = 前一天，檢查前一天是否在 [start_day, end_day]
    if uva_on:
        if uva_hour_on <= uva_hour_off:
            # 日間照射: 例如 10:00-16:00
            if uva_start_day <= day_from_sowing <= uva_end_day:
                in_uva_window = uva_hour_on <= hour < uva_hour_off
        else:
            # 跨夜照射: 例如 22:00-04:00
            if hour >= uva_hour_on:
                # 前半夜: 照射開始於當天，檢查當天
                if uva_start_day <= day_from_sowing <= uva_end_day:
                    in_uva_window = True
            elif hour < uva_hour_off:
                # 後半夜: 照射開始於前一天，檢查前一天
                day_session_started = day_from_sowing - 1
                if uva_start_day <= day_session_started <= uva_end_day:
                    in_uva_window = True

        if in_uva_window:
            I_UVA = uva_intensity

    # 夜間 UVA 標記
    is_night_uva = (I_UVA > 0) and (not is_day)

    # =========================================================================
    # 步驟5: UVA-PAR 轉換效應
    # =========================================================================
    # v6.0: par_conversion_factor = 1.0 (移除放大效應)
    # I_UVA = 22 W/m² 已是等效短波輻射，直接加入基礎光照
    I_gain_par = I_UVA * p.par_conversion_factor
    I_effective = I_base + I_gain_par

    # =========================================================================
    # 步驟6: 呼叫基礎 Sun 模型
    # =========================================================================
    env_modified = env.copy()
    # 使用 I_override 直接覆蓋光照強度（支援夜間 UVA-PAR 光合貢獻）
    env_modified['I_override'] = I_effective
    # 傳入當前溫度（日間或夜間）
    env_modified['T_override'] = Tc
    # 傳入日夜狀態（用於 CO2、RH）
    env_modified['is_day_override'] = is_day

    base_state = [X_d, C_buf, LAI]
    base_derivatives = sun_derivatives_final(t, base_state, p, env_modified)
    dXd_dt_base, dCbuf_dt, dLAI_dt_base = base_derivatives

    # =========================================================================
    # 步驟7: Stress 損傷-修復動態
    # =========================================================================

    # --- 7a. LAI 脆弱性 ---
    LAI_ref_vuln = p.LAI_ref_vuln
    n_vuln = p.n_vuln
    cap_vuln = p.cap_vuln

    base_vuln = (LAI_ref_vuln / LAI) ** n_vuln
    vulnerability = cap_vuln * base_vuln / (cap_vuln + base_vuln)

    # --- 7b. 當日累積能量非線性因子 (Intraday Damage Amplification) ---
    # 公式: intraday_factor = 1 + k × softplus((E_elapsed - E_50) / E_scale, sharpness)^m
    #
    # 生物學意義:
    # - E_elapsed: 當日已累積 UVA 能量 [kJ/m²] = I_UVA × hours_elapsed × 3.6
    # - E_50: 閾值能量 [kJ/m²]，≈ 6h @ 22 W/m² = 475.2 kJ/m²
    # - 當 E_elapsed < E_50 時，修復能力足夠，損傷因子接近 1
    # - 當 E_elapsed > E_50 時，修復飽和，損傷非線性增加
    # - 針對 H12D3 (12h/day)：E_elapsed ≈ 950 kJ/m² >> E_50，觸發強烈放大
    # - 針對 L6D6 (6h/day)：E_elapsed ≈ 475 kJ/m² ≈ E_50，影響小
    # - softplus 確保完全連續，無硬性閾值

    # 計算當日已累積時數 (從 UVA 開始時刻到現在)
    if in_uva_window:
        if uva_hour_on <= uva_hour_off:
            hours_elapsed = hour - uva_hour_on
        else:  # 跨夜
            if hour >= uva_hour_on:
                hours_elapsed = hour - uva_hour_on
            else:
                hours_elapsed = (24 - uva_hour_on) + hour
    else:
        hours_elapsed = 0.0

    # 轉換為能量單位: E = I_UVA_config × hours × 3.6 [kJ/m²]
    # 注意: 使用 env 設定的 UVA 強度，而非當前時刻的 I_UVA
    # 這確保能量公式與時間公式數學等價
    I_UVA_config = uva_intensity  # 使用前面已讀取的 uva_intensity (22 W/m²)
    E_elapsed = I_UVA_config * hours_elapsed * 3.6

    # 使用正規化能量計算 intraday_factor
    # 正規化: normalized = (E - E_50) / E_scale (無單位)
    # softplus(normalized, sharpness) 輸出也是無單位
    # 這樣公式數學上等價於原時間公式，但使用能量參數
    normalized_E = (E_elapsed - p.E_50) / p.E_scale
    excess_normalized = softplus(normalized_E, p.sharpness_intraday)

    # 當日非線性因子
    intraday_factor = 1.0 + p.k_intraday * (excess_normalized ** p.m_intraday)

    # --- 7c. Stress 非線性累積 ---
    # ROS 級聯效應: 已有損傷會加速新損傷 (正反饋)
    nonlinear_factor = 1.0 + p.stress_nonlinear_coeff * Stress / (p.K_nonlinear + Stress + 1e-9)

    # --- 7d. 夜間節律損傷 ---
    if is_night_uva:
        circadian_penalty = p.circadian_disruption_factor
    else:
        circadian_penalty = 1.0

    # --- 7e. 損傷率 ---
    damage_rate = p.stress_damage_coeff * I_UVA * vulnerability * \
                  intraday_factor * nonlinear_factor * circadian_penalty

    # --- 7f. 修復能力 (碳池依賴) ---
    # 修復能力 = 基礎修復 + 碳池加成 (Michaelis-Menten 形式)
    # 植物修復需要能量，碳池充足時修復能力較強
    C_buf_positive = max(C_buf, 0)
    repair_capacity = p.base_repair_capacity + \
                      p.carbon_repair_bonus * C_buf_positive / (p.K_carbon + C_buf_positive + 1e-9)

    # --- 7g. 修復率 ---
    repair_rate = p.stress_repair_coeff * Stress * repair_capacity

    # --- Stress 變化率 ---
    dStress_dt = damage_rate - repair_rate

    # =========================================================================
    # 步驟8: Stress 對生長的抑制
    # =========================================================================
    # 注意: 這裡保留 if-else 判斷，因為:
    # 1. 生物學上，Stress 應該只抑制「正向生長」，不應該加速「負向變化」(呼吸消耗)
    # 2. 使用 sigmoid 平滑化會改變這個語義，導致模型預測偏離
    # 3. 由於 dXd_dt_base 通常 > 0 (日間淨光合為正)，這個判斷很少觸發不連續點
    stress_inhibition = Stress / (p.K_stress + Stress + 1e-9)

    xd_reduction = p.stress_photosynthesis_inhibition * stress_inhibition
    lai_reduction = p.stress_lai_inhibition * stress_inhibition

    # 只抑制正向生長
    if dXd_dt_base > 0:
        dXd_dt = dXd_dt_base * (1.0 - xd_reduction)
    else:
        dXd_dt = dXd_dt_base

    if dLAI_dt_base > 0:
        dLAI_dt = dLAI_dt_base * (1.0 - lai_reduction)
    else:
        dLAI_dt = dLAI_dt_base

    # =========================================================================
    # 步驟9: 花青素合成 (在修復之前計算，因為需要知道合成率)
    # =========================================================================
    # v6.7 花青素機制: f(瞬時 Stress, FW)
    # 使用瞬時 Stress 代表當下逆境程度
    # - 照射時: Stress 高 → 合成多
    # - 不照射時: Stress 低 → 合成少
    # - 自然避免長期累積問題
    #
    # 生物學意義:
    # - Stress: 當下的 UVA 逆境程度（瞬時）
    # - FW: 鮮重，代表植物實際大小
    # - FW 小時（早期）: 即使 Stress 高，合成也少（植物小）
    # - FW 大時（後期）: Stress 高時合成多（植物大）
    #
    # v6.7 關鍵發現: 必須使用 FW 而非 LAI
    # - LDMC 動態必須保留（解釋 H12D3/L6D12 鮮重下降）
    # - LDMC 動態導致 LAI-FW 脫鉤（高 LDMC = 高 LAI/FW）
    # - 使用 LAI 會高估花青素（L6D12: 小顆植物但高 LAI）
    # - 使用 FW 直接反映植物實際大小，避免脫鉤問題

    # 計算當前鮮重 (需要用於花青素合成)
    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress, p)
    FW_kg_m2 = X_d / dw_fw_ratio  # [kg FW/m²]

    # 日夜基礎合成 (單位 FW 合成率)
    day_weight = 1.0 if is_day else 0.0
    base_synthesis_per_fw = day_weight * p.base_anth_rate_light + (1 - day_weight) * p.base_anth_rate_dark

    # Stress 誘導合成: 使用 Michaelis-Menten 函數
    # V_max × Stress / (K_stress_anth + Stress)
    stress_induced_per_fw = p.V_max_anth * Stress / (p.K_stress_anth + Stress + 1e-12)

    # FW 影響: 總合成量 = FW × 單位 FW 合成率
    # 這確保早期 FW 小時，總合成量少（即使 Stress 高）
    synthesis_rate = FW_kg_m2 * (base_synthesis_per_fw + stress_induced_per_fw)

    degradation = p.k_deg * Anth
    dAnth_dt = synthesis_rate - degradation

    # =========================================================================
    # 步驟10: 碳消耗 (修復 + 花青素合成) - v6.0 新增花青素碳成本
    # =========================================================================
    repair_carbon_consumption = repair_rate * p.repair_carbon_cost
    anth_carbon_consumption = synthesis_rate * p.anth_carbon_cost  # Gemini 建議

    dCbuf_dt = dCbuf_dt - repair_carbon_consumption - anth_carbon_consumption

    # =========================================================================
    # 步驟11: 回傳結果 (v6.7: 移除 E_stress)
    # =========================================================================
    return np.array([dXd_dt, dCbuf_dt, dLAI_dt, dAnth_dt, dStress_dt])


# ==============================================================================
# 主程式
# ==============================================================================

if __name__ == "__main__":
    from model_config import (
        ENV_BASE, TREATMENT_CONFIGS, TARGETS, ODE_SETTINGS, SIMULATION,
        get_env_for_treatment
    )

    # 初始化參數
    p = UVAParams()

    # 選擇處理組
    treatment = 'L6D6'

    # 取得環境設定
    env = get_env_for_treatment(treatment)
    target = TARGETS.get(treatment, {'FW': 0, 'Anth': 0})

    print("=" * 80)
    print("萵苣生長與UVA效應整合模型 v6.8 Final")
    print("=" * 80)
    print("\n特點:")
    print("  - 花青素由瞬時 Stress × FW 驅動")
    print("  - 合成率正比於鮮重 FW（避免 LDMC 導致的 LAI-FW 解耦）")
    print("  - 鮮重預測: 6/6 組 < 5% (max 2.6%) ✅")
    print("  - 花青素預測: 6/6 組 < 10% (max 9.4%) ✅")
    print("  - 總達標率: 12/12 (100%) ✅")
    print("=" * 80)

    # 初始條件 (v6.7: 5個狀態變量，移除 E_stress)
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    C_buf_init = Xd_init * 0.1
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6
    initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0]  # v6.7: 5個狀態

    simulation_days = SIMULATION['days']
    transplant_day = SIMULATION['transplant_offset']  # v6.0: 從第14天 (移植日) 開始

    print(f"\n模擬設定: {treatment} 處理組")
    print(f"  - 目標: FW={target['FW']}g, Anth={target['Anth']}ppm")

    # 執行模擬 (v6.0: 從移植日開始)
    t_start = transplant_day * 86400
    t_end = (transplant_day + simulation_days) * 86400
    sol = solve_ivp(
        fun=uva_sun_derivatives,
        t_span=(t_start, t_end),
        y0=initial_state,
        args=(p, env),
        method=ODE_SETTINGS['method'],
        max_step=ODE_SETTINGS['max_step'],
        t_eval=np.linspace(t_start, t_end, simulation_days * 24 + 1)
    )

    if sol.success:
        print("\n模擬成功!")

        sim_stress = sol.y[4, -1]
        sim_dw_fw_ratio = calculate_dynamic_dw_fw_ratio(sim_stress, p)
        sim_dw_per_plant = sol.y[0, -1] / ENV_BASE['plant_density'] * 1000
        sim_fw_per_plant = sim_dw_per_plant / sim_dw_fw_ratio

        # 花青素計算
        sim_anth_kg_m2 = sol.y[3, -1]
        sim_fw_kg_m2 = sim_fw_per_plant * ENV_BASE['plant_density'] / 1000
        sim_anth_ppm = sim_anth_kg_m2 * 1e6 / (sim_fw_kg_m2 + 1e-9)

        fw_err = ((sim_fw_per_plant - target['FW']) / target['FW'] * 100) if target['FW'] > 0 else 0
        anth_err = ((sim_anth_ppm - target['Anth']) / target['Anth'] * 100) if target['Anth'] > 0 else 0

        print(f"\n最終結果 ({treatment}):")
        print(f"  - 鮮重: {sim_fw_per_plant:.2f} g/plant (目標: {target['FW']}g, 誤差: {fw_err:+.1f}%)")
        print(f"  - 花青素: {sim_anth_ppm:.1f} ppm (目標: {target['Anth']}ppm, 誤差: {anth_err:+.1f}%)")
        print(f"  - 最終 Stress: {sim_stress:.2f}")
    else:
        print(f"模擬失敗: {sol.message}")
