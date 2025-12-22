"""
萵苣生長與UVA效應整合模型
===========================
版本: v5.9 (Stress 累積能量驅動花青素合成)
日期: 2025-12-22

目標：改善花青素預測精度 (VL3D12, L6D12)

修改紀錄:
- v5.1: 移除硬編碼參數到 UVAParams
- v5.2: 花青素校準 (從低值累積)
- v5.3: 程式碼審查 (softplus 閾值、保留生長抑制 if-else)
- v5.4: 移除指數型天數累積 (exp_nonlinear_beta, exp_nonlinear_D_ref)
- v5.5: 修正夜間 UVA-PAR 光合貢獻，使用 I_override 繞過 Sun 模型日夜判斷
- v5.6: 修正 repair_capacity，現在使用 base_repair_capacity 和 carbon_repair_bonus 參數
- v5.7: 將時間閾值改為能量單位公式 [kJ/m²]
- v5.8: 移除 par_conversion_factor 放大效應 (1.0)，恢復 c_alpha = 0.68
        理由: I_UVA = 22 W/m² 已是等效短波輻射，不需再乘以 3.0
- v5.9: 移除硬閾值 stress_threshold_anth，改用 Stress 累積能量 (E_stress)
        使用 Hill 函數描述花青素誘導，無硬編碼閾值，完全連續
        新增狀態變量 E_stress [Stress·day]
- v6.0: 移除 par_conversion_factor 放大效應，改為 1.0
        理由: Sun 原始模型已能正確模擬 CK 和 L6D6，不需額外放大
        UVA 22 W/m² 已是等效短波輻射，直接加入即可
"""

import numpy as np
from scipy.integrate import solve_ivp

# 導入基礎 Sun 模型
from lettuce_uva_carbon_complete_model import SunParams as BaseSunParams
from params_config import ALL_PARAMS
from lettuce_uva_carbon_complete_model import sun_derivatives_final


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
    UVA效應模型參數類別 (v6.0)

    v6.0 更新: 所有參數從 params_config.py 讀取，不再硬編碼
    """

    def __init__(self, params=None):
        super().__init__()

        # 使用提供的參數字典，或預設使用 ALL_PARAMS
        if params is None:
            params = ALL_PARAMS

        # ==================================================================
        # 1. 基礎模型校準參數
        # ==================================================================
        self.c_alpha = params['c_alpha']

        # ==================================================================
        # 2. UVA-PAR 轉換參數
        # ==================================================================
        self.par_conversion_factor = params['par_conversion_factor']

        # ==================================================================
        # 3. Stress 損傷與修復參數
        # ==================================================================
        self.stress_damage_coeff = params['stress_damage_coeff']
        self.stress_repair_coeff = params['stress_repair_coeff']
        self.stress_nonlinear_coeff = params['stress_nonlinear_coeff']
        self.K_nonlinear = params['K_nonlinear']

        # LAI 脆弱性
        self.LAI_ref_vuln = params['LAI_ref_vuln']
        self.n_vuln = params['n_vuln']
        self.cap_vuln = params['cap_vuln']

        # 日內能量非線性
        self.E_50 = params['E_50']
        self.E_scale = params['E_scale']
        self.k_intraday = params['k_intraday']
        self.m_intraday = params['m_intraday']
        self.sharpness_intraday = params['sharpness_intraday']

        # 夜間節律損傷
        self.circadian_disruption_factor = params['circadian_disruption_factor']

        # Stress 對生長的抑制
        self.stress_photosynthesis_inhibition = params['stress_photosynthesis_inhibition']
        self.stress_lai_inhibition = params['stress_lai_inhibition']
        self.K_stress = params['K_stress']

        # ==================================================================
        # 4. 碳修復參數
        # ==================================================================
        self.base_repair_capacity = params['base_repair_capacity']
        self.carbon_repair_bonus = params['carbon_repair_bonus']
        self.K_carbon = params['K_carbon']
        self.repair_carbon_cost = params['repair_carbon_cost']

        # ==================================================================
        # 5. 花青素參數
        # ==================================================================
        self.base_anth_rate_light = params['base_anth_rate_light']
        self.base_anth_rate_dark = params['base_anth_rate_dark']
        self.V_max_anth = params['V_max_anth']
        self.K_stress_anth = params['K_stress_anth']
        self.n_stress_anth = params['n_stress_anth']
        self.k_deg = params['k_deg']
        self.anth_carbon_cost = params['anth_carbon_cost']

        # ==================================================================
        # 6. LDMC 參數
        # ==================================================================
        self.dw_fw_ratio_base = params['dw_fw_ratio_base']
        self.ldmc_stress_sensitivity = params['ldmc_stress_sensitivity']
        self.K_ldmc = params['K_ldmc']
        self.dw_fw_ratio_max = params['dw_fw_ratio_max']

        # ==================================================================
        # 7. 其他參數
        # ==================================================================
        self.transplant_day = params['transplant_day']


# ==============================================================================
# 輔助函數
# ==============================================================================

def calculate_dynamic_dw_fw_ratio(Stress, p):
    """
    根據 Stress 計算動態 DW:FW 比例
    """
    stress_effect = p.ldmc_stress_sensitivity * Stress / (p.K_ldmc + Stress + 1e-9)
    ratio = p.dw_fw_ratio_base * (1.0 + stress_effect)
    return min(ratio, p.dw_fw_ratio_max)


def get_daily_hours(uva_hour_on, uva_hour_off):
    """
    計算每日UVA照射時數 (支援跨夜)
    """
    if uva_hour_on < uva_hour_off:
        return uva_hour_off - uva_hour_on
    else:
        return (24 - uva_hour_on) + uva_hour_off


# ==============================================================================
# 核心微分方程
# ==============================================================================

def uva_sun_derivatives(t, state, p, env):
    """
    UVA效應整合模型的微分方程 (v5.9 - Stress 累積能量驅動花青素)

    狀態變量 (6個):
    ---------------
    [X_d, C_buf, LAI, Anth, Stress, E_stress]

    X_d: 乾重碳量 [kg C/m²]
    C_buf: 碳緩衝池 [kg C/m²]
    LAI: 葉面積指數 [m²/m²]
    Anth: 花青素含量 [kg/m²]
    Stress: 瞬時脅迫程度 [無量綱]
    E_stress: Stress 累積能量 [Stress·day]，用於花青素誘導
    """

    # =========================================================================
    # 步驟1: 解包狀態變量
    # =========================================================================
    X_d, C_buf, LAI, Anth, Stress, E_stress = state

    # 邊界保護
    X_d = max(X_d, 1e-9)
    C_buf = max(C_buf, 0)
    LAI = max(LAI, 0.1)
    Anth = max(Anth, 0)
    Stress = max(Stress, 0)
    E_stress = max(E_stress, 0)

    # =========================================================================
    # 步驟2: 計算時間相關變量
    # =========================================================================
    hour = (t / 3600.0) % 24.0  # 當日小時 [0, 24)
    # v6.1 BUG FIX: t 是從播種開始的絕對時間 (t_start = transplant_day * 86400)
    # 所以 day_from_sowing = t / 86400.0 (直接轉換)
    day_from_sowing = t / 86400.0  # 播種後天數 (絕對時間)
    day_from_transplant = day_from_sowing - p.transplant_day  # 移植後天數 (相對時間)

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
        Xc_ppm = env['CO2_day']
        Xh = env['RH_day']
    else:
        I_base = 0.0
        Tc = env['T_night']
        Xc_ppm = env['CO2_night']
        Xh = env['RH_night']

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

    if uva_on and uva_start_day <= day_from_sowing < uva_end_day:
        if uva_hour_on <= uva_hour_off:
            in_uva_window = uva_hour_on <= hour < uva_hour_off
        else:  # 跨夜照射
            in_uva_window = hour >= uva_hour_on or hour < uva_hour_off

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
    # 公式: intraday_factor = 1 + k × softplus(E_elapsed - E_repair, sharpness)^m
    #
    # 生物學意義:
    # - E_elapsed: 當日已累積 UVA 能量 [kJ/m²] = I_UVA × hours × 3.6
    # - E_repair: 修復系統飽和能量 [kJ/m²]
    # - 當 E < E_repair 時，修復能力足夠，損傷因子接近 1
    # - 當 E > E_repair 時，修復飽和，損傷非線性增加
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
    I_UVA_config = env.get('I_UVA', 11.0)  # 從 env 取得設定的 UVA 強度
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
    # 日夜基礎合成
    day_weight = 1.0 if is_day else 0.0
    base_synthesis = day_weight * p.base_anth_rate_light + (1 - day_weight) * p.base_anth_rate_dark

    # UVA 誘導合成: v6.0 單一 Stress 驅動
    Stress_power_n = Stress ** p.n_stress_anth
    K_power_n = p.K_stress_anth ** p.n_stress_anth
    uva_induced = p.V_max_anth * Stress_power_n / (K_power_n + Stress_power_n + 1e-12)

    synthesis_rate = base_synthesis + uva_induced
    degradation = p.k_deg * Anth
    dAnth_dt = synthesis_rate - degradation

    # =========================================================================
    # 步驟10: 碳消耗 (修復 + 花青素合成) - v6.0 新增花青素碳成本
    # =========================================================================
    repair_carbon_consumption = repair_rate * p.repair_carbon_cost
    anth_carbon_consumption = synthesis_rate * p.anth_carbon_cost  # Gemini 建議

    dCbuf_dt = dCbuf_dt - repair_carbon_consumption - anth_carbon_consumption

    # =========================================================================
    # 步驟11: Stress 累積能量
    # =========================================================================
    # E_stress 是 Stress 對時間的積分 [Stress·day]
    dE_stress_dt = Stress / 86400.0

    # =========================================================================
    # 步驟12: 回傳結果
    # =========================================================================
    return np.array([dXd_dt, dCbuf_dt, dLAI_dt, dAnth_dt, dStress_dt, dE_stress_dt])


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
    print("萵苣生長與UVA效應整合模型 v6.0 (Stress驅動花青素，優化參數)")
    print("=" * 80)
    print("\n特點:")
    print("  - 花青素由瞬時 Stress 驅動 (非累積 E_stress)")
    print("  - 使用 Hill 函數，靈活參數化")
    print("  - 基於數據分析優化參數")
    print("  - 目標: 所有處理組花青素誤差 < 5%")
    print("=" * 80)

    # 初始條件 (v6.0: 修復 C_buf 初始化)
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    C_buf_init = Xd_init * 0.1  # v6.0 BUG FIX: C_buf 應該是 X_d * 0.1，不是 0
    LAI_init = (dw_init_g / 0.01) * 0.04
    # 花青素初始值設為很低 (約 5 ppm)，後續會累積
    # Anth_init [kg/m²] = 5 ppm × FW_total [kg/m²] / 1e6
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000  # kg/m²
    Anth_init = 5.0 * fw_total_init / 1e6  # 從約 5 ppm 開始
    initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]

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
        sim_E_stress = sol.y[5, -1]  # v5.9: 新增 E_stress 輸出
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
        print(f"  - 累積 E_stress: {sim_E_stress:.2f} Stress·day")
    else:
        print(f"模擬失敗: {sol.message}")
