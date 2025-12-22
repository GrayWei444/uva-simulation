#!/usr/bin/env python3
"""
統一的參數配置文件
將所有模型參數集中管理，避免散落在各個類別中
"""

# ==============================================================================
# 1. Sun 基礎模型參數
# ==============================================================================
SUN_PARAMS = {
    # 光合效率
    'c_alpha': 0.548,  # v6.0: 校準使 CK ≈ 87.0g (Sun 原始值 0.68)

    # 以下參數繼承自 Sun 原始模型，不修改
    # (在 SunParams 中定義，這裡不重複列出)
}

# ==============================================================================
# 2. UVA-PAR 轉換參數
# ==============================================================================
UVA_PAR_PARAMS = {
    'par_conversion_factor': 1.0,  # v6.0: 移除放大效應
    # 理由: I_UVA = 22 W/m² 已是等效短波輻射，直接加入即可
    # Sun 原始模型可以完美解釋 UVA 的鮮重促進效應
}

# ==============================================================================
# 3. Stress 損傷與修復參數
# ==============================================================================
STRESS_PARAMS = {
    # 損傷機制
    'stress_damage_coeff': 0.70e-6,  # 損傷係數 [1/(W/m²·s)] (v6.3: 降低以確保 L6D6 > CK)
    'stress_repair_coeff': 1.0e-5,   # 修復係數 [1/s]
    'stress_nonlinear_coeff': 8.0,   # Stress 累積非線性係數 (v6.3: 提高以改善 H12D3)
    'K_nonlinear': 0.8,              # Stress 非線性半飽和常數 (v6.3: 降低觸發點)

    # LAI 脆弱性 (幼苗更易受損)
    'LAI_ref_vuln': 6.5,   # 脆弱性參考 LAI
    'n_vuln': 8,           # 脆弱性指數 (v6.2: 保持高值，保護長期組)
    'cap_vuln': 100.0,     # 脆弱性上限

    # 日內能量非線性 (修復飽和)
    'E_50': 800.0,           # 半飽和能量 [kJ/m²] (v6.2: 保持高值)
    'E_scale': 150.0,        # 能量尺度 [kJ/m²]
    'k_intraday': 2.0,       # 非線性放大係數
    'm_intraday': 2.0,       # 指數
    'sharpness_intraday': 3.0,  # softplus 銳度

    # 夜間節律損傷
    'circadian_disruption_factor': 3.2,  # 夜間 UVA 損傷加成 (v6.3: 提高以調整 L6D6-N)

    # Stress 對生長的抑制
    'stress_photosynthesis_inhibition': 0.68,  # 光合抑制係數 (v6.2: 溫和提高)
    'stress_lai_inhibition': 0.68,             # LAI 抑制係數 (v6.2: 同步)
    'K_stress': 1.8,                           # 抑制半飽和常數 (v6.2: 適中)
}

# ==============================================================================
# 4. 碳修復參數
# ==============================================================================
CARBON_REPAIR_PARAMS = {
    'base_repair_capacity': 0.5,     # 基礎修復能力
    'carbon_repair_bonus': 0.5,      # 碳池加成
    'K_carbon': 0.001,               # 半飽和常數 [kg C/m²]
    'repair_carbon_cost': 1.0e-6,    # 每單位修復消耗的碳 [kg C / Stress]
}

# ==============================================================================
# 5. 花青素參數
# ==============================================================================
ANTHOCYANIN_PARAMS = {
    # 基礎合成 (無 UVA 時)
    'base_anth_rate_light': 4.0e-10,   # 日間基礎合成率 [kg/m²/s]
    'base_anth_rate_dark': 2.0e-10,    # 夜間基礎合成率

    # Stress 誘導合成 (Michaelis-Menten 型)
    'V_max_anth': 1.8e-10,      # 最大誘導合成率 [kg/m²/s]
    'K_stress_anth': 10.0,      # Stress 半飽和常數
    'n_stress_anth': 1.0,       # Hill 係數 (1 = Michaelis-Menten)

    # 降解
    'k_deg': 2.6e-6,             # 花青素降解率 [1/s] (半衰期 ~3天)

    # 碳成本 (v6.0 新增，但發現會嚴重抑制生長)
    'anth_carbon_cost': 0.0,    # 每單位花青素合成消耗的碳 [kg C / kg Anth]
    # 注意: 設為 0 是因為發現即使是小值也會過度抑制生長
    # 這可能表示基礎合成率設定有問題
}

# ==============================================================================
# 6. LDMC (葉乾物質含量) 參數
# ==============================================================================
LDMC_PARAMS = {
    'dw_fw_ratio_base': 0.05,           # 基礎 DW:FW 比例 (5%)
    'ldmc_stress_sensitivity': 1.0,     # LDMC 對 Stress 敏感度
    'K_ldmc': 50.0,                     # LDMC 半飽和 Stress 值
    'dw_fw_ratio_max': 0.12,            # 最大 DW:FW 比例 (12%)
}

# ==============================================================================
# 7. 其他參數
# ==============================================================================
OTHER_PARAMS = {
    'transplant_day': 14,  # 移植日 = 播種後第14天
}

# ==============================================================================
# 統一的參數字典 (整合所有參數)
# ==============================================================================
ALL_PARAMS = {}
ALL_PARAMS.update(SUN_PARAMS)
ALL_PARAMS.update(UVA_PAR_PARAMS)
ALL_PARAMS.update(STRESS_PARAMS)
ALL_PARAMS.update(CARBON_REPAIR_PARAMS)
ALL_PARAMS.update(ANTHOCYANIN_PARAMS)
ALL_PARAMS.update(LDMC_PARAMS)
ALL_PARAMS.update(OTHER_PARAMS)

# ==============================================================================
# 參數檢查與驗證
# ==============================================================================
def validate_params(params):
    """
    驗證參數的合理性
    """
    issues = []

    # 檢查花青素基礎合成率
    if params['base_anth_rate_light'] > 1e-9:
        issues.append(f"警告: base_anth_rate_light = {params['base_anth_rate_light']:.2e} 可能過高")

    # 檢查花青素碳成本
    if params['anth_carbon_cost'] > 0.1:
        issues.append(f"警告: anth_carbon_cost = {params['anth_carbon_cost']} 可能會嚴重抑制生長")

    # 檢查 c_alpha
    if params['c_alpha'] < 0.5 or params['c_alpha'] > 0.7:
        issues.append(f"警告: c_alpha = {params['c_alpha']} 超出合理範圍 [0.5, 0.7]")

    return issues

# ==============================================================================
# 打印參數摘要
# ==============================================================================
def print_params_summary():
    """打印參數摘要"""
    print('=' * 80)
    print('模型參數配置摘要')
    print('=' * 80)
    print()

    sections = [
        ('Sun 基礎模型', SUN_PARAMS),
        ('UVA-PAR 轉換', UVA_PAR_PARAMS),
        ('Stress 損傷與修復', STRESS_PARAMS),
        ('碳修復', CARBON_REPAIR_PARAMS),
        ('花青素', ANTHOCYANIN_PARAMS),
        ('LDMC', LDMC_PARAMS),
        ('其他', OTHER_PARAMS),
    ]

    for section_name, section_params in sections:
        print(f'{section_name}:')
        print('-' * 80)
        for key, value in section_params.items():
            if isinstance(value, float):
                if abs(value) < 1e-3 or abs(value) > 1e3:
                    print(f'  {key:40s}: {value:.3e}')
                else:
                    print(f'  {key:40s}: {value:.4f}')
            else:
                print(f'  {key:40s}: {value}')
        print()

    # 參數驗證
    issues = validate_params(ALL_PARAMS)
    if issues:
        print('參數問題:')
        print('-' * 80)
        for issue in issues:
            print(f'  ⚠️  {issue}')
        print()

if __name__ == "__main__":
    print_params_summary()
