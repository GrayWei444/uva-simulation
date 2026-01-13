# -*- coding: utf-8 -*-
"""
================================================================================
萵苣UVA模型共用設定 (Shared Configuration Module)
================================================================================
版本: v2.0
日期: 2025-12-05

功能說明:
---------
此模組集中管理所有腳本共用的設定，避免各檔案重複定義導致不同步。
包含:
  1. 環境基礎設定 (植物工廠環境)
  2. 處理組設定 (UVA照射時間和強度)
  3. 實驗目標數據 (觀測值)
  4. 資料集分割 (訓練/驗證/測試)
  5. ODE求解器設定
  6. 模擬設定
  7. 夜間抑制參數先驗區間

用法:
-----
    from model_config import ENV_BASE, TREATMENT_CONFIGS, TARGETS, ODE_SETTINGS

注意事項:
---------
  - ODE求解器設定 (RK45, max_step=60) 已校準，請勿修改!
  - 修改任何設定後請重新執行 run_validation.py 確認結果

依賴腳本:
---------
  - simulate_uva_model.py
  - run_validation.py
  - run_scenarios.py
  - optimize_parameters.py
  - sensitivity_analysis.py
================================================================================
"""

# ==============================================================================
# 第一部分: 環境基礎設定 (植物工廠環境)
# ==============================================================================
ENV_BASE = {
    # --- 光照設定 ---
    'light_on_hour': 6,       # 光照開始時間 (06:00)
    'light_off_hour': 22,     # 光照結束時間 (22:00) - 共16小時光照
    # PPFD 130 μmol/m²/s → PAR功率 = 130 × 0.219 = 28.5 W/m²
    # Sun模型假設 PAR = 0.5 × I，因此等效短波輻射 I = 28.5 / 0.5 = 57 W/m²
    'I_day': 57,              # 日間短波輻射 [W/m²] (物理計算值)

    # --- 溫度設定 ---
    'T_day': 25,              # 日間溫度 [°C]
    'T_night': 18,            # 夜間溫度 [°C]

    # --- CO2設定 ---
    'CO2_day': 1200,          # 日間CO2濃度 [ppm] (Sun模型校準值)
    'CO2_night': 1200,        # 夜間CO2濃度 [ppm]

    # --- 濕度設定 ---
    'RH_day': 0.70,           # 日間相對濕度 [0-1]
    'RH_night': 0.85,         # 夜間相對濕度 [0-1]

    # --- 種植密度 ---
    'plant_density': 36,      # 種植密度 [株/m²]
}


# ==============================================================================
# 第二部分: 處理組設定
# ==============================================================================
# 說明:
# - uva_start_day / uva_end_day: 播種後的天數
# - uva_hour_on / uva_hour_off: 每日UVA照射的開始/結束時間
# - 模擬從播種後第14天(移植日)開始，共21天，到播種後第35天結束

TREATMENT_CONFIGS = {
    # --- 對照組 ---
    'CK': {
        'uva_on': False,
        'description': '對照組 (無UVA)',
    },

    # --- 標準日間照射 (最佳處理) ---
    'L6D6': {
        'uva_on': True,
        'uva_intensity': 11.0,          # UVA強度 [W/m²] (v10.7: 實際LED功率)
        'uva_start_day': 29,            # 播種後第29天開始
        'uva_end_day': 35,              # 播種後第35天結束 (共6天)
        'uva_hour_on': 10,              # 10:00 開始
        'uva_hour_off': 16,             # 16:00 結束 (日間6小時)
        'description': '低劑量日間 (6h/day, 6天)',
    },

    # --- 夜間照射 (生理節律打斷效應) ---
    'L6D6-N': {
        'uva_on': True,
        'uva_intensity': 11.0,          # v10.7: 實際LED功率
        'uva_start_day': 29,
        'uva_end_day': 35,
        'uva_hour_on': 22,              # 22:00 開始
        'uva_hour_off': 4,              # 04:00 結束 (夜間跨夜6小時)
        'description': '低劑量夜間 (6h/night, 6天)',
    },

    # --- 高劑量脅迫 (損傷機制) ---
    'H12D3': {
        'uva_on': True,
        'uva_intensity': 11.0,          # v10.7: 實際LED功率
        'uva_start_day': 32,            # 播種後第32天開始
        'uva_end_day': 35,              # 播種後第35天結束 (共3天)
        'uva_hour_on': 6,               # 06:00 開始
        'uva_hour_off': 18,             # 18:00 結束 (12小時)
        'description': '高劑量脅迫 (12h/day, 3天)',
    },

    # --- 極低劑量長期 (適應效應) ---
    'VL3D12': {
        'uva_on': True,
        'uva_intensity': 11.0,          # v10.7: 實際LED功率
        'uva_start_day': 23,            # 播種後第23天開始
        'uva_end_day': 35,              # 播種後第35天結束 (共12天)
        'uva_hour_on': 10,              # 10:00 開始
        'uva_hour_off': 13,             # 13:00 結束 (3小時)
        'description': '極低劑量長期 (3h/day, 12天)',
    },

    # --- 低劑量長期 (形態效應主導) ---
    'L6D12': {
        'uva_on': True,
        'uva_intensity': 11.0,          # v10.7: 實際LED功率
        'uva_start_day': 23,
        'uva_end_day': 35,
        'uva_hour_on': 10,
        'uva_hour_off': 16,
        'description': '低劑量長期 (6h/day, 12天) - 形態效應累積',
    },

    # --- 正弦波漸進照射 (動態模式) ---
    # 照射強度隨時間呈正弦波變化，從0漸增到最大值再降至0
    # 週期12小時 (06:00-18:00)，總有效劑量相當於6小時連續照射
    'SIN': {
        'uva_on': True,
        'uva_mode': 'sinusoidal',        # 動態照射模式
        'uva_intensity': 22.0,           # 最大強度 [W/m²]
        'uva_start_day': 23,             # 播種後第23天開始
        'uva_end_day': 35,               # 播種後第35天結束 (共12天)
        'uva_hour_on': 6,                # 06:00 開始
        'uva_hour_off': 18,              # 18:00 結束 (12小時週期)
        'description': '正弦波漸進照射 (6h等效/day, 12天)',
    },

    # --- 間歇式照射 (動態模式) ---
    # 30分鐘開/30分鐘關，總計6小時有效照射分散在12小時內
    'INT': {
        'uva_on': True,
        'uva_mode': 'intermittent',      # 動態照射模式
        'uva_intensity': 22.0,           # 照射時強度 [W/m²]
        'uva_start_day': 23,             # 播種後第23天開始
        'uva_end_day': 35,               # 播種後第35天結束 (共12天)
        'uva_hour_on': 6,                # 06:00 開始
        'uva_hour_off': 18,              # 18:00 結束
        'intermittent_on_min': 30,       # 開啟時間 [分鐘]
        'intermittent_off_min': 30,      # 關閉時間 [分鐘]
        'description': '間歇式照射 (30min開/關, 6h等效/day, 12天)',
    },
}


# ==============================================================================
# 第三部分: 實驗目標數據 (觀測值)
# ==============================================================================
# 單位說明:
# - FW: 鮮重 [g/plant]
# - Anth: 花青素濃度 [ppm = mg/kg FW]

TARGETS = {
    # 數據更新 (2026-01-09)
    # 鮮重 (g/plant):
    #   Ref=87, 6h/6d(D)=91.4, 6h/6d(N)=80.8, 12h/3d=60.6, 3h/12d=67, 6h/12d=60.4
    # 花青素濃度 (ppm = mg/kg FW):
    #   單位已更正: 原 mg/100g → 現 mg/kg (×10)
    'CK': {'FW': 87.0, 'Anth': 433.0},       # Anth STD: 9.6
    'L6D6': {'FW': 91.4, 'Anth': 494.0},     # Anth STD: 6.8
    'L6D6-N': {'FW': 80.8, 'Anth': 493.0},   # Anth STD: 7.4
    'H12D3': {'FW': 60.6, 'Anth': 651.0},    # Anth STD: 14.6
    'VL3D12': {'FW': 67.0, 'Anth': 482.0},   # Anth STD: 2.7
    'L6D12': {'FW': 60.4, 'Anth': 518.0},    # Anth STD: 3.4
}


# ==============================================================================
# 第四部分: 資料集分割 (訓練/驗證/測試)
# ==============================================================================
# 說明:
# - 訓練集: 用於參數估計
# - 驗證集: 用於時間尺度泛化驗證
# - 測試集: SIN, INT (動態模式，尚未實現)

DATASET_SPLIT = {
    'train': ['CK', 'L6D6', 'L6D6-N', 'H12D3', 'VL3D12', 'L6D12'],  # 全部 6 組
    'validation': [],  # 無驗證集
    'test': [],
}


# ==============================================================================
# 第五部分: ODE求解器設定 (重要: 請勿修改!)
# ==============================================================================
# 說明:
# - method='RK45' 和 max_step=60 已校準，修改會導致結果不一致
# - LSODA 在不同步長下結果不穩定，故使用 RK45

ODE_SETTINGS = {
    'method': 'RK45',         # 求解器方法 (不要改!)
    'max_step': 60,           # 最大步長 [秒] (不要改!)
}


# ==============================================================================
# 第六部分: 模擬設定
# ==============================================================================
SIMULATION = {
    'days': 21,               # 模擬天數 (移植後0-21天)
    'transplant_offset': 14,  # 移植偏移 (播種後14天移植)
    'initial_fw_g': 10,       # 移植時初始鮮重 [g/plant]
}


# ==============================================================================
# 第七部分: 夜間抑制參數先驗區間 (基於文獻)
# ==============================================================================
# 參考文獻:
# - Bennie et al. (2016) J. Ecology DOI: 10.1111/1365-2745.12551
# - Deng et al. (2025) Biology DOI: 10.3390/biology14050571
# - Harmer (2009) Annu. Rev. Plant Biol.
# - Covington et al. (2008) Genome Biol.

NIGHT_INHIBITION_SCENARIOS = {
    'Low': {
        'night_uva_base_inhibition': 0.08,
        'circadian_inhibition_decay': 3.2e-5,  # 半衰期約6小時
        'description': '低抑制情境 (快速恢復)',
    },
    'Medium': {
        'night_uva_base_inhibition': 0.12,
        'circadian_inhibition_decay': 1.6e-5,  # 半衰期約12小時
        'description': '中等抑制情境 (文獻中值)',
    },
    'High': {
        'night_uva_base_inhibition': 0.18,
        'circadian_inhibition_decay': 8.0e-6,  # 半衰期約24小時
        'description': '高抑制情境 (慢速恢復)',
    },
}


# ==============================================================================
# 第八部分: 參數先驗範圍 (用於敏感度分析和優化)
# ==============================================================================
PARAM_PRIORS = {
    'night_uva_base_inhibition': {
        'range': (0.08, 0.18),
        'default': 0.12,
        'unit': '[-]',
        'source': 'Bennie et al. 2016; Deng et al. 2025',
    },
    'circadian_inhibition_decay': {
        'range': (8.0e-6, 3.2e-5),
        'default': 1.6e-5,
        'unit': '[1/s]',
        'source': 'Covington et al. 2008; Harmer 2009',
    },
    'adaptation_rate': {
        'range': (1e-7, 5e-7),
        'default': 3e-7,
        'unit': '[1/Ws]',
        'source': 'Hideg et al. 2013',
    },
    'max_adaptation': {
        'range': (0.4, 0.8),
        'default': 0.65,
        'unit': '[-]',
        'source': 'Jansen et al. 1998',
    },
}


# ==============================================================================
# 輔助函數
# ==============================================================================

def get_env_for_treatment(treatment: str) -> dict:
    """
    取得特定處理組的完整環境設定

    合併 ENV_BASE 和 TREATMENT_CONFIGS 中的設定。

    參數:
    -----
    treatment : str
        處理組代碼 (如 'CK', 'L6D6' 等)

    回傳:
    -----
    dict
        完整的環境設定字典
    """
    env = ENV_BASE.copy()
    if treatment in TREATMENT_CONFIGS:
        env.update(TREATMENT_CONFIGS[treatment])
    return env


def get_target(treatment: str) -> dict:
    """
    取得特定處理組的目標值

    參數:
    -----
    treatment : str
        處理組代碼

    回傳:
    -----
    dict
        目標值字典 {'FW': float, 'Anth': float}
    """
    return TARGETS.get(treatment, {'FW': 0, 'Anth': 0})


def get_treatments_by_dataset(dataset: str) -> list:
    """
    取得特定資料集的處理組列表

    參數:
    -----
    dataset : str
        資料集名稱 ('train', 'validation', 'test')

    回傳:
    -----
    list
        處理組代碼列表
    """
    return DATASET_SPLIT.get(dataset, [])


def print_config_summary():
    """
    印出設定摘要

    用於快速檢視當前模型設定。
    """
    print("=" * 70)
    print("萵苣UVA模型設定摘要")
    print("=" * 70)

    print("\n環境設定:")
    print(f"  光照: {ENV_BASE['light_on_hour']}:00-{ENV_BASE['light_off_hour']}:00, "
          f"{ENV_BASE['I_day']} μmol/m²/s")
    print(f"  溫度: 日間{ENV_BASE['T_day']}°C / 夜間{ENV_BASE['T_night']}°C")
    print(f"  CO2: 日間{ENV_BASE['CO2_day']} / 夜間{ENV_BASE['CO2_night']} ppm")
    print(f"  種植密度: {ENV_BASE['plant_density']} 株/m²")

    print("\n處理組:")
    for name, config in TREATMENT_CONFIGS.items():
        target = TARGETS.get(name, {})
        print(f"  {name}: {config.get('description', '')}")
        if target:
            print(f"       目標: FW={target['FW']}g, Anth={target['Anth']}ppm")

    print("\n資料集分割:")
    print(f"  訓練集: {DATASET_SPLIT['train']}")
    print(f"  驗證集: {DATASET_SPLIT['validation']}")
    print(f"  測試集: {DATASET_SPLIT['test']} (尚未實現)")

    print("\nODE設定:")
    print(f"  方法: {ODE_SETTINGS['method']}")
    print(f"  最大步長: {ODE_SETTINGS['max_step']} 秒")

    print("\n模擬設定:")
    print(f"  模擬天數: {SIMULATION['days']} 天")
    print(f"  初始鮮重: {SIMULATION['initial_fw_g']} g/plant")

    print("=" * 70)


if __name__ == '__main__':
    print_config_summary()
