# 參數審計摘要 (Parameter Audit Summary)

**日期**: 2025-12-23
**模型版本**: v6.7 (FW-based Anthocyanin)

---

## 核心結論

### ✅ 所有參數均已使用

- **總參數數**: 33個
- **已使用**: 33個 (100%)
- **未使用**: 0個
- **冗餘**: 0個

### 無需清理任何參數

---

## 快速分類

| 類別 | 參數數 | 全部使用 |
|------|--------|----------|
| Sun 基礎模型 | 1 | ✅ |
| UVA-PAR 轉換 | 1 | ✅ |
| Stress 損傷與修復 | 16 | ✅ |
| 碳修復 | 4 | ✅ |
| 花青素 | 6 | ✅ |
| LDMC | 4 | ✅ |
| 其他 | 1 | ✅ |
| **總計** | **33** | **✅** |

---

## 特殊參數說明

### 1. anth_carbon_cost = 0.0

- **狀態**: 設為 0（花青素合成無碳成本）
- **原因**: 非零值會過度抑制生長（已測試）
- **建議**: ✅ 保留（記錄設計決策）
- **位置**: params_config.py:93

### 2. par_conversion_factor = 1.0

- **狀態**: 設為 1.0（UVA-PAR 無放大效應）
- **歷史**: v6.0 移除放大效應（v5.x 為 3.0）
- **原因**: Sun 原始模型足以解釋 UVA 鮮重促進效應
- **建議**: ✅ 保留（未來調整彈性）
- **位置**: params_config.py:22

---

## v6.7 升級確認

### ✅ 花青素機制已正確升級

- **舊機制 (v6.6)**: LAI-based (synthesis_rate = LAI × ...)
- **新機制 (v6.7)**: FW-based (synthesis_rate = FW × ...)
- **狀態變量**: 5個 (X_d, C_buf, LAI, Anth, Stress)

### ✅ E_stress 已正確移除

- **v5.9 - v6.5**: 使用 E_stress（Stress 累積能量）驅動花青素
- **v6.6**: 改用 LAI × Stress
- **v6.7**: 改用 FW × Stress
- **確認**: 無遺留 E_stress 相關參數或代碼

### ✅ 參數正確性驗證

所有花青素參數均正確使用在 FW-based 機制中：

| 參數 | 使用位置 | 用途 |
|------|----------|------|
| base_anth_rate_light | simulate_uva_model.py:421 | 日間基礎合成率 [kg/(kg FW·s)] |
| base_anth_rate_dark | simulate_uva_model.py:421 | 夜間基礎合成率 [kg/(kg FW·s)] |
| V_max_anth | simulate_uva_model.py:425 | Stress 誘導最大合成率 [kg/(kg FW·s)] |
| K_stress_anth | simulate_uva_model.py:425 | Stress 半飽和常數 |
| k_deg | simulate_uva_model.py:431 | 降解速率 [1/s] |
| anth_carbon_cost | simulate_uva_model.py:438 | 合成碳成本（= 0.0）|

---

## 完整參數列表

### 1. Sun 基礎模型 (1個)
- c_alpha = 0.548

### 2. UVA-PAR 轉換 (1個)
- par_conversion_factor = 1.0

### 3. Stress 損傷與修復 (16個)

**損傷機制 (4個)**:
- stress_damage_coeff = 6.6e-7
- stress_repair_coeff = 1.0e-5
- stress_nonlinear_coeff = 8.0
- K_nonlinear = 0.8

**LAI 脆弱性 (3個)**:
- LAI_ref_vuln = 6.5
- n_vuln = 8
- cap_vuln = 100.0

**日內能量非線性 (5個)**:
- E_50 = 475.2
- E_scale = 237.6
- k_intraday = 49.0
- m_intraday = 2.0
- sharpness_intraday = 3.0

**節律與抑制 (4個)**:
- circadian_disruption_factor = 3.0
- stress_photosynthesis_inhibition = 0.66
- stress_lai_inhibition = 0.66
- K_stress = 1.9

### 4. 碳修復 (4個)
- base_repair_capacity = 0.5
- carbon_repair_bonus = 0.5
- K_carbon = 0.001
- repair_carbon_cost = 1.0e-6

### 5. 花青素 (6個)
- base_anth_rate_light = 2.0e-10
- base_anth_rate_dark = 1.0e-10
- V_max_anth = 2.35e-11
- K_stress_anth = 0.30
- k_deg = 3.02e-6
- anth_carbon_cost = 0.0

### 6. LDMC (4個)
- dw_fw_ratio_base = 0.05
- ldmc_stress_sensitivity = 1.0
- K_ldmc = 50.0
- dw_fw_ratio_max = 0.12

### 7. 其他 (1個)
- transplant_day = 14

---

## 檢查清單

- [x] 所有 params_config.py 中的參數均已在 UVAParams.__init__ 中載入
- [x] 所有載入的參數均在 uva_sun_derivatives 函數中使用
- [x] 無遺留的 E_stress 相關參數
- [x] 花青素機制已正確改為 FW-based (v6.7)
- [x] 無未使用或冗餘的參數
- [x] 特殊參數（anth_carbon_cost, par_conversion_factor）均有充分理由

---

## 建議

### ✅ 當前狀態良好，無需修改

所有參數均有明確用途，代碼結構清晰，參數管理規範。

### 文件參考

- 完整審計報告: `/home/kasm-user/projects/uva-simulation/PARAMETER_AUDIT.md`
- 參數配置: `/home/kasm-user/projects/uva-simulation/params_config.py`
- 模型主程式: `/home/kasm-user/projects/uva-simulation/simulate_uva_model.py`

---

**審計完成**: 2025-12-23
**審計工具**: Claude Code
**結論**: ✅ 無需清理，所有參數均正確使用
