# 參數使用審計報告 (Parameter Usage Audit)

**日期**: 2025-12-23
**版本**: v6.7 (FW-based Anthocyanin)

---

## 執行摘要

本報告審計了 `params_config.py` 和 `simulate_uva_model.py` 中定義的所有參數，檢查是否有未使用或冗餘的參數。

### 總體結果

- **總參數數**: 33個
- **已使用參數**: 33個 (100%)
- **未使用參數**: 0個 ✅
- **設為0的參數**: 1個 (anth_carbon_cost)

### 分類統計

| 類別 | 參數數 |
|------|--------|
| Sun 基礎模型 | 1 |
| UVA-PAR 轉換 | 1 |
| Stress 損傷與修復 | 16 |
| 碳修復 | 4 |
| 花青素 | 6 |
| LDMC | 4 |
| 其他 | 1 |
| **總計** | **33** |

---

## 詳細分析

### 1. Sun 基礎模型參數

| 參數 | 值 | 定義位置 | 使用位置 | 狀態 |
|------|-----|----------|----------|------|
| c_alpha | 0.548 | params_config.py:12 | UVAParams.__init__:84 | ✅ 已使用 |

**說明**: `c_alpha` 在 `UVAParams` 初始化時賦值給 `self.c_alpha`，並由父類 `BaseSunParams` 在 `sun_derivatives_final` 中使用（光合效率參數）。

---

### 2. UVA-PAR 轉換參數

| 參數 | 值 | 定義位置 | 使用位置 | 狀態 |
|------|-----|----------|----------|------|
| par_conversion_factor | 1.0 | params_config.py:22 | simulate_uva_model.py:274 | ✅ 已使用 |

**使用代碼**:
```python
I_gain_par = I_UVA * p.par_conversion_factor
```

**說明**: v6.0 移除放大效應，設為 1.0。雖然值為 1.0 使其效果等同於不乘，但保留此參數可提供未來調整的彈性。

**建議**: 保留（設計彈性）

---

### 3. Stress 損傷與修復參數

#### 3.1 損傷機制

| 參數 | 值 | 定義位置 | 使用位置 | 狀態 |
|------|-----|----------|----------|------|
| stress_damage_coeff | 0.66e-6 | params_config.py:32 | simulate_uva_model.py:355 | ✅ 已使用 |
| stress_repair_coeff | 1.0e-5 | params_config.py:33 | simulate_uva_model.py:366 | ✅ 已使用 |
| stress_nonlinear_coeff | 8.0 | params_config.py:34 | simulate_uva_model.py:346 | ✅ 已使用 |
| K_nonlinear | 0.8 | params_config.py:35 | simulate_uva_model.py:346 | ✅ 已使用 |

**使用代碼**:
```python
# 損傷率計算
damage_rate = p.stress_damage_coeff * I_UVA * vulnerability * \
              intraday_factor * nonlinear_factor * circadian_penalty

# Stress 非線性累積
nonlinear_factor = 1.0 + p.stress_nonlinear_coeff * Stress / (p.K_nonlinear + Stress + 1e-9)

# 修復率計算
repair_rate = p.stress_repair_coeff * Stress * repair_capacity
```

#### 3.2 LAI 脆弱性

| 參數 | 值 | 定義位置 | 使用位置 | 狀態 |
|------|-----|----------|----------|------|
| LAI_ref_vuln | 6.5 | params_config.py:38 | simulate_uva_model.py:297 | ✅ 已使用 |
| n_vuln | 8 | params_config.py:39 | simulate_uva_model.py:298 | ✅ 已使用 |
| cap_vuln | 100.0 | params_config.py:40 | simulate_uva_model.py:299 | ✅ 已使用 |

**使用代碼**:
```python
LAI_ref_vuln = p.LAI_ref_vuln
n_vuln = p.n_vuln
cap_vuln = p.cap_vuln

base_vuln = (LAI_ref_vuln / LAI) ** n_vuln
vulnerability = cap_vuln * base_vuln / (cap_vuln + base_vuln)
```

#### 3.3 日內能量非線性

| 參數 | 值 | 定義位置 | 使用位置 | 狀態 |
|------|-----|----------|----------|------|
| E_50 | 475.2 | params_config.py:44 | simulate_uva_model.py:338 | ✅ 已使用 |
| E_scale | 237.6 | params_config.py:45 | simulate_uva_model.py:338 | ✅ 已使用 |
| k_intraday | 49.0 | params_config.py:46 | simulate_uva_model.py:342 | ✅ 已使用 |
| m_intraday | 2.0 | params_config.py:47 | simulate_uva_model.py:342 | ✅ 已使用 |
| sharpness_intraday | 3.0 | params_config.py:48 | simulate_uva_model.py:339 | ✅ 已使用 |

**使用代碼**:
```python
normalized_E = (E_elapsed - p.E_50) / p.E_scale
excess_normalized = softplus(normalized_E, p.sharpness_intraday)
intraday_factor = 1.0 + p.k_intraday * (excess_normalized ** p.m_intraday)
```

#### 3.4 夜間節律損傷

| 參數 | 值 | 定義位置 | 使用位置 | 狀態 |
|------|-----|----------|----------|------|
| circadian_disruption_factor | 3.0 | params_config.py:51 | simulate_uva_model.py:350 | ✅ 已使用 |

**使用代碼**:
```python
if is_night_uva:
    circadian_penalty = p.circadian_disruption_factor
else:
    circadian_penalty = 1.0
```

#### 3.5 Stress 對生長的抑制

| 參數 | 值 | 定義位置 | 使用位置 | 狀態 |
|------|-----|----------|----------|------|
| stress_photosynthesis_inhibition | 0.66 | params_config.py:54 | simulate_uva_model.py:380 | ✅ 已使用 |
| stress_lai_inhibition | 0.66 | params_config.py:55 | simulate_uva_model.py:381 | ✅ 已使用 |
| K_stress | 1.9 | params_config.py:56 | simulate_uva_model.py:378 | ✅ 已使用 |

**使用代碼**:
```python
stress_inhibition = Stress / (p.K_stress + Stress + 1e-9)
xd_reduction = p.stress_photosynthesis_inhibition * stress_inhibition
lai_reduction = p.stress_lai_inhibition * stress_inhibition
```

---

### 4. 碳修復參數

| 參數 | 值 | 定義位置 | 使用位置 | 狀態 |
|------|-----|----------|----------|------|
| base_repair_capacity | 0.5 | params_config.py:63 | simulate_uva_model.py:362 | ✅ 已使用 |
| carbon_repair_bonus | 0.5 | params_config.py:64 | simulate_uva_model.py:363 | ✅ 已使用 |
| K_carbon | 0.001 | params_config.py:65 | simulate_uva_model.py:363 | ✅ 已使用 |
| repair_carbon_cost | 1.0e-6 | params_config.py:66 | simulate_uva_model.py:437 | ✅ 已使用 |

**使用代碼**:
```python
# 修復能力計算
repair_capacity = p.base_repair_capacity + \
                  p.carbon_repair_bonus * C_buf_positive / (p.K_carbon + C_buf_positive + 1e-9)

# 碳消耗計算
repair_carbon_consumption = repair_rate * p.repair_carbon_cost
```

---

### 5. 花青素參數

| 參數 | 值 | 定義位置 | 使用位置 | 狀態 |
|------|-----|----------|----------|------|
| base_anth_rate_light | 2.0e-10 | params_config.py:74 | simulate_uva_model.py:421 | ✅ 已使用 |
| base_anth_rate_dark | 1.0e-10 | params_config.py:75 | simulate_uva_model.py:421 | ✅ 已使用 |
| V_max_anth | 2.35e-11 | params_config.py:86 | simulate_uva_model.py:425 | ✅ 已使用 |
| K_stress_anth | 0.30 | params_config.py:87 | simulate_uva_model.py:425 | ✅ 已使用 |
| k_deg | 3.02e-6 | params_config.py:90 | simulate_uva_model.py:431 | ✅ 已使用 |
| anth_carbon_cost | 0.0 | params_config.py:93 | simulate_uva_model.py:438 | ⚠️ 設為 0 |

**使用代碼**:
```python
# 基礎合成
base_synthesis_per_fw = day_weight * p.base_anth_rate_light + (1 - day_weight) * p.base_anth_rate_dark

# Stress 誘導合成
stress_induced_per_fw = p.V_max_anth * Stress / (p.K_stress_anth + Stress + 1e-12)

# 降解
degradation = p.k_deg * Anth

# 碳成本
anth_carbon_consumption = synthesis_rate * p.anth_carbon_cost
```

**說明**: `anth_carbon_cost` 設為 0.0 是因為發現即使是小值也會過度抑制生長（見 params_config.py:93-95 註解）。雖然值為 0 使其效果等同於無成本，但保留此參數可提供未來調整的彈性。

**建議**: 保留（設計彈性 + 文件記錄）

---

### 6. LDMC 參數

| 參數 | 值 | 定義位置 | 使用位置 | 狀態 |
|------|-----|----------|----------|------|
| dw_fw_ratio_base | 0.05 | params_config.py:102 | simulate_uva_model.py:167, 480 | ✅ 已使用 |
| ldmc_stress_sensitivity | 1.0 | params_config.py:103 | simulate_uva_model.py:166 | ✅ 已使用 |
| K_ldmc | 50.0 | params_config.py:104 | simulate_uva_model.py:166 | ✅ 已使用 |
| dw_fw_ratio_max | 0.12 | params_config.py:105 | simulate_uva_model.py:168 | ✅ 已使用 |

**使用代碼**:
```python
def calculate_dynamic_dw_fw_ratio(Stress, p):
    stress_effect = p.ldmc_stress_sensitivity * Stress / (p.K_ldmc + Stress + 1e-9)
    ratio = p.dw_fw_ratio_base * (1.0 + stress_effect)
    return min(ratio, p.dw_fw_ratio_max)
```

---

### 7. 其他參數

| 參數 | 值 | 定義位置 | 使用位置 | 狀態 |
|------|-----|----------|----------|------|
| transplant_day | 14 | params_config.py:112 | simulate_uva_model.py:219 | ✅ 已使用 |

**使用代碼**:
```python
day_from_transplant = day_from_sowing - p.transplant_day
```

---

## v6.7 特定檢查: 花青素機制變更

### v6.7 移除的狀態變量

**E_stress** (Stress 累積能量):
- **狀態**: ✅ 已移除
- **確認**:
  - `uva_sun_derivatives` 只有 5 個狀態變量 (simulate_uva_model.py:189-197)
  - 初始條件只有 5 個元素 (simulate_uva_model.py:486)
  - 沒有 `dE_stress_dt` 計算

### v6.7 花青素機制驗證

**當前機制**: FW-based (v6.7)
- **驅動變量**: 瞬時 Stress + FW
- **合成公式**: `synthesis_rate = FW_kg_m2 * (base_synthesis_per_fw + stress_induced_per_fw)`
- **使用參數**:
  - `base_anth_rate_light` ✅
  - `base_anth_rate_dark` ✅
  - `V_max_anth` ✅
  - `K_stress_anth` ✅ (v6.7 改用 Stress 半飽和常數)
  - `k_deg` ✅

**未發現遺留 E_stress 相關參數** ✅

---

## 潛在問題與建議

### 1. 設為 0 的參數

#### anth_carbon_cost = 0.0
- **影響**: 花青素合成無碳成本
- **原因**: 非零值會過度抑制生長（已測試）
- **建議**:
  - ✅ **保留**：此參數記錄了設計決策（嘗試過但不可行）
  - 保留註解說明原因（已完成，見 params_config.py:93-95）

### 2. 值為 1.0 的參數

#### par_conversion_factor = 1.0
- **影響**: UVA-PAR 轉換無放大效應 (1:1)
- **歷史**: v6.0 移除放大效應（v5.x 為 3.0）
- **建議**:
  - ✅ **保留**：提供未來調整彈性
  - v6.0 結論: Sun 原始模型足以解釋 UVA 鮮重促進效應

#### ldmc_stress_sensitivity = 1.0
- **影響**: LDMC 對 Stress 的敏感度基準值
- **建議**:
  - ✅ **保留**：正常的校準參數

---

## 結論

### 參數清理狀態: ✅ 無需清理

所有 28 個參數均已在程式碼中使用。沒有冗餘或遺留參數。

### v6.7 升級確認: ✅ 完成

- E_stress 狀態變量已正確移除
- 花青素機制已正確改為 FW-based
- 所有花青素相關參數均正確使用
- 無遺留的 LAI-based 或 E_stress-based 機制

### 特殊參數說明

兩個參數設為"無效"值（0.0 或 1.0），但均有充分的設計理由：

1. **anth_carbon_cost = 0.0**: 記錄設計決策（嘗試過但不可行）
2. **par_conversion_factor = 1.0**: 保留未來調整彈性

**建議**: 兩者均保留，維持當前設計。

---

## 附錄: 參數使用統計

### 按類別統計

| 類別 | 參數數 | 已使用 | 未使用 |
|------|--------|--------|--------|
| Sun 基礎模型 | 1 | 1 | 0 |
| UVA-PAR 轉換 | 1 | 1 | 0 |
| Stress 損傷與修復 | 16 | 16 | 0 |
| 碳修復 | 4 | 4 | 0 |
| 花青素 | 6 | 6 | 0 |
| LDMC | 4 | 4 | 0 |
| 其他 | 1 | 1 | 0 |
| **總計** | **33** | **33** | **0** |

### 按使用頻率統計

| 使用次數 | 參數數 | 參數列表 |
|----------|--------|----------|
| 1次 | 32 | 大部分參數 |
| 2次 | 1 | dw_fw_ratio_base |
| 3次以上 | 0 | - |

**說明**: `dw_fw_ratio_base` 在兩處使用：
1. `calculate_dynamic_dw_fw_ratio` 函數（計算動態 DW:FW）
2. 初始條件設定（計算初始乾重）

---

**審計完成**: 2025-12-23
**審計者**: Claude Code
**模型版本**: v6.7 (FW-based Anthocyanin)
