# 參數使用映射表 (Parameter Usage Map)

**日期**: 2025-12-23
**目的**: 快速查找每個參數在代碼中的定義和使用位置

---

## 使用說明

本文件提供每個參數的完整追蹤路徑：

1. **定義位置**: params_config.py 中的定義行
2. **載入位置**: simulate_uva_model.py 中 UVAParams.__init__ 的載入行
3. **使用位置**: simulate_uva_model.py 中實際使用該參數的行號

---

## 參數映射表

### 1. Sun 基礎模型

#### c_alpha
- **定義**: params_config.py:12
- **載入**: simulate_uva_model.py:84
- **使用**: 繼承自 BaseSunParams，在 sun_derivatives_final 中使用
- **用途**: 光合效率係數（校準至 CK ≈ 87.0g）
- **值**: 0.548

---

### 2. UVA-PAR 轉換

#### par_conversion_factor
- **定義**: params_config.py:22
- **載入**: simulate_uva_model.py:89
- **使用**: simulate_uva_model.py:274
  ```python
  I_gain_par = I_UVA * p.par_conversion_factor
  ```
- **用途**: UVA 轉換為等效 PAR 的係數
- **值**: 1.0 (v6.0: 移除放大效應)

---

### 3. Stress 損傷與修復

#### stress_damage_coeff
- **定義**: params_config.py:32
- **載入**: simulate_uva_model.py:94
- **使用**: simulate_uva_model.py:355-356
  ```python
  damage_rate = p.stress_damage_coeff * I_UVA * vulnerability * \
                intraday_factor * nonlinear_factor * circadian_penalty
  ```
- **用途**: Stress 損傷係數 [1/(W/m²·s)]
- **值**: 6.6e-7

#### stress_repair_coeff
- **定義**: params_config.py:33
- **載入**: simulate_uva_model.py:95
- **使用**: simulate_uva_model.py:366
  ```python
  repair_rate = p.stress_repair_coeff * Stress * repair_capacity
  ```
- **用途**: Stress 修復係數 [1/s]
- **值**: 1.0e-5

#### stress_nonlinear_coeff
- **定義**: params_config.py:34
- **載入**: simulate_uva_model.py:96
- **使用**: simulate_uva_model.py:346
  ```python
  nonlinear_factor = 1.0 + p.stress_nonlinear_coeff * Stress / (p.K_nonlinear + Stress + 1e-9)
  ```
- **用途**: Stress 累積非線性係數（ROS 級聯效應）
- **值**: 8.0

#### K_nonlinear
- **定義**: params_config.py:35
- **載入**: simulate_uva_model.py:97
- **使用**: simulate_uva_model.py:346
- **用途**: Stress 非線性半飽和常數
- **值**: 0.8

---

### 4. LAI 脆弱性

#### LAI_ref_vuln
- **定義**: params_config.py:38
- **載入**: simulate_uva_model.py:100
- **使用**: simulate_uva_model.py:297, 301
  ```python
  LAI_ref_vuln = p.LAI_ref_vuln
  base_vuln = (LAI_ref_vuln / LAI) ** n_vuln
  ```
- **用途**: 脆弱性參考 LAI（幼苗更易受損）
- **值**: 6.5

#### n_vuln
- **定義**: params_config.py:39
- **載入**: simulate_uva_model.py:101
- **使用**: simulate_uva_model.py:298, 301
  ```python
  n_vuln = p.n_vuln
  base_vuln = (LAI_ref_vuln / LAI) ** n_vuln
  ```
- **用途**: 脆弱性指數
- **值**: 8

#### cap_vuln
- **定義**: params_config.py:40
- **載入**: simulate_uva_model.py:102
- **使用**: simulate_uva_model.py:299, 302
  ```python
  cap_vuln = p.cap_vuln
  vulnerability = cap_vuln * base_vuln / (cap_vuln + base_vuln)
  ```
- **用途**: 脆弱性上限
- **值**: 100.0

---

### 5. 日內能量非線性

#### E_50
- **定義**: params_config.py:44
- **載入**: simulate_uva_model.py:105
- **使用**: simulate_uva_model.py:338
  ```python
  normalized_E = (E_elapsed - p.E_50) / p.E_scale
  ```
- **用途**: 半飽和能量 [kJ/m²] (≈ 6h @ 22 W/m²)
- **值**: 475.2

#### E_scale
- **定義**: params_config.py:45
- **載入**: simulate_uva_model.py:106
- **使用**: simulate_uva_model.py:338
  ```python
  normalized_E = (E_elapsed - p.E_50) / p.E_scale
  ```
- **用途**: 能量尺度 [kJ/m²] (≈ 3h @ 22 W/m²)
- **值**: 237.6

#### k_intraday
- **定義**: params_config.py:46
- **載入**: simulate_uva_model.py:107
- **使用**: simulate_uva_model.py:342
  ```python
  intraday_factor = 1.0 + p.k_intraday * (excess_normalized ** p.m_intraday)
  ```
- **用途**: 日內非線性放大係數（H12D3 關鍵參數）
- **值**: 49.0

#### m_intraday
- **定義**: params_config.py:47
- **載入**: simulate_uva_model.py:108
- **使用**: simulate_uva_model.py:342
  ```python
  intraday_factor = 1.0 + p.k_intraday * (excess_normalized ** p.m_intraday)
  ```
- **用途**: 日內非線性指數
- **值**: 2.0

#### sharpness_intraday
- **定義**: params_config.py:48
- **載入**: simulate_uva_model.py:109
- **使用**: simulate_uva_model.py:339
  ```python
  excess_normalized = softplus(normalized_E, p.sharpness_intraday)
  ```
- **用途**: softplus 銳度（控制閾值平滑度）
- **值**: 3.0

---

### 6. 夜間節律損傷

#### circadian_disruption_factor
- **定義**: params_config.py:51
- **載入**: simulate_uva_model.py:112
- **使用**: simulate_uva_model.py:350
  ```python
  circadian_penalty = p.circadian_disruption_factor
  ```
- **用途**: 夜間 UVA 損傷加成（節律破壞）
- **值**: 3.0

---

### 7. Stress 對生長的抑制

#### stress_photosynthesis_inhibition
- **定義**: params_config.py:54
- **載入**: simulate_uva_model.py:115
- **使用**: simulate_uva_model.py:380
  ```python
  xd_reduction = p.stress_photosynthesis_inhibition * stress_inhibition
  ```
- **用途**: Stress 對光合作用的最大抑制比例
- **值**: 0.66

#### stress_lai_inhibition
- **定義**: params_config.py:55
- **載入**: simulate_uva_model.py:116
- **使用**: simulate_uva_model.py:381
  ```python
  lai_reduction = p.stress_lai_inhibition * stress_inhibition
  ```
- **用途**: Stress 對 LAI 生長的最大抑制比例
- **值**: 0.66

#### K_stress
- **定義**: params_config.py:56
- **載入**: simulate_uva_model.py:117
- **使用**: simulate_uva_model.py:378
  ```python
  stress_inhibition = Stress / (p.K_stress + Stress + 1e-9)
  ```
- **用途**: 生長抑制半飽和 Stress 值
- **值**: 1.9

---

### 8. 碳修復

#### base_repair_capacity
- **定義**: params_config.py:63
- **載入**: simulate_uva_model.py:122
- **使用**: simulate_uva_model.py:362
  ```python
  repair_capacity = p.base_repair_capacity + \
                    p.carbon_repair_bonus * C_buf_positive / (p.K_carbon + C_buf_positive + 1e-9)
  ```
- **用途**: 基礎修復能力
- **值**: 0.5

#### carbon_repair_bonus
- **定義**: params_config.py:64
- **載入**: simulate_uva_model.py:123
- **使用**: simulate_uva_model.py:363
- **用途**: 碳池對修復能力的加成
- **值**: 0.5

#### K_carbon
- **定義**: params_config.py:65
- **載入**: simulate_uva_model.py:124
- **使用**: simulate_uva_model.py:363
- **用途**: 碳修復半飽和常數 [kg C/m²]
- **值**: 0.001

#### repair_carbon_cost
- **定義**: params_config.py:66
- **載入**: simulate_uva_model.py:125
- **使用**: simulate_uva_model.py:437
  ```python
  repair_carbon_consumption = repair_rate * p.repair_carbon_cost
  ```
- **用途**: 每單位修復消耗的碳 [kg C / Stress]
- **值**: 1.0e-6

---

### 9. 花青素 (v6.7: FW-based 機制)

#### base_anth_rate_light
- **定義**: params_config.py:74
- **載入**: simulate_uva_model.py:130
- **使用**: simulate_uva_model.py:421
  ```python
  base_synthesis_per_fw = day_weight * p.base_anth_rate_light + (1 - day_weight) * p.base_anth_rate_dark
  ```
- **用途**: 日間基礎合成率 [kg Anth / (kg FW · s)]
- **值**: 2.0e-10

#### base_anth_rate_dark
- **定義**: params_config.py:75
- **載入**: simulate_uva_model.py:131
- **使用**: simulate_uva_model.py:421
- **用途**: 夜間基礎合成率 [kg Anth / (kg FW · s)]
- **值**: 1.0e-10

#### V_max_anth
- **定義**: params_config.py:86
- **載入**: simulate_uva_model.py:132
- **使用**: simulate_uva_model.py:425
  ```python
  stress_induced_per_fw = p.V_max_anth * Stress / (p.K_stress_anth + Stress + 1e-12)
  ```
- **用途**: 最大 Stress 誘導合成率 [kg Anth / (kg FW · s)]
- **值**: 2.35e-11

#### K_stress_anth
- **定義**: params_config.py:87
- **載入**: simulate_uva_model.py:133
- **使用**: simulate_uva_model.py:425
- **用途**: 花青素合成的 Stress 半飽和常數
- **值**: 0.30
- **備註**: v6.7 改用 Stress (而非 E_stress)

#### k_deg
- **定義**: params_config.py:90
- **載入**: simulate_uva_model.py:134
- **使用**: simulate_uva_model.py:431
  ```python
  degradation = p.k_deg * Anth
  ```
- **用途**: 花青素降解速率 [1/s]
- **值**: 3.02e-6

#### anth_carbon_cost
- **定義**: params_config.py:93
- **載入**: simulate_uva_model.py:135
- **使用**: simulate_uva_model.py:438
  ```python
  anth_carbon_consumption = synthesis_rate * p.anth_carbon_cost
  ```
- **用途**: 每單位花青素合成消耗的碳 [kg C / kg Anth]
- **值**: 0.0
- **備註**: 設為 0 因為非零值會過度抑制生長

---

### 10. LDMC (葉乾物質含量)

#### dw_fw_ratio_base
- **定義**: params_config.py:102
- **載入**: simulate_uva_model.py:140
- **使用**:
  - simulate_uva_model.py:167 (calculate_dynamic_dw_fw_ratio)
  - simulate_uva_model.py:480 (初始條件設定)
  ```python
  ratio = p.dw_fw_ratio_base * (1.0 + stress_effect)
  dw_init_g = fw_init_g * p.dw_fw_ratio_base
  ```
- **用途**: 基礎 DW:FW 比例 (5%)
- **值**: 0.05

#### ldmc_stress_sensitivity
- **定義**: params_config.py:103
- **載入**: simulate_uva_model.py:141
- **使用**: simulate_uva_model.py:166
  ```python
  stress_effect = p.ldmc_stress_sensitivity * Stress / (p.K_ldmc + Stress + 1e-9)
  ```
- **用途**: LDMC 對 Stress 的敏感度
- **值**: 1.0

#### K_ldmc
- **定義**: params_config.py:104
- **載入**: simulate_uva_model.py:142
- **使用**: simulate_uva_model.py:166
- **用途**: LDMC 半飽和 Stress 值
- **值**: 50.0

#### dw_fw_ratio_max
- **定義**: params_config.py:105
- **載入**: simulate_uva_model.py:143
- **使用**: simulate_uva_model.py:168
  ```python
  return min(ratio, p.dw_fw_ratio_max)
  ```
- **用途**: 最大 DW:FW 比例 (12%)
- **值**: 0.12

---

### 11. 其他

#### transplant_day
- **定義**: params_config.py:112
- **載入**: simulate_uva_model.py:148
- **使用**: simulate_uva_model.py:219
  ```python
  day_from_transplant = day_from_sowing - p.transplant_day
  ```
- **用途**: 移植日 (播種後第14天)
- **值**: 14

---

## 參數使用統計

### 按使用次數分類

| 使用次數 | 參數列表 | 數量 |
|----------|----------|------|
| 1次 | 大部分參數 | 32 |
| 2次 | dw_fw_ratio_base | 1 |
| **總計** | | **33** |

### 多次使用的參數

#### dw_fw_ratio_base (2次使用)

1. **calculate_dynamic_dw_fw_ratio** (simulate_uva_model.py:167)
   - 用途: 計算動態 DW:FW 比例（基準值）

2. **初始條件設定** (simulate_uva_model.py:480)
   - 用途: 計算初始乾重

---

## v6.7 機制驗證

### 花青素機制 (FW-based)

**合成公式**:
```python
# 第 1 步: 計算單位 FW 合成率
base_synthesis_per_fw = day_weight * p.base_anth_rate_light + (1 - day_weight) * p.base_anth_rate_dark
stress_induced_per_fw = p.V_max_anth * Stress / (p.K_stress_anth + Stress + 1e-12)

# 第 2 步: 乘以總 FW
synthesis_rate = FW_kg_m2 * (base_synthesis_per_fw + stress_induced_per_fw)

# 第 3 步: 扣除降解
degradation = p.k_deg * Anth
dAnth_dt = synthesis_rate - degradation
```

**使用參數**:
- base_anth_rate_light ✅
- base_anth_rate_dark ✅
- V_max_anth ✅
- K_stress_anth ✅
- k_deg ✅

**移除項目**:
- E_stress 狀態變量 ✅ (v6.7 移除)
- LAI-based 合成 ✅ (v6.6 → v6.7 改為 FW-based)

---

## 快速查找索引

### 按字母順序

| 參數名稱 | 定義行 | 載入行 | 使用行 |
|----------|--------|--------|--------|
| anth_carbon_cost | 93 | 135 | 438 |
| base_anth_rate_dark | 75 | 131 | 421 |
| base_anth_rate_light | 74 | 130 | 421 |
| base_repair_capacity | 63 | 122 | 362 |
| c_alpha | 12 | 84 | (繼承) |
| cap_vuln | 40 | 102 | 299, 302 |
| carbon_repair_bonus | 64 | 123 | 363 |
| circadian_disruption_factor | 51 | 112 | 350 |
| dw_fw_ratio_base | 102 | 140 | 167, 480 |
| dw_fw_ratio_max | 105 | 143 | 168 |
| E_50 | 44 | 105 | 338 |
| E_scale | 45 | 106 | 338 |
| K_carbon | 65 | 124 | 363 |
| K_ldmc | 104 | 142 | 166 |
| K_nonlinear | 35 | 97 | 346 |
| K_stress | 56 | 117 | 378 |
| K_stress_anth | 87 | 133 | 425 |
| k_deg | 90 | 134 | 431 |
| k_intraday | 46 | 107 | 342 |
| LAI_ref_vuln | 38 | 100 | 297, 301 |
| ldmc_stress_sensitivity | 103 | 141 | 166 |
| m_intraday | 47 | 108 | 342 |
| n_vuln | 39 | 101 | 298, 301 |
| par_conversion_factor | 22 | 89 | 274 |
| repair_carbon_cost | 66 | 125 | 437 |
| sharpness_intraday | 48 | 109 | 339 |
| stress_damage_coeff | 32 | 94 | 355 |
| stress_lai_inhibition | 55 | 116 | 381 |
| stress_nonlinear_coeff | 34 | 96 | 346 |
| stress_photosynthesis_inhibition | 54 | 115 | 380 |
| stress_repair_coeff | 33 | 95 | 366 |
| transplant_day | 112 | 148 | 219 |
| V_max_anth | 86 | 132 | 425 |

---

## 結論

✅ **所有 33 個參數均已正確映射**

- 每個參數都有明確的定義、載入和使用位置
- 無遺漏、無冗餘
- v6.7 花青素機制升級完成
- 代碼結構清晰，易於維護

---

**文件更新**: 2025-12-23
**審計工具**: Claude Code
**版本**: v6.7
