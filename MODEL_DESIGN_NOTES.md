# 萵苣UVA模型設計筆記 (Model Design Notes)

**重要：每次新聊天請先閱讀此文檔**
**最後更新**: 2025-12-11

---

## 最新驗證結果 (v5.5 - 夜間 UVA-PAR 修正)

### FW 預測結果
```
Treatment | FW_sim | FW_exp | FW_Err | Stress
---------------------------------------------
       CK |   87.6 |   87.0 |  +0.7% |    0.0
     L6D6 |   88.3 |   91.4 |  -3.3% |    2.0
   L6D6-N |   81.4 |   80.8 |  +0.8% |    4.8  <- 大幅改善！
    H12D3 |   62.5 |   60.6 |  +3.1% |   36.1
   VL3D12 |   63.6 |   67.0 |  -5.0% |    6.5
    L6D12 |   61.0 |   60.4 |  +0.9% |   11.2

Mean |FW_Err|: 2.3%  |  Max: 5.0%
```

### Anth 預測結果
```
Treatment | Anth_sim | Anth_exp | Anth_Err
------------------------------------------
       CK |    45.7  |    43.3  |   +5.6%
     L6D6 |    45.3  |    49.4  |   -8.3%
   L6D6-N |    49.2  |    49.3  |   -0.1%  <- 大幅改善！
    H12D3 |    68.7  |    65.1  |   +5.7%
   VL3D12 |    64.6  |    48.2  |  +34.1%
    L6D12 |    72.5  |    51.8  |  +40.0%

Mean |Anth_Err|: 15.6%  |  Max: 40.0%
```

**備註**:
- L6D6-N 誤差從 4.4%/0.5% 改善至 0.8%/0.1%（FW/Anth）
- 夜間 UVA 現在正確貢獻光合作用
- VL3D12、L6D12 花青素仍偏高，為模型結構性限制

---

## 1. 模型概述

本模型整合了 Sun et al. (2025) 的萵苣生長模型與 UVA 光處理效應，用於預測不同 UVA 處理方案對萵苣鮮重 (FW) 和花青素含量 (Anth) 的影響。

### 1.1 狀態變量 (5個) - v5.0 簡化版

| 變量 | 符號 | 單位 | 描述 |
|------|------|------|------|
| 乾物質密度 | X_d | kg/m² | 植物乾重 |
| 碳緩衝池 | C_buf | kg C/m² | 光合產物暫存 |
| 葉面積指數 | LAI | - | 葉面積/地面積 |
| 花青素含量 | Anth | g/m² | 花青素總量 |
| 脅迫程度 | Stress | - | 統一的損傷-修復狀態 |

---

## 2. 五個微分方程 (ODE)

### 2.1 ODE 1: 乾物質密度 (dX_d/dt)

**來源**: Sun et al. (2025) 基礎模型 + Stress 抑制

```
dX_d/dt = dX_d_base × (1 - β_photo × Stress/(K_stress + Stress))
```

- `dX_d_base`: 來自 Sun 模型的基礎光合同化
- `β_photo = 0.70`: Stress 對光合的最大抑制比例
- `K_stress = 5.0`: 半飽和 Stress 值

### 2.2 ODE 2: 碳緩衝池 (dC_buf/dt)

**來源**: Sun 模型 - 修復碳消耗

```
dC_buf/dt = dC_buf_base - repair_rate × repair_carbon_cost
```

- 修復損傷需要消耗碳池中的碳
- `repair_carbon_cost = 1e-6`: 每單位修復消耗的碳

### 2.3 ODE 3: 葉面積指數 (dLAI/dt)

**來源**: Sun 模型 + Stress 抑制

```
dLAI/dt = dLAI_base × (1 - β_LAI × Stress/(K_stress + Stress))
```

- `β_LAI = 0.70`: Stress 對 LAI 生長的最大抑制比例

### 2.4 ODE 4: 花青素 (dAnth/dt) - v5.2 更新

**來源**: 自創機制 + 閾值 + Michaelis-Menten 飽和

```
effective_stress = max(0, Stress - threshold)
uva_induced = V_max × effective_stress / (K_m + effective_stress)
dAnth/dt = base_synthesis + uva_induced - k_deg × Anth
```

| 參數 | 值 | 說明 |
|------|-----|------|
| base_anth_rate_light | 3.85e-10 kg/m²/s | 日間基礎合成率 |
| base_anth_rate_dark | 1.92e-10 kg/m²/s | 夜間基礎合成率 (50%) |
| V_max_anth | 2.5e-10 kg/m²/s | 最大 Stress 誘導合成率 |
| K_m_anth | 30.0 | 半飽和 Stress 值 |
| stress_threshold_anth | 15.0 | Stress 閾值 |
| k_deg | 2.5e-6 /s | 降解速率 (半衰期 ~3.2 天) |

**機制說明**:
1. 只有 Stress 超過閾值 (15) 時才會誘導額外花青素
2. 使用飽和函數避免高 Stress 過度響應
3. 這解釋了 H12D3 (高 Stress) 花青素顯著增加，而低 Stress 處理組維持基礎水平

### 2.5 ODE 5: Stress (dStress/dt) - 核心自創機制

```
dStress/dt = damage_rate - repair_rate

損傷項:
damage = k_damage × I_UVA × vulnerability × intraday_factor × nonlinear × circadian

修復項:
repair = k_repair × Stress × carbon_availability
```

| 參數 | 值 | 說明 |
|------|-----|------|
| k_damage | 3.5e-6 | 損傷係數 (校準值) |
| k_repair | 1.0e-5 | 修復係數 |
| carbon_availability | 0.5~1.0 | 碳依賴修復能力 |

---

## 3. 自創函數介紹

### 3.1 LAI 脆弱性函數 (vulnerability)

**位置**: `simulate_uva_model.py` 第 504-517 行

```
vulnerability = cap × (LAI_ref/LAI)^n / (cap + (LAI_ref/LAI)^n)
```

| 參數 | 值 | 說明 |
|------|-----|------|
| LAI_ref | 7.5 | 參考 LAI (L6D6 開始時) |
| n_vuln | 7 | 脆弱性指數 |
| cap_vuln | 100 | 上限值 |

**機制說明**: 低 LAI (幼嫩植物) 對 UVA 更脆弱，這解釋了 VL3D12/L6D12（早期照射）比 L6D6（晚期照射）受損更嚴重。

**典型值**:
- LAI=5.21 (day 23): vulnerability ≈ 5.3
- LAI=7.49 (day 29): vulnerability ≈ 0.84
- LAI=9.50 (day 35): vulnerability ≈ 0.17

### 3.2 當日照射時數非線性因子 (intraday_factor)

**位置**: `simulate_uva_model.py` 第 519-556 行

```
intraday_factor = 1 + k × softplus(hours_elapsed - 6)^m

softplus(x) = ln(1 + e^(sharpness × x)) / sharpness
```

| 參數 | 值 | 說明 |
|------|-----|------|
| k_intraday | 1.5 | 非線性放大係數 (校準值) |
| m_intraday | 2.0 | 指數 |
| sharpness | 3.0 | softplus 銳度 |

**機制說明**: 超過 6 小時後，修復機制飽和，損傷效率急劇上升。這解釋了 H12D3（每天12h）造成更嚴重的損傷。

**典型值**:
- 3h/day: intraday_factor ≈ 1.0
- 6h/day: intraday_factor ≈ 1.0
- 12h/day: intraday_factor ≈ 6.4

### 3.3 Stress 非線性累積因子 (nonlinear_factor)

**位置**: `simulate_uva_model.py` 第 295-297 行

```
nonlinear_factor = 1 + α × Stress / (K_nonlinear + Stress)
```

| 參數 | 值 | 說明 |
|------|-----|------|
| α (stress_nonlinear_coeff) | 1.5 | 非線性放大係數 |
| K_nonlinear | 3.0 | 半飽和常數 |

**機制說明**: ROS 級聯效應 - 已有損傷會加速新損傷（正反饋）。

**注意**: v5.4 版已移除「指數型天數累積」功能。D12 組的額外損傷主要由 LAI 脆弱性機制解釋（早期照射時 LAI 較低，vulnerability 較高）。

### 3.4 夜間節律抑制 (circadian_penalty) - v5.5 更新

**位置**: `simulate_uva_model.py` 第 304-308 行

```
circadian_penalty = 2.0 (v5.5 校準值: 補償夜間 UVA-PAR 光合貢獻)
```

**機制說明**:
- v5.5 修正了夜間 UVA-PAR 光合貢獻機制（使用 I_override 繞過 Sun 模型日夜判斷）
- 夜間 UVA 現在可以貢獻光合作用，但會增加 Stress
- `circadian_disruption_factor = 2.0` 補償夜間光合帶來的額外生長

**校準結果**:
| CDF | L6D6-N FW | 誤差 | Stress |
|-----|-----------|------|--------|
| 1.0 | 90.7g | +12.3% | 1.6 |
| 2.0 | 81.4g | +0.8% | 4.8 |
| 3.0 | 74.3g | -8.1% | 9.0 |

### 3.5 動態 DW:FW 比例 (calculate_dynamic_dw_fw_ratio)

**位置**: `simulate_uva_model.py` 第 304-336 行

```
ratio = base × (1 + sensitivity × Stress / (K_ldmc + Stress))
```

| 參數 | 值 | 說明 |
|------|-----|------|
| base | 0.05 (5%) | 基礎 DW:FW |
| sensitivity | 1.0 | LDMC 敏感度 |
| K_ldmc | 50.0 | 半飽和 Stress |
| max | 0.12 (12%) | 最大 DW:FW |

**機制說明**: UV 脅迫使葉片變厚變乾（LDMC 效應），來自 Qian et al. (2021)。

---

## 4. 核心機制總結

### 4.1 UVA-PAR 增益效應 (即時)

```
I_effective = I_base + I_UVA × par_conversion_factor
```

- `par_conversion_factor = 4.3`
- 只在日間 LED 開啟時有效
- 解釋了 L6D6 的 FW 能接近或超過 CK

### 4.2 Stress 對 FW 的影響路徑

1. **生長抑制**: Stress → 抑制 dXd_dt 和 dLAI_dt
2. **LDMC 效應**: Stress → DW/FW 比例升高 → FW 降低

---

## 5. 已校準參數配置 (v5.5 最終版)

```python
# simulate_uva_model.py 中的關鍵參數

# 基礎模型
c_alpha = 0.54  # 校準值，調整基礎光合以符合 CK

# UVA-PAR 轉換
par_conversion_factor = 3.0  # UVA 對 PAR 的增益

# Stress 損傷
stress_damage_coeff = 3.5e-6  # 損傷係數 (校準值)
stress_repair_coeff = 1.0e-5  # 修復係數

# 夜間節律 (v5.5 更新)
circadian_disruption_factor = 2.0  # 夜間損傷加成 (v5.5 校準: 補償夜間 UVA-PAR 光合貢獻)

# LAI 脆弱性 (解釋 D12 組早期照射的額外損傷)
LAI_ref_vuln = 7.5
n_vuln = 7
cap_vuln = 100.0  # 上限

# 當日照射時數非線性 (解釋 H12D3 的額外損傷)
k_intraday = 1.5  # 校準值 (修正 H12D3)
m_intraday = 2.0
sharpness_intraday = 3.0

# Stress 非線性累積 (ROS 級聯效應)
stress_nonlinear_coeff = 1.5
K_nonlinear = 3.0

# Stress 對生長抑制
stress_photosynthesis_inhibition = 0.70
stress_lai_inhibition = 0.70
K_stress = 5.0

# LDMC
K_ldmc = 50.0
ldmc_stress_sensitivity = 1.0
```

---

## 6. 實驗處理組設計

| 處理代碼 | 描述 | UVA時數/天 | 總天數 | 開始日 |
|----------|------|------------|--------|--------|
| CK | 對照組 | 0 | - | - |
| L6D6 | 低劑量日間 | 6h (10:00-16:00) | 6天 | Day 29 |
| L6D6-N | 低劑量夜間 | 6h (22:00-04:00) | 6天 | Day 29 |
| H12D3 | 高劑量脅迫 | 12h (06:00-18:00) | 3天 | Day 32 |
| VL3D12 | 極低劑量 | 3h (10:00-13:00) | 12天 | Day 23 |
| L6D12 | 低劑量長期 | 6h (10:00-16:00) | 12天 | Day 23 |

---

## 7. ODE求解器設定 (重要!)

```python
method = 'RK45'  # 不要改! LSODA 在不同步長下結果不穩定
max_step = 60    # 秒
```

---

## 8. 文件結構

| 文件 | 說明 |
|------|------|
| `simulate_uva_model.py` | 主模型 v5.0 (5狀態變量 + Stress 機制) - 已校準 |
| `model_config.py` | 共用設定模組 |
| `lettuce_uva_carbon_complete_model.py` | 基礎 Sun 模型 |
| `CLAUDE.md` | Claude 工作守則 |
| `MODEL_DESIGN_NOTES.md` | 本文件 |
| `HANDOFF_STATUS.md` | 交接狀態 |

---

## 9. 參考文獻

### 基礎模型
- Sun et al. (2025) - 萵苣生長模型基礎

### UVA-PAR 效應
- Verdaguer et al. (2017) - UVA 對光合作用的促進效應
- Kataria et al. (2014) - UV-A effects on photosynthesis

### Stress 與損傷
- Hideg et al. (2013) Trends Plant Sci. - UV 損傷與 ROS
- Frohnmeyer & Staiger (2003) - UV signaling and oxidative stress
- Foyer & Noctor (2005) - Oxidant and antioxidant signalling

### 修復機制
- Zhu et al. (2018) - Repair mechanisms require carbon/energy

### 節律抑制
- Harmer (2009) Annu. Rev. Plant Biol. - Circadian rhythm
- Covington et al. (2008) Genome Biol. - Night light disruption

### LDMC 效應
- Qian et al. (2021) Plant Physiol. DOI:10.1093/plphys/kiab262

### 花青素
- Winkel-Shirley (2002) - Anthocyanin biosynthesis
- Gould (2004) - Anthocyanin as stress protectant

---

**模型校準完成 (v5.5)！Mean |FW_Err| = 2.3%，Max = 5.0%，L6D6-N 誤差大幅改善至 <1%**
