# 交接狀態 (HANDOFF_STATUS)

**最後更新**: 2025-12-22 (v6.3 - 參數外部化 + PAR修正 - 5/6組達標✅)

---

## 當前進度 (v6.3)

### ✅ v6.3 校準完成 - 5/6 組達標，L6D6 > CK 確認！

**核心成就**:
1. ✅ **參數完全外部化** - 所有參數移至 [params_config.py](params_config.py)
2. ✅ **修正 PAR 計算** - 移除 par_conversion_factor 放大 (3.0 → 1.0)
3. ✅ **修正時間計算 Bug** - day_from_sowing 計算錯誤 (v6.1)
4. ✅ **L6D6 > CK** - 符合物理預期 (UVA 提供額外 PAR)
5. ✅ **5/6 組 < 5%** - 83.3% 達標率

### ⚠️ v5.9 基準結果 (已被 v6.0-v6.3 取代)

**版本演進**:
- v5.9: 花青素機制改用 Stress 累積能量 (E_stress)
- v6.0: 移除 PAR 放大效應 (par_conversion_factor: 3.0 → 1.0)
- v6.1: 修復時間計算 Bug (day_from_sowing)
- v6.2: 提高抑制靈敏度
- v6.3: L6D6 優先穩定 + 機制分離調整

**當前狀態**: ✅ v6.3 校準完成

**v6.3 最終預測結果**:

```
Treatment   FW_sim   FW_exp  FW_Err   狀態
-------------------------------------------
CK          87.8g    87.0g   +1.0%    ✅ < 5%
L6D6        89.3g    91.4g   -2.3%    ✅ < 5% ⭐ > CK
L6D6-N      81.1g    80.8g   +0.3%    ✅ < 5% (幾乎完美!)
VL3D12      64.9g    67.0g   -3.1%    ✅ < 5%
L6D12       59.3g    60.4g   -1.9%    ✅ < 5%
H12D3       90.8g    60.6g  +49.9%    ❌ > 10% (模型固有限制)
```

**統計摘要**:
- **< 5% 組數**: 5/6 (83.3%) ✅
- **平均誤差**: 9.7% ✅
- **中位數誤差**: 1.9% ✅
- **L6D6 > CK**: ✅ 符合物理預期

---

## v6.0-v6.3 主要變更 (2025-12-22)

### v6.0: 參數外部化 + PAR 修正

**關鍵變更**:
1. **參數完全外部化** - 創建 [params_config.py](params_config.py)
   - 所有參數集中管理，不再散落在類別中
   - 便於批次調整和版本控制

2. **移除 PAR 放大效應**
   - `par_conversion_factor: 3.0 → 1.0`
   - 理由: I_UVA = 22 W/m² 已是等效短波輻射，不需再放大
   - Sun 原始模型已能正確模擬 UVA 效應

3. **校準 c_alpha**
   - `c_alpha: 0.68 → 0.548`
   - 確保 CK ≈ 87.0g

### v6.1: 修復時間計算 Bug ⭐⭐⭐

**發現的 Bug**:
```python
# v6.0 (錯誤)
day_from_transplant = t / 86400.0
day_from_sowing = day_from_transplant + p.transplant_day  # BUG!

# v6.1 (正確)
day_from_sowing = t / 86400.0  # 播種後天數 (絕對時間)
day_from_transplant = day_from_sowing - p.transplant_day  # 移植後天數
```

**影響**:
- L6D6 在 Day 15 就被誤判為 Day 29 (UVA 照射開始!)
- 導致 Stress 提前累積，生長被抑制
- L6D6 誤差從 -24.3% 改善到 -2.5% (+21.8%) ✅

### v6.2: 提高抑制靈敏度

**策略**: 提高 Stress 對生長的抑制效果
- `stress_photosynthesis_inhibition: 0.60 → 0.68`
- `K_stress: 2.5 → 1.8`

**結果**: 4/6 組 < 5%

### v6.3: L6D6 優先穩定 + 機制分離調整

**核心策略**: 確保 **L6D6 > CK** (因為 UVA 提供額外 22 W/m² PAR)

**參數調整**:
```python
'stress_damage_coeff': 1.45e-6 → 0.70e-6,  # -52% (降低 L6D6 Stress)
'circadian_disruption_factor': 2.0 → 3.2,  # +60% (僅調整 L6D6-N)
'stress_nonlinear_coeff': 5.5 → 8.0,       # +45% (改善長期組)
'K_nonlinear': 1.0 → 0.8,                  # -20% (更早觸發)
```

**結果**: ✅ **5/6 組 < 5%**，L6D6 > CK 確認

---

## H12D3 問題分析 (v6.3)

### 為什麼 H12D3 無法達標？

**實驗配置**:
- H12D3: 22 W/m², 12h/day, **只照 3 天** (Day 32-35)
- 目標: 60.6g (vs CK 87.0g，需要 **-30% 抑制**)
- 預測: 90.8g (**+49.9%**)

**根本限制**:

1. **時間太短 (3 天)**
   - Stress 累積模型需要時間建立正反饋
   - 即使 12h/day 照射，3 天不足以讓 Stress 達到所需水平

2. **LAI 已經很高 (Day 32 時 LAI ≈ 4.5)**
   - 成熟植株，脆弱性低 (vulnerability ≈ 14 vs 目標 100)
   - 基礎損傷率太低

3. **數學限制**
   - 非線性累積是 Michaelis-Menten 型飽和函數
   - 只有當 Stress >> K_nonlinear 時才有強放大
   - H12D3 只有 3 天，Stress 無法進入快速累積階段

### 唯一解決方案：即時損傷機制 (Acute Damage)

**概念**: 高強度長時照射 (> 8-10 h/day) 直接抑制光合作用，不依賴 Stress 累積

**實現建議**:
```python
def calculate_growth_inhibition(I_UVA, Stress, hours_per_day, p):
    # 累積效應 (原有機制)
    cumulative_inhibition = p.stress_photosynthesis_inhibition * Stress / (p.K_stress + Stress)

    # 即時效應 (新增機制) - 只在長時照射時觸發
    if hours_per_day > 8:
        acute_dose = I_UVA * (hours_per_day - 8) / 100
        acute_inhibition = p.acute_inhibition_coeff * acute_dose / (p.K_acute + acute_dose)
    else:
        acute_inhibition = 0.0

    # 總抑制 (取較大值)
    total_inhibition = max(cumulative_inhibition, acute_inhibition)
    return total_inhibition
```

**優點**:
- ✅ 能夠捕捉短期高強度效應 (H12D3: 12h/day × 3 days)
- ✅ 不影響其他處理組 (L6D6: 6h/day, L6D6-N: 6h/day, VL3D12: 3h/day)
- ✅ 物理意義清晰 (光合機構飽和/抑制)

**缺點**:
- 增加 2 個參數 (acute_inhibition_coeff, K_acute)
- 需要重新校準

---

## 下一步待辦

1. **選項 A: 接受 v6.3 結果** (推薦 ⭐⭐⭐)
   - 5/6 組 < 5% (83.3%)
   - L6D6 > CK 確認
   - H12D3 為模型固有限制 (短期高強度問題)

2. **選項 B: 實現即時損傷機制** (推薦 ⭐⭐⭐⭐)
   - 能夠解決 H12D3 問題
   - 預期所有組都能 < 5%
   - 需要增加 2 個參數和重新校準

3. **花青素預測改善** (可選)
   - v6.3 主要關注鮮重預測
   - 花青素預測可能需要重新校準

---

## 歷史修改歷程 (v5.0-v5.9)

### 1. 將時間閾值改為能量單位公式

**需求**: 用戶希望將時間閾值 (6 小時) 改為能量單位 [kJ/m²]，使公式更具通用性。

**原公式 (時間基礎)**:
```python
intraday_factor = 1 + 1.5 × softplus(hours - 6, 3.0)²
```

**新公式 (能量基礎)**:
```python
intraday_factor = 1 + k × softplus((E - E_50) / E_scale, sharpness)^m

其中:
  E = I_UVA_config × hours_elapsed × 3.6 [kJ/m²]
  E_50 = 237.6 kJ/m² (≈6h @ 11 W/m²)
  E_scale = 39.6 kJ/m² (≈1h @ 11 W/m²)
  k = 1.5
  m = 2.0
  sharpness = 3.0
```

**關鍵修正**:
- 在計算 E_elapsed 時，使用 `env['I_UVA']` (設定的 UVA 強度) 而非 `I_UVA` (當前時刻的強度)
- 這確保能量公式與時間公式數學完全等價

**通用性**:
- 參數 E_50 和 E_scale 都是能量單位 [kJ/m²]
- 適用於任何 UVA 強度，只需根據不同實驗條件調整 E_50 和 E_scale
- softplus 確保完全連續，無硬性閾值

### 2. 修改的檔案位置

`simulate_uva_model.py`:
- 參數定義: 第 78-110 行 (新增 E_50, E_scale 參數)
- 計算邏輯: 第 337-351 行 (E_elapsed 和 intraday_factor 計算)

---

## 已完成事項 (v5.0-v6.3)

### v5.x 系列 (花青素機制改進)
1. ✅ 移除指數型天數累積 (v5.4)
2. ✅ 修正夜間 UVA-PAR 光合貢獻 (v5.5)
3. ✅ 將時間閾值改為能量單位公式 (v5.7)
4. ✅ 移除硬閾值，引入 Stress 累積能量 (v5.9)
5. ✅ 花青素機制改用 Hill 函數 (v5.9)

### v6.x 系列 (參數外部化 + PAR 修正)
6. ✅ **參數完全外部化** (v6.0) - 創建 params_config.py ⭐⭐⭐
7. ✅ **移除 PAR 放大效應** (v6.0) - par_conversion_factor: 3.0 → 1.0 ⭐⭐⭐
8. ✅ **修復時間計算 Bug** (v6.1) - day_from_sowing 計算錯誤 ⭐⭐⭐⭐⭐
9. ✅ 提高抑制靈敏度 (v6.2)
10. ✅ **L6D6 優先穩定** (v6.3) - 確保 L6D6 > CK ⭐⭐⭐⭐
11. ✅ **5/6 組 < 5%** (v6.3) - 83.3% 達標率 ⭐⭐⭐

### 其他
12. ✅ 完成敏感度分析 (25個參數 × 4處理組)
13. ✅ 更新論文敏感度分析章節
14. ✅ 更新專案超級總結文件

---

## 擴展敏感度分析結果摘要 (2025-12-17)

**分析範圍**: 25 個自創參數 × 4 處理組
- UVA_Stress: 13 個參數
- Anthocyanin: 7 個參數
- Carbon_Repair: 3 個參數
- LDMC: 2 個參數

### H12D3 處理組 (高劑量，最具代表性)

| 排名 | FW 敏感參數 | S_FW | Stress 敏感參數 | S_Stress | Anth 敏感參數 | S_Anth |
|------|-------------|------|-----------------|----------|---------------|--------|
| 1 | LAI_ref_vuln | -1.28 | LAI_ref_vuln | +13.2 | LAI_ref_vuln | +1.60 |
| 2 | E_50 | +0.61 | E_50 | -3.35 | k_deg | -0.95 |
| 3 | E_scale | +0.39 | E_scale | -2.57 | E_50 | -0.85 |

### L6D6-N 處理組 (夜間照射)

| 排名 | FW 敏感參數 | S_FW | Stress 敏感參數 | S_Stress |
|------|-------------|------|-----------------|----------|
| 1 | LAI_ref_vuln | -0.95 | LAI_ref_vuln | +19.5 |
| 2 | stress_damage_coeff | -0.20 | stress_damage_coeff | +1.56 |
| 3 | circadian_disruption_factor | -0.20 | circadian_disruption_factor | +1.53 |

### 關鍵發現

**UVA 逆境機制**:
- LAI_ref_vuln 是模型最敏感參數 (FW: -0.7~-1.3, Stress: +13~23)
- 能量閾值 E_50, E_scale 在高劑量處理中更重要
- circadian_disruption_factor 只在夜間處理 (L6D6-N) 有效

**花青素機制**:
- k_deg (降解率) 對花青素濃度影響最大 (S ≈ -1.0)
- base_anth_rate_light (日間合成) 次之 (S ≈ +0.8)
- V_max_anth (Stress誘導合成) 只在高劑量處理有輕微影響

**碳依賴修復**:
- base_repair_capacity 與 stress_repair_coeff 敏感度相近
- carbon_repair_bonus 和 K_carbon 影響極小 (碳池通常不限制修復)

**LDMC**:
- ldmc_stress_sensitivity 在高劑量處理有中等影響 (FW: -0.3)

**CK 對照組**:
- 所有 UVA 相關參數敏感度為 0 (符合預期)
- 只有花青素基礎合成率和降解率有影響

---

## 下一步待辦

1. ~~**敏感度分析**: 測試關鍵參數的敏感度~~ ✅ 已完成
2. **考慮移除 model_config.py 中未使用的參數先驗** (PARAM_PRIORS, NIGHT_INHIBITION_SCENARIOS)

### 閾值校準結果 (供參考)

網格搜尋結果顯示最佳閾值為 5.5h (Mean Error 2.0%)，但 6.0h (Mean Error 2.3%) 非常接近。
考慮到 6.0h 是整數且差異僅 0.3%，決定維持 E_50 = 237.6 kJ/m² (≈6h @ 11 W/m²)。

```
Threshold | Mean Error | Max Error
----------|------------|----------
  5.5h    |   2.0%     |   5.0%    <-- 數學最佳
  6.0h    |   2.3%     |   5.0%    <-- 採用值
```

---

## 檔案結構

### 核心模型文件
| 檔案 | 版本 | 說明 |
|------|------|------|
| **`simulate_uva_model.py`** | **v6.3** | **主模型 - 參數外部化版本** |
| **`params_config.py`** | **v6.3** | **參數配置文件 (NEW!)** |
| `lettuce_uva_carbon_complete_model.py` | v7.1 | Sun 基礎模型 |

### 文檔文件
| 檔案 | 說明 |
|------|------|
| `CLAUDE.md` | Claude 工作守則 |
| `MODEL_DESIGN_NOTES.md` | 模型設計筆記 (需更新至 v6.3) |
| `HANDOFF_STATUS.md` | 本文件 |

### 報告文件 (v6.x 系列)
| 檔案 | 版本 | 說明 |
|------|------|------|
| **`V63_FINAL_REPORT.md`** | **v6.3** | **v6.3 最終校準報告** ⭐ |
| `V62_FINAL_REPORT.md` | v6.2 | v6.2 校準報告 |
| `V61_BUG_FIX_REPORT.md` | v6.1 | v6.1 Bug 修復報告 |
| `V61_RECALIBRATION_REPORT.md` | v6.1 | v6.1 重新校準報告 |

### 分析腳本
| 檔案 | 說明 |
|------|------|
| `generate_paper_figures.py` | 論文圖表生成 |
| `sensitivity_analysis_extended.py` | 擴展敏感度分析 (25參數) |

---

**模型校準完成 v6.3！5/6 組 < 5%，L6D6 > CK 確認。FW: Mean 9.7%, Median 1.9%** ✅
