# v6.0 Debug 完整總結

**日期**: 2025-12-22
**狀態**: ✅ 主要 Bug 已修復，參數需要重新校準

---

## 🎯 任務回顧

### 用戶的原始要求

1. **移除 `par_conversion_factor` 放大效應**
   - 原值: 3.0 → 新值: 1.0
   - 理由: LED 57 W/m² + UVA 22 W/m² = 79 W/m²，Sun 原始模型已經能解釋鮮重增加

2. **增加 `c_alpha` 匹配 Sun 模型**
   - 原值: 0.54 → 新值: 0.68 (後校準為 0.57)

3. **加入花青素碳成本機制** (Gemini 建議)
   - 花青素合成會消耗碳

4. **參數重構**
   - 將所有參數移到統一配置文件

---

## 🔍 Debug 過程

### 階段 1: 初步測試

**發現問題**:
- CK 預測 102.7g (實驗 87.0g) → 誤差 +18.1%
- L6D6 預測 55.2g (實驗 91.4g) → 誤差 -39.5%

**初步分析**:
- `c_alpha = 0.68` 太高，需要校準

### 階段 2: 校準 c_alpha

**調整**: `c_alpha = 0.68 → 0.57`

**結果**:
- CK 用 Sun 原始模型: 86.7g ✅
- L6D6 用 Sun 原始模型: 90.6g ✅
- **Sun 模型完美解釋 UVA 效應！**

**BUT**:
- UVA 模型的 CK: 仍然很高
- UVA 模型的 L6D6: 仍然很低

### 階段 3: 追蹤 dC_buf/dt 異常

**發現**:
- Sun 模型: dC_buf/dt = +6.69e-10 (正值，碳增加)
- UVA 模型: dC_buf/dt = -7.21e-08 (負值，碳減少)
- 差異: **-10,880%** ❌

**懷疑**:
- 花青素碳消耗過大？
- 修復碳消耗過大？
- Sun 模型調用有問題？

**測試**:
- 花青素碳成本設為 0 → 問題仍存在
- Stress ≈ 0 → 修復消耗極小
- Sun 模型調用測試 → 正確

**結論**: 問題更深層

### 階段 4: 深度追蹤 (Breakthrough!)

**創建**: `debug_cbuf_trace.py` 逐步追蹤 dC_buf/dt

**發現**:
- 當 Stress = 0 時，UVA 模型與 Sun 模型**完全一致** ✅
- 這表示問題不在 UVA 模型的計算邏輯
- 問題在於: **為什麼實際模擬中 Stress ≠ 0？**

### 階段 5: 檢查實際 Stress 值

**創建**: `check_stress_values.py`

**震驚發現**:
- Day 28 (UVA 前): Stress = **1.90** ❌❌❌
- Day 30 (UVA 中): Stress = **0.27**
- Day 35 (結束): Stress = 0.03
- **C_buf 變成負值**: Day 31: C_buf = -0.00000004

**問題**: 為什麼在 UVA 開始前 (Day 28) 就有巨大的 Stress？

### 階段 6: 檢查 CK 組

**創建**: `check_ck_stress.py`

**發現**:
- CK 組 Stress = 0 ✅ (正確，因為沒有 UVA)
- **BUT**: CK 組預測 FW = **3120g** ❌❌❌
  - 實驗值: 87g
  - 誤差: **+3500%**

**震驚**: 即使 Stress = 0，CK 預測也完全錯誤！

### 階段 7: 對比單點計算

**創建**: `compare_derivatives_direct.py`

**發現**:
- 在單個時間點 (t, y) 測試
- Sun 模型和 UVA 模型輸出**完全一致** ✅
- 這表示 `uva_sun_derivatives` 函數本身沒問題

**結論**: 問題在**初始條件**或**時間範圍**！

### 階段 8: 找到 Bug! 🎉

**檢查 `simulate_uva_model.py` 主程式**:

**Bug 1**: Line 510
```python
initial_state = [Xd_init, 0.0, LAI_init, ...]
#                        ^^^
# C_buf 初始化為 0，而不是 X_d * 0.1！
```

**Bug 2**: Line 520
```python
t_span=(0, simulation_days * 86400)
#      ^
# 從 Day 0 開始，而不是 Day 14 (移植日)！
```

**影響**:
- C_buf = 0 導致初期碳飢餓，後期爆炸式生長
- t = 0 導致時間軸錯亂

### 階段 9: 修復並驗證

**修復**:
1. C_buf_init = X_d_init * 0.1
2. t_start = transplant_day * 86400

**結果**:
- CK: 3120g → **90.5g** ✅ (誤差 +4.1%)
- 巨大進步！

**BUT**:
- UVA 處理組仍被低估 18-34%
- 原因: Stress 累積過多

---

## 📊 當前結果

### 修復後的預測 (v6.0 Bug 已修復)

| 處理組 | 預測 FW | 目標 FW | 誤差 | 預測 Anth | 目標 Anth | 誤差 |
|--------|---------|---------|------|-----------|-----------|------|
| CK | 90.5g | 87.0g | +4.1% ✅ | 39.3 ppm | 43.3 ppm | -9.3% |
| L6D6 | 60.0g | 91.4g | -34.3% ⚠️ | 63.5 ppm | 49.4 ppm | +28.5% |
| L6D6-N | 57.1g | 80.8g | -29.3% ⚠️ | 67.6 ppm | 49.3 ppm | +37.3% |
| H12D3 | 48.0g | 60.6g | -20.9% ⚠️ | 96.7 ppm | 65.1 ppm | +48.7% |
| VL3D12 | 58.7g | 67.0g | -12.4% ⚠️ | 64.2 ppm | 48.2 ppm | +33.2% |
| L6D12 | 55.9g | 60.4g | -7.4% ✓ | 68.3 ppm | 51.8 ppm | +31.9% |

**統計**:
- 鮮重誤差: 平均 **18.1%**, 最大 34.3%
- 花青素誤差: 平均 **31.5%**, 最大 48.7%

---

## 🔬 根本原因分析

### 為什麼 UVA 組被低估？

**關鍵觀察**:
1. CK (無 UVA): Stress = 0, 預測正確 (+4.1%)
2. L6D6 (22 W/m² UVA): Stress ≈ 0.27, 預測偏低 (-34.3%)
3. 實驗中 L6D6 反而促進生長 (+5.1%)

**結論**: **Stress 損傷機制過於敏感**

### 物理解釋

UVA 強度 22 W/m² 是相對溫和的:
- 實驗顯示促進生長 (L6D6: 91.4g > CK: 87.0g)
- 模型預測抑制生長 (L6D6: 60.0g < CK: 90.5g)

**問題參數**:
- `stress_damage_coeff = 3.5e-6` → 太高
- `stress_photosynthesis_inhibition = 0.70` → 太強
- `stress_lai_inhibition = 0.70` → 太強

---

## ✅ 已完成的工作

1. ✅ **參數重構**: 創建 `params_config.py`
2. ✅ **確認 par_conversion_factor=1.0 正確**: Sun 模型測試證明
3. ✅ **找到 dC_buf/dt 異常根本原因**: C_buf 初始化 bug
4. ✅ **修復兩個致命 Bug**:
   - C_buf 初始化錯誤
   - 模擬起始時間錯誤
5. ✅ **CK 預測精準** (+4.1%)

---

## 📋 待完成的工作

### 優先級 1: 調整 Stress 參數 ⚠️

**問題**: Stress 損傷過於敏感，導致 UVA 組低估 18-34%

**建議**:
```python
# 當前值
stress_damage_coeff = 3.5e-6
stress_photosynthesis_inhibition = 0.70
stress_lai_inhibition = 0.70

# 建議值 (減少 50%)
stress_damage_coeff = 1.75e-6
stress_photosynthesis_inhibition = 0.35
stress_lai_inhibition = 0.35
```

**目標**: L6D6 能夠預測促進效應 (≈ 91.4g)

### 優先級 2: 校準花青素參數

**問題**: 基礎合成率過高 22倍

**建議**:
```python
# 當前值
base_anth_rate_light = 4.0e-10

# 建議值
base_anth_rate_light = 1.8e-11  # 原值的 4.5%
```

### 優先級 3: 全面驗證

**目標**:
- 所有處理組鮮重誤差 < 10%
- 所有處理組花青素誤差 < 20%
- Stress 動態生理合理

---

## 💡 重要結論

### ✅ v6.0 的核心理念是正確的

1. **par_conversion_factor = 1.0** 是合理的
   - Sun 模型測試證明: 額外 22 W/m² 可以解釋鮮重增加 4.6%
   - 實驗觀察: 鮮重增加 5.1%
   - 比例: 0.90x (非常接近)

2. **c_alpha = 0.57** 校準正確
   - CK 預測: 90.5g (誤差 +4.1%)
   - 證明基礎模型運作正常

3. **Sun 原始模型可以解釋 UVA 效應**
   - 不需要人為放大 PAR 轉換
   - UVA 提供的額外輻射本身就足夠

### ⚠️ 需要調整的是 Stress 機制參數

**問題不在概念，而在參數校準**:
- 當前 Stress 參數針對高劑量 UVA 校準
- 但對低劑量 (22 W/m²) 過於敏感
- 需要重新平衡損傷和修復的速率

---

## 📁 相關文件

### Debug 腳本
- [compare_cbuf_dynamics.py](compare_cbuf_dynamics.py) - 對比碳緩衝池動態
- [debug_cbuf_trace.py](debug_cbuf_trace.py) - 深度追蹤 dC_buf/dt
- [check_stress_values.py](check_stress_values.py) - 檢查 Stress 累積
- [check_ck_stress.py](check_ck_stress.py) - 檢查 CK 組 Stress
- [compare_derivatives_direct.py](compare_derivatives_direct.py) - 直接對比導數
- [test_all_treatments_v60fixed.py](test_all_treatments_v60fixed.py) - 全處理組測試

### 驗證腳本
- [test_sun_model_uva_effect.py](test_sun_model_uva_effect.py) - Sun 模型 UVA 效應測試

### 文檔
- [V60_DEBUG_LOG.md](V60_DEBUG_LOG.md) - Debug 日誌
- [V60_BUG_FIX_REPORT.md](V60_BUG_FIX_REPORT.md) - Bug 修復報告
- [params_config.py](params_config.py) - 統一參數配置

---

## 🎓 學到的教訓

1. **初始條件極為重要**
   - C_buf = 0 vs C_buf = X_d * 0.1 造成巨大差異
   - 必須與基礎模型保持一致

2. **時間軸必須正確**
   - 模擬起始時間必須匹配實驗設計
   - t = 0 vs t = 14 天會影響所有時間相關的判斷

3. **單點測試 vs 完整模擬**
   - 函數邏輯可能正確，但初始條件錯誤
   - 兩種測試都需要

4. **Debug 策略**
   - 從簡單到複雜
   - 對比 Sun 原始模型
   - 逐步縮小範圍
   - 最後找到根本原因

---

**v6.0 已經非常接近成功，只需要重新校準 Stress 相關參數即可。**
