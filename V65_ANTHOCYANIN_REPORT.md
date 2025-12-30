# v6.5 花青素優化報告

**日期**: 2025-12-22  
**版本**: v6.5  
**主要變更**: 花青素改用 E_stress (累積能量) 驅動 + 註解修正

---

## 核心變更

### 1. 花青素機制修正

**問題**: v6.4 使用瞬時 Stress，但註解和設計文檔說應該用 E_stress

**修正** (simulate_uva_model.py:386-391):
```python
# v6.4 (錯誤)
Stress_power_n = Stress ** p.n_stress_anth
uva_induced = p.V_max_anth * Stress_power_n / (K_power_n + Stress_power_n)

# v6.5 (正確)
E_stress_power_n = E_stress ** p.n_stress_anth
uva_induced = p.V_max_anth * E_stress_power_n / (K_power_n + E_stress_power_n)
```

**理由**:
- 花青素基因表達需要累積信號
- E_stress = ∫ Stress dt 代表總暴露劑量
- 符合 v5.9 原始設計意圖

### 2. 參數調整

**花青素參數** (params_config.py):
```python
# v6.4 → v6.5
'V_max_anth': 1.8e-10 → 2.0e-10,        # 最大誘導合成率
'K_stress_anth': 10.0 → 1.0,            # E_stress 半飽和常數 [Stress·day]
'k_deg': 2.6e-6 → 3.0e-6,               # 降解率 (提高以控制長期組累積)
```

**變更原因**:
- `K_stress_anth` 從 10.0 (Stress 單位) → 1.0 (Stress·day 單位)
- 長期組 (VL3D12, L6D12) 累積過多 E_stress，需提高降解率

### 3. 註解修正

修正了容易誤會或錯誤的註解：

1. **intraday_factor 註解** (simulate_uva_model.py:294-304):
   - 明確說明 E_50 = 475.2 kJ/m² (6h @ 22 W/m²)
   - 解釋 H12D3 vs L6D6 的差異

2. **E_stress 註解** (simulate_uva_model.py:411-416):
   - 明確物理意義：總暴露劑量
   - 說明用途：驅動花青素合成

3. **模型版本號更新**:
   - v6.4 → v6.5
   - 所有相關輸出和類別註解

---

## v6.5 預測結果

### 鮮重 (FW) - 保持 100% 達標 ✅

```
Treatment    FW_sim   FW_exp   FW_Err   狀態
----------------------------------------------
CK            88.6g    87.0g   +1.9%    ✅ < 5%
L6D6          89.5g    91.4g   -2.0%    ✅ < 5%
L6D6-N        79.7g    80.8g   -1.4%    ✅ < 5%
VL3D12        67.4g    67.0g   +0.6%    ✅ < 5%
L6D12         58.7g    60.4g   -2.9%    ✅ < 5%
H12D3         62.0g    60.6g   +2.3%    ✅ < 5%
```

### 花青素 (Anth) - 部分改善

```
Treatment   Anth_sim  Anth_exp  Anth_Err  v6.4      改善
-----------------------------------------------------------
CK           35.9ppm    43.3ppm  -17.0%    N/A       N/A
L6D6         42.8ppm    49.4ppm  -13.3%   -17.4%    +4.1% ✅
L6D6-N       55.2ppm    49.3ppm  +12.0%    N/A       N/A
VL3D12       72.6ppm    48.2ppm  +50.6%    N/A       N/A
L6D12        83.9ppm    51.8ppm  +62.0%    N/A       N/A
H12D3        63.8ppm    65.1ppm   -2.0%    N/A      ✅
```

**統計**:
- **H12D3**: -2.0% ✅ (優秀!)
- **L6D6**: -13.3% (改善 4.1%)
- **問題組**: VL3D12, L6D12 過高 (長期累積問題)

---

## 問題分析

### 長期組 (VL3D12, L6D12) 花青素過高

**原因**:
1. **E_stress 持續累積**: 12 天照射導致 E_stress ≈ 5-6 Stress·day
2. **飽和不足**: K_stress_anth = 1.0，長期組達到 ~85% 飽和
3. **降解不足**: 即使提高 k_deg，仍無法抵消累積

**物理解釋**:
- E_stress 是積分量，只會增加不會減少
- 長期低劑量暴露會讓 E_stress 累積很高
- 與實驗觀察不符：長期組花青素應該較低

### 可能的解決方案

#### 選項 A: 使用 E_stress 但加入衰減
```python
# 修改 E_stress 微分方程
dE_stress_dt = Stress / 86400.0 - decay_rate * E_stress
```
讓 E_stress 有「記憶衰退」，避免無限累積。

#### 選項 B: 混合機制
```python
# 短期響應 + 長期累積
instant_term = Stress / (K_instant + Stress)
cumulative_term = E_stress / (K_cumulative + E_stress)
uva_induced = V_max * (w1 * instant_term + w2 * cumulative_term)
```

#### 選項 C: 回到瞬時 Stress
重新校準 V_max_anth 和 K_stress_anth。

---

## 版本對比

| 指標 | v6.4 (Stress) | v6.5 (E_stress) | 變化 |
|------|--------------|----------------|------|
| **FW 達標** | 6/6 (100%) | 6/6 (100%) | 持平 ✅ |
| **L6D6 Anth** | -17.4% | -13.3% | +4.1% ✅ |
| **H12D3 Anth** | N/A | -2.0% | ✅ 優秀 |
| **VL3D12 Anth** | N/A | +50.6% | ❌ 問題 |
| **L6D12 Anth** | N/A | +62.0% | ❌ 問題 |

---

## 修改文件清單

1. **simulate_uva_model.py**:
   - Line 4: 版本號 v6.4 → v6.5
   - Line 31-33: 新增 v6.5 修改紀錄
   - Line 386-391: 花青素改用 E_stress
   - Line 294-304: 修正 intraday_factor 註解
   - Line 411-416: 修正 E_stress 註解
   - Line 441-448: 更新輸出說明

2. **params_config.py**:
   - Line 77-80: 花青素參數調整
   - Line 83: 降解率提高

---

## 下一步建議

### 短期 (推薦)
1. ✅ **接受 v6.5 鮮重結果** (6/6 組 < 5%)
2. ⚠️  **花青素需要進一步工作**:
   - 考慮選項 A (E_stress 衰減)
   - 或選項 B (混合機制)

### 長期
1. 系統性分析花青素機制
2. 收集更多實驗數據
3. 考慮其他環境因子 (溫度、光照等)

---

**版本**: v6.5  
**日期**: 2025-12-22  
**狀態**: 
- ✅ FW: 生產就緒 (6/6 組 < 5%)
- ⚠️  Anth: 需要進一步優化

**主要貢獻**: 修正花青素機制使用 E_stress，符合原始設計意圖。
