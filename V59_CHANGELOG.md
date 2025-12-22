# v5.9 版本變更說明

## 變更日期
2025-12-22

## 變更動機

### 問題
v5.7 的花青素預測誤差過大，特別是：
- VL3D12: 實驗 48.2 ppm, 預測 64.6 ppm (**誤差 +34.1%**)
- L6D12: 實驗 51.8 ppm, 預測 72.5 ppm (**誤差 +40.0%**)

### 根本原因

v5.7 使用硬閾值機制：
```python
# v5.7 舊機制
effective_stress = softplus(Stress - stress_threshold_anth, sharpness)
# stress_threshold_anth = 15.0  <-- 硬編碼閾值
```

**問題**:
1. **硬閾值 15.0 是硬編碼的**，不符合"所有參數可調、無硬編碼"的設計理念
2. **無法正確描述中等 Stress × 長時間的情境**
   - VL3D12: Stress=6.5 (< 15)，但照射 12 天
   - L6D12: Stress=11.2 (< 15)，但照射 12 天
   - 這兩組的瞬時 Stress 雖低，但累積效應應該誘導花青素
3. **不符合生物學**：花青素是對「累積脅迫」的響應，而非瞬時脅迫

## 解決方案

### 核心概念：Stress 累積能量 (E_stress)

引入新的狀態變量 `E_stress`，表示 Stress 對時間的積分：

```
E_stress = ∫ Stress(t) dt  [單位: Stress·day]
```

**生物學意義**: 植物對 UV 脅迫的「記憶」，累積的脅迫會持續誘導花青素合成。

### 新機制：Hill 函數

```python
# v5.9 新機制
E_power_n = E_stress ** n_hill_anth
K_power_n = K_E_anth ** n_hill_anth
uva_induced = V_max_anth * E_power_n / (K_power_n + E_power_n)
```

**特性**:
- ✅ **完全連續**：無硬閾值，符合微分方程連續性要求
- ✅ **可參數化**：所有參數都可調整
- ✅ **符合生物學**：累積脅迫驅動花青素合成
- ✅ **有飽和效應**：高 E_stress 時不會無限增長
- ✅ **通用性**：適用於其他實驗條件和作物

## 詳細變更

### 1. 新增狀態變量

```python
# 狀態變量從 5 個增加為 6 個
[X_d, C_buf, LAI, Anth, Stress, E_stress]  # v5.9
[X_d, C_buf, LAI, Anth, Stress]            # v5.7

# E_stress 微分方程
dE_stress_dt = Stress / 86400.0  # 轉換為 [Stress·day/s]
```

### 2. 修改花青素參數

**移除的參數** (v5.7):
```python
stress_threshold_anth = 15.0         # 硬閾值 (已移除)
anth_threshold_sharpness = 0.5       # 閾值平滑度 (已移除)
K_m_anth = 30.0                      # Michaelis-Menten 常數 (已移除)
```

**新增的參數** (v5.9):
```python
K_E_anth = 50.0          # 半飽和 Stress 累積能量 [Stress·day]
n_hill_anth = 2.0        # Hill 係數 (控制響應曲線陡度)
```

**保留的參數**:
```python
base_anth_rate_light = 4.11e-10   # 日間基礎合成率 [kg/m²/s]
base_anth_rate_dark = 2.05e-10    # 夜間基礎合成率
V_max_anth = 3.0e-10              # 最大誘導合成率 (值已調整)
k_deg = 2.5e-6                    # 降解率
```

### 3. 修改花青素計算邏輯

**v5.7 舊邏輯**:
```python
effective_stress = softplus(Stress - p.stress_threshold_anth, p.anth_threshold_sharpness)
uva_induced = p.V_max_anth * effective_stress / (p.K_m_anth + effective_stress)
```

**v5.9 新邏輯**:
```python
E_power_n = E_stress ** p.n_hill_anth
K_power_n = p.K_E_anth ** p.n_hill_anth
uva_induced = p.V_max_anth * E_power_n / (K_power_n + E_power_n + 1e-12)
```

## Hill 函數特性

### 數學形式
```
f(E) = E^n / (K^n + E^n)
```

### 關鍵點
- E = 0: f(0) = 0 (無誘導)
- E = K: f(K) = 0.5 (半飽和)
- E >> K: f(E) → 1 (飽和)

### n (Hill 係數) 的影響
- n = 1: Michaelis-Menten 型 (較平滑)
- n = 2: 中等陡度 (協同效應)
- n > 2: 更陡峭 (強協同效應)

### 參數物理意義

| 參數 | 意義 | 預期值範圍 |
|------|------|-----------|
| V_max_anth | 最大誘導合成率 | 1e-10 ~ 1e-8 kg/m²/s |
| K_E_anth | 達到一半最大響應所需的累積能量 | 20 ~ 100 Stress·day |
| n_hill_anth | 響應曲線陡度 | 1.5 ~ 4.0 |

## 預期效果

### 對各處理組的影響

| Treatment | Stress (瞬時) | 照射天數 | E_stress (估) | 預期改善 |
|-----------|--------------|---------|--------------|---------|
| CK        | 0            | 0       | 0            | 無影響 (基礎合成) |
| L6D6      | 2.0          | 6       | ~6           | 輕微影響 |
| L6D6-N    | 4.8          | 6       | ~14          | 輕微影響 |
| H12D3     | 36.1         | 3       | ~54          | 保持高誘導 |
| VL3D12    | 6.5          | 12      | ~39          | **降低預測值** ✅ |
| L6D12     | 11.2         | 12      | ~67          | **降低預測值** ✅ |

**關鍵改善**:
- VL3D12 和 L6D12 的 E_stress 雖然較高，但使用 Hill 函數後，響應會**飽和**
- 這符合實驗觀察：L6D12 (E_stress ≈ 67) 的花青素只比 VL3D12 (E_stress ≈ 39) 高一點點

## 校準計畫

### 階段 1: 手動測試 (test_anthocyanin_v59.py)
1. 使用建議參數測試
2. 觀察 E_stress 的實際值
3. 手動調整參數

### 階段 2: 自動優化 (calibrate_anthocyanin_v59.py)
1. 使用差分進化法全域優化
2. 給予 VL3D12 和 L6D12 更高權重
3. 同時考慮 FW 和 Anth 誤差

## 相容性

### 向後相容性
- ❌ **不相容**: 狀態變量從 5 個增加為 6 個
- 需要更新所有使用 `simulate_uva_model.py` 的腳本

### 需要更新的檔案
- [x] `simulate_uva_model.py` - 主模型
- [ ] `generate_paper_figures.py` - 圖表生成
- [ ] `run_scenarios.py` - 場景測試
- [ ] `optimize_parameters.py` - 參數優化
- [ ] 其他使用模型的腳本

## 驗證標準

**成功標準**:
- FW: Mean |Err| < 3%, Max < 5% (與 v5.7 相當)
- Anth: Mean |Err| < 15%, Max < 25% (改善 VL3D12, L6D12)

**特別關注**:
- VL3D12: 誤差從 +34% 降至 < +20%
- L6D12: 誤差從 +40% 降至 < +25%

## 參考文獻

### 花青素與累積脅迫
- Gould, K.S. (2004). *Anthocyanin as stress protectant.* J. Theor. Biol.
- Chalker-Scott, L. (1999). *Environmental significance of anthocyanins in plant stress responses.* Photochem. Photobiol.

### Hill 函數在生物學中的應用
- Hill, A.V. (1910). *The possible effects of the aggregation of the molecules of haemoglobin.* J. Physiol.
- Weiss, J.N. (1997). *The Hill equation revisited.* FASEB J.

## 附註

此版本是對 v5.7 的重要改進，移除了硬編碼閾值，使模型更具通用性和生物學合理性。

---

**版本**: v5.9
**作者**: Gray (with Claude)
**日期**: 2025-12-22
