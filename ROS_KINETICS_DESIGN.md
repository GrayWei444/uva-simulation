# ROS 動力學模型設計文檔 v8.0

**日期**: 2025-12-29
**狀態**: 探索完成，結論總結
**目標**: 用 ROS 動力學推導損傷-修復機制，取代 v7.0 的時間冪函數

---

## 執行摘要

經過深入探索，我們發現：
1. **v7.0 的 n=7 無法從 ROS 動力學推導** - 反應動力學最高約 3-4 級
2. **需要 n≈6.6 才能產生 12h/6h = 95 的比值** - 超出物理合理範圍
3. **v7.0 的高冪次是經驗擬合，不是第一原理推導**

建議：保留 v7.0 模型，但在論文中用 ROS 動力學 **解釋** 而非 **推導** 該機制。

---

## 1. 問題分析

### v7.0 的「軟閾值」問題

v7.0 使用純冪函數：
```
day_factor = 1 + k_day × hours^n_day
# k_day = 1e-5, n_day = 7
# 6h: factor = 3.8
# 12h: factor = 359 (128 倍)
```

**問題**：
- n_day = 7 是預設形狀，不是從第一原理推導
- 12^7 / 6^7 = 128 是人為設定，不是物理/生化原理的結果
- 無法解釋「為什麼是 7 次方？」

### 目標：基於 ROS 動力學推導

用酶動力學和能量積分，讓模型形狀從生化原理 **自然浮現**。

---

## 2. ROS 動力學基礎

### 2.1 關鍵酶參數 (文獻值)

| 酶 | 反應 | Km | 催化常數 | 來源 |
|---|---|---|---|---|
| CAT | 2H₂O₂ → O₂ + 2H₂O | 40-600 mM | 一級反應，無飽和 | Mhamdi et al. 2010 |
| SOD | 2O₂⁻ → H₂O₂ + O₂ | - | 2×10⁹ M⁻¹s⁻¹ | Wikipedia/PMC |
| APX | H₂O₂ + ASC → 2H₂O | 33-76 µM | 高親和力 | PMC |

### 2.2 關鍵生化洞察

1. **APX 先飽和**：Km_APX (~50 µM) << Km_CAT (~200 mM)
   - 低 ROS：APX 高效清除
   - 高 ROS：APX 飽和，CAT 接手
   - 極高 ROS：兩者都飽和 → 損傷累積

2. **ROS 產生與 UVA 能量成正比**
   - 光敏劑（葉綠素）吸收 UVA → 產生 ¹O₂
   - 產生率 ∝ I_UVA

3. **修復有上限**
   - 酶活性有限（Vmax）
   - Ascorbate/GSH 再生需要 ATP

---

## 3. 新模型設計：ROS 動態平衡

### 3.1 狀態變量

在原有 5 個狀態變量基礎上，新增：

```
[X_d, C_buf, LAI, Anth, Stress, ROS]
```

或者將 Stress 重新定義為「ROS 累積損傷」。

### 3.2 ROS 動態方程

```python
dROS/dt = production - scavenging_APX - scavenging_CAT

# ROS 產生 (能量積分)
production = k_ros × I_UVA × α_chlorophyll

# APX 清除 (高親和力，低 Km，容易飽和)
scavenging_APX = V_max_APX × ROS / (Km_APX + ROS) × [ASC]

# CAT 清除 (低親和力，高 Km，不易飽和)
scavenging_CAT = k_CAT × ROS  # 近似一級反應
```

### 3.3 損傷累積

```python
dDamage/dt = k_damage × max(0, ROS - ROS_threshold)

# 當 ROS 超過清除能力時，開始累積損傷
```

### 3.4 修復動態

```python
dRepair/dt = V_repair × C_buf / (K_repair + C_buf) - k_deg_repair × Repair

# 修復能力與碳資源相關
# 有降解/損耗
```

### 3.5 Stress 重新定義

```python
Stress = Damage - Repair

# 淨損傷 = 累積損傷 - 已修復
```

---

## 4. 簡化模型：不新增狀態變量

考慮到校準複雜度，可以用 **準穩態假設** 簡化。

### 4.1 核心假設

假設 ROS 處於準穩態 (dROS/dt ≈ 0)：

```python
production ≈ scavenging_total

k_ros × I_UVA = V_max_total × ROS_ss / (Km_eff + ROS_ss)
```

### 4.2 解準穩態 ROS

```python
ROS_ss = Km_eff × k_ros × I_UVA / (V_max_total - k_ros × I_UVA)
```

當 `k_ros × I_UVA` 接近 `V_max_total` 時，ROS_ss → ∞（清除能力飽和）

### 4.3 能量積分

關鍵創新：用 **累積能量** 而非 **照射時間**

```python
# 累積能量（J/m²）
E_accumulated += I_UVA × dt

# 清除能力隨能量累積而下降（酶/抗氧化劑耗竭）
V_eff = V_max × exp(-E_accumulated / E_scale)
```

### 4.4 損傷因子（取代 day_factor）

```python
# 能量依賴的損傷因子
E = I_UVA × hours_today × 3600  # J/m²

# Michaelis-Menten 型飽和
saturation = E / (Km_energy + E)

# 損傷因子
damage_factor = 1 + k_damage_max × saturation
```

**關鍵優勢**：
- `Km_energy` 對應酶的 Km（可從文獻值推導）
- 形狀是 Michaelis-Menten 而非任意冪函數
- 物理意義清晰：能量積分 → 酶飽和 → 損傷加速

---

## 5. 實作方案

### 5.1 方案 A：替換 day_factor（最小改動）

```python
# v7.0 原版
day_factor = 1 + k_day × hours^n_day

# v8.0 新版 (能量積分 + Michaelis-Menten)
E_today = I_UVA × hours_today × 3600  # J/m²
Km_energy = 200000  # J/m² (待校準，對應酶 Km)
saturation = E_today / (Km_energy + E_today)
day_factor = 1 + k_damage_max × saturation

# 或更完整的雙酶模型
sat_APX = E_today / (Km_APX_energy + E_today)  # APX 先飽和
sat_CAT = E_today / (Km_CAT_energy + E_today)  # CAT 後飽和
day_factor = 1 + k_APX × sat_APX + k_CAT × sat_CAT
```

### 5.2 方案 B：新增 ROS 狀態變量（完整動力學）

```python
# 狀態變量: [X_d, C_buf, LAI, Anth, Stress, ROS]

def dROS_dt(ROS, I_UVA, p):
    production = p.k_ros_prod × I_UVA
    scav_APX = p.V_APX × ROS / (p.Km_APX + ROS)
    scav_CAT = p.k_CAT × ROS
    return production - scav_APX - scav_CAT

# Stress 從 ROS 累積推導
def dStress_dt(Stress, ROS, p):
    damage = p.k_damage × max(0, ROS - p.ROS_threshold)
    repair = p.k_repair × Stress
    return damage - repair
```

---

## 6. 參數推導

### 6.1 從文獻值推導 Km_energy

APX 的 Km_H2O2 ≈ 50 µM = 5×10⁻⁵ M

假設：
- ROS 產生率 = k × I_UVA
- 當 ROS 達到 Km 時，APX 處於半飽和

需要估算 k（ROS 產生係數）：
- 葉綠素含量約 0.5 mg/g FW
- UVA 吸收截面 ~10⁻² m²/g chlorophyll
- 量子產率 ~0.01-0.1

### 6.2 Km_energy 估計

對於 I_UVA = 11 W/m²：
- 6 小時能量 = 11 × 6 × 3600 = 237,600 J/m²
- 12 小時能量 = 11 × 12 × 3600 = 475,200 J/m²

如果設 Km_energy = 200,000 J/m²：
- 6h: saturation = 237600 / (200000 + 237600) = 54%
- 12h: saturation = 475200 / (200000 + 475200) = 70%
- 比值 = 70%/54% = 1.3

如果設 Km_energy = 50,000 J/m²（更接近 APX 飽和）：
- 6h: saturation = 237600 / (50000 + 237600) = 83%
- 12h: saturation = 475200 / (50000 + 475200) = 90%
- 比值 = 1.1

**問題**：單一 Michaelis-Menten 函數無法產生 v7.0 的 128 倍差距

### 6.3 解決方案：多級酶飽和

生物現實：APX、CAT、GPX 依序飽和

```python
# 三級飽和模型
# Level 1: APX (最先飽和)
sat1 = E / (Km1 + E)  # Km1 = 50,000 J/m²

# Level 2: CAT (次飽和)
sat2 = E / (Km2 + E)  # Km2 = 200,000 J/m²

# Level 3: GSH/其他 (最後飽和)
sat3 = E / (Km3 + E)  # Km3 = 500,000 J/m²

# 損傷因子 = 所有層級的累積效應
day_factor = 1 + k1×sat1 + k2×sat1×sat2 + k3×sat1×sat2×sat3
```

這樣可以產生更陡峭的響應曲線。

---

## 7. 夜間節律的 ROS 解釋

### 7.1 生理機制

夜間 UVA 造成更大損傷的 ROS 解釋：
1. 夜間抗氧化酶活性下降（節律調控）
2. Ascorbate 再生需要光合作用的 ATP/NADPH
3. 夜間 GSH/GSSG 比例下降

### 7.2 模型實作

```python
# 夜間酶活性下降
V_APX_night = V_APX_day × (1 - k_circadian × hours_in_dark / (K_circ + hours_in_dark))

# 或用節律振盪函數
phase = (hour - 6) / 24 × 2π  # 假設清晨 6 點為高峰
V_APX = V_APX_max × (0.5 + 0.5 × cos(phase))
```

---

## 8. 實作計畫

### Phase 1：驗證概念
1. 實作方案 A（替換 day_factor）
2. 用雙酶 Michaelis-Menten 模型
3. 校準 Km 值使結果匹配 v7.0

### Phase 2：完整動力學
1. 新增 ROS 狀態變量
2. 實作完整的產生-清除動態
3. 驗證是否能捕捉實驗現象

### Phase 3：文獻驗證
1. 比較模型參數與文獻酶動力學值
2. 撰寫論文中的機制說明

---

## 9. 參考文獻

1. Mhamdi A et al. (2010). Catalase function in plants. J Exp Bot 61(15):4197-4220
2. SOD reaction rate: ~2×10⁹ M⁻¹s⁻¹ (Wikipedia/PMC)
3. APX Km values: 33-76 µM (PMC)
4. CAT Km values: 40-600 mM (Chelikani et al. 2004)

---

## 10. 探索結論 (2025-12-29)

### 10.1 測試結果摘要

| 方案 | 機制 | 12h/6h 比值 | 結果 |
|------|------|------------|------|
| v7.0 | 1 + k × hours^7 | 94.6 | ✓ 12/12 達標 |
| v8.0a (閾值) | 1 + k × (E - E_th)^3.29 | 84.1 | ✗ L6D6-N, VL3D12 失敗 |
| v8.0b (軟閾值) | softplus + 冪律 | 61.5 | ✗ VL3D12 Stress=0 |
| v8.0c (無閾值) | 1 + k × E^2.5 | 2.4 | ✗ 全部失敗 |

### 10.2 核心發現

1. **冪次要求分析**
   - v7.0 的 12h/6h = 94.6 比值來自 2^7 = 128
   - 要用能量冪次達到相同效果：2^n = 95 → n ≈ 6.6
   - 反應動力學最高約 3-4 級（二級反應 + 級聯效應）
   - **結論：無法從反應動力學推導 n=7**

2. **閾值模型問題**
   - 閾值設計導致短時間照射 (VL3D12, 3h/day) 完全沒有損傷
   - 軟閾值（softplus）在 E << E_threshold 時效果仍太小
   - **結論：閾值模型不適合此應用**

3. **Michaelis-Menten 限制**
   - 酶動力學會趨近飽和 (factor → constant)
   - 無法產生 v7.0 那種持續陡峭增長的曲線
   - **結論：單純酶動力學不足以解釋實驗現象**

### 10.3 物理意義解釋 (建議用於論文)

雖然無法從第一原理推導 n=7，但可以這樣解釋：

> UVA 照射導致 ROS 累積，當 ROS 產生率超過抗氧化系統清除能力時，
> 未清除的 ROS 引發一系列級聯反應：
> - 蛋白質氧化 → 酶失活 → 清除能力進一步下降
> - 脂質過氧化 → 膜損傷 → 離子洩漏
> - DNA 損傷 → 轉錄受阻
>
> 這種正反饋級聯導致損傷與照射時間呈高度非線性關係，
> 在模型中用冪函數 (hours^n, n≈7) 經驗擬合。
> 高冪次反映了多層級聯效應的累積，而非單一反應動力學。

### 10.4 建議

1. **保留 v7.0 模型**
   - 已達到 12/12 目標
   - 模型結構合理，預測準確

2. **論文中的表述**
   - 用 ROS 動力學 **解釋** 高冪次的生物學意義
   - 不聲稱從第一原理 **推導** 出 n=7
   - 強調這是「經驗擬合反映多層級聯效應」

3. **未來研究方向**
   - 實驗測量 ROS 動態以驗證級聯假說
   - 測量抗氧化酶活性隨照射時間的變化
   - 建立更詳細的多層 ROS 模型（研究用途，非預測用途）

---

## 附錄：ROS 動力學文獻值

| 酶 | 反應 | Km | 催化常數 | 來源 |
|---|---|---|---|---|
| CAT | 2H₂O₂ → O₂ + 2H₂O | 40-600 mM | 一級反應 | Mhamdi et al. 2010 |
| SOD | 2O₂⁻ → H₂O₂ + O₂ | - | 2×10⁹ M⁻¹s⁻¹ | PMC |
| APX | H₂O₂ + ASC → 2H₂O | 33-76 µM | 高親和力 | PMC |

來源:
- [Catalase in plants](https://academic.oup.com/jxb/article/61/15/4197/436494)
- [SOD kinetics](https://en.wikipedia.org/wiki/Superoxide_dismutase)
- [APX review](https://pmc.ncbi.nlm.nih.gov/articles/PMC10333675/)
