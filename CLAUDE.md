# Claude Code 工作守則 v2.0

**核心原則: 延續性優先，避免重複工作**

---

## 0. 最重要原則：禁止硬閾值

**所有機制必須使用連續函數，禁止硬閾值 (hard threshold)**

### 禁止的寫法：
```python
# 錯誤！硬閾值
if hours > 6:
    damage = high_damage
else:
    damage = low_damage

# 錯誤！硬閾值
circ = 3.8 if is_night else 1.0

# 錯誤！預設特定時間為閾值
E_threshold = 6 * 3600 * 22 / 1000  # "6小時是閾值"
```

### 正確的寫法：
```python
# 正確！連續冪函數
damage_factor = 1 + k * (hours ** n)

# 正確！連續 sigmoid
circ = 1 + (circ_max - 1) * sigmoid((hour - transition) / scale)

# 正確！連續 softplus，參數由校準決定
factor = 1 + k * softplus((E - E_ref) / E_scale) ** n
# E_ref 是校準結果，不是預設閾值
```

### 設計原則：
1. **所有參數皆可校準**：沒有「預設的最佳時間」
2. **連續可微分**：函數在所有點都平滑
3. **預測性**：模型可以預測任意照射時間的結果，不只是實驗中的 3h、6h、12h
4. **校準結果 vs 預設值**：E_ref=475.2 是「校準出來恰好對應 6 小時」，不是「預設 6 小時為閾值」

---

## 1. 開始新聊天時

**必須先讀取以下文件:**
1. `CLAUDE.md` (本文件) - 了解工作守則
2. `MODEL_DESIGN_NOTES.md` - 了解模型當前狀態和已校準參數
3. `HANDOFF_STATUS.md` - 了解上一個聊天的進度和待辦事項

---

## 2. 工作守則

1. **禁止硬閾值**: 所有機制必須是連續函數
2. **不要重複已完成的工作**: 先讀取文件確認已校準的參數
3. **繼承而非重來**: 使用已確認的 ODE 設定 (RK45 + max_step=60)
4. **記錄決策理由**: 任何參數修改都要記錄在 `MODEL_DESIGN_NOTES.md`
5. **驗證先於整合**: 修改模型後必須運行驗證

---

## 3. Token 快用完或結束聊天時

**必須執行自動交接:**
更新 `HANDOFF_STATUS.md` 文件

---

## 4. 參數校準策略

**校準順序（必須遵守）：**

| 順序 | 處理組 | 主要調整機制 | 目標 |
|------|--------|-------------|------|
| 1 | **L6D6** | 基準組 | Stress ≈ 0，幾乎無損傷 |
| 2 | **L6D6-N** | 夜間節律損傷 | 比 L6D6 損傷大 |
| 3 | **VL3D12** | LAI 脆弱性 | 早期照射損傷 |
| 4 | **L6D12** | LAI 脆弱性 | 早期照射損傷 |
| 5 | **H12D3** | ROS 非線性放大 | 長時間照射損傷 |

**關鍵原則：**
- L6D6 是基準，6h/day 日間照射應該 **幾乎無損傷**
- ROS 非線性放大 **只針對 H12D3**（12h/day）
- D12 組損傷 **靠 LAI 脆弱性**，不是 ROS 放大
- 先讓 L6D6 達標，再依序調整其他組

---

## 5. 已確認的核心機制 (v8.0)

| 機制 | 公式 | 說明 |
|------|------|------|
| ROS 動態 | `dROS/dt = k_prod × I_UVA × amp - k_clear × ROS` | 當日氧化壓力 |
| ROS 非線性放大 | `amp = 1 + k × hours^n` | 針對 H12D3 |
| LAI 脆弱性 | `vuln = A × exp(-k × LAI) + 1` | 針對 D12 組 |
| 夜間節律損傷 | `circ = k × I_UVA × hours_dark^n` | 針對 L6D6-N |
| 損傷公式 | `damage = k × ROS × vuln × (1 - protection)` | 簡化版 |

---

## 6. 關鍵文件結構

| 文件 | 用途 |
|------|------|
| `simulate_uva_model_v8.py` | v8.0 ROS 動態模型版 |
| `simulate_uva_model_v72.py` | v7.x 舊版（備份） |
| `MODEL_DESIGN_NOTES.md` | 模型設計筆記 |
| `HANDOFF_STATUS.md` | 聊天交接狀態 |
| `CLAUDE.md` | 本文件 |

---

**嚴格遵守以上守則，特別是：**
1. **禁止硬閾值**
2. **L6D6 優先達標，Stress ≈ 0**
3. **各組別用各自機制調整**
