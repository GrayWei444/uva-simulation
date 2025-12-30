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

## 4. 已確認的核心機制

| 機制 | 類型 | 說明 |
|------|------|------|
| 日內逆境累積 | 連續冪函數 | `1 + k × (hours/h_ref)^n` |
| 夜間節律干擾 | **待修改** | 需改為連續函數 |
| LAI 脆弱性 | 連續 sigmoid | `cap × (LAI_ref/LAI)^n / (cap + ...)` |
| ROS 清除 | 連續 Michaelis-Menten | `Vmax × S / (Km + S)` |

---

## 5. 關鍵文件結構

| 文件 | 用途 |
|------|------|
| `simulate_uva_model.py` | 主模擬腳本 (v6.9) |
| `simulate_uva_model_ros.py` | ROS 版本 |
| `MODEL_DESIGN_NOTES.md` | 模型設計筆記 |
| `HANDOFF_STATUS.md` | 聊天交接狀態 |
| `CLAUDE.md` | 本文件 |

---

**嚴格遵守以上守則，特別是「禁止硬閾值」原則！**
