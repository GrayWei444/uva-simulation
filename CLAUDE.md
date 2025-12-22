# Claude Code 工作守則 v1.0
_NOTES.md` - 了解模型當前狀態和已校準參數
**核心原則: 延續性優先，避免重複工作**

## 1. 開始新聊天時

**必須先讀取以下文件:**
1. `CLAUDE.md` (本文件) - 了解工作守則
2. `MODEL_DESIGN
3. `HANDOFF_STATUS.md` - 了解上一個聊天的進度和待辦事項

## 2. 工作守則

1. **不要重複已完成的工作**: 在進行任何模擬或參數調整之前，先讀取 `MODEL_DESIGN_NOTES.md` 確認已校準的參數。
2. **繼承而非重來**: 使用已確認的ODE設定 (RK45 + max_step=60)，不要再測試其他求解器。
3. **記錄決策理由**: 任何參數修改都要記錄在 `MODEL_DESIGN_NOTES.md` 中，包含為什麼這樣改。
4. **驗證先於整合**: 修改模型後必須運行驗證，確保CK和L6D6的結果仍在誤差範圍內。

## 3. Token 快用完或結束聊天時

**必須執行自動交接:**
更新 `HANDOFF_STATUS.md` 文件，包含以下內容：

```markdown
# 交接狀態 (HANDOFF_STATUS)
最後更新: [日期時間]

## 當前進度
- [描述已完成的工作]

## 下一步待辦
- [描述需要繼續的工作]

## 重要發現/決策
- [記錄重要的發現或決策]

## 當前模型狀態
- CK: FW=?g, Anth=?ppm
- L6D6: FW=?g, Anth=?ppm
```

## 4. 已確認的核心機制 (不要改變!)

| 機制 | 參數 | 值 |
|------|------|-----|
| ODE求解器 | method | RK45 |
| ODE求解器 | max_step | 60秒 |
| 光合效率 | c_alpha | 0.54 |
| UVA強度 | I_UVA | 11 W/m² |
| PAR轉換係數 | par_conversion_factor | 3.0 |
| LDMC係數 | dw_fw_ratio_reduction_per_dose | 1.0e-9 |

## 5. 關鍵文件結構

| 文件 | 用途 |
|------|------|
| `simulate_uva_model.py` | 主模擬腳本 (已校準) |
| `MODEL_DESIGN_NOTES.md` | 模型設計筆記和已校準參數 |
| `HANDOFF_STATUS.md` | 聊天交接狀態 |
| `CLAUDE.md` | 本文件 - Claude工作守則 |

---
**嚴格遵守以上守則，確保工作延續性。**
