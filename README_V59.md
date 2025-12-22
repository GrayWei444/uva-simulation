# 萵苣 UVA 模型 v5.9 - 快速開始

## 🚀 立即開始

### 選項 1: 測試當前參數 (推薦先執行)

```bash
cd /home/coder/projects/uva-simulation
python3 test_anthocyanin_v59.py
```

這個腳本會：
- 顯示所有 6 個處理組的預測結果
- 顯示實際的 E_stress 累積值
- 測試參數敏感度

### 選項 2: 自動校準參數

```bash
cd /home/coder/projects/uva-simulation
python3 calibrate_anthocyanin_v59.py
```

這個腳本會：
- 使用差分進化法自動尋找最佳參數
- 給予 VL3D12 和 L6D12 更高權重
- 顯示校準後的參數和結果

---

## 📁 重要檔案

| 檔案 | 用途 |
|------|------|
| `simulate_uva_model.py` | ⭐ 主模型 (v5.9) |
| `V59_SUMMARY.md` | ⭐ 快速總結 |
| `V59_USAGE_GUIDE.md` | ⭐ 使用指南 |
| `V59_CHANGELOG.md` | 詳細變更說明 |
| `E_stress_analysis.md` | 理論分析 |
| `HANDOFF_STATUS.md` | 當前進度 |
| `test_anthocyanin_v59.py` | 測試腳本 |
| `calibrate_anthocyanin_v59.py` | 校準腳本 |

---

## 🎯 v5.9 核心改進

### 問題
v5.7 的花青素預測誤差過大：
- **VL3D12**: 實驗 48.2 ppm, 預測 64.6 ppm (**誤差 +34%**)
- **L6D12**: 實驗 51.8 ppm, 預測 72.5 ppm (**誤差 +40%**)

### 解決方案
移除硬閾值，使用 **Stress 累積能量** (E_stress) 驅動花青素合成：

```python
# ❌ v5.7 舊機制 (硬閾值)
effective_stress = softplus(Stress - 15.0, sharpness)  # 硬編碼!

# ✅ v5.9 新機制 (能量累積)
E_stress = ∫ Stress(t) dt  # [Stress·day]
uva_induced = V_max × E^n / (K_E^n + E^n)  # Hill 函數
```

### 優勢
✅ 完全連續，無硬編碼閾值
✅ 符合生物學：植物對累積脅迫的響應
✅ 有飽和效應：避免長期低劑量過度響應
✅ 通用性高：適用於其他實驗條件

---

## 🔧 快速調整參數

### 位置
`simulate_uva_model.py` 第 166-171 行

### 參數說明

```python
# 1. 基礎合成率 (影響所有組，包括CK)
self.base_anth_rate_light = 4.11e-10   # 日間
self.base_anth_rate_dark = 2.05e-10    # 夜間

# 2. UVA 誘導參數 (只影響UVA組)
self.V_max_anth = 3.0e-10    # ⬆️增加→所有UVA組花青素↑
self.K_E_anth = 50.0         # ⬆️增加→需要更高E_stress才會誘導
self.n_hill_anth = 2.0       # ⬆️增加→響應曲線更陡

# 3. 降解率 (影響所有組)
self.k_deg = 2.5e-6          # ⬆️增加→所有組花青素↓
```

### 調整策略

| 問題 | 解決方案 |
|------|---------|
| 所有組花青素都偏高 | 增加 `k_deg` 或減少 `V_max_anth` |
| 所有組花青素都偏低 | 減少 `k_deg` 或增加 `V_max_anth` |
| VL3D12/L6D12 偏高 | 增加 `K_E_anth` (降低低能量響應) |
| H12D3 偏低 | 減少 `K_E_anth` 或增加 `n_hill_anth` |

---

## 📊 驗證標準

| 指標 | v5.7 基準 | v5.9 目標 |
|------|----------|----------|
| FW Mean | 2.3% | < 3% ✅ |
| FW Max | 5.0% | < 5% ✅ |
| **VL3D12 Anth** | **+34%** | **< +20%** ⭐ |
| **L6D12 Anth** | **+40%** | **< +25%** ⭐ |
| Anth Mean | 15.6% | < 12% |
| Anth Max | 40.0% | < 25% |

---

## 🆘 常見問題

### Q: 模型跑不動？
```bash
# 檢查 Python 套件
python3 -c "import numpy, scipy; print('OK')"

# 如果失敗，重新安裝
sudo apt-get install python3-numpy python3-scipy python3-matplotlib python3-pandas
```

### Q: E_stress 是什麼？
E_stress (Stress 累積能量) 是 Stress 對時間的積分，單位 [Stress·day]。代表植物對 UV 脅迫的「記憶」。

### Q: 為什麼不直接用 Stress？
Stress 會修復（有正有負），但植物對脅迫的記憶不會消失。花青素是對「累積脅迫歷史」的響應。

### Q: Hill 係數 n 的意義？
n > 1 表示協同效應：累積脅迫需達一定程度後，花青素才會快速合成（類似酶的協同調控）。

---

## 📞 技術支援

遇到問題？請檢查：
1. ✅ Python 環境是否正確安裝
2. ✅ 所有檔案是否在正確位置
3. ✅ 閱讀 `V59_USAGE_GUIDE.md` 使用指南
4. ✅ 檢查 `HANDOFF_STATUS.md` 當前進度

---

**版本**: v5.9
**日期**: 2025-12-22
**作者**: Gray (with Claude)
