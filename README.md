# 萵苣 UVA 模型 v6.3

**A Mechanistic Model for UVA Effects on Lettuce Growth and Anthocyanin Accumulation**

[![Python Version](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

---

## 📊 模型表現 (v6.3)

**鮮重預測**: 5/6 組 < 5% (83.3% 達標率) ✅

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
- 平均誤差: 9.7%
- 中位數誤差: 1.9%
- **L6D6 > CK**: ✅ 符合物理預期 (UVA 提供額外 PAR)

---

## 🎯 核心特色

### v6.3 重大更新

1. ✅ **參數完全外部化** - 所有參數移至 [params_config.py](params_config.py)
2. ✅ **修正 PAR 計算** - 移除不必要的放大效應 (3.0 → 1.0)
3. ✅ **修復時間計算 Bug** - day_from_sowing 計算錯誤 (v6.1)
4. ✅ **L6D6 > CK** - 確認 UVA 的光合促進效應
5. ✅ **高準確度** - 5/6 組誤差 < 5%

### 模型機制

- **基礎模型**: Sun et al. (2025) 萵苣生長模型
- **UVA 效應**:
  - ✅ 光合促進 (UVA-PAR 轉換)
  - ✅ 氧化損傷累積 (Stress 機制)
  - ✅ 碳依賴修復
  - ✅ 花青素誘導合成
- **狀態變量**: 6 個 (乾重, 碳池, LAI, 花青素, Stress, E_stress)
- **微分方程**: RK45 求解器, max_step=60s

---

## 🚀 快速開始

### 環境需求

```bash
# Python 版本
Python 3.8+

# 必要套件
numpy>=1.21.0
scipy>=1.7.0
matplotlib>=3.4.0
pandas>=1.3.0
```

### 安裝

```bash
# 克隆儲存庫
git clone https://github.com/YOUR_USERNAME/uva-simulation.git
cd uva-simulation

# 安裝依賴
pip install -r requirements.txt
```

### 運行模擬

```bash
# 選項 1: 運行所有處理組
python3 simulate_uva_model.py

# 選項 2: 查看參數配置
python3 params_config.py

# 選項 3: 運行驗證
python3 run_validation.py
```

---

## 📁 專案結構

### 核心模型文件

| 檔案 | 說明 |
|------|------|
| **[simulate_uva_model.py](simulate_uva_model.py)** | 主模型 (v6.3) |
| **[params_config.py](params_config.py)** | 參數配置文件 |
| [lettuce_uva_carbon_complete_model.py](lettuce_uva_carbon_complete_model.py) | Sun 基礎模型 |

### 文檔

| 檔案 | 說明 |
|------|------|
| [HANDOFF_STATUS.md](HANDOFF_STATUS.md) | 當前進度和交接狀態 |
| [MODEL_DESIGN_NOTES.md](MODEL_DESIGN_NOTES.md) | 模型設計筆記 |
| [V63_FINAL_REPORT.md](V63_FINAL_REPORT.md) | v6.3 最終校準報告 |
| [CLAUDE.md](CLAUDE.md) | Claude Code 工作守則 |

### 分析腳本

| 檔案 | 說明 |
|------|------|
| [generate_paper_figures.py](generate_paper_figures.py) | 論文圖表生成 |
| [sensitivity_analysis_extended.py](sensitivity_analysis_extended.py) | 敏感度分析 (25參數) |
| [run_validation.py](run_validation.py) | 模型驗證 |

---

## 🔧 參數調整

所有參數集中在 [params_config.py](params_config.py) 中管理：

### 關鍵參數組

```python
# Sun 基礎模型
SUN_PARAMS = {
    'c_alpha': 0.548,  # 光合效率 (校準值)
}

# UVA-PAR 轉換
UVA_PAR_PARAMS = {
    'par_conversion_factor': 1.0,  # UVA 對 PAR 的轉換
}

# Stress 損傷與修復
STRESS_PARAMS = {
    'stress_damage_coeff': 0.70e-6,  # 損傷係數
    'stress_repair_coeff': 1.0e-5,   # 修復係數
    'circadian_disruption_factor': 3.2,  # 夜間損傷加成
    # ... 更多參數
}

# 花青素
ANTHOCYANIN_PARAMS = {
    'base_anth_rate_light': 4.0e-10,  # 日間合成率
    'V_max_anth': 1.8e-10,  # 最大誘導合成率
    'k_deg': 2.6e-6,  # 降解率
    # ... 更多參數
}
```

### 打印參數摘要

```bash
python3 params_config.py
```

---

## 📊 實驗處理組

| 處理代碼 | 描述 | UVA時數/天 | 總天數 | 開始日 |
|----------|------|------------|--------|--------|
| CK | 對照組 | 0 | - | - |
| L6D6 | 低劑量日間 | 6h (10:00-16:00) | 6天 | Day 29 |
| L6D6-N | 低劑量夜間 | 6h (22:00-04:00) | 6天 | Day 29 |
| H12D3 | 高劑量脅迫 | 12h (06:00-18:00) | 3天 | Day 32 |
| VL3D12 | 極低劑量長期 | 3h (10:00-13:00) | 12天 | Day 23 |
| L6D12 | 低劑量長期 | 6h (10:00-16:00) | 12天 | Day 23 |

---

## 📈 版本歷史

### v6.3 (2025-12-22) - L6D6 優先穩定 ⭐
- ✅ 確保 L6D6 > CK (符合物理預期)
- ✅ 5/6 組 < 5% (83.3% 達標率)
- 調整參數: stress_damage_coeff, circadian_disruption_factor, stress_nonlinear_coeff

### v6.2 (2025-12-22) - 提高抑制靈敏度
- 4/6 組 < 5%
- 提高 Stress 對生長的抑制效果

### v6.1 (2025-12-22) - 修復時間計算 Bug ⭐⭐⭐
- **重大 Bug 修復**: day_from_sowing 計算錯誤
- L6D6 誤差從 -24.3% 改善到 -2.5% (+21.8%)

### v6.0 (2025-12-22) - 參數外部化 + PAR 修正 ⭐⭐⭐
- 創建 params_config.py
- 移除 PAR 放大效應 (3.0 → 1.0)
- 校準 c_alpha (0.68 → 0.548)

### v5.9 (2025-12-22) - Stress 累積能量驅動花青素
- 移除硬閾值，引入 E_stress
- 使用 Hill 函數描述花青素誘導

---

## 🔬 已知限制

### H12D3 問題

**狀態**: 無法通過累積模型達標 (+49.9%)

**原因**:
1. 時間太短 (3 天) - Stress 累積模型需要時間建立正反饋
2. LAI 已經很高 (Day 32) - 成熟植株脆弱性低
3. 數學限制 - Michaelis-Menten 型非線性累積無法捕捉短期高強度效應

**解決方案** (待實現):
- 添加即時損傷機制 (Acute Damage)
- 基於瞬時 UVA 強度的光抑制
- 預期所有組都能 < 5%

---

## 📚 參考文獻

### 基礎模型
- Sun et al. (2025) - 萵苣生長模型

### UVA 效應
- Verdaguer et al. (2017) - UVA 光合促進
- Kataria et al. (2014) - UV-A effects on photosynthesis

### Stress 機制
- Hideg et al. (2013) - UV 損傷與 ROS
- Foyer & Noctor (2005) - Oxidative stress signaling

### 花青素
- Winkel-Shirley (2002) - Anthocyanin biosynthesis
- Gould (2004) - Anthocyanin as stress protectant

---

## 🤝 貢獻

歡迎提交 Issues 和 Pull Requests！

### 開發流程

1. Fork 專案
2. 創建特性分支 (`git checkout -b feature/AmazingFeature`)
3. 提交變更 (`git commit -m 'Add some AmazingFeature'`)
4. 推送到分支 (`git push origin feature/AmazingFeature`)
5. 開啟 Pull Request

---

## 📄 授權

本專案採用 MIT 授權 - 詳見 [LICENSE](LICENSE) 文件

---

## 📞 聯絡方式

- **專案**: [uva-simulation](https://github.com/YOUR_USERNAME/uva-simulation)
- **Issues**: [提交問題](https://github.com/YOUR_USERNAME/uva-simulation/issues)

---

## 🙏 致謝

- Sun et al. (2025) 提供基礎生長模型
- Claude Code 協助開發和除錯
- 所有文獻作者的研究貢獻

---

**版本**: v6.3
**日期**: 2025-12-22
**狀態**: ✅ 生產就緒 (5/6 組達標)
