# 萵苣 UVA 模型 v10.39

**A Mechanistic Model for UVA Effects on Lettuce Growth and Anthocyanin Accumulation**

[![Python Version](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Status](https://img.shields.io/badge/status-Production%20Ready-brightgreen.svg)]()

---

## 模型表現 (v10.39)

**完美達標**: 12/12 目標達成 (訓練組 6/6 + 驗證組 6/6)

### 訓練組 (允收誤差 5%)

| 處理組 | FW預測 | FW觀測 | FW誤差 | Anth預測 | Anth觀測 | Anth誤差 |
|--------|--------|--------|--------|----------|----------|----------|
| CK | 86.5g | 87.0g | -0.5% | 439 | 433 | +1.3% |
| L6D6 | 92.5g | 91.4g | +1.2% | 474 | 494 | -4.0% |
| L6D6-N | 84.0g | 80.8g | +3.9% | 475 | 493 | -3.6% |
| VL3D12 | 69.4g | 67.0g | +3.6% | 492 | 482 | +2.0% |
| L6D12 | 58.9g | 60.4g | -2.5% | 496 | 518 | -4.3% |
| H12D3 | 61.3g | 60.6g | +1.2% | 651 | 651 | +0.0% |

### 驗證組 (允收誤差 10%)

| 處理組 | 時數 | FW預測 | FW觀測 | FW誤差 | Anth預測 | Anth觀測 | Anth誤差 |
|--------|------|--------|--------|--------|----------|----------|----------|
| CK | 0h | 86.5g | 85.2g | +1.6% | 439 | 413 | +6.2% |
| VL3D3 | 3h | 88.4g | 89.0g | -0.8% | 457 | 437 | +4.5% |
| L6D3 | 6h | 89.9g | 92.2g | -2.5% | 473 | 468 | +1.1% |
| M9D3 | 9h | 87.8g | 83.8g | +4.8% | 589 | 539 | +9.2% |
| H12D3 | 12h | 61.3g | 62.2g | -1.4% | 651 | 657 | -0.9% |
| VH15D3 | 15h | 51.2g | 51.3g | +0.0% | 532 | 578 | -7.9% |

---

## 核心參數 (v10.39)

### Gompertz 非線性因子
- threshold: 10.5 hours
- max_factor: 250.0
- steepness: 0.5

| 每日時數 | nonlinear_factor |
|----------|------------------|
| 3h | 1.0 |
| 6h | 1.0 |
| 9h | 31.1 |
| 12h | 156.9 |
| 15h | 226.0 |

### 花青素效率抑制 (Hill 函數)
- K = 800.0
- n = 1.5
- 公式: efficiency = 1 / (1 + (nonlin/K)^n)

---

## 核心特色

### 模型架構

- **基礎模型**: Sun et al. (2025) 萵苣生長模型
- **狀態變量**: 6 個 `[X_d, C_buf, LAI, Anth, Stress, ROS]`
- **ODE求解器**: RK45, max_step=300s

### UVA 效應機制

1. **UVA 形態效應** - UVA 促進 SLA 和 LAI 增長 (不直接疊加 PAR)
2. **ROS 動態** - UVA 產生 ROS，抗氧化系統清除
3. **Stress 損傷-衰減** - 累積損傷與自然衰減平衡
4. **LAI 脆弱性** - 幼嫩植株更易受損
5. **Gompertz 非線性** - 長時間照射觸發抗氧化系統崩潰
6. **夜間節律損傷** - 夜間照射造成額外壓力
7. **花青素誘導** - Stress 誘導 + UV 直接誘導
8. **水分抑制** - 極端逆境下花青素合成效率下降
9. **Hill 效率抑制** - 單調遞減抑制高 nonlinear_factor

---

## 檔案結構

```
.
├── simulate_uva_model_v10.py     # 主模型程式 (v10.39)
├── model_config.py               # 處理組設定與目標值
├── CLAUDE.md                     # 開發守則 (v3.0)
├── HANDOFF_STATUS.md             # 交接狀態
├── MODEL_DESIGN_NOTES.md         # 模型設計筆記
├── generate_paper_figures.py     # 論文圖表生成
├── optimize_uva_strategy.py      # 優化腳本 (uva_intensity=11 W/m²)
└── 紅葉萵苣UVA...論文_v5.txt     # 中文論文
```

---

## 快速開始

### 安裝

```bash
pip install numpy scipy matplotlib
```

### 運行模擬

```bash
python3 simulate_uva_model_v10.py
```

---

## 版本歷史

- **v10.39**: 單調遞減 Hill 效率函數 (K=800, n=1.5)
- **v10.37**: Gompertz threshold 9.5→10.5
- **v10.33**: 連續非對稱高斯 (已移除)
- **v10.9**: 水分抑制機制
- **v10.0**: UVA 形態效應取代直接 PAR 疊加

---

## 參考文獻

詳見 `紅葉萵苣UVA誘導花青素機制模型論文_中文版_v5.txt`
