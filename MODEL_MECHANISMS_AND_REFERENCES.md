# 萵苣UVA模型機制與參考文獻整理

**最後更新**: 2025-12-04
**用途**: 論文撰寫參考

---

## 1. 模型機制總覽

### 1.1 已實現的機制

| 機制 | 狀態變量 | 用於解釋處理組 | 參考文獻 |
|------|----------|----------------|----------|
| 基礎生長模型 | X_d, C_buf, LAI | CK (對照組) | Sun et al. (2025) |
| UVA-PAR增益 | - | L6D6 (日間UVA) | Yao et al. (2020) |
| LDMC效應 | D_UVA | L6D6, L6D6-N | Qian et al. (2021) |
| 花青素基礎合成 | Anth | CK | Petroni & Tonelli (2011) |
| UVA誘導花青素 | Anth | L6D6 | Lee et al. (2014) |
| 生理節律抑制 | CircInhib | L6D6-N (夜間UVA) | Bennie et al. (2016), Deng et al. (2025) |
| 適應效應 | Adapt | VL3D12, L6D12 (長期照射) | Hideg et al. (2013) |
| 損傷機制 | Stress | H12D3 (高劑量) | Kakani et al. (2003) |

### 1.2 狀態變量定義

```
X_d      : 乾物質密度 [kg/m²]
C_buf    : 碳緩衝池 [kg C/m²]
LAI      : 葉面積指數 [-]
Anth     : 花青素含量 [g/m²]
D_UVA    : 累積UVA劑量 [W·s/m²]
Stress   : 脅迫程度 [-]
Adapt    : 適應程度 [0-0.65]
CircInhib: 生理節律抑制程度 [0-0.95]
```

---

## 2. 各機制詳細說明與文獻支持

### 2.1 基礎生長模型 (Sun Model)

**機制**: 三狀態變量ODE模型，描述乾物質累積、碳緩衝和葉面積擴展

**方程式**:
```
dX_d/dt = Growth_rate - Respiration
dC_buf/dt = Assimilation - Growth - Maintenance
dLAI/dt = LAI_growth_rate
```

**參考文獻**:
- **Sun, J. et al. (2025).** [論文待確認] 萵苣植物工廠碳模型

---

### 2.2 UVA-PAR增益效應 (UVA-PAR Enhancement)

**機制**: UVA光被葉片感知，可增加光合作用效率

**生物學基礎**:
- 萵苣光飽和點約 250 μmol/m²/s
- 植物工廠環境 (130 μmol/m²/s) 遠未飽和
- UVA可被光受體吸收，轉換為等效PAR

**方程式**:
```
I_effective = I_PAR + I_UVA × par_conversion_factor
```
其中 `par_conversion_factor ≈ 14.0 μmol/m²/s per W/m²`

**參考文獻**:
- **Yao, Y. et al. (2020).** Effects of ultraviolet-A on photosynthesis and photoinhibition in lettuce. *J. Plant Physiol.* DOI: 待確認
- **Viršilė, A. et al. (2020).** LED Light Supplementation Effects on Photosynthesis, Growth and Phenolic Content of Lettuce. *Scientia Horticulturae* 262: 109022.

**驗證**: L6D6 (日間UVA + PAR) FW比L6D6-N (夜間UVA) 高7%

---

### 2.3 LDMC效應 (Leaf Dry Matter Content)

**機制**: UV脅迫使葉片變厚變乾，DW:FW比例升高

**生物學基礎**:
- UV光誘導細胞壁增厚
- 植物增加葉片厚度以減少UV穿透
- 這是植物的防禦反應

**方程式**:
```
dw_fw_ratio = dw_fw_ratio_base + dw_fw_ratio_increase_per_dose × D_UVA

# 非線性效應 (高劑量時)
if D_UVA > threshold:
    dw_fw_ratio += additional_increase
```

**參考文獻**:
- **Qian, M. et al. (2021).** UVA and UVB light independently affect photosynthetic electron transport, membrane potential and leaf morphology in Arabidopsis. *Plant Physiol.* DOI: 10.1093/plphys/kiab262
- **Wargent, J.J. & Jordan, B.R. (2013).** From ozone depletion to agriculture: understanding the role of UV radiation in sustainable crop production. *New Phytologist* 197: 1058-1076.

**驗證**: L6D6-N (有UVA但無PAR增益) FW較CK低2%

---

### 2.4 花青素基礎合成 (Basal Anthocyanin Synthesis)

**機制**: LED藍光和溫度效應誘導花青素基礎合成

**生物學基礎**:
- 藍光 (450-500nm) 活化 CRY 和 phototropin
- 光受體活化 MYB-bHLH-WD40 轉錄因子複合物
- 溫度調控花青素穩定性

**方程式**:
```
base_synthesis = base_anth_rate_light (日間)
              或 base_anth_rate_dark  (夜間)
```

**參考文獻**:
- **Petroni, K. & Tonelli, C. (2011).** Recent advances on the regulation of anthocyanin synthesis in reproductive organs. *Plant Science* 181: 219-229.
- **Albert, N.W. et al. (2014).** Light-induced vegetative anthocyanin pigmentation in Petunia. *J. Exp. Bot.* 65: 4609-4619.

**驗證**: CK組在無UVA下仍有 43 ppm 花青素

---

### 2.5 UVA誘導花青素合成 (UVA-Induced Anthocyanin)

**機制**: UVA劑量依賴的花青素額外誘導

**生物學基礎**:
- UVR8 光受體感知 UV (主要是 UVB，但 UVA 也有效)
- 活化 HY5 轉錄因子
- 增強 PAL, CHS, DFR, ANS 等生物合成酶表達

**方程式** (Hill 方程):
```
uva_induced = S_max × (D^n) / (K_m^n + D^n)

其中:
  S_max = 最大誘導速率
  K_m   = 半飽和常數
  n     = Hill 係數
  D     = 累積UVA劑量
```

**參考文獻**:
- **Lee, M.J. et al. (2014).** Ultraviolet-A radiation stimulates anthocyanin accumulation in some red leaf lettuce cultivars. *J. Amer. Soc. Hort. Sci.* 139: 484-489.
- **Rizzini, L. et al. (2011).** Perception of UV-B by the Arabidopsis UVR8 protein. *Science* 332: 103-106.
- **Jenkins, G.I. (2017).** Photomorphogenic responses to ultraviolet-B light. *Plant Cell Environ.* 40: 2544-2557.

**驗證**: L6D6 比 CK 多約 6-7 ppm 花青素

---

### 2.6 生理節律抑制效應 (Circadian Disruption)

**機制**: 夜間UVA照射打斷植物暗期恢復，抑制花青素合成

**生物學基礎**:
- 植物具有內源性生理節律時鐘
- 花青素合成酶 (如 PAL, CHS) 受生理節律調控
- 夜間光照打斷暗期恢復，抑制相關酶的表達
- 夜間脅迫造成額外的氧化損傷

**方程式**:
```
# 夜間UVA照射時累積抑制
if is_night_uva:
    dCircInhib/dt = inhibition_rate × (1 - CircInhib)
else:
    dCircInhib/dt = -decay_rate × CircInhib

# 抑制效應
base_synthesis = base_synthesis × (1 - min(CircInhib, 0.95))
```

**參考文獻**:
- **Bennie, J. et al. (2016).** Ecological effects of artificial light at night on wild plants. *J. Ecology* 104: 611-620. DOI: 10.1111/1365-2745.12551
- **Deng, Y. et al. (2025).** Effects of light supplementation at night on plant growth and physiological responses. *Biology* 14: 571. DOI: 10.3390/biology14050571
- **Harmer, S.L. (2009).** The circadian system in higher plants. *Annu. Rev. Plant Biol.* 60: 357-377.
- **Covington, M.F. et al. (2008).** Global transcriptome analysis reveals circadian regulation of key pathways in plant growth and development. *Genome Biol.* 9: R130.

**驗證**: L6D6-N (夜間) 花青素比 L6D6 (日間) 低約 5%

---

### 2.7 適應效應 (Adaptation/Acclimation)

**機制**: 長期UVA照射使植物產生適應，減少效應強度

**生物學基礎**:
- 長期UV曝露後，植物上調保護機制
- 抗氧化酶系統增強
- UV吸收化合物 (如黃酮類) 累積
- 對後續UV的敏感度降低

**方程式**:
```
if I_UVA > 0:
    dAdapt/dt = adaptation_rate × I_UVA × (1 - Adapt/max_adaptation)
else:
    dAdapt/dt = -decay_rate × Adapt

# 適應降低效應強度
effective_response = response × (1 - Adapt)
```

**參考文獻**:
- **Hideg, E. et al. (2013).** UV-B exposure, ROS, and stress: inseparable companions or loosely linked associates? *Trends Plant Sci.* 18: 107-115.
- **Jansen, M.A.K. et al. (1998).** Higher plants and UV-B radiation: balancing damage, repair and acclimation. *Trends Plant Sci.* 3: 131-135.
- **Brosché, M. & Strid, A. (2003).** Molecular events following perception of ultraviolet-B radiation by plants. *Physiol. Plant.* 117: 1-10.

**驗證**: VL3D12 (12天照射) 比 L6D6 (6天) 每日效應更弱

---

### 2.8 損傷機制 (Damage Mechanism)

**機制**: 高劑量UVA造成光氧化損傷，抑制生長

**生物學基礎**:
- 過量UV產生活性氧 (ROS)
- ROS 損傷葉綠體和光合系統
- 蛋白質、DNA、膜脂質受損
- 植物資源從生長轉向修復

**方程式**:
```
dStress/dt = stress_accumulation × I_UVA - recovery × Stress

# 損傷降低生長
growth_rate = base_growth × (1 - damage_sensitivity × Stress)
```

**參考文獻**:
- **Kakani, V.G. et al. (2003).** Effects of ultraviolet-B radiation on cotton (Gossypium hirsutum L.) morphology and anatomy. *Ann. Bot.* 91: 817-826.
- **Frohnmeyer, H. & Staiger, D. (2003).** Ultraviolet-B radiation-mediated responses in plants: balancing damage and protection. *Plant Physiol.* 133: 1420-1428.
- **A-H-Mackerness, S. et al. (2001).** Ultraviolet-B-induced stress and changes in gene expression in Arabidopsis thaliana: role of signalling pathways controlled by jasmonic acid, ethylene and reactive oxygen species. *Plant Cell Environ.* 24: 1373-1381.

**用途**: 用於 H12D3 (12h/day × 3天) 高劑量組，該組 FW 下降 33%

---

## 3. 驗證策略

### 3.1 資料集分割

| 資料集 | 處理組 | 目的 |
|--------|--------|------|
| 訓練集 (Training) | CK, L6D6, L6D6-N | 參數估計 |
| 驗證集 (Validation) | VL3D12, L6D12 | 驗證時間效應泛化 |
| 測試集 (Test) | H12D3, SIN, INT | 驗證極端條件和動態模式 |

### 3.2 各處理組驗證的機制

| 處理組 | 主要驗證機制 |
|--------|--------------|
| CK | 基礎生長模型、基礎花青素合成 |
| L6D6 | UVA-PAR增益、UVA誘導花青素、LDMC效應 |
| L6D6-N | 生理節律抑制、無PAR增益 (夜間) |
| VL3D12 | 適應效應、低日劑量效率 |
| L6D12 | 累積脅迫、長期LDMC效應 |
| H12D3 | 損傷機制、超補償花青素合成 |
| SIN | 動態照射模式 (正弦波) |
| INT | 間歇式照射 (恢復時間效應) |

---

## 4. 論文結構建議

### 4.1 Materials and Methods

1. **Plant Material and Growth Conditions**
2. **UVA Treatment Design**
3. **Mathematical Model**
   - 3.1 Base Growth Model
   - 3.2 UVA-PAR Enhancement
   - 3.3 LDMC Response
   - 3.4 Anthocyanin Synthesis
   - 3.5 Circadian Disruption (Night UVA)
   - 3.6 Adaptation/Acclimation
   - 3.7 Damage Mechanism
4. **Parameter Estimation**
   - Training/Validation Split
   - Weighted Least Squares Optimization
   - Bootstrap Confidence Intervals
5. **Model Validation**

### 4.2 需要的圖表

1. **Fig 1**: 模型結構示意圖 (狀態變量和機制關係)
2. **Fig 2**: 訓練集擬合結果 (CK, L6D6, L6D6-N)
3. **Fig 3**: 驗證集預測結果 (VL3D12, L6D12)
4. **Fig 4**: 參數敏感度分析
5. **Fig 5**: 模型預測 vs 實驗數據 (所有處理組)

### 4.3 表格

1. **Table 1**: 模型參數估計值 (含95%置信區間)
2. **Table 2**: 模型擬合統計量 (R², RMSE, AIC, BIC)
3. **Table 3**: 各處理組預測誤差

---

## 5. 參考文獻列表 (完整格式)

### 核心引用

1. Bennie, J., Davies, T.W., Cruse, D., & Gaston, K.J. (2016). Ecological effects of artificial light at night on wild plants. *Journal of Ecology*, 104(3), 611-620. https://doi.org/10.1111/1365-2745.12551

2. Deng, Y., et al. (2025). Effects of light supplementation at night on plant growth and physiological responses. *Biology*, 14(5), 571. https://doi.org/10.3390/biology14050571

3. Hideg, E., Jansen, M.A.K., & Strid, A. (2013). UV-B exposure, ROS, and stress: inseparable companions or loosely linked associates? *Trends in Plant Science*, 18(2), 107-115.

4. Lee, M.J., Son, J.E., & Oh, M.M. (2014). Ultraviolet-A radiation stimulates anthocyanin accumulation in some red leaf lettuce cultivars. *Journal of the American Society for Horticultural Science*, 139(4), 484-489.

5. Petroni, K., & Tonelli, C. (2011). Recent advances on the regulation of anthocyanin synthesis in reproductive organs. *Plant Science*, 181(3), 219-229.

6. Qian, M., et al. (2021). UVA and UVB light independently affect photosynthetic electron transport, membrane potential and leaf morphology in Arabidopsis. *Plant Physiology*, kiab262. https://doi.org/10.1093/plphys/kiab262

7. Rizzini, L., et al. (2011). Perception of UV-B by the Arabidopsis UVR8 protein. *Science*, 332(6025), 103-106.

8. Wargent, J.J., & Jordan, B.R. (2013). From ozone depletion to agriculture: understanding the role of UV radiation in sustainable crop production. *New Phytologist*, 197(4), 1058-1076.

---

**此文件應定期更新，確保參考文獻完整性。**
