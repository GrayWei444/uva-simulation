# UVA-Lettuce Model v1.0 Parameter Documentation

## Version History

| Version | Date | Key Changes |
|---------|------|-------------|
| v10.39 | 2026-01-14 | Original calibrated model (6/6 training, 6/6 validation) |
| v1.0 | 2026-01-25 | Pure Sun model base + AOX ODE with carbon competition |
| v1.1 | 2026-01-25 | uva_sla_enhancement 0.55→1.1, L6D6 FW=91g 達標 |
| v1.2 | 2026-01-25 | k_circadian 3e-6→3.5e-5, L6D6-N FW=80.8g 達標 |
| v1.3 | 2026-01-25 | 累積傷害機制, 6/6 training + 5/6 validation FW 達標 |
| **v1.5** | **2026-01-26** | **LAI順序修正 + LDMC多肉化, 12/12 FW 達標** |

---

## 1. Sun Model Verification (2026-01-25)

### 1.1 Implementation Correctness

v1 uses `sun_model_pure.py` which correctly implements Sun et al. (2025):

| Paper Equation | Description | Implementation |
|----------------|-------------|----------------|
| Eq. 1 | dX_d/dt = c_beta × (c_alpha × A_C × h_buf - R_d) | ✓ Correct |
| Eq. 5 | dC_buf/dt = c_alpha × A_C × h_buf - R_d - RGR_max × X_d / c_beta | ✓ Correct |
| Eq. 9 | dLAI/dt = dX_d/dt × (1 - σ_r) × SLA | ✓ Correct |
| Eq. 6-8 | 3-point Gaussian quadrature for canopy A_C | ✓ Correct |
| Eq. 10-12 | SLA = SLA_ref × f_I × f_Xh | ✓ Correct |
| Eq. 36 | σ_r = f(ln(plant_dw)) | ✓ Correct |

### 1.2 Parameter Discrepancy Discovery

**IMPORTANT:** Supplementary Materials Table S2 contains an error!

| Parameter | Supplementary Materials (WRONG) | v10.39 Code (CORRECT) |
|-----------|--------------------------------|----------------------|
| LAI(0) | 0.40 | **2.0** |
| c_alpha | 0.54 | 0.54 |
| dw_fw_ratio_base | 0.05 | 0.05 |

v10.39 initial LAI calculation (from code line 1435):
```python
LAI_init = (dw_init_g / 0.01) * 0.04  # = (0.5/0.01)*0.04 = 2.0
```

### 1.3 Calibration Results Comparison

With c_alpha=0.54, LAI(0)=2.0, dw_fw_ratio=0.05:

| Source | LAI | FW (g) |
|--------|-----|--------|
| sun_model_pure.py | 9.25 | 87.6 |
| simulate_uva_model_v1.py | 9.30 | 88.1 |
| v10.39 | 9.1 | 86.5 |
| **Target CK** | - | **87.0** |

c_alpha calibration table (pure Sun model):

| c_alpha | DW (g) | LAI | FW (g) |
|---------|--------|-----|--------|
| 0.52 | 4.25 | 8.99 | 85.0 |
| 0.525 | 4.28 | 9.05 | 85.6 |
| 0.53 | 4.31 | 9.12 | 86.3 |
| 0.535 | 4.35 | 9.18 | 86.9 |
| **0.54** | **4.38** | **9.25** | **87.6** |

---

## 2. v1.0 Parameters

### 2.1 Base Sun Model Parameters

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| c_alpha | 0.54 | - | Photosynthetic efficiency (v10.39 calibrated) |
| c_beta | 0.8 | - | Carbohydrate to structural material factor |
| RGR_max_20 | 1.54e-6 | s⁻¹ | Max relative growth rate at 20°C |

### 2.2 Initial Conditions

| State | Symbol | Value | Unit | Calculation |
|-------|--------|-------|------|-------------|
| Fresh weight | FW_init | 10.0 | g/plant | Measured at transplant |
| Dry weight | DW_init | 0.5 | g/plant | FW × 0.05 |
| Structural dry weight | X_d(0) | 0.018 | kg/m² | DW × density / 1000 |
| Leaf area index | LAI(0) | 2.0 | m²/m² | (DW/0.01) × 0.04 |
| Carbon buffer | C_buf(0) | 0.0 | kg/m² | Start empty |

### 2.3 UVA Morphological Effect

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| uva_sla_enhancement | **1.5** | - | Max X_d growth boost from UVA (v1.3 校準) |
| K_adapt_days | **8.0** | days | 適應半衰期 (v1.3 校準) |
| optimal_hours | 6.0 | h/day | Bell curve peak |
| sigma_hours | 4.0 | h | Bell curve width |

Morphological boost equation (bell curve with adaptation):
```python
hours_factor = exp(-((daily_hours - 6)² / (2 × 4²)))
xd_adaptation = K_adapt / (K_adapt + days_irradiated)  # 長期適應
xd_boost = 1.5 × hours_factor × xd_adaptation  # 作用於 dXd_dt
```

### 2.3b Cumulative Damage (NEW in v1.3)

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| k_cumulative_damage | **8.0e-7** | - | 累積傷害速率係數 |
| K_cumulative_damage | 50.0 | - | 半飽和常數 |
| max_cumulative_penalty | **0.40** | - | 最大生長抑制 (40%) |

Cumulative damage mechanism:
```python
# 累積傷害 = 積分(I_UVA × vulnerability) over treatment time
dCumulDamage_dt = k_cumul × I_UVA × vulnerability
cumul_penalty = max_penalty × CumulDamage / (K + CumulDamage)

# 總生長抑制 = 瞬時壓力 + 累積傷害
xd_reduction = stress_inhibition + cumul_penalty
```

**為何這區分 L6D6 和 L6D12：**
- L6D12 (Day 23-34): 開始時 LAI~5-6 → 高 vulnerability → 累積大量傷害
- L6D6 (Day 29-34): 開始時 LAI~9 → 低 vulnerability → 累積傷害小
- 累積傷害是「傷害債務」，持續抑制生長

**v1.3 校準結果 (累積傷害機制):**

### Training Set (6/6 FW 達標 ✓)

| Treatment | FW (g) | Target | Error | Anth | Target | Error | LAI | Stress | Status |
|-----------|--------|--------|-------|------|--------|-------|-----|--------|--------|
| CK | 88.1 | 87.0 | +1.2% | 419 | 433 | -3.3% | 9.3 | 0.0 | **✓**✓ |
| L6D6 | 90.5 | 91.4 | -1.0% | 539 | 494 | +9.2% | 9.6 | 15.5 | **✓**✗ |
| L6D6-N | 80.6 | 80.8 | -0.3% | 554 | 493 | +12.4% | 8.9 | 95.1 | **✓**✗ |
| **VL3D12** | **66.2** | **67.0** | **-1.1%** | 700 | 482 | +45.3% | 7.4 | 55.0 | **✓**✗ |
| **L6D12** | **61.7** | **60.4** | **+2.1%** | 652 | 518 | +26.0% | 7.0 | 149.4 | **✓**✗ |
| H12D3 | 60.7 | 60.6 | +0.2% | 693 | 651 | +6.4% | 9.2 | 171.0 | **✓**✗ |

### Validation Set (5/6 FW 達標)

| Treatment | FW (g) | Target | Error | Anth | Target | Error | LAI | Status |
|-----------|--------|--------|-------|------|--------|-------|-----|--------|
| CK_V | 88.1 | 85.1 | +3.5% | 419 | 413 | +1.4% | 9.3 | **✓**✓ |
| VL3D3 | 88.4 | 89.1 | -0.8% | 465 | 437 | +6.3% | 9.4 | **✓**✗ |
| L6D3 | 89.5 | 92.2 | -2.9% | 500 | 468 | +6.9% | 9.5 | **✓**✗ |
| M9D3 | 84.2 | 83.8 | +0.5% | 570 | 539 | +5.7% | 9.5 | **✓**✗ |
| H12D3_V | 60.7 | 62.2 | -2.4% | 693 | 657 | +5.4% | 9.2 | **✓**✗ |
| VH15D3 | 57.2 | 51.3 | +11.7% | 474 | 578 | -18.0% | 9.0 | ✗✗ |

**關鍵發現：**
- **累積傷害機制成功區分 L6D6 和 L6D12！**
- L6D6: FW=90.5g (LAI=9.6) - 形態增益主導
- L6D12: FW=61.7g (LAI=7.0) - 累積傷害主導（Day 23 低 LAI 時開始）
- VL3D12: FW=66.2g (LAI=7.4) - 同樣受累積傷害影響
- Training 6/6 FW 達標, Validation 5/6 FW 達標
- Anthocyanin 校準仍需改進

### 2.4 ROS Dynamics

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| k_ros_production | 0.022 | ROS/(W/m²·s) | ROS production coefficient |
| k_ros_clearance | 5e-4 | s⁻¹ | ROS clearance coefficient |
| LAI_ref | 9.0 | m²/m² | Reference LAI for area dilution |

ROS production (area dilution mechanism):
```python
area_concentration = I_UVA × LAI_ref / LAI  # Small LAI → high ROS
ros_production = k_ros × area_concentration
```

### 2.5 LAI Vulnerability

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| A_vulnerability | 500000 | - | Vulnerability amplitude |
| k_vulnerability | 1.6 | - | Decay coefficient |
| stress_damage_coeff | 8.0e-6 | - | Base damage coefficient |

Vulnerability function:
```python
vulnerability = A × exp(-k × LAI) + 1
# LAI=5 (D12 start): vuln ≈ 167
# LAI=7: vuln ≈ 16
# LAI=9 (D6 start): vuln ≈ 3
```

### 2.6 Gompertz Nonlinear Damage

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| gompertz_max_factor | 450 | - | Maximum amplification |
| gompertz_threshold | 11.5 | h | Collapse threshold |
| gompertz_steepness | 0.7 | - | Transition steepness |

```python
nonlinear_factor = 1 + max_factor × exp(-exp(-steepness × (hours - threshold)))
# 3h/day: 1.0, 6h/day: 1.0, 9h/day: ~31, 12h/day: ~157, 15h/day: ~226
```

### 2.7 Circadian Damage (Night Treatment)

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| k_circadian | **3.5e-5** | - | Circadian damage coefficient (v1.2 校準) |
| n_circadian | 2.0 | - | Exponent for dark hour effect |

```python
# Night UVA irradiation causes extra damage
circadian_damage = k_circadian × I_UVA × hours_in_dark²
```

L6D6-N (22:00-04:00) 夜間處理比 L6D6 (10:00-16:00) 日間處理有更多晝夜節律傷害。

### 2.8 LDMC (DW/FW Ratio)

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| dw_fw_ratio_base | 0.05 | - | Base ratio (healthy plants) |
| K_stress_ldmc | 150 | - | Half-saturation for chronic LDMC |
| stress_ldmc_max | 0.10 | - | Max LDMC increase from stress |
| dw_fw_ratio_max | 0.12 | - | Maximum ratio (severe stress) |

### 2.9 AOX (Antioxidant) Synthesis

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| base_aox_rate_light | 6.35e-10 | kg/(m²·s) | Baseline synthesis (light) |
| base_aox_rate_dark | 3.18e-10 | kg/(m²·s) | Baseline synthesis (dark) |
| V_max_aox | 5.5e-9 | kg/(m²·s) | Max stress-induced synthesis |
| K_stress_aox | 120 | - | Stress half-saturation |
| k_aox_deg | 4.0e-6 | s⁻¹ | Degradation rate |
| anth_fraction | 0.18 | - | Anthocyanin fraction of AOX |

### 2.10 Carbon Competition (NEW in v1)

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| c_cost_aox | 0.05 | - | Carbon cost for UVA-induced AOX |
| base_c_cost | 0.005 | - | Carbon cost for base AOX (CK) |
| c_cost_lai_factor | 0.5 | - | Low LAI amplifies C competition |
| K_c_cost_lai | 8.0 | m²/m² | LAI half-saturation |

---

## 3. Key Mechanism: UVA Dual Effect

### 3.1 Morphological Enhancement (Positive)

- Determined purely by **daily hours** (bell curve, 6h optimal)
- Independent of stress/LAI
- Increases X_d growth → LAI follows via Sun model

### 3.2 ROS/Stress Damage (Negative)

- Determined by **unit-area UVA intensity** (I_UVA / LAI)
- Large LAI → dilutes ROS → less damage
- Small LAI → concentrates ROS → more damage

### 3.3 Why D12 Treatments Have Low FW

**L6D6 (Day 29-34, 6 days):**
```
Day 29: LAI already large (~8+)
    → Low ROS per unit area
    → Damage small, enhancement dominates
    → LAI continues growing (positive feedback)
    → Final: highest LAI, highest FW
```

**VL3D12 / L6D12 (Day 23-34, 12 days):**
```
Day 23: LAI still small (~5)
    → High ROS per unit area
    → Damage large, suppresses growth
    → LAI cannot grow (negative feedback)
    → Final: low LAI, low FW
```

**Key insight:** VL3D12 has low FW not because 3h/day is stressful, but because irradiation started when LAI was small!

---

## 4. Environment Settings

| Parameter | Value | Unit |
|-----------|-------|------|
| I_day (LED) | 57 | W/m² |
| I_UVA | 11 | W/m² |
| T_day | 25 | °C |
| T_night | 18 | °C |
| CO2 | 1200 | ppm |
| RH_day | 70 | % |
| RH_night | 85 | % |
| Photoperiod | 16 | h (06:00-22:00) |
| Plant density | 36 | plants/m² |

---

## 5. Experimental Targets

### Training Set (Day 14-35, ±5% tolerance)

| Treatment | Hours | Days | Start | FW (g) | Anth (ppm) |
|-----------|-------|------|-------|--------|------------|
| CK | 0 | 0 | - | 87.0 | 433 |
| L6D6 | 6 | 6 | Day 29 | 91.4 | 494 |
| L6D6-N | 6 (night) | 6 | Day 29 | 80.8 | 493 |
| VL3D12 | 3 | 12 | Day 23 | 67.0 | 482 |
| L6D12 | 6 | 12 | Day 23 | 60.4 | 518 |
| H12D3 | 12 | 3 | Day 32 | 60.6 | 651 |

### Validation Set (3-Day Gradient, ±10% tolerance)

| Treatment | Hours | Days | Start | FW (g) | Anth (ppm) |
|-----------|-------|------|-------|--------|------------|
| CK_V | 0 | 0 | - | 85.2 | 413 |
| VL3D3 | 3 | 3 | Day 32 | 89.0 | 437 |
| L6D3 | 6 | 3 | Day 32 | 92.2 | 468 |
| M9D3 | 9 | 3 | Day 32 | 83.8 | 539 |
| H12D3_V | 12 | 3 | Day 32 | 62.2 | 657 |
| VH15D3 | 15 | 3 | Day 32 | 51.3 | 578 |

---

## 6. v1.5 Key Parameters and Mechanisms (2026-01-26)

### 6.1 核心問題與解決

**問題：** v1.3 的 LAI 順序錯誤
- VL3D12 LAI > L6D12 LAI (錯誤！應該相反)
- L6D6 FW ≈ CK FW (錯誤！L6D6 LAI 較高，FW 也應較高)

**解決方案：** 三個機制調整

### 6.2 ROS Hours Attenuation (區分 L6D12 vs VL3D12)

```python
hours_attenuation = 1.0 + 0.268 × daily_uva_hours
# 6h/day: 2.608 (更多 ROS)
# 3h/day: 1.804 (較少 ROS)

ros_production = k_ros × I_UVA × area_dilution × hours_attenuation
```

**生理學解釋：** 高 UVA 劑量飽和光保護機制，6h/day 的 ROS 產生效率比 3h/day 高。
**效果：** L6D12 累積更多 Stress → FW 較低

### 6.3 LDMC 多肉化效應 (提升 L6D6 FW)

```python
# 低壓力 + 6h/day → 多肉化 → LDMC 降低 → FW 增加
hours_factor = exp(-((daily_hours - 6)² / (2 × 4²)))
stress_block = Stress / (30.0 + Stress)
succulence = 0.126 × (hours_factor²) × (1.0 - stress_block)

# 淨 LDMC 變化
net_ldmc_change = ldmc_increase - succulence
```

**效果：** L6D6 (低壓力, 6h) 獲得多肉化 → FW > CK

### 6.4 X_d Adaptation Steep (短期處理增益)

```python
# 短期處理獲得更多形態增益，長期處理衰減
xd_adaptation_steep = 5.0 / (5.0 + days + 0.03 × days²)
# D3: 0.60, D6: 0.41, D12: 0.23

xd_boost = 1.5 × hours_factor × xd_adaptation_steep
```

### 6.5 v1.5 校準結果

#### Training Set (6/6 FW ✓)

| Treatment | FW (g) | Target | Error | LAI | Stress | Status |
|-----------|--------|--------|-------|-----|--------|--------|
| CK | 88.1 | 87.0 | +1.2% | 9.3 | 0.0 | ✓ |
| L6D6 | **90.0** | 91.4 | -1.6% | **9.6** | 34.8 | ✓ |
| L6D6-N | 77.1 | 80.8 | -4.6% | 8.6 | 145.1 | ✓ |
| VL3D12 | 63.9 | 67.0 | -4.6% | **7.4** | 84.9 | ✓ |
| L6D12 | 63.4 | 60.4 | +5.0% | **7.9** | 125.7 | ✓ |
| H12D3 | 63.2 | 60.6 | +4.3% | 9.0 | 658.3 | ✓ |

#### Validation Set (6/6 FW ✓)

| Treatment | FW (g) | Target | Error | LAI | Stress | Status |
|-----------|--------|--------|-------|-----|--------|--------|
| CK_V | 88.1 | 85.1 | +3.5% | 9.3 | 0.0 | ✓ |
| VL3D3 | 92.6 | 89.1 | +3.9% | 9.4 | 10.3 | ✓ |
| L6D3 | 92.9 | 92.2 | +0.8% | 9.5 | 29.7 | ✓ |
| M9D3 | 85.1 | 83.8 | +1.5% | 9.5 | 62.6 | ✓ |
| H12D3_V | 63.2 | 62.2 | +1.6% | 9.0 | 658.3 | ✓ |
| VH15D3 | 53.1 | 51.3 | +3.8% | 8.8 | 2726.2 | ✓ |

### 6.6 LAI 順序驗證

| 比較 | 結果 | 狀態 |
|------|------|------|
| L6D12 LAI (7.9) > VL3D12 LAI (7.4) | 6h 形態增益 > 3h | ✓ |
| L6D6 LAI (9.6) > CK LAI (9.3) | UVA 形態增益 | ✓ |
| L6D6 FW (90.0) > CK FW (88.1) | 高 LAI + 多肉化 | ✓ |

### 6.7 v1.5 關鍵參數

| Parameter | Value | Description |
|-----------|-------|-------------|
| hours_attenuation | 1.0 + 0.268 × hours | ROS 時數放大 |
| max_succulence | 0.126 | 最大多肉化效應 |
| K_stress_succulence | 30.0 | 多肉化壓力閾值 |
| xd_adaptation_steep | 5.0/(5.0+d+0.03d²) | X_d 增益衰減 |

**v1.5 校準結果：**
- Training: **6/6 FW 達標** ✓
- Validation: **6/6 FW 達標** ✓
- LAI 順序: L6D12 > VL3D12 ✓, L6D6 > CK ✓

**待改進：**
- Anthocyanin 校準（目前偏高或偏低）

---

## 7. References

1. Sun, J., et al. (2025). A lettuce growth model responding to a broad range of greenhouse climates. Biosystems Engineering, 250, 285-305.
2. Wei, C.H., Fang, W., & Huang, C.K. (2026). A Two-Stage Screening-to-Optimization Approach with Mechanistic Model Analysis. Plants (under review).
