# Handoff Status

**Last Updated**: 2026-01-17
**Current Version**: v10.39

---

## Core Formulas (100% Consistent with Code)

### 1. Gompertz Nonlinear Factor

```python
nonlinear_factor = 1 + max_factor * np.exp(-np.exp(-steepness * (hours - threshold)))
```

Parameters:
- `threshold = 10.5` hours
- `max_factor = 250.0`
- `steepness = 0.5`

### 2. Hill Efficiency Inhibition

```python
efficiency = 1.0 / (1.0 + (nonlinear_factor / K) ** n)
```

Parameters:
- `K = 800.0`
- `n = 1.5`

### 3. LDMC (Leaf Dry Matter Content)

```python
stress_effect = ldmc_sensitivity * Stress / (K_ldmc + Stress)
x = softplus((nonlin - acute_center) / acute_scale) * acute_scale
acute_factor = 1 + acute_k * x^acute_n / (acute_K^acute_n + x^acute_n)
LDMC = base * (1 + stress_effect * acute_factor)
```

Parameters:
- `dw_fw_ratio_base = 0.05`
- `ldmc_stress_sensitivity = 0.45`
- `K_ldmc = 1400`
- `dw_fw_ratio_max = 0.080`
- `acute_center = 50`, `acute_scale = 10`
- `acute_k = 9`, `acute_K = 120`, `acute_n = 2`

---

## Simulation Results (v10.39)

### Training Set (Tolerance: 5%)

| Treatment | FW Obs | FW Pred | FW Error | Anth Obs | Anth Pred | Anth Error |
|-----------|--------|---------|----------|----------|-----------|------------|
| CK | 87.0 | 86.5 | -0.5% ✓ | 433 | 439 | +1.3% ✓ |
| L6D6 | 91.4 | 92.5 | +1.2% ✓ | 494 | 474 | -4.0% ✓ |
| L6D6-N | 80.8 | 84.0 | +3.9% ✓ | 493 | 475 | -3.6% ✓ |
| VL3D12 | 67.0 | 69.4 | +3.6% ✓ | 482 | 492 | +2.0% ✓ |
| L6D12 | 60.4 | 58.9 | -2.5% ✓ | 518 | 496 | -4.3% ✓ |
| H12D3 | 60.6 | 61.3 | +1.2% ✓ | 651 | 651 | +0.0% ✓ |

**Training: FW 6/6, Anth 6/6, Total 12/12**

### Validation Set (Tolerance: 10%)

| Treatment | Hours | FW Obs | FW Pred | FW Error | Anth Obs | Anth Pred | Anth Error |
|-----------|-------|--------|---------|----------|----------|-----------|------------|
| CK | 0h | 85.2 | 86.5 | +1.6% ✓ | 413 | 439 | +6.2% ✓ |
| VL3D3 | 3h | 89.0 | 88.4 | -0.8% ✓ | 437 | 457 | +4.5% ✓ |
| L6D3 | 6h | 92.2 | 89.9 | -2.5% ✓ | 468 | 473 | +1.1% ✓ |
| M9D3 | 9h | 83.8 | 87.8 | +4.8% ✓ | 539 | 589 | +9.2% ✓ |
| H12D3 | 12h | 62.2 | 61.3 | -1.4% ✓ | 657 | 651 | -0.9% ✓ |
| VH15D3 | 15h | 51.3 | 51.2 | +0.0% ✓ | 578 | 532 | -7.9% ✓ |

**Validation: FW 6/6 (<5%), Anth 6/6 (<10%)**

---

## Nonlinear Factor Reference

| Hours/Day | nonlinear_factor | Efficiency |
|-----------|------------------|------------|
| 3h | 1.0 | 100.0% |
| 6h | 1.0 | 100.0% |
| 9h | 31.1 | 99.2% |
| 12h | 156.9 | 92.0% |
| 15h | 226.0 | 86.9% |

---

## File Synchronization Status

| File | Status |
|------|--------|
| simulate_uva_model_v10.py | v10.39 ✓ |
| model_config.py | Synchronized ✓ |
| lettuce_uva_carbon_complete_model.py | v7.1 ✓ |
| CLAUDE.md | v3.0 English ✓ |
| HANDOFF_STATUS.md | Synchronized ✓ |
| MODEL_DESIGN_NOTES.md | Synchronized ✓ |
| README.md | English ✓ |
| generate_paper_figures.py | v10.39 ✓ |
| optimize_uva_strategy.py | 11 W/m² ✓ |

---

## Key Parameters Summary

### UVA Morphological Effect
- `uva_sla_enhancement = 5.0`
- `K_uva_sla = 7.5` W/m²
- `uva_lai_boost = 1.70`
- `K_uva_lai = 7.5` W/m²

### ROS Dynamics
- `k_ros_production = 0.010`
- `k_ros_clearance = 5e-4`

### Stress Damage
- `stress_damage_coeff = 1.6e-7`
- `A_vulnerability = 8.5e7`
- `k_vulnerability = 2.0`
- `k_circadian = 3.0e-6`
- `n_circadian = 2.0`

### Stress Decay & Growth Inhibition
- `k_stress_decay = 2.14e-5` (half-life ≈ 9 hours)
- `stress_photosynthesis_inhibition = 0.85`
- `K_stress = 50.0`

### Anthocyanin Synthesis
- `base_anth_rate_light = 6.35e-10`
- `V_max_anth = 2.75e-9`
- `K_stress_anth = 100.0`
- `k_deg = 3.02e-6`

### Environment
- `I_day = 57` W/m² (PPFD 130 μmol/m²/s equivalent)
- `I_UVA = 11` W/m²
- `plant_density = 36` plants/m²

---

## Update History

### 2026-01-17
- Converted all code comments to English
- Updated README.md, CLAUDE.md, HANDOFF_STATUS.md to English
- Removed unnecessary test/tune files
- Prepared for GitHub commit

### 2026-01-16
- Fixed DOI bibliography errors (Jiménez → Budkowska)
- Added complete A_canopy formulas to paper
- Expanded equation (8) with 11 sub-formulas
- Added acute factor parameters to Table 4

### 2026-01-14
- v10.39: Hill efficiency function (K=800, n=1.5)
- Replaced sigmoid with monotonically decreasing Hill function
- All 12/12 targets achieved
