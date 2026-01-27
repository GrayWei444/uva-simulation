# UVA Lettuce Model v2.0 - Parameter Documentation

**Version:** v2.0 (Carbon Competition + AOX Framework)
**Date:** 2026-01-26
**Base:** v10.39 (all mechanisms have literature support)

---

## Core Innovation

v2.0 introduces **real carbon competition** between growth and defense:

1. **State variable change:** Anthocyanin (Anth) → AOX (Total Antioxidants)
2. **Output:** Anthocyanin = AOX × 18%
3. **Carbon competition:** AOX synthesis consumes C_buf (10% per timestep max)
4. **Growth-defense tradeoff:** Stress-induced defense diverts carbon from growth

---

## Literature Support

| Mechanism | Reference | DOI |
|-----------|-----------|-----|
| Growth-Differentiation Balance | Herms & Mattson (1992) | 10.1086/285343 |
| Resource allocation tradeoffs | Monson et al. (2022) | 10.1111/nph.17773 |
| Phenylpropanoid carbon flux | Vogt (2010) | 10.1093/mp/ssp106 |
| Metabolic costs | Gershenzon (1994) | 10.1007/BF02059810 |

---

## State Variables (6 total)

| Variable | Description | Unit |
|----------|-------------|------|
| X_d | Dry weight | kg/m² |
| C_buf | Carbon buffer pool | kg/m² |
| LAI | Leaf Area Index | m²/m² |
| AOX | Antioxidant content | kg/m² |
| Stress | Cumulative stress index | - |
| ROS | Reactive Oxygen Species | - |

---

## Model Parameters

### 1. Base Photosynthesis

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| c_alpha | 0.54 | - | Photosynthesis efficiency coefficient |

### 2. UVA Morphological Effect

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| uva_sla_enhancement | 5.00 | - | UVA max enhancement on SLA |
| K_uva_sla | 7.5 | W/m² | SLA enhancement half-saturation |
| uva_lai_boost | 1.70 | - | UVA max enhancement on LAI growth |
| K_uva_lai | 7.5 | W/m² | LAI enhancement half-saturation |

### 3. ROS Dynamics

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| k_ros_production | 0.010 | ROS/(W/m²·s) | ROS production coefficient |
| k_ros_clearance | 5e-4 | 1/s | ROS clearance coefficient |

### 4. Stress Damage

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| stress_damage_coeff | 1.6e-7 | - | Base damage coefficient |
| A_vulnerability | 8.5e7 | - | Vulnerability amplitude |
| k_vulnerability | 2.0 | 1/(m²/m²) | Vulnerability decay coefficient |
| gompertz_max_factor | 250.0 | - | Gompertz maximum factor |
| gompertz_threshold | 10.5 | h | Gompertz threshold (daily hours) |
| gompertz_steepness | 0.5 | - | Gompertz steepness |
| alpha_aox_protection | 0.5 | - | Maximum AOX protection efficiency |
| K_aox_protection | 2.78e-5 | kg/m² | AOX protection half-saturation |
| k_circadian | 3.0e-6 | - | Circadian damage coefficient |
| n_circadian | 2.0 | - | Circadian damage exponent |
| k_nonlinear_stress | 5.0e-6 | - | Nonlinear stress coefficient |

### 5. Stress Decay

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| k_stress_decay | 2.14e-5 | 1/s | Stress decay (half-life ~9 hours) |

### 6. Stress Inhibition on Growth

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| stress_photosynthesis_inhibition | 0.85 | - | Max photosynthesis inhibition |
| stress_lai_inhibition | 0.80 | - | Max LAI growth inhibition |
| K_stress | 50.0 | - | Stress half-saturation |

### 7. AOX Synthesis

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| base_aox_rate_light | 3.53e-9 | kg/m²/s | Base synthesis (day) |
| base_aox_rate_dark | 1.77e-9 | kg/m²/s | Base synthesis (night) |
| V_max_aox | 1.45e-8 | kg/m²/s | Max stress-induced synthesis |
| K_stress_aox | 100.0 | - | Stress half-saturation for synthesis |
| k_aox_deg | 3.02e-6 | 1/s | AOX degradation rate |

### 8. Carbon Competition

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| aox_carbon_cost | 1.0 | kg C/kg AOX | Carbon cost per AOX synthesized |
| carbon_competition_K | 1e-8 | kg/m²/s | Half-saturation for AOX carbon effect |
| stress_competition_K | 21.0 | - | Half-saturation for stress effect |
| stress_competition_max | 0.225 | - | Max stress-based competition |
| carbon_competition_max | 0.30 | - | Max AOX-based competition |
| max_cbuf_consumption | 0.10 | - | Max fraction of C_buf consumed per timestep |

**Carbon Competition Formula:**
```
stress_carbon_effect = 0.225 × Stress / (21 + Stress)
aox_carbon_effect = stress_aox_carbon_demand / (1e-8 + stress_aox_carbon_demand)
carbon_competition_effect = aox_carbon_effect × 0.30 + stress_carbon_effect
growth_penalty = 1.0 - carbon_competition_effect
```

**Real C_buf Consumption:**
```
max_consumption = C_buf × 0.10  (10% per timestep)
aox_carbon_consumption = min(aox_carbon_demand, max_consumption)
dC_buf/dt = dC_buf/dt - aox_carbon_consumption
```

### 9. AOX Water Inhibition

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| water_aox_threshold | 0.055 | - | DW/FW threshold for inhibition |
| water_aox_K | 0.020 | - | Water inhibition half-saturation |
| water_aox_max_inhib | 0.50 | - | Maximum water inhibition |
| water_n | 2.0 | - | Hill coefficient for water inhibition |

### 10. AOX Stress Inhibition

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| K_stress_inhib | 150.0 | - | Stress inhibition half-saturation |
| n_stress_inhib | 2.0 | - | Hill coefficient |
| max_stress_inhib | 0.80 | - | Maximum synthesis inhibition |

### 11. AOX Adaptation

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| K_adapt_days | 4.0 | days | Adaptation half-saturation |

### 12. UV-Induced AOX Synthesis

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| k_uv_aox | 7.78e-11 | kg/m²/s | UV induction coefficient |
| K_uv_hours | 30.0 | h | Half-saturation for UV hours effect |

### 13. LAI Efficiency

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| LAI_healthy | 9.0 | m²/m² | Reference LAI for healthy plant |
| n_LAI_eff | 2.0 | - | Hill coefficient for LAI efficiency |

### 14. Night Irradiation Efficiency

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| night_stress_efficiency | 0.4 | - | Efficiency multiplier for night irradiation |

### 15. Nonlinear AOX Efficiency

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| K_nonlin_aox | 800.0 | - | Half-saturation for nonlinear efficiency |
| n_nonlin_aox | 1.5 | - | Hill coefficient for nonlinear efficiency |

### 16. AOX ROS Consumption

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| k_aox_consumption | 1.8e-7 | - | AOX-ROS consumption rate |
| K_ros_consumption | 500.0 | - | ROS half-saturation for consumption |
| n_ros_consumption | 2.0 | - | Hill coefficient for ROS consumption |
| cons_amp_center | 200.0 | - | Softplus center for consumption amplification |
| cons_amp_scale | 15.0 | - | Softplus scale |
| cons_amp_k | 12.0 | - | Amplification coefficient |
| cons_amp_K | 20.0 | - | Amplification half-saturation |

### 17. LDMC (Leaf Dry Matter Content)

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| dw_fw_ratio_base | 0.05 | - | Base DW:FW ratio |
| ldmc_stress_sensitivity | 0.45 | - | LDMC stress sensitivity |
| K_ldmc | 1400.0 | - | LDMC half-saturation |
| dw_fw_ratio_max | 0.080 | - | Maximum DW:FW ratio |
| acute_center | 50.0 | - | Softplus center for acute LDMC |
| acute_scale | 10.0 | - | Softplus scale |
| acute_k | 9.0 | - | Acute LDMC coefficient |
| acute_K | 120.0 | - | Acute LDMC half-saturation |
| acute_n | 2.0 | - | Acute LDMC Hill coefficient |

### 18. Anthocyanin Fraction

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| anthocyanin_fraction | 0.18 | - | Anthocyanin = AOX × 18% |

---

## Key Equations

### Gompertz Nonlinear Damage Factor
```
factor = 1 + 250 × exp(-exp(-0.5 × (hours - 10.5)))
```

**Note:** The nonlinear damage factor is calculated using `hours_today` (current exposure progress)
for progressive damage accumulation during the day. The table below shows FINAL daily values:

| Daily Hours | Factor |
|-------------|--------|
| 3h | 1.0 |
| 6h | 1.0 |
| 9h | 31.1 |
| 12h | 156.9 |

### AOX Synthesis
```
stress_induced = V_max × Stress / (K + Stress) × night_eff × LAI_eff
uv_induced = k_uv_aox × total_hours / (K_uv_hours + total_hours)
aox_synthesis = LAI × (base + uv_induced + stress_induced × adaptation × nonlin_eff)
                × stress_eff × water_eff
```

### Total UVA Hours Calculation
```
# Uses integer day counting for accurate irradiation tracking
day_int = int(day_from_sowing)
days_irradiated = min(day_int - uva_start_day + 1, uva_end_day - uva_start_day + 1)
total_uva_hours = max(0, days_irradiated - 1) × daily_hours + hours_today
```

### Carbon Competition Growth Penalty
```
growth_penalty = 1.0 - (aox_carbon_effect × carbon_competition_max + stress_competition_max × Stress/(stress_competition_K + Stress))
dX_d/dt = dX_d/dt × growth_penalty
```

---

## Simulation Results

### Training Dataset (6 treatments)

| Treatment | Target FW | Sim FW | FW Error | Target Anth | Sim Anth | Anth Error | LAI | Stress |
|-----------|-----------|--------|----------|-------------|----------|------------|-----|--------|
| CK | 87.0 | 87.9 | +1.0% PASS | 433 | 440 | +1.6% PASS | 9.3 | 0.0 |
| L6D6 | 91.4 | 90.2 | -1.3% PASS | 494 | 487 | -1.4% PASS | 9.7 | 7.2 |
| L6D6-N | 80.8 | 83.5 | +3.4% PASS | 493 | 481 | -2.5% PASS | 9.1 | 10.5 |
| VL3D12 | 67.0 | 70.7 | +5.6% **FAIL** | 482 | 509 | +5.6% **FAIL** | 8.3 | 61.0 |
| L6D12 | 60.4 | 61.5 | +1.8% PASS | 518 | 533 | +2.9% PASS | 7.6 | 145.5 |
| H12D3 | 60.6 | 60.9 | +0.6% PASS | 651 | 655 | +0.6% PASS | 8.9 | 262.0 |

**Training Summary:** FW 5/6, Anth 5/6, **Total 10/12**

### Validation Dataset (6 treatments, Day 32-35)

| Treatment | Hours/Day | Target FW | Sim FW | FW Error | Target Anth | Sim Anth | Anth Error |
|-----------|-----------|-----------|--------|----------|-------------|----------|------------|
| CK | 0 | 85.2 | 87.9 | +3.2% PASS | 413 | 440 | +6.5% WARN |
| VL3D3 | 3 | 89.1 | 88.5 | -0.7% PASS | 437 | 460 | +5.4% WARN |
| L6D3 | 6 | 92.1 | 89.1 | -3.4% PASS | 468 | 480 | +2.5% PASS |
| M9D3 | 9 | 83.8 | 85.9 | +2.5% PASS | 539 | 596 | +10.6% FAIL |
| H12D3 | 12 | 62.3 | 60.9 | -2.0% PASS | 657 | 655 | -0.3% PASS |
| VH15D3 | 15 | 51.2 | 51.2 | +0.0% PASS | 578 | 605 | +4.8% PASS |

**Validation Summary:** FW <5%: 6/6, Anth <5%: 3/6, <10%: 5/6

---

## Model Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                    Environment Inputs                           │
│  I_PAR, T, CO2, UVA intensity, UVA schedule                    │
└─────────────────────┬───────────────────────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────────────────────┐
│                    Sun Model Base                               │
│  dX_d/dt, dC_buf/dt, dLAI/dt (from lettuce_uva_carbon_complete) │
└─────────────────────┬───────────────────────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────────────────────┐
│                    UVA Effects                                  │
│  ┌─────────────────────┐    ┌─────────────────────┐            │
│  │ Morphological Boost │    │ ROS Production      │            │
│  │ (SLA, LAI boost)    │    │ dROS/dt = k×I_UVA   │            │
│  └─────────────────────┘    └─────────────────────┘            │
└─────────────────────┬───────────────────────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────────────────────┐
│                    Stress Dynamics                              │
│  Vulnerability = A × exp(-k × LAI) + 1                         │
│  Damage = stress_coeff × ROS × vulnerability × (1 - AOX_prot)  │
│  dStress/dt = Damage - Decay                                   │
└─────────────────────┬───────────────────────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────────────────────┐
│                    AOX Dynamics (v2.0)                          │
│  synthesis = LAI × (base + UV + stress_induced) × efficiencies │
│  dAOX/dt = synthesis - degradation - ROS_consumption           │
└─────────────────────┬───────────────────────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────────────────────┐
│              Carbon Competition (v2.0 NEW)                      │
│  growth_penalty = 1 - carbon_competition_effect                │
│  dX_d/dt = dX_d/dt × growth_penalty                            │
│  dC_buf/dt = dC_buf/dt - AOX_carbon_consumption                │
└─────────────────────┬───────────────────────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────────────────────┐
│                    Output Conversion                            │
│  Anthocyanin = AOX × 0.18                                      │
│  FW = X_d / (DW:FW ratio)                                      │
│  Anth_ppm = Anthocyanin / FW × 1e6                             │
└─────────────────────────────────────────────────────────────────┘
```

---

## Changes from v10.39

| Aspect | v10.39 | v2.0 |
|--------|--------|------|
| State variable | Anthocyanin (Anth) | AOX (total antioxidants) |
| Output | Direct Anth | Anth = AOX × 18% |
| C_buf consumption | None | Real (10% per timestep max) |
| Growth penalty | Stress-based only | Stress + Carbon competition |
| Hardcoded constants | Many | All moved to ALL_PARAMS |
| days_irradiated | Float calculation | Integer for accurate counting |
| Unused parameters | k_repair_lai, k_days_accumulation, power function params | Removed |

---

## Notes

1. **VL3D12 is the hardest case** - Low stress (avgS=61) but needs FW reduction comparable to L6D12 (avgS=146). The stress-based penalty naturally affects L6D12 more.

2. **Carbon consumption is now realistic** - Changed from `C_buf / 300` (0.33%) to `C_buf × 0.10` (10%) per timestep.

3. **AOX scaling** - Since Anth = AOX × 0.18, all AOX rates are scaled by 1/0.18 = 5.56× from original Anth rates.

4. **All parameters now in ALL_PARAMS** - Previously hardcoded constants (LAI efficiency, night efficiency, UV induction, ROS consumption, LDMC acute) are now configurable.

5. **days_irradiated uses integer counting** - Ensures `total_uva_hours` only counts actual irradiation time, not fractional days.
