# UVA-Lettuce Model v11.0 Parameter Tables

## Version History

| Version | Date | Key Changes |
|---------|------|-------------|
| v10.39 | 2026-01-14 | Monotonic decreasing efficiency function |
| **v11.0** | **2026-01-20** | **Carbon Allocation Competition + Coupled Photosynthesis (Sun et al. 2025 framework)** |

---

## NEW Parameters in v11.0

### 1. FvCB Photosynthesis Parameters

| Parameter | Symbol | Value | Unit | Description | Reference |
|-----------|--------|-------|------|-------------|-----------|
| V_cmax25 | $V_{c,max}^{25}$ | 120.0 | µmol CO₂/m²/s | Maximum carboxylation rate at 25°C | Farquhar et al. (1980) |
| J_max25 | $J_{max}^{25}$ | 200.0 | µmol e⁻/m²/s | Maximum electron transport rate at 25°C | von Caemmerer (2000) |
| Rd25 | $R_d^{25}$ | 1.5 | µmol CO₂/m²/s | Dark respiration at 25°C | |
| Gamma_star25 | $\Gamma^{*25}$ | 42.75 | µmol/mol | CO₂ compensation point at 25°C | |
| K_c25 | $K_c^{25}$ | 404.0 | µmol/mol | Michaelis constant for CO₂ at 25°C | |
| K_o25 | $K_o^{25}$ | 278.0 | mmol/mol | Michaelis constant for O₂ at 25°C | |
| O_i | $O_i$ | 210.0 | mmol/mol | O₂ partial pressure | |

### 2. Temperature Response Parameters (Arrhenius)

| Parameter | Symbol | Value | Unit | Description |
|-----------|--------|-------|------|-------------|
| E_Vcmax | $E_{V_{c,max}}$ | 65330 | J/mol | Activation energy for Vcmax |
| E_Jmax | $E_{J_{max}}$ | 43900 | J/mol | Activation energy for Jmax |
| E_Rd | $E_{R_d}$ | 46390 | J/mol | Activation energy for Rd |
| E_Kc | $E_{K_c}$ | 79430 | J/mol | Activation energy for Kc |
| E_Ko | $E_{K_o}$ | 36380 | J/mol | Activation energy for Ko |
| E_Gamma | $E_{\Gamma^*}$ | 37830 | J/mol | Activation energy for Γ* |

### 3. Stomatal Conductance Parameters (Medlyn model)

| Parameter | Symbol | Value | Unit | Description | Reference |
|-----------|--------|-------|------|-------------|-----------|
| g0 | $g_0$ | 0.01 | mol H₂O/m²/s | Minimum stomatal conductance | Medlyn et al. (2011) |
| g1 | $g_1$ | 4.0 | - | Stomatal slope parameter | |

### 4. Carbon Allocation Competition Parameters (CORE CHANGE)

| Parameter | Symbol | Value | Unit | Description | Equation |
|-----------|--------|-------|------|-------------|----------|
| K_C_growth | $K_{C,growth}$ | 0.005 | kg C/m² | Carbon affinity for structural growth | Eq. 4 |
| K_C_anth | $K_{C,anth}$ | 0.0001 | kg C/m² | Carbon affinity for anthocyanin synthesis | Eq. 8 |
| mu_max | $\mu_{max}$ | 1.0e-6 | 1/s | Maximum relative growth rate | Eq. 4 |

**Note:** These parameters are UNCALIBRATED and require optimization to match experimental data.

---

## Key Equations (Modified in v11.0)

### Equation 4: Structural Growth Dynamics

$$\frac{dX_d}{dt} = \mu_{max} \cdot X_d \cdot \underbrace{\left( \frac{C_{buf}}{C_{buf} + K_{C,growth}} \right)}_{\text{Carbon Limitation}} \cdot (1 - f_{inhib}(Stress))$$

Where:
- $X_d$ = Structural dry biomass [kg/m²]
- $C_{buf}$ = Carbon buffer pool [kg/m²]
- $K_{C,growth}$ = Carbon affinity for growth (lower = higher affinity)
- $f_{inhib}(Stress)$ = Stress inhibition function

### Equation 8: Anthocyanin Dynamics

$$\frac{dAnth}{dt} = \underbrace{V_{max,anth} \cdot f_{stress}(ROS) \cdot \left( \frac{C_{buf}}{C_{buf} + K_{C,anth}} \right)}_{\text{Synthesis (Carbon-Limited)}} - k_{deg} \cdot Anth$$

Where:
- $K_{C,anth}$ = Carbon affinity for anthocyanin synthesis
- This creates "resource-dependent defense" mechanism

---

## Parameters Inherited from v10

### UVA Morphological Effect

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| uva_sla_enhancement | 5.00 | - | UVA max enhancement on SLA |
| K_uva_sla | 7.5 | W/m² | SLA enhancement half-saturation |
| uva_lai_boost | 1.70 | - | UVA max enhancement on LAI |
| K_uva_lai | 7.5 | W/m² | LAI enhancement half-saturation |

### ROS Dynamics

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| k_ros_production | 0.010 | ROS/(W/m²·s) | ROS production coefficient |
| k_ros_clearance | 5e-4 | 1/s | ROS clearance coefficient |

### Stress Damage

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| stress_damage_coeff | 1.6e-7 | Stress/(ROS·s) | Base damage coefficient |
| A_vulnerability | 8.5e7 | - | LAI vulnerability amplitude |
| k_vulnerability | 2.0 | 1/(m²/m²) | Vulnerability decay coefficient |
| gompertz_threshold | 10.5 | hours | Antioxidant collapse threshold |
| gompertz_max_factor | 250 | - | Maximum damage multiplier |
| gompertz_steepness | 0.5 | - | Collapse rate |

### Anthocyanin Synthesis

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| V_max_anth | 2.0e-8 | kg/(m²·s) | Max stress-induced synthesis (increased for v11) |
| K_stress_anth | 100.0 | - | Stress half-saturation |
| k_deg | 3.02e-6 | 1/s | Degradation rate |

### LDMC

| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| dw_fw_ratio_base | 0.05 | - | Base DW/FW ratio |
| ldmc_stress_sensitivity | 0.45 | - | Stress sensitivity |
| K_ldmc | 1400.0 | - | LDMC half-saturation |
| dw_fw_ratio_max | 0.080 | - | Maximum DW/FW ratio |

---

## Initial State Variables (v11)

| State | Symbol | Initial Value | Unit | Description |
|-------|--------|---------------|------|-------------|
| Dry weight | $X_d$ | 0.0126 | kg/m² | Structural biomass |
| Carbon buffer | $C_{buf}$ | 0.00126 | kg/m² | Non-structural carbon pool |
| LAI | LAI | 4.0 | m²/m² | Leaf area index |
| Anthocyanin | Anth | 6.3e-6 | kg/m² | Anthocyanin content |
| Stress | Stress | 0 | - | Cumulative stress |
| ROS | ROS | 0 | - | Reactive oxygen species |

### Calculation basis:
- Transplant at Day 14, initial FW = 10 g/plant
- Plant density = 36 plants/m²
- DW/FW ratio (base) = 0.05
- $X_d = 10 \times 0.05 / 1000 \times 36 = 0.0126$ kg/m²
- $C_{buf} = X_d \times 0.1 = 0.00126$ kg/m²
- $LAI = (10 / 0.01) \times 0.04 = 4.0$ m²/m²

---

## Expected Model Behavior Changes (v11 vs v10)

### v10 Behavior
- If UV high → ROS high → Anthocyanin keeps increasing regardless of carbon status

### v11 Behavior (More Realistic)

| Scenario | Photosynthesis | C_buf | Growth | Anthocyanin |
|----------|----------------|-------|--------|-------------|
| High light + Moderate UV | Strong | High | ✓ High | ✓ High |
| Extreme UV (stomata close) | ↓ Drops | ↓ Depleted | ↓ Slows | ↓ Throttled |
| Very high stress | ↓↓ Severely reduced | ↓↓ Near zero | ↓↓ Minimal | ↓↓ Limited by carbon |

### Key Insight: Hormesis Explanation
At very high UV doses, anthocyanin doesn't keep increasing because:
1. Stomata close (Medlyn model response to high VPD/stress)
2. Photosynthesis drops ($A_n$ decreases)
3. Carbon buffer depletes ($C_{buf}$ drops)
4. Anthocyanin synthesis is throttled by carbon limitation (Eq. 8)

This mechanistically explains the hormesis curve observed in VH15D3 (15h/day UV) where anthocyanin is lower than H12D3 (12h/day UV).

---

## Computational Considerations

| Aspect | v10 | v11 |
|--------|-----|-----|
| Photosynthesis calculation | Direct $O(1)$ | Iterative solver $O(10-20)$ per step |
| Estimated slowdown | - | 10-50x |
| Parameter sensitivity | - | K_C parameters may be insensitive if C_buf always high |
| Fitting strategy | Direct | Use simplified mode first, then full model |

### Solver Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| max_iterations | 20 | Max iterations for T_leaf, C_i solver |
| tol | 0.01 | Convergence tolerance |
| use_simplified_photosynthesis | False | Fast mode (skip iterative solver) |

---

## References

1. Farquhar, G.D., von Caemmerer, S., Berry, J.A. (1980). A biochemical model of photosynthetic CO2 assimilation in leaves of C3 species. Planta, 149, 78-90.
2. von Caemmerer, S. (2000). Biochemical Models of Leaf Photosynthesis. CSIRO Publishing.
3. Medlyn, B.E., et al. (2011). Reconciling the optimal and empirical approaches to modelling stomatal conductance. Global Change Biology, 17, 2134-2144.
4. Sun, J., et al. (2025). High-fidelity plant growth modeling framework.
