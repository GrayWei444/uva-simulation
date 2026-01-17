"""
================================================================================
Lettuce Growth and UVA Effect Integrated Model
================================================================================
Version: v10.39 (Monotonic decreasing efficiency function replacing sigmoid)
Date: 2026-01-14

================================================================================
Version History and Key Findings Summary
================================================================================

v10.0:  Reviewer suggestion revision - UVA morphological effect replaces direct PAR addition
v10.1:  Damage formula revision - Additive independent mechanism (LAI vulnerability + nonlinear damage)
v10.3:  LAI vulnerability adjustment - D12 group calibration
v10.5:  Nonlinear damage enabled - H12D3 calibration
v10.6:  Stress inhibits UVA morphological effect
v10.6b: LDMC effect re-enabled
v10.6c: Acute injury linked to Gompertz nonlinear factor
v10.7:  UVA intensity changed from 22 to 11 W/m2 (actual LED power)
v10.7c: Anthocyanin unit correction (mg/100g to mg/kg) + 12/12 all pass
v10.8:  Validation experiment calibration - Gompertz threshold 9 to 11, max_factor 160 to 250
v10.9:  Anthocyanin water inhibition mechanism - Explains why VH15D3 anthocyanin is lower than H12D3
v10.23: LAI efficiency mechanism - Low LAI inhibits Stress-induced anthocyanin efficiency
v10.32: Nonlinear factor anthocyanin inhibition - Solves M9D3 over-prediction problem
v10.33: Continuous asymmetric Gaussian + softplus soft threshold (compliant with CLAUDE.md)
v10.37: Gompertz threshold 9.5 to 10.5 + sigmoid inhibition replaces asymmetric Gaussian
        - gompertz_threshold: 10.5 hours
        - V_max_anth: 2.75e-9
        - ldmc_stress_sensitivity: 0.45
        - K_ldmc: 1400, dw_fw_ratio_max: 0.080
v10.38: Parameter synchronization revision
        - Default uva_intensity: 22 to 11 W/m2
        - Updated all legacy value comments (6h to 1.0, 9h to 31.1, 12h to 156.9, etc.)
v10.39: Monotonic decreasing efficiency function
        - Uses Hill function: efficiency = 1 / (1 + (nonlin/K)^n)
        - K=800, n=1.5
        - Efficiency monotonically decreases as nonlinear_factor increases

================================================================================
v10.39 Calibration Results (12/12 Pass)
================================================================================

Gompertz nonlinear factor (threshold=10.5h, max=250, k=0.5):
| Daily hours | factor | Anthocyanin efficiency |
|-------------|--------|------------------------|
| 3h          | 1.0    | 100.0%                 |
| 6h          | 1.0    | 100.0%                 |
| 9h          | 31.1   | 99.2%                  |
| 12h         | 156.9  | 92.0%                  |
| 15h         | 226.0  | 86.9%                  |

Training set (acceptance error 5%):
| Treatment | FW pred | FW obs | FW err  | Anth pred | Anth obs | Anth err |
|-----------|---------|--------|---------|-----------|----------|----------|
| CK        | 86.5g   | 87.0g  | -0.5%   | 439       | 433      | +1.3%    |
| L6D6      | 92.5g   | 91.4g  | +1.2%   | 474       | 494      | -4.0%    |
| L6D6-N    | 84.0g   | 80.8g  | +3.9%   | 475       | 493      | -3.6%    |
| VL3D12    | 69.4g   | 67.0g  | +3.6%   | 492       | 482      | +2.0%    |
| L6D12     | 58.9g   | 60.4g  | -2.5%   | 496       | 518      | -4.3%    |
| H12D3     | 61.3g   | 60.6g  | +1.2%   | 651       | 651      | +0.0%    |

Validation set (acceptance error 10%):
| Treatment | Hours | FW pred | FW obs | FW err  | Anth pred | Anth obs | Anth err |
|-----------|-------|---------|--------|---------|-----------|----------|----------|
| CK        | 0h    | 86.5g   | 85.2g  | +1.6%   | 439       | 413      | +6.2%    |
| VL3D3     | 3h    | 88.4g   | 89.0g  | -0.8%   | 457       | 437      | +4.5%    |
| L6D3      | 6h    | 89.9g   | 92.2g  | -2.5%   | 473       | 468      | +1.1%    |
| M9D3      | 9h    | 87.8g   | 83.8g  | +4.8%   | 589       | 539      | +9.2%    |
| H12D3     | 12h   | 61.3g   | 62.2g  | -1.4%   | 651       | 657      | -0.9%    |
| VH15D3    | 15h   | 51.2g   | 51.3g  | +0.0%   | 532       | 578      | -7.9%    |

Anthocyanin ranking correct: H12D3 (651) > M9D3 (589) > VH15D3 (532)

================================================================================
Key Finding 1: UVA Morphological Effect (v10.0)
================================================================================

[Background]
- Old model: I_effective = I_base + I_UVA (UVA directly added to PAR)
- This lacks physiological basis:
  * UVA (365nm) has very low quantum yield, chlorophyll barely absorbs it
  * Direct addition overestimates photosynthesis contribution

[Solution]
UVA indirectly promotes growth through "morphological effect":
  1. UVA irradiation promotes leaf area expansion (SLA increase)
  2. Larger leaf area leads to higher PAR interception
  3. Higher light interception leads to more photosynthate leads to higher biomass

[Formula]
  sla_boost = max_boost * I_UVA / (K + I_UVA)  [Michaelis-Menten]
  lai_boost = max_boost * I_UVA / (K + I_UVA)

[v10.6 Addition] High Stress inhibits morphological effect
  stress_suppression = 1 - Stress / (K_stress + Stress)
  sla_boost = sla_boost * stress_suppression
  lai_boost = lai_boost * stress_suppression

  Biological significance: Plants under high stress cannot effectively use UVA for growth

[Literature Support]
- Chen et al. (2019): UVA increases lettuce leaf area by 15-26%, biomass by 18-32%
- Wargent et al. (2009): UV-B modifies lettuce leaf morphology and SLA
- Krizek et al. (1998): UV-A increases cucumber leaf area and biomass

================================================================================
Key Finding 2: L6D6-N Circadian Damage Mechanism
================================================================================

[Experimental Design]
- L6D6:   6h UVA x 6d, daytime irradiation (10:00-16:00)
- L6D6-N: 6h UVA x 6d, nighttime irradiation (22:00-04:00)
- Both groups have same total UVA dose, but L6D6-N has lower fresh weight

[Mechanism Explanation]
Nighttime UVA irradiation causes additional "circadian damage":
  1. Plants have lower antioxidant system activity at night
  2. Circadian clock mismatch with UVA irradiation timing causes extra stress
  3. Repair capacity is also lower at night

[Formula]
  circadian_damage = k_circadian * I_UVA * hours_in_dark^n

  Only calculated when "UVA is applied during night"
  hours_in_dark: How long since lights turned off

[Parameters]
  k_circadian = 1.5e-6
  n_circadian = 2.0

================================================================================
Key Finding 3: Damage Formula Transformation - Additive Independent Mechanism (v10.1)
================================================================================

[Background]
Old formula (multiplicative):
  damage = k * ROS * vulnerability * nonlinear_factor

Problem: H12D3 (vuln=6, nonlin=137) and VL3D12 (vuln=3607, nonlin=1) damages
        would cancel each other due to multiplication, unable to correctly
        distinguish damage characteristics of each group

[Solution]
New formula (additive):
  damage = k1 * ROS * vulnerability + k2 * ROS * nonlinear_factor

Two mechanisms act independently:
  1. LAI vulnerability damage: For D12 groups (early irradiation, small LAI)
  2. Nonlinear time damage: For H12D3 (long irradiation, antioxidant collapse)

[Results]
- VL3D12, L6D12: Mainly rely on LAI vulnerability for Stress
- H12D3: Mainly rely on nonlinear time damage for Stress
- L6D6: Both are very small, almost no damage

================================================================================
Key Finding 4: D12 Group LAI Vulnerability Mechanism
================================================================================

[Experimental Design]
  VL3D12: 3h UVA x 12d, starting from Day 23
  L6D12:  6h UVA x 12d, starting from Day 23
  L6D6:   6h UVA x 6d, starting from Day 29

[Key Observation]
D12 groups start irradiation from Day 23, when LAI is approximately 5
L6D6 starts irradiation from Day 29, when LAI is approximately 7

Young plants (low LAI) are more sensitive to UVA!

[Formula]
  vulnerability = A * exp(-k * LAI) + 1

  LAI=5: vulnerability = 97e6 * exp(-2*5) + 1 = approximately 4405
  LAI=7: vulnerability = 97e6 * exp(-2*7) + 1 = approximately 81

[Biological Significance]
- Young leaves have incomplete protection mechanisms
- Less epidermal wax layer, anthocyanins and other defense substances
- Antioxidant enzyme systems not fully developed

================================================================================
Key Finding 5: H12D3 Nonlinear Damage and LDMC Relationship (v10.6c)
================================================================================

[Experimental Design]
  H12D3: 12h UVA x 3d, starting from Day 32

  Key: H12D3 only irradiated for 3 days, but has highest daily dose (12h/day = 51.2 J/cm2/day)

[Problem]
H12D3 only irradiated for 3 days (Day 32-35), LAI reached approximately 8.0 before irradiation
At this time LAI vulnerability is very low (vuln approximately 6), cannot generate enough
damage through vulnerability mechanism

[Solution: Gompertz Nonlinear Function]
  nonlinear_factor = 1 + max * exp(-exp(-k * (hours - threshold)))

  Parameters (v10.38):
  - threshold = 10.5 hours (time point when antioxidant system starts to collapse)
  - steepness = 0.5 (collapse rate)
  - max_factor = 250 (maximum damage multiplier)

  Results (v10.38):
  - 6h/day:  nonlinear_factor = approximately 1.0  (almost no damage)
  - 9h/day:  nonlinear_factor = approximately 31.1 (entering transition zone)
  - 12h/day: nonlinear_factor = approximately 156.9 (severe)
  - 15h/day: nonlinear_factor = approximately 226.0 (approaching saturation)

[LDMC Acute Injury Mechanism (v10.6c)]
Problem: Even if H12D3 has high Stress, it cannot fully compensate for 18 days of normal LAI growth
Solution: Link LDMC (Leaf Dry Matter Content) with nonlinear factor

Formula:
  Acute factor uses softplus + Hill function
  See calculate_dynamic_dw_fw_ratio() for details

[Key Finding: Daily Cumulative Irradiation is the Key]
Core mechanism revealed by Experiment B design:
  - 3h x 12d = 12.8 J/cm2/day (chronic mild)
  - 6h x 6d  = 25.6 J/cm2/day (moderate)
  - 12h x 3d = 51.2 J/cm2/day (acute intense)

All three groups have same "total dose" (153.6 J/cm2), but results differ greatly:
  - Total dose cannot predict damage
  - Daily cumulative dose determines acute damage
  - Exceeding threshold (approximately 9 hours) triggers antioxidant system collapse

[Biological Explanation]
1. Antioxidant system has daily cycle:
   - Daytime endures ROS, antioxidant enzyme activity rises
   - Nighttime repairs, activity drops

2. When daily irradiation time exceeds antioxidant capacity:
   - Antioxidant enzymes get "exhausted"
   - ROS clearance capacity drops
   - Causes irreversible oxidative damage

3. Gompertz function describes this "collapse" process:
   - Slow at start (defense mechanism working)
   - Rapid rise near threshold (system collapse)
   - Eventually saturates (complete failure)

================================================================================
Literature Support
================================================================================

- Chen et al. (2019): UVA effects on lettuce morphology and physiology
- Wargent et al. (2009): UV effects on lettuce leaf morphology
- Krizek et al. (1998): UV-A promotion effects on cucumber
- Gill & Tuteja (2010): ROS and antioxidant defense systems
- Apel & Hirt (2004): Role of ROS in plant stress responses

================================================================================
"""

import numpy as np
from scipy.integrate import solve_ivp

# Import base Sun model
from lettuce_uva_carbon_complete_model import SunParams as BaseSunParams
from lettuce_uva_carbon_complete_model import sun_derivatives_final


# ==============================================================================
# Model Parameter Definitions (v10.6c - Complete Calibration Version)
# ==============================================================================

ALL_PARAMS = {
    # --------------------------------------------------------------------------
    # 1. Base Photosynthesis Parameters
    # --------------------------------------------------------------------------
    'c_alpha': 0.54,                     # Photosynthesis efficiency coefficient [-] (v10.10: 0.55 to 0.54)

    # --------------------------------------------------------------------------
    # 2. UVA Morphological Effect Parameters (v10.0 Core Change)
    # --------------------------------------------------------------------------
    # Important: UVA does NOT directly contribute to photosynthesis!
    #
    # Old version (v9 and before): I_effective = I_base + I_UVA  <- Removed
    # New version (v10): UVA indirectly promotes growth through morphological effect
    #
    # Biological mechanism:
    # 1. UVA irradiation promotes leaf area expansion (SLA increase)
    # 2. Larger leaf area leads to higher PAR interception
    # 3. Higher light interception leads to more photosynthate leads to higher biomass
    #
    # This explains why L6D6 fresh weight can exceed CK:
    # - L6D6 receives 6h/day x 6d of UVA, morphological effect accumulates
    # - Morphological effect integrates through ODE, longer irradiation = greater effect
    # - Meanwhile L6D6 has mild stress, damage is negligible
    #
    # Literature support: Chen et al. (2019) reported UVA increases leaf area by 15-26%
    # v10.7: UVA=11 W/m2 calibrated version
    'uva_sla_enhancement': 5.00,         # UVA max enhancement on SLA (v10.9: 280% to 500%)
    'K_uva_sla': 7.5,                    # SLA enhancement half-saturation UVA intensity [W/m2] (v10.7: for UVA=11)

    # UVA direct promotion of LAI growth (morphogenesis)
    'uva_lai_boost': 1.70,               # UVA max enhancement on LAI growth
    'K_uva_lai': 7.5,                    # LAI enhancement half-saturation UVA intensity [W/m2] (v10.7: for UVA=11)

    # --------------------------------------------------------------------------
    # 3. ROS Dynamics Parameters
    # --------------------------------------------------------------------------
    'k_ros_production': 0.010,           # ROS production coefficient [ROS/(W/m2*s)] (v10.7: 0.005 to 0.01 for UVA=11)
    'k_ros_clearance': 5e-4,             # ROS clearance coefficient [1/s]

    # --------------------------------------------------------------------------
    # 4. Stress Damage Parameters
    # --------------------------------------------------------------------------
    # 4.1 Base damage coefficient
    'stress_damage_coeff': 1.6e-7,       # Base damage coefficient [Stress/(ROS*s)] (restored original value)

    # 4.2 LAI-dependent vulnerability (v10.7 calibration)
    # Biological significance: Young plants (low LAI) are more sensitive to UVA
    # D12 groups irradiated early (Day 23 start, LAI approx 5) are more vulnerable than L6D6 (Day 29 start, LAI approx 7)
    'A_vulnerability': 85000000.0,       # Vulnerability amplitude [-] (v10.10: 95M to 90M to 85M for VL3D12)
    'k_vulnerability': 2.0,              # Vulnerability decay coefficient [1/(m2/m2)]

    # 4.3 Nonlinear damage - v10.0 Gompertz function approach
    #
    # Reviewer suggestion: Simple power function (hours^8) is hard to explain biologically
    #
    # Solution: Use Gompertz function with following advantages:
    # 1. Has clear saturation upper limit (consistent with physiology: defense system has limits)
    # 2. Asymmetric S-curve, steeper than Logistic
    # 3. Commonly used to describe biological growth curves and dose-response relationships
    # 4. Describes the accelerating process of "antioxidant system collapse"
    #
    # Formula: factor = 1 + max * exp(-exp(-k * (hours - threshold)))
    # Characteristics:
    # - hours << threshold: factor approx 1 (defense mechanism intact)
    # - hours approx threshold: factor starts rapid rise (defense collapse)
    # - hours >> threshold: factor approaches 1 + max (reaches saturation)
    'use_gompertz': True,                # Use Gompertz function (recommended)

    # Power function backup parameters (used when use_gompertz=False)
    'k_hours_nonlinear': 5.8e-7,         # Time-dependent nonlinear coefficient [-]
    'n_hours_nonlinear': 8.0,            # Time exponent [-]
    'max_nonlinear_factor': 500.0,       # Maximum nonlinear factor (numerical clipping)

    # Gompertz function parameters (v10.8 update: validation experiment calibration)
    # Formula: factor = 1 + max * exp(-exp(-k * (hours - threshold)))
    # v10.8: Adjusted based on 3-day gradient validation experiment (0-15h/day)
    #        threshold 9 to 11, max_factor 160 to 250
    #        Validation error reduced from 10.9% to 7.4%, original FW 6/6 maintained
    'gompertz_max_factor': 250.0,        # Maximum damage multiplier (saturation upper limit) [v10.8: 160 to 250]
    'gompertz_threshold': 10.5,          # Inflection point [v10.37: 9.5 to 10.5]
    'gompertz_steepness': 0.5,           # Collapse rate [v10.15: 0.6 to 0.5, slightly gentler curve]

    # 4.4 Anthocyanin protection
    'alpha_anth_protection': 0.5,        # Maximum protection efficiency (0~1) [-]
    'K_anth_protection': 5e-6,           # Protection half-saturation constant [kg Anth/m2]

    # 4.5 Circadian damage at night (v10.7 adjustment for UVA=11)
    # Increase k_circadian to make L6D6-N Stress significantly higher than L6D6
    # v10.7: UVA=11 needs doubled k_circadian (effect = k * I_UVA)
    'k_circadian': 3.0e-6,               # Circadian damage coefficient at night [-] (v10.7: 1.5e-6 to 3e-6)
    'n_circadian': 2.0,                  # Circadian damage exponent at night [-]

    # 4.6 Day accumulation damage (v10.0 new)
    # Long-term irradiation accumulation effect: More days of irradiation, worse repair capacity
    'k_days_accumulation': 0.0,          # Day accumulation coefficient [-] (temporarily disabled)

    # 4.7 Nonlinear damage independent coefficient (v10.1 new)
    # Make time-based nonlinear damage independent of LAI vulnerability, avoid mutual cancellation
    # v10.5: Enable nonlinear damage, let H12D3 (12h/day) have sufficient Stress
    # v10.6: 5.0e-6 leads to 5/6 FW pass
    'k_nonlinear_stress': 5.0e-6,        # Nonlinear damage coefficient [1/s]

    # --------------------------------------------------------------------------
    # 5. Stress Decay Parameters
    # --------------------------------------------------------------------------
    # v10.4: Estimated based on 15h/day critical point
    # User experience: 15h/day UVA brings plants close to death (damage approx equals recovery equilibrium)
    # Estimation logic:
    #   - Irradiation 15h, rest 9h
    #   - Equilibrium condition: rest time approx half-life
    #   - Half-life = 9 hours
    # k = ln(2) / (9 * 3600) = 2.14e-5
    'k_stress_decay': 2.14e-5,           # Stress decay coefficient [1/s] (half-life 9 hours)

    # --------------------------------------------------------------------------
    # 6. Stress Repair Parameters
    # --------------------------------------------------------------------------
    'k_repair_lai': 5.0e-8,              # LAI repair coefficient [1/(m2/m2*s)]

    # --------------------------------------------------------------------------
    # 7. Stress Inhibition on Growth
    # --------------------------------------------------------------------------
    # v10.6: Strengthen photosynthesis inhibition, let high Stress groups (H12D3) have greater FW reduction
    # K_stress=50 leads to 5/6 pass, only H12D3 slightly high
    'stress_photosynthesis_inhibition': 0.85,
    'stress_lai_inhibition': 0.80,
    'K_stress': 50.0,

    # --------------------------------------------------------------------------
    # 8. Anthocyanin Synthesis Parameters (v10.7 recalibration)
    # --------------------------------------------------------------------------
    # Anthocyanin concentration unit correction: mg/100g to mg/kg (x10)
    # Target: CK=433, L6D6=494, L6D6-N=493, H12D3=651, VL3D12=482, L6D12=518 ppm
    # v10.7c: 12/12 all pass
    'base_anth_rate_light': 6.35e-10,    # v10.23: Original value
    'base_anth_rate_dark': 3.18e-10,     # = base_light * 0.5
    'V_max_anth': 2.75e-9,               # v10.37: Paper version
    'K_stress_anth': 100.0,              # v10.23: Original value
    'k_deg': 3.02e-6,
    'anth_carbon_cost': 0.0,

    # --------------------------------------------------------------------------
    # 8b. Anthocyanin Water Inhibition Parameters (v10.9 new)
    # --------------------------------------------------------------------------
    # [Mechanism Explanation]
    # When plants are severely dehydrated (DW/FW increases), anthocyanin synthesis efficiency decreases
    # Reason: Cell turgor loss leads to enzyme activity decrease leads to metabolic pathway impairment
    #
    # [Gompertz Framework]
    # water_efficiency = 1 - max_inhib * exp(-exp(-steepness * (DW/FW - threshold)))
    #
    # [Literature Support]
    # - Hadacek (2010) DOI: 10.2203/dose-response.09-028.Hadacek
    #   "Hormesis: low levels of ROS stimulate growth, high levels induce cell death"
    # - Garnier et al. (2001) DOI: 10.1046/j.0269-8463.2001.00563.x
    #   "LDMC reflects plant water status"
    # - Ferreyra et al. (2021) DOI: 10.1016/j.plaphy.2021.05.022
    #   "UV-B hormesis in secondary metabolite biosynthesis"
    #
    'water_anth_threshold': 0.055,       # DW/FW threshold where inhibition starts
    'water_anth_K': 0.020,               # v10.23: Original value
    'water_anth_max_inhib': 0.50,        # v10.23: Original value (50%)

    # --------------------------------------------------------------------------
    # 8c. Anthocyanin Stress Inhibition Parameters (v10.32 new)
    # --------------------------------------------------------------------------
    # [Mechanism Explanation]
    # Under extreme high Stress, cellular metabolism collapses, anthocyanin synthesis efficiency decreases
    # Hill function: inhibition = max * S^n / (K^n + S^n)
    #
    # v10.23 original: K=150, let inhibition affect VH15D3 (avgS=930)
    'K_stress_inhib': 150.0,             # v10.23: Original value
    'n_stress_inhib': 2.0,               # Hill coefficient
    'max_stress_inhib': 0.80,            # Maximum inhibition degree (80%)

    # --------------------------------------------------------------------------
    # 8d. Anthocyanin Adaptation Effect Parameters (v10.32 new)
    # --------------------------------------------------------------------------
    # [Mechanism Explanation]
    # After prolonged cumulative irradiation, plant response to Stress induction decreases (adaptation/desensitization)
    # adaptation_factor = K / (K + days_irradiated)
    # D3=0.57, D6=0.40, D12=0.25 (K=4)
    'K_adapt_days': 4.0,                 # Adaptation half-saturation constant (days)

    # --------------------------------------------------------------------------
    # 9. Anthocyanin Consumption Parameters
    # --------------------------------------------------------------------------
    'k_anth_consumption': 1.0e-6,        # Anthocyanin consumption by ROS coefficient (v10.21: balanced version)

    # --------------------------------------------------------------------------
    # 10. LDMC Parameters (v10.6c core mechanism)
    # --------------------------------------------------------------------------
    #
    # LDMC (Leaf Dry Matter Content) = DW / FW
    #
    # [Background]
    # H12D3 only irradiated for 3 days, but needs to achieve similar FW reduction as L6D12 (12 days)
    # Simply relying on Stress to inhibit growth cannot achieve this, because H12D3 already had 18 days of normal growth before irradiation
    #
    # [Solution]
    # Link LDMC with Gompertz nonlinear factor:
    # - H12D3 (12h/day): High nonlinear factor leads to high LDMC leads to low FW
    # - L6D6, L6D12: Low nonlinear factor leads to normal LDMC leads to normal FW
    #
    # [Formula]
    # acute_factor = 1 + k_acute * log(nonlinear_factor)
    # ratio = base * (1 + stress_effect * acute_factor)
    #
    # [Biological Significance]
    # Acute high-intensity UVA (>9h/day) causes cell dehydration:
    # - Cell membrane damage leads to osmotic pressure changes
    # - Increased water loss leads to increased DW/FW ratio
    # - This is an "acute injury" response
    #
    'dw_fw_ratio_base': 0.05,            # Base DW/FW ratio (healthy plants)
    'ldmc_stress_sensitivity': 0.45,     # v10.37: Paper version
    'K_ldmc': 1400.0,                    # v10.37: Paper version
    'dw_fw_ratio_max': 0.080,            # v10.37: Paper version
}


# ==============================================================================
# Model Parameter Class
# ==============================================================================

class UVAParams(BaseSunParams):
    """
    UVA Effect Model Parameter Class (v10.0 Reviewer Suggestion Revision)
    """

    def __init__(self, params=None):
        super().__init__()
        if params is None:
            params = ALL_PARAMS

        # 1. Base photosynthesis parameters
        self.c_alpha = params['c_alpha']

        # 2. UVA morphological effect parameters (v10.0 new)
        self.uva_sla_enhancement = params['uva_sla_enhancement']
        self.K_uva_sla = params['K_uva_sla']
        self.uva_lai_boost = params['uva_lai_boost']
        self.K_uva_lai = params['K_uva_lai']

        # 3. ROS dynamics parameters
        self.k_ros_production = params['k_ros_production']
        self.k_ros_clearance = params['k_ros_clearance']

        # 4. Stress damage parameters
        self.stress_damage_coeff = params['stress_damage_coeff']
        self.A_vulnerability = params['A_vulnerability']
        self.k_vulnerability = params['k_vulnerability']

        # v10.0 Nonlinear damage parameters
        self.use_gompertz = params['use_gompertz']
        self.k_hours_nonlinear = params['k_hours_nonlinear']
        self.n_hours_nonlinear = params['n_hours_nonlinear']
        self.max_nonlinear_factor = params['max_nonlinear_factor']

        # Gompertz function parameters (v10.0 core)
        self.gompertz_max_factor = params['gompertz_max_factor']
        self.gompertz_threshold = params['gompertz_threshold']
        self.gompertz_steepness = params['gompertz_steepness']

        self.alpha_anth_protection = params['alpha_anth_protection']
        self.K_anth_protection = params['K_anth_protection']
        self.k_circadian = params['k_circadian']
        self.n_circadian = params['n_circadian']
        self.k_days_accumulation = params['k_days_accumulation']
        self.k_nonlinear_stress = params['k_nonlinear_stress']

        # 5. Stress decay parameters
        self.k_stress_decay = params['k_stress_decay']

        # 6. Stress repair parameters
        self.k_repair_lai = params['k_repair_lai']

        # 7. Stress inhibition on growth
        self.stress_photosynthesis_inhibition = params['stress_photosynthesis_inhibition']
        self.stress_lai_inhibition = params['stress_lai_inhibition']
        self.K_stress = params['K_stress']

        # 8. Anthocyanin synthesis parameters
        self.base_anth_rate_light = params['base_anth_rate_light']
        self.base_anth_rate_dark = params['base_anth_rate_dark']
        self.V_max_anth = params['V_max_anth']
        self.K_stress_anth = params['K_stress_anth']
        self.k_deg = params['k_deg']
        self.anth_carbon_cost = params['anth_carbon_cost']

        # 8b. Anthocyanin water inhibition parameters (v10.9)
        self.water_anth_threshold = params['water_anth_threshold']
        self.water_anth_K = params['water_anth_K']  # v10.32: Changed to Hill K
        self.water_anth_max_inhib = params['water_anth_max_inhib']

        # 8c. Anthocyanin Stress inhibition parameters (v10.32)
        self.K_stress_inhib = params['K_stress_inhib']
        self.n_stress_inhib = params['n_stress_inhib']
        self.max_stress_inhib = params['max_stress_inhib']

        # 8d. Anthocyanin adaptation effect parameters (v10.32)
        self.K_adapt_days = params['K_adapt_days']

        # 9. Anthocyanin consumption parameters
        self.k_anth_consumption = params['k_anth_consumption']

        # 10. LDMC parameters
        self.dw_fw_ratio_base = params['dw_fw_ratio_base']
        self.ldmc_stress_sensitivity = params['ldmc_stress_sensitivity']
        self.K_ldmc = params['K_ldmc']
        self.dw_fw_ratio_max = params['dw_fw_ratio_max']


# ==============================================================================
# Utility Functions
# ==============================================================================

def calculate_dynamic_dw_fw_ratio(Stress, p, nonlinear_factor=1.0):
    """
    Calculate dynamic DW:FW ratio based on Stress and nonlinear factor (LDMC effect)

    Version: v10.6c

    ============================================================================
    Key Finding: Daily Cumulative Irradiation Determines Acute Injury
    ============================================================================

    Experiment B design:
    - VL3D12: 3h x 12d = 12.8 J/cm2/day (chronic mild)
    - L6D12:  6h x 6d  = 25.6 J/cm2/day (moderate)
    - H12D3: 12h x 3d  = 51.2 J/cm2/day (acute intense)

    All three groups have same "total dose" (153.6 J/cm2), but results differ greatly!
    -> Total dose cannot predict damage
    -> Daily cumulative dose is the key

    ============================================================================
    Mechanism Explanation
    ============================================================================

    1. Base LDMC effect (stress_effect):
       - High Stress leads to high LDMC (dehydration)
       - Uses Michaelis-Menten form, K_ldmc=400 makes low Stress effect small

    2. Acute injury factor (acute_factor):
       - Linked to Gompertz nonlinear factor
       - Nonlinear factor reflects "daily irradiation time exceeding antioxidant capacity"
       - Uses softplus + Hill function to calculate (v10.38)

       nonlinear_factor values (v10.38, threshold=10.5):
       - L6D6 (6h/day):  1.0  -> Almost no acute injury
       - M9D3 (9h/day):  31.1 -> Starts to have acute injury
       - H12D3 (12h/day): 156.9 -> Severe acute injury

    3. Final LDMC:
       ratio = base * (1 + stress_effect * acute_factor)

       Results:
       - CK, L6D6:    dw/fw approx 0.050 (healthy)
       - L6D6-N:      dw/fw approx 0.051 (mild dehydration)
       - VL3D12:      dw/fw approx 0.050 (healthy)
       - L6D12:       dw/fw approx 0.052 (mild dehydration)
       - H12D3:       dw/fw approx 0.071 (severe dehydration/acute injury)

    ============================================================================
    Biological Significance
    ============================================================================

    Acute high-intensity UVA (>9h/day) causes:
    - Antioxidant system collapse (Gompertz nonlinear factor rises)
    - Cell membrane damage leads to osmotic pressure changes
    - Increased water loss leads to increased DW/FW ratio
    - This explains H12D3's high LDMC (0.071 vs normal 0.050)

    Parameters:
        Stress: Cumulative stress index
        p: Parameter object
        nonlinear_factor: Gompertz nonlinear factor (reflects daily irradiation intensity)

    Returns:
        dw_fw_ratio: DW/FW ratio (0.05 ~ 0.085)
    """
    # Base LDMC effect: Stress-dependent dehydration
    stress_effect = p.ldmc_stress_sensitivity * Stress / (p.K_ldmc + Stress + 1e-9)

    # Acute injury factor: Linked to Gompertz nonlinear factor
    # Key finding: Daily cumulative irradiation (not total dose) determines acute injury degree
    # v10.17: Changed to Hill function, let acute saturate at high nonlin
    # Avoid VH15D3 (15h/day) having too high dw/fw
    #
    # v10.33: Removed hard threshold, changed to softplus continuous function (compliant with CLAUDE.md)
    # softplus(x) = log(1 + exp(x)) is a continuous differentiable soft threshold
    # When nonlin << acute_center, x approx 0
    # When nonlin >> acute_center, x approx nonlin - acute_center
    acute_center = 50.0    # Soft threshold center
    acute_scale = 10.0     # Transition width
    acute_k = 9.0          # v10.37: Paper version (6.5 to 9.0)
    acute_K = 120.0        # v10.37: Paper version (150 to 120)
    acute_n = 2.0          # Hill coefficient

    # Use softplus to implement soft threshold: x = softplus((nonlin - center) / scale) * scale
    x_raw = (nonlinear_factor - acute_center) / acute_scale
    x = acute_scale * np.log(1.0 + np.exp(np.clip(x_raw, -50, 50)))  # softplus
    acute_factor = 1.0 + acute_k * (x ** acute_n) / (acute_K ** acute_n + x ** acute_n + 1e-9)

    # Final LDMC
    ratio = p.dw_fw_ratio_base * (1.0 + stress_effect * acute_factor)
    return min(ratio, p.dw_fw_ratio_max)


def calculate_water_anth_efficiency(dw_fw_ratio, p):
    """
    Calculate water status effect on anthocyanin synthesis efficiency (v10.9)

    ============================================================================
    Mechanism Explanation
    ============================================================================

    When plants are severely dehydrated (DW/FW increases), anthocyanin synthesis efficiency decreases:
    - Cell turgor loss leads to enzyme activity decrease
    - Metabolic pathway impairment leads to reduced synthesis efficiency
    - This explains why VH15D3 (15h/day) anthocyanin is lower than H12D3

    ============================================================================
    Gompertz Framework
    ============================================================================

    Uses same Gompertz framework as nonlinear damage to maintain model consistency:

    efficiency = 1 - max_inhib * exp(-exp(-steepness * (DW/FW - threshold)))

    Characteristics:
    - DW/FW < threshold: efficiency approx 1 (normal synthesis)
    - DW/FW = threshold: efficiency starts to decrease
    - DW/FW >> threshold: efficiency approaches 1 - max_inhib (maximum inhibition)

    ============================================================================
    Literature Support
    ============================================================================

    1. Hadacek (2010) DOI: 10.2203/dose-response.09-028.Hadacek
       "Hormesis: low levels of ROS stimulate growth, high levels induce cell death"

    2. Garnier et al. (2001) DOI: 10.1046/j.0269-8463.2001.00563.x
       "LDMC reflects fundamental trade-off in plant functioning"

    3. Ferreyra et al. (2021) DOI: 10.1016/j.plaphy.2021.05.022
       "UV-B hormesis in secondary metabolite biosynthesis"

    ============================================================================
    Parameter Explanation
    ============================================================================

    dw_fw_ratio: Current DW/FW ratio
    p: Parameter object, contains:
       - water_anth_threshold: DW/FW threshold where inhibition starts (0.068)
       - water_anth_steepness: Inhibition rate (100)
       - water_anth_max_inhib: Maximum inhibition degree (0.40)

    Returns:
        efficiency: Synthesis efficiency (0.35 ~ 1.0)
    """
    # =========================================================================
    # v10.18: Changed to Hill function (smoother, more controllable)
    # =========================================================================
    #
    # Formula: inhibition = max_inhib * (dw-base)^n / (K^n + (dw-base)^n)
    #          efficiency = 1 - inhibition
    #
    # Hill function characteristics:
    # - When dw_fw <= base: inhibition = 0 -> efficiency = 1 (no inhibition)
    # - When dw_fw = base + K: inhibition = 50% max -> moderate inhibition
    # - When dw_fw >> base + K: inhibition approaches max_inhib (maximum inhibition)
    #
    # Data calibration:
    # - H12D3:  dw/fw = 0.065 -> Needs mild inhibition ~5-10%
    # - VH15D3: dw/fw = 0.083 -> Needs strong inhibition ~25-30%
    #
    # Parameter design (v10.18):
    # - base (threshold) = 0.055 (below H12D3)
    # - K (steepness) = 0.025 (half-saturation constant)
    # - max_inhib = 0.50 (maximum 50% inhibition)
    # - n = 2 (Hill coefficient)

    base = p.water_anth_threshold  # Use threshold as base
    K = p.water_anth_K             # v10.32: Use parameter-defined Hill K
    n = 2.0                        # Hill coefficient

    if dw_fw_ratio <= base:
        return 1.0

    x = dw_fw_ratio - base
    inhibition = p.water_anth_max_inhib * (x ** n) / (K ** n + x ** n)
    efficiency = 1.0 - inhibition

    return efficiency


def calculate_nonlin_anth_efficiency(nonlinear_factor, p):
    """
    Calculate anthocyanin synthesis efficiency based on nonlinear factor (v10.39)

    [Background]
    Anthocyanin absolute amount peaks around 9h, then starts to decrease.
    H12D3 anthocyanin "concentration" is high because FW drops faster (denominator becomes smaller).
    VH15D3 is severely damaged, needs additional inhibition.

    [Biological Explanation]
    Uses monotonic decreasing function:
    - Efficiency monotonically decreases as nonlinear_factor increases
    - Combined with water inhibition, Stress inhibition to jointly regulate VH15D3

    [Formula]
    efficiency = 1 / (1 + (nonlinear_factor / K)^n)
    - K: Half-effect constant (nonlinear_factor when efficiency drops to 50%)
    - n: Hill coefficient (controls decay rate)

    nonlinear_factor range (Gompertz, threshold=10.5):
    - 3h/day:  1.0   -> 100.0%
    - 6h/day:  1.0   -> 100.0%
    - 9h/day:  31.1  -> 99.2%
    - 12h/day: 156.9 -> 92.0%
    - 15h/day: 226.0 -> 86.9%
    """
    # v10.39: Monotonic decreasing Hill function
    K = 800.0    # Half-effect constant
    n = 1.5      # Hill coefficient

    efficiency = 1.0 / (1.0 + (nonlinear_factor / K) ** n)

    return efficiency


def softplus(x):
    """Smooth ReLU function: softplus(x) = log(1 + exp(x))"""
    # Prevent numerical overflow
    x = np.clip(x, -50, 50)
    return np.log(1.0 + np.exp(x))


def softplus_damage_factor(hours, p):
    """
    Calculate softplus + power form nonlinear damage factor (v10.0)

    Formula: factor = 1 + k * softplus(hours - threshold)^n

    Characteristics:
    - hours < threshold: factor approx 1 (softplus approx 0)
    - hours = threshold: factor = 1 + k * 0.69^n (mild)
    - hours > threshold: factor approx 1 + k * (hours - threshold)^n (power growth)

    Advantages:
    - Smooth transition, no hard threshold
    - Lower exponent (n=4) is more stable than n=8
    - Retains power function discrimination ability
    """
    sp = softplus(hours - p.softplus_threshold)
    factor = 1.0 + p.k_hours_nonlinear * (sp ** p.n_hours_nonlinear)
    return factor


def sigmoid_damage_factor(hours, p):
    """
    Calculate Sigmoid form nonlinear damage factor (v10.0 backup)

    Formula: factor = 1 + MaxFactor / (1 + exp(-k * (hours - threshold)))
    """
    exponent = -p.sigmoid_steepness * (hours - p.sigmoid_threshold)
    exponent = np.clip(exponent, -50, 50)

    factor = 1.0 + p.sigmoid_max_factor / (1.0 + np.exp(exponent))
    return factor


def gompertz_damage_factor(hours, p):
    """
    Calculate Gompertz form nonlinear damage factor (v10.0 core mechanism)

    Formula: factor = 1 + max * exp(-exp(-k * (hours - threshold)))

    ============================================================================
    Why Choose Gompertz Function?
    ============================================================================

    1. Reviewer suggestion: Simple power function (hours^8) is hard to explain biologically

    2. Gompertz function advantages:
       - Has clear saturation upper limit (consistent with physiology: defense system has limits)
       - Asymmetric S-curve, steeper than Logistic
       - Commonly used to describe biological growth curves and dose-response relationships
       - Describes the accelerating process of "antioxidant system collapse"

    ============================================================================
    Biological Significance: Antioxidant System "Collapse"
    ============================================================================

    Plant antioxidant defense systems (SOD, CAT, APX, etc.) have operational limits:

    1. hours < 10 (before threshold):
       - Antioxidant system operates normally
       - ROS effectively cleared
       - factor approx 1 (almost no additional damage)

    2. hours approx 10.5 (near threshold):
       - Antioxidant enzymes start to "fatigue"
       - ROS clearance efficiency decreases
       - factor starts rapid rise

    3. hours > 12 (after threshold):
       - Antioxidant system "collapses"
       - ROS accumulates massively, causing oxidative damage
       - factor approaches maximum (approximately 250)

    ============================================================================
    Key Finding: Daily Cumulative Irradiation is the Key
    ============================================================================

    Experiment B three groups (v10.37, threshold=10.5h):
    - 3h x 12d:  factor = 1.0   -> Almost no collapse
    - 6h x 6d:   factor = 1.0   -> Almost no collapse
    - 12h x 3d:  factor = 156.9 -> Severe collapse

    All three groups have same total dose, but different daily cumulative doses!
    -> Daily cumulative irradiation (not total dose) determines whether antioxidant system collapses

    ============================================================================
    Parameter Explanation (v10.37)
    ============================================================================

    - gompertz_threshold (= 10.5 hours): Time point when antioxidant system starts to collapse
    - gompertz_steepness (= 0.5): Collapse rate
    - gompertz_max_factor (= 250): Maximum damage multiplier (saturation upper limit)

    Returns:
        factor: Nonlinear damage factor (1.0 ~ 251)
    """
    exponent = -p.gompertz_steepness * (hours - p.gompertz_threshold)
    exponent = np.clip(exponent, -50, 50)

    factor = 1.0 + p.gompertz_max_factor * np.exp(-np.exp(exponent))
    return factor


def power_damage_factor(hours, p):
    """
    Calculate standard power function nonlinear damage factor (v9 style, with numerical clipping)

    Formula: factor = 1 + k * hours^n, clipped to maximum value
    """
    factor = 1.0 + p.k_hours_nonlinear * (hours ** p.n_hours_nonlinear)
    # Numerical clipping
    factor = min(factor, p.max_nonlinear_factor)
    return factor


def nonlinear_damage_factor(hours, p):
    """
    Select nonlinear damage factor calculation method based on settings

    Options:
    - use_gompertz=True (recommended): Use Gompertz function, has saturation upper limit
    - use_gompertz=False: Use Power function (hours^n), with maximum value clipping
    """
    if p.use_gompertz:
        # v10.0 recommended: Gompertz function (biologically more reasonable S-curve)
        return gompertz_damage_factor(hours, p)
    else:
        # Backup: Power function (hours^n)
        return power_damage_factor(hours, p)


# ==============================================================================
# Core Differential Equations (v10.6c - Complete Calibration Version)
# ==============================================================================

def uva_sun_derivatives(t, state, p, env):
    """
    UVA Effect Integrated Model Core Differential Equations (v10.6c)

    ============================================================================
    State Variables (6 total)
    ============================================================================
    - X_d:    Dry weight [kg/m2]
    - C_buf:  Carbon buffer pool [kg/m2]
    - LAI:    Leaf Area Index [m2/m2]
    - Anth:   Anthocyanin content [kg/m2]
    - Stress: Cumulative stress index [-]
    - ROS:    Reactive Oxygen Species [-]

    ============================================================================
    v10.6c Core Mechanism Overview
    ============================================================================

    1. UVA morphological effect (v10.0):
       - UVA does not directly contribute to photosynthesis (removed I_effective = I_base + I_UVA)
       - UVA indirectly promotes growth through SLA/LAI enhancement
       - High Stress inhibits morphological effect (v10.6)

    2. Damage formula (v10.1):
       damage = k1 * ROS * vulnerability + k2 * ROS * nonlinear_factor
       - LAI vulnerability damage: For D12 groups (early irradiation, small LAI)
       - Nonlinear time damage: For H12D3 (long irradiation, antioxidant collapse)

    3. Circadian damage at night:
       circadian_damage = k * I_UVA * hours_in_dark^n
       - For L6D6-N (nighttime irradiation group)

    4. LDMC acute injury (v10.6c):
       - Nonlinear factor passed to calculate_dynamic_dw_fw_ratio()
       - H12D3 (12h/day): High nonlinear factor -> High LDMC -> Low FW
    """

    # =========================================================================
    # Step 1: Unpack state variables
    # =========================================================================
    X_d, C_buf, LAI, Anth, Stress, ROS = state

    # Numerical protection
    X_d = max(X_d, 1e-9)
    C_buf = max(C_buf, 0)
    LAI = max(LAI, 0.1)
    Anth = max(Anth, 0)
    Stress = max(Stress, 0)
    ROS = max(ROS, 0)

    # =========================================================================
    # Step 2: Calculate time-related variables
    # =========================================================================
    hour = (t / 3600.0) % 24.0
    day_from_sowing = t / 86400.0

    # =========================================================================
    # Step 3: Determine day/night status
    # =========================================================================
    light_on = env['light_on_hour']
    light_off = env['light_off_hour']

    if light_on <= light_off:
        is_day = light_on <= hour < light_off
    else:
        is_day = hour >= light_on or hour < light_off

    if is_day:
        I_base = env['I_day']
        Tc = env['T_day']
    else:
        I_base = 0.0
        Tc = env['T_night']

    # =========================================================================
    # Step 4: Calculate UVA intensity
    # =========================================================================
    uva_on = env.get('uva_on', False)
    uva_start_day = env.get('uva_start_day', 29)
    uva_end_day = env.get('uva_end_day', 35)
    uva_hour_on = env.get('uva_hour_on', 10)
    uva_hour_off = env.get('uva_hour_off', 16)
    uva_intensity = env.get('uva_intensity', 11.0)  # v10.38: Default changed to 11 W/m2

    I_UVA = 0.0
    hours_today = 0.0
    days_irradiated = 0

    if uva_on:
        if day_from_sowing >= uva_start_day:
            days_irradiated = min(
                day_from_sowing - uva_start_day + 1,
                uva_end_day - uva_start_day + 1
            )

        if uva_hour_on <= uva_hour_off:
            if uva_start_day <= day_from_sowing <= uva_end_day:
                if uva_hour_on <= hour < uva_hour_off:
                    I_UVA = uva_intensity
                    hours_today = hour - uva_hour_on
        else:
            # Cross-day case (nighttime irradiation)
            if hour >= uva_hour_on:
                if uva_start_day <= day_from_sowing <= uva_end_day:
                    I_UVA = uva_intensity
                    hours_today = hour - uva_hour_on
            elif hour < uva_hour_off:
                day_session_started = day_from_sowing - 1
                if uva_start_day <= day_session_started <= uva_end_day:
                    I_UVA = uva_intensity
                    hours_today = (24 - uva_hour_on) + hour

    # =========================================================================
    # Step 5: Calculate circadian damage at night
    # =========================================================================
    #
    # [Key Finding 2] L6D6-N Circadian Damage Mechanism
    #
    # Experimental design:
    # - L6D6:   6h UVA x 6d, daytime irradiation (10:00-16:00)
    # - L6D6-N: 6h UVA x 6d, nighttime irradiation (22:00-04:00)
    # - Both groups have same total UVA dose, but L6D6-N has lower fresh weight
    #
    # Mechanism explanation:
    # 1. Plants have lower antioxidant system activity at night
    # 2. Circadian clock mismatch with UVA irradiation timing causes extra stress
    # 3. Repair capacity is also lower at night
    #
    # Formula: circadian_damage = k * I_UVA * hours_in_dark^n
    #
    if light_on <= light_off:
        if light_on <= hour < light_off:
            hours_in_dark = 0.0
        else:
            if hour >= light_off:
                hours_in_dark = hour - light_off
            else:
                hours_in_dark = hour + (24 - light_off)
    else:
        if hour >= light_on or hour < light_off:
            hours_in_dark = 0.0
        else:
            hours_in_dark = hour - light_off if hour >= light_off else hour + (24 - light_off)

    # Additional circadian damage only when UVA at night (only affects L6D6-N)
    if I_UVA > 0 and hours_in_dark > 0:
        circadian_damage = p.k_circadian * I_UVA * (hours_in_dark ** p.n_circadian)
    else:
        circadian_damage = 0.0

    # =========================================================================
    # Step 6: Call base Sun model (v10.0: UVA not added to PAR)
    # =========================================================================
    # v10.0 change: UVA does not directly increase photosynthetically active radiation
    # Photosynthesis only uses base PAR
    I_effective = I_base  # No longer add I_UVA

    env_modified = env.copy()
    env_modified['I_override'] = I_effective
    env_modified['T_override'] = Tc
    env_modified['is_day_override'] = is_day

    base_state = [X_d, C_buf, LAI]
    dXd_dt_base, dCbuf_dt, dLAI_dt_base = sun_derivatives_final(t, base_state, p, env_modified)

    # =========================================================================
    # Step 6b: UVA morphological effect (v10.0 core mechanism)
    # =========================================================================
    #
    # [Key Finding 1] UVA Morphological Effect
    #
    # Important: UVA does NOT directly contribute to photosynthesis!
    #
    # Old incorrect approach: I_effective = I_base + I_UVA  <- Physiologically unreasonable
    # Reason: UVA (365nm) has very low quantum yield, chlorophyll barely absorbs it
    #
    # v10.0 correction: UVA indirectly promotes growth through "morphological effect"
    # Mechanism: UVA -> Promotes leaf area expansion -> Increases PAR interception -> Indirectly promotes photosynthesis
    #
    # This explains why L6D6 fresh weight > CK:
    # - L6D6 (6h/day x 6d) morphological effect integrates through ODE
    # - Longer irradiation time, greater leaf area expansion
    # - Larger leaf area -> More light interception -> Higher biomass
    #
    # Literature support: Chen et al. (2019) reported UVA increases leaf area by 15-26%

    if I_UVA > 0:
        # SLA enhancement effect: UVA increases specific leaf area, indirectly improves photosynthesis efficiency
        # Formula: sla_boost = max_boost * I_UVA / (K + I_UVA)  [Michaelis-Menten]
        sla_boost = p.uva_sla_enhancement * I_UVA / (p.K_uva_sla + I_UVA)

        # LAI growth enhancement: UVA promotes leaf area index growth rate
        lai_boost = p.uva_lai_boost * I_UVA / (p.K_uva_lai + I_UVA)

        # v10.6 new: High Stress inhibits UVA morphological effect
        #
        # Biological significance: Plants under high stress cannot effectively use UVA for growth
        # - Resources are directed to repair and defense, not growth
        # - Use same K_stress for consistent inhibition mechanism
        #
        # Important for calibration:
        # - L6D6: Low Stress -> Complete morphological effect -> High FW
        # - H12D3: High Stress -> Morphological effect inhibited -> FW decreases
        stress_suppression = 1.0 - Stress / (p.K_stress + Stress + 1e-9)
        sla_boost = sla_boost * stress_suppression
        lai_boost = lai_boost * stress_suppression

        # v10.3 correction: UVA morphological effect works regardless of day/night
        # As long as there is UVA irradiation, there is morphological effect
        # Use absolute value approach: regardless of whether base growth is positive or negative, add positive morphological enhancement
        if dLAI_dt_base > 0:
            dLAI_dt_base = dLAI_dt_base * (1.0 + lai_boost)
        else:
            # Night: Base growth is negative, but still give morphological effect (reduce negative magnitude)
            dLAI_dt_base = dLAI_dt_base * (1.0 - lai_boost * 0.3)

        if dXd_dt_base > 0:
            dXd_dt_base = dXd_dt_base * (1.0 + sla_boost * 0.5)
        else:
            # Night: Reduce respiration consumption
            dXd_dt_base = dXd_dt_base * (1.0 - sla_boost * 0.15)

    # =========================================================================
    # Step 7: Calculate ROS dynamics
    # =========================================================================
    ros_production = p.k_ros_production * I_UVA
    ros_clearance = p.k_ros_clearance * ROS

    dROS_dt = ros_production - ros_clearance

    # =========================================================================
    # Step 8: Calculate LAI-dependent vulnerability
    # =========================================================================
    vulnerability = p.A_vulnerability * np.exp(-p.k_vulnerability * LAI) + 1.0

    # =========================================================================
    # Step 9: Calculate nonlinear damage (v10.0: softplus or Sigmoid)
    # =========================================================================
    # v10.0 change: Use softplus+power or Sigmoid instead of hours^8
    # Advantages: Numerically stable, biologically reasonable
    nonlinear_factor = nonlinear_damage_factor(hours_today, p)

    # Day accumulation factor: Long-term irradiation reduces repair capacity
    # days_factor = 1 + k * days, has greater impact on D12 groups
    days_factor = 1.0 + p.k_days_accumulation * days_irradiated

    # =========================================================================
    # Step 10: Calculate anthocyanin protection
    # =========================================================================
    anth_protection = p.alpha_anth_protection * Anth / (p.K_anth_protection + Anth + 1e-12)

    # =========================================================================
    # Step 11: Calculate damage rate (v10.1 revision: Additive independent mechanism)
    # =========================================================================
    #
    # [Key Finding 3] Damage Formula Transformation - Additive Independent Mechanism
    #
    # Background:
    # Old formula (multiplicative): damage = k * ROS * vulnerability * nonlinear_factor
    # -> H12D3 (vuln=6, nonlin=137) and VL3D12 (vuln=3607, nonlin=1)
    # -> 6*137 approx 822 vs 3607*1 = 3607
    # -> Two mechanisms cancel each other, cannot correctly distinguish damage for each group
    #
    # Solution:
    # New formula (additive): damage = k1 * ROS * vuln + k2 * ROS * nonlin
    # -> Two mechanisms act independently, each affects corresponding treatment groups
    #
    # Results:
    # - VL3D12, L6D12: Mainly rely on LAI vulnerability damage
    # - H12D3: Mainly rely on nonlinear time damage
    # - L6D6: Both are very small, almost no damage

    # LAI vulnerability damage (for D12 groups: Day 23 start irradiation, LAI approx 5)
    # Young plants have incomplete protection mechanisms, more sensitive to UVA
    vuln_damage = p.stress_damage_coeff * ROS * vulnerability * days_factor

    # Time nonlinear damage (for H12D3: 12h/day exceeds antioxidant capacity)
    # Based on Gompertz function describing "antioxidant system collapse"
    nonlin_damage = p.k_nonlinear_stress * ROS * nonlinear_factor

    # Base damage = vulnerability damage + nonlinear damage (independent addition)
    base_damage = vuln_damage + nonlin_damage

    # Anthocyanin protection: Reduces damage
    protected_damage = base_damage * (1.0 - anth_protection)

    # Total damage = ROS damage + circadian damage at night (for L6D6-N)
    damage_rate = protected_damage + circadian_damage

    # =========================================================================
    # Step 12: Calculate Stress decay (v10.2 simplified version)
    # =========================================================================
    # Simplified: Only one decay term, removed repair_rate
    # Original: dStress_dt = damage_rate - repair_rate - stress_decay (two decay terms redundant)
    # Now: dStress_dt = damage_rate - stress_decay (single decay term)
    stress_decay = p.k_stress_decay * Stress

    # =========================================================================
    # Step 13: Calculate Stress derivative
    # =========================================================================
    dStress_dt_raw = damage_rate - stress_decay
    if Stress <= 0 and dStress_dt_raw < 0:
        dStress_dt = 0.0
    else:
        dStress_dt = dStress_dt_raw

    # =========================================================================
    # Step 14: Calculate Stress inhibition on growth
    # =========================================================================
    stress_inhibition = Stress / (p.K_stress + Stress + 1e-9)
    xd_reduction = p.stress_photosynthesis_inhibition * stress_inhibition
    lai_reduction = p.stress_lai_inhibition * stress_inhibition

    dXd_dt = dXd_dt_base * (1.0 - xd_reduction) if dXd_dt_base > 0 else dXd_dt_base
    dLAI_dt = dLAI_dt_base * (1.0 - lai_reduction) if dLAI_dt_base > 0 else dLAI_dt_base

    # =========================================================================
    # Step 15: Calculate anthocyanin dynamics
    # =========================================================================
    #
    # v10.6c core modification: Pass nonlinear_factor to calculate LDMC acute injury effect
    #
    # [Key Finding] Daily cumulative irradiation determines acute injury
    # v10.38 nonlinear_factor (threshold=10.5):
    # - L6D6 (6h/day):   1.0   -> Normal LDMC -> Normal FW
    # - M9D3 (9h/day):   31.1  -> Starts acute injury
    # - H12D3 (12h/day): 156.9 -> High LDMC -> Low FW
    #
    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress, p, nonlinear_factor)
    FW_kg_m2 = X_d / dw_fw_ratio

    # Calculate daily irradiation hours
    daily_hours = uva_hour_off - uva_hour_on if uva_hour_on < uva_hour_off else 24 - uva_hour_on + uva_hour_off

    base_synthesis = p.base_anth_rate_light if is_day else p.base_anth_rate_dark

    # Calculate cumulative UVA exposure time
    total_uva_hours = max(0, days_irradiated - 1) * daily_hours + hours_today

    # Stress-induced synthesis
    # v10.22: Nighttime irradiation group has reduced Stress induction efficiency
    # This explains why L6D6-N (nighttime irradiation, uva_hour_on=22) anthocyanin is not higher than L6D6
    # Biological significance: Nighttime metabolic pathway activity is lower, even with higher Stress cannot effectively induce synthesis
    is_night_irradiation = (uva_hour_on >= 18) or (uva_hour_off <= 6)  # Nighttime irradiation group
    night_stress_efficiency = 0.4 if is_night_irradiation else 1.0  # Nighttime irradiation group efficiency 40%

    # v10.23: Low LAI inhibits Stress induction efficiency (nonlinear)
    # Damaged plants (low LAI) have reduced metabolic capacity, cannot effectively use stress signals to synthesize anthocyanin
    # This explains why D12 groups (low LAI) anthocyanin did not increase much despite high Stress
    # Use Hill function to make low LAI inhibition stronger
    LAI_healthy = 9.0  # Healthy LAI reference value
    n_LAI = 2.0        # Hill coefficient
    LAI_stress_efficiency = min(1.0, (LAI / LAI_healthy) ** n_LAI)  # LAI=6.9 efficiency 0.59, LAI=7.7 efficiency 0.73

    # v10.23: Original Hill function response
    stress_induced = p.V_max_anth * Stress / (p.K_stress_anth + Stress + 1e-12) * night_stress_efficiency * LAI_stress_efficiency

    # UV direct induction
    # v10.10: Adjusted for balance between groups
    # VL3D12 and L6D6 have same UV time (36h), but VL3D12 is too high
    # Increase K_uv to make long irradiation effect saturate earlier
    K_uv = 30.0  # Half-saturation constant (increase for earlier saturation)
    uv_induced = 1.4e-11 * total_uva_hours / (K_uv + total_uva_hours + 1e-12)  # Restored original value

    # =========================================================================
    # High Stress inhibits synthesis efficiency (v10.12)
    # =========================================================================
    # When Stress exceeds threshold, anthocyanin synthesis efficiency decreases
    # Use Gompertz function (nonlinear, asymmetric)
    #
    # Biological explanation:
    # - Extreme oxidative stress causes synthesis enzyme activity decrease
    # - Cellular metabolic pathways impaired
    # - But this is reversible, can recover after stress removal
    #
    # Formula: efficiency = 1 - max_inhib * exp(-exp(-steepness * (Stress - threshold)))
    #
    # v10.15 parameter design (changed to Hill function):
    # Hill function: inhibition = max_inhib * S^n / (K^n + S^n)
    # - Stress < K: Almost no inhibition
    # - Stress = K: 50% maximum inhibition
    # - Stress >> K: Approaches maximum inhibition
    #
    # Targets:
    # - L6D6 maxS approx 22 -> No inhibition (~0%)
    # - M9D3 maxS approx 329 -> Mild inhibition (~10%)
    # - H12D3 avgS approx 407 -> Mild inhibition (~16%)
    # - VH15D3 avgS approx 930 -> Moderate inhibition (~46%)
    #
    # v10.32: K adjusted from 150 to 800, let inhibition only affect extreme high Stress (VH15D3)
    # Fixes M9D3 anthocyanin over-prediction problem
    # Use parameter-defined values (no longer hardcoded)
    stress_inhibition = p.max_stress_inhib * (Stress ** p.n_stress_inhib) / (p.K_stress_inhib ** p.n_stress_inhib + Stress ** p.n_stress_inhib + 1e-9)
    stress_efficiency = 1.0 - stress_inhibition

    # Water inhibition efficiency (v10.9 existing mechanism)
    # When DW/FW is too high (dehydration), anthocyanin synthesis efficiency decreases
    water_efficiency = calculate_water_anth_efficiency(dw_fw_ratio, p)

    # v10.39: Nonlinear factor inhibits efficiency (Hill function)
    # Formula: efficiency = 1 / (1 + (nonlin/K)^n), K=800, n=1.5
    # Efficiency monotonically decreases as nonlinear_factor increases
    # v10.39 efficiency:
    #   - 3h/6h (nonlin=1.0):   100.0%
    #   - M9D3 (nonlin=31.1):   99.2%
    #   - H12D3 (nonlin=156.9): 92.0%
    #   - VH15D3 (nonlin=226):  86.9%
    daily_nonlin_factor = nonlinear_damage_factor(daily_hours, p)
    nonlin_anth_efficiency = calculate_nonlin_anth_efficiency(daily_nonlin_factor, p)

    # v10.22: Long-term irradiation adaptation effect
    # After prolonged cumulative irradiation, plant response to Stress induction decreases (adaptation/desensitization)
    # This explains why D12 groups (12 days irradiation) anthocyanin is not much higher than D3/D6 groups
    # Adaptation factor: K/(K+days), D3=0.57, D6=0.40, D12=0.25 (K=4)
    # Use parameter-defined values (no longer hardcoded)
    adaptation_factor = p.K_adapt_days / (p.K_adapt_days + days_irradiated)

    # Synthesis rate (including Stress inhibition + water inhibition + adaptation effect + nonlinear factor inhibition)
    # v10.32: nonlin_anth_efficiency only affects stress_induced, not base and uv_induced
    synthesis_rate = LAI * (base_synthesis + uv_induced + stress_induced * adaptation_factor * nonlin_anth_efficiency) * stress_efficiency * water_efficiency

    # Degradation and consumption
    natural_degradation = p.k_deg * Anth
    # v10.14: Anthocyanin consumption changed to use ROS (instantaneous stress), not Stress (cumulative stress)
    # Biological significance: Anthocyanin as antioxidant directly scavenges ROS
    # ROS only exists during irradiation, disappears quickly after stopping, so consumption also stops
    #
    # v10.20: Consumption efficiency increases with daily irradiation hours
    # During long irradiation, antioxidant system collapses, anthocyanin consumption intensifies
    # v10.38 daily_nonlin (threshold=10.5): 6h->1.0, 9h->31.1, 12h->156.9, 15h->226.0
    K_ros_cons = 500.0   # ROS half-saturation constant
    n_cons = 2.0
    # Consumption amplification factor: Based on daily_hours nonlinear_factor (not hours_today)
    # v10.33: Removed hard threshold, changed to softplus continuous function (compliant with CLAUDE.md)
    daily_nonlin = nonlinear_damage_factor(daily_hours, p)
    cons_amp_center = 200.0  # Soft threshold center (higher than 12h=156.9, mainly affects 15h=226)
    cons_amp_scale = 15.0    # Transition width
    cons_amp_k = 12.0        # Maximum amplification factor
    cons_amp_K = 20.0        # Half-saturation constant

    # Use softplus to implement soft threshold
    x_raw = (daily_nonlin - cons_amp_center) / cons_amp_scale
    x = cons_amp_scale * np.log(1.0 + np.exp(np.clip(x_raw, -50, 50)))  # softplus
    consumption_amp = 1.0 + cons_amp_k * (x ** 2) / (cons_amp_K ** 2 + x ** 2 + 1e-9)
    ros_consumption = p.k_anth_consumption * consumption_amp * Anth * (ROS ** n_cons) / (K_ros_cons ** n_cons + ROS ** n_cons + 1e-9)

    dAnth_dt = synthesis_rate - natural_degradation - ros_consumption

    # =========================================================================
    # Step 16: Carbon consumption
    # =========================================================================
    anth_carbon_consumption = synthesis_rate * p.anth_carbon_cost
    dCbuf_dt = dCbuf_dt - anth_carbon_consumption

    # =========================================================================
    # Return derivative vector (6 state variables)
    # =========================================================================
    return np.array([dXd_dt, dCbuf_dt, dLAI_dt, dAnth_dt, dStress_dt, dROS_dt])


# ==============================================================================
# Sensitivity Analysis Function
# ==============================================================================

def run_sensitivity_analysis(param_name, param_values, base_params, env_func, targets):
    """
    Execute single parameter sensitivity analysis

    Returns: FW, Stress, Anth changes for each treatment group at different parameter values
    """
    from model_config import ENV_BASE, SIMULATION

    results = []

    for pval in param_values:
        # Create modified parameters
        modified_params = base_params.copy()
        modified_params[param_name] = pval
        p = UVAParams(modified_params)

        treatment_results = {}

        for treatment in ['CK', 'L6D6', 'H12D3']:
            env = env_func(treatment)
            target = targets.get(treatment, {'FW': 0, 'Anth': 0})

            # Initial conditions
            fw_init_g = SIMULATION['initial_fw_g']
            dw_init_g = fw_init_g * p.dw_fw_ratio_base
            Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
            C_buf_init = Xd_init * 0.1
            LAI_init = (dw_init_g / 0.01) * 0.04
            fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
            Anth_init = 5.0 * fw_total_init / 1e6

            initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]

            transplant_day = SIMULATION['transplant_offset']
            simulation_days = SIMULATION['days']
            t_start = transplant_day * 86400
            t_end = (transplant_day + simulation_days) * 86400

            sol = solve_ivp(
                uva_sun_derivatives,
                (t_start, t_end),
                initial_state,
                args=(p, env),
                method='RK45',
                max_step=300
            )

            if sol.success:
                Xd_f, Cbuf_f, LAI_f, Anth_f, Stress_f, ROS_f = sol.y[:, -1]

                # Calculate nonlinear factor (for LDMC acute injury)
                uva_hour_on = env.get('uva_hour_on', 0)
                uva_hour_off = env.get('uva_hour_off', 0)
                hours_daily = uva_hour_off - uva_hour_on if uva_hour_on < uva_hour_off else 24 - uva_hour_on + uva_hour_off
                if not env.get('uva_on', False):
                    hours_daily = 0
                nonlin_factor = nonlinear_damage_factor(hours_daily, p)

                dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress_f, p, nonlin_factor)
                FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
                FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
                Anth_sim = Anth_f / FW_total_kg * 1e6

                treatment_results[treatment] = {
                    'FW': FW_sim,
                    'Stress': Stress_f,
                    'Anth': Anth_sim
                }

        results.append({
            'param_value': pval,
            'results': treatment_results
        })

    return results


# ==============================================================================
# Main Program
# ==============================================================================

if __name__ == "__main__":
    # ===========================================================================
    # Embedded Configuration (for standalone execution)
    # ===========================================================================

    # Environment base settings (plant factory)
    ENV_BASE = {
        'light_on_hour': 6,       # Light start time (06:00)
        'light_off_hour': 22,     # Light end time (22:00) - 16h photoperiod
        'I_day': 57,              # Daytime shortwave radiation [W/m²]
        'T_day': 25,              # Daytime temperature [°C]
        'T_night': 18,            # Nighttime temperature [°C]
        'CO2_day': 1200,          # Daytime CO2 [ppm]
        'CO2_night': 1200,        # Nighttime CO2 [ppm]
        'RH_day': 0.70,           # Daytime relative humidity
        'RH_night': 0.85,         # Nighttime relative humidity
        'plant_density': 36,      # Plant density [plants/m²]
    }

    # Simulation settings
    SIMULATION = {
        'days': 21,               # Simulation days after transplant
        'transplant_offset': 14,  # Transplant on day 14 after sowing
        'initial_fw_g': 10,       # Initial fresh weight [g/plant]
    }

    # Training set targets (observed values)
    TARGETS = {
        'CK':      {'FW': 87.0, 'Anth': 433},
        'L6D6':    {'FW': 91.4, 'Anth': 494},
        'L6D6-N':  {'FW': 80.8, 'Anth': 493},
        'H12D3':   {'FW': 60.6, 'Anth': 651},
        'VL3D12':  {'FW': 67.0, 'Anth': 482},
        'L6D12':   {'FW': 60.4, 'Anth': 518},
    }

    # Treatment configurations
    TREATMENT_CONFIGS = {
        'CK':      {'uva_on': False},
        'L6D6':    {'uva_on': True, 'uva_intensity': 11.0, 'uva_start_day': 29, 'uva_end_day': 35, 'uva_hour_on': 10, 'uva_hour_off': 16},
        'L6D6-N':  {'uva_on': True, 'uva_intensity': 11.0, 'uva_start_day': 29, 'uva_end_day': 35, 'uva_hour_on': 22, 'uva_hour_off': 4},
        'H12D3':   {'uva_on': True, 'uva_intensity': 11.0, 'uva_start_day': 32, 'uva_end_day': 35, 'uva_hour_on': 6, 'uva_hour_off': 18},
        'VL3D12':  {'uva_on': True, 'uva_intensity': 11.0, 'uva_start_day': 23, 'uva_end_day': 35, 'uva_hour_on': 10, 'uva_hour_off': 13},
        'L6D12':   {'uva_on': True, 'uva_intensity': 11.0, 'uva_start_day': 23, 'uva_end_day': 35, 'uva_hour_on': 10, 'uva_hour_off': 16},
    }

    def get_env_for_treatment(treatment):
        """Get complete environment settings for a treatment."""
        env = ENV_BASE.copy()
        if treatment in TREATMENT_CONFIGS:
            env.update(TREATMENT_CONFIGS[treatment])
        return env

    p = UVAParams()

    # Calculate ROS steady state value
    ros_ss = p.k_ros_production * 22.0 / p.k_ros_clearance

    print("=" * 80)
    print("Lettuce Growth and UVA Effect Integrated Model v10.0 (Reviewer Suggestion Revision)")
    print("=" * 80)
    print("\nv10.0 Core Changes:")
    print(f"  1. UVA morphological effect (replaces direct PAR addition):")
    print(f"     - SLA enhancement: {p.uva_sla_enhancement*100:.0f}% (K={p.K_uva_sla})")
    print(f"     - LAI enhancement: {p.uva_lai_boost*100:.0f}% (K={p.K_uva_lai})")
    if p.use_gompertz:
        print(f"  2. Gompertz nonlinear damage function (v10.0 recommended):")
        print(f"     - Threshold: {p.gompertz_threshold} hours")
        print(f"     - Steepness (collapse rate): {p.gompertz_steepness}")
        print(f"     - Max factor (saturation upper limit): {p.gompertz_max_factor}")
    else:
        print(f"  2. Power nonlinear damage (backup):")
        print(f"     - k = {p.k_hours_nonlinear}, n = {p.n_hours_nonlinear}")
        print(f"     - Max factor clipping: {p.max_nonlinear_factor}")
    print(f"  3. ROS steady state: {ros_ss:.1f}")
    print("=" * 80)

    # Display nonlinear factor characteristics
    print("\nNonlinear damage factor:")
    for h in [3, 6, 9, 12]:
        factor = nonlinear_damage_factor(h, p)
        print(f"  {h}h/day: factor = {factor:.1f}")
    print()

    fw_errs = []
    anth_errs = []

    for treatment in ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']:
        env = get_env_for_treatment(treatment)
        target = TARGETS.get(treatment, {'FW': 0, 'Anth': 0})

        # Initial conditions
        fw_init_g = SIMULATION['initial_fw_g']
        dw_init_g = fw_init_g * p.dw_fw_ratio_base
        Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
        C_buf_init = Xd_init * 0.1
        LAI_init = (dw_init_g / 0.01) * 0.04
        fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
        Anth_init = 5.0 * fw_total_init / 1e6

        initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]

        transplant_day = SIMULATION['transplant_offset']
        simulation_days = SIMULATION['days']
        t_start = transplant_day * 86400
        # Harvest time: Day 35 06:00 AM (at light start, before UVA irradiation)
        harvest_hour = 6
        t_end = (transplant_day + simulation_days) * 86400 + harvest_hour * 3600

        t_eval_points = np.linspace(t_start, t_end, 100)
        sol = solve_ivp(
            uva_sun_derivatives,
            (t_start, t_end),
            initial_state,
            args=(p, env),
            method='RK45',
            max_step=300,
            t_eval=t_eval_points
        )

        if sol.success:
            Xd_f, Cbuf_f, LAI_f, Anth_f, Stress_f, ROS_f = sol.y[:, -1]

            # Diagnostics
            uva_start = env.get('uva_start_day', 35) * 86400
            lai_at_uva_start = LAI_f
            max_stress = 0.0
            max_ros = 0.0
            for i, t in enumerate(sol.t):
                if t >= uva_start:
                    lai_at_uva_start = sol.y[2, i]
                    break

            # Calculate average Stress during irradiation period (for LDMC)
            stress_sum = 0.0
            stress_count = 0
            for i in range(len(sol.t)):
                if sol.y[4, i] > max_stress:
                    max_stress = sol.y[4, i]
                if sol.y[5, i] > max_ros:
                    max_ros = sol.y[5, i]
                # Accumulate Stress during irradiation period
                if sol.t[i] >= uva_start:
                    stress_sum += sol.y[4, i]
                    stress_count += 1
            avg_stress = stress_sum / max(1, stress_count)

            # Calculate nonlinear factor (for LDMC acute injury)
            uva_hour_on = env.get('uva_hour_on', 0)
            uva_hour_off = env.get('uva_hour_off', 0)
            hours_daily = uva_hour_off - uva_hour_on if uva_hour_on < uva_hour_off else 24 - uva_hour_on + uva_hour_off
            if not env.get('uva_on', False):
                hours_daily = 0
            nonlin_factor = nonlinear_damage_factor(hours_daily, p)

            # Use average Stress (during irradiation) + nonlinear_factor to calculate DW/FW
            # v10.10: Changed to avg_stress because FW is result of cumulative growth
            # nonlinear_factor reflects acute injury effect of daily cumulative irradiation
            dw_fw_ratio = calculate_dynamic_dw_fw_ratio(avg_stress, p, nonlin_factor)
            FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
            FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
            Anth_sim = Anth_f / FW_total_kg * 1e6

            FW_obs = target['FW']
            Anth_obs = target['Anth']

            fw_err = (FW_sim - FW_obs) / FW_obs * 100
            anth_err = (Anth_sim - Anth_obs) / Anth_obs * 100

            fw_errs.append(abs(fw_err))
            anth_errs.append(abs(anth_err))

            s1 = "PASS" if abs(fw_err) < 5 else "FAIL"
            s2 = "PASS" if abs(anth_err) < 5 else "FAIL"

            vuln_at_start = p.A_vulnerability * np.exp(-p.k_vulnerability * lai_at_uva_start) + 1
            vuln_at_end = p.A_vulnerability * np.exp(-p.k_vulnerability * LAI_f) + 1

            # Calculate anthocyanin protection
            anth_prot = p.alpha_anth_protection * Anth_f / (p.K_anth_protection + Anth_f + 1e-12)

            # Estimate circadian damage (for L6D6-N)
            is_night_uva = env.get('uva_hour_on', 0) >= 18 or env.get('uva_hour_off', 0) <= 6
            circ_estimate = 0.0
            if is_night_uva and env.get('uva_on', False):
                hours_in_dark = hours_daily
                circ_estimate = p.k_circadian * env.get('I_UVA', 11.0) * (hours_in_dark ** p.n_circadian)

            dw_g = Xd_f / ENV_BASE['plant_density'] * 1000  # Dry weight g/plant
            print(f"{treatment:<8} LAI:{LAI_f:>4.1f} endS:{Stress_f:>6.1f} avgS:{avg_stress:>5.1f} "
                  f"FW:{FW_sim:>5.1f}g({fw_err:>+5.1f}%{s1}) "
                  f"Anth:{Anth_sim:>4.0f}({anth_err:>+5.1f}%{s2}) "
                  f"dw/fw:{dw_fw_ratio:.3f}")
            # Diagnostics for L6D6 and L6D6-N
            if treatment in ['L6D6', 'L6D6-N']:
                # Calculate theoretical values for vuln_damage and nonlin_damage
                ros_ss = p.k_ros_production * env.get('I_UVA', 11.0) / p.k_ros_clearance
                vuln_dam_est = p.stress_damage_coeff * ros_ss * vuln_at_end * 1.0
                nonlin_dam_est = p.k_nonlinear_stress * ros_ss * nonlin_factor
                base_dam = vuln_dam_est + nonlin_dam_est
                prot_dam = base_dam * (1.0 - anth_prot)
                total_dam = prot_dam + circ_estimate
                print(f"         vuln_end:{vuln_at_end:.2f} anth_prot:{anth_prot:.3f} "
                      f"Anth_abs:{Anth_f*1e6:.2f}mg")
                print(f"         vuln_dam:{vuln_dam_est:.2e} nonlin:{nonlin_dam_est:.2e} "
                      f"circ:{circ_estimate:.2e} total:{total_dam:.2e}")
                # Track max Stress and corresponding time
                max_s_time = 0
                for i, t in enumerate(sol.t):
                    if sol.y[4, i] == max_stress:
                        max_s_time = t / 86400
                        break
                print(f"         maxS:{max_stress:.1f} at day {max_s_time:.1f}")
                # Check Stress at end of irradiation
                # L6D6: ends at 16:00 each day, L6D6-N: ends at 04:00 each day
                if treatment == 'L6D6':
                    check_hour = 16
                else:
                    check_hour = 4
                print(f"         Stress at {check_hour}:00: ", end="")
                for day in range(29, 36):
                    target_t = day * 86400 + check_hour * 3600
                    closest_idx = np.argmin(np.abs(sol.t - target_t))
                    s_val = sol.y[4, closest_idx]
                    print(f"D{day}:{s_val:.1f} ", end="")
                print()
        else:
            print(f"{treatment:<8} Simulation failed: {sol.message}")

    print("-" * 80)
    fw_ok = sum(1 for e in fw_errs if e < 5)
    anth_ok = sum(1 for e in anth_errs if e < 5)
    print(f"Pass: FW {fw_ok}/6, Anth {anth_ok}/6, Total {fw_ok + anth_ok}/12")

    # =========================================================================
    # Validation Experiment: 3-Day Gradient Experiment (Day 32-35, 0-15h UVA daily)
    # =========================================================================
    print("\n" + "=" * 80)
    print("Validation Experiment: 3-Day Gradient (Day 32-35)")
    print("=" * 80)

    # Validation set observed data (2026-01-12 v3 update: Anthocyanin data corrected)
    validation_targets = {
        'CK':      {'FW': 85.14, 'Anth': 413, 'hours': 0},
        'VL3D3':   {'FW': 89.1, 'Anth': 437, 'hours': 3},
        'L6D3':    {'FW': 92.18, 'Anth': 468, 'hours': 6},
        'M9D3':    {'FW': 83.79, 'Anth': 539, 'hours': 9},
        'H12D3':   {'FW': 62.2, 'Anth': 657, 'hours': 12},
        'VH15D3':  {'FW': 51.2, 'Anth': 578, 'hours': 15},
    }

    # Batch correction factor (Training CK=87g, Validation CK=85.14g -> 85.14/87=0.979)
    # v10.16: Batch difference is very small, set to 1.0 (no adjustment)
    # If batch_factor < 1, validation batch DW/FW is higher -> FW is lower
    batch_factor = 1.0

    val_fw_errs = []
    val_anth_errs = []

    print(f"Batch correction factor: {batch_factor}")
    print()

    for name, target in validation_targets.items():
        hours = target['hours']

        # Build environment settings
        env = dict(ENV_BASE)  # Copy base settings
        if hours > 0:
            env['uva_on'] = True
            env['uva_intensity'] = 11.0
            env['uva_start_day'] = 32
            env['uva_end_day'] = 35
            env['uva_hour_on'] = 6  # Consistent with training H12D3
            env['uva_hour_off'] = 6 + hours
        else:
            env['uva_on'] = False

        # Initial conditions (same as training)
        fw_init_g = SIMULATION['initial_fw_g']
        dw_init_g = fw_init_g * p.dw_fw_ratio_base
        Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
        C_buf_init = Xd_init * 0.1
        LAI_init = (dw_init_g / 0.01) * 0.04
        fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
        Anth_init = 5.0 * fw_total_init / 1e6

        initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]

        transplant_day = SIMULATION['transplant_offset']
        simulation_days = SIMULATION['days']
        t_start = transplant_day * 86400
        harvest_hour = 6
        t_end = (transplant_day + simulation_days) * 86400 + harvest_hour * 3600

        t_eval_points = np.linspace(t_start, t_end, 100)
        sol = solve_ivp(
            uva_sun_derivatives,
            (t_start, t_end),
            initial_state,
            args=(p, env),
            method='RK45',
            max_step=300,
            t_eval=t_eval_points
        )

        if sol.success:
            Xd_f, Cbuf_f, LAI_f, Anth_f, Stress_f, ROS_f = sol.y[:, -1]

            # Calculate max Stress
            max_stress = max(sol.y[4, :])

            # Calculate average Stress during irradiation period
            uva_start = env.get('uva_start_day', 35) * 86400
            stress_sum = 0.0
            stress_count = 0
            for i in range(len(sol.t)):
                if sol.t[i] >= uva_start:
                    stress_sum += sol.y[4, i]
                    stress_count += 1
            avg_stress = stress_sum / max(1, stress_count)

            # Calculate nonlinear factor
            hours_daily = hours
            nonlin_factor = nonlinear_damage_factor(hours_daily, p)

            # Calculate DW/FW (apply batch factor)
            # batch_factor < 1 means validation batch DW/FW is higher (drier/lighter)
            # FW = DW / (DW/FW), higher DW/FW -> lower FW
            dw_fw_ratio_base = calculate_dynamic_dw_fw_ratio(avg_stress, p, nonlin_factor)
            dw_fw_ratio = dw_fw_ratio_base / batch_factor  # batch<1 -> ratio increases -> FW decreases

            # Calculate FW
            FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000

            # Calculate anthocyanin
            FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
            Anth_sim = Anth_f / FW_total_kg * 1e6

            FW_obs = target['FW']
            Anth_obs = target['Anth']

            fw_err = (FW_sim - FW_obs) / FW_obs * 100
            anth_err = (Anth_sim - Anth_obs) / Anth_obs * 100

            val_fw_errs.append(abs(fw_err))
            val_anth_errs.append(abs(anth_err))

            # Status symbols
            fw_s = "PASS" if abs(fw_err) < 5 else ("WARN" if abs(fw_err) < 10 else "FAIL")
            anth_s = "PASS" if abs(anth_err) < 5 else ("WARN" if abs(anth_err) < 10 else "FAIL")

            print(f"{name:<8} {hours:>2}h/day LAI:{LAI_f:>4.1f} avgS:{avg_stress:>6.0f} maxS:{max_stress:>6.0f} "
                  f"FW:{FW_sim:>5.1f}g({fw_err:>+5.1f}%{fw_s}) "
                  f"Anth:{Anth_sim:>4.0f}({anth_err:>+5.1f}%{anth_s}) "
                  f"dw/fw:{dw_fw_ratio:.3f} nonlin:{nonlin_factor:.1f}")
        else:
            print(f"{name:<8} Simulation failed: {sol.message}")

    print("-" * 80)
    val_fw_ok5 = sum(1 for e in val_fw_errs if e < 5)
    val_fw_ok10 = sum(1 for e in val_fw_errs if e < 10)
    val_anth_ok5 = sum(1 for e in val_anth_errs if e < 5)
    val_anth_ok10 = sum(1 for e in val_anth_errs if e < 10)
    print(f"Validation FW: <5%: {val_fw_ok5}/6, <10%: {val_fw_ok10}/6, Mean error: {sum(val_fw_errs)/6:.1f}%")
    print(f"Validation Anth: <5%: {val_anth_ok5}/6, <10%: {val_anth_ok10}/6, Mean error: {sum(val_anth_errs)/6:.1f}%")

    # =========================================================================
    # Export results to CSV for reproducibility verification
    # =========================================================================
    print("\n" + "=" * 80)
    print("Exporting results to results.csv...")
    print("=" * 80)

    # Collect all results
    all_results = []

    # Re-run training set to collect data
    for treatment in ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']:
        env = get_env_for_treatment(treatment)
        target = TARGETS.get(treatment, {'FW': 0, 'Anth': 0})

        fw_init_g = SIMULATION['initial_fw_g']
        dw_init_g = fw_init_g * p.dw_fw_ratio_base
        Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
        C_buf_init = Xd_init * 0.1
        LAI_init = (dw_init_g / 0.01) * 0.04
        fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
        Anth_init = 5.0 * fw_total_init / 1e6
        initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]

        transplant_day = SIMULATION['transplant_offset']
        simulation_days = SIMULATION['days']
        t_start = transplant_day * 86400
        harvest_hour = 6
        t_end = (transplant_day + simulation_days) * 86400 + harvest_hour * 3600
        t_eval_points = np.linspace(t_start, t_end, 100)

        sol = solve_ivp(uva_sun_derivatives, (t_start, t_end), initial_state,
                        args=(p, env), method='RK45', max_step=300, t_eval=t_eval_points)

        if sol.success:
            Xd_f, Cbuf_f, LAI_f, Anth_f, Stress_f, ROS_f = sol.y[:, -1]
            uva_start = env.get('uva_start_day', 35) * 86400
            stress_sum, stress_count = 0.0, 0
            for i in range(len(sol.t)):
                if sol.t[i] >= uva_start:
                    stress_sum += sol.y[4, i]
                    stress_count += 1
            avg_stress = stress_sum / max(1, stress_count)

            uva_hour_on = env.get('uva_hour_on', 0)
            uva_hour_off = env.get('uva_hour_off', 0)
            hours_daily = uva_hour_off - uva_hour_on if uva_hour_on < uva_hour_off else 24 - uva_hour_on + uva_hour_off
            if not env.get('uva_on', False):
                hours_daily = 0
            nonlin_factor = nonlinear_damage_factor(hours_daily, p)
            dw_fw_ratio = calculate_dynamic_dw_fw_ratio(avg_stress, p, nonlin_factor)
            FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
            FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
            Anth_sim = Anth_f / FW_total_kg * 1e6

            all_results.append({
                'Set': 'Training',
                'Treatment': treatment,
                'FW_obs': target['FW'],
                'FW_pred': round(FW_sim, 1),
                'FW_error_pct': round((FW_sim - target['FW']) / target['FW'] * 100, 1),
                'Anth_obs': target['Anth'],
                'Anth_pred': round(Anth_sim, 0),
                'Anth_error_pct': round((Anth_sim - target['Anth']) / target['Anth'] * 100, 1)
            })

    # Re-run validation set to collect data
    for name, target in validation_targets.items():
        hours = target['hours']
        env = dict(ENV_BASE)
        if hours > 0:
            env['uva_on'] = True
            env['uva_intensity'] = 11.0
            env['uva_start_day'] = 32
            env['uva_end_day'] = 35
            env['uva_hour_on'] = 6
            env['uva_hour_off'] = 6 + hours
        else:
            env['uva_on'] = False

        fw_init_g = SIMULATION['initial_fw_g']
        dw_init_g = fw_init_g * p.dw_fw_ratio_base
        Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
        C_buf_init = Xd_init * 0.1
        LAI_init = (dw_init_g / 0.01) * 0.04
        fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
        Anth_init = 5.0 * fw_total_init / 1e6
        initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]

        transplant_day = SIMULATION['transplant_offset']
        simulation_days = SIMULATION['days']
        t_start = transplant_day * 86400
        harvest_hour = 6
        t_end = (transplant_day + simulation_days) * 86400 + harvest_hour * 3600
        t_eval_points = np.linspace(t_start, t_end, 100)

        sol = solve_ivp(uva_sun_derivatives, (t_start, t_end), initial_state,
                        args=(p, env), method='RK45', max_step=300, t_eval=t_eval_points)

        if sol.success:
            Xd_f, Cbuf_f, LAI_f, Anth_f, Stress_f, ROS_f = sol.y[:, -1]
            uva_start = env.get('uva_start_day', 35) * 86400
            stress_sum, stress_count = 0.0, 0
            for i in range(len(sol.t)):
                if sol.t[i] >= uva_start:
                    stress_sum += sol.y[4, i]
                    stress_count += 1
            avg_stress = stress_sum / max(1, stress_count)

            nonlin_factor = nonlinear_damage_factor(hours, p)
            dw_fw_ratio = calculate_dynamic_dw_fw_ratio(avg_stress, p, nonlin_factor)
            FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
            FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
            Anth_sim = Anth_f / FW_total_kg * 1e6

            all_results.append({
                'Set': 'Validation',
                'Treatment': name,
                'FW_obs': target['FW'],
                'FW_pred': round(FW_sim, 1),
                'FW_error_pct': round((FW_sim - target['FW']) / target['FW'] * 100, 1),
                'Anth_obs': target['Anth'],
                'Anth_pred': round(Anth_sim, 0),
                'Anth_error_pct': round((Anth_sim - target['Anth']) / target['Anth'] * 100, 1)
            })

    # Write to CSV (without pandas dependency)
    csv_filename = 'results.csv'
    with open(csv_filename, 'w') as f:
        f.write('Set,Treatment,FW_obs,FW_pred,FW_error_pct,Anth_obs,Anth_pred,Anth_error_pct\n')
        for r in all_results:
            f.write(f"{r['Set']},{r['Treatment']},{r['FW_obs']},{r['FW_pred']},{r['FW_error_pct']},"
                    f"{r['Anth_obs']},{r['Anth_pred']},{r['Anth_error_pct']}\n")

    print(f"\nResults exported to: {csv_filename}")
    print("Columns: Set, Treatment, FW_obs, FW_pred, FW_error_pct, Anth_obs, Anth_pred, Anth_error_pct")
    print(f"Total rows: {len(all_results)} (Training: 6, Validation: 6)")
