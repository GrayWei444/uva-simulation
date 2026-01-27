"""
================================================================================
Lettuce Growth and UVA Effect Integrated Model
================================================================================

Carbon Competition + AOX/Anthocyanin Framework

================================================================================
Core Innovation: Carbon Competition Mechanism
================================================================================

Key concept: AOX (Antioxidant) synthesis consumes carbon buffer (C_buf)

1. AOX is the state variable for antioxidants
   - Includes all antioxidant compounds: anthocyanins, flavonoids, ascorbate, etc.
   - Anthocyanin is approximately 18% of total AOX

2. Carbon competition:
   - dC_buf/dt = photosynthesis - respiration - growth - AOX_synthesis
   - Under stress, plant allocates more carbon to defense (AOX) vs growth (X_d)
   - This creates a physiologically meaningful trade-off

3. Literature support for carbon competition:
   - Herms & Mattson (1992): Growth-Differentiation Balance Hypothesis
     DOI: 10.1086/285343
   - Monson et al. (2022): Coordinated resource allocation to growth-defense tradeoffs
     DOI: 10.1111/nph.17773
   - Vogt (2010): Phenylpropanoid Biosynthesis - ~20% photosynthate to phenylpropanoids
     DOI: 10.1093/mp/ssp106
   - Gershenzon (1994): Metabolic costs of terpenoid accumulation
     DOI: 10.1007/BF02059810

================================================================================
State Variables (6 total)
================================================================================
- X_d:    Dry weight [kg/m2]
- C_buf:  Carbon buffer pool [kg/m2]
- LAI:    Leaf Area Index [m2/m2]
- AOX:    Antioxidant content [kg/m2] (includes anthocyanins, flavonoids, etc.)
- Stress: Cumulative stress index [-]
- ROS:    Reactive Oxygen Species [-]

Output:
- Anthocyanin = AOX * 0.18 (18% of total antioxidants)
================================================================================
"""

import numpy as np
from scipy.integrate import solve_ivp

# Import base Sun model
from lettuce_uva_carbon_complete_model import SunParams as BaseSunParams
from lettuce_uva_carbon_complete_model import sun_derivatives_final


# ==============================================================================
# Model Parameter Definitions
# ==============================================================================

ALL_PARAMS = {
    # --------------------------------------------------------------------------
    # 1. Base Photosynthesis Parameters
    # --------------------------------------------------------------------------
    'c_alpha': 0.54,                     # Photosynthesis efficiency coefficient [-]

    # --------------------------------------------------------------------------
    # 2. UVA Morphological Effect Parameters
    # --------------------------------------------------------------------------
    'uva_sla_enhancement': 5.00,         # UVA max enhancement on SLA
    'K_uva_sla': 7.5,                    # SLA enhancement half-saturation UVA intensity [W/m2]
    'uva_lai_boost': 1.70,               # UVA max enhancement on LAI growth
    'K_uva_lai': 7.5,                    # LAI enhancement half-saturation UVA intensity [W/m2]

    # --------------------------------------------------------------------------
    # 3. ROS Dynamics Parameters
    # --------------------------------------------------------------------------
    'k_ros_production': 0.010,           # ROS production coefficient [ROS/(W/m2*s)]
    'k_ros_clearance': 5e-4,             # ROS clearance coefficient [1/s]

    # --------------------------------------------------------------------------
    # 4. Stress Damage Parameters
    # --------------------------------------------------------------------------
    'stress_damage_coeff': 1.6e-7,       # Base damage coefficient
    'A_vulnerability': 85000000.0,       # Vulnerability amplitude [-]
    'k_vulnerability': 2.0,              # Vulnerability decay coefficient [1/(m2/m2)]

    # Gompertz function parameters for nonlinear damage
    'gompertz_max_factor': 250.0,
    'gompertz_threshold': 10.5,
    'gompertz_steepness': 0.5,

    # Anthocyanin/AOX protection
    # K must be scaled by 1/0.18 since AOX is 5.56x larger than Anth
    'alpha_aox_protection': 0.5,         # Maximum protection efficiency (0~1) [-]
    'K_aox_protection': 2.78e-5,         # = 5e-6 / 0.18 (scaled for AOX)

    # Circadian damage
    'k_circadian': 3.0e-6,
    'n_circadian': 2.0,

    # Nonlinear damage
    'k_nonlinear_stress': 5.0e-6,

    # --------------------------------------------------------------------------
    # 5. Stress Decay Parameters
    # --------------------------------------------------------------------------
    'k_stress_decay': 2.14e-5,           # Stress decay coefficient [1/s] (half-life 9 hours)

    # --------------------------------------------------------------------------
    # 6. Stress Inhibition on Growth
    # --------------------------------------------------------------------------
    'stress_photosynthesis_inhibition': 0.85,
    'stress_lai_inhibition': 0.80,
    'K_stress': 50.0,

    # --------------------------------------------------------------------------
    # 7. AOX Synthesis Parameters
    # --------------------------------------------------------------------------
    # AOX = Total antioxidants (anthocyanins, flavonoids, ascorbate, etc.)
    # Anthocyanin ~ 18% of AOX
    #
    # IMPORTANT: Since Anth = AOX * 0.18, AOX rates must be 1/0.18 = 5.56x of Anth rates
    # to achieve the same final anthocyanin concentration
    #
    'base_aox_rate_light': 3.53e-9,      # = 6.35e-10 / 0.18 (scaled for AOX)
    'base_aox_rate_dark': 1.77e-9,       # = 3.18e-10 / 0.18 (scaled for AOX)
    'V_max_aox': 1.45e-8,                # Slightly reduced for better Anth balance
    'K_stress_aox': 100.0,               # Half-saturation for stress-induced synthesis
    'k_aox_deg': 3.02e-6,                # AOX degradation rate (same as Anth)

    # --------------------------------------------------------------------------
    # 7b. Carbon Competition (Growth-Defense Trade-off)
    # --------------------------------------------------------------------------
    # AOX synthesis consumes carbon from C_buf
    # This creates growth vs defense trade-off
    #
    # Literature support:
    # - Herms & Mattson (1992): Growth-Differentiation Balance Hypothesis
    #   DOI: 10.1086/285343
    # - Monson et al. (2022): Coordinated resource allocation
    #   DOI: 10.1111/nph.17773
    #
    # Carbon cost estimation:
    # - Anthocyanins are phenylpropanoids: ~40% carbon by mass
    # - Synthesis requires additional ATP/NADPH (metabolic cost ~2x substrate)
    # - Estimated total carbon cost: 0.8-1.2 kg C per kg AOX synthesized
    #
    'aox_carbon_cost': 1.0,              # kg C_buf consumed per kg AOX synthesized
    'carbon_competition_K': 1e-8,        # Half-saturation for AOX carbon effect
    'stress_competition_K': 21.0,        # Half-saturation for stress competition effect
    'stress_competition_max': 0.225,     # Maximum stress-based competition
    'carbon_competition_max': 0.30,      # Maximum AOX-based competition
    'max_cbuf_consumption': 0.10,        # Max fraction of C_buf consumed per timestep

    # --------------------------------------------------------------------------
    # 7c. AOX Water Inhibition Parameters
    # --------------------------------------------------------------------------
    'water_aox_threshold': 0.055,
    'water_aox_K': 0.020,
    'water_aox_max_inhib': 0.50,
    'water_n': 2.0,                      # Hill coefficient for water inhibition

    # --------------------------------------------------------------------------
    # 7d. AOX Stress Inhibition Parameters
    # --------------------------------------------------------------------------
    'K_stress_inhib': 150.0,
    'n_stress_inhib': 2.0,
    'max_stress_inhib': 0.80,

    # --------------------------------------------------------------------------
    # 7e. AOX Adaptation Effect Parameters
    # --------------------------------------------------------------------------
    'K_adapt_days': 4.0,

    # --------------------------------------------------------------------------
    # 7f. UV-Induced AOX Synthesis Parameters
    # --------------------------------------------------------------------------
    'k_uv_aox': 7.78e-11,                # UV induction coefficient (= 1.4e-11 / 0.18)
    'K_uv_hours': 30.0,                  # Half-saturation for UV hours effect

    # --------------------------------------------------------------------------
    # 7g. LAI Efficiency Parameters
    # --------------------------------------------------------------------------
    'LAI_healthy': 9.0,                  # Reference LAI for healthy plant
    'n_LAI_eff': 2.0,                    # Hill coefficient for LAI efficiency

    # --------------------------------------------------------------------------
    # 7h. Night Irradiation Efficiency
    # --------------------------------------------------------------------------
    'night_stress_efficiency': 0.4,      # Efficiency multiplier for night irradiation

    # --------------------------------------------------------------------------
    # 7i. Nonlinear AOX Efficiency Parameters
    # --------------------------------------------------------------------------
    'K_nonlin_aox': 800.0,               # Half-saturation for nonlinear efficiency
    'n_nonlin_aox': 1.5,                 # Hill coefficient for nonlinear efficiency

    # --------------------------------------------------------------------------
    # 8. AOX Consumption Parameters (ROS scavenging)
    # --------------------------------------------------------------------------
    # Since AOX = Anth / 0.18, consumption rate should be k_anth * 0.18
    # to consume the same absolute amount
    'k_aox_consumption': 1.8e-7,         # = 1.0e-6 * 0.18
    'K_ros_consumption': 500.0,          # ROS half-saturation for consumption
    'n_ros_consumption': 2.0,            # Hill coefficient for ROS consumption
    # Consumption amplification for extreme daily hours
    'cons_amp_center': 200.0,            # Center for softplus activation
    'cons_amp_scale': 15.0,              # Scale for softplus
    'cons_amp_k': 12.0,                  # Amplification coefficient
    'cons_amp_K': 20.0,                  # Amplification half-saturation

    # --------------------------------------------------------------------------
    # 9. LDMC Parameters
    # --------------------------------------------------------------------------
    'dw_fw_ratio_base': 0.05,
    'ldmc_stress_sensitivity': 0.45,
    'K_ldmc': 1400.0,
    'dw_fw_ratio_max': 0.080,
    # Acute LDMC parameters (for high nonlinear factor)
    'acute_center': 50.0,                # Softplus center for acute LDMC
    'acute_scale': 10.0,                 # Softplus scale
    'acute_k': 9.0,                      # Acute LDMC coefficient
    'acute_K': 120.0,                    # Acute LDMC half-saturation
    'acute_n': 2.0,                      # Acute LDMC Hill coefficient

    # --------------------------------------------------------------------------
    # 10. Anthocyanin Fraction
    # --------------------------------------------------------------------------
    # Anthocyanin is a fraction of total AOX
    'anthocyanin_fraction': 0.18,        # 18% of AOX is anthocyanin
}


# ==============================================================================
# Model Parameter Class
# ==============================================================================

class UVAParams(BaseSunParams):
    """
    UVA Effect Model Parameter Class (Carbon Competition)
    """

    def __init__(self, params=None):
        super().__init__()
        if params is None:
            params = ALL_PARAMS

        # 1. Base photosynthesis parameters
        self.c_alpha = params['c_alpha']

        # 2. UVA morphological effect parameters
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

        # Gompertz parameters
        self.gompertz_max_factor = params['gompertz_max_factor']
        self.gompertz_threshold = params['gompertz_threshold']
        self.gompertz_steepness = params['gompertz_steepness']

        self.alpha_aox_protection = params['alpha_aox_protection']
        self.K_aox_protection = params['K_aox_protection']
        self.k_circadian = params['k_circadian']
        self.n_circadian = params['n_circadian']
        self.k_nonlinear_stress = params['k_nonlinear_stress']

        # 5. Stress decay parameters
        self.k_stress_decay = params['k_stress_decay']

        # 6. Stress inhibition on growth
        self.stress_photosynthesis_inhibition = params['stress_photosynthesis_inhibition']
        self.stress_lai_inhibition = params['stress_lai_inhibition']
        self.K_stress = params['K_stress']

        # 7. AOX synthesis parameters
        self.base_aox_rate_light = params['base_aox_rate_light']
        self.base_aox_rate_dark = params['base_aox_rate_dark']
        self.V_max_aox = params['V_max_aox']
        self.K_stress_aox = params['K_stress_aox']
        self.k_aox_deg = params['k_aox_deg']

        # 7b. Carbon competition
        self.aox_carbon_cost = params['aox_carbon_cost']
        self.carbon_competition_K = params['carbon_competition_K']
        self.stress_competition_K = params['stress_competition_K']
        self.stress_competition_max = params['stress_competition_max']
        self.carbon_competition_max = params['carbon_competition_max']
        self.max_cbuf_consumption = params['max_cbuf_consumption']

        # 7c. AOX water inhibition
        self.water_aox_threshold = params['water_aox_threshold']
        self.water_aox_K = params['water_aox_K']
        self.water_aox_max_inhib = params['water_aox_max_inhib']
        self.water_n = params['water_n']

        # 7f. UV-induced AOX synthesis
        self.k_uv_aox = params['k_uv_aox']
        self.K_uv_hours = params['K_uv_hours']

        # 7g. LAI efficiency
        self.LAI_healthy = params['LAI_healthy']
        self.n_LAI_eff = params['n_LAI_eff']

        # 7h. Night irradiation efficiency
        self.night_stress_efficiency = params['night_stress_efficiency']

        # 7i. Nonlinear AOX efficiency
        self.K_nonlin_aox = params['K_nonlin_aox']
        self.n_nonlin_aox = params['n_nonlin_aox']

        # 7d. AOX stress inhibition
        self.K_stress_inhib = params['K_stress_inhib']
        self.n_stress_inhib = params['n_stress_inhib']
        self.max_stress_inhib = params['max_stress_inhib']

        # 7e. AOX adaptation
        self.K_adapt_days = params['K_adapt_days']

        # 8. AOX consumption
        self.k_aox_consumption = params['k_aox_consumption']
        self.K_ros_consumption = params['K_ros_consumption']
        self.n_ros_consumption = params['n_ros_consumption']
        self.cons_amp_center = params['cons_amp_center']
        self.cons_amp_scale = params['cons_amp_scale']
        self.cons_amp_k = params['cons_amp_k']
        self.cons_amp_K = params['cons_amp_K']

        # 9. LDMC parameters
        self.dw_fw_ratio_base = params['dw_fw_ratio_base']
        self.ldmc_stress_sensitivity = params['ldmc_stress_sensitivity']
        self.K_ldmc = params['K_ldmc']
        self.dw_fw_ratio_max = params['dw_fw_ratio_max']
        self.acute_center = params['acute_center']
        self.acute_scale = params['acute_scale']
        self.acute_k = params['acute_k']
        self.acute_K = params['acute_K']
        self.acute_n = params['acute_n']

        # 10. Anthocyanin fraction
        self.anthocyanin_fraction = params['anthocyanin_fraction']


# ==============================================================================
# Utility Functions
# ==============================================================================

def calculate_dynamic_dw_fw_ratio(Stress, p, nonlinear_factor=1.0):
    """
    Calculate dynamic DW:FW ratio based on Stress and nonlinear factor
    """
    stress_effect = p.ldmc_stress_sensitivity * Stress / (p.K_ldmc + Stress + 1e-9)

    x_raw = (nonlinear_factor - p.acute_center) / p.acute_scale
    x = p.acute_scale * np.log(1.0 + np.exp(np.clip(x_raw, -50, 50)))
    acute_factor = 1.0 + p.acute_k * (x ** p.acute_n) / (p.acute_K ** p.acute_n + x ** p.acute_n + 1e-9)

    ratio = p.dw_fw_ratio_base * (1.0 + stress_effect * acute_factor)
    return min(ratio, p.dw_fw_ratio_max)


def calculate_water_aox_efficiency(dw_fw_ratio, p):
    """
    Calculate water status effect on AOX synthesis efficiency
    """
    base = p.water_aox_threshold
    K = p.water_aox_K
    n = p.water_n

    if dw_fw_ratio <= base:
        return 1.0

    x = dw_fw_ratio - base
    inhibition = p.water_aox_max_inhib * (x ** n) / (K ** n + x ** n)
    efficiency = 1.0 - inhibition

    return efficiency


def calculate_nonlin_aox_efficiency(nonlinear_factor, p):
    """
    Calculate AOX synthesis efficiency based on nonlinear factor
    Uses Hill function: efficiency = 1 / (1 + (nonlin/K)^n)
    """
    efficiency = 1.0 / (1.0 + (nonlinear_factor / p.K_nonlin_aox) ** p.n_nonlin_aox)
    return efficiency


def nonlinear_damage_factor(hours, p):
    """
    Calculate Gompertz form nonlinear damage factor
    Formula: factor = 1 + max * exp(-exp(-k * (hours - threshold)))

    Literature: Gompertz function commonly used for growth/damage modeling
    """
    exponent = -p.gompertz_steepness * (hours - p.gompertz_threshold)
    exponent = np.clip(exponent, -50, 50)
    factor = 1.0 + p.gompertz_max_factor * np.exp(-np.exp(exponent))
    return factor


# ==============================================================================
# Core Differential Equations
# ==============================================================================

def uva_sun_derivatives(t, state, p, env):
    """
    UVA Effect Integrated Model Core Differential Equations

    State Variables (6 total):
    - X_d:    Dry weight [kg/m2]
    - C_buf:  Carbon buffer pool [kg/m2]
    - LAI:    Leaf Area Index [m2/m2]
    - AOX:    Antioxidant content [kg/m2]
    - Stress: Cumulative stress index [-]
    - ROS:    Reactive Oxygen Species [-]

    Key Innovation: Carbon Competition
    - AOX synthesis consumes C_buf
    - dC_buf/dt = photosynthesis - respiration - growth - AOX_synthesis*carbon_cost
    """

    # =========================================================================
    # Step 1: Unpack state variables
    # =========================================================================
    X_d, C_buf, LAI, AOX, Stress, ROS = state

    # Numerical protection
    X_d = max(X_d, 1e-9)
    C_buf = max(C_buf, 0)
    LAI = max(LAI, 0.1)
    AOX = max(AOX, 0)
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
    uva_intensity = env.get('uva_intensity', 11.0)

    I_UVA = 0.0
    hours_today = 0.0
    days_irradiated = 0

    if uva_on:
        # Use integer day for counting completed irradiation days
        # This ensures total_uva_hours only counts actual irradiation time
        day_int = int(day_from_sowing)
        if day_int >= uva_start_day:
            days_irradiated = min(
                day_int - uva_start_day + 1,
                uva_end_day - uva_start_day + 1
            )

        if uva_hour_on <= uva_hour_off:
            if uva_start_day <= day_from_sowing <= uva_end_day:
                if uva_hour_on <= hour < uva_hour_off:
                    I_UVA = uva_intensity
                    hours_today = hour - uva_hour_on
        else:
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

    if I_UVA > 0 and hours_in_dark > 0:
        circadian_damage = p.k_circadian * I_UVA * (hours_in_dark ** p.n_circadian)
    else:
        circadian_damage = 0.0

    # =========================================================================
    # Step 6: Call base Sun model
    # =========================================================================
    I_effective = I_base

    env_modified = env.copy()
    env_modified['I_override'] = I_effective
    env_modified['T_override'] = Tc
    env_modified['is_day_override'] = is_day

    base_state = [X_d, C_buf, LAI]
    dXd_dt_base, dCbuf_dt, dLAI_dt_base = sun_derivatives_final(t, base_state, p, env_modified)

    # =========================================================================
    # Step 6b: UVA morphological effect
    # =========================================================================
    if I_UVA > 0:
        sla_boost = p.uva_sla_enhancement * I_UVA / (p.K_uva_sla + I_UVA)
        lai_boost = p.uva_lai_boost * I_UVA / (p.K_uva_lai + I_UVA)

        stress_suppression = 1.0 - Stress / (p.K_stress + Stress + 1e-9)
        sla_boost = sla_boost * stress_suppression
        lai_boost = lai_boost * stress_suppression

        if dLAI_dt_base > 0:
            dLAI_dt_base = dLAI_dt_base * (1.0 + lai_boost)
        else:
            dLAI_dt_base = dLAI_dt_base * (1.0 - lai_boost * 0.3)

        if dXd_dt_base > 0:
            dXd_dt_base = dXd_dt_base * (1.0 + sla_boost * 0.5)
        else:
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
    # Step 9: Calculate nonlinear damage factor
    # =========================================================================
    # NOTE: Using hours_today (current progress) for progressive damage accumulation
    # This is biologically realistic - damage accumulates over the exposure period
    # The documentation table shows FINAL daily values for reference
    # Calculate scheduled daily_hours for use in other calculations
    daily_hours = uva_hour_off - uva_hour_on if uva_hour_on < uva_hour_off else 24 - uva_hour_on + uva_hour_off
    if not uva_on:
        daily_hours = 0
    # nonlinear_factor based on current exposure progress (hours_today)
    nonlinear_factor = nonlinear_damage_factor(hours_today, p)

    # =========================================================================
    # Step 10: Calculate AOX protection
    # =========================================================================
    aox_protection = p.alpha_aox_protection * AOX / (p.K_aox_protection + AOX + 1e-12)

    # =========================================================================
    # Step 11: Calculate damage rate
    # =========================================================================
    vuln_damage = p.stress_damage_coeff * ROS * vulnerability
    nonlin_damage = p.k_nonlinear_stress * ROS * nonlinear_factor
    base_damage = vuln_damage + nonlin_damage
    protected_damage = base_damage * (1.0 - aox_protection)
    damage_rate = protected_damage + circadian_damage

    # =========================================================================
    # Step 12: Calculate Stress decay
    # =========================================================================
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
    # Step 15: Calculate AOX dynamics
    # =========================================================================
    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress, p, nonlinear_factor)

    base_synthesis = p.base_aox_rate_light if is_day else p.base_aox_rate_dark

    # total_uva_hours: only count completed days + current session progress
    # hours_today is 0 when UVA is off, so this correctly tracks actual irradiation
    total_uva_hours = max(0, days_irradiated - 1) * daily_hours + hours_today

    # Nighttime irradiation efficiency
    is_night_irradiation = (uva_hour_on >= 18) or (uva_hour_off <= 6)
    night_eff = p.night_stress_efficiency if is_night_irradiation else 1.0

    # LAI efficiency
    LAI_stress_efficiency = min(1.0, (LAI / p.LAI_healthy) ** p.n_LAI_eff)

    # Stress-induced synthesis
    stress_induced = p.V_max_aox * Stress / (p.K_stress_aox + Stress + 1e-12) * night_eff * LAI_stress_efficiency

    # UV direct induction
    uv_induced = p.k_uv_aox * total_uva_hours / (p.K_uv_hours + total_uva_hours + 1e-12)

    # Stress inhibition on synthesis
    stress_inhibition_synth = p.max_stress_inhib * (Stress ** p.n_stress_inhib) / (p.K_stress_inhib ** p.n_stress_inhib + Stress ** p.n_stress_inhib + 1e-9)
    stress_efficiency = 1.0 - stress_inhibition_synth

    # Water inhibition
    water_efficiency = calculate_water_aox_efficiency(dw_fw_ratio, p)

    # Nonlinear factor efficiency
    daily_nonlin_factor = nonlinear_damage_factor(daily_hours, p)
    nonlin_aox_efficiency = calculate_nonlin_aox_efficiency(daily_nonlin_factor, p)

    # Adaptation factor
    adaptation_factor = p.K_adapt_days / (p.K_adapt_days + days_irradiated)

    # Total AOX synthesis rate
    aox_synthesis_rate = LAI * (base_synthesis + uv_induced + stress_induced * adaptation_factor * nonlin_aox_efficiency) * stress_efficiency * water_efficiency

    # AOX degradation
    natural_degradation = p.k_aox_deg * AOX

    # AOX consumption by ROS
    daily_nonlin = nonlinear_damage_factor(daily_hours, p)

    # Consumption amplification for extreme daily hours (softplus activation)
    x_raw = (daily_nonlin - p.cons_amp_center) / p.cons_amp_scale
    x = p.cons_amp_scale * np.log(1.0 + np.exp(np.clip(x_raw, -50, 50)))
    consumption_amp = 1.0 + p.cons_amp_k * (x ** 2) / (p.cons_amp_K ** 2 + x ** 2 + 1e-9)

    ros_consumption = p.k_aox_consumption * consumption_amp * AOX * (ROS ** p.n_ros_consumption) / (p.K_ros_consumption ** p.n_ros_consumption + ROS ** p.n_ros_consumption + 1e-9)

    # Store synthesis rate for carbon competition calculation
    aox_synthesis_rate_base = aox_synthesis_rate

    # =========================================================================
    # Step 16: Carbon Competition (Growth-Defense Trade-off)
    # =========================================================================
    # AOX synthesis competes with growth for carbon resources
    # Key: Only STRESS-INDUCED AOX synthesis competes with growth
    # Base synthesis (constitutive) does not affect growth
    #
    # This reflects the biological reality:
    # - Constitutive defense (base AOX) is part of normal metabolism
    # - Stress-induced defense diverts resources from growth
    #

    # Calculate stress-induced AOX synthesis (the component that competes)
    stress_induced_aox = LAI * stress_induced * adaptation_factor * nonlin_aox_efficiency * stress_efficiency * water_efficiency
    stress_aox_carbon_demand = stress_induced_aox * p.aox_carbon_cost

    # Carbon competition from stress-induced synthesis AND cumulative stress
    # D12 groups have high cumulative stress, should have stronger competition
    #
    # Combined effect:
    # 1. Stress-induced AOX synthesis diverts carbon
    # 2. High cumulative stress indicates sustained defense allocation
    #
    # Literature: Monson et al. (2022) DOI: 10.1111/nph.17773
    aox_carbon_effect = stress_aox_carbon_demand / (p.carbon_competition_K + stress_aox_carbon_demand + 1e-12)

    # Additional competition from cumulative stress (for D12 groups)
    # VL3D12 has avgS~60, L6D12 has avgS~150
    stress_carbon_effect = p.stress_competition_max * Stress / (p.stress_competition_K + Stress + 1e-9)

    carbon_competition_effect = aox_carbon_effect * p.carbon_competition_max + stress_carbon_effect

    # Apply carbon competition penalty to growth
    growth_penalty = 1.0 - carbon_competition_effect
    if dXd_dt > 0:
        dXd_dt = dXd_dt * growth_penalty

    # Also reduce AOX synthesis rate when carbon is limited
    # (partial effect - 20% of growth penalty applies to synthesis)
    aox_synthesis_penalty = 1.0 - 0.20 * carbon_competition_effect
    aox_synthesis_rate = aox_synthesis_rate_base * aox_synthesis_penalty

    # Recalculate dAOX_dt with reduced synthesis
    dAOX_dt = aox_synthesis_rate - natural_degradation - ros_consumption

    # =========================================================================
    # Step 17: Real Carbon Consumption from C_buf
    # =========================================================================
    # AOX synthesis consumes carbon from C_buf
    #
    # Literature:
    # - Vogt (2010) DOI: 10.1093/mp/ssp106: ~20% photosynthate to phenylpropanoids
    # - Gershenzon (1994) DOI: 10.1007/BF02059810: metabolic cost 1.5-3x substrate
    #
    # Max consumption limited by C_buf availability for numerical stability
    aox_carbon_demand = aox_synthesis_rate * p.aox_carbon_cost
    if C_buf > 0:
        max_consumption = C_buf * p.max_cbuf_consumption
        aox_carbon_consumption = min(aox_carbon_demand, max_consumption)
    else:
        aox_carbon_consumption = 0.0

    dCbuf_dt = dCbuf_dt - aox_carbon_consumption

    # =========================================================================
    # Return derivative vector (6 state variables)
    # =========================================================================
    return np.array([dXd_dt, dCbuf_dt, dLAI_dt, dAOX_dt, dStress_dt, dROS_dt])


# ==============================================================================
# Output Conversion Functions
# ==============================================================================

def aox_to_anthocyanin(AOX, p):
    """
    Convert AOX to Anthocyanin
    Anthocyanin is approximately 18% of total antioxidants
    """
    return AOX * p.anthocyanin_fraction


def calculate_anthocyanin_ppm(AOX, FW_total_kg, p):
    """
    Calculate anthocyanin concentration in ppm (mg/kg FW)
    """
    anth_kg = aox_to_anthocyanin(AOX, p)
    anth_ppm = anth_kg / FW_total_kg * 1e6
    return anth_ppm


# ==============================================================================
# Main Program
# ==============================================================================

if __name__ == "__main__":

    # Environment base settings
    ENV_BASE = {
        'light_on_hour': 6,
        'light_off_hour': 22,
        'I_day': 57,
        'T_day': 25,
        'T_night': 18,
        'CO2_day': 1200,
        'CO2_night': 1200,
        'RH_day': 0.70,
        'RH_night': 0.85,
        'plant_density': 36,
    }

    # Simulation settings
    SIMULATION = {
        'days': 21,
        'transplant_offset': 14,
        'initial_fw_g': 10,
    }

    # Training set targets
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
        env = ENV_BASE.copy()
        if treatment in TREATMENT_CONFIGS:
            env.update(TREATMENT_CONFIGS[treatment])
        return env

    p = UVAParams()

    print("=" * 80)
    print("Lettuce Growth and UVA Effect Integrated Model")
    print("(Carbon Competition + AOX/Anthocyanin Framework)")
    print("=" * 80)
    print("\nCore Features:")
    print(f"  1. State variable: AOX (total antioxidants)")
    print(f"  2. Anthocyanin = AOX * {p.anthocyanin_fraction:.0%}")
    print(f"  3. Carbon competition: AOX synthesis consumes C_buf")
    print(f"     - Carbon cost: {p.aox_carbon_cost} kg C per kg AOX")
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
        # Initial AOX (adjust for new framework: AOX = Anth / 0.18)
        Anth_init_ppm = 5.0  # Initial anthocyanin concentration
        Anth_init = Anth_init_ppm * fw_total_init / 1e6
        AOX_init = Anth_init / p.anthocyanin_fraction

        initial_state = [Xd_init, C_buf_init, LAI_init, AOX_init, 0.0, 0.0]

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
            Xd_f, Cbuf_f, LAI_f, AOX_f, _, _ = sol.y[:, -1]

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
            uva_hour_on = env.get('uva_hour_on', 0)
            uva_hour_off = env.get('uva_hour_off', 0)
            hours_daily = uva_hour_off - uva_hour_on if uva_hour_on < uva_hour_off else 24 - uva_hour_on + uva_hour_off
            if not env.get('uva_on', False):
                hours_daily = 0
            nonlin_factor = nonlinear_damage_factor(hours_daily, p)

            # Calculate FW and Anthocyanin
            dw_fw_ratio = calculate_dynamic_dw_fw_ratio(avg_stress, p, nonlin_factor)
            FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
            FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']

            # Convert AOX to Anthocyanin
            Anth_sim = calculate_anthocyanin_ppm(AOX_f, FW_total_kg, p)

            FW_obs = target['FW']
            Anth_obs = target['Anth']

            fw_err = (FW_sim - FW_obs) / FW_obs * 100
            anth_err = (Anth_sim - Anth_obs) / Anth_obs * 100

            fw_errs.append(abs(fw_err))
            anth_errs.append(abs(anth_err))

            s1 = "PASS" if abs(fw_err) < 5 else "FAIL"
            s2 = "PASS" if abs(anth_err) < 5 else "FAIL"

            print(f"{treatment:<8} LAI:{LAI_f:>4.1f} avgS:{avg_stress:>5.1f} "
                  f"FW:{FW_sim:>5.1f}g({fw_err:>+5.1f}%{s1}) "
                  f"Anth:{Anth_sim:>4.0f}({anth_err:>+5.1f}%{s2}) "
                  f"dw/fw:{dw_fw_ratio:.3f} AOX:{AOX_f*1e6:.2f}mg C_buf:{Cbuf_f*1e3:.2f}mg")
        else:
            print(f"{treatment:<8} Simulation failed: {sol.message}")

    print("-" * 80)
    fw_ok = sum(1 for e in fw_errs if e < 5)
    anth_ok = sum(1 for e in anth_errs if e < 5)
    print(f"Pass: FW {fw_ok}/6, Anth {anth_ok}/6, Total {fw_ok + anth_ok}/12")

    # =========================================================================
    # Validation Experiment
    # =========================================================================
    print("\n" + "=" * 80)
    print("Validation Experiment: 3-Day Gradient (Day 32-35)")
    print("=" * 80)

    validation_targets = {
        'CK':      {'FW': 85.14, 'Anth': 413, 'hours': 0},
        'VL3D3':   {'FW': 89.1, 'Anth': 437, 'hours': 3},
        'L6D3':    {'FW': 92.18, 'Anth': 468, 'hours': 6},
        'M9D3':    {'FW': 83.79, 'Anth': 539, 'hours': 9},
        'H12D3':   {'FW': 62.2, 'Anth': 657, 'hours': 12},
        'VH15D3':  {'FW': 51.2, 'Anth': 578, 'hours': 15},
    }

    val_fw_errs = []
    val_anth_errs = []

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
        AOX_init = Anth_init / p.anthocyanin_fraction

        initial_state = [Xd_init, C_buf_init, LAI_init, AOX_init, 0.0, 0.0]

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
            Xd_f, Cbuf_f, LAI_f, AOX_f, _, _ = sol.y[:, -1]

            uva_start = env.get('uva_start_day', 35) * 86400
            stress_sum = 0.0
            stress_count = 0
            for i in range(len(sol.t)):
                if sol.t[i] >= uva_start:
                    stress_sum += sol.y[4, i]
                    stress_count += 1
            avg_stress = stress_sum / max(1, stress_count)

            hours_daily = hours
            nonlin_factor = nonlinear_damage_factor(hours_daily, p)

            dw_fw_ratio = calculate_dynamic_dw_fw_ratio(avg_stress, p, nonlin_factor)
            FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
            FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
            Anth_sim = calculate_anthocyanin_ppm(AOX_f, FW_total_kg, p)

            FW_obs = target['FW']
            Anth_obs = target['Anth']

            fw_err = (FW_sim - FW_obs) / FW_obs * 100
            anth_err = (Anth_sim - Anth_obs) / Anth_obs * 100

            val_fw_errs.append(abs(fw_err))
            val_anth_errs.append(abs(anth_err))

            fw_s = "PASS" if abs(fw_err) < 5 else ("WARN" if abs(fw_err) < 10 else "FAIL")
            anth_s = "PASS" if abs(anth_err) < 5 else ("WARN" if abs(anth_err) < 10 else "FAIL")

            print(f"{name:<8} {hours:>2}h/day LAI:{LAI_f:>4.1f} avgS:{avg_stress:>6.0f} "
                  f"FW:{FW_sim:>5.1f}g({fw_err:>+5.1f}%{fw_s}) "
                  f"Anth:{Anth_sim:>4.0f}({anth_err:>+5.1f}%{anth_s}) "
                  f"dw/fw:{dw_fw_ratio:.3f}")
        else:
            print(f"{name:<8} Simulation failed: {sol.message}")

    print("-" * 80)
    val_fw_ok5 = sum(1 for e in val_fw_errs if e < 5)
    val_fw_ok10 = sum(1 for e in val_fw_errs if e < 10)
    val_anth_ok5 = sum(1 for e in val_anth_errs if e < 5)
    val_anth_ok10 = sum(1 for e in val_anth_errs if e < 10)
    print(f"Validation FW: <5%: {val_fw_ok5}/6, <10%: {val_fw_ok10}/6")
    print(f"Validation Anth: <5%: {val_anth_ok5}/6, <10%: {val_anth_ok10}/6")
