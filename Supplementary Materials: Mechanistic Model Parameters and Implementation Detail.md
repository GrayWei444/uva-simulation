
Supplementary Materials: Mechanistic Model Parameters and Implementation Details

Supplementary Table S1. Complete parameter set for the mechanistic ODE model (Carbon Competition + AOX Framework)
This table lists all calibrated parameters, consistent with the repository implementation (lettuce_uva_model.py; ALL_PARAMS).

Category	Parameter	Description	Value	Unit
A	c_alpha	Photosynthetic efficiency scaling factor	0.54	-
B	uva_sla_enhancement	Maximum UV-A enhancement on SLA	5.0	-
B	K_uva_sla	Half-saturation UV-A intensity for SLA enhancement	7.5	W m^-2
B	uva_lai_boost	Maximum UV-A boost on LAI growth	1.70	-
B	K_uva_lai	Half-saturation UV-A intensity for LAI boost	7.5	W m^-2
C	k_ros_production	ROS production coefficient	0.010	ROS (W m^-2 s)^-1
C	k_ros_clearance	ROS clearance coefficient	5x10^-4	s^-1
D	stress_damage_coeff	Baseline damage coefficient	1.6x10^-7	Stress (ROS s)^-1
D	A_vulnerability	Amplitude of LAI-dependent vulnerability	8.5x10^7	-
D	k_vulnerability	Decay of vulnerability with LAI	2.0	-
D	gompertz_max_factor	Maximum damage amplification (upper bound)	250.0	-
D	gompertz_threshold	Gompertz inflection (collapse time)	10.5	h
D	gompertz_steepness	Gompertz steepness	0.5	h^-1
D	alpha_aox_protection	Maximum AOX protection efficiency	0.5	-
D	K_aox_protection	Half-saturation for AOX protection	2.78x10^-5	kg m^-2
D	k_circadian	Circadian damage coefficient (night irradiation)	3.0x10^-6	-
D	n_circadian	Exponent for circadian damage	2.0	-
D	k_nonlinear_stress	Nonlinear stress coefficient	5.0x10^-6	-
E	k_stress_decay	Stress decay rate (half-life ~ 9 h)	2.14x10^-5	s^-1
F	stress_photosynthesis_inhibition	Maximum photosynthesis inhibition by stress	0.85	-
F	stress_lai_inhibition	Maximum LAI-growth inhibition by stress	0.80	-
F	K_stress	Half-saturation stress level for growth inhibition	50.0	-
G	base_aox_rate_light	Baseline AOX synthesis rate (light)	3.53x10^-9	kg (m^2 LAI s)^-1
G	base_aox_rate_dark	Baseline AOX synthesis rate (dark)	1.77x10^-9	kg (m^2 LAI s)^-1
G	V_max_aox	Maximum stress-induced AOX synthesis rate	1.45x10^-8	kg (m^2 LAI s)^-1
G	K_stress_aox	Half-saturation for stress-induced synthesis	100.0	-
G	k_aox_deg	AOX degradation rate	3.02x10^-6	s^-1
G	k_aox_consumption	AOX consumption coefficient (ROS scavenging)	1.8x10^-7	s^-1
G	k_uv_aox	UV-induced AOX synthesis coefficient	7.78x10^-11	kg m^-2 s^-1
G	K_uv_hours	Half-saturation for UV hours effect	30.0	h
H	water_aox_threshold	LDMC threshold triggering water inhibition (DW/FW)	0.055	-
H	water_aox_K	Half-saturation for water inhibition	0.020	-
H	water_aox_max_inhib	Maximum water inhibition on synthesis	0.50	-
H	water_n	Hill coefficient for water inhibition	2.0	-
I	K_stress_inhib	Half-saturation for stress inhibition on synthesis	150.0	-
I	n_stress_inhib	Hill exponent for stress inhibition	2.0	-
I	max_stress_inhib	Maximum stress inhibition on synthesis	0.80	-
J	K_nonlin_aox	Half-effect constant for nonlinear efficiency	800.0	-
J	n_nonlin_aox	Hill exponent for nonlinear efficiency	1.5	-
K	K_adapt_days	Half-saturation days for adaptation effect	4.0	days
L	dw_fw_ratio_base	Baseline DW/FW ratio (healthy plants)	0.05	-
L	ldmc_stress_sensitivity	Stress sensitivity of LDMC	0.45	-
L	K_ldmc	Half-saturation for LDMC effect	1400.0	-
L	dw_fw_ratio_max	Maximum DW/FW ratio (severe stress)	0.080	-
L	acute_center	Softplus center for acute LDMC	50.0	-
L	acute_scale	Softplus scale for acute LDMC	10.0	-
L	acute_k	Acute LDMC coefficient	9.0	-
L	acute_K	Acute LDMC half-saturation	120.0	-
L	acute_n	Acute LDMC Hill coefficient	2.0	-
M	LAI_healthy	Reference LAI for healthy plant	9.0	m^2 m^-2
M	n_LAI_eff	Hill coefficient for LAI efficiency	2.0	-
M	night_stress_efficiency	Night irradiation stress efficiency	0.4	-
N	K_ros_consumption	ROS half-saturation for AOX consumption	500.0	-
N	n_ros_consumption	Hill coefficient for ROS consumption	2.0	-
N	cons_amp_center	Softplus center for consumption amplification	200.0	-
N	cons_amp_scale	Softplus scale for consumption amplification	15.0	-
N	cons_amp_k	Amplification coefficient	12.0	-
N	cons_amp_K	Amplification half-saturation	20.0	-
O	aox_carbon_cost	Carbon cost per AOX synthesized	1.0	kg C / kg AOX
O	carbon_competition_K	Half-saturation for AOX carbon effect	1x10^-8	kg m^-2 s^-1
O	stress_competition_K	Half-saturation for stress competition	21.0	-
O	stress_competition_max	Maximum stress-based competition	0.225	-
O	carbon_competition_max	Maximum AOX-based competition	0.30	-
O	max_cbuf_consumption	Max fraction of C_buf consumed per timestep	0.10	-
P	anthocyanin_fraction	Fraction of AOX that is anthocyanin	0.18	-
Q	PPFD	Photosynthetic photon flux density	130	umol m^-2 s^-1
Q	I_base	Equivalent baseline shortwave irradiance (LED)	57	W m^-2
Q	I_UVA	UV-A irradiance (365 nm)	11	W m^-2
Q	plant_density	Planting density	36	plants m^-2

Supplementary Table S2. Initial conditions and output conversion
The simulation starts at 14 DAS (transplanting) and runs to 35 DAS. All rates are in seconds. Initial fresh weight at transplanting is assumed to be FW_init = 10 g plant^-1, with LDMC_base = 0.05.

Initial conditions:
X_d(0) = FW_init x LDMC_base x plant_density / 1000 = 0.018 kg m^-2
C_buf(0) = 0.1 x X_d(0) = 0.0018 kg m^-2
LAI(0) = 2.0 m^2 m^-2 (calculated as FW_init × LDMC_base / 0.01 × 0.04)
AOX(0) = 1.0x10^-5 kg m^-2 (equivalent to ~5 ppm anthocyanin)
Stress(0) = 0
ROS(0) = 0

Output conversion:
FW (g plant^-1) = [X_d* / (plant_density x LDMC)] x 1000
Anth_ppm (mg kg^-1) = (AOX* x 0.18 / FW_total) x 10^6, where FW_total = FW x plant_density / 1000 (kg m^-2)

Supplementary Table S3. Carbon competition mechanism parameters and formulas
The carbon competition mechanism implements the Growth-Differentiation Balance Hypothesis (Herms & Mattson, 1992).

Formula: Carbon competition penalty
stress_carbon_effect = stress_competition_max x Stress / (stress_competition_K + Stress)
                     = 0.225 x Stress / (21 + Stress)

aox_carbon_effect = stress_aox_carbon_demand / (carbon_competition_K + stress_aox_carbon_demand)

carbon_competition_effect = aox_carbon_effect x carbon_competition_max + stress_carbon_effect
                          = aox_carbon_effect x 0.30 + stress_carbon_effect

growth_penalty = 1.0 - carbon_competition_effect
dX_d/dt = dX_d/dt x growth_penalty

Formula: Real C_buf consumption
max_consumption = C_buf x max_cbuf_consumption = C_buf x 0.10
aox_carbon_consumption = min(aox_synthesis_rate x aox_carbon_cost, max_consumption)
dC_buf/dt = dC_buf/dt - aox_carbon_consumption

Note: aox_synthesis_rate includes all components (baseline + UV-induced + stress-induced), reflecting that all AOX biosynthesis requires carbon substrates from the phenylpropanoid pathway.

Formula: AOX synthesis penalty (feedback inhibition)
aox_synthesis_penalty = 1.0 - 0.20 x carbon_competition_effect
S_AOX = S_AOX_base x aox_synthesis_penalty

This feedback reflects reduced synthetic capacity when carbon is limiting, consistent with metabolic regulation of secondary metabolism under resource constraints.

Supplementary Table S4. Literature support for carbon competition
Mechanism	Reference	DOI
Growth-Differentiation Balance Hypothesis	Herms & Mattson (1992)	10.1086/285343
Coordinated resource allocation	Monson et al. (2022)	10.1111/nph.17773
Phenylpropanoid carbon flux (~20% photosynthate)	Vogt (2010)	10.1093/mp/ssp106
Metabolic costs of secondary metabolites	Gershenzon (1994)	10.1007/BF02059810

Supplementary Table S5. Nonlinear damage factor (Gompertz function)
Formula: factor = 1 + gompertz_max_factor x exp(-exp(-gompertz_steepness x (hours - gompertz_threshold)))
       = 1 + 250 x exp(-exp(-0.5 x (hours - 10.5)))

Daily Hours	Factor
3h	1.0
6h	1.0
9h	31.1
12h	156.9
15h	226.7

Note: The nonlinear damage factor is calculated using hours_today (current exposure progress within the day) for progressive damage accumulation. The table shows FINAL daily values for reference.

Supplementary Table S6. Stress dynamics implementation detail
The stress accumulation in Equation (7) is implemented as follows:

vuln_damage = stress_damage_coeff × ROS × vulnerability(LAI)
nonlin_damage = k_nonlinear_stress × ROS × Nonlinear_Factor
base_damage = vuln_damage + nonlin_damage
protected_damage = base_damage × (1 − aox_protection)
total_damage = protected_damage + circadian_damage
dStress/dt = total_damage − k_stress_decay × Stress

where:
- vuln_damage: LAI-dependent vulnerability effect
- nonlin_damage: Gompertz-modulated nonlinear damage amplification
- base_damage: combined damage before protection
- protected_damage: base damage reduced by AOX protection (protection applies to both vulnerability and nonlinear components)
- circadian_damage: additional damage for night irradiation (not protected by AOX)

This structure applies AOX protection uniformly to the combined damage from both vulnerability and nonlinear amplification, while circadian damage remains unprotected to reflect that circadian disruption operates through distinct physiological pathways.

Supplementary Table S7. Sun model parameters (handled internally)
The main text equations include parameters K_carbon (Equation 4) and k_senescence (Equation 5) for conceptual clarity. These parameters are handled internally by the base Sun model framework (lettuce_uva_carbon_complete_model.py) rather than being explicitly calibrated in the UV-A extension layer.

Parameter	Description	Handling
K_carbon	Half-saturation constant for carbon-limited growth	Embedded in Sun model's photosynthesis and respiration calculations; the carbon buffer dynamics implicitly regulate growth through the C_buf state variable
k_senescence	LAI senescence/turnover rate	Incorporated in Sun model's leaf area dynamics; accounts for natural leaf turnover and age-dependent senescence

The conceptual equations in the main text (Equations 4-5) present a simplified view for clarity. The actual implementation inherits the full Sun model framework, which provides mechanistically-grounded carbon and LAI dynamics calibrated for lettuce growth under controlled environment conditions.

Code availability
The full implementation, configuration files and validation scripts are provided in the public GitHub repository (URL to be inserted by the authors). The main model file is lettuce_uva_model.py. We recommend citing a tagged release and archiving it via Zenodo to obtain a DOI for peer-review reproducibility.

Key model files:
- lettuce_uva_model.py: Main simulation code with carbon competition mechanism
- lettuce_uva_carbon_complete_model.py: Base Sun model (dependency)
- parameters.md: Complete parameter documentation
