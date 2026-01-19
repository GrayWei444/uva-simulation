
Supplementary Materials: Mechanistic Model Parameters and Implementation Details

Supplementary Table S1. Complete parameter set for the mechanistic ODE model (v10.39)
This table lists all calibrated parameters, consistent with the repository implementation (simulate_uva_model_v10.py; ALL_PARAMS).
Category	Parameter	Description	Value	Unit
A	c_alpha	Photosynthetic efficiency scaling factor	0.54	-
B	uva_sla_enhancement	Maximum UV-A enhancement on SLA	5.0	-
B	K_uva_sla	Half-saturation UV-A intensity for SLA enhancement	7.5	W m−2
B	uva_lai_boost	Maximum UV-A boost on LAI growth	1.70	-
B	K_uva_lai	Half-saturation UV-A intensity for LAI boost	7.5	W m−2
C	k_ros_production	ROS production coefficient	0.010	ROS (W m−2 s)−1
C	k_ros_clearance	ROS clearance coefficient	5×10−4	s−1
D	stress_damage_coeff	Baseline damage coefficient	1.6×10−7	Stress (ROS·s)−1
D	A_vulnerability	Amplitude of LAI-dependent vulnerability	8.5×10^7	-
D	k_vulnerability	Decay of vulnerability with LAI	2.0	-
D	gompertz_max_factor	Maximum damage amplification (upper bound)	250.0	-
D	gompertz_threshold	Gompertz inflection (collapse time)	10.5	h
D	gompertz_steepness	Gompertz steepness	0.5	h−1
D	alpha_anth_protection	Maximum anthocyanin protection efficiency	0.5	-
D	K_anth_protection	Half-saturation for anthocyanin protection	5×10−6	kg m−2
D	k_circadian	Circadian damage coefficient (night irradiation)	3.0×10−6	-
D	n_circadian	Exponent for circadian damage	2.0	-
E	k_repair_lai	LAI repair coefficient	5.0×10−8	(m2 m−2·s)−1
E	k_stress_decay	Stress decay rate (half-life ≈ 9 h)	2.14×10−5	s−1
F	stress_photosynthesis_inhibition	Maximum photosynthesis inhibition by stress	0.85	-
F	stress_lai_inhibition	Maximum LAI-growth inhibition by stress	0.80	-
F	K_stress	Half-saturation stress level for growth inhibition	50.0	-
G	base_anth_rate_light	Baseline anthocyanin synthesis rate (light)	6.35×10−10	kg (m2 LAI·s)−1
G	base_anth_rate_dark	Baseline anthocyanin synthesis rate (dark)	3.18×10−10	kg (m2 LAI·s)−1
G	V_max_anth	Maximum stress-induced synthesis rate	2.75×10−9	kg (m2 LAI·s)−1
G	K_stress_anth	Half-saturation for stress-induced synthesis	100.0	-
G	k_deg	Anthocyanin degradation rate	3.02×10−6	s−1
G	k_anth_consumption	Anthocyanin consumption coefficient (ROS usage)	1.0×10−6	s−1
H	water_anth_threshold	LDMC threshold triggering water inhibition (DW/FW)	0.055	-
H	water_anth_K	Half-saturation for water inhibition	0.020	-
H	water_anth_max_inhib	Maximum water inhibition on synthesis	0.50	-
I	K_stress_inhib	Half-saturation for stress inhibition on synthesis	150.0	-
I	n_stress_inhib	Hill exponent for stress inhibition	2.0	-
I	max_stress_inhib	Maximum stress inhibition on synthesis	0.80	-
J	hill_K	Half-effect constant for efficiency inhibition	800.0	-
J	hill_n	Hill exponent for efficiency inhibition	1.5	-
K	K_adapt_days	Half-saturation days for adaptation effect	4.0	days
L	dw_fw_ratio_base	Baseline DW/FW ratio (healthy plants)	0.05	-
L	ldmc_stress_sensitivity	Stress sensitivity of LDMC	0.45	-
L	K_ldmc	Half-saturation for LDMC effect	1400.0	-
L	dw_fw_ratio_max	Maximum DW/FW ratio (severe stress)	0.080	-
M	PPFD	Photosynthetic photon flux density	130	µmol m−2 s−1
M	I_base	Equivalent baseline shortwave irradiance (LED)	57	W m−2
M	I_UVA	UV-A irradiance (365 nm)	11	W m−2
M	plant_density	Planting density	36	plants m−2

Supplementary Table S2. Initial conditions and output conversion
The simulation starts at 14 DAS (transplanting) and runs to 35 DAS. All rates are in seconds. Initial fresh weight at transplanting is assumed to be FW_init = 10 g plant−1, with LDMC_base = 0.05.
Initial conditions:
X_d(0) = FW_init × LDMC_base × plant_density / 1000 = 0.018 kg m−2
C_buf(0) = 0.1 × X_d(0) = 0.0018 kg m−2
LAI(0) = 0.40 m2 m−2
Anth(0) = 1.8×10−6 kg m−2
Stress(0) = 0
ROS(0) = 0
Output conversion:
FW (g plant−1) = [X_d* / (plant_density × LDMC)] × 1000
Anth_ppm (mg kg−1) = (Anth* / FW_total) × 10^6, where FW_total = FW × plant_density / 1000 (kg m−2)
Code availability
The full implementation, configuration files and validation scripts are provided in the public GitHub repository (URL to be inserted by the authors). We recommend citing a tagged release and archiving it via Zenodo to obtain a DOI for peer-review reproducibility.
