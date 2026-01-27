A Two-Stage Screening-to-Optimization Approach with Mechanistic Model Analysis: Enhancing Anthocyanin in Lettuce Without Yield Loss
Zhi-Hao Wei¹, Wei Fang²*, Zhen-Kang Huang³*
1 Department of Bio-Industrial Mechatronics Engineering, National Taiwan University, Taipei 10617, Taiwan
2 Professor, Department of Bio-Industrial Mechatronics Engineering, National Taiwan University, Taipei 10617, Taiwan
3 Professor, Department of Bio-Industrial Mechatronics Engineering, National Taiwan University, Taipei 10617, Taiwan
* Correspondence: (to be filled)

Abstract
Enhancing anthocyanin accumulation in red-leaf lettuce grown in plant factories often incurs yield penalties. Here we propose a two-stage screening-to-optimization framework integrated with mechanistic modeling to resolve this tradeoff. In Stage 1, comparative experiments confirmed that UV-A is more compatible with growth and pigmentation than UV-B, and identified 'Lollo Rosso' as a highly responsive cultivar. In Stage 2, optimization experiments showed that L6D6 (6 h day⁻¹ for 6 days) increased the total anthocyanin per plant by 15.4% while maintaining fresh weight. Motivated by observed nonlinear phenomena including biomass overcompensation, circadian disruption under night irradiation, and ontogeny-dependent vulnerability, we developed a six-state ordinary differential equation (ODE) model that integrates reactive oxygen species (ROS) dynamics with stress damage–repair processes. A key innovation is the explicit representation of carbon competition between growth and antioxidant defense, where AOX synthesis consumes carbon from the buffer pool, creating a physiologically meaningful growth–defense tradeoff supported by the Growth-Differentiation Balance Hypothesis [1]. The model achieved high accuracy in an independent validation set that included extreme doses (errors ≤ 10.6%, with 11 of 12 metrics < 10%), supporting the physiological necessity of the introduced mechanisms. Global optimization based on the calibrated model predicted that 9 h day⁻¹ for 4 days is the theoretical optimum, potentially increasing total anthocyanin by 38.3% with minimal fresh weight reduction (−2.4%), substantially outperforming the best experimental treatment. This quantitative mechanistic framework provides a scientific basis for designing precise stress-light recipes in controlled-environment agriculture.
Keywords: UV-A radiation; red-leaf lettuce; anthocyanin biosynthesis; controlled-environment agriculture; LED supplementation; screening-to-optimization methodology; mechanistic modeling; ordinary differential equations; carbon competition; growth-defense tradeoff

1. Introduction
Lettuce (Lactuca sativa L.) is widely consumed, and red-leaf lettuce is distinguished by its high anthocyanin content. Dietary supplementation with red-pigmented lettuce has been reported to improve lipid profiles and antioxidant status in a mouse model, suggesting potential cardiovascular benefits [2]. Light is a key signal regulating plant development and metabolism across the life cycle [3]. Plants perceive light cues through photoreceptors such as phytochromes, cryptochromes and UVR8 [4], which activate transcription factors including HY5 to induce the expression of anthocyanin biosynthetic genes [5,6]. In subtropical plant factories, however, artificial lighting typically lacks ultraviolet bands, which may limit the normal pigmentation of red-leaf lettuce.
Pigmentation and the accumulation of secondary metabolites are largely defensive responses to environmental stresses such as UV radiation [7,8]. Although UV supplementation can promote anthocyanin accumulation [9], it is a double-edged sword and may simultaneously suppress growth [10]. This growth–defense tradeoff has been extensively studied in ecological literature. The Growth-Differentiation Balance Hypothesis [1] posits that plants allocate limited carbon resources between growth and defense, creating an inherent tradeoff. Recent work has further elucidated the coordinated resource allocation mechanisms underlying this tradeoff [11]. Understanding the carbon costs of defense compound synthesis [12,13] is crucial for optimizing light recipes that balance yield and quality.
Distinguishing band-specific effects is crucial. Weiland et al. [14] reported that UV-B can increase flavonoids but its higher photon energy easily leads to growth inhibition. In contrast, UV-A is considered a milder regulatory tool: Li and Kubota [15] showed that UV-A increased lettuce phytochemicals, and Chen et al. [16] further observed that UV-A can even increase biomass. Responses can be cultivar specific [17].
Regarding crop growth models, existing lettuce models [18,19] focus primarily on dry matter accumulation under non-stress conditions and do not explicitly represent the dynamic impact of UV stress on yield and quality. There is still a lack of a quantitative mechanistic model that captures UV-induced stress, damage–repair dynamics, carbon competition between growth and defense, and their interactive effects on biomass and anthocyanin. Therefore, this study aims to develop and validate a two-stage screening-to-optimization approach integrating experiments and mechanistic modeling: (1) Screening stage: compare UV-A and UV-B across cultivars to identify a safer band and a responsive cultivar; (2) Optimization stage: test multiple UV-A recipes for the selected cultivar to maximize anthocyanin without yield loss; and (3) Mechanistic integration: develop a continuous, physiologically interpretable model that integrates stress damage–repair dynamics, carbon competition, ontogeny-dependent vulnerability and UV-A morphological effects, and calibrate parameters using data-driven fitting to enable quantitative prediction.

2. Results
2.1. Experiment 1: cultivar screening under UV-A and UV-B
UV-B caused visible necrotic lesions and significant fresh-weight reduction in all cultivars, indicating severe photodamage. In contrast, UV-A induced red pigmentation without visible injury, and fresh weight was not statistically different from controls across cultivars. Representative phenotypes are shown in Figure 1, and fresh-weight responses are summarized in Figure 2.
Figure 1. Morphology of four lettuce cultivars at 35 DAS under control (no UV), UV-A and UV-B. (a) Whole-plant representative images; (b) close-up of the central leaves. UV-B caused clear necrosis, whereas UV-A produced uniform pigmentation without injury.
Figure 2. Effects of UV-A and UV-B on shoot fresh weight for four lettuce cultivars. Data are mean ± SD (n = 5). Different letters indicate significant differences (p < 0.05).
2.1.1. Color parameters and anthocyanin content
Both UV-A and UV-B significantly decreased L* (lightness) and increased a* (redness). Although UV-B produced the highest anthocyanin concentration, it was accompanied by strong biomass penalties. Considering the growth-quality tradeoff, 'Lollo Rosso' showed robust growth and strong pigmentation response and was selected for optimization. Leaf color parameters are reported in Figure 3, and anthocyanin concentration is shown in Figure 4.
Figure 3. Effects of UV-A and UV-B on leaf color parameters (L*, a*, b*) across cultivars. UV treatments increased a* values indicating enhanced redness.
Figure 4. Effects of UV-A and UV-B on anthocyanin concentration across cultivars. UV-B induced the highest concentration but with strong growth inhibition.
2.2. Experiment 2: optimization of UV-A recipes for 'Lollo Rosso'
2.2.1. Morphology and fresh weight
All UV-A treatments enhanced red coloration. Fresh weight was maintained or slightly increased only in UV-A-L6D6 (91.4 g) compared with CK (87.0 g). Other recipes (H12D3, VL3D12, L6D12) significantly reduced fresh weight, indicating that higher daily dose or longer cumulative exposure is detrimental to growth. Representative phenotypes are shown in Figure 5, and fresh-weight responses are summarized in Figure 6.
Figure 5. Morphology of 'Lollo Rosso' at 35 DAS under different UV-A optimization recipes. L6D6 produced compact plants with strong pigmentation.
Figure 6. Effects of different UV-A recipes on shoot fresh weight of 'Lollo Rosso'. Data are mean ± SD (n = 5). Different letters indicate significant differences (p < 0.05).
2.2.2. Color parameters and anthocyanin
All UV-A treatments significantly increased a* values. H12D3 produced the highest anthocyanin concentration (651 ppm), but due to low biomass its per-plant total anthocyanin was not optimal. UV-A-L6D6 combined higher biomass with elevated concentration (494 ppm), resulting in the highest anthocyanin amount per plant among tested recipes. Leaf color parameters are shown in Figure 7, and anthocyanin concentration and per-plant anthocyanin are shown in Figure 8.
Figure 7. Effects of UV-A recipes on leaf color parameters (L*, a*, b*) of 'Lollo Rosso'.
Figure 8. Effects of UV-A recipes on anthocyanin of 'Lollo Rosso': (a) concentration (ppm); (b) total anthocyanin per plant (mg), calculated as concentration × fresh weight.
Summary (anthocyanin in ppm = mg kg⁻¹):
CK: FW 87.0 g, Anth 433 ppm
UV-A-L6D6: FW 91.4 g, Anth 494 ppm (best balance)
UV-A-H12D3: FW 60.6 g, Anth 651 ppm (high-stress cost)
2.3. Mechanistic modeling results
2.3.1. Training-set calibration
The model fit the six training treatments with low errors (Table 1). Mean absolute errors were 2.3% for fresh weight (max 5.6%) and 2.6% for anthocyanin (max 5.6%), meeting the <5% target across 10 of 12 metrics. The carbon competition mechanism successfully captured the growth–defense tradeoff, with VL3D12 and L6D12 showing elevated stress levels and corresponding growth penalties.
Treatment	FW obs (g)	FW pred (g)	FW error	Anth obs (ppm)	Anth pred (ppm)	Anth error	Avg Stress
CK	87.0	87.9	+1.0%	433	440	+1.6%	0.0
L6D6	91.4	90.2	-1.3%	494	487	-1.4%	7.2
L6D6-N	80.8	83.5	+3.4%	493	481	-2.5%	10.5
VL3D12	67.0	70.7	+5.6%	482	509	+5.6%	61.0
L6D12	60.4	61.5	+1.8%	518	533	+2.9%	145.5
H12D3	60.6	60.9	+0.6%	651	655	+0.6%	262.0

Table 1. Training-set calibration results (Carbon Competition Model).

2.3.2. Independent validation (3-day gradient)
An independent 3-day gradient dataset was used to test generalization (Table 2). The validation control (CK_val) was from a different experimental batch, leading to slight differences from the training control. The model achieved <10% errors for fresh weight across all validation treatments and <10% errors for anthocyanin in 5 of 6 treatments. Fresh-weight and anthocyanin response curves are shown in Figures 9 and 10, and prediction–observation parity is shown in Figure 11.
Treatment	Hours	FW obs (g)	FW pred (g)	FW error	Anth obs (ppm)	Anth pred (ppm)	Anth error
CK_val	0	85.2	87.9	+3.2%	413	440	+6.5%
VL3D3	3	89.0	88.5	-0.7%	437	460	+5.4%
L6D3	6	92.2	89.1	-3.4%	468	480	+2.5%
M9D3	9	83.8	85.9	+2.5%	539	596	+10.6%
H12D3_val	12	62.2	60.9	-2.0%	657	655	-0.3%
VH15D3	15	51.3	51.2	+0.0%	578	605	+4.8%

Table 2. Validation-set prediction results (Carbon Competition Model).

The model reproduced the non-monotonic anthocyanin response: 12 h day⁻¹ produced the peak concentration (657 ppm) while 15 h day⁻¹ declined (578 ppm), consistent with hormesis [20,21].
2.3.3. Mechanistic diagnostics
Model diagnostics indicated that cumulative stress (∫Stress dt) better reflects total stress burden than endpoint stress. Early irradiation (Day 23; VL3D12 and L6D12) incurred high vulnerability at low LAI, causing large cumulative stress despite partial recovery later. Night irradiation increased stress relative to day irradiation, consistent with circadian regulation of ROS homeostasis [22–24].
The carbon competition mechanism explains why VL3D12 (low stress) and L6D12 (high stress) both show FW reduction: VL3D12's stress-based carbon competition (Stress=61) creates a growth penalty of ~17%, while L6D12's higher stress (Stress=145.5) results in ~19% penalty. This demonstrates that the tradeoff operates continuously across stress levels, not just at extremes.
Mechanistic functions and diagnostics are summarized in Figure 12 (LAI vulnerability), Figure 13 (nonlinear damage amplification), Figure 14 (training-set parity plots), Figure 15 (stress time series), Figure 9 (validation fresh-weight response), Figure 10 (validation anthocyanin response), Figure 11 (validation parity plots), Figure 16 (Hill-type inhibition), Figure S1 (system dynamics block diagram with carbon competition), and Figure S2 (hormesis response surface).
Treatment	Pattern	Total hours	Avg Stress	Carbon Competition	Dominant mechanisms
CK	none	0	0.0	0%	control
L6D6	6 h×6 d	36	7.2	~7%	low daytime dose, minimal penalty
L6D6-N	6 h×6 d (night)	36	10.5	~10%	circadian disruption
VL3D12	3 h×12 d	36	61.0	~17%	early-stage vulnerability, carbon competition
L6D12	6 h×12 d	72	145.5	~19%	cumulative stress + carbon competition
H12D3	12 h×3 d	36	262.0	~20%	nonlinear damage + acute LDMC

Table 3. Stress diagnostics, carbon competition effects, and dominant mechanisms.

Figure 9. Validation set: fresh weight response vs daily UV-A hours (0–15 h day⁻¹).
Figure 10. Validation set: anthocyanin response vs daily UV-A hours (0–15 h day⁻¹), showing hormesis at the highest dose.
Figure 11. Validation set: model prediction vs observation with ±5% and ±10% error bands.
Figure 12. LAI vulnerability function v(LAI) = A·exp(−k·LAI) + 1 and the LAI at treatment onset for each recipe.
Figure 13. Nonlinear damage amplification as a function of daily irradiation hours modeled by a Gompertz function.
Figure 14. Training set: model prediction vs observation for (a) fresh weight and (b) anthocyanin. The diagonal is the 1:1 line; shaded bands indicate ±5% error.
Figure 15. Stress time series for each treatment from Day 14 to Day 35. Vertical dashed lines indicate UV-A onset times.
Figure 16. Monotonically decreasing Hill-type inhibition of AOX synthesis efficiency vs nonlinear damage factor (K = 800, n = 1.5).
Figure 17. Sensitivity analysis plots: (a) parameter sensitivity curves; (b) sensitivity heatmap across treatments and outputs.
Figure 18. Optimization heatmaps: (a) fresh-weight change vs daily hours and total days; (b) total anthocyanin change vs daily hours and total days. The best strategy (9 h × 4 d) maximizes total anthocyanin while maintaining FW within −5% of control.
Figure S1. System dynamics block diagram of the six-state model showing carbon competition between growth (X_d) and defense (AOX), with UV-A dual effects through growth promotion (LAI) and stress induction (ROS→Stress).
Figure S2. Hormesis response surface of anthocyanin concentration across daily hours (1–12 h day⁻¹) and total days (2–12 d) predicted by the model.
Figure S3. Carbon competition mechanism visualization showing carbon flow, C_buf dynamics, and AOX accumulation under different treatments.

3. Discussion
Experiment 1 confirmed that UV-B is high risk: it can strongly induce pigmentation but readily causes photodamage and yield loss, whereas UV-A provides a safer window that enhances coloration while preserving biomass, consistent with Weiland et al. [14]. Cultivar differences highlight the genotype-specific nature of light-stress recipes [17].
Experiment 2 demonstrated a classical growth–defense tradeoff [25]. While H12D3 maximized anthocyanin concentration, it reduced biomass by ~30%. From a production perspective, total anthocyanin per plant (concentration × fresh weight) is a more relevant metric. The L6D6 recipe achieved the best balance and even slightly increased fresh weight, consistent with UV-A-associated morphological benefits reported by Chen et al. [16].
The mechanistic model was built by mapping distinct experimental observations to physiologically motivated terms: ontogeny-dependent vulnerability captured the stronger damage when irradiation began early [26,27]; a Gompertz-type nonlinear factor represented the collapse of protective capacity under prolonged daily exposure [28,29]; and a circadian damage term captured the inferior performance of night irradiation, consistent with clock regulation of ROS homeostasis [22–24]. The model's successful prediction of the 15 h day⁻¹ hormesis reversal supports the inclusion of smooth efficiency inhibition under severe stress [20,21].
A key innovation of this model is the explicit representation of carbon competition between growth and antioxidant defense, implementing the Growth-Differentiation Balance Hypothesis [1]. This mechanism provides a physiologically meaningful explanation for the observed growth–defence tradeoff: under UV-A stress, plants allocate carbon from the buffer pool to AOX synthesis, reducing the carbon available for structural growth. The carbon cost of antioxidant synthesis (approximately 1.0 kg C per kg AOX) reflects the substantial metabolic investment in phenylpropanoid biosynthesis, consistent with estimates that up to 20% of photosynthate can be channeled to phenylpropanoids under stress conditions [13]. Recent work on coordinated resource allocation [11] further supports this mechanistic framework, showing that growth–defense tradeoffs arise from fundamental carbon allocation constraints.
Sensitivity analysis (Table 4) indicated that parameters governing early-stage vulnerability and the nonlinear damage threshold are most influential for biomass and stress, whereas AOX-specific parameters primarily affect anthocyanin output, suggesting reasonable decoupling between growth and pigment submodules. The carbon competition parameters (stress_competition_K, stress_competition_max) showed moderate sensitivity for both FW and Anth, reflecting their role in coupling growth and defense. The sensitivity plots are provided in Figure 17.
Finally, global search using the calibrated model suggested that 9 h day⁻¹ for 4 days is the best theoretical strategy that maximizes total anthocyanin while keeping fresh weight within acceptable limits (≥−5%). This optimum lies below the nonlinear damage threshold region and balances sufficient stress induction with limited growth penalty and carbon competition. The predicted hormesis response surface is shown in Figure S2, and the optimization heatmaps are shown in Figure 18. Re-running the released simulation code with default numerical tolerances reproduced the same optimum (9 h × 4 d: FW 84.9 g; Anth 613 ppm).
Parameter	Category	S_FW	S_Stress	S_Anth	Sensitivity
A_vulnerability	LAI vulnerability	-1.28	+13.2	+1.60	Very high
gompertz_threshold	Nonlinear threshold	-0.85	+8.5	+0.92	High
V_max_aox	Stress-induced synthesis	~0	~0	+0.85	High
k_aox_deg	AOX degradation	~0	~0	-0.95	High
base_aox_rate_light	Baseline synthesis (light)	~0	~0	+0.75	High
stress_competition_K	Carbon competition	-0.35	+1.2	-0.25	Medium
K_nonlin_aox	Nonlinear inhibition	~0	~0	+0.30	Medium
k_aox_consumption	AOX consumption	~0	~0	-0.60	Medium
k_stress_decay	Stress decay	+0.45	-2.8	-0.35	Medium

Table 4. Sensitivity analysis of key parameters (Carbon Competition Model).

Rank	Hours per day	Total days	Start day	FW (g)	FW change	Anth (ppm)	Anth change	Total anth change
1	9	4	Day 31	84.9	−2.4%	613	+41.6%	+38.3%
2	9	5	Day 30	84.0	−3.5%	620	+43.1%	+38.1%
3	9	6	Day 29	82.8	−4.9%	620	+43.3%	+36.3%
4	10	2	Day 33	84.2	−3.2%	610	+40.9%	+36.3%
5	9	3	Day 32	85.9	−1.3%	596	+37.7%	+36.0%

Table 5. Top 5 safe UV-A strategies predicted by global search (FW ≥ −5% vs control). Percent changes are relative to the experimental control at harvest (CK: FW = 87.0 g; Anth = 433 ppm).

4. Materials and Methods
4.1. Cultivation environment and equipment
All experiments were conducted in controlled-environment growth chambers at the Department of Bio-Industrial Mechatronics Engineering, National Taiwan University. The hydroponic system consisted of four stacked layers sharing a recirculating nutrient solution to ensure consistent nutrition. Environmental setpoints were: day/night temperature 25/18 ± 1 °C, relative humidity 60 ± 10%, and CO2 concentration 1000 ± 200 ppm. Nutrient solution (HuaBao No. 1, Taiwan) was maintained at EC 1.2 mS cm⁻¹ and pH 6.5. Baseline lighting was provided by T8 LED tubes (Epistar, Taiwan) at a 16 h/8 h photoperiod with PPFD 130 μmol m⁻² s⁻¹.
LED type	Red (%)	Green (%)	Blue (%)	PPFD (μmol m⁻² s⁻¹)
LED 4000K	30	45	25	130

Table 6. Spectral composition of the primary LED lighting used for cultivation.

4.2. Experiment 1: Cultivar screening under UV-A and UV-B
4.2.1. Plant materials
Four red-leaf lettuce cultivars were used: 'Lollo Rosso' (LS-006, Known-You Seed, Taiwan), 'Purple Romaine' (LE4705, Jiase, Taiwan), 'Champagne Red Flame' (Dacheng, Taiwan) and 'Oakleaf Purple' (LE4806, Jiase, Taiwan).
4.2.2. Cultivation and UV treatments
UV treatments started at 32 days after sowing (DAS) and plants were harvested at 35 DAS (n = 5 per treatment for each cultivar). Treatment conditions are summarized in Table 7.
Code	Description	Lamp type	Peak wavelength (nm)	Irradiance (μW cm⁻²)	Hours per day (h)	Duration (days)
UV-A	UV-A irradiation	Fluorescent insect-attracting lamp	365	1104	4	3
UV-B	UV-B irradiation	Germicidal lamp	302	1510	4	3
CK	Control (no UV)	-	-	-	-	-

Table 7. Treatment conditions for the cultivar screening experiment (Experiment 1).

4.3. Experiment 2: UV-A recipe optimization for the selected cultivar
Based on Experiment 1, 'Lollo Rosso' was selected for optimization.
4.3.1. Cultivation procedure
Seeds were soaked for 5 h and sown in sponge cubes. Germination was performed under PPFD 100 μmol m⁻² s⁻¹. Seedlings were transplanted at 14 DAS and harvested at 35 DAS.
4.3.2. Daily UV-A energy calculation
To quantify UV-A dose across recipes, daily irradiation energy was calculated as:
E_daily = I_UVA × t_daily
where I_UVA = 11 W m⁻² (= 1100 μW cm⁻²) and t_daily is the daily irradiation time. Table 8a lists the calculated daily energies for different irradiation durations.
For ecological relevance, winter outdoor UV-A is on the order of 1500–3000 μW cm⁻² during effective sun hours [30]. Assuming a representative value of 2250 μW cm⁻² over 6 h (09:30–15:30), the estimated outdoor daily UV-A energy is:
E_outdoor = 2250 × 21600 / 10000 = 48.6 J cm⁻²
This estimate indicates that the 12 h day⁻¹ treatment (47.5 J cm⁻²) is comparable to winter sunlight, supporting the dose range used here [10].
Daily duration (h)	Seconds per day (s)	Daily energy (J cm⁻²)
3	10,800	11.9
6	21,600	23.8
9	32,400	35.6
12	43,200	47.5
15	54,000	59.4

Table 8a. Daily UV-A energy calculation for I_UVA = 11 W m⁻².

4.3.3. Training-set treatments (optimization recipes)
UV-A irradiance was fixed at 11 W m⁻² and different combinations of daily duration and total days were designed (Table 8b).
Code	Description	I_UVA (W m⁻²)	Hours per day (h)	Total days (d)	Start day
CK	Control (no UV-A)	0	0	0	-
UV-A-L6D6	Low daily dose	11	6	6	Day 29
UV-A-L6D6-N	Night irradiation	11	6 (night)	6	Day 29
UV-A-H12D3	High daily dose	11	12	3	Day 32
UV-A-VL3D12	Very low daily dose (long)	11	3	12	Day 23
UV-A-L6D12	Low daily dose (long)	11	6	12	Day 23

Table 8b. Training-set UV-A recipes for 'Lollo Rosso' (Experiment 2).

4.3.4. Validation-set treatments (3-day gradient)
To evaluate generalization, a 3-day gradient experiment was conducted with fixed duration (3 days) and varying daily irradiation time (0–15 h day⁻¹), as summarized in Table 8c.
Code	Description	I_UVA (W m⁻²)	Hours per day (h)	Total days (d)	Start day	Daily energy (J cm⁻²)
CK	Control	0	0	0	-	0
VL3D3	Very low daily dose	11	3	3	Day 32	11.9
L6D3	Low daily dose	11	6	3	Day 32	23.8
M9D3	Medium daily dose	11	9	3	Day 32	35.6
H12D3	High daily dose	11	12	3	Day 32	47.5
VH15D3	Very high daily dose	11	15	3	Day 32	59.4

Table 8c. Validation-set UV-A recipes (3-day gradient).

4.4. Measurements and data analysis
Fresh weight (FW) of shoots was measured at harvest. Leaf color parameters (L*, a*, b*) were measured using a colorimeter (NE-4000, Nippon Denshoku, Japan). Anthocyanin content was measured following a modified protocol of Hung et al. [31]: 0.5 g of leaf tissue was homogenized in potassium phosphate buffer, centrifuged, and the absorbance of the supernatant was measured at 600 nm.
4.5. Mechanistic model of UV-A effects
A six-state ODE model was developed based on the lettuce growth modeling framework of Sun et al. [19], and extended to represent the dual effects of UV-A: (i) UV-A-induced ROS accumulation leading to stress and growth inhibition, and (ii) UV-A-driven morphological changes that increase leaf area and thereby indirectly enhance PAR interception. A key innovation is the explicit representation of carbon competition between growth and antioxidant defense, implementing the Growth-Differentiation Balance Hypothesis [1].
The model comprises six state variables: structural dry weight (X_d), non-structural carbon buffer (C_buf), leaf area index (LAI), total antioxidants (AOX), a dimensionless stress state (Stress), and a dimensionless ROS variable (ROS). Anthocyanin concentration is derived as 18% of total AOX, reflecting the proportion of anthocyanins among total antioxidant compounds (flavonoids, ascorbate, etc.).
For reproducibility and to reduce main-text length, the full parameter list, initial conditions, and implementation details are provided in the Supplementary Materials.
4.5.1. State variables
The six state variables are: X_d (structural dry weight, kg m⁻²), C_buf (non-structural carbon buffer, kg m⁻²), LAI (m² m⁻²), AOX (total antioxidants, kg m⁻²), Stress (dimensionless), and ROS (dimensionless). Anthocyanin is calculated as Anth = AOX × 0.18.
4.5.2. Photosynthesis and energy conversion
The baseline growth model follows Sun et al. [19] and retains their energy balance and stomatal conductance formulation. Canopy photosynthesis is based on the Farquhar biochemical model [32] and integrates canopy layers using three-point Gaussian quadrature. Rubisco-limited and RuBP-regeneration-limited assimilation rates are:
A_c = Vc,max (Ci − Γ*) / (Ci + Kc (1 + O/Ko))  (1)
A_j = J (Ci − Γ*) / (4 Ci + 8 Γ*)  (2)
where J is the electron transport rate dependent on irradiance and temperature.
To harmonize units for artificial lighting, PPFD was converted to PAR power using a spectrum-weighted coefficient κ based on Table 6:
κ = 0.30×0.181 + 0.45×0.217 + 0.25×0.266 = 0.219 J µmol⁻¹
Thus, PPFD = 130 µmol m⁻² s⁻¹ corresponds to PAR ≈ 28.5 W m⁻². Because the Sun model assumes PAR = σ_PAR × I with σ_PAR = 0.5, we used an equivalent shortwave irradiance I_LED = 57 W m⁻² for the baseline light.
UV-A at 365 nm does not directly contribute to photosynthesis and is introduced as a driver of ROS production and as a trigger for morphological effects. The measured UV-A intensity was I_UVA = 11 W m⁻².
4.5.3. Core ODE system
(3) Carbon buffer dynamics with carbon competition:
dC_buf/dt = A_canopy − (r_g · dX_d/dt) − r_m · X_d − S_AOX · C_cost  (3)
where S_AOX is the AOX synthesis rate and C_cost = 1.0 kg C per kg AOX represents the carbon cost of antioxidant synthesis. This formulation explicitly implements the carbon competition between growth and defense [1,11].
(4) Structural growth dynamics with carbon competition penalty:
dX_d/dt = [C_buf / (C_buf + K_carbon)] · RGR_max · X_d · (1 − Stress_inhib) · (1 − CC_penalty)  (4)
The carbon competition penalty is calculated as:
CC_penalty = (AOX_carbon_effect × 0.30) + (0.225 × Stress / (21 + Stress))
where AOX_carbon_effect is a Hill-type saturation function of stress-induced AOX synthesis carbon demand:
AOX_carbon_effect = D_AOX / (K_cc + D_AOX)
D_AOX = S_AOX,stress × C_cost
Here, S_AOX,stress is the stress-induced component of AOX synthesis rate (kg m⁻² s⁻¹), C_cost = 1.0 kg C per kg AOX is the carbon cost of synthesis, and K_cc = 1×10⁻⁸ kg m⁻² s⁻¹ is the half-saturation constant. This Hill function ensures smooth transitions and bounded effects (0 ≤ AOX_carbon_effect ≤ 1). The first term (×0.30) represents the direct carbon diversion to AOX synthesis, while the second term represents the indirect effect of cumulative stress on resource allocation [11].
(5) LAI dynamics:
dLAI/dt = (dX_d/dt) · SLA_new − k_senescence · LAI  (5)
SLA_new = SLA_base × [1 + (k_LAI · I_UVA)/(K_LAI + I_UVA)]
(6) ROS dynamics:
dROS/dt = k_prod · I_UVA − k_clear · ROS  (6)
(7) Stress dynamics:
dStress/dt = (Damage + Circadian_Damage) − Decay  (7)
Damage = stress_coeff · ROS · (A_vuln · exp(−k_vuln · LAI) + 1) · Nonlinear_Factor · (1 − AOX_protection)
Nonlinear_Factor = 1 + F_max · exp(−exp(−k_steep (Hours − H_threshold)))
Circadian_Damage = k_circ · I_UVA · Hours_in_Dark^n_circ  (night-only)
(8) AOX dynamics:
dAOX/dt = Synthesis − Degradation − Consumption  (8)
AOX synthesis includes baseline, stress-induced and UV-induced components [33,34]. Anthocyanins constitute approximately 18% of total antioxidants, so Anth = AOX × 0.18. AOX consumption represents antioxidant usage for ROS scavenging under high oxidative pressure [35–37]. Additional smooth inhibition terms represent reduced synthetic efficiency under severe dehydration [21,38] and reduced carbon supply under low LAI [39,40].
4.5.4. Carbon competition mechanism
The carbon competition mechanism is a key innovation of this model, grounded in the Growth-Differentiation Balance Hypothesis [1]. Under UV-A stress, plants face a tradeoff: allocate carbon to growth (X_d) or to defense (AOX synthesis). This is implemented through two pathways:
1. Direct carbon consumption: AOX synthesis consumes carbon from C_buf at a rate proportional to synthesis rate × carbon cost (1.0 kg C per kg AOX). The carbon cost reflects the metabolic investment in phenylpropanoid biosynthesis, which can consume up to 20% of photosynthate [13].
2. Growth penalty: High stress and AOX synthesis rates reduce growth efficiency through the CC_penalty term, reflecting resource allocation to defense at the expense of growth [12].
The carbon consumption is capped at 10% of available C_buf per timestep for numerical stability:
AOX_carbon_consumption = min(S_AOX × C_cost, C_buf × 0.10)
4.6. Statistical analysis
Experimental data are presented as mean ± SD (n = 5). One-way ANOVA followed by Tukey's HSD test was performed in SPSS (p < 0.05).

5. Conclusions
We developed a two-stage framework integrating experiments and mechanistic modeling to enhance anthocyanin in red-leaf lettuce without yield loss. Experiments showed that UV-A is safer than UV-B and identified 'Lollo Rosso' as a highly responsive cultivar. Among tested recipes, L6D6 (6 h day⁻¹ for 6 days) maintained fresh weight while increasing anthocyanin concentration and total anthocyanin per plant.
The six-state mechanistic ODE model, featuring explicit carbon competition between growth and antioxidant defense based on the Growth-Differentiation Balance Hypothesis [1], reproduced training data with <5% error in 10 of 12 metrics and predicted an independent validation dataset with <10% error in 11 of 12 metrics. The model captures key nonlinear features including ontogeny-dependent vulnerability, nonlinear damage amplification, circadian disruption, carbon competition, and hormesis at extreme doses. Global optimization predicted that 9 h day⁻¹ for 4 days is the theoretical optimum, increasing total anthocyanin by ~38% (relative to control) with minimal yield reduction (−2.4%), and warrants future experimental verification.
Overall, the workflow of "experimental observation → mechanistic hypothesis → mathematical formulation → data-driven calibration → model-based optimization" offers a general paradigm for resolving yield–quality tradeoffs in controlled-environment agriculture. The carbon competition framework provides a mechanistic basis for understanding and optimizing the growth–defense balance in crops.

References
1. Herms, D.A.; Mattson, W.J. The dilemma of plants: to grow or defend. Q. Rev. Biol. 1992, 67, 283–335. https://doi.org/10.1086/285343
2. Lee, J.H.; Felipe, P.; Yang, Y.H.; Kim, M.Y.; Kwon, O.Y.; Sok, D.E.; Kim, H.C.; Kim, M.R. Effects of dietary supplementation with red-pigmented leafy lettuce on lipid profiles and antioxidant status in C57BL/6J mice fed a high-fat high-cholesterol diet. Br. J. Nutr. 2009, 101, 1246–1254. https://doi.org/10.1017/S0007114508073650
3. Paik, I.; Huq, E. Plant photoreceptors: multifunctional sensory proteins and their signaling networks. Semin. Cell Dev. Biol. 2019, 92, 114–121. https://doi.org/10.1016/j.semcdb.2019.03.007
4. Chaves, I.; Pokorny, R.; Byrdin, M.; Hoang, N.; Ritz, T.; Brettel, K.; Essen, L.O.; van der Horst, G.T.; Batschauer, A.; Ahmad, M. The cryptochromes: blue light photoreceptors in plants and animals. Annu. Rev. Plant Biol. 2011, 62, 335–364. https://doi.org/10.1146/annurev-arplant-042110-103759
5. Chen, X.; Ding, Y.; Yang, Y.; Song, C.; Wang, B.; Yang, S.; Guo, Y.; Gong, Z. HY5: A pivotal regulator of light-dependent development in higher plants. Front. Plant Sci. 2021, 12, 800989. https://doi.org/10.3389/fpls.2021.800989
6. Nguyen, N.H. HY5, an integrator of light and temperature signals in the regulation of anthocyanins biosynthesis in Arabidopsis. AIMS Mol. Sci. 2020, 7, 126–141. https://doi.org/10.3934/molsci.2020005
7. Stapleton, A.E. Ultraviolet radiation and plants: burning questions. Plant Cell 1992, 4, 1353–1358. https://doi.org/10.1105/tpc.4.11.1353
8. Tevini, M.; Teramura, A.H. UV-B effects on terrestrial plants. Photochem. Photobiol. 1989, 50, 479–487. https://doi.org/10.1111/j.1751-1097.1989.tb05552.x
9. Tsormpatsidis, E.; Henbest, R.G.C.; Davis, F.J.; Battey, N.H.; Hadley, P.; Wagstaffe, A. The influence of ultraviolet radiation on growth, photosynthesis and phenolic levels of green and red lettuce. Ann. Appl. Biol. 2010, 156, 357–366. https://doi.org/10.1111/j.1744-7348.2010.00393.x
10. Krizek, D.T.; Mirecki, R.M.; Britz, S.J. Inhibitory effects of ambient levels of solar UV-A and UV-B radiation on growth of cucumber. Physiol. Plant. 1997, 100, 886–893. https://doi.org/10.1111/j.1399-3054.1997.tb00014.x
11. Monson, R.K.; Trowbridge, A.M.; Lindroth, R.L.; Lerdau, M.T. Coordinated resource allocation to plant growth–defense tradeoffs. New Phytol. 2022, 233, 1051–1066. https://doi.org/10.1111/nph.17773
12. Gershenzon, J. Metabolic costs of terpenoid accumulation in higher plants. J. Chem. Ecol. 1994, 20, 1281–1328. https://doi.org/10.1007/BF02059810
13. Vogt, T. Phenylpropanoid biosynthesis. Mol. Plant 2010, 3, 2–20. https://doi.org/10.1093/mp/ssp106
14. Weiland, M.; Smink, H.S.; Kannangara, R.M.; Bak, S. A comparison of consistent UV treatment versus inconsistent UV treatment in horticultural production of lettuce. Photochem. Photobiol. Sci. 2023, 22, 1611–1624. https://doi.org/10.1007/s43630-023-00402-8
15. Li, Q.; Kubota, C. Effects of supplemental light quality on growth and phytochemicals of baby leaf lettuce. Environ. Exp. Bot. 2009, 67, 59–64. https://doi.org/10.1016/j.envexpbot.2009.06.011
16. Chen, Y.; Li, T.; Yang, Q.; Zhang, Y.; Zou, J.; Bian, Z.; Wen, X. UVA radiation is beneficial for yield and quality of indoor cultivated lettuce. Front. Plant Sci. 2019, 10, 1563. https://doi.org/10.3389/fpls.2019.01563
17. García-Macías, P.; Ordidge, M.; Vysini, E.; Warber, S.; Baker, A.; Hadley, P. Changes in the flavonoid and phenolic acid contents and antioxidant activity of red leaf lettuce (Lollo Rosso) due to cultivation under plastic films varying in ultraviolet transparency. J. Agric. Food Chem. 2007, 55, 10168–10172. https://doi.org/10.1021/jf071570m
18. van Henten, E.J. Validation of a dynamic lettuce growth model for greenhouse climate control. Agric. Syst. 1994, 45, 55–72. https://doi.org/10.1016/S0308-521X(94)90280-1
19. Sun, W.; Weiss, A.; Zhou, G. A lettuce growth model responding to a broad range of greenhouse climates. Biosyst. Eng. 2025, 250, 285–305. https://doi.org/10.1016/j.biosystemseng.2025.01.008
20. Hadacek, F.; Kreutzer, S.; Göhler, T. Plant secondary metabolites as defense systems against abiotic stresses: a comprehensive analysis of low-dose effects. Dose-Response 2010, 9, 60–75. https://doi.org/10.2203/dose-response.09-028.Hadacek
21. Ferreyra, M.L.F.; Muhlemann, J.K.; Grotewold, E. Ultraviolet-B radiation hormesis in plants: perspectives for the enhancement of secondary metabolite biosynthesis. Plant Physiol. Biochem. 2021, 165, 95–105. https://doi.org/10.1016/j.plaphy.2021.05.022
22. Lai, A.G.; Doherty, C.J.; Mueller-Roeber, B.; Kay, S.A.; Schippers, J.H.M.; Dijkwel, P.P. CIRCADIAN CLOCK-ASSOCIATED 1 regulates ROS homeostasis and oxidative stress responses. Proc. Natl. Acad. Sci. USA 2012, 109, 17129–17134. https://doi.org/10.1073/pnas.1209148109
23. Jiménez, A.; Sevilla, F.; Martí, M.C. The influence of circadian rhythm on the activity of oxidative stress enzymes. Int. J. Mol. Sci. 2022, 23, 14275. https://doi.org/10.3390/ijms232214275
24. Nitschke, S.; Cortleven, A.; Iven, T.; Feussner, I.; Havaux, M.; Riefler, M.; Schmülling, T. Circadian stress regimes affect the circadian clock and cause jasmonic acid-dependent cell death in cytokinin-deficient Arabidopsis plants. Plant Cell 2016, 28, 1616–1639. https://doi.org/10.1105/tpc.16.00016
25. Huot, B.; Yao, J.; Montgomery, B.L.; He, S.Y. Growth–defense tradeoffs in plants: a balancing act to optimize fitness. Mol. Plant 2014, 7, 1267–1287. https://doi.org/10.1093/mp/ssu049
26. Láposi, R.; Veres, S.; Lakatos, G.; Oláh, V.; Fieldsend, A.; Mészáros, I. Species-specific and leaf-age dependent effects of ultraviolet radiation on two Brassicaceae. Photochem. Photobiol. Sci. 2007, 6, 105–112. https://doi.org/10.1039/b612590f
27. Agati, G.; Stefano, G.; Biricolti, S.; Tattini, M. Mesophyll distribution of 'antioxidant' flavonoid glycosides in Ligustrum vulgare leaves under contrasting sunlight irradiance. Ann. Bot. 2011, 107, 383–390. https://doi.org/10.1093/aob/mcq237
28. Jenkins, G.I. Signal transduction in responses to UV-B radiation. Annu. Rev. Plant Biol. 2009, 60, 407–431. https://doi.org/10.1146/annurev.arplant.59.032607.092953
29. Czégény, G.; Mátai, A.; Hideg, É. UV-B effects on leaves—oxidative stress and acclimation in controlled environments. Plant Sci. 2016, 248, 57–63. https://doi.org/10.1016/j.plantsci.2016.04.013
30. Verdaguer, D.; Jansen, M.A.K.; Llorens, L.; Morales, L.O.; Neugart, S. UV-A radiation effects on higher plants: exploring the known unknown. Plant Sci. 2017, 255, 72–81. https://doi.org/10.1016/j.plantsci.2016.11.014
31. Hung, K.T.; Cheng, D.G.; Hsu, Y.T.; Kao, C.H. Abscisic acid-induced hydrogen peroxide is required for anthocyanin accumulation in leaves of rice seedlings. J. Plant Physiol. 2008, 165, 1280–1287. https://doi.org/10.1016/j.jplph.2007.10.008
32. Farquhar, G.D.; von Caemmerer, S.; Berry, J.A. A biochemical model of photosynthetic CO2 assimilation in leaves of C3 species. Planta 1980, 149, 78–90. https://doi.org/10.1007/BF00386231
33. Winkel-Shirley, B. Biosynthesis of flavonoids and effects of stress. Curr. Opin. Plant Biol. 2002, 5, 218–223. https://doi.org/10.1016/S1369-5266(02)00256-X
34. An, J.P.; Zhang, X.W.; Liu, Y.J.; Wang, X.F.; You, C.X.; Hao, Y.J. The photomorphogenic transcription factor PpHY5 regulates anthocyanin accumulation in response to UV-A and UV-B irradiation. Front. Plant Sci. 2020, 11, 603178. https://doi.org/10.3389/fpls.2020.603178
35. Gould, K.S. Nature's Swiss Army Knife: the diverse protective roles of anthocyanins in leaves. J. Biomed. Biotechnol. 2004, 2004, 314–320. https://doi.org/10.1155/S1110724304406147
36. Xu, Z.; Rothstein, S.J. ROS-induced anthocyanin production provides feedback protection by scavenging ROS and maintaining photosynthetic capacity in Arabidopsis. Plant Signal. Behav. 2018, 13, e1451708. https://doi.org/10.1080/15592324.2018.1451708
37. Cerqueira, J.V.A.; Silveira, J.A.G.; Carvalho, F.E.L.; Cunha, J.R.; Neto, M.C.L. Anthocyanins and reactive oxygen species: a team of rivals regulating plant development? Plant Mol. Biol. 2023, 112, 213–223. https://doi.org/10.1007/s11103-023-01362-4
38. Garnier, E.; Shipley, B.; Roumet, C.; Laurent, G. A standardized protocol for the determination of specific leaf area and leaf dry matter content. Funct. Ecol. 2001, 15, 688–695. https://doi.org/10.1046/j.0269-8463.2001.00563.x
39. Wang, L.; Wei, L.; Chen, W.; Zhang, T. Regulation of anthocyanin and sugar accumulation in grape berry through carbon limitation and exogenous ABA application. Food Res. Int. 2022, 160, 111478. https://doi.org/10.1016/j.foodres.2022.111478
40. Zhao, S.; Gao, H.; Luo, J.; Wang, H.; Dong, Q.; Wang, Y.; Yang, K.; Mao, K.; Ma, F. Anthocyanin accumulation provides protection against high light stress while reducing photosynthesis in apple leaves. Int. J. Mol. Sci. 2022, 23, 12616. https://doi.org/10.3390/ijms232012616
