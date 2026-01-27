"""
Lettuce Carbon Allocation Model (Base Growth Module)
=====================================================
Based on: Sun et al. (2025) "A mechanistic model for simulating lettuce
          growth and resource allocation"

This module implements the base lettuce growth model as a three-state ODE system:
- X_d: Structural dry weight [kg/m²]
- C_buf: Carbon buffer (non-structural carbohydrates) [kg/m²]
- LAI: Leaf area index [m²/m²]

Used as the foundation for the UVA effect model in:
Wei, C.H., Fang, W., & Huang, C.K. (2026). A Two-Stage Screening-to-Optimization
Approach with Mechanistic Model Analysis: Enhancing Anthocyanin in Lettuce
Without Yield Loss. Plants (under review).
"""
import numpy as np

class SunParams:
    def __init__(self):
        # Parameters same as previous version...
        self.M_CO2=44e-3; self.Rg=8.314; self.T0_K=273.15; self.c_alpha=0.68
        self.c_beta=0.8; self.sigma_buf=0.2; self.RGR_max_20=1.54e-6  # Corrected: Consistent with Sun original model
        self.Q10_gr=1.6; self.T_c_RGR=25.0; self.c_Rd_25_sh=3.47e-7
        self.c_Rd_25_r=1.16e-7; self.Q10_Rd=2.0; self.cr_I=0.22; self.cr_PAR=0.07
        self.kI=0.48; self.kPAR=0.9; self.sigma_PAR=0.5; self.SLA_ref=47.93
        self.Ia_L_ref=50.3; self.Xh_ref=0.75; self.beta_I=-4.74e-3
        self.beta_Xh=0.912; self.c_sigma_r_1=-0.026; self.c_sigma_r_2=-0.076
        self.epsilon_0=17e-9; self.Gamma_T20=40.0; self.Q10_Gamma=2.0
        self.Jmax_25=210.15; self.EJ=3.7e4; self.cH=2.2e5; self.cS=710.0
        self.T25_K=298.15; self.c_zeta=1.6; self.r_H2O_min=82.0; self.Le=1.47
        self.lf=0.1; self.va=0.09; self.rt=50.0; self.c_rc_1=0.315
        self.c_rc_2=-27.35; self.c_rc_3=790.7; self.rho_CO2_T0=1.98

def sun_derivatives_final(t, state, p, env):
    """
    Sun Model Differential Equations

    Parameters:
    -----
    t : float - Time [seconds]
    state : array - [X_d, C_buf, LAI]
    p : SunParams - Parameter object
    env : dict - Environment settings, may include:
        - I_override: If provided, use this irradiance directly (overrides day/night logic)
        - T_override: If provided, use this temperature directly
        - is_day_override: If provided, force this day/night state (for temperature, CO2, RH)
    """
    # 1. Unpack three state variables
    X_d, C_buf, LAI = state
    X_d, C_buf, LAI = max(X_d, 1e-9), max(C_buf, 0), max(LAI, 1e-9)

    # 2. Get environmental conditions
    hour = (t / 3600) % 24

    # Support day/night determination across midnight
    light_on = env['light_on_hour']
    light_off = env['light_off_hour']
    if light_on <= light_off:
        is_day_internal = light_on <= hour < light_off
    else:
        is_day_internal = hour >= light_on or hour < light_off

    # Allow external override of day/night state (for temperature, CO2, RH)
    is_day = env.get('is_day_override', is_day_internal)

    # Irradiance: prioritize I_override (supports nighttime UVA-PAR photosynthetic contribution)
    if 'I_override' in env:
        I = env['I_override']
    else:
        I = env['I_day'] if is_day else 0.0

    # Temperature: prioritize T_override
    if 'T_override' in env:
        Tc = env['T_override']
    else:
        Tc = env['T_day'] if is_day else env['T_night']

    # CO2 and RH
    Xc_ppm = env['CO2_day'] if is_day else env['CO2_night']
    Xh = env['RH_day'] if is_day else env['RH_night']

    Tc_K = Tc + p.T0_K

    # 3. Calculate auxiliary variables
    plant_dw = X_d / env['plant_density']; sr_val = np.clip(p.c_sigma_r_1*np.log(plant_dw+1e-9)+p.c_sigma_r_2,0.05,0.35)
    I_a = (1 - p.cr_I) * I * (1 - np.exp(-p.kI * LAI))
    Ia_pl = I_a / (LAI + 1e-9)
    f_I_SLA = 1 / (1 + p.beta_I * (p.Ia_L_ref - Ia_pl)); f_Xh_SLA = 1 / (1 + p.beta_Xh * (p.Xh_ref - Xh)); SLA = p.SLA_ref * f_I_SLA * f_Xh_SLA

    # 4. Total photosynthesis rate A_C
    Gamma = p.Gamma_T20 * (p.Q10_Gamma**((Tc - 20) / 10)); eps = p.epsilon_0 * (Xc_ppm - Gamma) / (Xc_ppm + 2 * Gamma + 1e-9)
    Jmax = p.Jmax_25 * np.exp(p.EJ * (Tc_K - p.T25_K) / (Tc_K * p.Rg * p.T25_K)) * (1 + np.exp((p.cS * p.T25_K - p.cH) / (p.Rg * p.T25_K))) / (1 + np.exp((p.cS * Tc_K - p.cH) / (p.Rg * Tc_K)) + 1e-9)
    AL_mm = p.M_CO2 * Jmax / 4.0 * 1e-6; rc = max((p.c_rc_1 * Tc**2 + p.c_rc_2 * Tc + p.c_rc_3), 10.0); rb = p.Le**0.67 * 1174 * p.lf**0.5 / ((p.lf * abs(Tc - Tc) + 207 * p.va**2)**0.25 + 1e-9)
    es = 10**(2.7857 + 7.5 * Tc / (237.3 + Tc)); ec_a = es * (1 - Xh); fXh_s = 4.0 / ((1 + 255 * np.exp(-0.54e-2 * ec_a))**0.25 + 1e-9)
    fXc_s = 1 + 6.1e-7 * (Xc_ppm - 200)**2 if I > 3 and Xc_ppm < 1100 else (1.5 if I > 3 else 1.0); fTc_s = 1 + 0.5e-2 * (Tc - 33.6)**2 if I <= 3 else 1 + 2.3e-2 * (Tc - 24.5)**2
    fI_s = (I_a / (2 * LAI + 1e-9) + 4.3) / (I_a / (2 * LAI + 1e-9) + 0.54); rs = p.c_zeta * p.r_H2O_min * fI_s * fTc_s * fXc_s * fXh_s; rCO2 = rs + rb + rc + p.rt
    rho = p.rho_CO2_T0 * p.T0_K / (Tc_K + 1e-9); AL_cn = max(rho * (Xc_ppm - Gamma) / (rCO2 + 1e-9) * 1e-6, 0.0); AL_sat_n = min(AL_cn, AL_mm)
    R_d_for_A_L_sat = (p.c_Rd_25_sh * (1 - sr_val) + p.c_Rd_25_r * sr_val) * X_d * (p.Q10_Rd**((Tc - 25) / 10)); AL_sat = max(AL_sat_n + (R_d_for_A_L_sat / (LAI + 1e-9)) / p.c_alpha, 0.0)
    # Correction: Use 3-point Gaussian integration for canopy photosynthesis (Eq. 7-8)
    l_1 = (0.5 - np.sqrt(0.15)) * LAI; l_2 = 0.5 * LAI; l_3 = (0.5 + np.sqrt(0.15)) * LAI
    PARa_1 = p.kPAR * (1 - p.cr_PAR) * I * p.sigma_PAR * np.exp(-p.kPAR * l_1)
    PARa_2 = p.kPAR * (1 - p.cr_PAR) * I * p.sigma_PAR * np.exp(-p.kPAR * l_2)
    PARa_3 = p.kPAR * (1 - p.cr_PAR) * I * p.sigma_PAR * np.exp(-p.kPAR * l_3)
    A_L_1 = max(AL_sat * (1 - np.exp(-eps * PARa_1 / (AL_sat + 1e-9))), 0.0)
    A_L_2 = max(AL_sat * (1 - np.exp(-eps * PARa_2 / (AL_sat + 1e-9))), 0.0)
    A_L_3 = max(AL_sat * (1 - np.exp(-eps * PARa_3 / (AL_sat + 1e-9))), 0.0)
    A_L_C = (A_L_1 + 1.6 * A_L_2 + A_L_3) / 3.6; A_C = A_L_C * LAI

    # 5. Calculate carbon fluxes
    R_d = (p.c_Rd_25_sh * (1 - sr_val) + p.c_Rd_25_r * sr_val) * X_d * (p.Q10_Rd**((Tc - 25) / 10))
    C_buf_max = p.sigma_buf * X_d
    RGR_max = p.RGR_max_20 * (p.Q10_gr**(((Tc - 20) if Tc <= p.T_c_RGR else -(Tc - 20)) / 10))
    h_buf = 1.0
    if C_buf >= C_buf_max: h_buf = min((R_d + (RGR_max * X_d / p.c_beta)) / (p.c_alpha * A_C + 1e-9), 1.0)

    # 6. Calculate differential equations
    # Original Sun et al. model equations (restored)
    net_assimilation_term = p.c_alpha * A_C * h_buf - R_d

    dXd_dt = p.c_beta * net_assimilation_term # Eq. 1

    growth_consumption = RGR_max * X_d
    dCbuf_dt = net_assimilation_term - (growth_consumption / p.c_beta) # Eq. 5

    dLAI_dt = dXd_dt * (1 - sr_val) * SLA # Eq. 9

    # Boundary conditions for C_buf
    # v11.0 fix: Enforce non-negativity of C_buf
    # Problem: ODE solvers can overshoot even with derivative protection
    # Solution: When C_buf is very small or negative, prevent further decrease
    C_buf_min_threshold = 1e-5  # Small positive threshold
    if C_buf <= C_buf_min_threshold:
        if dCbuf_dt < 0:
            # Strongly suppress negative derivative when C_buf is low
            # Use exponential damping to smoothly approach zero
            dCbuf_dt = dCbuf_dt * max(0.0, C_buf / C_buf_min_threshold) ** 2
    if C_buf >= C_buf_max and dCbuf_dt > 0:
        dCbuf_dt = 0
    if X_d < (0.03 / 1000 * env['plant_density']) and dXd_dt < 0: dXd_dt = 0
    if LAI < 0.01 and dLAI_dt < 0: dLAI_dt = 0 # Added protection for LAI

    return np.array([dXd_dt, dCbuf_dt, dLAI_dt])
