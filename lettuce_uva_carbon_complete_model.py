"""
莴苣UVA-碳分配完整动态模型
================================
版本: v7.1 - 支援外部光照覆蓋
作者: Roo & 研究团队
日期: 2025-12-11

核心进展:
1. 严格遵循论文，将模型重构为 [X_d, C_buf, LAI] 三状态变量系统。
2. dLAI/dt 作为独立的微分方程实现，移除所有近似计算。
3. 此文件现在是一个纯粹的、与论文100%对应的可导入模型。

修改紀錄:
- v7.1 (2025-12-11): 新增 I_override 參數，允許外部覆蓋光照強度（支援夜間 UVA-PAR 光合貢獻）
"""
import numpy as np

class SunParams:
    def __init__(self):
        # 参数与之前版本相同...
        self.M_CO2=44e-3; self.Rg=8.314; self.T0_K=273.15; self.c_alpha=0.68
        self.c_beta=0.8; self.sigma_buf=0.2; self.RGR_max_20=1.54e-6  # 修正：與 Sun 原始模型一致
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
    Sun 模型微分方程

    參數:
    -----
    t : float - 時間 [秒]
    state : array - [X_d, C_buf, LAI]
    p : SunParams - 參數物件
    env : dict - 環境設定，可包含:
        - I_override: 若提供，直接使用此光照強度（覆蓋日夜判斷）
        - T_override: 若提供，直接使用此溫度
        - is_day_override: 若提供，強制使用此日夜狀態（用於溫度、CO2、RH）
    """
    # 1. 解包三状态变量
    X_d, C_buf, LAI = state
    X_d, C_buf, LAI = max(X_d, 1e-9), max(C_buf, 0), max(LAI, 1e-9)

    # 2. 获取环境条件
    hour = (t / 3600) % 24

    # 支援跨夜的日夜判斷
    light_on = env['light_on_hour']
    light_off = env['light_off_hour']
    if light_on <= light_off:
        is_day_internal = light_on <= hour < light_off
    else:
        is_day_internal = hour >= light_on or hour < light_off

    # 允許外部覆蓋日夜狀態（用於決定溫度、CO2、RH）
    is_day = env.get('is_day_override', is_day_internal)

    # 光照強度: 優先使用 I_override（支援夜間 UVA-PAR 光合貢獻）
    if 'I_override' in env:
        I = env['I_override']
    else:
        I = env['I_day'] if is_day else 0.0

    # 溫度: 優先使用 T_override
    if 'T_override' in env:
        Tc = env['T_override']
    else:
        Tc = env['T_day'] if is_day else env['T_night']

    # CO2 和 RH
    Xc_ppm = env['CO2_day'] if is_day else env['CO2_night']
    Xh = env['RH_day'] if is_day else env['RH_night']

    Tc_K = Tc + p.T0_K
    
    # 3. 计算辅助变量
    plant_dw = X_d / env['plant_density']; sr_val = np.clip(p.c_sigma_r_1*np.log(plant_dw+1e-9)+p.c_sigma_r_2,0.05,0.35)
    I_a = (1 - p.cr_I) * I * (1 - np.exp(-p.kI * LAI))
    Ia_pl = I_a / (LAI + 1e-9)
    f_I_SLA = 1 / (1 + p.beta_I * (p.Ia_L_ref - Ia_pl)); f_Xh_SLA = 1 / (1 + p.beta_Xh * (p.Xh_ref - Xh)); SLA = p.SLA_ref * f_I_SLA * f_Xh_SLA
    
    # 4. 总光合速率 A_C
    Gamma = p.Gamma_T20 * (p.Q10_Gamma**((Tc - 20) / 10)); eps = p.epsilon_0 * (Xc_ppm - Gamma) / (Xc_ppm + 2 * Gamma + 1e-9)
    Jmax = p.Jmax_25 * np.exp(p.EJ * (Tc_K - p.T25_K) / (Tc_K * p.Rg * p.T25_K)) * (1 + np.exp((p.cS * p.T25_K - p.cH) / (p.Rg * p.T25_K))) / (1 + np.exp((p.cS * Tc_K - p.cH) / (p.Rg * Tc_K)) + 1e-9)
    AL_mm = p.M_CO2 * Jmax / 4.0 * 1e-6; rc = max((p.c_rc_1 * Tc**2 + p.c_rc_2 * Tc + p.c_rc_3), 10.0); rb = p.Le**0.67 * 1174 * p.lf**0.5 / ((p.lf * abs(Tc - Tc) + 207 * p.va**2)**0.25 + 1e-9)
    es = 10**(2.7857 + 7.5 * Tc / (237.3 + Tc)); ec_a = es * (1 - Xh); fXh_s = 4.0 / ((1 + 255 * np.exp(-0.54e-2 * ec_a))**0.25 + 1e-9)
    fXc_s = 1 + 6.1e-7 * (Xc_ppm - 200)**2 if I > 3 and Xc_ppm < 1100 else (1.5 if I > 3 else 1.0); fTc_s = 1 + 0.5e-2 * (Tc - 33.6)**2 if I <= 3 else 1 + 2.3e-2 * (Tc - 24.5)**2
    fI_s = (I_a / (2 * LAI + 1e-9) + 4.3) / (I_a / (2 * LAI + 1e-9) + 0.54); rs = p.c_zeta * p.r_H2O_min * fI_s * fTc_s * fXc_s * fXh_s; rCO2 = rs + rb + rc + p.rt
    rho = p.rho_CO2_T0 * p.T0_K / (Tc_K + 1e-9); AL_cn = max(rho * (Xc_ppm - Gamma) / (rCO2 + 1e-9) * 1e-6, 0.0); AL_sat_n = min(AL_cn, AL_mm)
    R_d_for_A_L_sat = (p.c_Rd_25_sh * (1 - sr_val) + p.c_Rd_25_r * sr_val) * X_d * (p.Q10_Rd**((Tc - 25) / 10)); AL_sat = max(AL_sat_n + (R_d_for_A_L_sat / (LAI + 1e-9)) / p.c_alpha, 0.0)
    # 修正：使用 3 點高斯積分計算冠層光合 (Eq. 7-8)
    l_1 = (0.5 - np.sqrt(0.15)) * LAI; l_2 = 0.5 * LAI; l_3 = (0.5 + np.sqrt(0.15)) * LAI
    PARa_1 = p.kPAR * (1 - p.cr_PAR) * I * p.sigma_PAR * np.exp(-p.kPAR * l_1)
    PARa_2 = p.kPAR * (1 - p.cr_PAR) * I * p.sigma_PAR * np.exp(-p.kPAR * l_2)
    PARa_3 = p.kPAR * (1 - p.cr_PAR) * I * p.sigma_PAR * np.exp(-p.kPAR * l_3)
    A_L_1 = max(AL_sat * (1 - np.exp(-eps * PARa_1 / (AL_sat + 1e-9))), 0.0)
    A_L_2 = max(AL_sat * (1 - np.exp(-eps * PARa_2 / (AL_sat + 1e-9))), 0.0)
    A_L_3 = max(AL_sat * (1 - np.exp(-eps * PARa_3 / (AL_sat + 1e-9))), 0.0)
    A_L_C = (A_L_1 + 1.6 * A_L_2 + A_L_3) / 3.6; A_C = A_L_C * LAI

    # 5. 计算碳流
    R_d = (p.c_Rd_25_sh * (1 - sr_val) + p.c_Rd_25_r * sr_val) * X_d * (p.Q10_Rd**((Tc - 25) / 10))
    C_buf_max = p.sigma_buf * X_d
    RGR_max = p.RGR_max_20 * (p.Q10_gr**(((Tc - 20) if Tc <= p.T_c_RGR else -(Tc - 20)) / 10))
    h_buf = 1.0
    if C_buf >= C_buf_max: h_buf = min((R_d + (RGR_max * X_d / p.c_beta)) / (p.c_alpha * A_C + 1e-9), 1.0)

    # 6. 计算微分方程
    net_assimilation_term = p.c_alpha * A_C * h_buf - R_d
    
    dXd_dt = p.c_beta * net_assimilation_term # Eq. 1
    
    growth_consumption = RGR_max * X_d
    dCbuf_dt = net_assimilation_term - (growth_consumption / p.c_beta) # Eq. 5
    
    dLAI_dt = dXd_dt * (1 - sr_val) * SLA # Eq. 9

    # 边界条件
    if (C_buf <= 0 and dCbuf_dt < 0) or (C_buf >= C_buf_max and dCbuf_dt > 0): dCbuf_dt = 0
    if X_d < (0.03 / 1000 * env['plant_density']) and dXd_dt < 0: dXd_dt = 0
    if LAI < 0.01 and dLAI_dt < 0: dLAI_dt = 0 # 增加对LAI的保护
        
    return np.array([dXd_dt, dCbuf_dt, dLAI_dt])
