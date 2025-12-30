"""純連續版本 v6 - 測試 n=7, k=1e-5 的完整結果"""
import numpy as np
from scipy.integrate import solve_ivp
import sys
sys.path.insert(0, '/home/kasm-user/projects/uva-simulation')

from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment
from lettuce_uva_carbon_complete_model import SunParams as BaseSunParams
from lettuce_uva_carbon_complete_model import sun_derivatives_final


class P(BaseSunParams):
    def __init__(self, params):
        super().__init__()
        for k, v in params.items():
            setattr(self, k, v)


def calc_dw_fw(S, p):
    e = p.ldmc_stress_sensitivity * S / (p.K_ldmc + S + 1e-9)
    return min(p.dw_fw_ratio_base * (1 + e), p.dw_fw_ratio_max)


def deriv(t, state, p, env):
    X_d, C_buf, LAI, Anth, Stress = state
    X_d = max(X_d, 1e-9)
    C_buf = max(C_buf, 0)
    LAI = max(LAI, 0.1)
    Anth = max(Anth, 0)
    Stress = max(Stress, 0)

    hour = (t / 3600) % 24
    day = t / 86400
    light_on = env['light_on_hour']
    light_off = env['light_off_hour']
    is_day = light_on <= hour < light_off
    I_base = env['I_day'] if is_day else 0
    Tc = env['T_day'] if is_day else env['T_night']

    I_UVA = 0
    hours_today = 0
    uva_hour_on = env.get('uva_hour_on', 10)
    uva_hour_off = env.get('uva_hour_off', 16)

    if env.get('uva_on'):
        uva_start = env['uva_start_day']
        uva_end = env['uva_end_day']
        if uva_hour_on <= uva_hour_off:
            if uva_start <= day <= uva_end and uva_hour_on <= hour < uva_hour_off:
                I_UVA = env['uva_intensity']
                hours_today = hour - uva_hour_on
        else:
            if hour >= uva_hour_on:
                if uva_start <= day <= uva_end:
                    I_UVA = env['uva_intensity']
                    hours_today = hour - uva_hour_on
            elif hour < uva_hour_off:
                if uva_start <= day - 1 <= uva_end:
                    I_UVA = env['uva_intensity']
                    hours_today = (24 - uva_hour_on) + hour

    # ===== 連續機制 1: 日內逆境 (冪函數) =====
    # factor = 1 + k × hours^n
    # n=7 讓長時間照射 (12h) 產生遠大於短時間 (6h) 的效果
    # 12^7 / 6^7 = 128 倍
    day_factor = 1.0 + p.k_day * (hours_today ** p.n_day)

    # ===== 連續機制 2: 夜間節律 =====
    if light_on <= hour < light_off:
        hours_in_dark = 0
    else:
        hours_in_dark = hour - light_off if hour >= light_off else hour + (24 - light_off)
    circ = 1.0 + p.k_night * (hours_in_dark ** p.n_night) if (I_UVA > 0 and hours_in_dark > 0) else 1.0

    env_mod = env.copy()
    env_mod['I_override'] = I_base + I_UVA * p.par_conversion_factor
    env_mod['T_override'] = Tc
    env_mod['is_day_override'] = is_day

    dXd_base, dCbuf, dLAI_base = sun_derivatives_final(t, [X_d, C_buf, LAI], p, env_mod)

    # ===== 連續機制 3: LAI 脆弱性 =====
    vuln = p.cap_vuln * (p.LAI_ref_vuln / LAI) ** p.n_vuln / (p.cap_vuln + (p.LAI_ref_vuln / LAI) ** p.n_vuln)

    # ===== 連續機制 4: Stress 非線性 =====
    nl = 1 + p.stress_nonlinear_coeff * Stress / (p.K_nonlinear + Stress + 1e-9)

    damage = p.stress_damage_coeff * I_UVA * vuln * nl * circ * day_factor
    repair_cap = p.base_repair_capacity + p.carbon_repair_bonus * C_buf / (p.K_carbon + C_buf + 1e-9)
    repair = p.stress_repair_coeff * Stress * repair_cap
    dStress = damage - repair

    inh = Stress / (p.K_stress + Stress + 1e-9)
    dXd = dXd_base * (1 - p.stress_photosynthesis_inhibition * inh) if dXd_base > 0 else dXd_base
    dLAI = dLAI_base * (1 - p.stress_lai_inhibition * inh) if dLAI_base > 0 else dLAI_base

    dw_fw = calc_dw_fw(Stress, p)
    FW = X_d / dw_fw
    synth = FW * ((p.base_anth_rate_light if is_day else p.base_anth_rate_dark) +
                  p.V_max_anth * Stress / (p.K_stress_anth + Stress + 1e-12))
    dAnth = synth - p.k_deg * Anth
    dCbuf = dCbuf - repair * p.repair_carbon_cost

    return np.array([dXd, dCbuf, dLAI, dAnth, dStress])


def test_all(params, label=""):
    print(f"\n{label}")
    print("-" * 80)

    results = {}
    for treat in ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']:
        p = P(params)
        env = get_env_for_treatment(treat)
        fw_init = SIMULATION['initial_fw_g']
        dw_init = fw_init * p.dw_fw_ratio_base
        Xd_init = dw_init / 1000 * ENV_BASE['plant_density']
        y0 = [Xd_init, Xd_init * 0.1, (dw_init / 0.01) * 0.04,
              5.0 * fw_init * ENV_BASE['plant_density'] / 1000 / 1e6, 0.0]

        t0 = SIMULATION['transplant_offset'] * 86400
        t1 = (SIMULATION['transplant_offset'] + SIMULATION['days']) * 86400

        sol = solve_ivp(deriv, (t0, t1), y0, args=(p, env), method='RK45', max_step=300, t_eval=[t1])

        if sol.success:
            Xd, Cbuf, LAI, Anth, Stress = sol.y[:, -1]
            dw_fw = calc_dw_fw(Stress, p)
            FW = Xd / ENV_BASE['plant_density'] / dw_fw * 1000
            FW_kg = FW / 1000 * ENV_BASE['plant_density']
            Anth_ppm = Anth / FW_kg * 1e6
            results[treat] = (FW, Anth_ppm, LAI, Stress)

    fw_errs = []
    anth_errs = []
    for t in ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']:
        if t in results:
            FW, Anth, LAI, Stress = results[t]
            target = TARGETS[t]
            fw_err = (FW - target['FW']) / target['FW'] * 100
            anth_err = (Anth - target['Anth']) / target['Anth'] * 100
            fw_errs.append(abs(fw_err))
            anth_errs.append(abs(anth_err))
            s1 = 'ok' if abs(fw_err) < 5 else 'NG'
            s2 = 'ok' if abs(anth_err) < 10 else 'NG'
            print(f'{t:<8} FW:{FW:>5.1f}g({fw_err:>+5.1f}%{s1}) Anth:{Anth:>5.1f}({anth_err:>+5.1f}%{s2}) LAI:{LAI:>4.1f} S:{Stress:>5.1f}')

    fw_ok = sum(1 for e in fw_errs if e < 5)
    anth_ok = sum(1 for e in anth_errs if e < 10)
    print("-" * 80)
    print(f'達標: FW {fw_ok}/6, Anth {anth_ok}/6, Total {fw_ok + anth_ok}/12')
    return fw_ok + anth_ok


BASE = {
    'c_alpha': 0.555, 'par_conversion_factor': 1.0,
    'stress_damage_coeff': 0.66e-6, 'stress_repair_coeff': 1e-5,
    'stress_nonlinear_coeff': 8.0, 'K_nonlinear': 0.8,
    'LAI_ref_vuln': 6.5, 'n_vuln': 8, 'cap_vuln': 100.0,
    'stress_photosynthesis_inhibition': 0.66, 'stress_lai_inhibition': 0.66,
    'K_stress': 1.9, 'base_repair_capacity': 0.5, 'carbon_repair_bonus': 0.5,
    'K_carbon': 0.001, 'repair_carbon_cost': 1e-6,
    'base_anth_rate_light': 2e-10, 'base_anth_rate_dark': 1e-10,
    'V_max_anth': 2.35e-11, 'K_stress_anth': 0.30, 'k_deg': 3.02e-6,
    'dw_fw_ratio_base': 0.05, 'ldmc_stress_sensitivity': 1.0,
    'K_ldmc': 50.0, 'dw_fw_ratio_max': 0.12,
    # 連續因子參數
    'k_day': 1e-5,     # 日內逆境係數
    'n_day': 7.0,      # 日內逆境冪次
    'k_night': 0.1,    # 夜間節律係數
    'n_night': 2.0,    # 夜間節律冪次
}

print('=' * 80)
print('純連續版本 v6 - 基於 n=7, k=1e-5')
print('=' * 80)
print()
print('連續機制:')
print('  1. 日內逆境: factor = 1 + 1e-5 × hours^7')
print('     - 6h: 1 + 1e-5 × 279936 = 3.8')
print('     - 12h: 1 + 1e-5 × 35831808 = 359.3')
print('  2. 夜間節律: circ = 1 + 0.1 × hours_in_dark^2')
print('  3. LAI 脆弱性: sigmoid 型')
print('  4. Stress 非線性: Michaelis-Menten 型')

test_all(BASE, "基礎設定 (n=7, k=1e-5)")

# 調整夜間節律參數
print('\n調整夜間節律參數以改善 L6D6-N:')
for k_night in [0.05, 0.1, 0.15, 0.2, 0.3]:
    params = BASE.copy()
    params['k_night'] = k_night
    test_all(params, f"k_night={k_night}")

# 微調 k_day
print('\n微調 k_day:')
for k_day in [0.5e-5, 0.8e-5, 1e-5, 1.2e-5, 1.5e-5]:
    params = BASE.copy()
    params['k_day'] = k_day
    test_all(params, f"k_day={k_day:.1e}")
