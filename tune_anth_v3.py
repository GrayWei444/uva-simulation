"""系統性花青素調整 v3"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

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

    day_factor = 1.0 + p.k_day * (hours_today ** p.n_day)

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

    vuln = p.cap_vuln * (p.LAI_ref_vuln / LAI) ** p.n_vuln / (p.cap_vuln + (p.LAI_ref_vuln / LAI) ** p.n_vuln)
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
    'k_day': 1e-5, 'n_day': 7.0,
    'k_night': 0.15, 'n_night': 2.0,
}

ANTH_TARGETS = {'CK': 43.3, 'L6D6': 49.4, 'L6D6-N': 49.3, 'VL3D12': 48.2, 'L6D12': 51.8, 'H12D3': 65.1}


def test_one(treat, params):
    p = P(params)
    env = get_env_for_treatment(treat)
    fw_init = SIMULATION['initial_fw_g']
    dw_init = fw_init * p.dw_fw_ratio_base
    Xd_init = dw_init / 1000 * ENV_BASE['plant_density']
    y0 = [Xd_init, Xd_init * 0.1, (dw_init / 0.01) * 0.04,
          5.0 * fw_init * ENV_BASE['plant_density'] / 1000 / 1e6, 0.0]
    t0 = SIMULATION['transplant_offset'] * 86400
    t1 = (SIMULATION['transplant_offset'] + SIMULATION['days']) * 86400
    sol = solve_ivp(deriv, (t0, t1), y0, args=(p, env), method='RK45', max_step=1800, t_eval=[t1])
    if sol.success:
        Xd, Cbuf, LAI, Anth, Stress = sol.y[:, -1]
        dw_fw = calc_dw_fw(Stress, p)
        FW = Xd / ENV_BASE['plant_density'] / dw_fw * 1000
        FW_kg = FW / 1000 * ENV_BASE['plant_density']
        Anth_ppm = Anth / FW_kg * 1e6
        return FW, Anth_ppm, Stress
    return None


print("=" * 70)
print("系統性花青素調整 v3")
print("=" * 70)

# 策略：
# 1. 用較低的 base_rate (1.84e-10) 讓 CK 接近目標
# 2. 調整 V_max 和 K_stress_anth 讓其他組達標

# 問題分析:
# - CK (S=0): 只有基礎合成 → 需要降低 base
# - L6D6 (S=0.29): 幾乎達標
# - L6D6-N (S=1.85): +6.7% → 需要降低中等 Stress 的誘導
# - VL3D12 (S=1.69): +9.3% → 同上
# - L6D12 (S=6.73): +2.3% → 達標
# - H12D3 (S=36.6): -2.0% → 需要提高高 Stress 的誘導

# 解決方案: 提高 K_stress_anth
# - K 越大，需要更高的 Stress 才能達到飽和
# - 這會降低中等 Stress 組的誘導，提高相對於高 Stress 組的差異

print("\n策略: 提高 K_stress_anth 來拉開 Stress 誘導的差異")
print()

# 細緻搜索
best_max_err = 100
best_params = None
best_results = None

for base in np.linspace(1.7e-10, 1.9e-10, 5):
    for v_max in np.linspace(2.0e-11, 3.5e-11, 8):
        for k_anth in np.linspace(0.3, 1.5, 10):
            params = BASE.copy()
            params['base_anth_rate_light'] = base
            params['base_anth_rate_dark'] = base / 2
            params['V_max_anth'] = v_max
            params['K_stress_anth'] = k_anth

            errs = []
            results = {}
            for t in ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']:
                r = test_one(t, params)
                if r:
                    err = (r[1] - ANTH_TARGETS[t]) / ANTH_TARGETS[t] * 100
                    errs.append(abs(err))
                    results[t] = r

            if errs:
                max_err = max(errs)
                if max_err < best_max_err:
                    best_max_err = max_err
                    best_params = params.copy()
                    best_results = results.copy()
                    if max_err < 5:
                        print(f"找到! base={base:.2e}, V_max={v_max:.2e}, K={k_anth:.2f}: max_err={max_err:.1f}%")

print(f"\n最佳 max_err: {best_max_err:.1f}%")
if best_params:
    print(f"參數:")
    print(f"  base_anth_rate_light = {best_params['base_anth_rate_light']:.2e}")
    print(f"  V_max_anth = {best_params['V_max_anth']:.2e}")
    print(f"  K_stress_anth = {best_params['K_stress_anth']:.2f}")

    print("\n詳細結果:")
    print("-" * 70)
    for t in ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']:
        if t in best_results:
            FW, Anth, Stress = best_results[t]
            fw_target = TARGETS[t]['FW']
            anth_target = ANTH_TARGETS[t]
            fw_err = (FW - fw_target) / fw_target * 100
            anth_err = (Anth - anth_target) / anth_target * 100
            s1 = "✓" if abs(fw_err) < 5 else "✗"
            s2 = "✓" if abs(anth_err) < 5 else "✗"
            print(f"{t:<8} FW:{FW:>5.1f}({fw_err:>+5.1f}%{s1}) Anth:{Anth:>5.1f}({anth_err:>+5.1f}%{s2}) S:{Stress:>5.1f}")
