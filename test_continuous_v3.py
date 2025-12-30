"""純連續版本 v3 - 系統性參數搜索

目標：找到連續函數參數使所有處理達標
機制：
1. 日內逆境：factor = 1 + k_day × (hours)^n_day
2. 夜間節律：circ = 1 + k_night × (hours_in_dark)^n_night
3. LAI 脆弱性：vuln = cap × (LAI_ref/LAI)^n / (cap + ...)
4. Stress 非線性：nl = 1 + k × Stress / (K + Stress)
"""
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

    # UVA 排程
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
            # 跨夜
            if hour >= uva_hour_on:
                if uva_start <= day <= uva_end:
                    I_UVA = env['uva_intensity']
                    hours_today = hour - uva_hour_on
            elif hour < uva_hour_off:
                if uva_start <= day - 1 <= uva_end:
                    I_UVA = env['uva_intensity']
                    hours_today = (24 - uva_hour_on) + hour

    # ===== 連續機制 1: 日內逆境 =====
    # factor = 1 + k × (hours)^n
    # 用較小的 k 和較高的 n 來讓短時間照射影響小，長時間照射影響大
    day_factor = 1.0 + p.k_day * (hours_today ** p.n_day)

    # ===== 連續機制 2: 夜間節律 =====
    # 計算進入暗期多久
    if light_on <= hour < light_off:
        hours_in_dark = 0
    else:
        if hour >= light_off:
            hours_in_dark = hour - light_off
        else:
            hours_in_dark = hour + (24 - light_off)

    # circ = 1 + k × hours_in_dark^n (只在有 UVA 時生效)
    if I_UVA > 0 and hours_in_dark > 0:
        circ = 1.0 + p.k_night * (hours_in_dark ** p.n_night)
    else:
        circ = 1.0

    # 環境
    env_mod = env.copy()
    env_mod['I_override'] = I_base + I_UVA * p.par_conversion_factor
    env_mod['T_override'] = Tc
    env_mod['is_day_override'] = is_day

    # Sun 基礎
    dXd_base, dCbuf, dLAI_base = sun_derivatives_final(t, [X_d, C_buf, LAI], p, env_mod)

    # ===== 連續機制 3: LAI 脆弱性 =====
    vuln = p.cap_vuln * (p.LAI_ref_vuln / LAI) ** p.n_vuln / (p.cap_vuln + (p.LAI_ref_vuln / LAI) ** p.n_vuln)

    # ===== 連續機制 4: Stress 非線性 =====
    nl = 1 + p.stress_nonlinear_coeff * Stress / (p.K_nonlinear + Stress + 1e-9)

    # 損傷 (所有因子相乘)
    damage = p.stress_damage_coeff * I_UVA * vuln * nl * circ * day_factor

    # 修復
    repair_cap = p.base_repair_capacity + p.carbon_repair_bonus * C_buf / (p.K_carbon + C_buf + 1e-9)
    repair = p.stress_repair_coeff * Stress * repair_cap
    dStress = damage - repair

    # Stress 抑制
    inh = Stress / (p.K_stress + Stress + 1e-9)
    dXd = dXd_base * (1 - p.stress_photosynthesis_inhibition * inh) if dXd_base > 0 else dXd_base
    dLAI = dLAI_base * (1 - p.stress_lai_inhibition * inh) if dLAI_base > 0 else dLAI_base

    # 花青素
    dw_fw = calc_dw_fw(Stress, p)
    FW = X_d / dw_fw
    synth = FW * ((p.base_anth_rate_light if is_day else p.base_anth_rate_dark) +
                  p.V_max_anth * Stress / (p.K_stress_anth + Stress + 1e-12))
    dAnth = synth - p.k_deg * Anth

    dCbuf = dCbuf - repair * p.repair_carbon_cost

    return np.array([dXd, dCbuf, dLAI, dAnth, dStress])


def test(params, verbose=False):
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

    # 計算誤差
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
            if verbose:
                s1 = 'ok' if abs(fw_err) < 5 else 'NG'
                s2 = 'ok' if abs(anth_err) < 10 else 'NG'
                print(f'{t:<8} FW:{FW:>5.1f}({fw_err:>+5.1f}%{s1}) Anth:{Anth:>5.1f}({anth_err:>+5.1f}%{s2}) LAI:{LAI:>4.1f} S:{Stress:>5.1f}')
        else:
            fw_errs.append(100)
            anth_errs.append(100)

    fw_ok = sum(1 for e in fw_errs if e < 5)
    anth_ok = sum(1 for e in anth_errs if e < 10)
    return fw_ok + anth_ok, results


# 基礎參數
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
    # 日內逆境 (連續)
    'k_day': 0.001,   # 很小的係數
    'n_day': 3.0,     # 冪次
    # 夜間節律 (連續)
    'k_night': 0.1,   # 係數
    'n_night': 2.0,   # 冪次
}

print('=' * 90)
print('純連續版本 v3 - 系統性參數搜索')
print('=' * 90)
print()
print('目標：')
print('  - L6D6 (6h/天×6天): FW ~91.5g, 輕微損傷')
print('  - H12D3 (12h/天×3天): FW ~60.6g, 較大損傷')
print('  - L6D6-N (夜間): 需要夜間節律來增加損傷')
print('  - L6D12 (6h/天×12天): 需要日數累積來增加損傷')
print()

# 分析問題：
# 12^n vs 6^n 的比例決定 H12D3 vs L6D6 的差異
# 如果 n=2: 144/36 = 4x
# 如果 n=3: 1728/216 = 8x
# 如果 n=4: 20736/1296 = 16x

print('分析日內因子比例 (12h vs 6h):')
for n in [1, 2, 3, 4, 5, 6]:
    ratio = (12**n) / (6**n)
    print(f'  n={n}: 12^{n}/6^{n} = {ratio:.1f}x')
print()

# 但還要考慮天數：H12D3 只有 3 天，L6D6 有 6 天
# 總累積：L6D6 = 6 × (1 + k×6^n), H12D3 = 3 × (1 + k×12^n)
print('考慮天數後的總累積因子:')
print('  L6D6: 6天 × day_factor(6h)')
print('  H12D3: 3天 × day_factor(12h)')
print()

# 測試不同 k_day 和 n_day 組合
print('搜索 k_day 和 n_day:')
print('-' * 90)

best_score = 0
best_params = None

for n_day in [2.0, 2.5, 3.0, 3.5, 4.0]:
    for k_day in [0.0001, 0.0005, 0.001, 0.002, 0.005, 0.01]:
        params = BASE.copy()
        params['k_day'] = k_day
        params['n_day'] = n_day
        score, results = test(params)

        # 計算關鍵誤差
        if 'L6D6' in results and 'H12D3' in results:
            fw_l6d6 = results['L6D6'][0]
            fw_h12d3 = results['H12D3'][0]
            err_l6d6 = (fw_l6d6 - 91.5) / 91.5 * 100
            err_h12d3 = (fw_h12d3 - 60.6) / 60.6 * 100
            print(f'n={n_day}, k={k_day:.4f}: {score}/12  L6D6:{err_l6d6:>+5.1f}% H12D3:{err_h12d3:>+5.1f}%')

            if score > best_score:
                best_score = score
                best_params = params.copy()

print()
print('=' * 90)
print(f'最佳: {best_score}/12')
if best_params:
    print(f"參數: k_day={best_params['k_day']}, n_day={best_params['n_day']}")
    print()
    print('詳細結果:')
    print('-' * 90)
    test(best_params, verbose=True)
