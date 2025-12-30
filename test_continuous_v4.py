"""純連續版本 v4 - 增加日內累積損傷機制

問題：H12D3 (12h×3天) 損傷不夠
原因：day_factor = 1 + k × hours^n 只考慮當前時刻的照射時間
     但 H12D3 每天照射 12 小時，累積能量比 L6D6 (6h×6天) 更集中

解決方案：
1. 日內損傷累積：隨著當天照射時間增加，損傷效率上升
   - 不是用 hours^n，而是累積能量的效應
   - 模擬 ROS 累積超過清除能力的效應

2. 每小時的瞬時損傷 = base × (1 + 累積放大係數)
   累積放大係數 = k × ∫[0,hours] f(h) dh

更簡單的實現：用 hours^(n+1) 模擬累積效應
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


def softplus(x, k=1.0):
    """平滑的 max(0, x) 函數"""
    return np.log(1 + np.exp(np.clip(k * x, -500, 500))) / k


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
            if hour >= uva_hour_on:
                if uva_start <= day <= uva_end:
                    I_UVA = env['uva_intensity']
                    hours_today = hour - uva_hour_on
            elif hour < uva_hour_off:
                if uva_start <= day - 1 <= uva_end:
                    I_UVA = env['uva_intensity']
                    hours_today = (24 - uva_hour_on) + hour

    # ===== 連續機制 1: 日內累積逆境 =====
    # 問題：12h×3天 的總照射時間 = 6h×6天，但實驗顯示 12h×3天 損傷更大
    # 說明：每天長時間照射會造成累積毒性（ROS 超過清除能力）
    #
    # 模型：瞬時損傷 = base × (1 + k × hours^n)
    # 當 hours=6: factor = 1 + k × 6^n
    # 當 hours=12: factor = 1 + k × 12^n
    #
    # 但這只是瞬時值。真正的日累積損傷是積分：
    # ∫[0,H] (1 + k × h^n) dh = H + k × H^(n+1)/(n+1)
    #
    # 對於 L6D6 (6天): 總 = 6 × [6 + k × 6^(n+1)/(n+1)]
    # 對於 H12D3 (3天): 總 = 3 × [12 + k × 12^(n+1)/(n+1)]
    #
    # 比較 (假設 n=2, k=0.01):
    # L6D6: 6 × [6 + 0.01 × 216/3] = 6 × [6 + 0.72] = 40.32
    # H12D3: 3 × [12 + 0.01 × 1728/3] = 3 × [12 + 5.76] = 53.28
    # 比例: 53.28/40.32 = 1.32x
    #
    # 如果用 n=3:
    # L6D6: 6 × [6 + 0.01 × 1296/4] = 6 × [6 + 3.24] = 55.44
    # H12D3: 3 × [12 + 0.01 × 20736/4] = 3 × [12 + 51.84] = 191.52
    # 比例: 191.52/55.44 = 3.45x

    # 使用瞬時值：factor = 1 + k × hours^n
    # 這會透過 ODE 積分自然累積
    day_factor = 1.0 + p.k_day * (hours_today ** p.n_day)

    # ===== 連續機制 2: 夜間節律 =====
    if light_on <= hour < light_off:
        hours_in_dark = 0
    else:
        if hour >= light_off:
            hours_in_dark = hour - light_off
        else:
            hours_in_dark = hour + (24 - light_off)

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

    # 損傷
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
    'k_day': 0.01,
    'n_day': 3.0,
    'k_night': 0.1,
    'n_night': 2.0,
}

print('=' * 90)
print('純連續版本 v4 - 調整 n_day 以增加 H12D3 損傷')
print('=' * 90)
print()

# 計算理論上需要的 n_day
# 目標：
# L6D6: FW 91.5g (幾乎無損傷)
# H12D3: FW 60.6g (約 30% 損傷)
#
# 日積分損傷比例：
# H12D3/L6D6 = 3 × ∫[0,12] (1+k×h^n) dh / 6 × ∫[0,6] (1+k×h^n) dh
#            = 3 × [12 + k×12^(n+1)/(n+1)] / 6 × [6 + k×6^(n+1)/(n+1)]
#            = [12 + k×12^(n+1)/(n+1)] / [12 + 2k×6^(n+1)/(n+1)]

print('理論日積分損傷比例 (H12D3/L6D6) vs n_day:')
k = 0.01
for n in [2, 3, 4, 5, 6]:
    int_12 = 12 + k * (12**(n+1)) / (n+1)
    int_6 = 6 + k * (6**(n+1)) / (n+1)
    # 天數加權
    total_h12d3 = 3 * int_12
    total_l6d6 = 6 * int_6
    ratio = total_h12d3 / total_l6d6
    print(f'  n={n}: H12D3積分={total_h12d3:.1f}, L6D6積分={total_l6d6:.1f}, 比例={ratio:.2f}x')

print()
print('搜索 n_day (高冪次讓長時間照射損傷更大):')
print('-' * 90)

best_score = 0
best_params = None

for n_day in [4.0, 5.0, 6.0, 7.0, 8.0]:
    for k_day in [0.0001, 0.0005, 0.001, 0.002, 0.003, 0.005]:
        params = BASE.copy()
        params['k_day'] = k_day
        params['n_day'] = n_day
        score, results = test(params)

        if 'L6D6' in results and 'H12D3' in results:
            fw_l6d6 = results['L6D6'][0]
            fw_h12d3 = results['H12D3'][0]
            err_l6d6 = (fw_l6d6 - 91.5) / 91.5 * 100
            err_h12d3 = (fw_h12d3 - 60.6) / 60.6 * 100
            mark = "***" if score >= 11 else ""
            print(f'n={n_day}, k={k_day:.4f}: {score:>2}/12  L6D6:{err_l6d6:>+6.1f}% H12D3:{err_h12d3:>+6.1f}% {mark}')

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
