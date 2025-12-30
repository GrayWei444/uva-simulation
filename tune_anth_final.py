"""花青素參數最終微調 - 目標: 降低 CK 和 L6D6 誤差"""

import numpy as np
from scipy.integrate import solve_ivp
from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment
from simulate_uva_model_v7 import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio, ALL_PARAMS


def test_all(base_light, v_max, k_anth):
    params = ALL_PARAMS.copy()
    params['base_anth_rate_light'] = base_light
    params['base_anth_rate_dark'] = base_light / 2
    params['V_max_anth'] = v_max
    params['K_stress_anth'] = k_anth
    p = UVAParams(params)
    results = []

    for treatment in ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']:
        env = get_env_for_treatment(treatment)
        target = TARGETS[treatment]
        fw_init_g = SIMULATION['initial_fw_g']
        dw_init_g = fw_init_g * p.dw_fw_ratio_base
        Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
        initial_state = [Xd_init, Xd_init * 0.1, (dw_init_g / 0.01) * 0.04,
                         5.0 * fw_init_g * ENV_BASE['plant_density'] / 1000 / 1e6, 0.0]
        t_start = SIMULATION['transplant_offset'] * 86400
        t_end = (SIMULATION['transplant_offset'] + SIMULATION['days']) * 86400

        sol = solve_ivp(uva_sun_derivatives, (t_start, t_end), initial_state,
                        args=(p, env), method='RK45', max_step=300, t_eval=[t_end])

        if sol.success:
            Xd_f, _, _, Anth_f, Stress_f = sol.y[:, -1]
            dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress_f, p)
            FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
            FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
            Anth_sim = Anth_f / FW_total_kg * 1e6
            anth_err = (Anth_sim - target['Anth']) / target['Anth'] * 100
            results.append((treatment, Anth_sim, anth_err, Stress_f))
        else:
            results.append((treatment, 0, 100, 0))

    return results


print("花青素參數最終微調")
print("=" * 60)

# 當前: CK +5.2%, L6D6 -5.9%
# CK 太高 → 降低 base
# L6D6 太低 → 提高 V 或降低 K

# 問題: CK 和 L6D6 的誤差方向相反！
# CK (Stress=0): Anth = base 效應
# L6D6 (Stress=0.3): Anth = base + V*0.3/(K+0.3)
#
# 降低 base → CK 降低 ✓，L6D6 也降低 ✗
# 提高 V → CK 不變，L6D6 提高 ✓
# 降低 K → CK 不變，L6D6 提高 ✓

# 策略: 微降 base，提高 V 或降低 K

print("\n當前參數: base=1.94e-10, V=2.71e-11, K=1.38")
results = test_all(1.94e-10, 2.71e-11, 1.38)
for t, anth, err, s in results:
    print(f"  {t:<8}: {err:+.1f}%")

print("\n微調搜索...")
best = (100, None, None)

# 在當前最佳點附近搜索
for b in np.linspace(1.88e-10, 1.96e-10, 9):
    for v in np.linspace(2.5e-11, 3.2e-11, 9):
        for k in np.linspace(0.8, 1.6, 9):
            results = test_all(b, v, k)
            max_err = max(abs(r[2]) for r in results)
            if max_err < best[0]:
                best = (max_err, (b, v, k), results)

print(f"\n最佳 max_err: {best[0]:.2f}%")

if best[1]:
    b, v, k = best[1]
    print(f"\n最佳參數:")
    print(f"  'base_anth_rate_light': {b:.2e},")
    print(f"  'base_anth_rate_dark': {b/2:.2e},")
    print(f"  'V_max_anth': {v:.2e},")
    print(f"  'K_stress_anth': {k:.2f},")

    print(f"\n詳細結果:")
    results = test_all(b, v, k)
    n_ok = 0
    for t, anth, err, s in results:
        target = TARGETS[t]['Anth']
        c = "✓" if abs(err) < 5 else "✗"
        if abs(err) < 5:
            n_ok += 1
        print(f"  {t:<8}: 目標={target:.1f} 模擬={anth:.1f} ({err:+.1f}%{c})")

    print(f"\n達標: {n_ok}/6 (<5%)")
    if n_ok == 6:
        print("✅ 完美！所有誤差 < 5%！")
