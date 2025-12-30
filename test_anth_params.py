"""測試不同花青素參數"""

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

    return results


print("=" * 70)
print("測試不同花青素參數")
print("=" * 70)

# 測試 v6.9 原始參數
print("\n1. v6.9 原始參數 (base=2e-10, V=2.35e-11, K=0.30):")
results = test_all(2.0e-10, 2.35e-11, 0.30)
for t, anth, err, s in results:
    c = "✓" if abs(err) < 5 else ("△" if abs(err) < 10 else "✗")
    print(f"  {t:<8}: {anth:.1f} ({err:+.1f}%{c}) S={s:.1f}")

# 測試目前 v7.0 優化參數
print("\n2. v7.0 優化參數 (base=1.94e-10, V=2.71e-11, K=1.38):")
results = test_all(1.94e-10, 2.71e-11, 1.38)
for t, anth, err, s in results:
    c = "✓" if abs(err) < 5 else ("△" if abs(err) < 10 else "✗")
    print(f"  {t:<8}: {anth:.1f} ({err:+.1f}%{c}) S={s:.1f}")

# 策略: 先讓 CK 和 L6D6 準確
# CK 目標 43.3, L6D6 目標 49.4
# 差距 = 49.4 - 43.3 = 6.1 ppm (+14%)
# L6D6 的 Stress = 0.3
# 誘導 = V × 0.3 / (K + 0.3)
# 需要讓誘導產生 6.1 ppm 的增量

print("\n3. 策略: 先確保 CK 和 L6D6 準確")
print("   CK 只有 base，L6D6 有 base + Stress 誘導")
print("   目標差距: 49.4 - 43.3 = 6.1 ppm")

# 如果 CK 用 base=2e-10 得到 47.0 (+8.5%)
# 需要 base × 0.921 = 1.84e-10 來讓 CK = 43.3
# 但這會讓 L6D6 也下降...

# 讓我測試一下 v6.9 的 Stress 值
print("\n   v6.9 的 Stress 值:")
for t, anth, err, s in test_all(2.0e-10, 2.35e-11, 0.30):
    print(f"     {t}: S={s:.2f}")
