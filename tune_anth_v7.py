"""調整 v7.0 花青素參數，目標 <5% 誤差"""
import numpy as np
from scipy.integrate import solve_ivp

from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment
from simulate_uva_model_v7 import ALL_PARAMS, UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio


def test_params(params_override):
    params = ALL_PARAMS.copy()
    params.update(params_override)

    class P(UVAParams):
        def __init__(self):
            super().__init__(params)

    p = P()
    results = []

    for treatment in ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']:
        env = get_env_for_treatment(treatment)
        target = TARGETS.get(treatment, {'FW': 0, 'Anth': 0})

        fw_init_g = SIMULATION['initial_fw_g']
        dw_init_g = fw_init_g * p.dw_fw_ratio_base
        Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
        C_buf_init = Xd_init * 0.1
        LAI_init = (dw_init_g / 0.01) * 0.04
        fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
        Anth_init = 5.0 * fw_total_init / 1e6
        initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0]

        t_start = SIMULATION['transplant_offset'] * 86400
        t_end = (SIMULATION['transplant_offset'] + SIMULATION['days']) * 86400

        sol = solve_ivp(
            uva_sun_derivatives,
            (t_start, t_end),
            initial_state,
            args=(p, env),
            method='RK45',
            max_step=300,
            t_eval=np.array([t_end])
        )

        if sol.success:
            Xd_f, Cbuf_f, LAI_f, Anth_f, Stress_f = sol.y[:, -1]
            dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress_f, p)
            FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
            FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
            Anth_sim = Anth_f / FW_total_kg * 1e6

            fw_err = (FW_sim - target['FW']) / target['FW'] * 100
            anth_err = (Anth_sim - target['Anth']) / target['Anth'] * 100

            results.append({
                'treatment': treatment,
                'FW': FW_sim,
                'Anth': Anth_sim,
                'Stress': Stress_f,
                'fw_err': fw_err,
                'anth_err': anth_err,
                'target_fw': target['FW'],
                'target_anth': target['Anth']
            })

    return results


def score_results(results):
    fw_ok = sum(1 for r in results if abs(r['fw_err']) < 5)
    anth_ok = sum(1 for r in results if abs(r['anth_err']) < 5)
    max_fw = max(abs(r['fw_err']) for r in results)
    max_anth = max(abs(r['anth_err']) for r in results)
    return fw_ok, anth_ok, max_fw, max_anth


def print_results(results, label=""):
    if label:
        print(f"\n{label}")
        print("-" * 90)

    for r in results:
        t = r['treatment']
        fw = r['FW']
        anth = r['Anth']
        fe = r['fw_err']
        ae = r['anth_err']
        s = r['Stress']
        s1 = "✓" if abs(fe) < 5 else "✗"
        s2 = "✓" if abs(ae) < 5 else "✗"
        print(f"{t:<8} FW:{fw:>5.1f}({fe:>+5.1f}%{s1}) Anth:{anth:>5.1f}({ae:>+5.1f}%{s2}) S:{s:>5.1f}")

    fw_ok, anth_ok, max_fw, max_anth = score_results(results)
    print(f"達標: FW {fw_ok}/6 (<5%), Anth {anth_ok}/6 (<5%)")
    print(f"最大誤差: FW {max_fw:.1f}%, Anth {max_anth:.1f}%")


print("=" * 90)
print("調整 v7.0 花青素參數")
print("目標: 所有花青素誤差 < 5%")
print("=" * 90)

# 目前結果
print("\n當前參數:")
results = test_params({})
print_results(results)

# 問題分析:
# CK: +8.5% (過高) - 需要降低基礎合成率
# VL3D12: +9.4% (過高) - 需要降低 Stress 誘導合成
# L6D6-N: +6.8% (過高)

# 策略 1: 降低基礎合成率
print("\n" + "=" * 90)
print("策略 1: 降低基礎合成率 (影響 CK)")
print("=" * 90)

for base_light in [1.8e-10, 1.6e-10, 1.4e-10]:
    base_dark = base_light / 2
    results = test_params({
        'base_anth_rate_light': base_light,
        'base_anth_rate_dark': base_dark
    })
    fw_ok, anth_ok, max_fw, max_anth = score_results(results)
    ck_anth = results[0]['anth_err']
    print(f"base_light={base_light:.1e}: CK={ck_anth:+.1f}%, max_anth={max_anth:.1f}%, anth_ok={anth_ok}/6")

# 策略 2: 調整 V_max_anth
print("\n" + "=" * 90)
print("策略 2: 調整 V_max_anth (影響 UVA 處理組)")
print("=" * 90)

for v_max in [2.0e-11, 2.2e-11, 2.35e-11, 2.5e-11]:
    results = test_params({
        'V_max_anth': v_max
    })
    fw_ok, anth_ok, max_fw, max_anth = score_results(results)
    h12d3_anth = results[5]['anth_err']
    vl3d12_anth = results[3]['anth_err']
    print(f"V_max={v_max:.2e}: H12D3={h12d3_anth:+.1f}%, VL3D12={vl3d12_anth:+.1f}%, max_anth={max_anth:.1f}%")

# 策略 3: 調整 K_stress_anth
print("\n" + "=" * 90)
print("策略 3: 調整 K_stress_anth (影響 Stress 誘導靈敏度)")
print("=" * 90)

for k_anth in [0.2, 0.25, 0.30, 0.35, 0.4]:
    results = test_params({
        'K_stress_anth': k_anth
    })
    fw_ok, anth_ok, max_fw, max_anth = score_results(results)
    print(f"K_stress_anth={k_anth}: max_anth={max_anth:.1f}%, anth_ok={anth_ok}/6")

# 策略 4: 組合調整
print("\n" + "=" * 90)
print("策略 4: 組合調整")
print("=" * 90)

best_score = 0
best_params = {}

for base_light in [1.4e-10, 1.5e-10, 1.6e-10, 1.7e-10]:
    for v_max in [2.0e-11, 2.1e-11, 2.2e-11, 2.3e-11]:
        for k_anth in [0.25, 0.30, 0.35]:
            params = {
                'base_anth_rate_light': base_light,
                'base_anth_rate_dark': base_light / 2,
                'V_max_anth': v_max,
                'K_stress_anth': k_anth
            }
            results = test_params(params)
            fw_ok, anth_ok, max_fw, max_anth = score_results(results)

            # 目標: 全部 <5%
            if anth_ok == 6 and max_anth < 5:
                score = 100 - max_anth
                if score > best_score:
                    best_score = score
                    best_params = params.copy()
                    print(f"找到! base={base_light:.1e}, V_max={v_max:.1e}, K={k_anth}: max_anth={max_anth:.1f}%")

if best_params:
    print("\n" + "=" * 90)
    print("最佳參數組合:")
    print("=" * 90)
    for k, v in best_params.items():
        print(f"  {k}: {v}")
    results = test_params(best_params)
    print_results(results, "最佳結果")
else:
    print("\n未找到全部 <5% 的組合，顯示最接近的結果:")
    # 嘗試更廣的範圍
    for base_light in [1.2e-10, 1.3e-10, 1.4e-10, 1.5e-10]:
        for v_max in [1.8e-11, 1.9e-11, 2.0e-11, 2.1e-11, 2.2e-11]:
            for k_anth in [0.20, 0.25, 0.30, 0.35, 0.40]:
                params = {
                    'base_anth_rate_light': base_light,
                    'base_anth_rate_dark': base_light / 2,
                    'V_max_anth': v_max,
                    'K_stress_anth': k_anth
                }
                results = test_params(params)
                fw_ok, anth_ok, max_fw, max_anth = score_results(results)

                if anth_ok >= 5 and max_anth < 6:
                    print(f"base={base_light:.1e}, V_max={v_max:.1e}, K={k_anth}: anth_ok={anth_ok}/6, max={max_anth:.1f}%")
