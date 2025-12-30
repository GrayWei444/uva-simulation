"""調校 v6.9 - 同時調整損傷係數和耗竭參數"""
import numpy as np
from scipy.integrate import solve_ivp
import sys
sys.path.insert(0, '/home/kasm-user/projects/uva-simulation')

from simulate_uva_model import ALL_PARAMS, UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TREATMENT_CONFIGS, TARGETS, SIMULATION

def test_params(params_override):
    """測試參數組合"""
    params = ALL_PARAMS.copy()
    params.update(params_override)
    p = UVAParams(params)

    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    C_buf_init = Xd_init * 0.1
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6
    initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0]

    transplant_day = SIMULATION['transplant_offset']
    simulation_days = SIMULATION['days']
    t_start = transplant_day * 86400
    t_end = (transplant_day + simulation_days) * 86400

    results = {}
    for treatment in ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']:
        config = TREATMENT_CONFIGS[treatment]
        env = ENV_BASE.copy()
        env.update(config)

        sol = solve_ivp(
            uva_sun_derivatives,
            (t_start, t_end),
            initial_state,
            args=(p, env),
            method='RK45',
            max_step=60,
            t_eval=np.array([t_end])
        )

        if sol.success:
            Xd_f, Cbuf_f, LAI_f, Anth_f, Stress_f = sol.y[:, -1]
            dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress_f, p)
            X_total = Xd_f / ENV_BASE['plant_density']
            FW_sim = X_total / dw_fw_ratio * 1000
            FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
            Anth_sim = (Anth_f / FW_total_kg) * 1e6

            FW_obs = TARGETS[treatment]['FW']
            Anth_obs = TARGETS[treatment].get('Anth', 0)

            results[treatment] = {
                'FW_sim': FW_sim, 'FW_obs': FW_obs,
                'Anth_sim': Anth_sim, 'Anth_obs': Anth_obs,
                'Stress': Stress_f
            }
    return results

def evaluate(results):
    """評估結果"""
    fw_errs = []
    anth_errs = []
    for t in ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']:
        r = results[t]
        fw_err = abs((r['FW_sim'] - r['FW_obs']) / r['FW_obs']) * 100
        anth_err = abs((r['Anth_sim'] - r['Anth_obs']) / r['Anth_obs']) * 100 if r['Anth_obs'] > 0 else 0
        fw_errs.append(fw_err)
        anth_errs.append(anth_err)
    return fw_errs, anth_errs

def print_results(results, label, show_full=True):
    print(f"\n{label}")
    print("-" * 100)
    fw_errs, anth_errs = evaluate(results)
    for i, t in enumerate(['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']):
        r = results[t]
        fw_err = (r['FW_sim'] - r['FW_obs']) / r['FW_obs'] * 100
        anth_err = (r['Anth_sim'] - r['Anth_obs']) / r['Anth_obs'] * 100 if r['Anth_obs'] > 0 else 0
        if show_full:
            print(f"{t:<8} FW:{r['FW_sim']:>6.1f}g({fw_err:>+5.1f}%)  Anth:{r['Anth_sim']:>5.1f}({anth_err:>+5.1f}%)  S:{r['Stress']:>5.2f}")
    fw_ok = sum(1 for e in fw_errs if e < 5)
    anth_ok = sum(1 for e in anth_errs if e < 10)
    print(f"FW <5%: {fw_ok}/6, Anth <10%: {anth_ok}/6, 總達標: {fw_ok + anth_ok}/12")
    return fw_ok + anth_ok

# 網格搜索
print("=" * 100)
print("v6.9 參數聯合調校")
print("=" * 100)

best_score = 0
best_params = None

# 策略: 增大 stress_damage_coeff，用較大的 K_repair_capacity
for damage_coeff in [0.5e-6, 0.66e-6, 0.8e-6, 1.0e-6, 1.2e-6, 1.5e-6]:
    for K in [200, 300, 400, 500, 600]:
        for n in [2.0, 2.5, 3.0, 3.5]:
            params = {
                'stress_damage_coeff': damage_coeff,
                'K_repair_capacity': K,
                'repair_depletion_power': n
            }
            try:
                results = test_params(params)
                fw_errs, anth_errs = evaluate(results)
                fw_ok = sum(1 for e in fw_errs if e < 5)
                anth_ok = sum(1 for e in anth_errs if e < 10)
                score = fw_ok + anth_ok

                if score > best_score:
                    best_score = score
                    best_params = params.copy()
                    print(f"d={damage_coeff:.2e}, K={K}, n={n:.1f} → FW:{fw_ok}/6, Anth:{anth_ok}/6, Score:{score}/12 *NEW BEST*")
            except:
                pass

print(f"\n最佳參數: {best_params}")
print(f"最佳分數: {best_score}/12")

if best_params:
    results = test_params(best_params)
    print_results(results, f"最佳結果: damage={best_params['stress_damage_coeff']:.2e}, K={best_params['K_repair_capacity']}, n={best_params['repair_depletion_power']}")
