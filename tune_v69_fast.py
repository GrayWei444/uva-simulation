"""快速調校 v6.9 - 精簡版"""
import numpy as np
from scipy.integrate import solve_ivp
import sys
sys.path.insert(0, '/home/kasm-user/projects/uva-simulation')

from simulate_uva_model import ALL_PARAMS, UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TREATMENT_CONFIGS, TARGETS, SIMULATION

def test_params(params_override):
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
    # 只測試 CK, L6D6, H12D3 三個關鍵組
    for treatment in ['CK', 'L6D6', 'H12D3']:
        config = TREATMENT_CONFIGS[treatment]
        env = ENV_BASE.copy()
        env.update(config)

        sol = solve_ivp(
            uva_sun_derivatives, (t_start, t_end), initial_state,
            args=(p, env), method='RK45', max_step=60, t_eval=np.array([t_end])
        )

        if sol.success:
            Xd_f, Cbuf_f, LAI_f, Anth_f, Stress_f = sol.y[:, -1]
            dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress_f, p)
            X_total = Xd_f / ENV_BASE['plant_density']
            FW_sim = X_total / dw_fw_ratio * 1000
            results[treatment] = {'FW_sim': FW_sim, 'Stress': Stress_f}
    return results

print("=" * 80)
print("v6.9 快速參數搜索 (只測試 CK, L6D6, H12D3)")
print("=" * 80)

# 目標: CK ~87g, L6D6 ~91g, H12D3 ~60g
# L6D6 應該比 CK 略高 (有益效果)
# H12D3 應該比 L6D6 低約 33%

best_score = float('inf')
best_params = None

for damage in [0.4e-6, 0.5e-6, 0.66e-6, 0.8e-6, 1.0e-6]:
    for K in [100, 150, 200, 300, 400, 500]:
        for n in [2.0, 3.0, 4.0]:
            params = {'stress_damage_coeff': damage, 'K_repair_capacity': K, 'repair_depletion_power': n}
            try:
                r = test_params(params)
                # 評分
                ck_err = abs(r['CK']['FW_sim'] - 87.0) / 87.0 * 100
                l6d6_err = abs(r['L6D6']['FW_sim'] - 91.4) / 91.4 * 100
                h12d3_err = abs(r['H12D3']['FW_sim'] - 60.6) / 60.6 * 100

                score = ck_err + l6d6_err + h12d3_err * 2  # H12D3 權重加倍

                if score < best_score:
                    best_score = score
                    best_params = params
                    print(f"d={damage:.1e} K={K:>3} n={n:.1f} | CK:{r['CK']['FW_sim']:>5.1f}({ck_err:>4.1f}%) "
                          f"L6D6:{r['L6D6']['FW_sim']:>5.1f}({l6d6_err:>4.1f}%) "
                          f"H12D3:{r['H12D3']['FW_sim']:>5.1f}({h12d3_err:>4.1f}%) Score:{score:>5.1f} *")
            except Exception as e:
                pass

print(f"\n最佳: {best_params}")
