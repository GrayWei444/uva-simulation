"""v6.9 參數調校腳本"""
import numpy as np
from scipy.integrate import solve_ivp
import sys
sys.path.insert(0, '/home/kasm-user/projects/uva-simulation')

from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio, ALL_PARAMS
from model_config import ENV_BASE, TREATMENT_CONFIGS, TARGETS, SIMULATION

def run_simulation(params_override=None):
    """運行模擬並返回結果"""
    # 合併參數
    params = ALL_PARAMS.copy()
    if params_override:
        params.update(params_override)

    # 初始條件
    fw_init_g = SIMULATION['initial_fw_g']
    p = UVAParams(params)
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    C_buf_init = Xd_init * 0.1
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6
    initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0]

    # 模擬時間
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
            Anth_conc_sim = (Anth_f / FW_total_kg) * 1e6

            FW_obs = TARGETS[treatment]['FW']
            Anth_obs = TARGETS[treatment].get('Anth', 0)

            results[treatment] = {
                'FW_sim': FW_sim,
                'FW_obs': FW_obs,
                'FW_err': abs((FW_sim - FW_obs) / FW_obs * 100),
                'Anth_sim': Anth_conc_sim,
                'Anth_obs': Anth_obs,
                'Anth_err': abs((Anth_conc_sim - Anth_obs) / Anth_obs * 100) if Anth_obs > 0 else 0,
                'Stress': Stress_f
            }
        else:
            results[treatment] = None

    return results

def evaluate_params(params_override):
    """評估參數組合"""
    results = run_simulation(params_override)
    if any(r is None for r in results.values()):
        return float('inf'), float('inf')

    fw_errors = [r['FW_err'] for r in results.values()]
    anth_errors = [r['Anth_err'] for r in results.values()]

    return np.mean(fw_errors), np.mean(anth_errors)

# 參數掃描
print("=" * 80)
print("v6.9 參數調校")
print("=" * 80)

# 首先測試不同的損傷係數和修復飽和常數
best_params = None
best_score = float('inf')

# 損傷係數範圍
damage_coefs = [1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5]
# 修復飽和常數範圍
k_saturations = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
# 能量放大係數
k_energies = [0.0005, 0.001, 0.002, 0.005, 0.01]

print("\n掃描損傷係數...")
for damage in damage_coefs:
    params = {'stress_damage_coeff': damage}
    fw_err, anth_err = evaluate_params(params)
    print(f"  damage={damage:.1e}: FW_err={fw_err:.1f}%, Anth_err={anth_err:.1f}%")

print("\n掃描修復飽和常數...")
for k_sat in k_saturations:
    params = {'K_repair_saturation': k_sat}
    fw_err, anth_err = evaluate_params(params)
    print(f"  K_sat={k_sat}: FW_err={fw_err:.1f}%, Anth_err={anth_err:.1f}%")

print("\n掃描能量放大係數...")
for k_e in k_energies:
    params = {'k_energy_amplify': k_e}
    fw_err, anth_err = evaluate_params(params)
    print(f"  k_energy={k_e}: FW_err={fw_err:.1f}%, Anth_err={anth_err:.1f}%")

# 測試組合
print("\n測試參數組合...")
for damage in [1e-5, 2e-5, 5e-5]:
    for k_sat in [0.5, 1.0, 2.0]:
        for k_e in [0.002, 0.005]:
            params = {
                'stress_damage_coeff': damage,
                'K_repair_saturation': k_sat,
                'k_energy_amplify': k_e
            }
            fw_err, anth_err = evaluate_params(params)
            score = fw_err + anth_err
            if score < best_score:
                best_score = score
                best_params = params.copy()
            if fw_err < 10 and anth_err < 15:
                print(f"  damage={damage:.1e}, K_sat={k_sat}, k_e={k_e}: "
                      f"FW={fw_err:.1f}%, Anth={anth_err:.1f}%")

print(f"\n最佳參數: {best_params}")
print(f"最佳分數: FW+Anth = {best_score:.1f}%")

# 詳細結果
if best_params:
    print("\n" + "=" * 80)
    print("最佳參數的詳細結果:")
    print("=" * 80)
    results = run_simulation(best_params)
    for treatment, r in results.items():
        if r:
            print(f"{treatment}: FW={r['FW_sim']:.1f}g ({r['FW_err']:+.1f}%), "
                  f"Anth={r['Anth_sim']:.1f} ({r['Anth_err']:+.1f}%), "
                  f"Stress={r['Stress']:.2f}")
