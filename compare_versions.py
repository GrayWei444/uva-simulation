"""比較 v6.9 和 v7.0 的花青素結果"""
import numpy as np
from scipy.integrate import solve_ivp

from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment

# v6.9
from simulate_uva_model import UVAParams as V69Params, uva_sun_derivatives as v69_deriv, calculate_dynamic_dw_fw_ratio as v69_dw_fw

# v7.0
from simulate_uva_model_v7 import UVAParams as V70Params, uva_sun_derivatives as v70_deriv, calculate_dynamic_dw_fw_ratio as v70_dw_fw


def run_model(deriv_func, params_class, dw_fw_func, label):
    p = params_class()
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
            deriv_func,
            (t_start, t_end),
            initial_state,
            args=(p, env),
            method='RK45',
            max_step=300,
            t_eval=np.array([t_end])
        )

        if sol.success:
            Xd_f, Cbuf_f, LAI_f, Anth_f, Stress_f = sol.y[:, -1]
            dw_fw_ratio = dw_fw_func(Stress_f, p)
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
                'target_anth': target['Anth']
            })

    return results


print("=" * 90)
print("v6.9 vs v7.0 花青素比較")
print("=" * 90)

v69_results = run_model(v69_deriv, V69Params, v69_dw_fw, "v6.9")
v70_results = run_model(v70_deriv, V70Params, v70_dw_fw, "v7.0")

print(f"\n{'處理組':<10} {'目標':>6} │ {'v6.9 Anth':>10} {'誤差':>8} │ {'v7.0 Anth':>10} {'誤差':>8} │ {'差異':>8}")
print("-" * 90)

for r69, r70 in zip(v69_results, v70_results):
    t = r69['treatment']
    target = r69['target_anth']
    a69 = r69['Anth']
    e69 = r69['anth_err']
    a70 = r70['Anth']
    e70 = r70['anth_err']
    diff = a70 - a69

    s69 = "✓" if abs(e69) < 10 else "✗"
    s70 = "✓" if abs(e70) < 10 else "✗"

    print(f"{t:<10} {target:>6.1f} │ {a69:>8.1f} {e69:>+6.1f}%{s69} │ {a70:>8.1f} {e70:>+6.1f}%{s70} │ {diff:>+7.1f}")

print("-" * 90)

# 計算平均絕對誤差
mae_69 = np.mean([abs(r['anth_err']) for r in v69_results])
mae_70 = np.mean([abs(r['anth_err']) for r in v70_results])
print(f"\n平均絕對誤差 (MAE):")
print(f"  v6.9: {mae_69:.2f}%")
print(f"  v7.0: {mae_70:.2f}%")

# 也比較 Stress
print(f"\n{'處理組':<10} │ {'v6.9 Stress':>12} │ {'v7.0 Stress':>12} │ {'差異':>10}")
print("-" * 60)
for r69, r70 in zip(v69_results, v70_results):
    t = r69['treatment']
    s69 = r69['Stress']
    s70 = r70['Stress']
    diff = s70 - s69
    print(f"{t:<10} │ {s69:>12.2f} │ {s70:>12.2f} │ {diff:>+10.2f}")
