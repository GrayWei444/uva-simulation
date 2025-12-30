"""驗證 v7.0 結果"""
import numpy as np
from scipy.integrate import solve_ivp
from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment
from simulate_uva_model_v7 import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio

p = UVAParams()

print("v7.0 模型驗證")
print("=" * 80)

fw_errs = []
anth_errs = []

for treatment in ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']:
    env = get_env_for_treatment(treatment)
    target = TARGETS[treatment]

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
        max_step=300,  # 精確模擬
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

        fw_errs.append(abs(fw_err))
        anth_errs.append(abs(anth_err))

        s1 = "✓" if abs(fw_err) < 5 else "✗"
        s2 = "✓" if abs(anth_err) < 5 else ("△" if abs(anth_err) < 10 else "✗")
        print(f"{treatment:<8} FW: {FW_sim:>5.1f}g (目標{target['FW']:.1f}, {fw_err:>+5.1f}%{s1}) "
              f"Anth: {Anth_sim:>5.1f} (目標{target['Anth']:.1f}, {anth_err:>+5.1f}%{s2})")

print("-" * 80)
fw_ok_5 = sum(1 for e in fw_errs if e < 5)
anth_ok_5 = sum(1 for e in anth_errs if e < 5)
anth_ok_10 = sum(1 for e in anth_errs if e < 10)

print(f"\nFW 達標: {fw_ok_5}/6 (<5%)")
print(f"Anth 達標: {anth_ok_5}/6 (<5%), {anth_ok_10}/6 (<10%)")
print(f"\n平均誤差: FW={np.mean(fw_errs):.1f}%, Anth={np.mean(anth_errs):.1f}%")
print(f"最大誤差: FW={max(fw_errs):.1f}%, Anth={max(anth_errs):.1f}%")
