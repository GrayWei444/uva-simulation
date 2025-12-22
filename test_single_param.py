"""
測試單一參數組合
"""
import sys
import numpy as np
from scipy.integrate import solve_ivp

from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TARGETS, ODE_SETTINGS, SIMULATION, get_env_for_treatment

# 從命令列獲取參數
if len(sys.argv) != 4:
    print("用法: python3 test_single_param.py V_max K_E n_hill")
    sys.exit(1)

V_max = float(sys.argv[1])
K_E = float(sys.argv[2])
n_hill = float(sys.argv[3])

p = UVAParams()
p.V_max_anth = V_max
p.K_E_anth = K_E
p.n_hill_anth = n_hill

print(f"測試參數: V_max={V_max:.2e}, K_E={K_E:.1f}, n_hill={n_hill:.1f}\n")

fw_errors = []
anth_errors = []

for treatment in ['CK', 'L6D6', 'L6D6-N', 'H12D3', 'VL3D12', 'L6D12']:
    env = get_env_for_treatment(treatment)
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6
    initial_state = [Xd_init, 0.0, LAI_init, Anth_init, 0.0, 0.0]

    sol = solve_ivp(
        fun=uva_sun_derivatives,
        t_span=(0, SIMULATION['days'] * 86400),
        y0=initial_state,
        args=(p, env),
        method=ODE_SETTINGS['method'],
        max_step=ODE_SETTINGS['max_step']
    )

    if sol.success:
        target = TARGETS[treatment]
        sim_stress = sol.y[4, -1]
        sim_E_stress = sol.y[5, -1]
        sim_dw_fw_ratio = calculate_dynamic_dw_fw_ratio(sim_stress, p)
        sim_dw_per_plant = sol.y[0, -1] / ENV_BASE['plant_density'] * 1000
        sim_fw_per_plant = sim_dw_per_plant / sim_dw_fw_ratio
        sim_anth_kg_m2 = sol.y[3, -1]
        sim_fw_kg_m2 = sim_fw_per_plant * ENV_BASE['plant_density'] / 1000
        sim_anth_ppm = sim_anth_kg_m2 * 1e6 / (sim_fw_kg_m2 + 1e-9)

        fw_err = (sim_fw_per_plant - target['FW']) / target['FW'] * 100
        anth_err = (sim_anth_ppm - target['Anth']) / target['Anth'] * 100

        fw_errors.append(abs(fw_err))
        anth_errors.append(abs(anth_err))

        print(f"{treatment:<10} FW:{sim_fw_per_plant:6.1f}g ({fw_err:+6.1f}%)  "
              f"Anth:{sim_anth_ppm:6.1f}ppm ({anth_err:+6.1f}%)  "
              f"Stress:{sim_stress:5.1f}  E_stress:{sim_E_stress:6.1f}")

print(f"\nMean |Err|: FW={np.mean(fw_errors):.1f}%, Anth={np.mean(anth_errors):.1f}%")
print(f"Max |Err|:  FW={np.max(fw_errors):.1f}%, Anth={np.max(anth_errors):.1f}%")
