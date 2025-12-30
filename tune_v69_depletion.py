"""調校 v6.9 修復容量耗竭參數"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import sys
sys.path.insert(0, '/home/kasm-user/projects/uva-simulation')

from simulate_uva_model import ALL_PARAMS, UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TREATMENT_CONFIGS, TARGETS, SIMULATION

# 初始條件設置函數
def get_initial_state(p):
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    C_buf_init = Xd_init * 0.1
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6
    return [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0]

def simulate_treatment(treatment, params_override=None):
    """模擬單一處理組"""
    params = ALL_PARAMS.copy()
    if params_override:
        params.update(params_override)

    p = UVAParams(params)
    initial_state = get_initial_state(p)

    transplant_day = SIMULATION['transplant_offset']
    simulation_days = SIMULATION['days']
    t_start = transplant_day * 86400
    t_end = (transplant_day + simulation_days) * 86400

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
        FW_sim = X_total / dw_fw_ratio * 1000  # g
        FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
        Anth_sim = (Anth_f / FW_total_kg) * 1e6  # ppm
        return FW_sim, Anth_sim, Stress_f
    return None, None, None

def objective(x):
    """優化目標函數"""
    K_repair_capacity, repair_depletion_power = x

    if K_repair_capacity < 10 or K_repair_capacity > 2000:
        return 1000
    if repair_depletion_power < 1.0 or repair_depletion_power > 5.0:
        return 1000

    params_override = {
        'K_repair_capacity': K_repair_capacity,
        'repair_depletion_power': repair_depletion_power
    }

    total_error = 0
    for treatment in ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']:
        FW_sim, Anth_sim, Stress_f = simulate_treatment(treatment, params_override)
        if FW_sim is None:
            return 1000

        FW_obs = TARGETS[treatment]['FW']
        Anth_obs = TARGETS[treatment].get('Anth', 0)

        fw_err = abs((FW_sim - FW_obs) / FW_obs) * 100
        anth_err = abs((Anth_sim - Anth_obs) / Anth_obs) * 100 if Anth_obs > 0 else 0

        # 給 H12D3 更高權重，因為它是最難的
        weight = 2.0 if treatment == 'H12D3' else 1.0
        total_error += weight * (fw_err + anth_err * 0.5)

    return total_error

# 網格搜索
print("=" * 80)
print("v6.9 修復容量耗竭參數調校")
print("=" * 80)

best_error = float('inf')
best_params = None

print("\n階段 1: 網格搜索")
print("-" * 80)

for K in [50, 100, 150, 200, 250, 300, 400, 500]:
    for n in [1.5, 2.0, 2.5, 3.0, 3.5, 4.0]:
        err = objective([K, n])
        if err < best_error:
            best_error = err
            best_params = [K, n]
            print(f"K={K:>4.0f}, n={n:.1f} → Error={err:>6.1f} *NEW BEST*")

print(f"\n最佳參數: K_repair_capacity={best_params[0]:.1f}, repair_depletion_power={best_params[1]:.2f}")

# 局部優化
print("\n階段 2: 局部優化")
print("-" * 80)

result = minimize(
    objective,
    best_params,
    method='Nelder-Mead',
    options={'maxiter': 200, 'xatol': 1.0, 'fatol': 0.1}
)

final_K = result.x[0]
final_n = result.x[1]
print(f"優化後: K_repair_capacity={final_K:.1f}, repair_depletion_power={final_n:.2f}")

# 最終驗證
print("\n" + "=" * 80)
print("最終驗證結果")
print("=" * 80)
print(f"{'Treatment':<12} {'FW_obs':>8} {'FW_sim':>8} {'FW_Err':>8}   {'Anth_obs':>9} {'Anth_sim':>9} {'Anth_Err':>9}   {'Stress':>8}")
print("-" * 80)

params_override = {
    'K_repair_capacity': final_K,
    'repair_depletion_power': final_n
}

fw_errors = []
anth_errors = []

for treatment in ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']:
    FW_sim, Anth_sim, Stress_f = simulate_treatment(treatment, params_override)
    FW_obs = TARGETS[treatment]['FW']
    Anth_obs = TARGETS[treatment].get('Anth', 0)

    fw_err = (FW_sim - FW_obs) / FW_obs * 100
    anth_err = (Anth_sim - Anth_obs) / Anth_obs * 100 if Anth_obs > 0 else 0

    fw_errors.append(abs(fw_err))
    anth_errors.append(abs(anth_err))

    print(f"{treatment:<12} {FW_obs:>7.1f}g {FW_sim:>7.1f}g {fw_err:>+6.1f}%   "
          f"{Anth_obs:>8.1f} {Anth_sim:>8.1f} {anth_err:>+7.1f}%   {Stress_f:>8.2f}")

print("=" * 80)
print(f"\n精度達標:")
print(f"{'✅' if sum(1 for e in fw_errors if e < 5) == 6 else '❌'} FW: {sum(1 for e in fw_errors if e < 5)}/6 < 5%")
print(f"{'✅' if sum(1 for e in anth_errors if e < 10) == 6 else '❌'} Anth: {sum(1 for e in anth_errors if e < 10)}/6 < 10%")

print(f"\n建議更新參數:")
print(f"'K_repair_capacity': {final_K:.1f},")
print(f"'repair_depletion_power': {final_n:.2f},")
