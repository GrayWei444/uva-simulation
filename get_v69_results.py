"""獲取 v6.9 完整結果 (FW, Anth, Stress)"""
import numpy as np
from scipy.integrate import solve_ivp
import sys
sys.path.insert(0, '/home/kasm-user/projects/uva-simulation')

from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TREATMENT_CONFIGS, TARGETS, SIMULATION

# 初始條件 (v6.9: 5個狀態變量)
fw_init_g = SIMULATION['initial_fw_g']
p = UVAParams()
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

print("=" * 120)
print("v6.9 完整驗證結果 (修復容量耗竭機制 - 無預設閾值)")
print("=" * 120)
print(f"{'Treatment':<12} {'FW_obs':>8} {'FW_sim':>8} {'FW_Err':>8}   {'Anth_obs':>9} {'Anth_sim':>9} {'Anth_Err':>9}   {'Stress':>8}")
print("-" * 120)

fw_errors = []
anth_errors = []

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
        FW_sim = X_total / dw_fw_ratio * 1000  # g

        # 計算花青素濃度 (ppm)
        FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
        Anth_conc_sim = (Anth_f / FW_total_kg) * 1e6  # ppm

        # 獲取實驗觀測值
        FW_obs = TARGETS[treatment]['FW']
        Anth_obs = TARGETS[treatment].get('Anth', 0)

        # 計算誤差
        fw_err = (FW_sim - FW_obs) / FW_obs * 100
        anth_err = (Anth_conc_sim - Anth_obs) / Anth_obs * 100 if Anth_obs > 0 else 0

        fw_errors.append(abs(fw_err))
        anth_errors.append(abs(anth_err))

        print(f"{treatment:<12} {FW_obs:>7.1f}g {FW_sim:>7.1f}g {fw_err:>+6.1f}%   "
              f"{Anth_obs:>8.1f} {Anth_conc_sim:>8.1f} {anth_err:>+7.1f}%   {Stress_f:>8.2f}")
    else:
        print(f"{treatment:<12} 模擬失敗")

print("=" * 120)
print(f"\n統計摘要:")
print(f"FW  - 平均絕對誤差: {np.mean(fw_errors):.1f}%, 最大絕對誤差: {np.max(fw_errors):.1f}%")
print(f"Anth - 平均絕對誤差: {np.mean(anth_errors):.1f}%, 最大絕對誤差: {np.max(anth_errors):.1f}%")
print(f"\n精度達標:")
print(f"{'✅' if sum(1 for e in fw_errors if e < 5) == 6 else '❌'} FW: {sum(1 for e in fw_errors if e < 5)}/6 < 5%")
print(f"{'✅' if sum(1 for e in anth_errors if e < 10) == 6 else '❌'} Anth: {sum(1 for e in anth_errors if e < 10)}/6 < 10%")
