#!/usr/bin/env python3
"""
測試 v6.6 所有處理組 (E_elapsed 機制 + 花青素)
"""
import numpy as np
from scipy.integrate import solve_ivp
import sys
sys.path.insert(0, '/home/kasm-user/projects/uva-simulation')

from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TREATMENT_CONFIGS, TARGETS, SIMULATION

# 初始條件 (v6.7: 5個狀態變量，移除 E_stress)
fw_init_g = SIMULATION['initial_fw_g']
p = UVAParams()
dw_init_g = fw_init_g * p.dw_fw_ratio_base
Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
C_buf_init = Xd_init * 0.1
LAI_init = (dw_init_g / 0.01) * 0.04
fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
Anth_init = 5.0 * fw_total_init / 1e6
initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0]  # 5個狀態

# 模擬時間 (從移植日開始)
transplant_day = SIMULATION['transplant_offset']
simulation_days = SIMULATION['days']
t_start = transplant_day * 86400
t_end = (transplant_day + simulation_days) * 86400

print("=" * 100)
print("v6.6 所有處理組驗證 (E_elapsed 機制)")
print("=" * 100)
print(f"{'Treatment':<10} {'FW_sim':>8} {'FW_exp':>8} {'FW_Err':>8} {'FW_OK':>6}  {'Anth_sim':>9} {'Anth_exp':>9} {'Anth_Err':>9} {'Anth_OK':>8}")
print("-" * 100)

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
        Xd_f, Cbuf_f, LAI_f, Anth_f, Stress_f = sol.y[:, -1]  # v6.7: 5個狀態

        # 計算鮮重
        dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress_f, p)
        X_total = Xd_f / ENV_BASE['plant_density']
        FW = X_total / dw_fw_ratio * 1000

        # 計算花青素 [ppm = mg/kg FW]
        Anth_total_kg = Anth_f / ENV_BASE['plant_density']  # kg/plant
        Anth_ppm = (Anth_total_kg / (FW / 1000)) * 1e6  # mg/kg

        # 鮮重誤差
        target_fw = TARGETS[treatment]['FW']
        fw_err = (FW - target_fw) / target_fw * 100
        fw_ok = "✅" if abs(fw_err) < 5.0 else "❌"
        fw_errors.append(abs(fw_err))

        # 花青素誤差
        target_anth = TARGETS[treatment]['Anth']
        anth_err = (Anth_ppm - target_anth) / target_anth * 100
        anth_ok = "✅" if abs(anth_err) < 10.0 else "❌"
        anth_errors.append(abs(anth_err))

        print(f"{treatment:<10} {FW:>7.1f}g {target_fw:>7.1f}g {fw_err:>+6.1f}% {fw_ok:>6}  {Anth_ppm:>8.1f}p {target_anth:>8.1f}p {anth_err:>+8.1f}% {anth_ok:>8}")
    else:
        print(f"{treatment:<10} 模擬失敗")

print("=" * 100)
print(f"FW 統計: 平均 {np.mean(fw_errors):.1f}%, 中位數 {np.median(fw_errors):.1f}%, 最大 {np.max(fw_errors):.1f}%")
print(f"Anth 統計: 平均 {np.mean(anth_errors):.1f}%, 中位數 {np.median(anth_errors):.1f}%, 最大 {np.max(anth_errors):.1f}%")
print(f"FW <5%: {sum(1 for e in fw_errors if e < 5)}/6  |  Anth <10%: {sum(1 for e in anth_errors if e < 10)}/6")
print("=" * 100)
