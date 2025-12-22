#!/usr/bin/env python3
"""
分析為什麼 L6D6 < CK
L6D6 有額外 22 W/m² PAR，理論上應該 > CK
"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TREATMENT_CONFIGS, TARGETS, SIMULATION, get_env_for_treatment

def simulate_and_compare(treatment_name):
    """模擬並詳細分析"""
    p = UVAParams()
    env = get_env_for_treatment(treatment_name)

    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    X_d_init = dw_init_g / 1000.0 * ENV_BASE['plant_density']
    C_buf_init = X_d_init * 0.1
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6

    y0 = [X_d_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]

    t_start = SIMULATION['transplant_offset'] * 86400
    t_end = (SIMULATION['transplant_offset'] + SIMULATION['days']) * 86400

    sol = solve_ivp(
        fun=lambda t, y: uva_sun_derivatives(t, y, p, env),
        t_span=(t_start, t_end),
        y0=y0,
        method='RK45',
        max_step=60,
        dense_output=True
    )

    y_final = sol.y[:, -1]
    X_d, C_buf, LAI, Anth, Stress, E_stress = y_final

    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress, p)
    dw_per_plant = X_d / ENV_BASE['plant_density'] * 1000
    fw_per_plant = dw_per_plant / dw_fw_ratio

    return {
        'fw': fw_per_plant,
        'X_d': X_d,
        'dw': dw_per_plant,
        'dw_fw_ratio': dw_fw_ratio,
        'stress': Stress,
        'LAI': LAI,
    }

print('=' * 80)
print('L6D6 vs CK 詳細對比')
print('=' * 80)
print()

ck = simulate_and_compare('CK')
l6d6 = simulate_and_compare('L6D6')

print('最終結果:')
print(f'  CK:   {ck["fw"]:.2f}g')
print(f'  L6D6: {l6d6["fw"]:.2f}g')
print(f'  差異: {l6d6["fw"] - ck["fw"]:+.2f}g ({(l6d6["fw"] - ck["fw"]) / ck["fw"] * 100:+.1f}%)')
print()

print('=' * 80)
print('問題診斷')
print('=' * 80)
print()

print('1. 乾重對比 (X_d):')
print(f'   CK:   {ck["X_d"]:.6f} kg/m²')
print(f'   L6D6: {l6d6["X_d"]:.6f} kg/m²')
print(f'   差異: {(l6d6["X_d"] - ck["X_d"]) / ck["X_d"] * 100:+.1f}%')
print()

print('2. 單株乾重對比:')
print(f'   CK:   {ck["dw"]:.2f}g')
print(f'   L6D6: {l6d6["dw"]:.2f}g')
print(f'   差異: {(l6d6["dw"] - ck["dw"]) / ck["dw"] * 100:+.1f}%')
print()

print('3. LDMC 對比 (dw_fw_ratio):')
print(f'   CK:   {ck["dw_fw_ratio"]:.4f}')
print(f'   L6D6: {l6d6["dw_fw_ratio"]:.4f}')
print(f'   差異: {(l6d6["dw_fw_ratio"] - ck["dw_fw_ratio"]) / ck["dw_fw_ratio"] * 100:+.1f}%')
print()

print('4. Stress 對比:')
print(f'   CK:   {ck["stress"]:.6f}')
print(f'   L6D6: {l6d6["stress"]:.6f}')
print()

print('5. LAI 對比:')
print(f'   CK:   {ck["LAI"]:.4f}')
print(f'   L6D6: {l6d6["LAI"]:.4f}')
print()

print('=' * 80)
print('分析')
print('=' * 80)
print()

# 理論計算
I_ck = 57  # W/m²
I_l6d6 = 57 + 22  # 79 W/m²
par_increase = (I_l6d6 - I_ck) / I_ck * 100

print(f'理論上:')
print(f'  CK PAR: {I_ck} W/m²')
print(f'  L6D6 PAR: {I_l6d6} W/m² (增加 {par_increase:.1f}%)')
print(f'  預期 L6D6 乾重應該增加 ~{par_increase:.1f}%')
print()

actual_increase = (l6d6["X_d"] - ck["X_d"]) / ck["X_d"] * 100

print(f'實際:')
print(f'  L6D6 乾重增加: {actual_increase:+.1f}%')
print()

if l6d6["X_d"] > ck["X_d"]:
    print('✓ 乾重符合預期 (L6D6 > CK)')
    print()
    print('問題在 LDMC:')
    if l6d6["dw_fw_ratio"] > ck["dw_fw_ratio"]:
        print(f'  L6D6 的 LDMC 更高 ({l6d6["dw_fw_ratio"]:.4f} vs {ck["dw_fw_ratio"]:.4f})')
        print(f'  這是因為 Stress = {l6d6["stress"]:.4f} > 0')
        print(f'  Stress 導致 LDMC 降低機制被觸發')
        print()
        print('解決方案:')
        print('  → 降低 L6D6 的 Stress 累積')
        print('  → 或者調整 LDMC 對 Stress 的敏感度')
else:
    print('✗ 乾重也比 CK 低！')
    print('  問題: Stress 導致生長抑制過強')
    print()
    print('解決方案:')
    print('  → 大幅降低 stress_damage_coeff')
    print('  → 或降低 stress_photosynthesis_inhibition')

print()
print('=' * 80)
print('檢查 LDMC 動態公式')
print('=' * 80)
print()

p = UVAParams()

print('當前 LDMC 參數:')
print(f'  dw_fw_ratio_base = {p.dw_fw_ratio_base}')
print(f'  dw_fw_ratio_reduction_per_dose = {p.dw_fw_ratio_reduction_per_dose}')
print()

# 計算 LDMC
ldmc_ck = calculate_dynamic_dw_fw_ratio(ck["stress"], p)
ldmc_l6d6 = calculate_dynamic_dw_fw_ratio(l6d6["stress"], p)

print('LDMC 計算:')
print(f'  CK (Stress={ck["stress"]:.4f}): LDMC = {ldmc_ck:.4f}')
print(f'  L6D6 (Stress={l6d6["stress"]:.4f}): LDMC = {ldmc_l6d6:.4f}')
print()

# E_stress 影響
print('E_stress 累積劑量:')
print(f'  CK: E_stress = {ck.get("E_stress", "N/A")}')
print(f'  L6D6: E_stress = {l6d6.get("E_stress", "N/A")}')

print()
print('=' * 80)
