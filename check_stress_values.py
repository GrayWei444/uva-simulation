#!/usr/bin/env python3
"""
檢查 L6D6 模擬中的 Stress 實際值
找出為什麼 Stress 不是 0
"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives
from model_config import ENV_BASE, TREATMENT_CONFIGS, SIMULATION

def simulate_with_logging(treatment_name):
    """模擬並記錄關鍵時間點的 Stress 值"""

    p = UVAParams()
    env = ENV_BASE.copy()
    env.update(TREATMENT_CONFIGS[treatment_name])

    plant_density = env['plant_density']

    # 初始條件
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * 0.05
    X_d_init = dw_init_g / 1000.0 * plant_density
    C_buf_init = X_d_init * 0.1
    LAI_init = 0.5

    y0 = [X_d_init, C_buf_init, LAI_init, 0.0, 0.0, 0.0]  # [X_d, C_buf, LAI, Anth, Stress, E_stress]

    # 模擬時間
    t_start = 14 * 86400
    t_end = 35 * 86400

    # 記錄時間點和狀態
    logged_times = []
    logged_states = []

    def event_logger(t, y):
        """記錄狀態"""
        logged_times.append(t)
        logged_states.append(y.copy())
        return 1  # 永遠不觸發事件

    # ODE 求解
    sol = solve_ivp(
        fun=lambda t, y: uva_sun_derivatives(t, y, p, env),
        t_span=(t_start, t_end),
        y0=y0,
        method='RK45',
        max_step=60,
        dense_output=True
    )

    return sol

print('=' * 80)
print('檢查 L6D6 模擬中的 Stress 值')
print('=' * 80)
print()

# 模擬 L6D6
print('模擬 L6D6...')
sol = simulate_with_logging('L6D6')
print('完成')
print()

# 關鍵時間點
key_times = {
    'Day 14 (開始)': 14 * 86400,
    'Day 28 (UVA 前)': 28 * 86400,
    'Day 29 12:00 (UVA 第1天中午)': (29 * 86400) + (12 * 3600),
    'Day 30 12:00 (UVA 第2天中午)': (30 * 86400) + (12 * 3600),
    'Day 31 12:00 (UVA 第3天中午)': (31 * 86400) + (12 * 3600),
    'Day 35 (結束)': 35 * 86400,
}

print('關鍵時間點的狀態:')
print('-' * 80)
print(f'{"時間點":<30s} {"X_d":>12s} {"C_buf":>12s} {"LAI":>8s} {"Anth":>12s} {"Stress":>12s}')
print('-' * 80)

for time_name, t in key_times.items():
    y = sol.sol(t)
    X_d, C_buf, LAI, Anth, Stress, E_stress = y
    print(f'{time_name:<30s} {X_d:12.6f} {C_buf:12.8f} {LAI:8.4f} {Anth:12.8e} {Stress:12.8f}')

print()

# 檢查 Day 30 12:00 的導數
print('=' * 80)
print('Day 30 12:00 的導數詳情')
print('=' * 80)
print()

t_test = (30 * 86400) + (12 * 3600)
y_test = sol.sol(t_test)

p = UVAParams()
env = ENV_BASE.copy()
env.update(TREATMENT_CONFIGS['L6D6'])

derivs = uva_sun_derivatives(t_test, y_test, p, env)

X_d, C_buf, LAI, Anth, Stress, E_stress = y_test
dXd_dt, dCbuf_dt, dLAI_dt, dAnth_dt, dStress_dt, dE_stress_dt = derivs

print(f'狀態:')
print(f'  X_d     = {X_d:.8f} kg/m²')
print(f'  C_buf   = {C_buf:.10f} kg C/m²')
print(f'  LAI     = {LAI:.6f}')
print(f'  Anth    = {Anth:.8e} kg/m²')
print(f'  Stress  = {Stress:.8f}  ← 關鍵！')
print(f'  E_stress= {E_stress:.6f}')
print()

print(f'導數:')
print(f'  dX_d/dt     = {dXd_dt:.8e} kg/m²/s')
print(f'  dC_buf/dt   = {dCbuf_dt:.8e} kg C/m²/s  ← 問題！')
print(f'  dLAI/dt     = {dLAI_dt:.8e} /s')
print(f'  dAnth/dt    = {dAnth_dt:.8e} kg/m²/s')
print(f'  dStress/dt  = {dStress_dt:.8e} /s')
print()

# 分析 Stress 的影響
if Stress > 0.1:
    print(f'⚠️  Stress = {Stress:.2f} 顯著大於 0！')
    print(f'   這會導致:')

    # 修復碳消耗
    C_buf_positive = max(C_buf, 0)
    repair_capacity = p.base_repair_capacity + \
                      p.carbon_repair_bonus * C_buf_positive / (p.K_carbon + C_buf_positive + 1e-9)
    repair_rate = p.stress_repair_coeff * Stress * repair_capacity
    repair_carbon_consumption = repair_rate * p.repair_carbon_cost

    print(f'   - 修復率: {repair_rate:.8e} /s')
    print(f'   - 修復碳消耗: {repair_carbon_consumption:.8e} kg C/m²/s')

    # Stress 抑制
    stress_inhibition = Stress / (p.K_stress + Stress + 1e-9)
    xd_reduction = p.stress_photosynthesis_inhibition * stress_inhibition

    print(f'   - Stress 抑制因子: {stress_inhibition:.4f}')
    print(f'   - 光合作用抑制: {xd_reduction:.2%}')
    print()

    print(f'根本問題: 為什麼 Stress 會累積到這麼高？')
    print(f'   UVA 只有 22 W/m²，按理說損傷應該很小')

else:
    print(f'✓ Stress = {Stress:.4f} 很小，不是問題')

print()
print('=' * 80)
