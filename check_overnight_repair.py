#!/usr/bin/env python3
"""
檢查 H12D3 的夜間修復
關鍵問題：Stress 峰值雖高 (213)，但夜間修復太強，隔天又降回 0？
"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TREATMENT_CONFIGS, TARGETS, SIMULATION, get_env_for_treatment
from params_config import ALL_PARAMS

def simulate_full_period(treatment_name, params):
    """完整模擬獲取詳細動態"""
    p = UVAParams(params)
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

    # 密集採樣，特別是 Day 32-35
    t_eval_list = []

    # Day 32-35: 每小時一個點
    for day in range(32, 36):
        for hour in range(24):
            t = (day * 86400) + (hour * 3600)
            if t >= t_start and t <= t_end:
                t_eval_list.append(t)

    t_eval = np.array(t_eval_list) if t_eval_list else None

    sol = solve_ivp(
        fun=lambda t, y: uva_sun_derivatives(t, y, p, env),
        t_span=(t_start, t_end),
        y0=y0,
        method='RK45',
        max_step=60,
        dense_output=True,
        t_eval=t_eval
    )

    return sol, p, env

print('=' * 80)
print('H12D3 夜間修復分析')
print('=' * 80)
print()

# 使用當前參數
params = ALL_PARAMS.copy()
sol, p, env = simulate_full_period('H12D3', params)

print('當前參數:')
print(f'  stress_repair_coeff = {p.stress_repair_coeff}')
print(f'  base_repair_capacity = {p.base_repair_capacity}')
print(f'  carbon_repair_bonus = {p.carbon_repair_bonus}')
print()

# 分析 Day 32-35 的詳細動態
print('=' * 80)
print('Day 32-35 每小時 Stress 動態')
print('=' * 80)
print()

print(f'{"Day":<6s} {"Hour":<6s} {"Stress":>10s} {"說明":<30s}')
print('-' * 60)

for i, t in enumerate(sol.t):
    day = t / 86400
    hour = (t / 3600) % 24

    if day < 32 or day >= 36:
        continue

    stress = sol.y[4, i]

    # 標註關鍵時刻
    note = ''
    if hour == 6 and abs(hour - round(hour)) < 0.1:
        note = 'UVA 開始'
    elif hour == 18 and abs(hour - round(hour)) < 0.1:
        note = 'UVA 結束，開始修復'
    elif hour == 0 and abs(hour - round(hour)) < 0.1:
        note = '隔天開始'

    # 只顯示關鍵時刻
    if note or (stress > 10 and hour in [9, 12, 15, 21]):
        print(f'{day:<6.2f} {hour:<6.0f} {stress:10.2f} {note:<30s}')

print()

# 分析每天的峰值和修復
print('=' * 80)
print('每日 Stress 峰值與夜間修復')
print('=' * 80)
print()

print(f'{"Day":<6s} {"峰值":>10s} {"峰值時刻":>10s} {"隔天6am":>10s} {"修復率":>10s}')
print('-' * 60)

for day in range(32, 35):
    # 找當天的 Stress 資料
    day_mask = (sol.t >= day * 86400) & (sol.t < (day + 1) * 86400)
    day_stress = sol.y[4, day_mask]
    day_times = sol.t[day_mask]

    if len(day_stress) == 0:
        continue

    # 峰值
    peak_stress = np.max(day_stress)
    peak_idx = np.argmax(day_stress)
    peak_time = day_times[peak_idx]
    peak_hour = (peak_time / 3600) % 24

    # 隔天早上 6am 的 Stress
    next_6am = (day + 1) * 86400 + 6 * 3600
    # 找最接近的時間點
    closest_idx = np.argmin(np.abs(sol.t - next_6am))
    next_morning_stress = sol.y[4, closest_idx]

    # 修復率
    repair_rate = (peak_stress - next_morning_stress) / peak_stress * 100 if peak_stress > 0 else 0

    print(f'{day:<6d} {peak_stress:10.2f} {peak_hour:10.1f}h {next_morning_stress:10.2f} {repair_rate:9.1f}%')

print()

# 最終結果
y_final = sol.y[:, -1]
X_d, C_buf, LAI, Anth, Stress_final, E_stress = y_final

dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress_final, p)
dw_per_plant = X_d / ENV_BASE['plant_density'] * 1000
fw_per_plant = dw_per_plant / dw_fw_ratio

target_fw = TARGETS['H12D3']['FW']
fw_err = (fw_per_plant - target_fw) / target_fw * 100

print('=' * 80)
print('最終結果')
print('=' * 80)
print(f'  最終 Stress: {Stress_final:.2f}')
print(f'  預測鮮重: {fw_per_plant:.2f}g')
print(f'  目標鮮重: {target_fw}g')
print(f'  誤差: {fw_err:+.1f}%')
print()

print('=' * 80)
print('診斷')
print('=' * 80)
print()

print('''
關鍵發現：

1. **Stress 峰值很高**（~200+）
   - 證明當日內累積機制正常工作
   - 非線性正反饋有效

2. **夜間修復效果如何？**
   - 查看上表的"修復率"欄
   - 如果修復率很高（>80%），說明夜間修復太強
   - Stress 無法跨夜累積

3. **可能的解決方案**：
   a) 降低修復率 (stress_repair_coeff)
   b) 降低修復容量 (base_repair_capacity)
   c) 或者：保留部分 Stress 不被修復（設定修復上限）

4. **為什麼最終 Stress 很低？**
   - 如果每天都修復掉大部分 Stress
   - 到 Day 35 時累積的總 Stress 就會很低
   - 導致抑制不足
''')

print()
print('=' * 80)
print('測試：降低修復率')
print('=' * 80)
print()

# 測試降低修復率
for repair_coeff in [1.0e-5, 8.0e-6, 6.0e-6, 4.0e-6]:
    test_params = ALL_PARAMS.copy()
    test_params['stress_repair_coeff'] = repair_coeff

    sol_test, p_test, env_test = simulate_full_period('H12D3', test_params)

    # 最終 Stress
    stress_final = sol_test.y[4, -1]

    # 預測鮮重
    X_d = sol_test.y[0, -1]
    dw_fw = calculate_dynamic_dw_fw_ratio(stress_final, p_test)
    dw_per_plant = X_d / ENV_BASE['plant_density'] * 1000
    fw = dw_per_plant / dw_fw

    err = (fw - 60.6) / 60.6 * 100

    print(f'  repair_coeff={repair_coeff:.1e}: 最終Stress={stress_final:6.2f}, FW={fw:6.2f}g, 誤差={err:+6.1f}%')

print()
print('=' * 80)
