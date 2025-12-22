#!/usr/bin/env python3
"""
深入分析 H12D3 非線性累積機制
為什麼提高 stress_nonlinear_coeff 無法解決 H12D3？
"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TREATMENT_CONFIGS, TARGETS, ODE_SETTINGS, SIMULATION, get_env_for_treatment
from params_config import ALL_PARAMS

def simulate_with_tracking(treatment_name, params):
    """模擬並追蹤完整動態"""
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

    # 密集採樣
    t_eval = np.linspace(t_start, t_end, 500)

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
print('H12D3 非線性累積機制深入分析')
print('=' * 80)
print()

print('H12D3 配置:')
print('  UVA: 22 W/m², 12h/day, Day 32-35 (3 days)')
print('  目標: 60.6g (vs CK 87.0g, 需要 -30% 抑制)')
print()

# 測試不同的非線性參數
test_configs = [
    {'name': '當前 (v6.1.1)', 'stress_nonlinear_coeff': 4.0, 'K_nonlinear': 1.2},
    {'name': '極端非線性', 'stress_nonlinear_coeff': 10.0, 'K_nonlinear': 0.5},
    {'name': '超極端', 'stress_nonlinear_coeff': 20.0, 'K_nonlinear': 0.3},
]

for config in test_configs:
    print('=' * 80)
    print(f'測試: {config["name"]}')
    print(f'  stress_nonlinear_coeff = {config["stress_nonlinear_coeff"]}')
    print(f'  K_nonlinear = {config["K_nonlinear"]}')
    print('=' * 80)
    print()

    params = ALL_PARAMS.copy()
    params.update({
        'stress_nonlinear_coeff': config['stress_nonlinear_coeff'],
        'K_nonlinear': config['K_nonlinear'],
    })

    sol, p, env = simulate_with_tracking('H12D3', params)

    # 分析 UVA 照射期間的動態
    uva_start = 32 * 86400
    uva_end = 35 * 86400

    # 找出 UVA 期間的索引
    uva_indices = np.where((sol.t >= uva_start) & (sol.t <= uva_end))[0]

    if len(uva_indices) > 0:
        print('UVA 照射期間 Stress 動態:')
        print(f'{"Day":<8s} {"Hour":<8s} {"Stress":>10s} {"nonlinear_factor":>18s} {"damage_rate":>15s}')
        print('-' * 80)

        # 顯示關鍵時間點
        for i in [0, len(uva_indices)//4, len(uva_indices)//2, 3*len(uva_indices)//4, -1]:
            if i >= len(uva_indices):
                continue

            idx = uva_indices[i]
            t = sol.t[idx]
            y = sol.y[:, idx]
            Stress = y[4]

            # 計算非線性因子
            nonlinear_factor = 1.0 + p.stress_nonlinear_coeff * Stress / (p.K_nonlinear + Stress)

            # 估算損傷率 (簡化)
            damage_rate = p.stress_damage_coeff * 22.0 * nonlinear_factor  # 簡化版本

            day = t / 86400
            hour = (t / 3600) % 24

            print(f'{day:<8.2f} {hour:<8.1f} {Stress:10.4f} {nonlinear_factor:18.4f} {damage_rate:15.6e}')

    # 最終結果
    y_final = sol.y[:, -1]
    X_d, C_buf, LAI, Anth, Stress_final, E_stress = y_final

    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress_final, p)
    dw_per_plant = X_d / ENV_BASE['plant_density'] * 1000
    fw_per_plant = dw_per_plant / dw_fw_ratio

    target_fw = TARGETS['H12D3']['FW']
    fw_err = (fw_per_plant - target_fw) / target_fw * 100

    print()
    print('最終結果:')
    print(f'  最終 Stress: {Stress_final:.4f}')
    print(f'  預測鮮重: {fw_per_plant:.2f}g')
    print(f'  目標鮮重: {target_fw}g')
    print(f'  誤差: {fw_err:+.1f}%')
    print()

print('=' * 80)
print('關鍵發現')
print('=' * 80)
print()

print('''
為什麼非線性累積對 H12D3 效果有限？

1. **非線性公式的數學特性**:

   nonlinear_factor = 1 + stress_nonlinear_coeff * Stress / (K_nonlinear + Stress)

   問題: 這是一個 Michaelis-Menten 型飽和函數
   - 當 Stress << K_nonlinear 時，nonlinear_factor ≈ 1 (幾乎無放大)
   - 只有當 Stress >> K_nonlinear 時，才接近最大放大

2. **H12D3 的 Stress 累積過低**:

   - H12D3 只照射 3 天
   - Stress 在 Day 32-35 從 0 增長到 ~1.0
   - 即使 K_nonlinear = 0.5，Stress = 1.0 時:
     nonlinear_factor = 1 + 10.0 * 1.0 / (0.5 + 1.0) = 1 + 6.67 = 7.67

   - 但 Stress 大部分時間 < 0.5，放大效果更弱

3. **需要的放大倍數**:

   - 目標: 60.6g (vs CK 87.0g, -30% 抑制)
   - 當前預測: ~90g (+49% 誤差)
   - 需要將 Stress 提高 ~10x 才能達到目標
   - 但非線性放大最多只能達到 7-8x (在 Stress 很高時)

4. **惡性循環**:

   - 初期 Stress 低 → 非線性因子小 → 累積慢
   - 累積慢 → Stress 仍然低 → 無法觸發更強的非線性
   - 需要更長的時間讓 Stress 進入快速累積階段
   - 但 H12D3 只有 3 天！

5. **為什麼長期組 (VL3D12, L6D12) 可以？**:

   - VL3D12: 12 天照射，Stress 有時間累積到 > 2.0
   - L6D12: 12 天照射，Stress 有時間累積
   - 足夠的時間讓 Stress 進入非線性快速累積階段

結論:
------
非線性累積機制需要**時間**讓 Stress 累積到臨界點。
H12D3 只有 3 天，不足以讓 Stress 進入快速累積階段。

解決方案:
--------
1. 提高基礎損傷率 (但會影響 L6D6)
2. 添加即時損傷機制 (acute damage)
3. 提高抑制靈敏度 (低 Stress 也能產生強抑制)
''')

print()
print('=' * 80)
print('測試方案 3: 提高抑制靈敏度')
print('=' * 80)
print()

print('當前抑制公式:')
print('  inhibition = stress_photosynthesis_inhibition * Stress / (K_stress + Stress)')
print()
print('問題: 當 Stress = 1.0, K_stress = 2.5 時:')
print('  inhibition = 0.60 * 1.0 / (2.5 + 1.0) = 0.171 (17.1% 抑制)')
print('  需要 ~30% 抑制才能達到目標')
print()

print('測試: 降低 K_stress，提高 inhibition 係數')

for k_stress in [1.0, 1.5, 2.0]:
    for inhibition in [0.70, 0.75, 0.80]:
        # 計算 Stress = 1.0 時的抑制
        inh_at_1 = inhibition * 1.0 / (k_stress + 1.0)

        print(f'  K_stress={k_stress}, inhibition={inhibition}: Stress=1.0時抑制={inh_at_1*100:.1f}%')

print()
print('=' * 80)
