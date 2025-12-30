"""
v8.0 完整校準腳本

問題：day_factor 的曲線形狀對了，但整體模型結果不對
原因：v7.0 的其他參數是針對 v7.0 的 day_factor 校準的

策略：
1. 保持 v8.0 的 ROS 動力學形狀 (E_threshold + power law)
2. 調整 stress_damage_coeff 使整體損傷水平正確
3. 可能需要調整 E_threshold 使 L6D6 和 H12D3 都達標
"""

import numpy as np
from scipy.integrate import solve_ivp
from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment
from simulate_uva_model_v8_ros import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio, ALL_PARAMS


def simulate_treatment(treatment, params):
    """模擬單一處理組"""
    p = UVAParams(params)
    env = get_env_for_treatment(treatment)

    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    initial_state = [Xd_init, Xd_init * 0.1, (dw_init_g / 0.01) * 0.04,
                     5.0 * fw_init_g * ENV_BASE['plant_density'] / 1000 / 1e6, 0.0]
    t_start = SIMULATION['transplant_offset'] * 86400
    t_end = (SIMULATION['transplant_offset'] + SIMULATION['days']) * 86400

    sol = solve_ivp(uva_sun_derivatives, (t_start, t_end), initial_state,
                    args=(p, env), method='RK45', max_step=3600, t_eval=[t_end])

    if sol.success:
        Xd_f, _, _, Anth_f, Stress_f = sol.y[:, -1]
        dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress_f, p)
        FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
        FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
        Anth_sim = Anth_f / FW_total_kg * 1e6
        return FW_sim, Anth_sim, Stress_f
    return None, None, None


def evaluate_params(params):
    """評估參數組合"""
    results = {}
    for t in ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']:
        FW, Anth, Stress = simulate_treatment(t, params)
        if FW is None:
            return None
        target = TARGETS[t]
        results[t] = {
            'FW_sim': FW, 'Anth_sim': Anth, 'Stress': Stress,
            'FW_err': (FW - target['FW']) / target['FW'] * 100,
            'Anth_err': (Anth - target['Anth']) / target['Anth'] * 100
        }
    return results


def print_results(results, label=""):
    print(f"\n{label}")
    print("-" * 80)
    for t in ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']:
        r = results[t]
        fw_ok = "✓" if abs(r['FW_err']) < 5 else "✗"
        anth_ok = "✓" if abs(r['Anth_err']) < 10 else "✗"
        print(f"{t:<8} FW:{r['FW_sim']:>5.1f}g({r['FW_err']:>+5.1f}%{fw_ok}) "
              f"Anth:{r['Anth_sim']:>5.1f}({r['Anth_err']:>+5.1f}%{anth_ok}) S:{r['Stress']:>5.1f}")

    fw_ok_count = sum(1 for t in results if abs(results[t]['FW_err']) < 5)
    anth_ok_count = sum(1 for t in results if abs(results[t]['Anth_err']) < 10)
    print(f"\n達標: FW {fw_ok_count}/6, Anth {anth_ok_count}/6, Total {fw_ok_count + anth_ok_count}/12")
    return fw_ok_count + anth_ok_count


print("=" * 80)
print("v8.0 完整校準")
print("=" * 80)

# 當前結果
print("\n當前參數:")
print(f"  E_threshold = {ALL_PARAMS['E_threshold']:.0f}")
print(f"  n_ros = {ALL_PARAMS['n_ros']:.2f}")
print(f"  k_ros_damage = {ALL_PARAMS['k_ros_damage']:.2e}")
print(f"  stress_damage_coeff = {ALL_PARAMS['stress_damage_coeff']:.2e}")

results = evaluate_params(ALL_PARAMS)
if results:
    print_results(results, "當前結果:")

# 分析：
# - L6D6 損傷太大 → 需要降低 stress_damage_coeff 或提高 E_threshold
# - H12D3 損傷也太大
# - VL3D12 其實還好

# 策略 1: 大幅降低 stress_damage_coeff，讓 day_factor 主導
print("\n" + "=" * 80)
print("策略 1: 降低 stress_damage_coeff")
print("=" * 80)

best_score = 0
best_params = None

for sdc in [1e-8, 5e-8, 1e-7, 2e-7, 5e-7, 1e-6, 2e-6]:
    params = ALL_PARAMS.copy()
    params['stress_damage_coeff'] = sdc
    results = evaluate_params(params)
    if results:
        fw_ok = sum(1 for t in results if abs(results[t]['FW_err']) < 5)
        anth_ok = sum(1 for t in results if abs(results[t]['Anth_err']) < 10)
        score = fw_ok + anth_ok
        if score >= best_score:
            print(f"sdc={sdc:.0e}: FW {fw_ok}/6, Anth {anth_ok}/6, Total {score}/12")
            if score > best_score:
                best_score = score
                best_params = params.copy()

if best_params:
    print(f"\n最佳 stress_damage_coeff = {best_params['stress_damage_coeff']:.2e}")
    results = evaluate_params(best_params)
    if results:
        print_results(results, "最佳結果:")

# 策略 2: 調整 E_threshold
print("\n" + "=" * 80)
print("策略 2: 調整 E_threshold + stress_damage_coeff")
print("=" * 80)

best_score = 0
best_params = None

for E_th in [150000, 180000, 200000, 220000]:
    for sdc in [5e-8, 1e-7, 2e-7, 5e-7, 1e-6]:
        params = ALL_PARAMS.copy()
        params['E_threshold'] = E_th
        params['stress_damage_coeff'] = sdc
        # 需要重新計算 k_ros_damage
        # 使 12h factor 大約等於 v7.0 的值
        I_UVA = 11.0
        E_12h = I_UVA * 12 * 3600
        E_excess = max(0, E_12h - E_th)
        target_factor = 359.0  # v7.0 的 12h factor
        if E_excess > 0:
            k = (target_factor - 1) / (E_excess ** params['n_ros'])
            params['k_ros_damage'] = k

        results = evaluate_params(params)
        if results:
            fw_ok = sum(1 for t in results if abs(results[t]['FW_err']) < 5)
            anth_ok = sum(1 for t in results if abs(results[t]['Anth_err']) < 10)
            score = fw_ok + anth_ok
            if score >= best_score:
                if score > best_score:
                    print(f"E_th={E_th:.0f}, sdc={sdc:.0e}, k={params['k_ros_damage']:.2e}: "
                          f"FW {fw_ok}/6, Anth {anth_ok}/6, Total {score}/12 *")
                    best_score = score
                    best_params = params.copy()

if best_params:
    print(f"\n最佳參數:")
    print(f"  E_threshold = {best_params['E_threshold']:.0f}")
    print(f"  stress_damage_coeff = {best_params['stress_damage_coeff']:.2e}")
    print(f"  k_ros_damage = {best_params['k_ros_damage']:.2e}")
    results = evaluate_params(best_params)
    if results:
        print_results(results, "最佳結果:")

# 策略 3: 更細的網格搜索
print("\n" + "=" * 80)
print("策略 3: 細網格搜索")
print("=" * 80)

if best_params:
    E_th_center = best_params['E_threshold']
    sdc_center = best_params['stress_damage_coeff']

    for E_th in np.linspace(E_th_center * 0.8, E_th_center * 1.2, 10):
        for sdc in np.linspace(sdc_center * 0.5, sdc_center * 2.0, 10):
            params = ALL_PARAMS.copy()
            params['E_threshold'] = E_th
            params['stress_damage_coeff'] = sdc
            # 重新計算 k_ros_damage
            I_UVA = 11.0
            E_12h = I_UVA * 12 * 3600
            E_excess = max(0, E_12h - E_th)
            if E_excess > 0:
                k = (359.0 - 1) / (E_excess ** params['n_ros'])
                params['k_ros_damage'] = k

            results = evaluate_params(params)
            if results:
                fw_ok = sum(1 for t in results if abs(results[t]['FW_err']) < 5)
                anth_ok = sum(1 for t in results if abs(results[t]['Anth_err']) < 10)
                score = fw_ok + anth_ok
                if score > best_score:
                    print(f"E_th={E_th:.0f}, sdc={sdc:.2e}: FW {fw_ok}/6, Anth {anth_ok}/6, Total {score}/12 *")
                    best_score = score
                    best_params = params.copy()

if best_params:
    print("\n" + "=" * 80)
    print("最終最佳參數")
    print("=" * 80)
    print(f"""
    'E_threshold': {best_params['E_threshold']:.0f},
    'n_ros': {best_params['n_ros']:.2f},
    'k_ros_damage': {best_params['k_ros_damage']:.2e},
    'stress_damage_coeff': {best_params['stress_damage_coeff']:.2e},
    """)
    results = evaluate_params(best_params)
    if results:
        print_results(results, "最終結果:")
