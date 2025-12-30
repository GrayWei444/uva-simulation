"""
測試：讓 n_day 和 k_day 自由校準

問題：v7.0 的 n_day=7 是「預設形狀」
目標：讓優化器自己找出最佳 n_day，看是否收斂到 ~7

策略：
1. 固定其他參數不變
2. 讓 n_day 在 [1, 15] 範圍內搜索
3. k_day 會根據 n_day 自動調整（使 6h 時 factor ≈ 3.8）
4. 評估 FW 達標率
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment
from simulate_uva_model_v7 import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio, ALL_PARAMS


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
    """評估參數組合，返回總誤差"""
    total_error = 0
    results = {}
    for t in ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']:
        FW, Anth, Stress = simulate_treatment(t, params)
        if FW is None:
            return 1e10, None
        target = TARGETS[t]
        fw_err = abs((FW - target['FW']) / target['FW'] * 100)
        anth_err = abs((Anth - target['Anth']) / target['Anth'] * 100)
        total_error += fw_err + anth_err * 0.5  # FW 權重較高
        results[t] = {
            'FW_sim': FW, 'Anth_sim': Anth, 'Stress': Stress,
            'FW_err': (FW - target['FW']) / target['FW'] * 100,
            'Anth_err': (Anth - target['Anth']) / target['Anth'] * 100
        }
    return total_error, results


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
print("測試：讓 n_day 和 k_day 自由校準")
print("=" * 80)

# 當前 v7.0 結果
print("\n當前 v7.0 參數 (k_day=1e-5, n_day=7):")
error, results = evaluate_params(ALL_PARAMS)
if results:
    score = print_results(results, "v7.0 結果:")

# 策略 1：網格搜索 n_day
print("\n" + "=" * 80)
print("策略 1：網格搜索 n_day (k_day 自動調整)")
print("=" * 80)
print("讓 k_day 自動調整使 6h 時 day_factor ≈ 3.8")
print()

best_score = 0
best_n_day = 7
best_k_day = 1e-5
best_results = None

for n_day in np.linspace(1, 15, 29):  # 從 1 到 15，每 0.5 測試一次
    # 自動計算 k_day 使 6h 時 factor = 3.8
    # factor = 1 + k × 6^n = 3.8
    # k = 2.8 / 6^n
    k_day = 2.8 / (6 ** n_day)

    params = ALL_PARAMS.copy()
    params['n_day'] = n_day
    params['k_day'] = k_day

    error, results = evaluate_params(params)
    if results:
        fw_ok = sum(1 for t in results if abs(results[t]['FW_err']) < 5)
        anth_ok = sum(1 for t in results if abs(results[t]['Anth_err']) < 10)
        score = fw_ok + anth_ok

        # 計算 12h factor
        factor_12h = 1 + k_day * (12 ** n_day)
        ratio_12_6 = factor_12h / 3.8  # 相對於 6h

        marker = " *" if score >= best_score else ""
        print(f"n_day={n_day:>5.2f}: k_day={k_day:.2e}, 12h_factor={factor_12h:>8.1f}, "
              f"12h/6h={ratio_12_6:>6.1f}x, FW {fw_ok}/6, Anth {anth_ok}/6 = {score}/12{marker}")

        if score > best_score:
            best_score = score
            best_n_day = n_day
            best_k_day = k_day
            best_results = results

print(f"\n最佳 n_day = {best_n_day:.2f}")
print(f"對應 k_day = {best_k_day:.2e}")
print(f"12h factor = {1 + best_k_day * (12 ** best_n_day):.1f}")
if best_results:
    print_results(best_results, "最佳結果:")

# 策略 2：更細的搜索
print("\n" + "=" * 80)
print("策略 2：在最佳 n_day 附近細搜索")
print("=" * 80)

for n_day in np.linspace(best_n_day - 1, best_n_day + 1, 21):
    k_day = 2.8 / (6 ** n_day)

    params = ALL_PARAMS.copy()
    params['n_day'] = n_day
    params['k_day'] = k_day

    error, results = evaluate_params(params)
    if results:
        fw_ok = sum(1 for t in results if abs(results[t]['FW_err']) < 5)
        anth_ok = sum(1 for t in results if abs(results[t]['Anth_err']) < 10)
        score = fw_ok + anth_ok

        if score > best_score:
            factor_12h = 1 + k_day * (12 ** n_day)
            print(f"n_day={n_day:>5.2f}: k_day={k_day:.2e}, 12h_factor={factor_12h:>8.1f}, {score}/12 *")
            best_score = score
            best_n_day = n_day
            best_k_day = k_day
            best_results = results

print("\n" + "=" * 80)
print("結論")
print("=" * 80)
print(f"""
優化器找到的最佳參數:
  n_day = {best_n_day:.2f}
  k_day = {best_k_day:.2e}
  12h factor = {1 + best_k_day * (12 ** best_n_day):.1f}

與 v7.0 原始設定比較:
  v7.0: n_day = 7.0, k_day = 1e-5
  優化: n_day = {best_n_day:.2f}, k_day = {best_k_day:.2e}

結論:
  {'優化器收斂到 n_day ≈ 7，證明這不是人為預設，而是數據驅動的結果' if abs(best_n_day - 7) < 1 else '優化器找到不同的 n_day，需要進一步分析'}
""")

# 策略 3：完全自由搜索（不固定 6h factor）
print("\n" + "=" * 80)
print("策略 3：完全自由搜索 (k_day, n_day)")
print("=" * 80)

best_score_free = 0
best_params_free = None

for n_day in np.linspace(3, 12, 19):
    for k_day_exp in np.linspace(-8, -2, 13):
        k_day = 10 ** k_day_exp

        params = ALL_PARAMS.copy()
        params['n_day'] = n_day
        params['k_day'] = k_day

        error, results = evaluate_params(params)
        if results:
            fw_ok = sum(1 for t in results if abs(results[t]['FW_err']) < 5)
            anth_ok = sum(1 for t in results if abs(results[t]['Anth_err']) < 10)
            score = fw_ok + anth_ok

            if score > best_score_free:
                factor_6h = 1 + k_day * (6 ** n_day)
                factor_12h = 1 + k_day * (12 ** n_day)
                print(f"n_day={n_day:>5.2f}, k_day={k_day:.2e}: "
                      f"6h={factor_6h:.1f}, 12h={factor_12h:.1f}, {score}/12 *")
                best_score_free = score
                best_params_free = (n_day, k_day, results)

if best_params_free:
    n_day, k_day, results = best_params_free
    print(f"\n完全自由搜索最佳結果:")
    print(f"  n_day = {n_day:.2f}")
    print(f"  k_day = {k_day:.2e}")
    print(f"  6h factor = {1 + k_day * (6 ** n_day):.1f}")
    print(f"  12h factor = {1 + k_day * (12 ** n_day):.1f}")
    print_results(results, "最佳結果:")
