"""
花青素參數優化 - 優先讓 CK 和 L6D6 達標

策略：
1. CK (Stress=0): Anth 主要來自 base_rate
2. L6D6 (Stress=0.3): Anth 來自 base + V*S/(K+S)

目標：先讓這兩組都 <5%，再調整其他組
"""

import numpy as np
from scipy.integrate import solve_ivp
from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment
from simulate_uva_model_v7 import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio, ALL_PARAMS


def simulate_treatment(treatment, params):
    """模擬單一處理組"""
    p = UVAParams(params)
    env = get_env_for_treatment(treatment)
    target = TARGETS[treatment]

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
        return {
            'FW_sim': FW_sim, 'Anth_sim': Anth_sim, 'Stress': Stress_f,
            'FW_err': (FW_sim - target['FW']) / target['FW'] * 100,
            'Anth_err': (Anth_sim - target['Anth']) / target['Anth'] * 100
        }
    return None


def test_params(base, v_max, k_stress):
    """測試一組參數"""
    params = ALL_PARAMS.copy()
    params['base_anth_rate_light'] = base
    params['base_anth_rate_dark'] = base / 2
    params['V_max_anth'] = v_max
    params['K_stress_anth'] = k_stress

    results = {}
    for t in ['CK', 'L6D6']:
        results[t] = simulate_treatment(t, params)
    return results


print("=" * 70)
print("花青素參數優化 - CK 和 L6D6 優先")
print("=" * 70)

# 目標值
print("\n目標值:")
print(f"  CK:   Anth = {TARGETS['CK']['Anth']:.1f} ppm")
print(f"  L6D6: Anth = {TARGETS['L6D6']['Anth']:.1f} ppm")

# 首先測試當前參數
print("\n當前 v7.0 參數:")
print(f"  base = {ALL_PARAMS['base_anth_rate_light']:.2e}")
print(f"  V_max = {ALL_PARAMS['V_max_anth']:.2e}")
print(f"  K_stress = {ALL_PARAMS['K_stress_anth']:.2f}")

results = test_params(
    ALL_PARAMS['base_anth_rate_light'],
    ALL_PARAMS['V_max_anth'],
    ALL_PARAMS['K_stress_anth']
)
print("\n結果:")
for t, r in results.items():
    if r:
        print(f"  {t}: Anth={r['Anth_sim']:.1f} ({r['Anth_err']:+.1f}%), Stress={r['Stress']:.2f}")

# 分析問題
print("\n" + "=" * 70)
print("問題分析:")
print("=" * 70)
print("""
CK (Stress=0):
  Anth 主要來自 base_rate × FW × time
  要讓 CK 準確，需要調整 base_rate

L6D6 (Stress≈0.3):
  Anth 來自 (base + V × S / (K + S)) × FW × time
  S/(K+S) ≈ 0.3/(1.38+0.3) ≈ 0.18
  誘導部分 ≈ V × 0.18

如果 CK 太高：降低 base → L6D6 也會降低
如果 L6D6 太低：提高 V 或降低 K → CK 不受影響

關鍵洞察：V 和 K 只影響有 Stress 的組，不影響 CK！
""")

# 策略：先找讓 CK 準確的 base，再調整 V 和 K 讓 L6D6 準確
print("=" * 70)
print("策略：先校準 CK (調整 base)，再校準 L6D6 (調整 V, K)")
print("=" * 70)

# 步驟 1: 找讓 CK 準確的 base
print("\n步驟 1: 搜索讓 CK 準確的 base...")
best_base = None
best_ck_err = 100

for base in np.linspace(1.5e-10, 2.5e-10, 21):
    r = test_params(base, 2.71e-11, 1.38)['CK']
    if r and abs(r['Anth_err']) < abs(best_ck_err):
        best_ck_err = r['Anth_err']
        best_base = base
        print(f"  base={base:.2e} → CK Anth_err={r['Anth_err']:+.1f}%")

print(f"\n最佳 base = {best_base:.2e} (CK 誤差 = {best_ck_err:+.1f}%)")

# 步驟 2: 固定 base，調整 V 和 K 讓 L6D6 準確
print("\n步驟 2: 固定 base，搜索 V 和 K 讓 L6D6 準確...")
print(f"  使用 base = {best_base:.2e}")

best_combo = None
best_l6d6_err = 100

for v in np.linspace(1.0e-11, 8.0e-11, 15):
    for k in np.linspace(0.1, 2.0, 15):
        results = test_params(best_base, v, k)
        if results['L6D6']:
            l6d6_err = results['L6D6']['Anth_err']
            ck_err = results['CK']['Anth_err']
            # 確保 CK 仍然準確
            if abs(ck_err) < 5 and abs(l6d6_err) < abs(best_l6d6_err):
                best_l6d6_err = l6d6_err
                best_combo = (v, k)
                print(f"  V={v:.2e}, K={k:.2f} → CK={ck_err:+.1f}%, L6D6={l6d6_err:+.1f}%")

if best_combo:
    print(f"\n最佳組合: V={best_combo[0]:.2e}, K={best_combo[1]:.2f}")
    print(f"  CK 誤差 < 5%, L6D6 誤差 = {best_l6d6_err:+.1f}%")

# 最終驗證
print("\n" + "=" * 70)
print("最終驗證 (所有處理組)")
print("=" * 70)

final_params = ALL_PARAMS.copy()
final_params['base_anth_rate_light'] = best_base
final_params['base_anth_rate_dark'] = best_base / 2
if best_combo:
    final_params['V_max_anth'] = best_combo[0]
    final_params['K_stress_anth'] = best_combo[1]

print(f"\n最終參數:")
print(f"  base_anth_rate_light = {final_params['base_anth_rate_light']:.2e}")
print(f"  base_anth_rate_dark = {final_params['base_anth_rate_dark']:.2e}")
print(f"  V_max_anth = {final_params['V_max_anth']:.2e}")
print(f"  K_stress_anth = {final_params['K_stress_anth']:.2f}")

print("\n所有處理組結果:")
anth_errs = []
for treatment in ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']:
    r = simulate_treatment(treatment, final_params)
    if r:
        target = TARGETS[treatment]['Anth']
        err = r['Anth_err']
        anth_errs.append(abs(err))
        s = "✓" if abs(err) < 5 else ("△" if abs(err) < 10 else "✗")
        print(f"  {treatment:<8}: 目標={target:.1f} 模擬={r['Anth_sim']:.1f} ({err:+.1f}%{s}) Stress={r['Stress']:.1f}")

print(f"\n達標: {sum(1 for e in anth_errs if e < 5)}/6 (<5%), {sum(1 for e in anth_errs if e < 10)}/6 (<10%)")
print(f"平均誤差: {np.mean(anth_errs):.1f}%")
print(f"最大誤差: {max(anth_errs):.1f}%")
