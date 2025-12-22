#!/usr/bin/env python3
"""
Debug: 為什麼 H12D3 Stress 累積這麼低？
峰值只有 0.4，遠低於預期的數百
"""

import numpy as np
from scipy.integrate import solve_ivp
from simulate_uva_model import UVAParams, uva_sun_derivatives
from model_config import ENV_BASE, TREATMENT_CONFIGS, SIMULATION, get_env_for_treatment

def trace_stress_dynamics(treatment_name, day_number):
    """追蹤單一天內的 Stress 動態，逐小時分析"""
    p = UVAParams()
    env = get_env_for_treatment(treatment_name)

    # 從前一天晚上的狀態開始
    # 簡化：使用典型的 Day 32 狀態
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    X_d_init = dw_init_g / 1000.0 * ENV_BASE['plant_density'] * 3.0  # Day 32 已經長大
    C_buf_init = X_d_init * 0.15
    LAI_init = 4.5  # Day 32 時 LAI 較高
    Anth_init = 1e-6
    Stress_init = 0.2  # 前一天遺留的 Stress

    y0 = [X_d_init, C_buf_init, LAI_init, Anth_init, Stress_init, 0.0]

    # 模擬一天
    t_start = day_number * 86400
    t_end = (day_number + 1) * 86400

    # 每 15 分鐘一個點
    t_eval = np.linspace(t_start, t_end, 96)

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
print('H12D3 Stress 累積 Debug')
print('=' * 80)
print()

sol, p, env = trace_stress_dynamics('H12D3', 32)

# 找幾個關鍵時刻，手動計算 dStress/dt
print('關鍵時刻的 Stress 動態分析:')
print()

# 選擇 12:00 (照射中段)
target_hour = 12
target_t = (32 * 86400) + (target_hour * 3600)
idx = np.argmin(np.abs(sol.t - target_t))

t = sol.t[idx]
y = sol.y[:, idx]
X_d, C_buf, LAI, Anth, Stress, E_stress = y

print(f'時刻: Day 32, {target_hour}:00')
print(f'狀態:')
print(f'  X_d = {X_d:.6f}')
print(f'  C_buf = {C_buf:.6f}')
print(f'  LAI = {LAI:.4f}')
print(f'  Stress = {Stress:.6f}')
print()

# 手動計算導數
derivs = uva_sun_derivatives(t, y, p, env)
dStress_dt = derivs[4]

print(f'dStress/dt = {dStress_dt:.8e}')
print()

# 分解各個因子
print('=' * 80)
print('損傷率分解 (手動計算)')
print('=' * 80)
print()

# 基礎參數
I_UVA = 22.0  # W/m²
damage_coeff = p.stress_damage_coeff

print(f'1. 基礎損傷率:')
print(f'   stress_damage_coeff = {damage_coeff:.2e}')
print(f'   I_UVA = {I_UVA} W/m²')
print(f'   基礎損傷 = {damage_coeff} × {I_UVA} = {damage_coeff * I_UVA:.2e}')
print()

# LAI 脆弱性
LAI_ref_vuln = p.LAI_ref_vuln
n_vuln = p.n_vuln
cap_vuln = p.cap_vuln

base_vuln = (LAI_ref_vuln / LAI) ** n_vuln
vulnerability = cap_vuln * base_vuln / (cap_vuln + base_vuln)

print(f'2. LAI 脆弱性:')
print(f'   LAI = {LAI:.4f}')
print(f'   LAI_ref_vuln = {LAI_ref_vuln}')
print(f'   n_vuln = {n_vuln}')
print(f'   base_vuln = ({LAI_ref_vuln}/{LAI})^{n_vuln} = {base_vuln:.4f}')
print(f'   vulnerability = {cap_vuln} × {base_vuln} / ({cap_vuln} + {base_vuln}) = {vulnerability:.4f}')
print()

# 日內非線性
hour = (t / 3600) % 24
I_UVA_config = env.get('uva_intensity', 22.0)
uva_hour_on = env.get('uva_hour_on', 6)

hours_elapsed = hour - uva_hour_on  # 12 - 6 = 6h
E_elapsed = I_UVA_config * hours_elapsed * 3.6  # 22 × 6 × 3.6 = 475.2 kJ/m²

normalized_E = (E_elapsed - p.E_50) / p.E_scale

def softplus(x, sharpness):
    if x > 10:
        return x
    return np.log(1.0 + np.exp(sharpness * x)) / sharpness

excess_normalized = softplus(normalized_E, p.sharpness_intraday)
intraday_factor = 1.0 + p.k_intraday * (excess_normalized ** p.m_intraday)

print(f'3. 日內能量非線性:')
print(f'   hours_elapsed = {hours_elapsed} h')
print(f'   E_elapsed = {I_UVA_config} × {hours_elapsed} × 3.6 = {E_elapsed:.1f} kJ/m²')
print(f'   normalized_E = ({E_elapsed} - {p.E_50}) / {p.E_scale} = {normalized_E:.3f}')
print(f'   intraday_factor = {intraday_factor:.4f}')
print()

# Stress 非線性累積
nonlinear_factor = 1.0 + p.stress_nonlinear_coeff * Stress / (p.K_nonlinear + Stress)

print(f'4. Stress 非線性累積:')
print(f'   Stress = {Stress:.6f}')
print(f'   stress_nonlinear_coeff = {p.stress_nonlinear_coeff}')
print(f'   K_nonlinear = {p.K_nonlinear}')
print(f'   nonlinear_factor = 1 + {p.stress_nonlinear_coeff} × {Stress:.6f} / ({p.K_nonlinear} + {Stress:.6f})')
print(f'                    = {nonlinear_factor:.4f}')
print()

# 總損傷率
damage_rate = damage_coeff * I_UVA * vulnerability * intraday_factor * nonlinear_factor

print(f'5. 總損傷率:')
print(f'   damage_rate = {damage_coeff:.2e} × {I_UVA} × {vulnerability:.4f} × {intraday_factor:.4f} × {nonlinear_factor:.4f}')
print(f'               = {damage_rate:.8e}')
print()

# 修復率
C_buf_positive = max(C_buf, 0)
repair_capacity = p.base_repair_capacity + p.carbon_repair_bonus * C_buf_positive / (p.K_carbon + C_buf_positive)
repair_rate = p.stress_repair_coeff * Stress * repair_capacity

print(f'6. 修復率:')
print(f'   stress_repair_coeff = {p.stress_repair_coeff:.2e}')
print(f'   repair_capacity = {repair_capacity:.4f}')
print(f'   repair_rate = {p.stress_repair_coeff:.2e} × {Stress:.6f} × {repair_capacity:.4f}')
print(f'               = {repair_rate:.8e}')
print()

# 淨變化率
net_rate = damage_rate - repair_rate

print(f'7. 淨變化率:')
print(f'   dStress/dt = damage_rate - repair_rate')
print(f'              = {damage_rate:.8e} - {repair_rate:.8e}')
print(f'              = {net_rate:.8e}')
print()

print(f'與模型計算對比: {dStress_dt:.8e}')
print()

print('=' * 80)
print('問題診斷')
print('=' * 80)
print()

print(f'''
關鍵發現:

1. **LAI 脆弱性極低**: {vulnerability:.4f}
   - LAI = {LAI:.2f} 已經很高（成熟植株）
   - LAI_ref_vuln = {LAI_ref_vuln}
   - 脆弱性被 LAI 大幅降低
   - 這是主要問題！

2. **日內非線性未觸發** (intraday_factor = {intraday_factor:.4f})
   - 照射 6h 時累積能量 = {E_elapsed:.1f} kJ/m²
   - E_50 = {p.E_50} kJ/m²
   - 還沒達到觸發點

3. **Stress 非線性也低** (nonlinear_factor = {nonlinear_factor:.4f})
   - Stress = {Stress:.6f} 太小
   - 無法觸發正反饋

4. **修復vs損傷平衡**:
   - 損傷率 = {damage_rate:.2e}
   - 修復率 = {repair_rate:.2e}
   - 淨累積 = {net_rate:.2e}
   - 修復率接近損傷率！

解決方案:
--------
**問題根源**: LAI 太高 → 脆弱性極低 → 損傷率低

方案 A: 降低 LAI_ref_vuln (讓高 LAI 也脆弱)
   LAI_ref_vuln: {LAI_ref_vuln} → 3.0 或更低

方案 B: 降低 n_vuln (減弱 LAI 對脆弱性的影響)
   n_vuln: {n_vuln} → 4 或更低

方案 C: 提高基礎損傷率（但會影響所有組）
   stress_damage_coeff: {damage_coeff:.2e} → 更高
''')

print()
print('=' * 80)
print('測試方案: 降低 LAI_ref_vuln')
print('=' * 80)
print()

for lai_ref_test in [6.5, 5.0, 4.0, 3.0, 2.5]:
    # 重新計算脆弱性
    base_vuln_test = (lai_ref_test / LAI) ** n_vuln
    vulnerability_test = cap_vuln * base_vuln_test / (cap_vuln + base_vuln_test)

    # 重新計算損傷率
    damage_rate_test = damage_coeff * I_UVA * vulnerability_test * intraday_factor * nonlinear_factor

    # 預估 Stress 峰值 (簡化計算)
    # 假設照射 12h，修復忽略
    stress_estimate = damage_rate_test * 12 * 3600

    print(f'  LAI_ref={lai_ref_test}: vulnerability={vulnerability_test:6.2f}, '
          f'damage_rate={damage_rate_test:.2e}, 估計峰值Stress≈{stress_estimate:.2f}')

print()
print('=' * 80)
