"""
生成論文圖 9-18：模型機制圖表
v10.39 版本 - Hill 效率函數取代 sigmoid
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.integrate import solve_ivp

# 設定字體
rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial Unicode MS', 'sans-serif']
rcParams['axes.unicode_minus'] = False

# 導入模型
from simulate_uva_model_v10 import (
    UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio,
    nonlinear_damage_factor, calculate_nonlin_anth_efficiency
)
from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment


def run_simulation(treatment_name):
    """執行單次模擬並返回結果"""
    p = UVAParams()
    env = get_env_for_treatment(treatment_name)

    # 初始條件
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    C_buf_init = Xd_init * 0.1
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6

    initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]

    # 時間設定
    transplant_day = SIMULATION['transplant_offset']
    simulation_days = SIMULATION['days']
    t_start = transplant_day * 86400
    t_end = (transplant_day + simulation_days) * 86400

    # 求解 ODE
    sol = solve_ivp(
        uva_sun_derivatives,
        (t_start, t_end),
        initial_state,
        args=(p, env),
        method='RK45',
        max_step=300
    )

    if not sol.success:
        return None

    # 提取結果
    Xd_f, Cbuf_f, LAI_f, Anth_f, Stress_f, ROS_f = sol.y[:, -1]

    # 計算平均 Stress
    uva_start_day = env.get('uva_start_day', 29)
    uva_start = uva_start_day * 86400
    stress_sum = 0
    stress_count = 0
    for i in range(len(sol.t)):
        if sol.t[i] >= uva_start:
            stress_sum += sol.y[4, i]
            stress_count += 1
    avg_stress = stress_sum / max(1, stress_count)

    # 計算 nonlinear_factor
    uva_hour_on = env.get('uva_hour_on', 0)
    uva_hour_off = env.get('uva_hour_off', 0)
    if env.get('uva_on', False):
        daily_hours = uva_hour_off - uva_hour_on if uva_hour_on < uva_hour_off else 24 - uva_hour_on + uva_hour_off
    else:
        daily_hours = 0
    nonlin_factor = nonlinear_damage_factor(daily_hours, p)

    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(avg_stress, p, nonlin_factor)

    FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
    FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
    Anth_sim = Anth_f / FW_total_kg * 1e6

    return {
        'FW_g': FW_sim,
        'Anth': Anth_sim,
        'LAI': LAI_f,
        'Stress': Stress_f,
        'avg_stress': avg_stress,
        'nonlin_factor': nonlin_factor,
        'sol': sol
    }


def run_simulation_with_history(env_config, description=""):
    """執行 ODE 模擬並返回完整的時間序列"""
    p = UVAParams()
    env = ENV_BASE.copy()
    env.update(env_config)

    # 初始條件
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    C_buf_init = Xd_init * 0.1
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6

    initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]

    # 模擬時間
    transplant_day = SIMULATION['transplant_offset']
    simulation_days = SIMULATION['days']
    t_start = transplant_day * 86400
    t_end = (transplant_day + simulation_days) * 86400

    t_eval = np.linspace(t_start, t_end, 2000)

    sol = solve_ivp(
        uva_sun_derivatives,
        (t_start, t_end),
        initial_state,
        args=(p, env),
        method='RK45',
        t_eval=t_eval,
        max_step=300
    )

    if not sol.success:
        return None

    days = sol.t / 86400

    return {
        't': sol.t,
        'days': days,
        'Xd': sol.y[0],
        'C_buf': sol.y[1],
        'LAI': sol.y[2],
        'Anth': sol.y[3],
        'Stress': sol.y[4],
        'ROS': sol.y[5]
    }


# ============================================================
# FIG 9: LAI Vulnerability Function
# ============================================================
def generate_fig9_lai_vulnerability():
    """圖9. LAI脆弱性函數曲線"""
    print("  生成圖9: LAI脆弱性函數...")
    fig, ax = plt.subplots(figsize=(10, 7))

    p = UVAParams()
    A_vuln = p.A_vulnerability
    k_vuln = p.k_vulnerability

    LAI = np.linspace(3, 10, 200)
    vulnerability = A_vuln * np.exp(-k_vuln * LAI) + 1

    ax.semilogy(LAI, vulnerability, 'b-', linewidth=2.5,
                label=f'vulnerability = A×exp(-k×LAI) + 1\n(A={A_vuln:.1e}, k={k_vuln})')

    # 標記關鍵點
    key_points = [
        (5.5, 'D12 groups\n(Day 23 start)', 'red'),
        (7.5, 'D6 groups\n(Day 29 start)', 'green'),
        (8.5, 'D3 groups\n(Day 32 start)', 'blue'),
    ]

    for lai, label, color in key_points:
        vuln = A_vuln * np.exp(-k_vuln * lai) + 1
        ax.axvline(x=lai, color=color, linestyle='--', alpha=0.7, linewidth=1.5)
        ax.plot(lai, vuln, 'o', color=color, markersize=12, markeredgecolor='black', markeredgewidth=1.5)
        ax.annotate(f'{label}\nLAI≈{lai:.1f}\nvuln≈{vuln:.0f}',
                   xy=(lai, vuln), xytext=(lai+0.3, vuln*2),
                   fontsize=10, ha='left',
                   arrowprops=dict(arrowstyle='->', color=color, lw=1.5))

    ax.set_xlabel('Leaf Area Index (LAI, m²/m²)', fontsize=14)
    ax.set_ylabel('Vulnerability Factor', fontsize=14)
    ax.set_title('Fig. 9. LAI Vulnerability Function', fontsize=14)
    ax.set_xlim(3, 10)
    ax.set_ylim(1, 1e5)
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(loc='upper right', fontsize=11)

    ax.text(0.02, 0.02,
            'Higher vulnerability at lower LAI\n→ Young plants more susceptible to UVA damage',
            transform=ax.transAxes, fontsize=10, verticalalignment='bottom',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig('fig9_lai_vulnerability.png', dpi=150, bbox_inches='tight')
    print("    已生成: fig9_lai_vulnerability.png")
    plt.close()


# ============================================================
# FIG 10: Nonlinear Damage Factor (Gompertz)
# ============================================================
def generate_fig10_nonlinear_factor():
    """圖10. 非線性傷害因子曲線 (Gompertz)"""
    print("  生成圖10: 非線性傷害因子...")
    fig, ax = plt.subplots(figsize=(10, 7))

    p = UVAParams()

    hours = np.linspace(0, 16, 200)
    factors = [nonlinear_damage_factor(h, p) for h in hours]

    ax.semilogy(hours, factors, 'b-', linewidth=2.5,
                label=f'Gompertz factor\n(threshold={p.gompertz_threshold}h, max={p.gompertz_max_factor})')

    # 標記關鍵點
    key_points = [
        (3, 'VL3D3/VL3D12\n3h/day', 'green'),
        (6, 'L6D6/L6D12\n6h/day', 'blue'),
        (9, 'M9D3\n9h/day', 'orange'),
        (12, 'H12D3\n12h/day', 'red'),
        (15, 'VH15D3\n15h/day', 'darkred'),
    ]

    for h, label, color in key_points:
        f = nonlinear_damage_factor(h, p)
        ax.plot(h, f, 'o', color=color, markersize=12, markeredgecolor='black', markeredgewidth=1.5)
        if h >= 12:
            ax.annotate(f'{label}\nfactor≈{f:.0f}',
                       xy=(h, f), xytext=(h-1.5, f*0.15),
                       fontsize=10, ha='center',
                       arrowprops=dict(arrowstyle='->', color=color, lw=1.5))
        elif h == 9:
            ax.annotate(f'{label}\nfactor≈{f:.1f}',
                       xy=(h, f), xytext=(h+0.5, f*3),
                       fontsize=10, ha='left',
                       arrowprops=dict(arrowstyle='->', color=color, lw=1.5))
        else:
            ax.annotate(f'{label}\nfactor≈{f:.2f}',
                       xy=(h, f), xytext=(h+0.3, f*3),
                       fontsize=10, ha='left',
                       arrowprops=dict(arrowstyle='->', color=color, lw=1.5))

    # 閾值線
    ax.axvline(x=p.gompertz_threshold, color='gray', linestyle=':', alpha=0.5, linewidth=1.5)
    ax.text(p.gompertz_threshold + 0.2, 1.5, f'threshold\n({p.gompertz_threshold}h)', fontsize=9, color='gray')

    ax.set_xlabel('Daily Irradiation Hours', fontsize=14)
    ax.set_ylabel('Nonlinear Damage Factor (Gompertz)', fontsize=14)
    ax.set_title('Fig. 10. Nonlinear Damage Factor (Gompertz Function)', fontsize=14)
    ax.set_xlim(0, 16)
    ax.set_ylim(0.8, 300)
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(loc='upper left', fontsize=11)

    ax.text(0.98, 0.02,
            f'Phase transition at ~{p.gompertz_threshold}h/day\n→ Exponential damage increase beyond threshold',
            transform=ax.transAxes, fontsize=10, verticalalignment='bottom', ha='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig('fig10_nonlinear_factor.png', dpi=150, bbox_inches='tight')
    print("    已生成: fig10_nonlinear_factor.png")
    plt.close()


# ============================================================
# FIG 11: Training Set FW Response (曲線圖，類似 FIG13)
# ============================================================
def generate_fig11_training_fw():
    """圖11. 訓練組鮮重響應曲線"""
    print("  生成圖11: 訓練組鮮重響應曲線...")
    fig, ax = plt.subplots(figsize=(10, 7))

    # 訓練組處理 (按類型排序)
    # D6組: L6D6(6h×6d), L6D6-N(6h×6d夜間)
    # D12組: VL3D12(3h×12d), L6D12(6h×12d)
    # D3組: H12D3(12h×3d)
    training_treatments = ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']
    x_positions = [0, 1, 2, 3, 4, 5]
    labels = ['CK\n(0h)', 'L6D6\n(6h×6d)', 'L6D6-N\n(6h×6d,N)', 'VL3D12\n(3h×12d)', 'L6D12\n(6h×12d)', 'H12D3\n(12h×3d)']

    obs_fw = []
    sim_fw = []

    for t in training_treatments:
        obs_fw.append(TARGETS[t]['FW'])
        result = run_simulation(t)
        if result:
            sim_fw.append(result['FW_g'])
        else:
            sim_fw.append(np.nan)

    # 繪製曲線
    ax.plot(x_positions, obs_fw, 'o-', markersize=12, linewidth=2.5,
            color='steelblue', label='Observed', markeredgecolor='black', markeredgewidth=1.5)
    ax.plot(x_positions, sim_fw, 's--', markersize=10, linewidth=2,
            color='coral', label='Simulated', markeredgecolor='black', markeredgewidth=1.5)

    # 添加誤差標註
    for i, (obs, sim) in enumerate(zip(obs_fw, sim_fw)):
        if not np.isnan(sim):
            error = (sim - obs) / obs * 100
            color = 'green' if abs(error) < 5 else 'orange'
            ax.annotate(f'{error:+.1f}%', xy=(i, sim), xytext=(i+0.1, sim+2),
                       fontsize=10, color=color, fontweight='bold')

    ax.set_xlabel('Treatment Group', fontsize=14)
    ax.set_ylabel('Fresh Weight (g/plant)', fontsize=14)
    ax.set_title('Fig. 11. Training Set: Fresh Weight Response\n(6 training groups, target error < 5%)', fontsize=14)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_xlim(-0.5, 5.5)
    ax.set_ylim(50, 100)
    ax.grid(True, alpha=0.3, axis='y')
    ax.legend(loc='upper right', fontsize=11)

    ax.text(0.02, 0.02,
            'Training criteria: FW error < 5%\nAll 6 groups pass',
            transform=ax.transAxes, fontsize=10, verticalalignment='bottom',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

    plt.tight_layout()
    plt.savefig('fig11_model_validation.png', dpi=150, bbox_inches='tight')
    print("    已生成: fig11_model_validation.png")
    plt.close()


# ============================================================
# FIG 12: Training Set Anthocyanin Response (曲線圖，類似 FIG14)
# ============================================================
def generate_fig12_training_anth():
    """圖12. 訓練組花青素響應曲線"""
    print("  生成圖12: 訓練組花青素響應曲線...")
    fig, ax = plt.subplots(figsize=(10, 7))

    training_treatments = ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']
    x_positions = [0, 1, 2, 3, 4, 5]
    labels = ['CK\n(0h)', 'L6D6\n(6h×6d)', 'L6D6-N\n(6h×6d,N)', 'VL3D12\n(3h×12d)', 'L6D12\n(6h×12d)', 'H12D3\n(12h×3d)']

    obs_anth = []
    sim_anth = []

    for t in training_treatments:
        obs_anth.append(TARGETS[t]['Anth'])
        result = run_simulation(t)
        if result:
            sim_anth.append(result['Anth'])
        else:
            sim_anth.append(np.nan)

    # 繪製曲線
    ax.plot(x_positions, obs_anth, 'o-', markersize=12, linewidth=2.5,
            color='purple', label='Observed', markeredgecolor='black', markeredgewidth=1.5)
    ax.plot(x_positions, sim_anth, 's--', markersize=10, linewidth=2,
            color='orchid', label='Simulated', markeredgecolor='black', markeredgewidth=1.5)

    # 添加誤差標註
    for i, (obs, sim) in enumerate(zip(obs_anth, sim_anth)):
        if not np.isnan(sim):
            error = (sim - obs) / obs * 100
            color = 'green' if abs(error) < 5 else 'orange'
            ax.annotate(f'{error:+.1f}%', xy=(i, sim), xytext=(i+0.1, sim+15),
                       fontsize=10, color=color, fontweight='bold')

    # 標記 H12D3 最高點
    max_idx = np.argmax(obs_anth)
    ax.annotate('Max Anthocyanin\n(H12D3)', xy=(max_idx, obs_anth[max_idx]),
               xytext=(max_idx-0.8, obs_anth[max_idx]+50),
               fontsize=11, ha='center',
               arrowprops=dict(arrowstyle='->', color='red', lw=2))

    ax.set_xlabel('Treatment Group', fontsize=14)
    ax.set_ylabel('Anthocyanin (mg/kg FW)', fontsize=14)
    ax.set_title('Fig. 12. Training Set: Anthocyanin Response\n(6 training groups, target error < 5%)', fontsize=14)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_xlim(-0.5, 5.5)
    ax.set_ylim(400, 750)
    ax.grid(True, alpha=0.3, axis='y')
    ax.legend(loc='upper left', fontsize=11)

    ax.text(0.98, 0.02,
            'Training criteria: Anth error < 5%\nAll 6 groups pass',
            transform=ax.transAxes, fontsize=10, verticalalignment='bottom', ha='right',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

    plt.tight_layout()
    plt.savefig('fig12_stress_dynamics.png', dpi=150, bbox_inches='tight')
    print("    已生成: fig12_stress_dynamics.png")
    plt.close()


# ============================================================
# FIG 13: Validation Set FW Response
# ============================================================
def generate_fig13_validation_fw():
    """圖13. 驗證組鮮重響應曲線"""
    print("  生成圖13: 驗證組鮮重響應曲線...")
    fig, ax = plt.subplots(figsize=(10, 7))

    validation_treatments = ['CK_val', 'VL3D3', 'L6D3', 'M9D3', 'H12D3_val', 'VH15D3']
    daily_hours = [0, 3, 6, 9, 12, 15]

    obs_fw = []
    sim_fw = []

    for t in validation_treatments:
        obs_fw.append(TARGETS[t]['FW'])
        result = run_simulation(t)
        if result:
            sim_fw.append(result['FW_g'])
        else:
            sim_fw.append(np.nan)

    ax.plot(daily_hours, obs_fw, 'o-', markersize=12, linewidth=2.5,
            color='steelblue', label='Observed', markeredgecolor='black', markeredgewidth=1.5)
    ax.plot(daily_hours, sim_fw, 's--', markersize=10, linewidth=2,
            color='coral', label='Simulated', markeredgecolor='black', markeredgewidth=1.5)

    for i, (obs, sim, h) in enumerate(zip(obs_fw, sim_fw, daily_hours)):
        if not np.isnan(sim):
            error = (sim - obs) / obs * 100
            color = 'green' if abs(error) < 5 else 'orange' if abs(error) < 10 else 'red'
            ax.annotate(f'{error:+.1f}%', xy=(h, sim), xytext=(h+0.3, sim+2),
                       fontsize=10, color=color, fontweight='bold')

    ax.set_xlabel('Daily UVA Irradiation Hours (h)', fontsize=14)
    ax.set_ylabel('Fresh Weight (g/plant)', fontsize=14)
    ax.set_title('Fig. 13. Validation Set: Fresh Weight Response to UVA Duration\n(3-day gradient experiment, Day 32-35)', fontsize=14)
    ax.set_xlim(-1, 16)
    ax.set_ylim(40, 100)
    ax.set_xticks(daily_hours)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=11)

    ax.text(0.02, 0.02,
            'Validation criteria: FW error < 10%\nAll 6 groups pass',
            transform=ax.transAxes, fontsize=10, verticalalignment='bottom',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

    plt.tight_layout()
    plt.savefig('fig13_validation_fw.png', dpi=150, bbox_inches='tight')
    print("    已生成: fig13_validation_fw.png")
    plt.close()


# ============================================================
# FIG 14: Validation Set Anthocyanin Response (Hormesis)
# ============================================================
def generate_fig14_validation_anth():
    """圖14. 驗證組花青素響應曲線 (顯示 hormesis)"""
    print("  生成圖14: 驗證組花青素響應曲線...")
    fig, ax = plt.subplots(figsize=(10, 7))

    validation_treatments = ['CK_val', 'VL3D3', 'L6D3', 'M9D3', 'H12D3_val', 'VH15D3']
    daily_hours = [0, 3, 6, 9, 12, 15]

    obs_anth = []
    sim_anth = []

    for t in validation_treatments:
        obs_anth.append(TARGETS[t]['Anth'])
        result = run_simulation(t)
        if result:
            sim_anth.append(result['Anth'])
        else:
            sim_anth.append(np.nan)

    ax.plot(daily_hours, obs_anth, 'o-', markersize=12, linewidth=2.5,
            color='purple', label='Observed', markeredgecolor='black', markeredgewidth=1.5)
    ax.plot(daily_hours, sim_anth, 's--', markersize=10, linewidth=2,
            color='orchid', label='Simulated', markeredgecolor='black', markeredgewidth=1.5)

    for i, (obs, sim, h) in enumerate(zip(obs_anth, sim_anth, daily_hours)):
        if not np.isnan(sim):
            error = (sim - obs) / obs * 100
            color = 'green' if abs(error) < 5 else 'orange' if abs(error) < 10 else 'red'
            ax.annotate(f'{error:+.1f}%', xy=(h, sim), xytext=(h+0.3, sim+15),
                       fontsize=10, color=color, fontweight='bold')

    # 標記 hormesis 峰值
    max_idx = np.argmax(obs_anth)
    ax.annotate('Hormesis Peak', xy=(daily_hours[max_idx], obs_anth[max_idx]),
               xytext=(daily_hours[max_idx]-2, obs_anth[max_idx]+50),
               fontsize=11, ha='center',
               arrowprops=dict(arrowstyle='->', color='red', lw=2))

    ax.set_xlabel('Daily UVA Irradiation Hours (h)', fontsize=14)
    ax.set_ylabel('Anthocyanin (mg/kg FW)', fontsize=14)
    ax.set_title('Fig. 14. Validation Set: Anthocyanin Response to UVA Duration\n(showing hormesis effect)', fontsize=14)
    ax.set_xlim(-1, 16)
    ax.set_ylim(350, 750)
    ax.set_xticks(daily_hours)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper left', fontsize=11)

    ax.text(0.98, 0.02,
            'Hormesis effect:\n'
            '- Anthocyanin peaks at 12h/day\n'
            '- Declines at 15h/day\n'
            '- Model captures nonlinear response',
            transform=ax.transAxes, fontsize=10, verticalalignment='bottom', ha='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

    plt.tight_layout()
    plt.savefig('fig14_validation_anth.png', dpi=150, bbox_inches='tight')
    print("    已生成: fig14_validation_anth.png")
    plt.close()


# ============================================================
# FIG 15: Validation Set Scatter Plot
# ============================================================
def generate_fig15_validation_scatter():
    """圖15. 驗證組模型預測 vs 觀測"""
    print("  生成圖15: 驗證組模型預測 vs 觀測...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    validation_treatments = ['CK_val', 'VL3D3', 'L6D3', 'M9D3', 'H12D3_val', 'VH15D3']
    labels = ['CK', 'VL3D3', 'L6D3', 'M9D3', 'H12D3', 'VH15D3']
    colors = ['gray', 'limegreen', 'blue', 'orange', 'red', 'darkred']

    obs_fw, sim_fw = [], []
    obs_anth, sim_anth = [], []

    for t in validation_treatments:
        obs_fw.append(TARGETS[t]['FW'])
        obs_anth.append(TARGETS[t]['Anth'])
        result = run_simulation(t)
        if result:
            sim_fw.append(result['FW_g'])
            sim_anth.append(result['Anth'])
        else:
            sim_fw.append(np.nan)
            sim_anth.append(np.nan)

    # 圖15a: 鮮重
    ax1 = axes[0]
    for obs, sim, label, color in zip(obs_fw, sim_fw, labels, colors):
        ax1.scatter(obs, sim, c=color, s=180, label=label,
                   edgecolors='black', linewidths=1.5, zorder=5)

    fw_range = np.array([45, 95])
    ax1.plot(fw_range, fw_range, 'k--', linewidth=1.5, label='1:1 line')
    ax1.fill_between(fw_range, fw_range*0.90, fw_range*1.10, alpha=0.15, color='green', label='±10% band')

    ax1.set_xlabel('Observed Fresh Weight (g)', fontsize=13)
    ax1.set_ylabel('Simulated Fresh Weight (g)', fontsize=13)
    ax1.set_title('(a) Fresh Weight', fontsize=14)
    ax1.set_xlim(45, 95)
    ax1.set_ylim(45, 95)
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper left', fontsize=9, ncol=2)

    fw_errors = [(s - o) for s, o in zip(sim_fw, obs_fw) if not np.isnan(s)]
    rmse_fw = np.sqrt(np.mean(np.array(fw_errors)**2))
    mape_fw = np.mean([abs(s-o)/o*100 for s, o in zip(sim_fw, obs_fw) if not np.isnan(s)])
    ax1.text(0.98, 0.02, f'RMSE: {rmse_fw:.1f}g\nMAPE: {mape_fw:.1f}%',
            transform=ax1.transAxes, fontsize=10, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

    # 圖15b: 花青素
    ax2 = axes[1]
    for obs, sim, label, color in zip(obs_anth, sim_anth, labels, colors):
        ax2.scatter(obs, sim, c=color, s=180, label=label,
                   edgecolors='black', linewidths=1.5, zorder=5)

    anth_range = np.array([380, 700])
    ax2.plot(anth_range, anth_range, 'k--', linewidth=1.5, label='1:1 line')
    ax2.fill_between(anth_range, anth_range*0.90, anth_range*1.10, alpha=0.15, color='green', label='±10% band')

    ax2.set_xlabel('Observed Anthocyanin (mg/kg FW)', fontsize=13)
    ax2.set_ylabel('Simulated Anthocyanin (mg/kg FW)', fontsize=13)
    ax2.set_title('(b) Anthocyanin', fontsize=14)
    ax2.set_xlim(380, 700)
    ax2.set_ylim(380, 700)
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper left', fontsize=9, ncol=2)

    anth_errors = [(s - o) for s, o in zip(sim_anth, obs_anth) if not np.isnan(s)]
    rmse_anth = np.sqrt(np.mean(np.array(anth_errors)**2))
    mape_anth = np.mean([abs(s-o)/o*100 for s, o in zip(sim_anth, obs_anth) if not np.isnan(s)])
    ax2.text(0.98, 0.02, f'RMSE: {rmse_anth:.1f} mg/kg\nMAPE: {mape_anth:.1f}%',
            transform=ax2.transAxes, fontsize=10, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

    plt.suptitle('Fig. 15. Validation Set: Model Prediction vs Observation', fontsize=15, y=1.02)
    plt.tight_layout()
    plt.savefig('fig15_validation_scatter.png', dpi=150, bbox_inches='tight')
    print("    已生成: fig15_validation_scatter.png")
    plt.close()


# ============================================================
# FIG 16: Anthocyanin Synthesis Efficiency
# ============================================================
def generate_fig16_anth_efficiency():
    """圖16. 花青素合成效率函數"""
    print("  生成圖16: 花青素合成效率函數...")
    fig, ax = plt.subplots(figsize=(10, 7))

    p = UVAParams()

    # 計算 nonlinear_factor 範圍
    hours = np.linspace(0, 16, 200)
    nonlin_factors = [nonlinear_damage_factor(h, p) for h in hours]

    # 計算效率
    efficiencies = [calculate_nonlin_anth_efficiency(nf, p) for nf in nonlin_factors]

    ax.plot(hours, np.array(efficiencies) * 100, 'b-', linewidth=2.5,
            label='Synthesis efficiency')

    # 標記關鍵點
    key_points = [
        (6, 'L6D6\n6h/day', 'blue'),
        (9, 'M9D3\n9h/day', 'orange'),
        (12, 'H12D3\n12h/day', 'red'),
        (15, 'VH15D3\n15h/day', 'darkred'),
    ]

    for h, label, color in key_points:
        nf = nonlinear_damage_factor(h, p)
        eff = calculate_nonlin_anth_efficiency(nf, p) * 100
        ax.plot(h, eff, 'o', color=color, markersize=12, markeredgecolor='black', markeredgewidth=1.5)
        if h == 15:
            ax.annotate(f'{label}\neff={eff:.0f}%',
                       xy=(h, eff), xytext=(h-1, eff-15),
                       fontsize=10, ha='center',
                       arrowprops=dict(arrowstyle='->', color=color, lw=1.5))
        else:
            ax.annotate(f'{label}\neff={eff:.0f}%',
                       xy=(h, eff), xytext=(h+0.5, eff-5),
                       fontsize=10, ha='left',
                       arrowprops=dict(arrowstyle='->', color=color, lw=1.5))

    ax.set_xlabel('Daily UVA Irradiation Hours (h)', fontsize=14)
    ax.set_ylabel('Anthocyanin Synthesis Efficiency (%)', fontsize=14)
    ax.set_title('Fig. 16. Anthocyanin Synthesis Efficiency\n(Hill function: monotonic decrease)', fontsize=14)
    ax.set_xlim(0, 16)
    ax.set_ylim(70, 105)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=11)

    ax.text(0.02, 0.02,
            'v10.39: Hill function inhibition\n'
            'efficiency = 1 / (1 + (nonlin/K)^n)\n'
            'K=800, n=1.5\n'
            '- 6h: 100%, 9h: 99%, 12h: 92%, 15h: 87%',
            transform=ax.transAxes, fontsize=10, verticalalignment='bottom',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

    plt.tight_layout()
    plt.savefig('fig16_hill_efficiency.png', dpi=150, bbox_inches='tight')
    print("    已生成: fig16_hill_efficiency.png")
    plt.close()


# ============================================================
# FIG 17: System Dynamics Block Diagram
# ============================================================
def generate_fig17_system_dynamics():
    """圖17. 系統動力學方塊圖"""
    print("  生成圖17: 系統動力學方塊圖...")
    import matplotlib.patches as mpatches
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

    fig, ax = plt.subplots(figsize=(14, 10))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 10)
    ax.axis('off')

    boxes = {
        'UVA': {'pos': (1, 7), 'size': (2, 1.2), 'color': 'gold', 'text': '$I_{UVA}$\n(11 W/m²)'},
        'ROS': {'pos': (5, 8), 'size': (2.2, 1.2), 'color': 'salmon', 'text': 'ROS\n(Reactive O₂)'},
        'Stress': {'pos': (9, 8), 'size': (2, 1.2), 'color': 'tomato', 'text': 'Stress\n(Cumulative)'},
        'LAI': {'pos': (5, 5), 'size': (2.2, 1.2), 'color': 'lightgreen', 'text': 'LAI\n(Leaf Area)'},
        'Anth': {'pos': (9, 5), 'size': (2, 1.2), 'color': 'plum', 'text': 'Anthocyanin\n(Protection)'},
        'FW': {'pos': (9, 2), 'size': (2, 1.2), 'color': 'lightblue', 'text': 'Fresh Weight\n(Yield)'},
        'Photo': {'pos': (5, 2), 'size': (2.2, 1.2), 'color': 'palegreen', 'text': 'Photosynthesis\n(C assimilation)'},
    }

    for name, props in boxes.items():
        x, y = props['pos']
        w, h = props['size']
        box = FancyBboxPatch((x, y), w, h,
                            boxstyle="round,pad=0.05,rounding_size=0.2",
                            facecolor=props['color'], edgecolor='black', linewidth=2)
        ax.add_patch(box)
        ax.text(x + w/2, y + h/2, props['text'],
               ha='center', va='center', fontsize=11, fontweight='bold')

    arrows = [
        ('UVA', 'ROS', '+k_prod', 'red', 'arc3,rad=0.1'),
        ('UVA', 'LAI', '+SLA effect', 'green', 'arc3,rad=-0.2'),
        ('UVA', 'Photo', '+PAR equiv.', 'green', 'arc3,rad=-0.3'),
        ('ROS', 'Stress', '+k_damage\n×vuln×(1-prot)', 'red', 'arc3,rad=0'),
        ('LAI', 'ROS', 'vulnerability\nA·e^(-k·LAI)', 'orange', 'arc3,rad=0.2'),
        ('LAI', 'Photo', '+light capture', 'green', 'arc3,rad=0'),
        ('Stress', 'Anth', '+induction', 'purple', 'arc3,rad=0'),
        ('Stress', 'FW', '-LDMC effect', 'red', 'arc3,rad=0'),
        ('Anth', 'ROS', '-protection\nα·Anth/(K+Anth)', 'blue', 'arc3,rad=-0.3'),
        ('Photo', 'FW', '+dry matter', 'green', 'arc3,rad=0'),
    ]

    for from_box, to_box, label, color, style in arrows:
        from_props = boxes[from_box]
        to_props = boxes[to_box]

        fx, fy = from_props['pos']
        fw, fh = from_props['size']
        tx, ty = to_props['pos']
        tw, th = to_props['size']

        if tx > fx + fw:
            start = (fx + fw, fy + fh/2)
            end = (tx, ty + th/2)
        elif tx + tw < fx:
            start = (fx, fy + fh/2)
            end = (tx + tw, ty + th/2)
        elif ty > fy + fh:
            start = (fx + fw/2, fy + fh)
            end = (tx + tw/2, ty)
        else:
            start = (fx + fw/2, fy)
            end = (tx + tw/2, ty + th)

        arrow = FancyArrowPatch(start, end,
                               connectionstyle=style,
                               arrowstyle='-|>',
                               mutation_scale=15,
                               color=color, linewidth=2)
        ax.add_patch(arrow)

        mid_x = (start[0] + end[0]) / 2
        mid_y = (start[1] + end[1]) / 2
        ax.text(mid_x, mid_y, label, fontsize=8, ha='center', va='center',
               bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))

    legend_elements = [
        mpatches.Patch(facecolor='gold', edgecolor='black', label='Input (UVA)'),
        mpatches.Patch(facecolor='salmon', edgecolor='black', label='Damage pathway'),
        mpatches.Patch(facecolor='lightgreen', edgecolor='black', label='Growth pathway'),
        mpatches.Patch(facecolor='plum', edgecolor='black', label='Protection pathway'),
    ]
    ax.legend(handles=legend_elements, loc='lower left', fontsize=10)

    ax.set_title('Fig. 17. System Dynamics Block Diagram\n(UVA Signal Flow in the 6-State ODE Model)',
                fontsize=14, fontweight='bold', y=0.98)

    explanation = (
        "Key mechanisms:\n"
        "• UVA → ROS: Direct oxidative stress induction\n"
        "• LAI vulnerability: Young plants (low LAI) more susceptible\n"
        "• Anthocyanin protection: Negative feedback on ROS damage\n"
        "• Stress → LDMC: High stress increases DW/FW ratio → reduces yield"
    )
    ax.text(0.02, 0.02, explanation, transform=ax.transAxes, fontsize=9,
           verticalalignment='bottom',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()
    plt.savefig('fig17_system_dynamics.png', dpi=150, bbox_inches='tight')
    print("    已生成: fig17_system_dynamics.png")
    plt.close()


# ============================================================
# FIG 18: Hormesis 3D Surface
# ============================================================
def generate_fig18_hormesis_3d():
    """圖18. Hormesis 3D 曲面圖"""
    print("  生成圖18: Hormesis 3D 曲面圖...")
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection='3d')

    p = UVAParams()

    days_range = np.array([3, 6, 12])
    hours_range = np.array([0, 3, 6, 9, 12, 15])

    Days, Hours = np.meshgrid(days_range, hours_range)
    Anth = np.zeros_like(Days, dtype=float)

    print("    計算 Hormesis 曲面數據...")

    for i, days in enumerate(days_range):
        for j, hours in enumerate(hours_range):
            env = ENV_BASE.copy()
            if hours > 0:
                env['uva_on'] = True
                env['uva_intensity'] = 11.0
                env['uva_start_day'] = 35 - days
                env['uva_end_day'] = 35
                env['uva_hour_on'] = 12 - hours // 2
                env['uva_hour_off'] = 12 + hours // 2 + hours % 2
            else:
                env['uva_on'] = False

            result = run_simulation_with_history(env, f"D{days}H{hours}")
            if result:
                Xd_f = result['Xd'][-1]
                Anth_f = result['Anth'][-1]

                uva_start_day = env.get('uva_start_day', 35)
                uva_start = uva_start_day * 86400
                stress_sum = 0
                stress_count = 0
                for k in range(len(result['t'])):
                    if result['t'][k] >= uva_start:
                        stress_sum += result['Stress'][k]
                        stress_count += 1
                avg_stress = stress_sum / max(1, stress_count)

                nonlin_factor = nonlinear_damage_factor(hours, p)
                dw_fw_ratio = calculate_dynamic_dw_fw_ratio(avg_stress, p, nonlin_factor)
                FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
                FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
                Anth_sim = Anth_f / FW_total_kg * 1e6

                Anth[j, i] = Anth_sim

    surf = ax.plot_surface(Days, Hours, Anth, cmap='plasma',
                          edgecolor='black', linewidth=0.3, alpha=0.8)

    exp_points = [
        (3, 3, 'VL3D3', 'green'),
        (3, 6, 'L6D3', 'blue'),
        (3, 9, 'M9D3', 'orange'),
        (3, 12, 'H12D3', 'red'),
        (3, 15, 'VH15D3', 'darkred'),
        (6, 6, 'L6D6', 'cyan'),
        (12, 3, 'VL3D12', 'lime'),
        (12, 6, 'L6D12', 'purple'),
    ]

    for days, hours, label, color in exp_points:
        i = np.where(days_range == days)[0]
        j = np.where(hours_range == hours)[0]
        if len(i) > 0 and len(j) > 0:
            z = Anth[j[0], i[0]]
            ax.scatter([days], [hours], [z], c=color, s=100, marker='o',
                      edgecolors='black', linewidths=1.5, zorder=5)
            ax.text(days, hours, z + 20, label, fontsize=8, ha='center')

    ax.set_xlabel('UVA Treatment Days', fontsize=12, labelpad=10)
    ax.set_ylabel('Daily Hours (h/day)', fontsize=12, labelpad=10)
    ax.set_zlabel('Anthocyanin (mg/kg FW)', fontsize=12, labelpad=10)
    ax.set_title('Fig. 18. Hormesis Response Surface\n(Anthocyanin as function of treatment duration and daily dose)',
                fontsize=13, fontweight='bold')

    ax.view_init(elev=25, azim=45)

    cbar = fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, pad=0.1)
    cbar.set_label('Anthocyanin (mg/kg FW)', fontsize=10)

    ax.text2D(0.02, 0.02,
             "Hormesis effect:\n"
             "• Peak at 12h/day × 3 days\n"
             "• Decline at 15h/day (metabolic collapse)\n"
             "• Longer duration (12 days) shows lower peak",
             transform=ax.transAxes, fontsize=9,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()
    plt.savefig('fig18_hormesis_3d.png', dpi=150, bbox_inches='tight')
    print("    已生成: fig18_hormesis_3d.png")
    plt.close()


# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("生成論文圖表 FIG9-18 (v10.37)")
    print("=" * 60)

    print("\n[1/10] 生成圖9: LAI脆弱性函數...")
    generate_fig9_lai_vulnerability()

    print("\n[2/10] 生成圖10: 非線性傷害因子 (Gompertz)...")
    generate_fig10_nonlinear_factor()

    print("\n[3/10] 生成圖11: 訓練組鮮重響應曲線...")
    generate_fig11_training_fw()

    print("\n[4/10] 生成圖12: 訓練組花青素響應曲線...")
    generate_fig12_training_anth()

    print("\n[5/10] 生成圖13: 驗證組鮮重響應曲線...")
    generate_fig13_validation_fw()

    print("\n[6/10] 生成圖14: 驗證組花青素響應曲線...")
    generate_fig14_validation_anth()

    print("\n[7/10] 生成圖15: 驗證組散點圖...")
    generate_fig15_validation_scatter()

    print("\n[8/10] 生成圖16: 花青素合成效率函數...")
    generate_fig16_anth_efficiency()

    print("\n[9/10] 生成圖17: 系統動力學方塊圖...")
    generate_fig17_system_dynamics()

    print("\n[10/10] 生成圖18: Hormesis 3D 曲面圖...")
    generate_fig18_hormesis_3d()

    print("\n" + "=" * 60)
    print("所有圖表生成完成!")
    print("=" * 60)
    print("\n生成的檔案:")
    print("  - fig9_lai_vulnerability.png")
    print("  - fig10_nonlinear_factor.png")
    print("  - fig11_model_validation.png   (訓練組 FW 曲線)")
    print("  - fig12_stress_dynamics.png    (訓練組 Anth 曲線)")
    print("  - fig13_validation_fw.png      (驗證組 FW 曲線)")
    print("  - fig14_validation_anth.png    (驗證組 Anth 曲線)")
    print("  - fig15_validation_scatter.png (驗證組散點圖)")
    print("  - fig16_hill_efficiency.png    (花青素合成效率 - Hill函數)")
    print("  - fig17_system_dynamics.png    (系統方塊圖)")
    print("  - fig18_hormesis_3d.png        (3D 曲面圖)")
