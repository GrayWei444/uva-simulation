"""
================================================================================
生成驗證組圖表 (圖13-16)
================================================================================
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

    # 計算平均 Stress (與主模型一致 - 只計算 UVA 照射期間)
    uva_start_day = env.get('uva_start_day', 29)
    uva_start = uva_start_day * 86400
    stress_sum = 0
    stress_count = 0
    for i in range(len(sol.t)):
        if sol.t[i] >= uva_start:
            stress_sum += sol.y[4, i]
            stress_count += 1
    avg_stress = stress_sum / max(1, stress_count)

    # 計算 nonlinear_factor (根據每日照射時數)
    uva_hour_on = env.get('uva_hour_on', 0)
    uva_hour_off = env.get('uva_hour_off', 0)
    if env.get('uva_on', False):
        daily_hours = uva_hour_off - uva_hour_on if uva_hour_on < uva_hour_off else 24 - uva_hour_on + uva_hour_off
    else:
        daily_hours = 0
    nonlin_factor = nonlinear_damage_factor(daily_hours, p)

    # 使用 avg_stress 和 nonlin_factor 計算 dw_fw_ratio (與主模型一致)
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
        'nonlin_factor': nonlin_factor
    }


def generate_fig13_validation_fw():
    """圖13. 驗證組鮮重響應曲線"""
    print("  生成圖13: 驗證組鮮重響應曲線...")

    fig, ax = plt.subplots(figsize=(10, 7))

    # 驗證組處理 (排序按每日照射時數)
    validation_treatments = ['CK_val', 'VL3D3', 'L6D3', 'M9D3', 'H12D3_val', 'VH15D3']
    daily_hours = [0, 3, 6, 9, 12, 15]

    # 收集數據
    obs_fw = []
    sim_fw = []

    for t in validation_treatments:
        # 觀測值
        obs_fw.append(TARGETS[t]['FW'])

        # 模擬值
        result = run_simulation(t)
        if result:
            sim_fw.append(result['FW_g'])
        else:
            sim_fw.append(np.nan)

    # 繪製曲線
    ax.plot(daily_hours, obs_fw, 'o-', markersize=12, linewidth=2.5,
            color='steelblue', label='Observed', markeredgecolor='black', markeredgewidth=1.5)
    ax.plot(daily_hours, sim_fw, 's--', markersize=10, linewidth=2,
            color='coral', label='Simulated', markeredgecolor='black', markeredgewidth=1.5)

    # 添加誤差標註
    for i, (obs, sim, h) in enumerate(zip(obs_fw, sim_fw, daily_hours)):
        if not np.isnan(sim):
            error = (sim - obs) / obs * 100
            color = 'green' if abs(error) < 5 else 'orange'
            ax.annotate(f'{error:+.1f}%', xy=(h, sim), xytext=(h+0.3, sim+2),
                       fontsize=9, color=color)

    # 標記訓練組參考 (H12D3)
    ax.axhline(y=TARGETS['H12D3']['FW'], color='red', linestyle=':', alpha=0.5, linewidth=1.5)
    ax.text(1, TARGETS['H12D3']['FW']+1, 'Training H12D3', fontsize=9, color='red', alpha=0.7)

    ax.set_xlabel('Daily UVA Irradiation Hours (h)', fontsize=14)
    ax.set_ylabel('Fresh Weight (g/plant)', fontsize=14)
    ax.set_title('Fig. 13. Validation Set: Fresh Weight Response to UVA Duration\n(3-day gradient experiment, Day 32-35)', fontsize=14)
    ax.set_xlim(-1, 16)
    ax.set_ylim(40, 100)
    ax.set_xticks(daily_hours)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=11)

    # 添加說明
    ax.text(0.02, 0.02,
            'Key observation:\n'
            '- FW decreases with increasing daily hours\n'
            '- Nonlinear damage threshold at ~9h/day',
            transform=ax.transAxes, fontsize=10, verticalalignment='bottom',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig('fig13_validation_fw.png', dpi=150, bbox_inches='tight')
    print("    已生成: fig13_validation_fw.png")
    plt.close()


def generate_fig14_validation_anth():
    """圖14. 驗證組花青素響應曲線 (顯示 hormesis)"""
    print("  生成圖14: 驗證組花青素響應曲線...")

    fig, ax = plt.subplots(figsize=(10, 7))

    # 驗證組處理
    validation_treatments = ['CK_val', 'VL3D3', 'L6D3', 'M9D3', 'H12D3_val', 'VH15D3']
    daily_hours = [0, 3, 6, 9, 12, 15]

    # 收集數據
    obs_anth = []
    sim_anth = []

    for t in validation_treatments:
        obs_anth.append(TARGETS[t]['Anth'])
        result = run_simulation(t)
        if result:
            sim_anth.append(result['Anth'])
        else:
            sim_anth.append(np.nan)

    # 繪製曲線
    ax.plot(daily_hours, obs_anth, 'o-', markersize=12, linewidth=2.5,
            color='purple', label='Observed', markeredgecolor='black', markeredgewidth=1.5)
    ax.plot(daily_hours, sim_anth, 's--', markersize=10, linewidth=2,
            color='orchid', label='Simulated', markeredgecolor='black', markeredgewidth=1.5)

    # 添加誤差標註
    for i, (obs, sim, h) in enumerate(zip(obs_anth, sim_anth, daily_hours)):
        if not np.isnan(sim):
            error = (sim - obs) / obs * 100
            color = 'green' if abs(error) < 10 else 'orange'
            ax.annotate(f'{error:+.1f}%', xy=(h, sim), xytext=(h+0.3, sim+10),
                       fontsize=9, color=color)

    # 標記 hormesis 峰值 (H12D3)
    max_idx = np.argmax(obs_anth)
    ax.annotate('Hormesis Peak', xy=(daily_hours[max_idx], obs_anth[max_idx]),
               xytext=(daily_hours[max_idx]-2, obs_anth[max_idx]+40),
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

    # 添加說明
    ax.text(0.98, 0.02,
            'Hormesis effect:\n'
            '- Anthocyanin peaks at 12h/day\n'
            '- Declines at 15h/day (metabolic collapse)\n'
            '- Model captures nonlinear response',
            transform=ax.transAxes, fontsize=10, verticalalignment='bottom', ha='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig('fig14_validation_anth.png', dpi=150, bbox_inches='tight')
    print("    已生成: fig14_validation_anth.png")
    plt.close()


def generate_fig15_validation_scatter():
    """圖15. 驗證組模型預測 vs 觀測"""
    print("  生成圖15: 驗證組模型預測 vs 觀測...")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # 驗證組處理
    validation_treatments = ['CK_val', 'VL3D3', 'L6D3', 'M9D3', 'H12D3_val', 'VH15D3']
    labels = ['CK', 'VL3D3', 'L6D3', 'M9D3', 'H12D3', 'VH15D3']
    colors = ['gray', 'green', 'blue', 'orange', 'red', 'darkred']

    # 收集數據
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
    for i, (obs, sim, label, color) in enumerate(zip(obs_fw, sim_fw, labels, colors)):
        ax1.scatter(obs, sim, c=color, s=150, label=label,
                   edgecolors='black', linewidths=1.5, zorder=5)

    # 1:1 線和誤差帶
    fw_range = np.array([45, 95])
    ax1.plot(fw_range, fw_range, 'k--', linewidth=1.5, label='1:1 line')
    ax1.fill_between(fw_range, fw_range*0.90, fw_range*1.10, alpha=0.15, color='gray', label='±10% band')

    ax1.set_xlabel('Observed Fresh Weight (g)', fontsize=13)
    ax1.set_ylabel('Simulated Fresh Weight (g)', fontsize=13)
    ax1.set_title('(a) Fresh Weight', fontsize=14)
    ax1.set_xlim(45, 95)
    ax1.set_ylim(45, 95)
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper left', fontsize=9, ncol=2)

    # 計算 RMSE
    fw_errors = [(s - o) for s, o in zip(sim_fw, obs_fw) if not np.isnan(s)]
    rmse_fw = np.sqrt(np.mean(np.array(fw_errors)**2))
    mape_fw = np.mean([abs(s-o)/o*100 for s, o in zip(sim_fw, obs_fw) if not np.isnan(s)])
    ax1.text(0.98, 0.02, f'RMSE: {rmse_fw:.1f}g\nMAPE: {mape_fw:.1f}%',
            transform=ax1.transAxes, fontsize=10, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

    # 圖15b: 花青素
    ax2 = axes[1]
    for i, (obs, sim, label, color) in enumerate(zip(obs_anth, sim_anth, labels, colors)):
        ax2.scatter(obs, sim, c=color, s=150, label=label,
                   edgecolors='black', linewidths=1.5, zorder=5)

    # 1:1 線和誤差帶
    anth_range = np.array([380, 700])
    ax2.plot(anth_range, anth_range, 'k--', linewidth=1.5, label='1:1 line')
    ax2.fill_between(anth_range, anth_range*0.90, anth_range*1.10, alpha=0.15, color='gray', label='±10% band')

    ax2.set_xlabel('Observed Anthocyanin (mg/kg FW)', fontsize=13)
    ax2.set_ylabel('Simulated Anthocyanin (mg/kg FW)', fontsize=13)
    ax2.set_title('(b) Anthocyanin', fontsize=14)
    ax2.set_xlim(380, 700)
    ax2.set_ylim(380, 700)
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper left', fontsize=9, ncol=2)

    # 計算 RMSE
    anth_errors = [(s - o) for s, o in zip(sim_anth, obs_anth) if not np.isnan(s)]
    rmse_anth = np.sqrt(np.mean(np.array(anth_errors)**2))
    mape_anth = np.mean([abs(s-o)/o*100 for s, o in zip(sim_anth, obs_anth) if not np.isnan(s)])
    ax2.text(0.98, 0.02, f'RMSE: {rmse_anth:.1f} ppm\nMAPE: {mape_anth:.1f}%',
            transform=ax2.transAxes, fontsize=10, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

    plt.suptitle('Fig. 15. Validation Set: Model Prediction vs Observation\n(Independent 3-day gradient experiment)', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig('fig15_validation_scatter.png', dpi=150, bbox_inches='tight')
    print("    已生成: fig15_validation_scatter.png")
    plt.close()


def generate_fig16_asymmetric_gaussian():
    """圖16. 非對稱高斯效率函數曲線"""
    print("  生成圖16: 非對稱高斯效率函數...")

    fig, ax = plt.subplots(figsize=(10, 7))

    p = UVAParams()

    # 計算非線性因子與效率的關係
    hours = np.linspace(0, 16, 200)
    nonlin_factors = np.array([nonlinear_damage_factor(h, p) for h in hours])
    efficiencies = np.array([calculate_nonlin_anth_efficiency(nf, p) for nf in nonlin_factors])

    # 繪製效率曲線
    ax.plot(hours, efficiencies * 100, 'b-', linewidth=2.5, label='Efficiency')

    # 標記關鍵點
    key_points = [
        (3, 'VL3D3', 'green'),
        (6, 'L6D3', 'blue'),
        (9, 'M9D3', 'orange'),
        (12, 'H12D3', 'red'),
        (15, 'VH15D3', 'darkred'),
    ]

    for h, label, color in key_points:
        nf = nonlinear_damage_factor(h, p)
        eff = calculate_nonlin_anth_efficiency(nf, p)
        ax.plot(h, eff * 100, 'o', color=color, markersize=14,
               markeredgecolor='black', markeredgewidth=2, zorder=5)

        # 標註
        if h == 9:
            ax.annotate(f'{label}\nnf={nf:.1f}\neff={eff*100:.0f}%',
                       xy=(h, eff*100), xytext=(h-2, eff*100-15),
                       fontsize=10, ha='center',
                       arrowprops=dict(arrowstyle='->', color=color, lw=1.5))
        else:
            ax.annotate(f'{label}\nnf={nf:.1f}\neff={eff*100:.0f}%',
                       xy=(h, eff*100), xytext=(h+0.5, eff*100+5),
                       fontsize=10, ha='left')

    # 標記抑制中心
    ax.axvline(x=9, color='orange', linestyle='--', alpha=0.5, linewidth=1.5)
    ax.axhline(y=50, color='orange', linestyle='--', alpha=0.5, linewidth=1.5)

    ax.set_xlabel('Daily UVA Irradiation Hours (h)', fontsize=14)
    ax.set_ylabel('Anthocyanin Synthesis Efficiency (%)', fontsize=14)
    ax.set_title('Fig. 16. Asymmetric Gaussian Efficiency Function\n(Stress-induced anthocyanin synthesis)', fontsize=14)
    ax.set_xlim(0, 16)
    ax.set_ylim(40, 110)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=11)

    # 添加公式說明
    formula_text = (
        r'$\sigma = \sigma_L \cdot \frac{1-t}{2} + \sigma_R \cdot \frac{1+t}{2}$' + '\n'
        r'$t = \tanh((nf - center) / scale)$' + '\n'
        r'$eff = 1 - max\_inhib \cdot \exp(-\frac{(nf-center)^2}{2\sigma^2})$' + '\n\n'
        r'Parameters:' + '\n'
        r'$center=70, \sigma_L=15, \sigma_R=50$' + '\n'
        r'$max\_inhib=0.50, scale=20$'
    )
    ax.text(0.02, 0.02, formula_text,
            transform=ax.transAxes, fontsize=9, verticalalignment='bottom',
            family='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

    # 添加生物學解釋
    ax.text(0.98, 0.98,
            'Biological interpretation:\n'
            '- Low hours (3-6h): Full efficiency\n'
            '- 9h/day: Maximum inhibition (50%)\n'
            '- 12-15h: Recovery (adaptation)',
            transform=ax.transAxes, fontsize=10, verticalalignment='top', ha='right',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7))

    plt.tight_layout()
    plt.savefig('fig16_asymmetric_gaussian.png', dpi=150, bbox_inches='tight')
    print("    已生成: fig16_asymmetric_gaussian.png")
    plt.close()


if __name__ == "__main__":
    print("=" * 60)
    print("生成驗證組圖表 (圖13-16)")
    print("=" * 60)

    print("\n[1/4] 圖13: 驗證組鮮重響應曲線")
    generate_fig13_validation_fw()

    print("\n[2/4] 圖14: 驗證組花青素響應曲線")
    generate_fig14_validation_anth()

    print("\n[3/4] 圖15: 驗證組模型預測 vs 觀測")
    generate_fig15_validation_scatter()

    print("\n[4/4] 圖16: 非對稱高斯效率函數")
    generate_fig16_asymmetric_gaussian()

    print("\n" + "=" * 60)
    print("所有驗證組圖表生成完成!")
    print("=" * 60)
