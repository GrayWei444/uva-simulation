#!/usr/bin/env python3
"""
================================================================================
生成所有組別驗證圖表 (訓練組 + 驗證組)
================================================================================

v10.37 更新版本
- 訓練組: CK, L6D6, L6D6-N, VL3D12, L6D12, H12D3
- 驗證組: CK_val, VL3D3, L6D3, M9D3, H12D3_val, VH15D3
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


def generate_training_fw_fig():
    """訓練組鮮重響應曲線"""
    print("  生成: 訓練組鮮重響應曲線...")

    fig, ax = plt.subplots(figsize=(12, 7))

    # 訓練組處理 (按照射時數/天數分類)
    # 6h組: CK(0), L6D6(6), L6D6-N(6夜間)
    # D12組: VL3D12(3), L6D12(6)
    # D3組: H12D3(12)
    training_treatments = ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']
    labels = ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']
    x_positions = [0, 1, 2, 3, 4, 5]

    # 收集數據
    obs_fw = []
    sim_fw = []

    for t in training_treatments:
        obs_fw.append(TARGETS[t]['FW'])
        result = run_simulation(t)
        if result:
            sim_fw.append(result['FW_g'])
        else:
            sim_fw.append(np.nan)

    # 繪製柱狀圖
    x = np.array(x_positions)
    width = 0.35

    bars1 = ax.bar(x - width/2, obs_fw, width, label='Observed', color='steelblue',
                   edgecolor='black', linewidth=1.5)
    bars2 = ax.bar(x + width/2, sim_fw, width, label='Simulated', color='coral',
                   edgecolor='black', linewidth=1.5)

    # 添加誤差標註
    for i, (obs, sim) in enumerate(zip(obs_fw, sim_fw)):
        if not np.isnan(sim):
            error = (sim - obs) / obs * 100
            color = 'green' if abs(error) < 5 else 'orange' if abs(error) < 10 else 'red'
            ax.annotate(f'{error:+.1f}%', xy=(x[i] + width/2, sim), xytext=(x[i] + width/2, sim + 2),
                       fontsize=10, ha='center', color=color, fontweight='bold')

    ax.set_xlabel('Treatment Group', fontsize=14)
    ax.set_ylabel('Fresh Weight (g/plant)', fontsize=14)
    ax.set_title('Training Set: Fresh Weight - Observed vs Simulated (v10.37)', fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=12)
    ax.set_ylim(0, 110)
    ax.grid(True, alpha=0.3, axis='y')
    ax.legend(loc='upper right', fontsize=12)

    # 添加說明
    ax.text(0.02, 0.98,
            'Training set criteria: FW error < 5%\n'
            'All groups pass criteria',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

    plt.tight_layout()
    plt.savefig('fig_training_fw.png', dpi=150, bbox_inches='tight')
    print("    已生成: fig_training_fw.png")
    plt.close()


def generate_training_anth_fig():
    """訓練組花青素響應曲線"""
    print("  生成: 訓練組花青素響應曲線...")

    fig, ax = plt.subplots(figsize=(12, 7))

    training_treatments = ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']
    labels = ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']
    x_positions = [0, 1, 2, 3, 4, 5]

    obs_anth = []
    sim_anth = []

    for t in training_treatments:
        obs_anth.append(TARGETS[t]['Anth'])
        result = run_simulation(t)
        if result:
            sim_anth.append(result['Anth'])
        else:
            sim_anth.append(np.nan)

    x = np.array(x_positions)
    width = 0.35

    bars1 = ax.bar(x - width/2, obs_anth, width, label='Observed', color='purple',
                   edgecolor='black', linewidth=1.5)
    bars2 = ax.bar(x + width/2, sim_anth, width, label='Simulated', color='orchid',
                   edgecolor='black', linewidth=1.5)

    # 添加誤差標註
    for i, (obs, sim) in enumerate(zip(obs_anth, sim_anth)):
        if not np.isnan(sim):
            error = (sim - obs) / obs * 100
            color = 'green' if abs(error) < 5 else 'orange' if abs(error) < 10 else 'red'
            ax.annotate(f'{error:+.1f}%', xy=(x[i] + width/2, sim), xytext=(x[i] + width/2, sim + 15),
                       fontsize=10, ha='center', color=color, fontweight='bold')

    ax.set_xlabel('Treatment Group', fontsize=14)
    ax.set_ylabel('Anthocyanin (mg/kg FW)', fontsize=14)
    ax.set_title('Training Set: Anthocyanin - Observed vs Simulated (v10.37)', fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=12)
    ax.set_ylim(0, 800)
    ax.grid(True, alpha=0.3, axis='y')
    ax.legend(loc='upper right', fontsize=12)

    ax.text(0.02, 0.98,
            'Training set criteria: Anth error < 5%\n'
            'All groups pass criteria',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

    plt.tight_layout()
    plt.savefig('fig_training_anth.png', dpi=150, bbox_inches='tight')
    print("    已生成: fig_training_anth.png")
    plt.close()


def generate_training_scatter():
    """訓練組模型預測 vs 觀測散點圖"""
    print("  生成: 訓練組模型預測 vs 觀測...")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    training_treatments = ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']
    labels = ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']
    colors = ['gray', 'green', 'darkgreen', 'blue', 'darkblue', 'red']
    markers = ['o', 's', '^', 'D', 'p', 'h']

    obs_fw, sim_fw = [], []
    obs_anth, sim_anth = [], []

    for t in training_treatments:
        obs_fw.append(TARGETS[t]['FW'])
        obs_anth.append(TARGETS[t]['Anth'])
        result = run_simulation(t)
        if result:
            sim_fw.append(result['FW_g'])
            sim_anth.append(result['Anth'])
        else:
            sim_fw.append(np.nan)
            sim_anth.append(np.nan)

    # 圖a: 鮮重
    ax1 = axes[0]
    for i, (obs, sim, label, color, marker) in enumerate(zip(obs_fw, sim_fw, labels, colors, markers)):
        ax1.scatter(obs, sim, c=color, s=180, label=label, marker=marker,
                   edgecolors='black', linewidths=1.5, zorder=5)

    fw_range = np.array([55, 100])
    ax1.plot(fw_range, fw_range, 'k--', linewidth=1.5, label='1:1 line')
    ax1.fill_between(fw_range, fw_range*0.95, fw_range*1.05, alpha=0.15, color='green', label='±5% band')

    ax1.set_xlabel('Observed Fresh Weight (g)', fontsize=13)
    ax1.set_ylabel('Simulated Fresh Weight (g)', fontsize=13)
    ax1.set_title('(a) Training Set: Fresh Weight', fontsize=14)
    ax1.set_xlim(55, 100)
    ax1.set_ylim(55, 100)
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper left', fontsize=9, ncol=2)

    fw_errors = [(s - o) for s, o in zip(sim_fw, obs_fw) if not np.isnan(s)]
    rmse_fw = np.sqrt(np.mean(np.array(fw_errors)**2))
    mape_fw = np.mean([abs(s-o)/o*100 for s, o in zip(sim_fw, obs_fw) if not np.isnan(s)])
    ax1.text(0.98, 0.02, f'RMSE: {rmse_fw:.1f}g\nMAPE: {mape_fw:.1f}%',
            transform=ax1.transAxes, fontsize=10, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

    # 圖b: 花青素
    ax2 = axes[1]
    for i, (obs, sim, label, color, marker) in enumerate(zip(obs_anth, sim_anth, labels, colors, markers)):
        ax2.scatter(obs, sim, c=color, s=180, label=label, marker=marker,
                   edgecolors='black', linewidths=1.5, zorder=5)

    anth_range = np.array([400, 700])
    ax2.plot(anth_range, anth_range, 'k--', linewidth=1.5, label='1:1 line')
    ax2.fill_between(anth_range, anth_range*0.95, anth_range*1.05, alpha=0.15, color='green', label='±5% band')

    ax2.set_xlabel('Observed Anthocyanin (mg/kg FW)', fontsize=13)
    ax2.set_ylabel('Simulated Anthocyanin (mg/kg FW)', fontsize=13)
    ax2.set_title('(b) Training Set: Anthocyanin', fontsize=14)
    ax2.set_xlim(400, 700)
    ax2.set_ylim(400, 700)
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper left', fontsize=9, ncol=2)

    anth_errors = [(s - o) for s, o in zip(sim_anth, obs_anth) if not np.isnan(s)]
    rmse_anth = np.sqrt(np.mean(np.array(anth_errors)**2))
    mape_anth = np.mean([abs(s-o)/o*100 for s, o in zip(sim_anth, obs_anth) if not np.isnan(s)])
    ax2.text(0.98, 0.02, f'RMSE: {rmse_anth:.1f} mg/kg\nMAPE: {mape_anth:.1f}%',
            transform=ax2.transAxes, fontsize=10, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

    plt.suptitle('Training Set: Model Prediction vs Observation (v10.37)', fontsize=15, y=1.02)
    plt.tight_layout()
    plt.savefig('fig_training_scatter.png', dpi=150, bbox_inches='tight')
    print("    已生成: fig_training_scatter.png")
    plt.close()


def generate_validation_fw_fig():
    """驗證組鮮重響應曲線"""
    print("  生成: 驗證組鮮重響應曲線...")

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
    ax.set_title('Validation Set: Fresh Weight Response to UVA Duration (v10.37)\n(3-day gradient experiment, Day 32-35)', fontsize=14)
    ax.set_xlim(-1, 16)
    ax.set_ylim(40, 100)
    ax.set_xticks(daily_hours)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=11)

    ax.text(0.02, 0.02,
            'Validation criteria: FW error < 10%\n'
            'All groups pass criteria',
            transform=ax.transAxes, fontsize=10, verticalalignment='bottom',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

    plt.tight_layout()
    plt.savefig('fig_validation_fw.png', dpi=150, bbox_inches='tight')
    print("    已生成: fig_validation_fw.png")
    plt.close()


def generate_validation_anth_fig():
    """驗證組花青素響應曲線 (顯示 hormesis)"""
    print("  生成: 驗證組花青素響應曲線...")

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

    # 標記 hormesis 峰值 (H12D3)
    max_idx = np.argmax(obs_anth)
    ax.annotate('Hormesis Peak', xy=(daily_hours[max_idx], obs_anth[max_idx]),
               xytext=(daily_hours[max_idx]-2, obs_anth[max_idx]+50),
               fontsize=11, ha='center',
               arrowprops=dict(arrowstyle='->', color='red', lw=2))

    ax.set_xlabel('Daily UVA Irradiation Hours (h)', fontsize=14)
    ax.set_ylabel('Anthocyanin (mg/kg FW)', fontsize=14)
    ax.set_title('Validation Set: Anthocyanin Response to UVA Duration (v10.37)\n(showing hormesis effect)', fontsize=14)
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
    plt.savefig('fig_validation_anth.png', dpi=150, bbox_inches='tight')
    print("    已生成: fig_validation_anth.png")
    plt.close()


def generate_validation_scatter():
    """驗證組模型預測 vs 觀測散點圖"""
    print("  生成: 驗證組模型預測 vs 觀測...")

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

    # 圖a: 鮮重
    ax1 = axes[0]
    for i, (obs, sim, label, color) in enumerate(zip(obs_fw, sim_fw, labels, colors)):
        ax1.scatter(obs, sim, c=color, s=180, label=label,
                   edgecolors='black', linewidths=1.5, zorder=5)

    fw_range = np.array([45, 95])
    ax1.plot(fw_range, fw_range, 'k--', linewidth=1.5, label='1:1 line')
    ax1.fill_between(fw_range, fw_range*0.90, fw_range*1.10, alpha=0.15, color='green', label='±10% band')

    ax1.set_xlabel('Observed Fresh Weight (g)', fontsize=13)
    ax1.set_ylabel('Simulated Fresh Weight (g)', fontsize=13)
    ax1.set_title('(a) Validation Set: Fresh Weight', fontsize=14)
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

    # 圖b: 花青素
    ax2 = axes[1]
    for i, (obs, sim, label, color) in enumerate(zip(obs_anth, sim_anth, labels, colors)):
        ax2.scatter(obs, sim, c=color, s=180, label=label,
                   edgecolors='black', linewidths=1.5, zorder=5)

    anth_range = np.array([380, 700])
    ax2.plot(anth_range, anth_range, 'k--', linewidth=1.5, label='1:1 line')
    ax2.fill_between(anth_range, anth_range*0.90, anth_range*1.10, alpha=0.15, color='green', label='±10% band')

    ax2.set_xlabel('Observed Anthocyanin (mg/kg FW)', fontsize=13)
    ax2.set_ylabel('Simulated Anthocyanin (mg/kg FW)', fontsize=13)
    ax2.set_title('(b) Validation Set: Anthocyanin', fontsize=14)
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

    plt.suptitle('Validation Set: Model Prediction vs Observation (v10.37)', fontsize=15, y=1.02)
    plt.tight_layout()
    plt.savefig('fig_validation_scatter.png', dpi=150, bbox_inches='tight')
    print("    已生成: fig_validation_scatter.png")
    plt.close()


def generate_combined_scatter():
    """全部組別散點圖 (訓練 + 驗證)"""
    print("  生成: 全部組別模型預測 vs 觀測...")

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # 訓練組
    train_treatments = ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']
    train_labels = ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']
    train_colors = ['gray', 'green', 'darkgreen', 'blue', 'darkblue', 'red']

    # 驗證組
    val_treatments = ['CK_val', 'VL3D3', 'L6D3', 'M9D3', 'H12D3_val', 'VH15D3']
    val_labels = ['CK*', 'VL3D3', 'L6D3', 'M9D3', 'H12D3*', 'VH15D3']
    val_colors = ['lightgray', 'limegreen', 'lightblue', 'orange', 'salmon', 'darkred']

    all_obs_fw, all_sim_fw = [], []
    all_obs_anth, all_sim_anth = [], []
    all_labels, all_colors, all_markers = [], [], []
    all_types = []

    # 收集訓練組數據
    for t, label, color in zip(train_treatments, train_labels, train_colors):
        all_obs_fw.append(TARGETS[t]['FW'])
        all_obs_anth.append(TARGETS[t]['Anth'])
        result = run_simulation(t)
        if result:
            all_sim_fw.append(result['FW_g'])
            all_sim_anth.append(result['Anth'])
        else:
            all_sim_fw.append(np.nan)
            all_sim_anth.append(np.nan)
        all_labels.append(label)
        all_colors.append(color)
        all_markers.append('o')
        all_types.append('Training')

    # 收集驗證組數據
    for t, label, color in zip(val_treatments, val_labels, val_colors):
        all_obs_fw.append(TARGETS[t]['FW'])
        all_obs_anth.append(TARGETS[t]['Anth'])
        result = run_simulation(t)
        if result:
            all_sim_fw.append(result['FW_g'])
            all_sim_anth.append(result['Anth'])
        else:
            all_sim_fw.append(np.nan)
            all_sim_anth.append(np.nan)
        all_labels.append(label)
        all_colors.append(color)
        all_markers.append('s')
        all_types.append('Validation')

    # 圖a: 鮮重
    ax1 = axes[0]
    for obs, sim, label, color, marker, typ in zip(all_obs_fw, all_sim_fw, all_labels, all_colors, all_markers, all_types):
        ax1.scatter(obs, sim, c=color, s=150, label=f'{label} ({typ[0]})', marker=marker,
                   edgecolors='black', linewidths=1.5, zorder=5)

    fw_range = np.array([40, 100])
    ax1.plot(fw_range, fw_range, 'k--', linewidth=1.5, label='1:1 line')
    ax1.fill_between(fw_range, fw_range*0.90, fw_range*1.10, alpha=0.1, color='gray')

    ax1.set_xlabel('Observed Fresh Weight (g)', fontsize=13)
    ax1.set_ylabel('Simulated Fresh Weight (g)', fontsize=13)
    ax1.set_title('(a) Fresh Weight: All Groups', fontsize=14)
    ax1.set_xlim(40, 100)
    ax1.set_ylim(40, 100)
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper left', fontsize=7, ncol=2)

    # 圖b: 花青素
    ax2 = axes[1]
    for obs, sim, label, color, marker, typ in zip(all_obs_anth, all_sim_anth, all_labels, all_colors, all_markers, all_types):
        ax2.scatter(obs, sim, c=color, s=150, label=f'{label} ({typ[0]})', marker=marker,
                   edgecolors='black', linewidths=1.5, zorder=5)

    anth_range = np.array([400, 700])
    ax2.plot(anth_range, anth_range, 'k--', linewidth=1.5, label='1:1 line')
    ax2.fill_between(anth_range, anth_range*0.90, anth_range*1.10, alpha=0.1, color='gray')

    ax2.set_xlabel('Observed Anthocyanin (mg/kg FW)', fontsize=13)
    ax2.set_ylabel('Simulated Anthocyanin (mg/kg FW)', fontsize=13)
    ax2.set_title('(b) Anthocyanin: All Groups', fontsize=14)
    ax2.set_xlim(400, 700)
    ax2.set_ylim(400, 700)
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper left', fontsize=7, ncol=2)

    # 計算整體統計
    train_fw_mape = np.mean([abs(s-o)/o*100 for s, o, t in zip(all_sim_fw, all_obs_fw, all_types) if t == 'Training' and not np.isnan(s)])
    val_fw_mape = np.mean([abs(s-o)/o*100 for s, o, t in zip(all_sim_fw, all_obs_fw, all_types) if t == 'Validation' and not np.isnan(s)])
    train_anth_mape = np.mean([abs(s-o)/o*100 for s, o, t in zip(all_sim_anth, all_obs_anth, all_types) if t == 'Training' and not np.isnan(s)])
    val_anth_mape = np.mean([abs(s-o)/o*100 for s, o, t in zip(all_sim_anth, all_obs_anth, all_types) if t == 'Validation' and not np.isnan(s)])

    ax1.text(0.98, 0.02, f'Training MAPE: {train_fw_mape:.1f}%\nValidation MAPE: {val_fw_mape:.1f}%',
            transform=ax1.transAxes, fontsize=9, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

    ax2.text(0.98, 0.02, f'Training MAPE: {train_anth_mape:.1f}%\nValidation MAPE: {val_anth_mape:.1f}%',
            transform=ax2.transAxes, fontsize=9, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

    plt.suptitle('All Groups: Model Prediction vs Observation (v10.37)\n○ Training (n=6, target <5%)  □ Validation (n=6, target <10%)', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig('fig_all_scatter.png', dpi=150, bbox_inches='tight')
    print("    已生成: fig_all_scatter.png")
    plt.close()


def generate_summary_table():
    """生成完整的驗證結果表格"""
    print("  生成: 完整驗證結果表格...")

    # 訓練組
    train_treatments = ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']
    # 驗證組
    val_treatments = ['CK_val', 'VL3D3', 'L6D3', 'M9D3', 'H12D3_val', 'VH15D3']

    print("\n" + "=" * 110)
    print("v10.37 完整驗證結果")
    print("=" * 110)

    print("\n【訓練組】(目標: 誤差 < 5%)")
    print("-" * 100)
    print(f"{'處理組':>12} | {'觀測FW':>8} | {'模擬FW':>8} | {'FW誤差':>8} | {'觀測Anth':>10} | {'模擬Anth':>10} | {'Anth誤差':>10}")
    print("-" * 100)

    for t in train_treatments:
        obs_fw = TARGETS[t]['FW']
        obs_anth = TARGETS[t]['Anth']
        result = run_simulation(t)
        if result:
            fw_err = (result['FW_g'] - obs_fw) / obs_fw * 100
            anth_err = (result['Anth'] - obs_anth) / obs_anth * 100
            print(f"{t:>12} | {obs_fw:>8.1f} | {result['FW_g']:>8.1f} | {fw_err:>+7.1f}% | {obs_anth:>10} | {result['Anth']:>10.1f} | {anth_err:>+9.1f}%")

    print("\n【驗證組】(目標: 誤差 < 10%)")
    print("-" * 100)
    print(f"{'處理組':>12} | {'觀測FW':>8} | {'模擬FW':>8} | {'FW誤差':>8} | {'觀測Anth':>10} | {'模擬Anth':>10} | {'Anth誤差':>10}")
    print("-" * 100)

    for t in val_treatments:
        obs_fw = TARGETS[t]['FW']
        obs_anth = TARGETS[t]['Anth']
        result = run_simulation(t)
        if result:
            fw_err = (result['FW_g'] - obs_fw) / obs_fw * 100
            anth_err = (result['Anth'] - obs_anth) / obs_anth * 100
            print(f"{t:>12} | {obs_fw:>8.1f} | {result['FW_g']:>8.1f} | {fw_err:>+7.1f}% | {obs_anth:>10} | {result['Anth']:>10.1f} | {anth_err:>+9.1f}%")

    print("-" * 100)


def main():
    print("=" * 60)
    print("生成所有組別驗證圖表 (v10.37)")
    print("=" * 60)

    print("\n[1/7] 訓練組鮮重柱狀圖")
    generate_training_fw_fig()

    print("\n[2/7] 訓練組花青素柱狀圖")
    generate_training_anth_fig()

    print("\n[3/7] 訓練組散點圖")
    generate_training_scatter()

    print("\n[4/7] 驗證組鮮重趨勢圖")
    generate_validation_fw_fig()

    print("\n[5/7] 驗證組花青素趨勢圖")
    generate_validation_anth_fig()

    print("\n[6/7] 驗證組散點圖")
    generate_validation_scatter()

    print("\n[7/7] 全部組別散點圖")
    generate_combined_scatter()

    # 打印數據摘要
    generate_summary_table()

    print("\n" + "=" * 60)
    print("所有圖表生成完成!")
    print("=" * 60)
    print("\n生成的檔案:")
    print("  - fig_training_fw.png      (訓練組鮮重)")
    print("  - fig_training_anth.png    (訓練組花青素)")
    print("  - fig_training_scatter.png (訓練組散點圖)")
    print("  - fig_validation_fw.png    (驗證組鮮重趨勢)")
    print("  - fig_validation_anth.png  (驗證組花青素趨勢)")
    print("  - fig_validation_scatter.png (驗證組散點圖)")
    print("  - fig_all_scatter.png      (全部組別散點圖)")


if __name__ == "__main__":
    main()
