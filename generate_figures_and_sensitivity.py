"""
================================================================================
UVA 萵苣模型 v10.34 圖表生成與敏感性分析
================================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.integrate import solve_ivp
import copy

# 設定字體
rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial Unicode MS', 'SimHei']
rcParams['axes.unicode_minus'] = False

# 導入模型 (更新為 v10)
from simulate_uva_model_v10 import (
    ALL_PARAMS, UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
)
from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment


def run_simulation(params_dict, treatment_name, return_trajectory=False):
    """運行單次模擬"""
    # 建立參數
    p = UVAParams(params_dict)
    env = get_env_for_treatment(treatment_name)

    # 初始條件
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    C_buf_init = Xd_init * 0.1
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6
    initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]  # v9: 6 states

    # 時間設定
    transplant_day = SIMULATION['transplant_offset']
    simulation_days = SIMULATION['days']
    t_start = transplant_day * 86400
    t_end = (transplant_day + simulation_days) * 86400

    # 設定評估時間點
    if return_trajectory:
        n_points = simulation_days * 24 + 1
        t_eval = np.linspace(t_start, t_end, n_points)
    else:
        t_eval = np.array([t_end])

    # 運行 ODE
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

    # 提取終點結果
    Xd_f, Cbuf_f, LAI_f, Anth_f, Stress_f, ROS_f = sol.y[:, -1]  # v9: 6 states
    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress_f, p)
    FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
    FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
    Anth_sim = Anth_f / FW_total_kg * 1e6  # ppm

    result = {
        'FW_g': FW_sim,
        'Anth': Anth_sim,
        'LAI': LAI_f,
        'Stress': Stress_f
    }

    if return_trajectory:
        result['time'] = (sol.t - t_start) / 86400  # 轉換為天數
        result['X_d'] = sol.y[0]
        result['C_buf'] = sol.y[1]
        result['LAI_t'] = sol.y[2]
        result['Anth_t'] = sol.y[3]
        result['Stress_t'] = sol.y[4]
        result['ROS_t'] = sol.y[5]  # v9: add ROS

    return result


def generate_comparison_figure():
    """生成模擬與觀測比較圖"""
    print("生成模擬與觀測比較圖...")

    treatments = ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']

    # 運行所有模擬
    results = {}
    for t in treatments:
        results[t] = run_simulation(ALL_PARAMS, t)

    # 準備數據
    sim_fw = [results[t]['FW_g'] for t in treatments]
    obs_fw = [TARGETS[t]['FW'] for t in treatments]
    sim_anth = [results[t]['Anth'] for t in treatments]
    obs_anth = [TARGETS[t]['Anth'] for t in treatments]

    # 繪圖
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    x = np.arange(len(treatments))
    width = 0.35

    # 鮮重比較
    ax1 = axes[0]
    bars1 = ax1.bar(x - width/2, obs_fw, width, label='Observed', color='steelblue', alpha=0.8)
    bars2 = ax1.bar(x + width/2, sim_fw, width, label='Simulated', color='coral', alpha=0.8)
    ax1.set_ylabel('Fresh Weight (g/plant)', fontsize=12)
    ax1.set_title('Fresh Weight: Observed vs Simulated', fontsize=14)
    ax1.set_xticks(x)
    ax1.set_xticklabels(treatments, fontsize=11)
    ax1.legend()
    ax1.set_ylim(0, 110)

    # 添加誤差標註
    for i, (obs, sim) in enumerate(zip(obs_fw, sim_fw)):
        error = (sim - obs) / obs * 100
        color = 'green' if abs(error) < 5 else 'red'
        ax1.annotate(f'{error:+.1f}%', xy=(x[i] + width/2, sim + 2),
                    ha='center', fontsize=9, color=color)

    # 花青素比較
    ax2 = axes[1]
    bars3 = ax2.bar(x - width/2, obs_anth, width, label='Observed', color='purple', alpha=0.8)
    bars4 = ax2.bar(x + width/2, sim_anth, width, label='Simulated', color='orchid', alpha=0.8)
    ax2.set_ylabel('Anthocyanin (mg/100g FW)', fontsize=12)
    ax2.set_title('Anthocyanin: Observed vs Simulated', fontsize=14)
    ax2.set_xticks(x)
    ax2.set_xticklabels(treatments, fontsize=11)
    ax2.legend()
    ax2.set_ylim(0, 80)

    # 添加誤差標註
    for i, (obs, sim) in enumerate(zip(obs_anth, sim_anth)):
        error = (sim - obs) / obs * 100
        color = 'green' if abs(error) < 10 else 'red'
        ax2.annotate(f'{error:+.1f}%', xy=(x[i] + width/2, sim + 2),
                    ha='center', fontsize=9, color=color)

    plt.tight_layout()
    plt.savefig('fig1_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> fig1_comparison.png")


def generate_dynamics_figure():
    """生成動態變化圖"""
    print("生成動態變化圖...")

    treatments = ['CK', 'L6D6', 'L6D6-N', 'H12D3']
    colors = {'CK': 'gray', 'L6D6': 'blue', 'L6D6-N': 'orange', 'H12D3': 'red'}

    results = {}
    for t in treatments:
        results[t] = run_simulation(ALL_PARAMS, t, return_trajectory=True)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 乾物質動態
    ax1 = axes[0, 0]
    for t in treatments:
        ax1.plot(results[t]['time'], results[t]['X_d'] * 1000 / ENV_BASE['plant_density'],
                label=t, color=colors[t], linewidth=2)
    ax1.set_xlabel('Days after transplant', fontsize=12)
    ax1.set_ylabel('Dry Weight (g/plant)', fontsize=12)
    ax1.set_title('Dry Matter Accumulation', fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # LAI 動態
    ax2 = axes[0, 1]
    for t in treatments:
        ax2.plot(results[t]['time'], results[t]['LAI_t'],
                label=t, color=colors[t], linewidth=2)
    ax2.set_xlabel('Days after transplant', fontsize=12)
    ax2.set_ylabel('LAI (m²/m²)', fontsize=12)
    ax2.set_title('Leaf Area Index', fontsize=14)
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Stress 動態
    ax3 = axes[1, 0]
    for t in treatments:
        ax3.plot(results[t]['time'], results[t]['Stress_t'],
                label=t, color=colors[t], linewidth=2)
    ax3.set_xlabel('Days after transplant', fontsize=12)
    ax3.set_ylabel('Stress Index', fontsize=12)
    ax3.set_title('Cumulative Stress', fontsize=14)
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # 花青素動態 (近似濃度)
    ax4 = axes[1, 1]
    for t in treatments:
        X_d = results[t]['X_d']
        Anth = results[t]['Anth_t']
        # 簡化：用固定 DW:FW 比例
        FW_kg_m2 = X_d / 0.05
        Anth_ppm = np.where(FW_kg_m2 > 1e-6, Anth / FW_kg_m2 * 1e6, 0)
        ax4.plot(results[t]['time'], Anth_ppm,
                label=t, color=colors[t], linewidth=2)
    ax4.set_xlabel('Days after transplant', fontsize=12)
    ax4.set_ylabel('Anthocyanin (ppm)', fontsize=12)
    ax4.set_title('Anthocyanin Content', fontsize=14)
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('fig2_dynamics.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> fig2_dynamics.png")


def sensitivity_analysis():
    """執行敏感性分析"""
    print("\n執行敏感性分析...")

    # 定義要分析的參數 (v9 參數名)
    params_to_analyze = [
        ('stress_damage_coeff', 'k_damage'),
        ('A_vulnerability', 'A_vuln'),
        ('k_hours_nonlinear', 'k_hours'),
        ('k_repair_lai', 'k_repair'),
        ('k_circadian', 'k_circ'),
        ('V_max_anth', 'V_max_anth'),
        ('K_stress_anth', 'K_anth'),
        ('alpha_anth_protection', 'α_protect'),
        ('stress_photosynthesis_inhibition', 'β_photo'),
    ]

    # 變化範圍：-30% 到 +30%
    variations = np.linspace(-0.3, 0.3, 7)

    # 目標處理組
    target_treatments = ['L6D6', 'L6D6-N', 'H12D3']

    results = {}

    for param_name, param_label in params_to_analyze:
        print(f"  {param_label}...", end="", flush=True)
        results[param_name] = {
            'label': param_label,
            'variations': variations,
            'FW': {t: [] for t in target_treatments},
            'Anth': {t: [] for t in target_treatments}
        }

        base_value = ALL_PARAMS[param_name]

        for var in variations:
            params = copy.deepcopy(ALL_PARAMS)
            params[param_name] = base_value * (1 + var)

            for t in target_treatments:
                try:
                    result = run_simulation(params, t)
                    if result:
                        results[param_name]['FW'][t].append(result['FW_g'])
                        results[param_name]['Anth'][t].append(result['Anth'])
                    else:
                        results[param_name]['FW'][t].append(np.nan)
                        results[param_name]['Anth'][t].append(np.nan)
                except Exception as e:
                    results[param_name]['FW'][t].append(np.nan)
                    results[param_name]['Anth'][t].append(np.nan)
        print(" done")

    return results, params_to_analyze, target_treatments, variations


def generate_sensitivity_figure(results, params_to_analyze, target_treatments, variations):
    """生成敏感性分析圖"""
    print("\n生成敏感性分析圖...")

    n_params = len(params_to_analyze)
    fig, axes = plt.subplots(n_params, 2, figsize=(14, 3 * n_params))

    colors = {'L6D6': 'blue', 'L6D6-N': 'orange', 'H12D3': 'red'}

    for i, (param_name, param_label) in enumerate(params_to_analyze):
        base_fw = {t: results[param_name]['FW'][t][3] for t in target_treatments}
        base_anth = {t: results[param_name]['Anth'][t][3] for t in target_treatments}

        # FW 敏感性
        ax1 = axes[i, 0]
        for t in target_treatments:
            if base_fw[t] and base_fw[t] > 0:
                fw_change = [(v - base_fw[t]) / base_fw[t] * 100
                            for v in results[param_name]['FW'][t]]
                ax1.plot(variations * 100, fw_change, 'o-',
                        label=t, color=colors[t], linewidth=2)
        ax1.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        ax1.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
        ax1.set_xlabel('Parameter Change (%)', fontsize=10)
        ax1.set_ylabel('FW Change (%)', fontsize=10)
        ax1.set_title(f'{param_label} → FW', fontsize=11)
        ax1.legend(fontsize=8)
        ax1.grid(True, alpha=0.3)
        ax1.set_xlim(-35, 35)

        # Anth 敏感性
        ax2 = axes[i, 1]
        for t in target_treatments:
            if base_anth[t] and base_anth[t] > 0:
                anth_change = [(v - base_anth[t]) / base_anth[t] * 100
                              for v in results[param_name]['Anth'][t]]
                ax2.plot(variations * 100, anth_change, 'o-',
                        label=t, color=colors[t], linewidth=2)
        ax2.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
        ax2.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
        ax2.set_xlabel('Parameter Change (%)', fontsize=10)
        ax2.set_ylabel('Anth Change (%)', fontsize=10)
        ax2.set_title(f'{param_label} → Anthocyanin', fontsize=11)
        ax2.legend(fontsize=8)
        ax2.grid(True, alpha=0.3)
        ax2.set_xlim(-35, 35)

    plt.tight_layout()
    plt.savefig('fig3_sensitivity.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> fig3_sensitivity.png")


def calculate_sensitivity_indices(results, params_to_analyze, target_treatments):
    """計算敏感性指數"""
    sensitivity = {}

    for param_name, param_label in params_to_analyze:
        sensitivity[param_name] = {'label': param_label, 'FW': {}, 'Anth': {}}

        for t in target_treatments:
            fw_low = results[param_name]['FW'][t][0]
            fw_high = results[param_name]['FW'][t][-1]
            fw_base = results[param_name]['FW'][t][3]

            anth_low = results[param_name]['Anth'][t][0]
            anth_high = results[param_name]['Anth'][t][-1]
            anth_base = results[param_name]['Anth'][t][3]

            # 敏感性指數 = (Output Change %) / (Input Change %)
            if fw_base and fw_base > 0:
                fw_sens = ((fw_high - fw_low) / fw_base * 100) / 60
            else:
                fw_sens = 0

            if anth_base and anth_base > 0:
                anth_sens = ((anth_high - anth_low) / anth_base * 100) / 60
            else:
                anth_sens = 0

            sensitivity[param_name]['FW'][t] = fw_sens
            sensitivity[param_name]['Anth'][t] = anth_sens

    return sensitivity


def print_sensitivity_table(sensitivity, params_to_analyze, target_treatments):
    """打印敏感性分析表格"""
    print("\n" + "=" * 75)
    print("敏感性分析結果")
    print("=" * 75)
    print("\n敏感性指數 S = (輸出變化%) / (參數變化%)")
    print("  |S| > 1: 高敏感 (輸出變化 > 輸入變化)")
    print("  |S| < 0.5: 低敏感")
    print()

    # FW 敏感性表
    print("-" * 75)
    print("鮮重 (FW) 敏感性指數")
    print("-" * 75)
    header = f"{'參數':<20} | " + " | ".join([f"{t:>8}" for t in target_treatments]) + " |  平均"
    print(header)
    print("-" * 75)

    for param_name, param_label in params_to_analyze:
        values = [sensitivity[param_name]['FW'][t] for t in target_treatments]
        avg = np.mean(values)
        row = f"{param_label:<20} | " + " | ".join([f"{v:>8.3f}" for v in values]) + f" | {avg:>6.3f}"
        print(row)

    # Anth 敏感性表
    print("\n" + "-" * 75)
    print("花青素 (Anth) 敏感性指數")
    print("-" * 75)
    print(header)
    print("-" * 75)

    for param_name, param_label in params_to_analyze:
        values = [sensitivity[param_name]['Anth'][t] for t in target_treatments]
        avg = np.mean(values)
        row = f"{param_label:<20} | " + " | ".join([f"{v:>8.3f}" for v in values]) + f" | {avg:>6.3f}"
        print(row)

    print("=" * 75)


def generate_sensitivity_heatmap(sensitivity, params_to_analyze, target_treatments):
    """生成敏感性熱力圖"""
    print("\n生成敏感性熱力圖...")

    fig, axes = plt.subplots(1, 2, figsize=(12, 8))

    param_labels = [p[1] for p in params_to_analyze]

    # FW 熱力圖
    fw_data = np.array([[sensitivity[p[0]]['FW'][t] for t in target_treatments]
                        for p in params_to_analyze])

    ax1 = axes[0]
    im1 = ax1.imshow(fw_data, cmap='RdBu_r', aspect='auto', vmin=-1.5, vmax=1.5)
    ax1.set_xticks(range(len(target_treatments)))
    ax1.set_xticklabels(target_treatments, fontsize=11)
    ax1.set_yticks(range(len(param_labels)))
    ax1.set_yticklabels(param_labels, fontsize=10)
    ax1.set_title('Fresh Weight Sensitivity', fontsize=14)

    for i in range(len(param_labels)):
        for j in range(len(target_treatments)):
            text = ax1.text(j, i, f'{fw_data[i, j]:.2f}',
                           ha='center', va='center', fontsize=9,
                           color='white' if abs(fw_data[i, j]) > 0.8 else 'black')

    plt.colorbar(im1, ax=ax1, label='Sensitivity Index')

    # Anth 熱力圖
    anth_data = np.array([[sensitivity[p[0]]['Anth'][t] for t in target_treatments]
                          for p in params_to_analyze])

    ax2 = axes[1]
    im2 = ax2.imshow(anth_data, cmap='RdBu_r', aspect='auto', vmin=-1.5, vmax=1.5)
    ax2.set_xticks(range(len(target_treatments)))
    ax2.set_xticklabels(target_treatments, fontsize=11)
    ax2.set_yticks(range(len(param_labels)))
    ax2.set_yticklabels(param_labels, fontsize=10)
    ax2.set_title('Anthocyanin Sensitivity', fontsize=14)

    for i in range(len(param_labels)):
        for j in range(len(target_treatments)):
            text = ax2.text(j, i, f'{anth_data[i, j]:.2f}',
                           ha='center', va='center', fontsize=9,
                           color='white' if abs(anth_data[i, j]) > 0.8 else 'black')

    plt.colorbar(im2, ax=ax2, label='Sensitivity Index')

    plt.tight_layout()
    plt.savefig('fig4_sensitivity_heatmap.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("  -> fig4_sensitivity_heatmap.png")


def main():
    """主程式"""
    print("=" * 75)
    print("UVA 萵苣模型 v7.3 圖表生成與敏感性分析")
    print("=" * 75)

    # 1. 生成比較圖
    generate_comparison_figure()

    # 2. 生成動態圖
    generate_dynamics_figure()

    # 3. 執行敏感性分析
    results, params_to_analyze, target_treatments, variations = sensitivity_analysis()

    # 4. 生成敏感性分析圖
    generate_sensitivity_figure(results, params_to_analyze, target_treatments, variations)

    # 5. 計算敏感性指數
    sensitivity = calculate_sensitivity_indices(results, params_to_analyze, target_treatments)

    # 6. 打印敏感性表格
    print_sensitivity_table(sensitivity, params_to_analyze, target_treatments)

    # 7. 生成敏感性熱力圖
    generate_sensitivity_heatmap(sensitivity, params_to_analyze, target_treatments)

    print("\n" + "=" * 75)
    print("完成！生成的圖表：")
    print("  - fig1_comparison.png: 模擬與觀測比較")
    print("  - fig2_dynamics.png: 動態變化圖")
    print("  - fig3_sensitivity.png: 敏感性分析曲線")
    print("  - fig4_sensitivity_heatmap.png: 敏感性熱力圖")
    print("=" * 75)


if __name__ == "__main__":
    main()
