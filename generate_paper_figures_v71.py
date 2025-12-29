"""
論文圖表生成腳本 (v7.1)
======================

生成以下圖表：
1. FW 模擬 vs 實驗比較圖
2. Anth 模擬 vs 實驗比較圖
3. 日內逆境因子曲線 (day_factor vs hours)
4. n_day 敏感性分析圖
5. 模型達標率 vs n_day 曲線
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment
from simulate_uva_model_v7 import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio, ALL_PARAMS

# 設定中文字體
plt.rcParams['font.family'] = ['DejaVu Sans', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False

# 處理組順序和顏色
TREATMENTS = ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']
COLORS = ['#2ecc71', '#3498db', '#9b59b6', '#f39c12', '#e74c3c', '#1abc9c']


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


def get_all_results(params):
    """獲取所有處理組的模擬結果"""
    results = {}
    for t in TREATMENTS:
        FW, Anth, Stress = simulate_treatment(t, params)
        if FW is not None:
            target = TARGETS[t]
            results[t] = {
                'FW_sim': FW, 'Anth_sim': Anth, 'Stress': Stress,
                'FW_obs': target['FW'], 'Anth_obs': target['Anth'],
                'FW_err': (FW - target['FW']) / target['FW'] * 100,
                'Anth_err': (Anth - target['Anth']) / target['Anth'] * 100
            }
    return results


def plot_fw_comparison(results, save_path='paper_figures_v71/fig1_fw_comparison.png'):
    """圖1: FW 模擬 vs 實驗比較"""
    fig, ax = plt.subplots(figsize=(10, 6))

    x = np.arange(len(TREATMENTS))
    width = 0.35

    fw_obs = [results[t]['FW_obs'] for t in TREATMENTS]
    fw_sim = [results[t]['FW_sim'] for t in TREATMENTS]
    fw_err = [results[t]['FW_err'] for t in TREATMENTS]

    bars1 = ax.bar(x - width/2, fw_obs, width, label='Observed', color='#3498db', alpha=0.8)
    bars2 = ax.bar(x + width/2, fw_sim, width, label='Simulated', color='#e74c3c', alpha=0.8)

    # 添加誤差標籤
    for i, (bar, err) in enumerate(zip(bars2, fw_err)):
        ax.annotate(f'{err:+.1f}%',
                    xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                    xytext=(0, 3), textcoords='offset points',
                    ha='center', va='bottom', fontsize=9)

    ax.set_xlabel('Treatment', fontsize=12)
    ax.set_ylabel('Fresh Weight (g/plant)', fontsize=12)
    ax.set_title('Fresh Weight: Observed vs Simulated (v7.1)', fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(TREATMENTS)
    ax.legend()
    ax.set_ylim(0, 110)

    # 添加 5% 誤差線
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")


def plot_anth_comparison(results, save_path='paper_figures_v71/fig2_anth_comparison.png'):
    """圖2: Anth 模擬 vs 實驗比較"""
    fig, ax = plt.subplots(figsize=(10, 6))

    x = np.arange(len(TREATMENTS))
    width = 0.35

    anth_obs = [results[t]['Anth_obs'] for t in TREATMENTS]
    anth_sim = [results[t]['Anth_sim'] for t in TREATMENTS]
    anth_err = [results[t]['Anth_err'] for t in TREATMENTS]

    bars1 = ax.bar(x - width/2, anth_obs, width, label='Observed', color='#9b59b6', alpha=0.8)
    bars2 = ax.bar(x + width/2, anth_sim, width, label='Simulated', color='#f39c12', alpha=0.8)

    # 添加誤差標籤
    for i, (bar, err) in enumerate(zip(bars2, anth_err)):
        ax.annotate(f'{err:+.1f}%',
                    xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                    xytext=(0, 3), textcoords='offset points',
                    ha='center', va='bottom', fontsize=9)

    ax.set_xlabel('Treatment', fontsize=12)
    ax.set_ylabel('Anthocyanin (mg/kg FW)', fontsize=12)
    ax.set_title('Anthocyanin: Observed vs Simulated (v7.1)', fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(TREATMENTS)
    ax.legend()
    ax.set_ylim(0, 80)

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")


def plot_day_factor_curve(params, save_path='paper_figures_v71/fig3_day_factor.png'):
    """圖3: 日內逆境因子曲線"""
    fig, ax = plt.subplots(figsize=(10, 6))

    hours = np.linspace(0, 12, 100)
    k_day = params['k_day']
    n_day = params['n_day']

    factor = 1 + k_day * (hours ** n_day)

    ax.plot(hours, factor, 'b-', linewidth=2, label=f'n = {n_day:.1f} (data-driven)')

    # 標記實驗點
    exp_hours = [3, 6, 12]
    exp_factors = [1 + k_day * (h ** n_day) for h in exp_hours]
    ax.scatter(exp_hours, exp_factors, c='red', s=100, zorder=5, label='Experimental timepoints')

    # 添加標籤
    for h, f in zip(exp_hours, exp_factors):
        ax.annotate(f'{h}h: {f:.1f}', xy=(h, f), xytext=(5, 10),
                    textcoords='offset points', fontsize=10)

    ax.set_xlabel('UVA Exposure Duration (hours)', fontsize=12)
    ax.set_ylabel('Day Factor (dimensionless)', fontsize=12)
    ax.set_title(f'Intraday Stress Factor: 1 + k × hours^n\n(k = {k_day:.2e}, n = {n_day:.1f})', fontsize=14)
    ax.legend()
    ax.set_xlim(0, 12)
    ax.set_ylim(0, max(factor) * 1.1)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")


def plot_n_day_sensitivity(save_path='paper_figures_v71/fig4_n_day_sensitivity.png'):
    """圖4: n_day 敏感性分析"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    n_values = np.linspace(2, 12, 21)
    scores = []
    fw_ok_list = []
    anth_ok_list = []

    factor_6h_target = 3.8

    for n_day in n_values:
        k_day = (factor_6h_target - 1) / (6 ** n_day)
        params = ALL_PARAMS.copy()
        params['n_day'] = n_day
        params['k_day'] = k_day

        results = get_all_results(params)
        if results:
            fw_ok = sum(1 for t in results if abs(results[t]['FW_err']) < 5)
            anth_ok = sum(1 for t in results if abs(results[t]['Anth_err']) < 10)
            scores.append(fw_ok + anth_ok)
            fw_ok_list.append(fw_ok)
            anth_ok_list.append(anth_ok)
        else:
            scores.append(0)
            fw_ok_list.append(0)
            anth_ok_list.append(0)

    # 左圖：達標率
    ax1 = axes[0]
    ax1.plot(n_values, scores, 'b-o', linewidth=2, markersize=6, label='Total (FW + Anth)')
    ax1.plot(n_values, fw_ok_list, 'g--s', linewidth=1.5, markersize=5, label='FW (<5%)')
    ax1.plot(n_values, anth_ok_list, 'r--^', linewidth=1.5, markersize=5, label='Anth (<10%)')

    # 標記最佳區間
    ax1.axvspan(6.5, 7.0, alpha=0.2, color='green', label='Optimal range')
    ax1.axvline(x=6.6, color='red', linestyle='--', alpha=0.7, label=f'Fitted n = 6.6')

    ax1.set_xlabel('n_day (power exponent)', fontsize=12)
    ax1.set_ylabel('Number of targets met', fontsize=12)
    ax1.set_title('Model Performance vs n_day', fontsize=14)
    ax1.legend(loc='lower right')
    ax1.set_xlim(2, 12)
    ax1.set_ylim(0, 13)
    ax1.grid(True, alpha=0.3)

    # 右圖：12h/6h 比值
    ax2 = axes[1]
    ratios = []
    for n_day in n_values:
        k_day = (factor_6h_target - 1) / (6 ** n_day)
        factor_12h = 1 + k_day * (12 ** n_day)
        ratios.append(factor_12h / factor_6h_target)

    ax2.semilogy(n_values, ratios, 'b-o', linewidth=2, markersize=6)
    ax2.axhline(y=73.8, color='red', linestyle='--', alpha=0.7, label='Fitted ratio = 73.8')
    ax2.axvline(x=6.6, color='red', linestyle='--', alpha=0.7, label='Fitted n = 6.6')

    ax2.set_xlabel('n_day (power exponent)', fontsize=12)
    ax2.set_ylabel('12h/6h Factor Ratio (log scale)', fontsize=12)
    ax2.set_title('12h/6h Ratio vs n_day', fontsize=14)
    ax2.legend()
    ax2.set_xlim(2, 12)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")


def plot_scatter_comparison(results, save_path='paper_figures_v71/fig5_scatter_comparison.png'):
    """圖5: 模擬 vs 觀測散點圖"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # FW 散點圖
    ax1 = axes[0]
    fw_obs = [results[t]['FW_obs'] for t in TREATMENTS]
    fw_sim = [results[t]['FW_sim'] for t in TREATMENTS]

    ax1.scatter(fw_obs, fw_sim, c=COLORS, s=100, zorder=5)
    for i, t in enumerate(TREATMENTS):
        ax1.annotate(t, (fw_obs[i], fw_sim[i]), xytext=(5, 5), textcoords='offset points')

    # 1:1 線
    lims = [min(min(fw_obs), min(fw_sim)) * 0.9, max(max(fw_obs), max(fw_sim)) * 1.1]
    ax1.plot(lims, lims, 'k--', alpha=0.5, label='1:1 line')

    # 計算 R²
    from scipy import stats
    slope, intercept, r_value, p_value, std_err = stats.linregress(fw_obs, fw_sim)
    ax1.plot(fw_obs, [slope * x + intercept for x in fw_obs], 'r-', alpha=0.7,
             label=f'R² = {r_value**2:.3f}')

    ax1.set_xlabel('Observed FW (g)', fontsize=12)
    ax1.set_ylabel('Simulated FW (g)', fontsize=12)
    ax1.set_title('Fresh Weight Validation', fontsize=14)
    ax1.legend()
    ax1.set_xlim(lims)
    ax1.set_ylim(lims)
    ax1.set_aspect('equal')

    # Anth 散點圖
    ax2 = axes[1]
    anth_obs = [results[t]['Anth_obs'] for t in TREATMENTS]
    anth_sim = [results[t]['Anth_sim'] for t in TREATMENTS]

    ax2.scatter(anth_obs, anth_sim, c=COLORS, s=100, zorder=5)
    for i, t in enumerate(TREATMENTS):
        ax2.annotate(t, (anth_obs[i], anth_sim[i]), xytext=(5, 5), textcoords='offset points')

    lims = [min(min(anth_obs), min(anth_sim)) * 0.9, max(max(anth_obs), max(anth_sim)) * 1.1]
    ax2.plot(lims, lims, 'k--', alpha=0.5, label='1:1 line')

    slope, intercept, r_value, p_value, std_err = stats.linregress(anth_obs, anth_sim)
    ax2.plot(anth_obs, [slope * x + intercept for x in anth_obs], 'r-', alpha=0.7,
             label=f'R² = {r_value**2:.3f}')

    ax2.set_xlabel('Observed Anthocyanin (mg/kg)', fontsize=12)
    ax2.set_ylabel('Simulated Anthocyanin (mg/kg)', fontsize=12)
    ax2.set_title('Anthocyanin Validation', fontsize=14)
    ax2.legend()
    ax2.set_xlim(lims)
    ax2.set_ylim(lims)
    ax2.set_aspect('equal')

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")


def generate_sensitivity_table(save_path='paper_figures_v71/sensitivity_table.csv'):
    """生成敏感性分析表格"""
    import csv

    n_values = np.linspace(2, 12, 21)
    factor_6h_target = 3.8

    with open(save_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['n_day', 'k_day', 'factor_6h', 'factor_12h', '12h/6h_ratio',
                        'FW_ok', 'Anth_ok', 'Total_ok'])

        for n_day in n_values:
            k_day = (factor_6h_target - 1) / (6 ** n_day)
            params = ALL_PARAMS.copy()
            params['n_day'] = n_day
            params['k_day'] = k_day

            factor_6h = 1 + k_day * (6 ** n_day)
            factor_12h = 1 + k_day * (12 ** n_day)
            ratio = factor_12h / factor_6h

            results = get_all_results(params)
            if results:
                fw_ok = sum(1 for t in results if abs(results[t]['FW_err']) < 5)
                anth_ok = sum(1 for t in results if abs(results[t]['Anth_err']) < 10)
                total = fw_ok + anth_ok
            else:
                fw_ok = anth_ok = total = 0

            writer.writerow([f'{n_day:.2f}', f'{k_day:.2e}', f'{factor_6h:.1f}',
                           f'{factor_12h:.1f}', f'{ratio:.1f}', fw_ok, anth_ok, total])

    print(f"Saved: {save_path}")


if __name__ == "__main__":
    import os

    # 創建輸出目錄
    os.makedirs('paper_figures_v71', exist_ok=True)

    print("=" * 60)
    print("論文圖表生成 (v7.1)")
    print("=" * 60)

    # 獲取模擬結果
    print("\n1. Running simulations...")
    results = get_all_results(ALL_PARAMS)

    # 打印結果摘要
    print("\nSimulation Results:")
    print("-" * 60)
    for t in TREATMENTS:
        r = results[t]
        fw_ok = "ok" if abs(r['FW_err']) < 5 else "NG"
        anth_ok = "ok" if abs(r['Anth_err']) < 10 else "NG"
        print(f"{t:<8} FW: {r['FW_sim']:.1f}g ({r['FW_err']:+.1f}% {fw_ok}) "
              f"Anth: {r['Anth_sim']:.1f} ({r['Anth_err']:+.1f}% {anth_ok})")

    # 生成圖表
    print("\n2. Generating figures...")
    plot_fw_comparison(results)
    plot_anth_comparison(results)
    plot_day_factor_curve(ALL_PARAMS)
    plot_n_day_sensitivity()
    plot_scatter_comparison(results)

    # 生成敏感性分析表格
    print("\n3. Generating sensitivity table...")
    generate_sensitivity_table()

    print("\n" + "=" * 60)
    print("All figures saved to paper_figures_v71/")
    print("=" * 60)
