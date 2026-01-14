"""
================================================================================
UVA 照射策略優化分析
================================================================================
目標: 找出最佳的 UVA 照射時間（每日小時數）和天數組合
評估指標:
  1. 花青素含量最大化
  2. 鮮重不減少（維持或增加）
  3. 綜合評分 = 花青素提升比例 × (1 + 鮮重變化比例)

================================================================================
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib import rcParams

# 設定中文字體
rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial Unicode MS', 'sans-serif']
rcParams['axes.unicode_minus'] = False

# 導入模型
from simulate_uva_model_v10 import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio, nonlinear_damage_factor
from model_config import ENV_BASE, SIMULATION

def simulate_treatment(hours_per_day, num_days, start_day=None, night_mode=False):
    """
    模擬指定的 UVA 處理方案

    參數:
    - hours_per_day: 每天照射小時數 (1-16)
    - num_days: 照射總天數 (1-15)
    - start_day: 開始照射的日期 (從播種算起)，預設為 35 - num_days
    - night_mode: 是否為夜間照射模式

    返回:
    - FW: 鮮重 (g)
    - Anth: 花青素濃度 (μg/g FW)
    - success: 是否成功
    """
    p = UVAParams()

    # 計算開始日期（確保在 Day 35 採收）
    harvest_day = 35
    if start_day is None:
        start_day = harvest_day - num_days

    end_day = start_day + num_days - 1

    # 設定照射時段
    if night_mode:
        # 夜間照射: 從 22:00 開始
        uva_hour_on = 22
        uva_hour_off = (22 + hours_per_day) % 24
    else:
        # 日間照射: 從 10:00 開始
        uva_hour_on = 10
        uva_hour_off = 10 + hours_per_day

    # 建立環境參數
    env = ENV_BASE.copy()
    env['uva_on'] = True
    env['uva_start_day'] = start_day
    env['uva_end_day'] = end_day
    env['uva_hour_on'] = uva_hour_on
    env['uva_hour_off'] = uva_hour_off
    env['uva_intensity'] = 11.0  # W/m² (與論文一致)

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
        return None, None, False

    Xd_f, Cbuf_f, LAI_f, Anth_f, Stress_f, ROS_f = sol.y[:, -1]

    # 計算平均 Stress (與主模型一致 - 只計算 UVA 照射期間)
    uva_start = start_day * 86400
    stress_sum = 0
    stress_count = 0
    for i in range(len(sol.t)):
        if sol.t[i] >= uva_start:
            stress_sum += sol.y[4, i]
            stress_count += 1
    avg_stress = stress_sum / max(1, stress_count)

    # 計算 nonlinear_factor
    nonlin_factor = nonlinear_damage_factor(hours_per_day, p)

    # 使用 avg_stress 和 nonlin_factor 計算 dw_fw_ratio (與主模型一致)
    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(avg_stress, p, nonlin_factor)

    FW_sim = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
    FW_total_kg = FW_sim / 1000 * ENV_BASE['plant_density']
    Anth_sim = Anth_f / FW_total_kg * 1e6

    return FW_sim, Anth_sim, True


def run_optimization():
    """執行優化搜尋"""

    print("=" * 80)
    print("UVA 照射策略優化分析")
    print("=" * 80)

    # 先模擬對照組
    print("\n[1] 模擬對照組 (無 UVA)...")

    p = UVAParams()
    env_ck = ENV_BASE.copy()
    env_ck['uva_on'] = False

    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    C_buf_init = Xd_init * 0.1
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6

    initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0, 0.0]

    transplant_day = SIMULATION['transplant_offset']
    simulation_days = SIMULATION['days']
    t_start = transplant_day * 86400
    t_end = (transplant_day + simulation_days) * 86400

    sol = solve_ivp(
        uva_sun_derivatives,
        (t_start, t_end),
        initial_state,
        args=(p, env_ck),
        method='RK45',
        max_step=300
    )

    Xd_f, Cbuf_f, LAI_f, Anth_f, Stress_f, ROS_f = sol.y[:, -1]
    dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress_f, p)
    FW_ck = Xd_f / ENV_BASE['plant_density'] / dw_fw_ratio * 1000
    FW_total_kg = FW_ck / 1000 * ENV_BASE['plant_density']
    Anth_ck = Anth_f / FW_total_kg * 1e6

    print(f"   對照組: FW = {FW_ck:.1f} g, Anth = {Anth_ck:.1f} ug/g FW")

    # 搜尋參數範圍
    hours_range = range(1, 13)  # 1-12 小時/天
    days_range = range(2, 13)   # 2-12 天

    print(f"\n[2] 搜尋範圍: {min(hours_range)}-{max(hours_range)} 小時/天, {min(days_range)}-{max(days_range)} 天")
    print("    (日間照射模式，從 Day 35 往前推算開始日)")

    # 儲存結果
    results = []

    print("\n[3] 執行優化搜尋...")
    total = len(hours_range) * len(days_range)
    count = 0

    for hours in hours_range:
        for days in days_range:
            count += 1
            if count % 20 == 0:
                print(f"    進度: {count}/{total} ({100*count/total:.0f}%)")

            FW, Anth, success = simulate_treatment(hours, days)

            if success:
                # 計算花青素總量 (μg/plant)
                anth_total = Anth * FW  # 濃度 × 鮮重
                anth_total_ck = Anth_ck * FW_ck

                # 計算相對變化
                fw_change = (FW - FW_ck) / FW_ck * 100
                anth_change = (Anth - Anth_ck) / Anth_ck * 100  # 濃度變化
                anth_total_change = (anth_total - anth_total_ck) / anth_total_ck * 100  # 總量變化

                # 新評分: 花青素總量變化 (同時考慮濃度和鮮重)
                score = anth_total_change

                # 安全評分: 鮮重不減超過5%時的花青素總量變化
                score_safe = anth_total_change if fw_change >= -5 else -999

                results.append({
                    'hours': hours,
                    'days': days,
                    'FW': FW,
                    'Anth': Anth,
                    'anth_total': anth_total,
                    'fw_change': fw_change,
                    'anth_change': anth_change,
                    'anth_total_change': anth_total_change,
                    'score': score,
                    'score_safe': score_safe,
                    'total_hours': hours * days
                })

    print(f"\n[4] 搜尋完成，共 {len(results)} 個有效組合")

    # 排序結果
    results_sorted_score = sorted(results, key=lambda x: x['score'], reverse=True)
    results_sorted_safe = sorted(results, key=lambda x: x['score_safe'], reverse=True)
    results_sorted_anth = sorted(results, key=lambda x: x['anth_change'], reverse=True)

    # 顯示結果
    print("\n" + "=" * 80)
    print("優化結果")
    print("=" * 80)

    print("\n【策略一】最高花青素總量 (濃度 × 鮮重):")
    print("-" * 100)
    print(f"{'排名':<4} {'小時/天':<8} {'天數':<6} {'總時數':<8} {'FW(g)':<10} {'Anth濃度':<10} {'花青素總量':<12} {'FW變化':<10} {'總量變化':<10}")
    print("-" * 100)
    for i, r in enumerate(results_sorted_score[:10]):
        print(f"{i+1:<4} {r['hours']:<8} {r['days']:<6} {r['total_hours']:<8} "
              f"{r['FW']:<10.1f} {r['Anth']:<10.1f} {r['anth_total']:<12.0f} {r['fw_change']:>+8.1f}% {r['anth_total_change']:>+8.1f}%")

    print("\n【策略二】鮮重不減少 (FW >= -5%) 下最高花青素總量:")
    print("-" * 100)
    print(f"{'排名':<4} {'小時/天':<8} {'天數':<6} {'總時數':<8} {'FW(g)':<10} {'Anth濃度':<10} {'花青素總量':<12} {'FW變化':<10} {'總量變化':<10}")
    print("-" * 100)
    safe_results = [r for r in results_sorted_safe if r['score_safe'] > -999]
    for i, r in enumerate(safe_results[:10]):
        print(f"{i+1:<4} {r['hours']:<8} {r['days']:<6} {r['total_hours']:<8} "
              f"{r['FW']:<10.1f} {r['Anth']:<10.1f} {r['anth_total']:<12.0f} {r['fw_change']:>+8.1f}% {r['anth_total_change']:>+8.1f}%")

    print("\n【策略三】最高花青素濃度 (不考慮鮮重):")
    print("-" * 100)
    print(f"{'排名':<4} {'小時/天':<8} {'天數':<6} {'總時數':<8} {'FW(g)':<10} {'Anth濃度':<10} {'花青素總量':<12} {'FW變化':<10} {'濃度變化':<10}")
    print("-" * 100)
    for i, r in enumerate(results_sorted_anth[:10]):
        print(f"{i+1:<4} {r['hours']:<8} {r['days']:<6} {r['total_hours']:<8} "
              f"{r['FW']:<10.1f} {r['Anth']:<10.1f} {r['anth_total']:<12.0f} {r['fw_change']:>+8.1f}% {r['anth_change']:>+8.1f}%")

    # 建立熱圖
    print("\n[5] 生成熱圖...")

    # 準備熱圖數據
    hours_list = sorted(set(r['hours'] for r in results))
    days_list = sorted(set(r['days'] for r in results))

    fw_matrix = np.zeros((len(days_list), len(hours_list)))
    anth_total_matrix = np.zeros((len(days_list), len(hours_list)))

    for r in results:
        i = days_list.index(r['days'])
        j = hours_list.index(r['hours'])
        fw_matrix[i, j] = r['fw_change']
        anth_total_matrix[i, j] = r['anth_total_change']

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # 鮮重變化熱圖
    im1 = axes[0].imshow(fw_matrix, cmap='RdYlGn', aspect='auto',
                         extent=[min(hours_list)-0.5, max(hours_list)+0.5,
                                max(days_list)+0.5, min(days_list)-0.5])
    axes[0].set_xlabel('Hours per Day')
    axes[0].set_ylabel('Number of Days')
    axes[0].set_title('Fresh Weight Change (%)')
    axes[0].set_xticks(hours_list)
    axes[0].set_yticks(days_list)
    plt.colorbar(im1, ax=axes[0])

    # 花青素總量變化熱圖
    im2 = axes[1].imshow(anth_total_matrix, cmap='RdYlGn', aspect='auto',
                         extent=[min(hours_list)-0.5, max(hours_list)+0.5,
                                max(days_list)+0.5, min(days_list)-0.5])
    axes[1].set_xlabel('Hours per Day')
    axes[1].set_ylabel('Number of Days')
    axes[1].set_title('Total Anthocyanin Change (%)\n(Concentration × Fresh Weight)')
    axes[1].set_xticks(hours_list)
    axes[1].set_yticks(days_list)
    plt.colorbar(im2, ax=axes[1])

    # 標記最佳點 (花青素總量最高)
    best = results_sorted_score[0]
    for ax in axes:
        ax.plot(best['hours'], best['days'], 'w*', markersize=15, markeredgecolor='black')

    # 標記安全最佳點 (FW >= -5% 下花青素總量最高)
    best_safe = safe_results[0] if safe_results else None
    if best_safe:
        for ax in axes:
            ax.plot(best_safe['hours'], best_safe['days'], 'y*', markersize=15, markeredgecolor='black')

    plt.tight_layout()
    plt.savefig('optimization_heatmap.png', dpi=150, bbox_inches='tight')
    print("    熱圖已儲存: optimization_heatmap.png")

    # 最終建議
    print("\n" + "=" * 80)
    print("最終建議")
    print("=" * 80)
    print(f"\n對照組 (CK): FW = {FW_ck:.1f} g, Anth = {Anth_ck:.1f} μg/g, 總量 = {Anth_ck * FW_ck:.0f} μg/plant")

    best_overall = results_sorted_score[0]
    best_safe = safe_results[0] if safe_results else None

    print(f"\n☆ 最佳花青素總量策略 (白星):")
    print(f"  - 每天照射: {best_overall['hours']} 小時")
    print(f"  - 持續天數: {best_overall['days']} 天")
    print(f"  - 開始日期: Day {35 - best_overall['days']}")
    print(f"  - 預期鮮重: {best_overall['FW']:.1f} g ({best_overall['fw_change']:+.1f}%)")
    print(f"  - 預期花青素濃度: {best_overall['Anth']:.1f} μg/g ({best_overall['anth_change']:+.1f}%)")
    print(f"  - 預期花青素總量: {best_overall['anth_total']:.0f} μg/plant ({best_overall['anth_total_change']:+.1f}%)")

    if best_safe:
        print(f"\n★ 最佳安全策略 (黃星, 鮮重不減超過5%):")
        print(f"  - 每天照射: {best_safe['hours']} 小時")
        print(f"  - 持續天數: {best_safe['days']} 天")
        print(f"  - 開始日期: Day {35 - best_safe['days']}")
        print(f"  - 預期鮮重: {best_safe['FW']:.1f} g ({best_safe['fw_change']:+.1f}%)")
        print(f"  - 預期花青素濃度: {best_safe['Anth']:.1f} μg/g ({best_safe['anth_change']:+.1f}%)")
        print(f"  - 預期花青素總量: {best_safe['anth_total']:.0f} μg/plant ({best_safe['anth_total_change']:+.1f}%)")

    return results


if __name__ == "__main__":
    results = run_optimization()
