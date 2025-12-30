"""
敏感度分析腳本 (Sensitivity Analysis v6.7)
==========================================
版本: v6.7 Final
日期: 2025-12-23

分析 v6.7 模型中所有自定義參數的敏感性
使用 H12D3 處理組作為基準 (最能展現參數影響)
"""

import sys
sys.stdout.reconfigure(encoding='utf-8')

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import pandas as pd

from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio
from model_config import ENV_BASE, TARGETS, SIMULATION, get_env_for_treatment, ODE_SETTINGS

# =============================================================================
# v6.7 參數列表 (基於 params_config.py)
# =============================================================================
PARAMS_TO_ANALYZE = {
    # UVA-PAR轉換
    'par_conversion_factor': {'base': 1.0, 'range': 0.30, 'description': 'UVA-PAR轉換係數'},

    # Stress 損傷參數
    'stress_damage_coeff': {'base': 0.66e-6, 'range': 0.30, 'description': '損傷係數'},
    'LAI_ref_vuln': {'base': 6.5, 'range': 0.30, 'description': 'LAI參考值 (脆弱性)'},
    'n_vuln': {'base': 8, 'range': 0.20, 'description': '脆弱性指數'},
    'E_50': {'base': 475.2, 'range': 0.30, 'description': '半效能量 [kJ/m²]'},
    'E_scale': {'base': 237.6, 'range': 0.30, 'description': '能量尺度 [kJ/m²]'},
    'k_intraday': {'base': 49.0, 'range': 0.30, 'description': '日內非線性係數'},
    'm_intraday': {'base': 2.0, 'range': 0.20, 'description': '日內非線性指數'},
    'stress_nonlinear_coeff': {'base': 8.0, 'range': 0.30, 'description': 'Stress非線性係數'},
    'K_nonlinear': {'base': 0.8, 'range': 0.30, 'description': '非線性半飽和常數'},
    'circadian_disruption_factor': {'base': 3.0, 'range': 0.30, 'description': '節律損傷因子'},

    # Stress 修復參數
    'stress_repair_coeff': {'base': 1.0e-5, 'range': 0.30, 'description': '修復係數'},
    'repair_carbon_cost': {'base': 1.0e-6, 'range': 0.30, 'description': '修復碳消耗'},
    'K_carbon': {'base': 0.001, 'range': 0.30, 'description': '碳依賴修復半飽和常數'},

    # Stress 對生長抑制
    'stress_photosynthesis_inhibition': {'base': 0.66, 'range': 0.20, 'description': '光合抑制係數'},
    'stress_lai_inhibition': {'base': 0.66, 'range': 0.20, 'description': 'LAI抑制係數'},
    'K_stress': {'base': 1.9, 'range': 0.30, 'description': '抑制半飽和常數'},

    # 花青素參數 (v6.7 FW-based)
    'base_anth_rate_light': {'base': 2.0e-10, 'range': 0.30, 'description': '日間基礎合成速率'},
    'base_anth_rate_dark': {'base': 1.0e-10, 'range': 0.30, 'description': '夜間基礎合成速率'},
    'V_max_anth': {'base': 2.35e-11, 'range': 0.30, 'description': 'Stress誘導最大速率'},
    'K_stress_anth': {'base': 0.30, 'range': 0.30, 'description': 'Stress半飽和常數 (花青素)'},
    'k_deg': {'base': 3.02e-6, 'range': 0.30, 'description': '降解速率'},

    # LDMC 參數
    'dw_fw_ratio_base': {'base': 0.05, 'range': 0.20, 'description': '基礎DW:FW比例'},
    'ldmc_stress_sensitivity': {'base': 1.0, 'range': 0.30, 'description': 'LDMC敏感度'},
    'K_ldmc': {'base': 50.0, 'range': 0.30, 'description': 'LDMC半飽和Stress值'},
}

# 測試處理組 (使用 H12D3 作為主要測試對象，因為最能展現參數影響)
PRIMARY_TREATMENT = 'H12D3'
SECONDARY_TREATMENTS = ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12']

# =============================================================================
# 模擬函數
# =============================================================================
def run_simulation(p, treatment):
    """運行單次模擬並返回結果"""
    env = get_env_for_treatment(treatment)

    # 初始條件
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    C_buf_init = Xd_init * 0.1
    LAI_init = (dw_init_g / 0.01) * 0.04
    fw_total_init = fw_init_g * ENV_BASE['plant_density'] / 1000
    Anth_init = 5.0 * fw_total_init / 1e6
    initial_state = [Xd_init, C_buf_init, LAI_init, Anth_init, 0.0]  # v6.7: 5 state variables

    simulation_days = SIMULATION['days']
    transplant_day = SIMULATION['transplant_offset']
    t_start = transplant_day * 86400
    t_end = (transplant_day + simulation_days) * 86400

    try:
        sol = solve_ivp(
            fun=uva_sun_derivatives,
            t_span=(t_start, t_end),
            y0=initial_state,
            args=(p, env),
            method=ODE_SETTINGS['method'],
            max_step=ODE_SETTINGS['max_step'],
        )

        if sol.success:
            Xd_f, Cbuf_f, LAI_f, Anth_f, Stress_f = sol.y[:, -1]
            dw_fw_ratio = calculate_dynamic_dw_fw_ratio(Stress_f, p)
            dw_per_plant = Xd_f / ENV_BASE['plant_density'] * 1000
            fw_per_plant = dw_per_plant / dw_fw_ratio

            # 計算花青素濃度
            fw_total_kg = fw_per_plant / 1000 * ENV_BASE['plant_density']
            anth_ppm = (Anth_f / fw_total_kg) * 1e6 if fw_total_kg > 0 else 0

            return {
                'FW': fw_per_plant,
                'Anth': anth_ppm,
                'Stress': Stress_f,
                'LAI': LAI_f,
                'success': True
            }
    except Exception as e:
        print(f"  ⚠️  模擬失敗 ({treatment}): {str(e)}")

    return {
        'FW': np.nan,
        'Anth': np.nan,
        'Stress': np.nan,
        'LAI': np.nan,
        'success': False
    }


# =============================================================================
# 敏感性分析函數
# =============================================================================
def calculate_sensitivity(param_name, param_info, treatment):
    """
    計算參數敏感性

    敏感性定義: S = (ΔY/Y) / (ΔX/X)
    其中 Y 是輸出 (FW, Anth, Stress)
    X 是參數值
    """
    base_value = param_info['base']
    variation = param_info['range']

    low_value = base_value * (1 - variation)
    high_value = base_value * (1 + variation)

    # 基準值
    p_base = UVAParams()
    base_result = run_simulation(p_base, treatment)

    # 低值
    p_low = UVAParams()
    setattr(p_low, param_name, low_value)
    low_result = run_simulation(p_low, treatment)

    # 高值
    p_high = UVAParams()
    setattr(p_high, param_name, high_value)
    high_result = run_simulation(p_high, treatment)

    if not all([base_result['success'], low_result['success'], high_result['success']]):
        return None

    # 計算敏感性
    results = {
        'param': param_name,
        'base_value': base_value,
        'variation': variation * 100,  # %
    }

    # FW 敏感性
    fw_base = base_result['FW']
    fw_low = low_result['FW']
    fw_high = high_result['FW']
    delta_fw_ratio = (fw_high - fw_low) / fw_base if fw_base > 0 else 0
    delta_x_ratio = 2 * variation
    s_fw = delta_fw_ratio / delta_x_ratio if delta_x_ratio > 0 else 0

    results['FW_base'] = fw_base
    results['FW_low'] = fw_low
    results['FW_high'] = fw_high
    results['S_FW'] = s_fw

    # Anth 敏感性
    anth_base = base_result['Anth']
    anth_low = low_result['Anth']
    anth_high = high_result['Anth']
    delta_anth_ratio = (anth_high - anth_low) / anth_base if anth_base > 0 else 0
    s_anth = delta_anth_ratio / delta_x_ratio if delta_x_ratio > 0 else 0

    results['Anth_base'] = anth_base
    results['Anth_low'] = anth_low
    results['Anth_high'] = anth_high
    results['S_Anth'] = s_anth

    # Stress 敏感性
    stress_base = base_result['Stress']
    stress_low = low_result['Stress']
    stress_high = high_result['Stress']
    delta_stress_ratio = (stress_high - stress_low) / (stress_base + 0.01)
    s_stress = delta_stress_ratio / delta_x_ratio if delta_x_ratio > 0 else 0

    results['Stress_base'] = stress_base
    results['Stress_low'] = stress_low
    results['Stress_high'] = stress_high
    results['S_Stress'] = s_stress

    return results


# =============================================================================
# 主程序
# =============================================================================
def main():
    print("=" * 80)
    print("v6.7 敏感性分析 (Sensitivity Analysis)")
    print("=" * 80)
    print(f"測試處理組: {PRIMARY_TREATMENT}")
    print(f"參數數量: {len(PARAMS_TO_ANALYZE)}")
    print(f"變化範圍: ±20-30% (視參數而定)")
    print()

    results_list = []

    for i, (param_name, param_info) in enumerate(PARAMS_TO_ANALYZE.items(), 1):
        print(f"[{i}/{len(PARAMS_TO_ANALYZE)}] 分析參數: {param_name} ({param_info['description']})")
        print(f"  基準值: {param_info['base']}, 變化: ±{param_info['range']*100:.0f}%")

        result = calculate_sensitivity(param_name, param_info, PRIMARY_TREATMENT)

        if result:
            results_list.append(result)
            print(f"  ✅ S_FW={result['S_FW']:+.3f}, S_Stress={result['S_Stress']:+.3f}, S_Anth={result['S_Anth']:+.3f}")
        else:
            print(f"  ❌ 模擬失敗")
        print()

    # 創建 DataFrame
    df = pd.DataFrame(results_list)

    # 按照 FW 敏感性絕對值排序
    df['abs_S_FW'] = df['S_FW'].abs()
    df['abs_S_Stress'] = df['S_Stress'].abs()
    df['abs_S_Anth'] = df['S_Anth'].abs()
    df = df.sort_values('abs_S_FW', ascending=False)

    # 保存結果
    output_file = 'sensitivity_analysis_v67_results.csv'
    df.to_csv(output_file, index=False, encoding='utf-8-sig')
    print(f"結果已保存至: {output_file}")

    # 打印摘要
    print("\n" + "=" * 80)
    print("敏感性分析摘要 (按 |S_FW| 降序)")
    print("=" * 80)
    print(f"{'參數':<30} {'S_FW':>8} {'S_Stress':>10} {'S_Anth':>9} {'敏感度等級'}")
    print("-" * 80)

    for _, row in df.iterrows():
        # 敏感度分級
        if row['abs_S_FW'] > 1.0:
            level = '極高 ⚠️'
        elif row['abs_S_FW'] > 0.5:
            level = '高'
        elif row['abs_S_FW'] > 0.2:
            level = '中'
        else:
            level = '低'

        print(f"{row['param']:<30} {row['S_FW']:>+7.3f} {row['S_Stress']:>+9.2f} {row['S_Anth']:>+8.3f}  {level}")

    print("=" * 80)

    # 生成視覺化
    plot_sensitivity_results(df)

    return df


def plot_sensitivity_results(df):
    """生成敏感性分析視覺化"""
    print("\n生成敏感性分析圖表...")

    fig, axes = plt.subplots(1, 3, figsize=(15, 6))

    # 取前10個最敏感的參數
    top_params = df.head(10)

    # FW 敏感性
    axes[0].barh(top_params['param'], top_params['S_FW'],
                 color=['red' if x < 0 else 'blue' for x in top_params['S_FW']])
    axes[0].set_xlabel('Sensitivity (S_FW)')
    axes[0].set_title('Fresh Weight Sensitivity\n(Top 10 Parameters)')
    axes[0].axvline(0, color='black', linewidth=0.5)
    axes[0].grid(True, alpha=0.3, axis='x')

    # Stress 敏感性
    axes[1].barh(top_params['param'], top_params['S_Stress'],
                 color=['red' if x < 0 else 'orange' for x in top_params['S_Stress']])
    axes[1].set_xlabel('Sensitivity (S_Stress)')
    axes[1].set_title('Stress Sensitivity\n(Top 10 Parameters)')
    axes[1].axvline(0, color='black', linewidth=0.5)
    axes[1].grid(True, alpha=0.3, axis='x')

    # Anth 敏感性
    axes[2].barh(top_params['param'], top_params['S_Anth'],
                 color=['red' if x < 0 else 'green' for x in top_params['S_Anth']])
    axes[2].set_xlabel('Sensitivity (S_Anth)')
    axes[2].set_title('Anthocyanin Sensitivity\n(Top 10 Parameters)')
    axes[2].axvline(0, color='black', linewidth=0.5)
    axes[2].grid(True, alpha=0.3, axis='x')

    plt.tight_layout()
    plt.savefig('sensitivity_analysis_v67_H12D3.png', dpi=300, bbox_inches='tight')
    plt.savefig('sensitivity_analysis_v67_H12D3.pdf', bbox_inches='tight')
    plt.close()

    print("  圖表已保存: sensitivity_analysis_v67_H12D3.png/pdf")


if __name__ == "__main__":
    df_results = main()
