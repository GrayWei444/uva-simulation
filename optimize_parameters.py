"""
================================================================================
萵苣UVA模型參數估計腳本 (Parameter Estimation for Lettuce UVA Model)
================================================================================
方法論 (適用於學術論文):

1. 資料分割策略:
   - 訓練集 (Training Set): CK, L6D6, L6D6-N
     * CK: 對照組，無UVA效應，用於校準基礎生長模型
     * L6D6: 日間UVA處理，用於校準UVA-PAR效應和花青素誘導
     * L6D6-N: 夜間UVA處理，用於校準生理節律效應

   - 驗證集 (Validation Set): VL3D12, L6D12
     * 用於測試模型對不同照射時長的泛化能力

   - 測試集 (Test Set): H12D3, SIN, INT
     * 用於評估模型對極端條件和動態照射的預測能力
     * 這些組別可能需要額外機制 (如損傷模型)

2. 優化方法:
   - 使用 scipy.optimize.differential_evolution 全局優化
   - 目標函數: 加權最小二乘法 (Weighted Least Squares)
   - 使用 Bootstrap 估計參數不確定性

3. 模型選擇:
   - 使用 AIC/BIC 評估模型複雜度與擬合度的權衡
   - 只保留統計顯著的機制

4. 結果報告:
   - 參數估計值 ± 標準誤
   - R² 值
   - RMSE
   - 殘差分析

================================================================================
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution, minimize
from scipy.stats import t as t_dist

# 使用共用設定模組
from model_config import (
    ENV_BASE, TREATMENT_CONFIGS, TARGETS, ODE_SETTINGS, SIMULATION,
    get_env_for_treatment
)
import warnings
import sys
sys.path.insert(0, '.')
from simulate_uva_model import UVAParams, uva_sun_derivatives, calculate_dynamic_dw_fw_ratio

# ==============================================================================
# 實驗數據 (來自論文)
# ==============================================================================
# 資料集分割策略:
# - 訓練集: CK, L6D6, L6D6-N, H12D3
#   * CK: 基礎生長模型
#   * L6D6: UVA-PAR增益、花青素誘導、LDMC效應
#   * L6D6-N: 生理節律抑制、夜間效應
#   * H12D3: 損傷機制 (高劑量脅迫)
# - 驗證集: VL3D12, L6D12 (驗證時間尺度泛化: 6天 vs 12天)
# - 測試集: SIN, INT (驗證動態照射模式: 連續 vs 非連續)

EXPERIMENTAL_DATA = {
    # 訓練集 (Training Set) - 涵蓋所有關鍵機制
    # 更新日期: 2025-12-05 (新實驗數據)
    'CK': {'FW': 87.0, 'FW_se': 3.5, 'Anth': 43.0, 'Anth_se': 2.0, 'set': 'train'},
    'L6D6': {'FW': 91.4, 'FW_se': 4.0, 'Anth': 49.40, 'Anth_se': 2.5, 'set': 'train'},
    'L6D6-N': {'FW': 80.8, 'FW_se': 3.8, 'Anth': 49.26, 'Anth_se': 2.2, 'set': 'train'},
    'H12D3': {'FW': 60.6, 'FW_se': 5.0, 'Anth': 65.06, 'Anth_se': 4.0, 'set': 'train'},  # 損傷機制必須在訓練集

    # 驗證集 (Validation Set) - 驗證時間尺度泛化
    # 重要: VL3D12, L6D12 的 FW 大幅下降!
    'VL3D12': {'FW': 67.0, 'FW_se': 3.6, 'Anth': 51.90, 'Anth_se': 2.1, 'set': 'validation'},
    'L6D12': {'FW': 60.4, 'FW_se': 4.2, 'Anth': 55.92, 'Anth_se': 3.0, 'set': 'validation'},

    # 動態照射模式
    'SIN': {'FW': 63.6, 'FW_se': 3.7, 'Anth': 53.26, 'Anth_se': 2.3, 'set': 'test'},
    'INT': {'FW': 62.8, 'FW_se': 3.9, 'Anth': 52.19, 'Anth_se': 2.4, 'set': 'test'},
}

# ==============================================================================
# 環境配置 (從 model_config 匯入)
# ==============================================================================
# ENV_BASE, TREATMENT_CONFIGS 已從 model_config 匯入

# ==============================================================================
# 待估計的參數 (基於物理/生物意義分組)
# ==============================================================================
# 核心參數 (從訓練集校準):
# - CK: c_alpha, base_anth_rate_light, base_anth_rate_dark
# - L6D6: par_conversion_factor, dw_fw_ratio_increase_per_dose, S_max_uva, K_m, n_hill
# - L6D6-N: night_uva_base_inhibition
# - H12D3: damage_sensitivity (損傷機制)

CORE_PARAMS = {
    # 基礎生長 (從CK校準)
    'c_alpha': (0.7, 0.95),  # 光合效率

    # UVA-PAR效應 (從L6D6校準)
    'par_conversion_factor': (8.0, 20.0),  # UVA到PAR轉換

    # LDMC效應
    'dw_fw_ratio_increase_per_dose': (1e-10, 5e-9),

    # 花青素基礎合成 (從CK校準)
    'base_anth_rate_light': (1.0e-7, 2.5e-7),
    'base_anth_rate_dark': (1.5e-8, 4.0e-8),

    # UVA誘導花青素 (從L6D6校準)
    'S_max_uva': (4.0e-8, 9.0e-8),
    'K_m': (80, 200),
    'n_hill': (1.0, 1.8),

    # 夜間效應 (從L6D6-N校準)
    'night_uva_base_inhibition': (0.05, 0.25),

    # 損傷機制 (從H12D3校準)
    # H12D3: FW下降33% (87->58.6g)，花青素上升51% (43->65.1ppm)
    'damage_sensitivity': (1e-8, 1e-6),  # 損傷敏感度

    # 脅迫恢復係數 (影響間歇照射時的 Stress 修復速度)
    'recovery_coeff': (1e-5, 5e-4),

    # 脅迫-花青素連續效應參數 (用於 SIN/INT)
    # stress_amplification = (Stress / stress_reference)^stress_exponent
    'stress_reference': (10.0, 50.0),   # 參考 Stress 值
    'stress_exponent': (0.2, 0.8),       # 指數 (控制非線性程度)
}

PARAM_NAMES = list(CORE_PARAMS.keys())

# ==============================================================================
# 模擬函數
# ==============================================================================
def run_simulation(treatment, params_dict, suppress_warnings=True):
    """執行單一處理組的模擬"""
    if suppress_warnings:
        warnings.filterwarnings('ignore')

    p = UVAParams()

    # 設定優化參數
    for name, value in params_dict.items():
        if hasattr(p, name):
            setattr(p, name, value)

    # 環境配置 (使用 model_config)
    env = get_env_for_treatment(treatment)

    # 初始條件
    fw_init_g = SIMULATION['initial_fw_g']
    dw_init_g = fw_init_g * p.dw_fw_ratio_base
    Xd_init = dw_init_g / 1000 * ENV_BASE['plant_density']
    LAI_init = (dw_init_g / 0.01) * 0.04
    Anth_init = 2.0 * (fw_init_g * ENV_BASE['plant_density'] / 1000) / 1000
    initial_state = [Xd_init, 0.0, LAI_init, Anth_init, 0.0, 0.0, 0.0, 0.0]

    try:
        sol = solve_ivp(
            fun=uva_sun_derivatives,
            t_span=(0, SIMULATION['days'] * 86400),
            y0=initial_state,
            args=(p, env),
            method=ODE_SETTINGS['method'],
            max_step=ODE_SETTINGS['max_step'],
            t_eval=np.array([SIMULATION['days'] * 86400])
        )

        if not sol.success:
            return None, None

        # 計算結果
        sim_D_UVA = sol.y[4, -1]
        sim_dw_fw_ratio = calculate_dynamic_dw_fw_ratio(sim_D_UVA, p)
        sim_dw_per_plant = sol.y[0, -1] / ENV_BASE['plant_density'] * 1000
        sim_fw = sim_dw_per_plant / sim_dw_fw_ratio

        sim_fw_kg_per_m2 = sim_fw * ENV_BASE['plant_density'] / 1000
        sim_anth_mg_per_m2 = sol.y[3, -1] * 1000
        sim_anth_ppm = sim_anth_mg_per_m2 / (sim_fw_kg_per_m2 + 1e-9)

        return sim_fw, sim_anth_ppm
    except Exception:
        return None, None

# ==============================================================================
# 目標函數 (加權最小二乘法)
# ==============================================================================
def weighted_residuals(x, dataset='train'):
    """計算加權殘差 (用於最小二乘法)"""
    params_dict = {name: x[i] for i, name in enumerate(PARAM_NAMES)}

    residuals = []
    for treatment, data in EXPERIMENTAL_DATA.items():
        if data['set'] != dataset:
            continue

        sim_fw, sim_anth = run_simulation(treatment, params_dict)

        if sim_fw is None:
            return np.array([1e6] * 6)  # 大殘差表示失敗

        # 加權殘差: 殘差 / 標準誤
        fw_residual = (sim_fw - data['FW']) / data['FW_se']
        anth_residual = (sim_anth - data['Anth']) / data['Anth_se']

        residuals.extend([fw_residual, anth_residual])

    return np.array(residuals)

def objective_wls(x):
    """加權最小二乘目標函數"""
    residuals = weighted_residuals(x, 'train')
    return np.sum(residuals**2)

def objective_rmse(x, dataset='train'):
    """RMSE 目標函數 (用於報告)"""
    params_dict = {name: x[i] for i, name in enumerate(PARAM_NAMES)}

    errors = []
    for treatment, data in EXPERIMENTAL_DATA.items():
        if data['set'] != dataset:
            continue

        sim_fw, sim_anth = run_simulation(treatment, params_dict)

        if sim_fw is None:
            return 1e10

        # 相對誤差
        fw_err = (sim_fw - data['FW']) / data['FW']
        anth_err = (sim_anth - data['Anth']) / data['Anth']

        errors.extend([fw_err**2, anth_err**2])

    return np.sqrt(np.mean(errors)) * 100  # 百分比

# ==============================================================================
# 統計分析函數
# ==============================================================================
def compute_statistics(best_params, x_opt):
    """計算統計量: R², RMSE, AIC, BIC"""
    n_params = len(PARAM_NAMES)

    # 計算各數據集的統計量
    results = {}

    for dataset in ['train', 'validation']:
        n_obs = sum(1 for d in EXPERIMENTAL_DATA.values() if d['set'] == dataset) * 2

        if n_obs == 0:
            continue

        # 計算殘差
        residuals = weighted_residuals(x_opt, dataset)
        ss_res = np.sum(residuals**2)

        # 計算 R² (相對於均值的改進)
        y_true = []
        y_pred = []
        for treatment, data in EXPERIMENTAL_DATA.items():
            if data['set'] != dataset:
                continue
            sim_fw, sim_anth = run_simulation(treatment, best_params)
            y_true.extend([data['FW'], data['Anth']])
            y_pred.extend([sim_fw, sim_anth])

        y_true = np.array(y_true)
        y_pred = np.array(y_pred)
        ss_tot = np.sum((y_true - np.mean(y_true))**2)
        r2 = 1 - np.sum((y_true - y_pred)**2) / ss_tot

        # RMSE
        rmse = objective_rmse(x_opt, dataset)

        # AIC / BIC (只對訓練集有意義)
        if dataset == 'train':
            log_likelihood = -0.5 * n_obs * np.log(ss_res / n_obs)
            aic = 2 * n_params - 2 * log_likelihood
            bic = n_params * np.log(n_obs) - 2 * log_likelihood
        else:
            aic = bic = None

        results[dataset] = {
            'n_obs': n_obs,
            'ss_res': ss_res,
            'r2': r2,
            'rmse': rmse,
            'aic': aic,
            'bic': bic,
        }

    return results

def bootstrap_confidence_intervals(x_opt, n_bootstrap=100, alpha=0.05):
    """使用 Bootstrap 估計參數置信區間"""
    print(f"\n正在進行 Bootstrap 分析 (n={n_bootstrap})...")

    # 收集訓練集數據
    train_data = [(t, d) for t, d in EXPERIMENTAL_DATA.items() if d['set'] == 'train']
    n_treatments = len(train_data)

    bootstrap_params = []

    for i in range(n_bootstrap):
        if (i + 1) % 20 == 0:
            print(f"  進度: {i+1}/{n_bootstrap}")

        # 重採樣 (對處理組進行重採樣)
        indices = np.random.choice(n_treatments, n_treatments, replace=True)

        # 創建重採樣的目標函數
        sampled_treatments = [train_data[j][0] for j in indices]

        def objective_bootstrap(x):
            params_dict = {name: x[k] for k, name in enumerate(PARAM_NAMES)}
            residuals = []
            for treatment in sampled_treatments:
                data = EXPERIMENTAL_DATA[treatment]
                sim_fw, sim_anth = run_simulation(treatment, params_dict)
                if sim_fw is None:
                    return 1e10
                fw_residual = (sim_fw - data['FW']) / data['FW_se']
                anth_residual = (sim_anth - data['Anth']) / data['Anth_se']
                residuals.extend([fw_residual**2, anth_residual**2])
            return np.sum(residuals)

        # 從原始最佳解出發進行局部優化
        bounds = [CORE_PARAMS[name] for name in PARAM_NAMES]
        result = minimize(
            objective_bootstrap,
            x_opt,
            method='L-BFGS-B',
            bounds=bounds,
            options={'maxiter': 50}
        )

        if result.success:
            bootstrap_params.append(result.x)

    bootstrap_params = np.array(bootstrap_params)

    # 計算置信區間
    ci_lower = np.percentile(bootstrap_params, 100 * alpha / 2, axis=0)
    ci_upper = np.percentile(bootstrap_params, 100 * (1 - alpha / 2), axis=0)
    std_error = np.std(bootstrap_params, axis=0)

    return {
        'params': bootstrap_params,
        'ci_lower': ci_lower,
        'ci_upper': ci_upper,
        'std_error': std_error,
    }

# ==============================================================================
# 主優化函數
# ==============================================================================
def run_optimization(maxiter=100, popsize=15, seed=42, do_bootstrap=True):
    """執行完整的參數估計流程"""
    print("="*80)
    print("萵苣UVA模型參數估計 (Parameter Estimation)")
    print("="*80)
    print("\n1. 資料集分割:")
    for dataset in ['train', 'validation', 'test']:
        treatments = [t for t, d in EXPERIMENTAL_DATA.items() if d['set'] == dataset]
        print(f"   {dataset.capitalize():12s}: {treatments}")

    print(f"\n2. 待估計參數 ({len(PARAM_NAMES)}個):")
    for name in PARAM_NAMES:
        bounds = CORE_PARAMS[name]
        print(f"   {name:35s}: [{bounds[0]:.2e}, {bounds[1]:.2e}]")

    print("\n3. 優化方法: Differential Evolution + L-BFGS-B polishing")
    print("   目標函數: 加權最小二乘法 (Weighted Least Squares)")
    print("="*80)

    # 執行全局優化
    print("\n正在進行全局優化...")
    bounds = [CORE_PARAMS[name] for name in PARAM_NAMES]

    # 進度回調函數
    generation_counter = [0]  # 使用列表以便在回調中修改

    def progress_callback(xk, convergence):
        """每代完成後輸出進度"""
        generation_counter[0] += 1
        gen = generation_counter[0]
        current_obj = objective_wls(xk)
        rmse = objective_rmse(xk, 'train')
        print(f"  [世代 {gen:3d}] 目標函數: {current_obj:10.4f}, 訓練RMSE: {rmse:6.2f}%, 收斂度: {convergence:.6f}")
        import sys
        sys.stdout.flush()
        return False  # 繼續優化

    result = differential_evolution(
        objective_wls,
        bounds,
        maxiter=maxiter,
        popsize=popsize,
        seed=seed,
        disp=False,  # 關閉內建輸出，使用我們的回調
        callback=progress_callback,
        workers=1,
        updating='deferred',
        polish=True,
        tol=1e-6,
    )

    x_opt = result.x
    best_params = {name: x_opt[i] for i, name in enumerate(PARAM_NAMES)}

    print("\n" + "="*80)
    print("4. 優化結果")
    print("="*80)

    # 計算統計量
    stats = compute_statistics(best_params, x_opt)

    print("\n參數估計值:")
    print("-"*70)
    print(f"{'參數名稱':<35s} {'估計值':>15s} {'範圍':>20s}")
    print("-"*70)
    for i, name in enumerate(PARAM_NAMES):
        bounds = CORE_PARAMS[name]
        print(f"{name:<35s} {x_opt[i]:>15.6e} [{bounds[0]:.1e}, {bounds[1]:.1e}]")

    # Bootstrap 置信區間
    if do_bootstrap:
        bootstrap = bootstrap_confidence_intervals(x_opt, n_bootstrap=50)

        print("\n參數估計值 (含95%置信區間):")
        print("-"*80)
        print(f"{'參數名稱':<35s} {'估計值':>12s} {'標準誤':>12s} {'95% CI':>20s}")
        print("-"*80)
        for i, name in enumerate(PARAM_NAMES):
            ci_str = f"[{bootstrap['ci_lower'][i]:.2e}, {bootstrap['ci_upper'][i]:.2e}]"
            print(f"{name:<35s} {x_opt[i]:>12.4e} {bootstrap['std_error'][i]:>12.4e} {ci_str:>20s}")

    # 模型擬合統計
    print("\n" + "="*80)
    print("5. 模型擬合統計")
    print("="*80)

    for dataset, stat in stats.items():
        print(f"\n{dataset.upper()} SET:")
        print(f"  觀測數 (n):    {stat['n_obs']}")
        print(f"  R-squared:     {stat['r2']:.4f}")
        print(f"  RMSE:          {stat['rmse']:.2f}%")
        if stat['aic'] is not None:
            print(f"  AIC:           {stat['aic']:.2f}")
            print(f"  BIC:           {stat['bic']:.2f}")

    # 詳細預測結果
    print("\n" + "="*80)
    print("6. 詳細預測結果")
    print("="*80)
    print(f"\n{'處理組':<10s} | {'資料集':<10s} | {'目標FW':>8s} | {'預測FW':>8s} | {'誤差':>8s} | {'目標Anth':>9s} | {'預測Anth':>9s} | {'誤差':>8s}")
    print("-"*95)

    for treatment in ['CK', 'L6D6', 'L6D6-N', 'H12D3', 'VL3D12', 'L6D12', 'SIN', 'INT']:
        data = EXPERIMENTAL_DATA[treatment]
        sim_fw, sim_anth = run_simulation(treatment, best_params)

        fw_err = (sim_fw - data['FW']) / data['FW'] * 100
        anth_err = (sim_anth - data['Anth']) / data['Anth'] * 100

        status = "✓" if abs(fw_err) < 5 and abs(anth_err) < 5 else "!"
        print(f"{treatment:<10s} | {data['set']:<10s} | {data['FW']:>8.1f} | {sim_fw:>8.2f} | {fw_err:>+7.1f}% | {data['Anth']:>9.1f} | {sim_anth:>9.2f} | {anth_err:>+7.1f}% {status}")

    print("-"*95)

    # 輸出可複製的參數
    print("\n" + "="*80)
    print("7. 可複製到 simulate_uva_model.py 的參數")
    print("="*80)
    for name, value in best_params.items():
        print(f"self.{name} = {value:.6e}")

    return best_params, stats, result

# ==============================================================================
# 主程序
# ==============================================================================
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='萵苣UVA模型參數估計')
    parser.add_argument('--maxiter', type=int, default=50, help='最大迭代次數')
    parser.add_argument('--popsize', type=int, default=10, help='種群大小')
    parser.add_argument('--seed', type=int, default=42, help='隨機種子')
    parser.add_argument('--no-bootstrap', action='store_true', help='跳過Bootstrap分析')
    parser.add_argument('--evaluate-only', action='store_true', help='只評估當前參數')
    args = parser.parse_args()

    if args.evaluate_only:
        # 使用當前模型參數評估
        p = UVAParams()
        current_params = {name: getattr(p, name) for name in PARAM_NAMES}
        print("\n當前參數:")
        for name, value in current_params.items():
            print(f"  {name}: {value:.6e}")

        x_current = np.array([current_params[name] for name in PARAM_NAMES])
        stats = compute_statistics(current_params, x_current)

        print("\n模型擬合統計:")
        for dataset, stat in stats.items():
            print(f"\n{dataset.upper()} SET:")
            print(f"  R-squared: {stat['r2']:.4f}")
            print(f"  RMSE:  {stat['rmse']:.2f}%")
    else:
        # 執行優化
        run_optimization(
            maxiter=args.maxiter,
            popsize=args.popsize,
            seed=args.seed,
            do_bootstrap=not args.no_bootstrap
        )
