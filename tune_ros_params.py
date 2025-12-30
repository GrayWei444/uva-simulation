"""
ROS 動力學參數校準腳本 (v8.0)

目標: 找到 E_threshold, n_ros, k_ros_damage 使模型結果匹配實驗數據

策略:
1. E_threshold: 抗氧化劑耗竭的臨界能量
   - 應該在 6h 照射前就開始耗竭 (讓 L6D6 有效果)
   - 但在 3-4h 左右 (讓 VL3D12 前半段安全)

2. n_ros: 損傷指數
   - 需要足夠大以產生 6h vs 12h 的顯著差異
   - 但要有物理意義 (反應動力學)

3. k_ros_damage: 校準係數
   - 調整整體損傷水平

目標值 (來自 v7.0):
- 6h: factor ≈ 3.8
- 12h: factor ≈ 359
- 12h/6h 比值 ≈ 95
"""

import numpy as np
from simulate_uva_model_v8_ros import calculate_ros_damage_factor, UVAParams, ALL_PARAMS

# v7.0 目標值
I_UVA = 11.0  # W/m²
E_6h = I_UVA * 6 * 3600    # 237,600 J/m²
E_12h = I_UVA * 12 * 3600  # 475,200 J/m²

target_factor_6h = 3.8
target_factor_12h = 359.0

print("=" * 70)
print("ROS 動力學參數校準")
print("=" * 70)
print(f"\n目標:")
print(f"  6h factor: {target_factor_6h}")
print(f"  12h factor: {target_factor_12h}")
print(f"  12h/6h 比值: {target_factor_12h/target_factor_6h:.1f}")

print(f"\n能量值:")
print(f"  E_6h = {E_6h} J/m²")
print(f"  E_12h = {E_12h} J/m²")

# 分析: 給定 E_threshold, n_ros, 可以計算 k_ros_damage
print("\n" + "=" * 70)
print("參數分析")
print("=" * 70)

print("""
公式: factor = 1 + k × (E - E_threshold)^n

給定 E_threshold 和 n，可以從兩個目標值計算 k:

factor_6h = 1 + k × (E_6h - E_threshold)^n
factor_12h = 1 + k × (E_12h - E_threshold)^n

解:
k = (factor - 1) / (E - E_threshold)^n

一致性檢查:
如果 k_from_6h ≈ k_from_12h，則參數合理
""")

def analyze_params(E_threshold, n):
    """分析給定 E_threshold 和 n 的一致性"""
    E_excess_6h = max(0, E_6h - E_threshold)
    E_excess_12h = max(0, E_12h - E_threshold)

    if E_excess_6h <= 0:
        return None, None, None

    k_from_6h = (target_factor_6h - 1) / (E_excess_6h ** n)
    k_from_12h = (target_factor_12h - 1) / (E_excess_12h ** n)

    ratio = k_from_12h / k_from_6h if k_from_6h > 0 else float('inf')

    return k_from_6h, k_from_12h, ratio

print("\n測試不同參數組合:")
print("-" * 70)
print(f"{'E_threshold':<12} {'n':<8} {'k_from_6h':<14} {'k_from_12h':<14} {'比值':<10}")
print("-" * 70)

for E_th in [100000, 150000, 180000, 200000]:
    for n in [2.0, 2.5, 3.0, 3.5, 4.0]:
        k_6h, k_12h, ratio = analyze_params(E_th, n)
        if k_6h is not None:
            print(f"{E_th:<12.0f} {n:<8.1f} {k_6h:<14.2e} {k_12h:<14.2e} {ratio:<10.2f}")

print("-" * 70)

# 關鍵洞察
print("""
關鍵洞察:
=========

v7.0 的 12h/6h = 95 比值意味著:
- 12h 的損傷是 6h 的 95 倍
- 這需要非常陡峭的曲線

對於純冪律 factor = 1 + k × E^n:
- (E_12h/E_6h)^n = (475200/237600)^n = 2^n
- 要達到 95 倍，需要 2^n ≈ 95，即 n ≈ 6.6

對於閾值冪律 factor = 1 + k × (E - E_threshold)^n:
- 比值 = (E_12h - E_th)^n / (E_6h - E_th)^n
- 當 E_th 接近 E_6h 時，分母很小，比值可以很大

結論: E_threshold 應該接近但不超過 E_6h
""")

# 搜索最佳參數
print("=" * 70)
print("最佳參數搜索")
print("=" * 70)

def objective(params):
    """目標函數: 最小化兩個目標的加權誤差"""
    E_th, n, log_k = params
    k = 10 ** log_k

    E_excess_6h = max(0, E_6h - E_th)
    E_excess_12h = max(0, E_12h - E_th)

    if E_excess_6h <= 0:
        return 1e10

    factor_6h = 1 + k * (E_excess_6h ** n)
    factor_12h = 1 + k * (E_excess_12h ** n)

    # 加權誤差 (對數尺度因為數值差異大)
    err_6h = (np.log(factor_6h) - np.log(target_factor_6h)) ** 2
    err_12h = (np.log(factor_12h) - np.log(target_factor_12h)) ** 2

    return err_6h + err_12h

# 網格搜索
best_error = float('inf')
best_params = None

for E_th in np.linspace(150000, 230000, 20):
    for n in np.linspace(2.0, 5.0, 15):
        for log_k in np.linspace(-25, -15, 20):
            err = objective([E_th, n, log_k])
            if err < best_error:
                best_error = err
                best_params = (E_th, n, 10 ** log_k)

E_th_opt, n_opt, k_opt = best_params
print(f"\n最佳參數:")
print(f"  E_threshold = {E_th_opt:.0f} J/m² ({E_th_opt/(I_UVA*3600):.2f}h @ {I_UVA} W/m²)")
print(f"  n_ros = {n_opt:.2f}")
print(f"  k_ros_damage = {k_opt:.2e}")

# 驗證
E_excess_6h = max(0, E_6h - E_th_opt)
E_excess_12h = max(0, E_12h - E_th_opt)
factor_6h_opt = 1 + k_opt * (E_excess_6h ** n_opt)
factor_12h_opt = 1 + k_opt * (E_excess_12h ** n_opt)

print(f"\n驗證:")
print(f"  6h factor: {factor_6h_opt:.1f} (目標: {target_factor_6h})")
print(f"  12h factor: {factor_12h_opt:.1f} (目標: {target_factor_12h})")
print(f"  12h/6h 比值: {factor_12h_opt/factor_6h_opt:.1f} (目標: {target_factor_12h/target_factor_6h:.1f})")

# 完整曲線比較
print("\n" + "=" * 70)
print("完整曲線比較 (使用最佳參數)")
print("=" * 70)

print(f"\n{'時間':<8} {'能量(J/m²)':<12} {'v7.0 factor':<14} {'v8.0 factor':<14}")
print("-" * 56)

k_day_v7 = 1.0e-5
n_day_v7 = 7.0

for hours in [0, 1, 2, 3, 4, 5, 6, 8, 10, 12]:
    E = I_UVA * hours * 3600

    # v7.0 factor
    factor_v7 = 1.0 + k_day_v7 * (hours ** n_day_v7)

    # v8.0 factor
    E_excess = max(0, E - E_th_opt)
    factor_v8 = 1.0 + k_opt * (E_excess ** n_opt) if E_excess > 0 else 1.0

    print(f"{hours}h{'':<6} {E:<12.0f} {factor_v7:<14.1f} {factor_v8:<14.1f}")

print("-" * 56)

# 建議的參數更新
print("\n" + "=" * 70)
print("建議的參數更新 (複製到 simulate_uva_model_v8_ros.py)")
print("=" * 70)
print(f"""
    'E_threshold': {E_th_opt:.0f},             # 抗氧化劑耗竭臨界能量 [J/m²]
    'n_ros': {n_opt:.2f},                        # ROS 損傷指數 (基於反應動力學)
    'k_ros_damage': {k_opt:.2e},               # ROS 損傷係數 [1/(J/m²)^n]
""")
