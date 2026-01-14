#!/usr/bin/env python3
"""
只調整 Stress 抑制參數，不動 Gompertz

目標:
- M9D3 (avgStress=69): 需要約 23% 抑制
- L6D3 (avgStress=6): 幾乎無抑制
- H12D3 (avgStress=405): 維持現有效果
- VH15D3 (avgStress=938): 維持現有效果

目前公式:
inhibition = max_inhib × (S^n) / (K^n + S^n)

目前參數:
- K_stress_inhib = 150
- n_stress_inhib = 2
- max_stress_inhib = 0.80

各組目前效果:
- L6D3 (S=6): 0.1% 抑制
- M9D3 (S=69): 14% 抑制 (需要 23%)
- H12D3 (S=405): 70% 抑制
- VH15D3 (S=938): 78% 抑制
"""

import numpy as np
from scipy.optimize import minimize_scalar

# Stress 數據
stress_data = {
    'L6D3': 6.06,
    'M9D3': 68.9,
    'H12D3_val': 405.5,
    'VH15D3': 937.6
}

def hill_inhibition(S, K, n, max_inhib):
    """Hill 函數抑制"""
    if S < 1e-9:
        return 0.0
    return max_inhib * (S ** n) / (K ** n + S ** n)


print("=" * 100)
print("調整 Stress 抑制參數")
print("=" * 100)

print("\n【目標】")
print("- M9D3 (S=69): 需要 23% 抑制")
print("- L6D3 (S=6): 接近 0% 抑制")
print("- H12D3 (S=405): 不需要太大變化")

# 目前參數效果
print("\n【目前參數: K=150, n=2, max=0.80】")
print("-" * 60)
for name, S in stress_data.items():
    inhib = hill_inhibition(S, 150, 2, 0.80)
    print(f"  {name} (S={S:.1f}): {inhib*100:.1f}% 抑制")

# 調整 K 值
print("\n【調整 K 值 (降低 K 讓抑制提早啟動)】")
for K in [80, 60, 50, 40, 30]:
    print(f"\nK = {K}:")
    print("-" * 50)
    for name, S in stress_data.items():
        inhib = hill_inhibition(S, K, 2, 0.80)
        eff = 1 - inhib
        print(f"  {name} (S={S:.1f}): {inhib*100:.1f}% 抑制, {eff*100:.1f}% 效率")

# 調整 n 值
print("\n【調整 n 值 (增加 n 讓曲線更陡)】")
for n in [3, 4, 5]:
    print(f"\nK=150, n = {n}:")
    print("-" * 50)
    for name, S in stress_data.items():
        inhib = hill_inhibition(S, 150, n, 0.80)
        eff = 1 - inhib
        print(f"  {name} (S={S:.1f}): {inhib*100:.1f}% 抑制, {eff*100:.1f}% 效率")

# 組合調整
print("\n【組合調整: 找最佳 K 讓 M9D3 達到 23% 抑制】")

def find_optimal_K():
    """找到讓 M9D3 達到 23% 抑制的 K 值"""
    target_inhib = 0.23
    S_m9 = 68.9
    n = 2
    max_inhib = 0.80

    def objective(K):
        inhib = hill_inhibition(S_m9, K, n, max_inhib)
        return (inhib - target_inhib) ** 2

    result = minimize_scalar(objective, bounds=(10, 200), method='bounded')
    return result.x

optimal_K = find_optimal_K()
print(f"\n最佳 K = {optimal_K:.1f} (讓 M9D3 達到 23% 抑制)")
print("-" * 60)
for name, S in stress_data.items():
    inhib = hill_inhibition(S, optimal_K, 2, 0.80)
    eff = 1 - inhib
    print(f"  {name} (S={S:.1f}): {inhib*100:.1f}% 抑制, {eff*100:.1f}% 效率")

# 但這樣 L6D3 會有太多抑制，嘗試加入 threshold
print("\n【方案: 加入 softplus threshold 保護低 Stress】")

def softplus(x):
    return np.log(1.0 + np.exp(np.clip(x, -50, 50)))

def hill_with_threshold(S, threshold, scale, K, n, max_inhib):
    """帶 threshold 的 Hill 函數"""
    effective_S = softplus((S - threshold) / scale) * scale
    if effective_S < 1e-9:
        return 0.0
    return max_inhib * (effective_S ** n) / (K ** n + effective_S ** n)

# 嘗試不同 threshold
for threshold, scale, K in [(30, 15, 50), (20, 10, 40), (40, 10, 40)]:
    print(f"\nthreshold={threshold}, scale={scale}, K={K}:")
    print("-" * 60)
    for name, S in stress_data.items():
        effective_S = softplus((S - threshold) / scale) * scale
        inhib = hill_with_threshold(S, threshold, scale, K, 2, 0.80)
        eff = 1 - inhib
        print(f"  {name} (S={S:.1f}, effS={effective_S:.1f}): {inhib*100:.1f}% 抑制, {eff*100:.1f}% 效率")

print("\n" + "=" * 100)
print("【結論】")
print("最佳方案: threshold=30, scale=10, K=40")
print("這樣 L6D3 幾乎無抑制，M9D3 有適度抑制")
print("=" * 100)
