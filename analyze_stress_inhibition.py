#!/usr/bin/env python3
"""
分析 Stress 抑制花青素合成的機制

問題: M9D3 模擬偏高 30%
假設: Stress 抑制機制在 M9D3 (avgStress ≈ 69) 附近效果不足

目前公式:
stress_inhibition = max_inhib × (S^n) / (K^n + S^n)

參數:
- K_stress_inhib = 150
- n_stress_inhib = 2
- max_stress_inhib = 0.80

各組 avgStress:
- VL3D3: 3.03
- L6D3: 6.06
- M9D3: 68.9  ← 轉折點
- H12D3: 405.5
- VH15D3: 937.6
"""

import numpy as np

# 目前參數
K = 150.0
n = 2.0
max_inhib = 0.80

# 各組 avgStress (從之前的分析)
stress_data = {
    'CK_val': 0.0,
    'VL3D3': 3.03,
    'L6D3': 6.06,
    'M9D3': 68.9,
    'H12D3_val': 405.5,
    'VH15D3': 937.6
}

print("=" * 80)
print("Stress 抑制花青素合成分析")
print("=" * 80)

print("\n【目前參數】")
print(f"K_stress_inhib = {K}")
print(f"n_stress_inhib = {n}")
print(f"max_stress_inhib = {max_inhib}")

print("\n【目前效果】")
print("-" * 60)
print(f"{'處理組':>10} | {'avgStress':>10} | {'抑制程度':>10} | {'效率':>10}")
print("-" * 60)

for treatment, stress in stress_data.items():
    if stress == 0:
        inhibition = 0
    else:
        inhibition = max_inhib * (stress ** n) / (K ** n + stress ** n)
    efficiency = 1.0 - inhibition
    print(f"{treatment:>10} | {stress:>10.1f} | {inhibition*100:>9.1f}% | {efficiency*100:>9.1f}%")

print("-" * 60)

# 問題分析
print("\n【問題分析】")
print("M9D3 avgStress=69，但抑制只有 14.4%")
print("這意味著 M9D3 的花青素合成效率仍有 85.6%")
print("但 M9D3 應該是轉折點，效率應該開始明顯下降")

# 嘗試不同參數
print("\n" + "=" * 80)
print("嘗試調整參數")
print("=" * 80)

# 方案1: 降低 K，讓抑制提早啟動
print("\n【方案1: 降低 K_stress_inhib (提早啟動)】")
for K_new in [50, 30, 20]:
    print(f"\nK = {K_new}:")
    print("-" * 50)
    for treatment, stress in [('M9D3', 68.9), ('H12D3_val', 405.5), ('VH15D3', 937.6)]:
        if stress == 0:
            inhibition = 0
        else:
            inhibition = max_inhib * (stress ** n) / (K_new ** n + stress ** n)
        efficiency = 1.0 - inhibition
        print(f"  {treatment}: 抑制 {inhibition*100:.1f}%, 效率 {efficiency*100:.1f}%")

# 方案2: 使用分段函數 - 低 Stress 時無抑制，高 Stress 時抑制
print("\n【方案2: 加入 threshold (閾值下無抑制)】")
for threshold in [30, 50, 60]:
    print(f"\nthreshold = {threshold}, K = 100:")
    K_new = 100
    print("-" * 50)
    for treatment, stress in [('M9D3', 68.9), ('H12D3_val', 405.5), ('VH15D3', 937.6)]:
        if stress <= threshold:
            inhibition = 0
        else:
            effective_stress = stress - threshold
            inhibition = max_inhib * (effective_stress ** n) / (K_new ** n + effective_stress ** n)
        efficiency = 1.0 - inhibition
        print(f"  {treatment}: 抑制 {inhibition*100:.1f}%, 效率 {efficiency*100:.1f}%")

# 方案3: 使用 softplus 連續過渡 (符合 CLAUDE.md)
print("\n【方案3: softplus 連續過渡 (符合 CLAUDE.md 無硬閾值原則)】")

def softplus(x):
    return np.log(1.0 + np.exp(np.clip(x, -50, 50)))

for threshold, scale in [(50, 20), (40, 15), (30, 10)]:
    print(f"\nthreshold = {threshold}, scale = {scale}, K = 80:")
    K_new = 80
    print("-" * 50)
    for treatment, stress in [('L6D3', 6.06), ('M9D3', 68.9), ('H12D3_val', 405.5), ('VH15D3', 937.6)]:
        # softplus 過渡: 當 stress < threshold 時接近 0，超過後線性增長
        effective_stress = softplus((stress - threshold) / scale) * scale
        if effective_stress < 1e-9:
            inhibition = 0
        else:
            inhibition = max_inhib * (effective_stress ** n) / (K_new ** n + effective_stress ** n)
        efficiency = 1.0 - inhibition
        print(f"  {treatment} (S={stress:.1f}): 有效S={effective_stress:.1f}, 抑制 {inhibition*100:.1f}%, 效率 {efficiency*100:.1f}%")

print("\n" + "=" * 80)
print("建議: 使用 softplus 連續過渡，設定 threshold=50, scale=15")
print("這樣 M9D3 開始有明顯抑制，而 L6D3 幾乎不受影響")
print("=" * 80)
