#!/usr/bin/env python3
"""
調整 Gompertz 參數，讓 M9D3 的轉折更明顯

目標:
- 9h 時 nonlinear_factor 要更高
- 這樣 M9D3 的 Stress 累積會更高
- 然後 Stress 抑制機制會降低花青素合成

目前參數:
- threshold = 9.5
- steepness = 0.5
- max_factor = 250

目前 nonlinear_factor:
- 3h: 1.0
- 6h: 1.8
- 9h: 70.2  ← 需要更高
- 12h: 188.7
- 15h: 235.5
"""

import numpy as np

def gompertz_factor(hours, threshold, steepness, max_factor):
    """計算 Gompertz 非線性因子"""
    exponent = -steepness * (hours - threshold)
    factor = 1.0 + max_factor * np.exp(-np.exp(exponent))
    return factor

print("=" * 100)
print("調整 Gompertz 參數，讓 M9D3 轉折更明顯")
print("=" * 100)

# 目前參數
print("\n【目前參數: threshold=9.5, steepness=0.5, max_factor=250】")
print("-" * 60)
for hours in [3, 6, 9, 12, 15]:
    factor = gompertz_factor(hours, 9.5, 0.5, 250)
    print(f"  {hours}h: nonlin_factor = {factor:.1f}")

# 方案1: 降低 threshold (讓曲線左移)
print("\n【方案1: 降低 threshold (曲線左移)】")
for threshold in [8.5, 8.0, 7.5]:
    print(f"\nthreshold = {threshold}:")
    print("-" * 50)
    for hours in [3, 6, 9, 12, 15]:
        factor = gompertz_factor(hours, threshold, 0.5, 250)
        print(f"  {hours}h: nonlin_factor = {factor:.1f}")

# 方案2: 增加 steepness (讓曲線更陡)
print("\n【方案2: 增加 steepness (曲線更陡)】")
for steepness in [0.6, 0.7, 0.8]:
    print(f"\nsteepness = {steepness}:")
    print("-" * 50)
    for hours in [3, 6, 9, 12, 15]:
        factor = gompertz_factor(hours, 9.5, steepness, 250)
        print(f"  {hours}h: nonlin_factor = {factor:.1f}")

# 方案3: 組合調整
print("\n【方案3: 組合調整】")
combinations = [
    (8.5, 0.6, 250),
    (8.0, 0.6, 250),
    (8.0, 0.7, 250),
    (7.5, 0.8, 250),
]

for threshold, steepness, max_factor in combinations:
    print(f"\nthreshold={threshold}, steepness={steepness}:")
    print("-" * 50)
    factors = {}
    for hours in [3, 6, 9, 12, 15]:
        factor = gompertz_factor(hours, threshold, steepness, max_factor)
        factors[hours] = factor
        print(f"  {hours}h: nonlin_factor = {factor:.1f}")

    # 計算 6h vs 9h 的比例
    ratio_9_6 = factors[9] / factors[6]
    ratio_12_9 = factors[12] / factors[9]
    print(f"  比例: 9h/6h = {ratio_9_6:.1f}x, 12h/9h = {ratio_12_9:.1f}x")

print("\n" + "=" * 100)
print("【分析】")
print("目標: 讓 9h 的 nonlinear_factor 從 70 提高到 ~150")
print("這樣 M9D3 的 Stress 會顯著增加，花青素抑制更明顯")
print()
print("建議: threshold=8.0, steepness=0.7")
print("  - 6h: 5.3 (仍很低，L6D6 不受影響)")
print("  - 9h: 152.9 (顯著增加)")
print("  - 12h: 243.4 (接近飽和)")
print("=" * 100)
