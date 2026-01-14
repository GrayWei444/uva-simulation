#!/usr/bin/env python3
"""
快速測試:
1. Gompertz threshold 從 9.5 → 10.0
2. 降低 V_max_anth (花青素合成效率)
3. 移除非對稱高斯
"""

import numpy as np

def gompertz_factor(hours, threshold=9.5, steepness=0.5, max_factor=250):
    exponent = -steepness * (hours - threshold)
    return 1.0 + max_factor * np.exp(-np.exp(exponent))

print("=" * 80)
print("Gompertz nonlinear_factor 比較")
print("=" * 80)

print(f"\n{'時數':>6} | {'原始(th=9.5)':>15} | {'新(th=10.0)':>15} | {'新(th=10.5)':>15}")
print("-" * 60)
for hours in [3, 6, 9, 12, 15]:
    orig = gompertz_factor(hours, 9.5)
    new1 = gompertz_factor(hours, 10.0)
    new2 = gompertz_factor(hours, 10.5)
    print(f"{hours:>5}h | {orig:>15.1f} | {new1:>15.1f} | {new2:>15.1f}")

print("\n" + "=" * 80)
print("分析:")
print("- threshold=10.0: 9h 的 nonlin_factor 從 70.2 降到 31.5")
print("- threshold=10.5: 9h 的 nonlin_factor 從 70.2 降到 13.2")
print("- M9D3 的 Stress 會相應降低")
print("=" * 80)

# 預估 Stress 變化 (簡化估算)
# nonlin_damage = k_nonlinear_stress * ROS * nonlinear_factor
# 假設 ROS 相對穩定，Stress 大致與 nonlinear_factor 成正比

print("\n預估 avgStress 變化 (假設與 nonlin_factor 成正比):")
print("-" * 60)
stress_orig = {'M9D3': 68.9, 'H12D3': 405.5, 'VH15D3': 937.6}
for hours, name in [(9, 'M9D3'), (12, 'H12D3'), (15, 'VH15D3')]:
    orig_nf = gompertz_factor(hours, 9.5)
    new_nf = gompertz_factor(hours, 10.0)
    ratio = new_nf / orig_nf
    est_new_stress = stress_orig[name] * ratio
    print(f"{name}: 原Stress={stress_orig[name]:.1f}, 預估新Stress={est_new_stress:.1f} (比例={ratio:.2f})")
