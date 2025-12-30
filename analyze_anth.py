"""分析花青素問題"""
import numpy as np

# 目標值
ANTH_TARGETS = {
    'CK': 43.3,      # 無 UVA
    'L6D6': 49.4,    # 6h×6天
    'L6D6-N': 49.3,  # 夜間 6h×6天
    'VL3D12': 48.2,  # 3h×12天
    'L6D12': 51.8,   # 6h×12天
    'H12D3': 65.1,   # 12h×3天
}

# 目前模擬結果 (v7.0)
ANTH_SIM = {
    'CK': 47.0,
    'L6D6': 49.5,
    'L6D6-N': 52.6,
    'VL3D12': 52.7,
    'L6D12': 53.0,
    'H12D3': 63.8,
}

# Stress 值 (v7.0)
STRESS = {
    'CK': 0.0,
    'L6D6': 0.29,
    'L6D6-N': 1.85,
    'VL3D12': 1.69,
    'L6D12': 6.73,
    'H12D3': 36.57,
}

print("=" * 70)
print("花青素分析")
print("=" * 70)

print("\n花青素 = 基礎合成 + Stress誘導合成")
print("synthesis = FW × (base_rate + V_max × Stress / (K + Stress))")
print()

print(f"{'處理組':<10} {'目標':>6} {'模擬':>6} {'誤差':>8} {'Stress':>8}")
print("-" * 50)
for t in ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']:
    target = ANTH_TARGETS[t]
    sim = ANTH_SIM[t]
    err = (sim - target) / target * 100
    stress = STRESS[t]
    s = "✓" if abs(err) < 5 else "✗"
    print(f"{t:<10} {target:>6.1f} {sim:>6.1f} {err:>+6.1f}%{s} {stress:>8.2f}")

print("\n問題分析:")
print("1. CK: +8.5% → 基礎合成率過高")
print("2. VL3D12: +9.4% → Stress 誘導合成過高 (Stress=1.69)")
print("3. H12D3: -1.9% → 接近目標")
print()

# 計算需要的調整
print("調整方向:")
print("- 降低 base_rate: 會降低所有組")
print("- 但 H12D3 已經 -1.9%，降低 base 會讓它更低")
print()

# 分析 Stress 誘導的相對貢獻
print("Stress 誘導貢獻分析 (假設 K=0.3):")
K = 0.30
for t in ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']:
    S = STRESS[t]
    induced = S / (K + S) if S > 0 else 0
    print(f"  {t:<10}: Stress={S:>5.2f}, 誘導因子={induced:.3f}")

print()
print("結論:")
print("- CK (Stress=0) 花青素完全來自基礎合成")
print("- 如果降低 base_rate 讓 CK 達標，需要提高 V_max 讓 H12D3 達標")
print("- 但這會讓中等 Stress 的組 (VL3D12, L6D12) 誤差更大")
print()

# 計算需要的 base_rate
print("計算需要的 base_rate (讓 CK 達標):")
# CK 只有基礎合成，Anth ∝ base_rate
# 目前 base=2e-10 → Anth=47.0
# 目標 Anth=43.3
# 需要 base = 2e-10 × (43.3/47.0) = 1.84e-10
new_base = 2e-10 * (43.3 / 47.0)
print(f"  需要 base_rate = {new_base:.2e}")
print()

# 如果用這個 base_rate，H12D3 會變成多少？
# H12D3 目前 Anth=63.8，其中基礎合成部分會降低
# 降低比例 = 43.3/47.0 = 0.921
print("影響分析 (假設基礎合成佔比):")
for t in ['CK', 'L6D6', 'L6D6-N', 'VL3D12', 'L6D12', 'H12D3']:
    target = ANTH_TARGETS[t]
    sim = ANTH_SIM[t]
    # 估計基礎部分 (假設 CK 全是基礎)
    base_part = 47.0  # CK 的值作為基礎
    induced_part = max(0, sim - base_part * 0.8)  # 粗估
    print(f"  {t:<10}: 目前={sim:.1f}, 目標={target:.1f}")
