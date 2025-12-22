# v5.9 使用指南

## 快速開始

### 1. 測試新機制

```bash
python3 test_anthocyanin_v59.py
```

這個腳本會：
- 顯示所有處理組的預測結果
- 顯示 E_stress 累積值
- 測試參數敏感度
- 比較不同參數設定的效果

### 2. 校準參數

```bash
python3 calibrate_anthocyanin_v59.py
```

這個腳本會：
- 使用差分進化法自動優化花青素參數
- 給予 VL3D12 和 L6D12 更高權重
- 顯示最佳參數和驗證結果

### 3. 單一處理組模擬

```bash
python3 simulate_uva_model.py
```

預設模擬 L6D6 處理組，可以修改腳本中的 `treatment` 變量。

## 參數調整指南

### 花青素參數位置

在 `simulate_uva_model.py` 的 `UVAParams` 類別中（第 166-171 行）：

```python
self.base_anth_rate_light = 4.11e-10   # 日間基礎合成率
self.base_anth_rate_dark = 2.05e-10    # 夜間基礎合成率
self.V_max_anth = 3.0e-10              # 最大誘導合成率
self.K_E_anth = 50.0                   # 半飽和 Stress 累積能量
self.n_hill_anth = 2.0                 # Hill 係數
self.k_deg = 2.5e-6                    # 降解率
```

### 參數調整建議

#### 1. V_max_anth (最大誘導合成率)
- **增加** → 所有 UVA 處理組的花青素都會增加
- **減少** → 所有 UVA 處理組的花青素都會減少
- **建議範圍**: 1e-10 ~ 1e-8 kg/m²/s

#### 2. K_E_anth (半飽和常數)
- **增加** → 花青素響應變慢，需要更高的 E_stress 才能誘導
- **減少** → 花青素響應變快，低 E_stress 也能誘導
- **建議範圍**: 20 ~ 100 Stress·day

**影響模式**:
```
K_E = 30: 低 E_stress 處理組花青素↑，高 E_stress 處理組花青素略↑
K_E = 60: 低 E_stress 處理組花青素↓，高 E_stress 處理組花青素略↓
```

#### 3. n_hill_anth (Hill 係數)
- **增加** → 響應曲線變陡，低/高 E_stress 之間差異更大
- **減少** → 響應曲線變平滑，低/高 E_stress 差異較小
- **建議範圍**: 1.5 ~ 4.0

**影響模式**:
```
n = 1.5: 較平滑的響應
n = 3.0: 較陡峭的響應（明顯的閾值效應）
```

#### 4. k_deg (降解率)
- **增加** → 所有處理組花青素都會減少
- **減少** → 所有處理組花青素都會增加
- **建議**: 保持 2.5e-6 (半衰期 ~3.2 天)

## 調整策略

### 問題：所有處理組花青素都偏高
**解決**: 增加 `k_deg` 或減少 `V_max_anth`

### 問題：所有處理組花青素都偏低
**解決**: 減少 `k_deg` 或增加 `V_max_anth`

### 問題：VL3D12 和 L6D12 偏高，其他正常
**解決**: 增加 `K_E_anth` (使低 E_stress 響應減弱)

### 問題：H12D3 偏低，其他正常
**解決**: 減少 `K_E_anth` 或增加 `n_hill_anth` (加強高 E_stress 響應)

### 問題：低/高 E_stress 處理組差異不明顯
**解決**: 增加 `n_hill_anth` (使曲線更陡)

## 預期 E_stress 範圍

基於 v5.7 的 Stress 值，粗估的 E_stress：

```
CK:      0 Stress·day
L6D6:    ~10 Stress·day
L6D6-N:  ~20 Stress·day
H12D3:   ~50 Stress·day
VL3D12:  ~40 Stress·day
L6D12:   ~70 Stress·day
```

**注意**: 實際值需要運行模擬才能得到。

## Hill 函數計算器

使用以下公式計算不同 E_stress 下的響應：

```python
def hill_response(E_stress, V_max, K_E, n):
    E_power_n = E_stress ** n
    K_power_n = K_E ** n
    response = V_max * E_power_n / (K_power_n + E_power_n)
    return response

# 範例
E_stress = 40  # Stress·day
V_max = 3.0e-10
K_E = 50.0
n = 2.0

response = hill_response(E_stress, V_max, K_E, n)
# response ≈ 1.7e-10 kg/m²/s
```

## 驗證檢查清單

校準完成後，檢查：

- [ ] FW Mean |Err| < 3%
- [ ] FW Max |Err| < 5%
- [ ] Anth Mean |Err| < 15%
- [ ] Anth Max |Err| < 25%
- [ ] VL3D12 Anth 誤差 < +20%
- [ ] L6D12 Anth 誤差 < +25%
- [ ] CK Anth 接近實驗值 (43.3 ppm)
- [ ] H12D3 Anth 仍然最高 (~65 ppm)

## 常見問題

### Q: E_stress 一直累積，不會降低嗎？
A: E_stress 是累積量，只增不減。這符合「植物記住脅迫」的概念。花青素的飽和效應由 Hill 函數控制。

### Q: 為什麼不用 Stress 本身而要用 E_stress？
A: Stress 會修復（有正有負），但植物對脅迫的「記憶」不會消失。花青素合成是對「累積脅迫歷史」的響應。

### Q: Hill 係數 n 的生物學意義？
A: n > 1 表示協同效應，即累積脅迫需要達到一定程度後，花青素才會快速合成。類似於酶的協同調控。

### Q: 如何選擇 K_E_anth 的初始值？
A: 看實驗中哪個處理組達到「中等花青素誘導」，其 E_stress 約為 K_E_anth 的合理估計。

## 進階使用

### 導出 E_stress 時間序列

在 `simulate_uva_model.py` 中，`sol.y[5, :]` 包含整個模擬期間的 E_stress 軌跡：

```python
import matplotlib.pyplot as plt

# 執行模擬後
time_days = sol.t / 86400
E_stress_trace = sol.y[5, :]

plt.plot(time_days, E_stress_trace)
plt.xlabel('Days from transplant')
plt.ylabel('E_stress [Stress·day]')
plt.title('Stress 累積能量軌跡')
plt.grid(True)
plt.show()
```

### 繪製 Hill 響應曲線

```python
import numpy as np
import matplotlib.pyplot as plt

E_range = np.linspace(0, 100, 200)
p = UVAParams()

responses = []
for E in E_range:
    E_power_n = E ** p.n_hill_anth
    K_power_n = p.K_E_anth ** p.n_hill_anth
    response = p.V_max_anth * E_power_n / (K_power_n + E_power_n)
    responses.append(response)

plt.plot(E_range, responses)
plt.axvline(p.K_E_anth, color='r', linestyle='--', label=f'K_E = {p.K_E_anth}')
plt.xlabel('E_stress [Stress·day]')
plt.ylabel('UVA-induced synthesis rate [kg/m²/s]')
plt.title(f'Hill Response (n={p.n_hill_anth})')
plt.legend()
plt.grid(True)
plt.show()
```

---

**版本**: v5.9
**日期**: 2025-12-22
