# Lettuce UVA Model

**A Mechanistic Model for UVA Effects on Lettuce Growth and Anthocyanin Accumulation**

[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

---

## About

This repository provides reproducible code for a mechanistic model that simulates UVA lighting effects on red-leaf lettuce (*Lactuca sativa* L.) grown in controlled-environment agriculture. The model predicts:

- **Fresh weight (FW)** - Biomass accumulation under various UVA treatments
- **Anthocyanin concentration** - Secondary metabolite induction by UVA stress

---

## Model Architecture

### Six-State ODE System

| State | Symbol | Description | Unit |
|-------|--------|-------------|------|
| Dry weight | X_d | Structural biomass | kg/m² |
| Carbon buffer | C_buf | Non-structural carbohydrates | kg/m² |
| Leaf area index | LAI | Canopy light interception | m²/m² |
| Anthocyanin | Anth | Secondary metabolite | kg/m² |
| Stress | Stress | Cumulative oxidative damage | - |
| ROS | ROS | Reactive oxygen species | - |

### Key Mechanisms

1. **Base growth model** - Sun et al. (2025) lettuce carbon allocation
2. **UVA morphological effect** - Enhanced SLA and LAI
3. **ROS dynamics** - Production and clearance balance
4. **Gompertz nonlinearity** - Antioxidant collapse under prolonged exposure
5. **Anthocyanin induction** - Stress-triggered biosynthesis

---

## Model Performance

### Training Set (n=6, Tolerance: 5%)

| Treatment | FW Obs | FW Pred | FW Error | Anth Obs | Anth Pred | Anth Error |
|-----------|--------|---------|----------|----------|-----------|------------|
| CK | 87.0g | 86.5g | -0.5% | 433 | 439 | +1.3% |
| L6D6 | 91.4g | 92.5g | +1.2% | 494 | 474 | -4.0% |
| L6D6-N | 80.8g | 84.0g | +3.9% | 493 | 475 | -3.6% |
| VL3D12 | 67.0g | 69.4g | +3.6% | 482 | 492 | +2.0% |
| L6D12 | 60.4g | 58.9g | -2.5% | 518 | 496 | -4.3% |
| H12D3 | 60.6g | 61.3g | +1.2% | 651 | 651 | +0.0% |

### Validation Set (n=6, Tolerance: 10%)

| Treatment | Hours | FW Obs | FW Pred | FW Error | Anth Obs | Anth Pred | Anth Error |
|-----------|-------|--------|---------|----------|----------|-----------|------------|
| CK | 0h | 85.2g | 86.5g | +1.6% | 413 | 439 | +6.2% |
| VL3D3 | 3h | 89.0g | 88.4g | -0.8% | 437 | 457 | +4.5% |
| L6D3 | 6h | 92.2g | 89.9g | -2.5% | 468 | 473 | +1.1% |
| M9D3 | 9h | 83.8g | 87.8g | +4.8% | 539 | 589 | +9.2% |
| H12D3 | 12h | 62.2g | 61.3g | -1.4% | 657 | 651 | -0.9% |
| VH15D3 | 15h | 51.3g | 51.2g | +0.0% | 578 | 532 | -7.9% |

**Result: 12/12 targets achieved** (Training 6/6 <5%, Validation 6/6 <10%)

---

## Installation

```bash
git clone https://github.com/GrayWei444/uva-simulation.git
cd uva-simulation
pip install -r requirements.txt
```

### Requirements

- Python 3.8+
- NumPy >= 1.21.0
- SciPy >= 1.7.0
- Matplotlib >= 3.4.0

---

## Usage

```bash
python simulate_uva_model_v10.py
```

This runs all training and validation simulations and outputs predicted vs. observed values.

---

## File Structure

```
.
├── simulate_uva_model_v10.py          # Main model (standalone)
├── lettuce_uva_carbon_complete_model.py  # Base Sun model
├── requirements.txt                   # Dependencies
├── LICENSE                            # MIT License
└── README.md                          # This file
```

---

## Key Parameters

### Gompertz Nonlinear Factor

```
nonlinear_factor = 1 + 250 × exp(-exp(-0.5 × (hours - 10.5)))
```

| Daily UVA | Factor | Description |
|-----------|--------|-------------|
| 3h | 1.0 | No amplification |
| 6h | 1.0 | No amplification |
| 9h | 31.1 | Transition zone |
| 12h | 156.9 | Severe stress |
| 15h | 226.0 | Near saturation |

---

## Citation

> Wei, C.H., Fang, W., & Huang, C.K. (2026). A Two-Stage Screening-to-Optimization Approach with Mechanistic Model Analysis: Enhancing Anthocyanin in Lettuce Without Yield Loss. *Plants* (under review).

---

## License

MIT License - See [LICENSE](LICENSE) for details.
