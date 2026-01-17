# Lettuce UVA Model

**A Mechanistic Model for UVA Effects on Lettuce Growth and Anthocyanin Accumulation**

[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

---

## About

This repository provides reproducible code for a mechanistic model that simulates UVA lighting effects on red-leaf lettuce (*Lactuca sativa* L.) grown in controlled-environment agriculture. The model predicts:

- **Fresh weight (FW)** - Biomass accumulation under various UVA treatments
- **Anthocyanin concentration** - Secondary metabolite induction by UVA stress

The goal is to **enhance anthocyanin accumulation without yield loss**.

---

## Model Architecture

### Six-State ODE System

| State Variable | Symbol | Description | Unit |
|----------------|--------|-------------|------|
| Dry weight | X_d | Structural biomass | kg/m² |
| Carbon buffer | C_buf | Non-structural carbohydrates | kg/m² |
| Leaf area index | LAI | Canopy light interception | m²/m² |
| Anthocyanin | Anth | Secondary metabolite content | kg/m² |
| Stress | Stress | Cumulative oxidative damage | - |
| ROS | ROS | Reactive oxygen species | - |

### Key Mechanisms

1. **Base growth model** - Sun et al. (2025) lettuce carbon allocation model
2. **UVA morphological effect** - Enhanced SLA and LAI development
3. **ROS dynamics** - Production under UVA, clearance by antioxidant system
4. **Stress accumulation** - Damage-decay balance with circadian modulation
5. **Gompertz nonlinearity** - Antioxidant collapse under prolonged exposure
6. **Anthocyanin induction** - Stress-triggered biosynthesis pathway

---

## Model Performance

### Training Set (n=6, Tolerance: 5%)

| Treatment | FW Obs (g) | FW Pred (g) | FW Error | Anth Obs | Anth Pred | Anth Error |
|-----------|------------|-------------|----------|----------|-----------|------------|
| CK (Control) | 87.0 | 86.5 | -0.5% | 433 | 439 | +1.3% |
| L6D6 (6h×6d) | 91.4 | 92.5 | +1.2% | 494 | 474 | -4.0% |
| L6D6-N (Night) | 80.8 | 84.0 | +3.9% | 493 | 475 | -3.6% |
| VL3D12 (3h×12d) | 67.0 | 69.4 | +3.6% | 482 | 492 | +2.0% |
| L6D12 (6h×12d) | 60.4 | 58.9 | -2.5% | 518 | 496 | -4.3% |
| H12D3 (12h×3d) | 60.6 | 61.3 | +1.2% | 651 | 651 | +0.0% |

### Validation Set (n=6, Tolerance: 10%)

| Treatment | Hours/day | FW Obs (g) | FW Pred (g) | FW Error | Anth Obs | Anth Pred | Anth Error |
|-----------|-----------|------------|-------------|----------|----------|-----------|------------|
| CK | 0h | 85.2 | 86.5 | +1.6% | 413 | 439 | +6.2% |
| VL3D3 | 3h | 89.0 | 88.4 | -0.8% | 437 | 457 | +4.5% |
| L6D3 | 6h | 92.2 | 89.9 | -2.5% | 468 | 473 | +1.1% |
| M9D3 | 9h | 83.8 | 87.8 | +4.8% | 539 | 589 | +9.2% |
| H12D3 | 12h | 62.2 | 61.3 | -1.4% | 657 | 651 | -0.9% |
| VH15D3 | 15h | 51.3 | 51.2 | +0.0% | 578 | 532 | -7.9% |

**Result: 12/12 targets achieved** (Training 6/6 + Validation 6/6)

---

## Installation

```bash
# Clone the repository
git clone https://github.com/GrayWei444/uva-simulation.git
cd uva-simulation

# Install dependencies
pip install -r requirements.txt
```

### Requirements

- Python 3.8+
- NumPy
- SciPy
- Matplotlib

---

## Usage

### Run Model Validation

```bash
python simulate_uva_model_v10.py
```

This runs all training and validation treatments and outputs predicted vs. observed values.

### Generate Paper Figures

```bash
python generate_paper_figures.py
```

Generates all figures used in the manuscript.

### Run UVA Optimization

```bash
python optimize_uva_strategy.py
```

Searches for optimal UVA treatment protocols that maximize anthocyanin while maintaining yield.

---

## File Structure

```
.
├── simulate_uva_model_v10.py          # Main model (6-state ODE system)
├── lettuce_uva_carbon_complete_model.py  # Base Sun model
├── model_config.py                    # Treatment configurations & targets
├── generate_paper_figures.py          # Figure generation script
├── generate_fig20_quick.py            # Quick optimization heatmap
├── optimize_uva_strategy.py           # Optimization search
├── requirements.txt                   # Python dependencies
└── README.md                          # This file
```

---

## Key Parameters

### Gompertz Nonlinear Damage Factor

```
nonlinear_factor = 1 + 250 × exp(-exp(-0.5 × (hours - 10.5)))
```

| Daily UVA (h) | Nonlinear Factor | Description |
|---------------|------------------|-------------|
| 3 | 1.0 | No amplification |
| 6 | 1.0 | No amplification |
| 9 | 31.1 | Transition zone |
| 12 | 156.9 | Severe stress |
| 15 | 226.0 | Near saturation |

### Hill Efficiency Inhibition

```
efficiency = 1 / (1 + (nonlinear_factor / 800)^1.5)
```

---

## Citation

If you use this code in your research, please cite:

> Wei, C.H., Fang, W., & Huang, C.K. (2026). A Two-Stage Screening-to-Optimization Approach with Mechanistic Model Analysis: Enhancing Anthocyanin in Lettuce Without Yield Loss. *Plants* (under review).

---

## License

MIT License - See [LICENSE](LICENSE) for details.

---

## Contact

For questions or collaborations, please open an issue on this repository.
