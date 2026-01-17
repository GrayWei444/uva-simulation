# Lettuce UVA Model v10.39

**A Mechanistic Model for UVA Effects on Lettuce Growth and Anthocyanin Accumulation**

[![Python Version](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Status](https://img.shields.io/badge/status-Production%20Ready-brightgreen.svg)]()

---

## About

This repository provides a reproducible pipeline and mechanistic modeling code to analyze and optimize UVA lighting recipes for red-leaf lettuce grown in controlled-environment agriculture. The goal is to **enhance anthocyanin accumulation without yield loss** (fresh weight).

The workflow consists of:

1. **Screening stage**: comparing UVA vs. UVB across multiple cultivars to identify a "safe band" and a responsive cultivar.

2. **Optimization stage**: exploring UVA photoperiod × duration combinations under a fixed irradiance and extrapolating responses via modeling and global search.

At the core is a **six-state ODE mechanistic model** integrating UVA-driven morphogenesis, ROS production/clearance, damage–repair dynamics, ontogeny-dependent vulnerability, nonlinear amplification under long exposure, and circadian effects. Parameters are calibrated in a data-driven manner to reproduce endpoint fresh weight and anthocyanin outcomes and to capture non-monotonic (hormetic) responses.

---

## What You Can Do With This Repo

- **Reproduce training/validation results** and summary tables
- **Generate response surfaces / heatmaps** for UVA recipes
- **Run constrained optimization** (e.g., FW ≥ CK − 5%) to identify candidate "best" recipes
- **Adapt the framework** to other cultivars or environments (with re-calibration)

---

## Model Performance (v10.39)

**Perfect Score**: 12/12 targets achieved (Training 6/6 + Validation 6/6)

### Training Set (Tolerance: 5%)

| Treatment | FW Pred | FW Obs | FW Error | Anth Pred | Anth Obs | Anth Error |
|-----------|---------|--------|----------|-----------|----------|------------|
| CK | 86.5g | 87.0g | -0.5% | 439 | 433 | +1.3% |
| L6D6 | 92.5g | 91.4g | +1.2% | 474 | 494 | -4.0% |
| L6D6-N | 84.0g | 80.8g | +3.9% | 475 | 493 | -3.6% |
| VL3D12 | 69.4g | 67.0g | +3.6% | 492 | 482 | +2.0% |
| L6D12 | 58.9g | 60.4g | -2.5% | 496 | 518 | -4.3% |
| H12D3 | 61.3g | 60.6g | +1.2% | 651 | 651 | +0.0% |

### Validation Set (Tolerance: 10%)

| Treatment | Hours | FW Pred | FW Obs | FW Error | Anth Pred | Anth Obs | Anth Error |
|-----------|-------|---------|--------|----------|-----------|----------|------------|
| CK | 0h | 86.5g | 85.2g | +1.6% | 439 | 413 | +6.2% |
| VL3D3 | 3h | 88.4g | 89.0g | -0.8% | 457 | 437 | +4.5% |
| L6D3 | 6h | 89.9g | 92.2g | -2.5% | 473 | 468 | +1.1% |
| M9D3 | 9h | 87.8g | 83.8g | +4.8% | 589 | 539 | +9.2% |
| H12D3 | 12h | 61.3g | 62.2g | -1.4% | 651 | 657 | -0.9% |
| VH15D3 | 15h | 51.2g | 51.3g | +0.0% | 532 | 578 | -7.9% |

---

## Model Architecture

### State Variables

Six-state ODE system: `[X_d, C_buf, LAI, Anth, Stress, ROS]`

| Variable | Description | Unit |
|----------|-------------|------|
| X_d | Dry weight biomass | kg/m² |
| C_buf | Carbon buffer | kg/m² |
| LAI | Leaf Area Index | m²/m² |
| Anth | Anthocyanin content | kg/m² |
| Stress | Cumulative stress | - |
| ROS | Reactive oxygen species | - |

### Base Model

- **Sun et al. (2025)** lettuce growth model
- **ODE Solver**: RK45, max_step=300s

### UVA Effect Mechanisms

1. **UVA Morphological Effect** - UVA promotes SLA and LAI growth (does not directly add to PAR)
2. **ROS Dynamics** - UVA generates ROS, antioxidant system clears it
3. **Stress Damage-Decay** - Balance between cumulative damage and natural decay
4. **LAI Vulnerability** - Young plants are more susceptible to damage
5. **Gompertz Nonlinearity** - Prolonged irradiation triggers antioxidant system collapse
6. **Circadian Damage** - Nighttime irradiation causes additional stress
7. **Anthocyanin Induction** - Stress-induced + UV direct induction
8. **Water Inhibition** - Anthocyanin synthesis efficiency decreases under extreme stress
9. **Hill Efficiency Inhibition** - Monotonically decreasing inhibition at high nonlinear_factor

---

## Core Parameters (v10.39)

### Gompertz Nonlinear Factor

```
nonlinear_factor = 1 + max_factor × exp(-exp(-steepness × (hours - threshold)))
```

| Parameter | Value |
|-----------|-------|
| threshold | 10.5 hours |
| max_factor | 250.0 |
| steepness | 0.5 |

| Hours/Day | nonlinear_factor |
|-----------|------------------|
| 3h | 1.0 |
| 6h | 1.0 |
| 9h | 31.1 |
| 12h | 156.9 |
| 15h | 226.0 |

### Anthocyanin Efficiency Inhibition (Hill Function)

```
efficiency = 1 / (1 + (nonlinear_factor / K)^n)
```

| Parameter | Value |
|-----------|-------|
| K | 800.0 |
| n | 1.5 |

---

## File Structure

```
.
├── simulate_uva_model_v10.py          # Main model (v10.39)
├── lettuce_uva_carbon_complete_model.py  # Sun model base
├── model_config.py                    # Treatment configurations and targets
├── generate_paper_figures.py          # Paper figure generation
├── generate_fig20_quick.py            # Quick optimization heatmap
├── optimize_uva_strategy.py           # Optimization script (11 W/m²)
├── CLAUDE.md                          # Development guidelines (v3.0)
├── HANDOFF_STATUS.md                  # Handoff status
└── MODEL_DESIGN_NOTES.md              # Model design notes
```

---

## Quick Start

### Installation

```bash
pip install numpy scipy matplotlib
```

### Run Simulation

```bash
python3 simulate_uva_model_v10.py
```

### Generate Paper Figures

```bash
python3 generate_paper_figures.py
```

### Run Optimization

```bash
python3 optimize_uva_strategy.py
```

---

## Version History

- **v10.39**: Monotonically decreasing Hill efficiency function (K=800, n=1.5)
- **v10.37**: Gompertz threshold 9.5→10.5
- **v10.33**: Continuous asymmetric Gaussian (removed)
- **v10.9**: Water inhibition mechanism
- **v10.0**: UVA morphological effect replaces direct PAR addition

---

## Citation

If you use this model in your research, please cite the associated manuscript and the repository:

> Wei, C.H., Fang, W., & Huang, C.K. (2026). A Two-Stage Screening-to-Optimization Approach with Mechanistic Model Analysis: Enhancing Anthocyanin in Lettuce Without Yield Loss. *Plants* (under review).

Repository: https://github.com/GrayWei444/uva-simulation

---

## License

This project is for academic research purposes.
