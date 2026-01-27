# Lettuce UVA Model v2.0

**A Mechanistic Model for UVA Effects on Lettuce Growth and Anthocyanin Accumulation with Carbon Competition**

[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

---

## About

This repository provides reproducible code for a mechanistic model that simulates UVA lighting effects on red-leaf lettuce (*Lactuca sativa* L.) grown in controlled-environment agriculture. The model predicts:

- **Fresh weight (FW)** - Biomass accumulation under various UVA treatments
- **Anthocyanin concentration** - Secondary metabolite induction by UVA stress

### v2.0 Key Innovation: Carbon Competition Mechanism

The v2.0 model introduces explicit **carbon competition between growth and defense**, implementing the Growth-Differentiation Balance Hypothesis (Herms & Mattson, 1992):

- AOX (antioxidant) synthesis consumes carbon from the buffer pool (C_buf)
- Creates a physiologically meaningful growth–defense tradeoff
- Supported by literature on carbon costs of phenylpropanoid biosynthesis (Vogt, 2010)

---

## Model Architecture

### Six-State ODE System

| State | Symbol | Description | Unit |
|-------|--------|-------------|------|
| Dry weight | X_d | Structural biomass | kg/m² |
| Carbon buffer | C_buf | Non-structural carbohydrates | kg/m² |
| Leaf area index | LAI | Canopy light interception | m²/m² |
| Antioxidants | AOX | Total antioxidants (Anth = AOX × 18%) | kg/m² |
| Stress | Stress | Cumulative oxidative damage | - |
| ROS | ROS | Reactive oxygen species | - |

### Key Mechanisms

1. **Base growth model** - Sun et al. (2025) lettuce carbon allocation
2. **UVA morphological effect** - Enhanced SLA and LAI
3. **ROS dynamics** - Production and clearance balance
4. **Gompertz nonlinearity** - Antioxidant collapse under prolonged exposure
5. **Carbon competition** - AOX synthesis consumes C_buf, creating growth–defense tradeoff
6. **Anthocyanin induction** - Stress-triggered biosynthesis (Anth = AOX × 0.18)

---

## Model Performance

### Training Set (n=6, Tolerance: 5%)

| Treatment | FW Obs | FW Pred | FW Error | Anth Obs | Anth Pred | Anth Error |
|-----------|--------|---------|----------|----------|-----------|------------|
| CK | 87.0g | 87.9g | +1.0% | 433 | 440 | +1.6% |
| L6D6 | 91.4g | 90.2g | -1.3% | 494 | 487 | -1.4% |
| L6D6-N | 80.8g | 83.5g | +3.4% | 493 | 481 | -2.5% |
| VL3D12 | 67.0g | 70.7g | +5.6% | 482 | 509 | +5.6% |
| L6D12 | 60.4g | 61.5g | +1.8% | 518 | 533 | +2.9% |
| H12D3 | 60.6g | 60.9g | +0.6% | 651 | 655 | +0.6% |

### Validation Set (n=6, Tolerance: 10%)

| Treatment | Hours | FW Obs | FW Pred | FW Error | Anth Obs | Anth Pred | Anth Error |
|-----------|-------|--------|---------|----------|----------|-----------|------------|
| CK | 0h | 85.2g | 87.9g | +3.2% | 413 | 440 | +6.5% |
| VL3D3 | 3h | 89.0g | 88.5g | -0.7% | 437 | 460 | +5.4% |
| L6D3 | 6h | 92.2g | 89.1g | -3.4% | 468 | 480 | +2.5% |
| M9D3 | 9h | 83.8g | 85.9g | +2.5% | 539 | 596 | +10.6% |
| H12D3 | 12h | 62.2g | 60.9g | -2.0% | 657 | 655 | -0.3% |
| VH15D3 | 15h | 51.3g | 51.2g | +0.0% | 578 | 605 | +4.8% |

**Result: 11/12 targets achieved** (Training 10/12 <5%, Validation 11/12 <10%)

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

### Run Simulation

```bash
python simulate_uva_model_v2.py
```

This runs all training and validation simulations and exports results to `results.csv`.

### Generate Paper Figures

```bash
python generate_paper_figures.py
```

Generates all figures (Fig 9-21) in the `paper_figures/` directory.

---

## File Structure

```
.
├── simulate_uva_model_v2.py          # Main model v2.0 (with carbon competition)
├── lettuce_uva_carbon_complete_model.py  # Base Sun model
├── generate_paper_figures.py         # Figure generation script
├── requirements.txt                   # Dependencies
├── LICENSE                            # MIT License
├── README.md                          # This file
├── paper_figures/                     # Generated figures (not tracked)
│   ├── Fig9_LAI_vulnerability.png
│   ├── Fig10_Gompertz_nonlinear.png
│   ├── Fig11_training_parity.png
│   ├── ...
│   └── Fig21_carbon_competition.png
└── results.csv                        # Simulation results
```

---

## Key Parameters

### Carbon Competition (v2.0 New)

| Parameter | Value | Description |
|-----------|-------|-------------|
| aox_carbon_cost | 1.0 kg C/kg AOX | Carbon cost of AOX synthesis |
| carbon_competition_max | 0.30 | Maximum growth penalty from AOX synthesis |
| stress_competition_K | 21.0 | Half-saturation for stress-based competition |
| stress_competition_max | 0.225 | Maximum stress-based growth penalty |

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

## Literature Support for Carbon Competition

1. **Herms & Mattson (1992)** - Growth-Differentiation Balance Hypothesis
   - DOI: [10.1086/285343](https://doi.org/10.1086/285343)

2. **Monson et al. (2022)** - Coordinated resource allocation to growth–defense tradeoffs
   - DOI: [10.1111/nph.17773](https://doi.org/10.1111/nph.17773)

3. **Vogt (2010)** - Phenylpropanoid Biosynthesis (~20% photosynthate to phenylpropanoids)
   - DOI: [10.1093/mp/ssp106](https://doi.org/10.1093/mp/ssp106)

4. **Gershenzon (1994)** - Metabolic costs of terpenoid accumulation
   - DOI: [10.1007/BF02059810](https://doi.org/10.1007/BF02059810)

---

## Citation

> Wei, C.H., Fang, W., & Huang, C.K. (2026). A Two-Stage Screening-to-Optimization Approach with Mechanistic Model Analysis: Enhancing Anthocyanin in Lettuce Without Yield Loss. *Plants* (under review).

---

## License

MIT License - See [LICENSE](LICENSE) for details.
