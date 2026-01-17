# Claude Code Guidelines v3.0

**Core Principle: Reproducibility & Consistency**

---

## 0. Most Important Principle: Reproducibility Requirements

### Strictly Prohibited:
1. **No result manipulation** - Do not arbitrarily modify parameters to match target values
2. **No hallucinations** - Documentation must be 100% consistent with code
3. **No hard thresholds** - All mechanisms must use continuous functions

### Reproducibility Verification Process:
After modifying any parameters or formulas, execute:
```bash
python3 simulate_uva_model_v10.py
```
Record the **complete output** in HANDOFF_STATUS.md

### Document Synchronization Rules:
1. Formulas in code = Formulas in HANDOFF_STATUS.md = Formulas in MODEL_DESIGN_NOTES.md
2. Parameter values in code = Parameter values in documentation
3. Simulation output values = Values recorded in documentation

**Violating reproducibility = Invalid modification**

---

## 1. When Starting a New Chat

**Must execute in order:**
1. Read `CLAUDE.md` - Understand guidelines
2. Read `HANDOFF_STATUS.md` - Understand current version, parameters, simulation results
3. Read `MODEL_DESIGN_NOTES.md` - Understand model mechanisms
4. **Run verification once**: `python3 simulate_uva_model_v10.py`
5. **Compare output** with documented records for consistency

If inconsistent, documentation must be corrected before other work.

---

## 2. Current Version Status (v10.39)

### Core Formulas

**Gompertz Nonlinear Factor:**
```python
# threshold=10.5, max_factor=250, steepness=0.5
nonlinear_factor = 1 + max_factor * np.exp(-np.exp(-steepness * (hours - threshold)))
```

**Anthocyanin Synthesis Efficiency Inhibition (Hill Function):**
```python
# K=800, n=1.5
efficiency = 1.0 / (1.0 + (nonlinear_factor / K) ** n)
```

### Efficiency Reference Table
| Hours/Day | nonlinear_factor | Efficiency |
|-----------|------------------|------------|
| 3h | 1.0 | 100.0% |
| 6h | 1.0 | 100.0% |
| 9h | 31.1 | 99.2% |
| 12h | 156.9 | 92.0% |
| 15h | 226.0 | 86.9% |

### Actual Simulation Results (Must match documentation)

**Training Set:**
| Treatment | FW Obs | FW Pred | FW Error | Anth Obs | Anth Pred | Anth Error |
|-----------|--------|---------|----------|----------|-----------|------------|
| CK | 87.0 | 86.5 | -0.5% | 433 | 439 | +1.3% |
| L6D6 | 91.4 | 92.5 | +1.2% | 494 | 474 | -4.0% |
| L6D6-N | 80.8 | 84.0 | +3.9% | 493 | 475 | -3.6% |
| VL3D12 | 67.0 | 69.4 | +3.6% | 482 | 492 | +2.0% |
| L6D12 | 60.4 | 58.9 | -2.5% | 518 | 496 | -4.3% |
| H12D3 | 60.6 | 61.3 | +1.2% | 651 | 651 | +0.0% |

**Validation Set:**
| Treatment | FW Obs | FW Pred | FW Error | Anth Obs | Anth Pred | Anth Error |
|-----------|--------|---------|----------|----------|-----------|------------|
| CK | 85.2 | 86.5 | +1.6% | 413 | 439 | +6.2% |
| VL3D3 | 89.0 | 88.4 | -0.8% | 437 | 457 | +4.5% |
| L6D3 | 92.2 | 89.9 | -2.5% | 468 | 473 | +1.1% |
| M9D3 | 83.8 | 87.8 | +4.8% | 539 | 589 | +9.2% |
| H12D3 | 62.2 | 61.3 | -1.4% | 657 | 651 | -0.9% |
| VH15D3 | 51.3 | 51.2 | +0.0% | 578 | 532 | -7.9% |

---

## 3. No Hard Thresholds

### Prohibited patterns:
```python
if hours > 6:
    damage = high_damage

circ = 3.8 if is_night else 1.0
```

### Correct patterns:
```python
# Continuous Hill function
efficiency = 1.0 / (1.0 + (x / K) ** n)

# Continuous Gompertz
factor = 1 + max_factor * exp(-exp(-k * (hours - threshold)))
```

---

## 4. When Token Limit Approaches or Ending Chat

**Must execute:**
1. Run `python3 simulate_uva_model_v10.py` to get complete output
2. Update `HANDOFF_STATUS.md` including:
   - Current version number
   - All formulas (exactly matching code)
   - All parameter values
   - Complete simulation results table
3. Synchronize `MODEL_DESIGN_NOTES.md`

---

## 5. File Structure

| File | Purpose |
|------|---------|
| `simulate_uva_model_v10.py` | Main model (v10.39) |
| `model_config.py` | Treatment configurations and target values |
| `CLAUDE.md` | Development guidelines (this file) |
| `HANDOFF_STATUS.md` | Version status, formulas, parameters, results |
| `MODEL_DESIGN_NOTES.md` | Model mechanism documentation |

---

## 6. Strict Compliance

1. **Reproducibility first** - All modifications must be verifiable and reproducible
2. **Document synchronization** - Code and documentation must be 100% consistent
3. **No hallucinations** - Do not fabricate data or formulas
4. **Verify first** - Must run verification after modifications
