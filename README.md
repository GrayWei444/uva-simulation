# Lettuce UVA Model v2.0

Mechanistic UVA–lettuce model with explicit carbon competition between growth and antioxidant defense.

## Contents

- `lettuce_uva_model.py` — main v2.0 model (AOX + carbon competition)
- `lettuce_uva_carbon_complete_model.py` — base Sun model dependency
- `simulate_uva_model_v2.py` — runs training + validation simulations and prints results
- `parameters_v2.md` — v2.0 parameter documentation

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python simulate_uva_model_v2.py
```

The script prints training/validation metrics to stdout.

## Notes

- Anthocyanin is computed as `Anth = AOX × 0.18`.
- Carbon competition reduces growth and slightly feeds back on AOX synthesis.

## License

MIT License — see `LICENSE`.
