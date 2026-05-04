# petrochar

A Python tool for continuous-distribution characterization of heavy petroleum fractions, producing discrete pseudo-components with PC-SAFT parameters ready for use in Aspen Plus or any PC-SAFT-capable simulator.

---

## What it does

Takes routinely-measured refinery data:
- Distillation curve (any ASTM method: D86, D1160 AET, D7169, TBP)
- Bulk specific gravity at 15°C
- Bulk molecular weight
- SARA wt%

And produces:
- A discrete pseudo-component table with mole fractions, MW, Tb, SG, K_W, γ, and PC-SAFT parameters (m, σ, ε/k)
- Aspen-Plus-paste-ready exports
- Validation diagnostics

The methodology follows Riazi (2005) MNL50 Chapter 4 for the continuous distribution and Gaussian-quadrature discretization, Panuganti (2012) for PC-SAFT parameter correlations on distillable pseudo-components, and Gonzalez (2007) for asphaltene defaults.

## Why

Most published PC-SAFT calibration workflows for petroleum require bench-scale n-heptane titration data on dead-oil samples. This is unavailable in operating refineries. petrochar produces a defensible PC-SAFT parameter set from data the refinery routinely produces.

## Status

Under active development. See `CURRENT_STATUS.md` for the latest phase and outstanding items.

## Citation

(Pending publication.)

## License

MIT.
