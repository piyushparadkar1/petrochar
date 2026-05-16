# petrochar

**Continuous-distribution characterization of heavy petroleum fractions for PC-SAFT phase-equilibrium modeling from routine refinery data.**

---

## What petrochar does

petrochar takes the distillation curve, bulk specific gravity, bulk molecular weight, and SARA composition of a heavy petroleum fraction and produces a discrete set of pseudo-components with PC-SAFT molecular parameters ready for use in Aspen Plus or any PC-SAFT-capable process simulator.

The methodology follows Riazi (2005) for distribution fitting and quadrature discretization, Panuganti et al. (2012) for PC-SAFT parameter correlations, and Gonzalez et al. (2007) for the discrete asphaltene component. The entire characterization runs from routine refinery assay data — no solvent-deasphalting data, no AOP measurements, and no tuning against process observables are required or used.

The outputs are: a pseudo-component table (z, M, Tb, SG, Watson K, aromaticity γ, m, σ, ε/k), a PC-SAFT parameter export formatted for direct paste into Aspen Plus, and a set of distribution diagnostics (PDF, CDF, fit-quality metrics, SARA closure check). A Streamlit web interface exposes all functionality and includes an in-app validation tab that reproduces Riazi textbook examples live.

---

## Installation

**From source (recommended during development):**

```bash
git clone https://github.com/piyushparadkar/petrochar.git
cd petrochar
pip install -e ".[dev]"
```

The `[dev]` extra installs pytest. The base install is sufficient for the Streamlit app.

**Requirements:** Python ≥ 3.10. Dependencies (installed automatically):
`numpy`, `scipy`, `pandas`, `matplotlib`, `streamlit`, `openpyxl`.

---

## Quick start

```bash
streamlit run app.py
```

1. **Tab 1 — Input Data:** The pre-populated template uses a synthetic VTB-like feed
   (D1160 AET, weight basis, bulk SG = 1.020, bulk MW = 700 g/mol, SARA 12/38/38/12 wt%).
   Click **Run pipeline**.

2. **Tab 4 — Pseudocomponents:** Inspect the 6-row pseudo-component table
   (5 distillable + 1 discrete asphaltene). Download CSV for further use.

3. **Tab 5 — PC-SAFT Export:** Copy the TSV tables directly into Aspen Plus property
   input sheets, or download as CSV.

4. **Tab 6 — Validation:** Run Riazi Examples 4.13 and 4.14 in-app. The tool reproduces
   all published values within the stated tolerances without leaving the browser.

---

## Methodology

Full methodology with equations, references, and architectural decisions:
[`docs/methodology.md`](docs/methodology.md)

The methodology document is structured as a draft for Paper 1 Section 2.

---

## Validation

Test suite (372 tests, all passing):

```bash
python -m pytest
```

Numerical validation report comparing tool output against Riazi MNL50 textbook
examples and the synthetic VTB-like test feed: [`docs/validation_report.md`](docs/validation_report.md).

---

## Architectural commitments

These rules are enforced at every level of the codebase and are non-negotiable:

- **No solvent-deasphalting data required.** All inputs are routine refinery assay
  measurements (distillation curve, bulk SG, bulk MW, SARA wt%).
- **No tuning against process observables.** Pseudo-component properties (Tb, SG,
  PC-SAFT parameters) are determined entirely by published correlations. There are
  no adjustable parameters fitted to asphaltene onset pressure data, bubble-point
  pressure data, or density measurements.
- **Asphaltenes always discrete.** The asphaltene component is a single discrete
  pseudo-component using Gonzalez (2007) nanoaggregate defaults (m = 33, σ = 4.3 Å,
  ε/k = 400 K). It does not enter the distillation-curve fitting domain.
- **Bulk properties as inputs, not targets.** Bulk MW and bulk SG are feed
  characterization inputs. The pseudo-component distribution is fitted to the
  distillation curve directly; the bulk properties serve as closure diagnostics,
  not as tuning targets.

---

## Citation

If you use petrochar in your research, please cite the companion paper once published:

```
Paradkar, P. (2026). petrochar: Continuous-distribution characterization of heavy
petroleum fractions for PC-SAFT phase-equilibrium modeling from routine refinery data.
Energy & Fuels (in preparation).
```

*This citation will be updated upon journal publication.*

---

## License

MIT License. See [LICENSE](LICENSE).

---

## Acknowledgments

petrochar implements methodology from the following primary sources:

- Riazi, M. R. (2005). *Characterization and Properties of Petroleum Fractions*.
  ASTM Manual MNL50. ISBN 0-8031-3361-8.
  (Distribution fitting §4.5.4.1, quadrature §4.6.1.1, Watson K §2.4.1.)
- Panuganti, S. R., Tavakkoli, M., Vargas, F. M., Gonzalez, D. L., Chapman, W. G.
  (2012). PC-SAFT characterization of crude oils and modeling of asphaltene phase
  behavior. *Fuel*, 93, 658–669.
  (PC-SAFT parameter correlations, Table 6.)
- Gonzalez, D. L., Hirasaki, G. J., Creek, J., Chapman, W. G. (2007). Modeling study
  of CO₂-induced asphaltene precipitation. *Energy & Fuels*, 21, 1230–1234.
  (Asphaltene nanoaggregate PC-SAFT defaults: m = 33, σ = 4.3 Å, ε/k = 400 K.)
- Watson, K. M. and Nelson, E. F. (1933). Improved methods for approximating critical
  and thermal properties of petroleum fractions. *Industrial & Engineering Chemistry*,
  25, 880.
  (Watson K factor, Eq. 2.13.)
