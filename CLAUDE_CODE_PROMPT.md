# CLAUDE CODE PROMPT — petrochar

## Read this entire file before writing any code.

You are starting a new project, **petrochar**: a standalone Python tool for characterizing heavy petroleum fractions into discrete pseudo-components with PC-SAFT parameters, suitable as input to Aspen Plus or any PC-SAFT-capable simulator.

This tool is the methodology contribution for an upcoming research paper:
**"Continuous-distribution characterization of heavy petroleum fractions for PC-SAFT phase-equilibrium modeling from routine refinery data"**, target journals: *Energy & Fuels* or *Fuel*.

This is **Paper 1** in a two-paper sequence. Paper 2 (a propane-deasphalting application) is in a separate, unrelated repository and you will not interact with it during this work.

---

## STRICT BOUNDARIES — READ FIRST

1. **Do NOT read or reference any file outside this repository.** Specifically, there is an older project at `~/projects/pda_pcsaft_tool/` (or similar location). You are forbidden from reading it. It contains a different architecture (a 7-component "L/H" framework) that this project deliberately moves away from. Reading it would corrupt your design decisions through pattern-matching on the wrong abstractions.

2. **Build everything from primary sources only.** The acceptable references are:
   - Riazi, M. R. (2005). *Characterization and Properties of Petroleum Fractions*. ASTM Manual MNL50. Especially Chapter 4 (Sections 4.5, 4.6).
   - Panuganti, S. R. et al. (2012). PC-SAFT characterization of crude oils and modeling of asphaltene phase behavior. *Fuel*, 93, 658-669.
   - Gonzalez, D. L. et al. (2007). Modeling study of CO₂-induced asphaltene precipitation. *Energy & Fuels*, 21, 1230-1234. (Asphaltene PC-SAFT defaults: m=33, σ=4.3 Å, ε/k=400 K.)
   - Watson, K. M. and Nelson, E. F. (1933). *Industrial & Engineering Chemistry*, 25, 880. (Watson K factor.)
   - Whitson, C. H. (1983). Characterizing hydrocarbon plus fractions. *SPE Journal*, 23(4), 683-694.
   - García Cárdenas, J. and Ancheyta, J. (2022). Modeling of the deasphalting process. *Industrial & Engineering Chemistry Research*, 61, 3383-3394.

3. **No DAO data. No process-side data. Ever.** This tool's contract is feed-side characterization only. The codebase must contain no DAO inputs, no DAO outputs, no DAO references. Validation is against textbook examples and bulk-property closures only.

4. **No phase-equilibrium computation.** No flashes. No LLE solvers. No iterative chemical-potential equality. This tool produces parameters; downstream simulators (Aspen Plus, etc.) compute phase splits.

5. **No tuning against process observables.** All pseudo-component properties (MW, Tb, SG, distribution shape) come from feed measurements and published correlations. Nothing is fitted against DAO yield or composition.

---

## OBJECTIVE AND CONTRACT

**Inputs (all routinely measured at any operating refinery):**
- Distillation curve (any ASTM method: D86, D1160 AET, D7169, TBP). Cumulative fraction (weight, volume, or mole) vs Tb in K or °C.
- Bulk specific gravity at 15°C.
- Bulk molecular weight (lab VPO or correlation estimate).
- SARA wt% (saturates, aromatics, resins, asphaltenes), summing to 100%.

**Outputs:**
- Discrete pseudo-component table: for each component, (z_i, M_i, Tb_K, SG, K_W, γ, m, σ, ε/k).
- Aspen-Plus-paste-ready TSV/CSV exports.
- Distribution diagnostics (PDFs, CDFs, fit-quality metrics, closure errors).
- Validation report comparing tool output against Riazi MNL50 textbook examples.

**Architecture commitment:**
- Discretization via Gaussian quadrature (Riazi §4.6.1.1), default 5 points + 1 discrete asphaltene component + propane = 7 components total. User-selectable 3-point variant.
- Distribution model: Riazi generalized model (Eq. 4.56) for Tb and SG. Gamma model (Eq. 4.31) available as a legacy option but documented as inappropriate for VTB-scale fractions per Riazi p. 178-179.
- Asphaltenes are always discrete with literature defaults (Gonzalez 2007). They never enter the distillation-curve fitting domain.
- PC-SAFT parameters via Panuganti 2012 Table 6 correlations for distillable pseudo-components; Gonzalez 2007 for asphaltenes; pure-component values for propane.
- Aromaticity factor γ derived per pseudo-component from its individual Watson K factor (computed from its own Tb_i and SG_i), via the linear clamp `γ = clamp((13.0 − K_W) / (13.0 − 9.5), 0, 1)`.

---

## DEVELOPMENT PROTOCOL — STRICT

1. **No shortcuts. No stubs.** Every function is fully implemented. No `pass`, no `NotImplementedError`, no commented-out fallbacks.
2. **Validation-first.** Each numerical module has a pass-gate test that reproduces a Riazi textbook example to within stated tolerance. **You may not proceed to the next phase until the current phase's pass-gate is green.** This is not negotiable.
3. **Citation density.** Every correlation, equation, and constant has a docstring citing the source: page and equation number for Riazi MNL50; section and equation number for Panuganti 2012 and Gonzalez 2007.
4. **Units explicit, always.** Every function signature documents input and output units. Tb in K internally; convert at I/O boundaries. SARA in wt% (not fractions). Distillation curves carry an explicit basis flag (cumulative weight, volume, or mole fraction).
5. **No proprietary names anywhere.** No client refinery, no plant numbers, no internal codes. UI says "VTB sample" or "heavy petroleum fraction".
6. **Session protocol.** At the start of every session, read `CURRENT_STATUS.md`. At the end of every session, update it: mark phase complete, log decisions, append session entry, commit.

---

## REPOSITORY LAYOUT

Create exactly this structure in Phase 0:

```
petrochar/
├── README.md
├── LICENSE                              # MIT
├── requirements.txt
├── runtime.txt                          # python-3.11
├── pyproject.toml                       # for later pip-installable packaging
├── CLAUDE.md                            # project rules, mirror of section "STRICT BOUNDARIES" + "DEVELOPMENT PROTOCOL"
├── CLAUDE_CODE_PROMPT.md                # this file (committed for traceability)
├── CURRENT_STATUS.md                    # session log
├── .gitignore
├── app.py                               # Streamlit entry; built last, in Phase 9
├── core/
│   ├── __init__.py
│   ├── correlations.py                  # Phase 1
│   ├── distillation.py                  # Phase 2
│   ├── distribution.py                  # Phase 3
│   ├── sg_distribution.py               # Phase 4
│   ├── mw_distribution.py               # Phase 4
│   ├── quadrature.py                    # Phase 5
│   ├── sara.py                          # Phase 6
│   ├── watson_k.py                      # Phase 7
│   └── pcsaft_params.py                 # Phase 7
├── tests/
│   ├── __init__.py
│   ├── test_phase1_correlations.py
│   ├── test_phase2_distillation.py
│   ├── test_phase3_riazi_example_4_13.py
│   ├── test_phase4_sg_mw.py
│   ├── test_phase5_riazi_example_4_14.py
│   ├── test_phase6_sara_closure.py
│   ├── test_phase7_pcsaft.py
│   └── test_phase8_pipeline.py
├── data/
│   └── riazi_reference/
│       ├── README.md                    # provenance and citation
│       ├── table_4_6_scn_groups.csv
│       ├── table_4_11_example_4_7.csv
│       ├── table_4_13_distribution_coeffs.csv
│       ├── table_4_22_quadrature_3pt.csv
│       └── table_4_23_lumping_methods.csv
├── tabs/                                # Phase 9, Streamlit UI modules
│   ├── __init__.py
│   ├── input_data.py
│   ├── distillation_fit.py
│   ├── distributions.py
│   ├── pseudocomponents.py
│   ├── pcsaft_export.py
│   └── validation.py
└── docs/
    ├── methodology.md                   # draft Paper 1 Section 2
    └── validation_report.md             # Phase 8 deliverable; draft Paper 1 Section 3
```

Initialize git in Phase 0. Confirm `pytest` runs (with no tests yet, exit code zero).

---

## PHASES

Each phase is a single Claude Code session (or two if needed). Phase pass-gates are non-negotiable.

### Phase 0 — Repository scaffold
Create directory structure above. Write minimal `README.md`, `LICENSE` (MIT), `requirements.txt` (numpy, scipy, pandas, matplotlib, pytest, streamlit, openpyxl), `runtime.txt` (`python-3.11`), `pyproject.toml` (project name `petrochar`, version `0.1.0`), `CLAUDE.md`, `CURRENT_STATUS.md`, `.gitignore`. Commit. Vendor Riazi reference CSVs into `data/riazi_reference/` per `data/riazi_reference/README.md` (see "Reference data" section below).

### Phase 1 — Shared correlations
`core/correlations.py`. Functions: `riazi_daubert_Tb(M, SG)`, `riazi_daubert_M(Tb_K, SG)`, `watson_K(Tb_K, SG)`, `gamma_function(x)`, `incomplete_gamma_upper(a, q)`, `aromaticity_to_gamma(K_W)`. Each function fully implemented with docstring citing Riazi MNL50 page/equation, Panuganti 2012, Watson 1933, etc.

**Pass gate (`tests/test_phase1_correlations.py`):**
- For each row of Riazi Table 4.11 (`data/riazi_reference/table_4_11_example_4_7.csv`), `riazi_daubert_Tb(M, SG)` reproduces tabulated Tb to ±2 K.
- `riazi_daubert_M(Tb_K, SG)` (numerical inverse via `scipy.optimize.brentq`) reproduces M to ±1.5%.
- `watson_K(Tb_K=350, SG=0.7597)` returns ~12.0 (Riazi Example 4.7 bulk values).
- `gamma_function(2.0) == 1.0` exactly (recurrence base case).
- `gamma_function(0.5)` reproduces √π to ±0.1%.

### Phase 2 — Distillation curve handling
`core/distillation.py`. Class `DistillationCurve` with method tags {D86, D1160_AET, D7169, TBP} and basis tags {weight, volume, mole}. Internal canonical form: cumulative weight fraction `x_cw` array vs `Tb_K` array. Implements:
- D86 → TBP conversion via Riazi §3.1.3 Eqs. 3.18-3.20 (Daubert procedure).
- D1160 AET: no temperature conversion (data is already AET); only unit handling.
- D7169, TBP: take as-is in K.
- Volume → weight conversion using bulk SG distribution (if available) or method-flag refusal.
- Validity flags: D86 unreliable above 370°C, D1160 capped at user-stated cracking limit, D7169 capped at ~720°C.
- `.measured_max_T_K`, `.cutoff_warnings`, `.to_TBP()` exposed.

**Pass gate:** unit test using a fabricated D86 curve with hand-computed expected D86→TBP conversion (from Riazi Eqs. 3.18-3.20). Tolerance ±5 K at any cumulative point.

### Phase 3 — Generalized distribution fit (Riazi Eq. 4.56)
`core/distribution.py`. Class `GeneralizedDistribution` implementing Riazi §4.5.4.1.

Mathematics: P* = [(A/B) ln(1/(1−x_c))]^(1/B), where P* = (P − P_o) / P_o.

Fitting via `scipy.optimize.least_squares` on linearized form (Eq. 4.57). Defaults per Riazi p. 174: `B_T ≈ 1.5` for Tb. Two-parameter fit (B fixed) and three-parameter fit (all of P_o, A, B free) both supported, user-selectable. PDF (Eq. 4.70) and CDF (Eq. 4.56) exposed. Average property via Eq. 4.74-4.75.

**Pass gate (`tests/test_phase3_riazi_example_4_13.py`):**
Reproduce Riazi Example 4.13 (Method B with TBP curve). Input data: Table 4.11 columns 3 (Tb,K) and 5 (x_cw). Reference output: Table 4.13 row "T_b" — T_o = 350 K, A_T = 0.1679, B_T = 1.2586 (three-parameter fit).
Tolerance: T_o ±5 K, A_T ±5%, B_T ±5%.

### Phase 4 — SG and MW distributions
`core/sg_distribution.py` and `core/mw_distribution.py`.

SG distribution: two methods, both implemented:
1. Constant Watson K (Whitson; Riazi p. 178). Compute K_W from bulk M and bulk SG; assume same K_W for all subfractions; SG_i derived from Tb_i and K_W.
2. Generalized fit via Eq. 4.56 with B_SG ≈ 3 default (Riazi p. 174, Table 4.13).

Default to constant Watson K for VTB-scale fractions. User-toggleable.

MW distribution: for each cut along x_cw, given Tb_i and SG_i, compute M_i via `riazi_daubert_M`. The bulk M closure (harmonic-mean integral, Riazi Eq. 4.32) is a verification, not a constraint.

**Pass gate (`tests/test_phase4_sg_mw.py`):**
- Reproduce Riazi Table 4.13 SG row: SG_o = 0.705, A_SG = 0.0232, B_SG = 1.811, tolerance ±5%.
- Bulk M_av and SG_av from Eqs. 4.74-4.76 within 1% of Riazi Example 4.7 measured values (M_7+ = 118.9, SG_7+ = 0.7597).

### Phase 5 — Gaussian quadrature discretization
`core/quadrature.py`. Reference: Riazi MNL50 §4.6.1.1, Table 4.21.

Hard-code quadrature points/weights for N ∈ {3, 5}. Implement `discretize_generalized(N, distribution)` using Eqs. 4.83-4.91. (Eqs. 4.91-4.92 with B = 1 shortcut as fallback for light-oil / gas-condensate cases; not the default.)

Output: list of `Pseudocomponent` records with fields `(z_i, M_i, Tb_K, SG, x_cw_lower, x_cw_upper)`.

**Pass gate (`tests/test_phase5_riazi_example_4_14.py`):**
Reproduce Riazi Example 4.14 (3-point Gaussian quadrature on Table 4.13 distribution). Reference output is Table 4.22:

| i | y_i  | w_i = z_i | M_i   |
|---|------|-----------|-------|
| 1 | 0.416 | 0.711    | 103.6 |
| 2 | 2.294 | 0.279    | 154.6 |
| 3 | 6.290 | 0.010    | 252.2 |
| Mixture | — | 1.000 | 119.4 |

Tolerance: ±1% on M_i, ±0.005 on z_i. Mixture M_av must match Riazi's stated 118.9 ±0.5%.

### Phase 6 — SARA closure and asphaltene handling
`core/sara.py`. 

Logic:
1. User supplies SARA wt% summing to 100.
2. ASP wt% is removed from the distillation-curve domain. Distillation curve is verified to represent only (1 − ASP_wt%/100) of total mass; if user-supplied curve includes ASP mass, renormalize and warn.
3. After Phase 5 produces N distillable pseudo-components, append a single discrete ASP component:
   - Mass fraction: ASP_wt%/100
   - MW: 1700 g/mol (Gonzalez 2007 nanoaggregate)
   - Tb: 800 °C (numerical convention, document explicitly)
   - SG: 1.15 (convention)
   - PC-SAFT params: Gonzalez 2007 defaults (handled in Phase 7)
4. K_W-bin closure check: aggregate distillable pseudo-components into nominal SAT/ARO/RES classes by K_W:
   - K_W ≥ 12.0 → SAT
   - 11.0 ≤ K_W < 12.0 → ARO
   - K_W < 11.0 → RES
   Compare to user-supplied SARA. Flag deviation > ±3 wt% per class as data-consistency warning. **This is a closure check, not a constraint — pseudo-components are not retuned to match SARA.**

K_W bin thresholds documented as conventions; user-overridable. Cite Riazi p. 75.

**Pass gate (`tests/test_phase6_sara_closure.py`):**
Synthetic test with a constructed three-component model fluid (paraffin, alkylbenzene, naphthenic) with known mass fractions. Run distillation curve through full pipeline; verify K_W binning recovers input class fractions to ±3 wt%.

### Phase 7 — PC-SAFT parameter generation
`core/watson_k.py` and `core/pcsaft_params.py`.

`watson_k.py`: `compute_K_W_per_pseudocomponent(pseudocomp_list)` — adds K_W and γ to each component using Phase 1's `aromaticity_to_gamma`.

`pcsaft_params.py`:
- `panuganti_alkyl_aromatic_resin_params(M, gamma)` — Panuganti 2012 *Fuel* 93, Table 6, alkyl-aromatic / resin family correlations.
- `panuganti_saturate_params(M)` — saturate end-member correlations from same source.
- `gonzalez_asphaltene_params()` — fixed: m = 33, σ = 4.3 Å, ε/k = 400 K. Cite Gonzalez 2007 *Energy & Fuels* 21, 1230.
- `propane_params()` — fixed pure-component: m = 2.002, σ = 3.6180 Å, ε/k = 208.11 K. **σ = 3.6180 Å, not 3.168 Å** — this is a known typo in some Aspen built-in databases.
- `generate_pcsaft_table(pseudocomp_list)` — top-level: returns pandas DataFrame, all components, all parameters, ready for Phase 9 export.

**Pass gate (`tests/test_phase7_pcsaft.py`):**
For three reference (M, γ) pairs from Panuganti 2012 Table 6 worked examples, verify (m, σ, ε/k) reproduce to ±2%. If exact reference values are unavailable in the published table, document the correlation formulas explicitly in code with equation numbers, and verify continuity properties: m monotonically increasing in M; σ in [3.5, 4.5] Å; ε/k in [200, 450] K for petroleum pseudocomponents.

### Phase 8 — End-to-end pipeline integration
`tests/test_phase8_pipeline.py` and `docs/validation_report.md`.

Run full pipeline (Phases 1-7) on a synthetic VTB-like input:
- Distillation curve: fabricated D1160 AET, IBP 280 °C, FBP 540 °C measured, 11 cumulative-weight points.
- Bulk SG = 1.020. Bulk MW = 700 g/mol. SARA = 12 / 38 / 38 / 12 wt%.

Verify:
- Distillation fit residual < 5 K RMS.
- Bulk MW closure < 2%.
- Bulk SG closure < 2%.
- SARA K_W-bin closure within ±3 wt% per class.
- Output: 5 quadrature points + 1 ASP = 6 pseudo-components (+ propane in Aspen export).
- All PC-SAFT parameters within physical bounds.

Generate `docs/validation_report.md` with all Phase 3, 4, 5, 6, 7, 8 pass-gate results in tables. **This document is the draft for Paper 1 Section 3 (Method validation).**

This is the moment the engine is declared correct. **Do not proceed to Phase 9 until this passes.**

### Phase 9 — Streamlit UI
`app.py` plus `tabs/*.py`. Six tabs:

1. **Input Data.** File upload (Excel/CSV) or manual entry. Distillation method/basis selectors. Bulk SG, bulk MW, SARA inputs. Live validity warnings.
2. **Distillation Fit.** Visual overlay: measured points, fitted Eq. 4.56 curve, extrapolated tail. Fit-quality metrics (RMS, AAD, R²). SG-distribution method selector.
3. **Distributions.** PDFs and CDFs for Tb, M, SG. Continuous distributions overlaid with quadrature point markers.
4. **Pseudocomponents.** Discrete table: z_i, M_i, Tb_i, SG_i, K_W, γ, m, σ, ε/k. ASP component highlighted. Quadrature-point selector (3 / 5).
5. **PC-SAFT Export.** Aspen-paste-ready tables with TSV/CSV downloads. Pure-component params, hypothetical properties, propane parameters, BIP matrix (default 0.010 for all C3-pseudocomponent pairs, 0 elsewhere).
6. **Validation.** Built-in test mode: run Riazi Examples 4.13, 4.14, 4.15 in-app and display computed-vs-published values side-by-side. This is a transparency feature for paper reviewers and end users.

UI must contain no proprietary names. Generic only.

### Phase 10 — Documentation and packaging
- `README.md` — full user-facing documentation: install, run, inputs, outputs, citation.
- `docs/methodology.md` — structured as draft Paper 1 Section 2 (Method) and Section 3 (Validation). Each subsection cites Riazi page/equation, Panuganti, Gonzalez. This document matures into the manuscript.
- `pyproject.toml` — finalize for `pip install petrochar`. So Paper 2's PDA tool can later import this as a versioned dependency.

---

## REFERENCE DATA — `data/riazi_reference/*.csv`

Vendor the following reference data points from Riazi MNL50 Chapter 4. These are factual data points (boiling points, densities, fitted distribution coefficients) and reproducing them as test fixtures is fair use. Each CSV must have a header row, an `__about__` row at the top documenting source citation, and units explicit in column names.

**`table_4_11_example_4_7.csv`** — North Sea gas condensate C7+ fraction, Riazi Example 4.7. Columns: `Fraction_No`, `Carbon_No`, `x_w`, `M_g_per_mol`, `Tb_K`, `SG`, `x_m`, `x_v`, `x_cm`, `x_cw`, `x_cv`. 12 rows + header.

**`table_4_13_distribution_coeffs.csv`** — Distribution coefficients for Example 4.7 system. Columns: `Property`, `Type_of_x_c`, `P_o`, `A`, `B`, `RMS`, `pct_AAD`, `R_squared`. Rows: M (mole basis), Tb (weight basis), SG (volume basis), SG (weight basis), and the fixed-B variants.

**`table_4_22_quadrature_3pt.csv`** — Three-point Gaussian quadrature output for Example 4.14. Columns: `i`, `y_i`, `w_i_eq_z_i`, `M_i`, `z_i_M_i`. 3 rows + mixture row + header.

**`table_4_23_lumping_methods.csv`** — Five-component lumping for Example 4.15 (North Sea oil). Columns: `Component_i`, `Method`, `mole_fraction`, `weight_fraction`, `M_i`, `SG_i`. Two methods × 5 components + 2 mixture rows.

**`table_4_6_scn_groups.csv`** — SCN group properties C5-C45 from Riazi Table 4.6. Columns: `SCN`, `Carbon_No`, `M_g_per_mol`, `Tb_K`, `SG`, `n_20`. 41 rows + header. (Used as fallback in M and SG correlations when distillation data absent.)

`data/riazi_reference/README.md` — citation, ASTM Manual MNL50 (2005) by M. R. Riazi, ISBN 0-8031-3361-8, ASTM International, West Conshohocken PA. Note that these tables are reproduced for software validation purposes and acknowledge the source on every CSV.

---

## CURRENT_STATUS.md TEMPLATE — initialize in Phase 0

```markdown
# petrochar — Current Status

## Project
Standalone Python tool for heavy petroleum fraction characterization (Paper 1).

## Current Phase
Phase 0 — Repository scaffold

## Phases Completed
(none yet)

## Decisions Made
1. Default distribution model: Riazi generalized (Eq. 4.56), not gamma. Gamma fails on heavy oils per Riazi p. 178-179.
2. Default quadrature: 5-point + 1 ASP discrete + propane = 7 components.
3. Asphaltenes always discrete with Gonzalez 2007 defaults (m=33, σ=4.3, ε/k=400). Tb=800°C and SG=1.15 are numerical conventions.
4. Aromaticity γ via per-component Watson K factor (linear clamp 9.5-13.0).
5. SARA: ASP wt% used as hard tail constraint; SAT/ARO/RES used as K_W-bin closure check (NOT a constraint, never tuned to).
6. No DAO data in this repository, ever.
7. No reading of `~/projects/pda_pcsaft_tool/` or any external project, ever.
8. Validation-first: each phase requires reproducing a Riazi textbook example before proceeding.

## Known Issues / Blockers
(none yet)

## Session Log
- YYYY-MM-DD | Phase 0 — repo scaffold + reference CSVs vendored | next: Phase 1 correlations
```

---

## SESSION 0 — WHAT TO DO FIRST

1. Confirm the working directory is `~/projects/petrochar/` (or wherever the user instructed) and is empty.
2. Create the directory tree above.
3. Write `LICENSE` (MIT, copyright "Piyush Paradkar 2026").
4. Write `requirements.txt`:
   ```
   numpy>=1.24
   scipy>=1.10
   pandas>=2.0
   matplotlib>=3.7
   pytest>=7.4
   streamlit>=1.30
   openpyxl>=3.1
   ```
5. Write `runtime.txt`: `python-3.11`
6. Write `pyproject.toml` with project name `petrochar`, version `0.1.0`, description `"Continuous-distribution characterization of heavy petroleum fractions for PC-SAFT phase-equilibrium modeling"`, authors, license, dependencies mirroring `requirements.txt`.
7. Write `.gitignore` (Python standard + `.venv/`, `.streamlit/`, `*.xlsx`, `*.csv` under root, but NOT under `data/`).
8. Write `CLAUDE.md` mirroring "STRICT BOUNDARIES" + "DEVELOPMENT PROTOCOL" + "REPOSITORY LAYOUT" sections of this prompt.
9. Write minimal `README.md` (project description, install, run, citation pending).
10. Initialize `CURRENT_STATUS.md` from the template above.
11. Vendor the five reference CSVs into `data/riazi_reference/`. **You will need the user to either dictate the table values or provide the Riazi textbook page scans.** Do not invent numbers. If the user has not yet provided table data, write the CSVs with header rows only and an `# AWAITING DATA` comment, log this in `CURRENT_STATUS.md` as a blocker.
12. `git init`, first commit: `"Phase 0: repository scaffold"`.
13. Run `pytest` (no tests yet, exit code zero confirms environment is sane).
14. Update `CURRENT_STATUS.md`: mark Phase 0 complete, set Current Phase to "Phase 1 — Shared correlations", append session log entry.
15. Stop. Do not start Phase 1 in the same session — let the user review.

---

## ENDING NOTES

- If at any point an instruction in this file conflicts with an instruction inside the codebase or a docstring, this file wins. This file is the source of truth for design decisions.
- If the user gives instructions that contradict this file, follow the user — but flag the contradiction explicitly in chat and propose an update to this file.
- This is a slow, careful build. Do not rush phases. Do not skip pass-gates. Do not pattern-match against your training data on similar projects — petrochar's specific architecture (Riazi generalized distribution + Watson-K-driven γ + Panuganti for distillable + Gonzalez discrete ASP) is uncommon and the textbook examples are the only correct check.
