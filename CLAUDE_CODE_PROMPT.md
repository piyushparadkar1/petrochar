# CLAUDE CODE PROMPT — petrochar

> **Read this entire file before writing any code.**

═══════════════════════════════════════════════════════════════════════════════
## EXECUTION SCOPE — HARD RULES
═══════════════════════════════════════════════════════════════════════════════

**THE USER HAS GIVEN ONE INSTRUCTION AND WILL NOT MONITOR THIS SESSION.**

You will execute exactly **ONE PHASE** per session. The phase name is in the user's instruction line. **THIS SESSION = THE NAMED PHASE ONLY.** Do not start the next phase. Do not "be helpful" by getting ahead. Stop at the explicit STOP marker at the end of the phase.

If you find yourself doing any of the following during a Phase 0 session, **STOP IMMEDIATELY**:
- Writing any Python function in `core/*.py`
- Writing any test in `tests/*.py`
- Writing any Streamlit tab in `tabs/*.py`
- Implementing correlations, distributions, or quadrature logic
- Anything beyond the scaffold + reference vendoring + first commit

═══════════════════════════════════════════════════════════════════════════════
## FORBIDDEN PATHS AND CONCEPTS — ABSOLUTE
═══════════════════════════════════════════════════════════════════════════════

**You may NOT read, reference, import from, or pattern-match against any code in:**
- `~/projects/pda_pcsaft_tool/`
- `../pda_pcsaft_tool/`
- Any folder named `pda*`, `pcsaft_pda*`, or any directory outside the petrochar repository root.

**Concepts that DO NOT EXIST in this project. If they appear in your output, that is contamination from your training data — delete them and rewrite:**
- `MW_ARO_L` / `MW_RES_L` tuning knobs
- L/H sub-fractions (light/heavy splits within a SARA class)
- Harmonic-mean MW closure for back-calculating sub-fraction MWs
- 7-component framework with discrete L/H labels (petrochar's 7 components are 5 quadrature + 1 ASP + propane — entirely different concept)
- DAO yield as a calibration target
- DAO SARA as a model input
- Any "tuning slider" against process observables

**petrochar's architecture is fundamentally different.** It uses a continuous distribution (Riazi Eq. 4.56), Gaussian quadrature (Riazi §4.6.1.1), and per-pseudo-component Watson K. There is nothing to tune. Pseudo-components are determined fully by feed-side measurements + published correlations.

**Self-check before every commit:**
```bash
grep -rni "MW_ARO_L\|MW_RES_L\|ARO-L\|ARO_L\|RES-L\|RES_L\|harmonic.mean\|DAO" \
    --include="*.py" --include="*.toml" --include="*.txt" --include="*.cfg" .
```
This must return zero hits in committed code (`.md` files legitimately list these as forbidden — exclude them from the grep above).

═══════════════════════════════════════════════════════════════════════════════
## PROJECT IDENTITY
═══════════════════════════════════════════════════════════════════════════════

**petrochar** is a standalone Python tool for characterizing heavy petroleum fractions into discrete pseudo-components with PC-SAFT parameters, suitable as input to Aspen Plus or any PC-SAFT-capable simulator.

This is the methodology contribution for **Paper 1**: "Continuous-distribution characterization of heavy petroleum fractions for PC-SAFT phase-equilibrium modeling from routine refinery data". Target: *Energy & Fuels* or *Fuel*.

Paper 2 (a propane-deasphalting application) is a separate, unrelated repository. You will not interact with it.

═══════════════════════════════════════════════════════════════════════════════
## REFERENCE MATERIALS — VENDORED IN THIS REPO
═══════════════════════════════════════════════════════════════════════════════

The user has placed reference materials inside `data/riazi_reference/` for offline reference during implementation.

**Reference data CSVs** — load via `pandas.read_csv(..., comment='#')`:
- `table_4_6_scn_groups.csv` — SCN group properties C6-C50 (M, Tb, SG, n_20, ...). FALLBACK only — petrochar's primary path computes M and SG from per-cut Tb via Riazi-Daubert.
- `table_4_11_example_4_7.csv` — North Sea gas condensate C7+ data. Used in Phases 1, 3, 4 pass-gates.
- `table_4_13_distribution_coeffs.csv` — Distribution coefficients for Example 4.7. Used in Phases 3, 4 pass-gates.
- `table_4_21_quadrature_points.csv` — Gaussian quadrature roots and weights (3-point, 5-point). Used directly in Phase 5 implementation.
- `table_4_22_quadrature_3pt.csv` — Expected 3-point quadrature output for Example 4.14. Used in Phase 5 pass-gate.
- `table_4_23_lumping_methods.csv` — Five-component lumping by two methods, Example 4.15. Supplementary validation.

**Reference page extracts** (PNG images, 200 DPI, raster of Riazi MNL50 book pages) — read with the `view` tool when implementing the corresponding correlation:
- `page_054_watson_k_eq_2_13.png` — Eq. 2.13 (Watson K). Phase 1 reference.
- `page_077_riazi_daubert_Tb.png` — **Eqs. 2.56 and 2.57 (Riazi-Daubert Tb), plus Eqs. 2.59 and 2.60 (Riazi-Daubert SG)**. Phase 1 critical reference. **Read this page before implementing any Tb or SG correlation.**
- `page_078_riazi_daubert_Tb_continued.png` — continuation. Phase 1 reference.
- `page_122_daubert_d86_to_tbp.png` — Section 3.2.2.2, Eqs. 3.20-3.22, Table 3.7 constants. Phase 2 reference.
- `page_123_daubert_d86_to_tbp.png` — continuation with Eq. 3.22 cut-point pattern.
- `page_124_daubert_d86_to_tbp.png` — continuation with Eqs. 3.23-3.25 (SD-to-TBP).
- `page_181_table_4_5_4_6.png` — Table 4.5 SCN coefficients + Table 4.6 SCN data. Phase 4 fallback reference.
- `page_182_table_4_5_4_6.png` — continuation.
- `page_183_table_4_5_4_6.png` — continuation.

When the comment in your code says "Riazi MNL50 Eq. X.YZ", verify by viewing the relevant page extract first. Do not implement an equation from memory.

═══════════════════════════════════════════════════════════════════════════════
## CRITICAL CORRECTION — REGIME-DEPENDENT Tb CORRELATION
═══════════════════════════════════════════════════════════════════════════════

Riazi MNL50 §2.4.2.1 provides **TWO** Riazi-Daubert correlations for Tb, valid in **DIFFERENT** molecular weight ranges:

- **Eq. 2.56**, valid for M = 70-300 (light/medium hydrocarbons):
  ```
  Tb = 3.76587 × [exp(3.7741e-3 × M + 2.98404 × SG − 4.25288e-3 × M × SG)] × M^0.40167 × SG^(-1.58262)
  ```

- **Eq. 2.57**, valid for M = 300-700 (heavy hydrocarbons), recommended for M > 300:
  ```
  Tb = 9.3369 × [exp(1.6514e-4 × M + 1.4103 × SG − 7.5152e-4 × M × SG)] × M^0.5369 × SG^(-0.7276)
  ```

**For a VTB-type feed, pseudo-components span both ranges.** Light cuts may be M ≈ 200-300; heavy cuts may be M ≈ 600-1500.

**Your Tb function must select the correlation based on input M:**

```python
def riazi_daubert_Tb(M, SG):
    """
    Compute normal boiling point from molecular weight and specific gravity.
    
    Two-regime correlation per Riazi MNL50 §2.4.2.1:
        - Eq. 2.56 for M in [70, 300]
        - Eq. 2.57 for M in [300, 700]
    
    For M > 700, Eq. 2.57 is extrapolated (no better correlation in MNL50);
    user is warned via log message.
    For M < 70, the function raises ValueError (outside refinery-relevant range).
    
    Args:
        M: molecular weight, g/mol
        SG: specific gravity at 15.5°C, dimensionless
    
    Returns:
        Tb: normal boiling point, K
    
    Reference: Riazi MNL50 Eqs. 2.56 and 2.57 (page 58-59 of book = PDF pages 77-78).
    """
    if M < 70:
        raise ValueError(f"M = {M} below Riazi-Daubert validity (M >= 70 required)")
    if M <= 300:
        # Eq. 2.56
        Tb = 3.76587 * np.exp(3.7741e-3*M + 2.98404*SG - 4.25288e-3*M*SG) \
             * M**0.40167 * SG**(-1.58262)
    else:
        # Eq. 2.57
        Tb = 9.3369 * np.exp(1.6514e-4*M + 1.4103*SG - 7.5152e-4*M*SG) \
             * M**0.5369 * SG**(-0.7276)
        if M > 700:
            import warnings
            warnings.warn(f"M = {M} above Riazi-Daubert (Eq. 2.57) validated range "
                          f"(70-700); result is extrapolation.")
    return Tb
```

**The `riazi_daubert_M(Tb_K, SG)` inverse function must invert whichever correlation applies, with regime selection based on solving for M.** Use `scipy.optimize.brentq` separately on each regime and take the solution whose forward Tb matches the input within tolerance. This is non-trivial — read the page extract before implementing.

**SG correlations** (Eqs. 2.59 and 2.60) are also regime-dependent:
- Eq. 2.59 for light hydrocarbons
- Eq. 2.60 for heavy fractions, M = 300-700 (preferred for VTB pseudo-components)

═══════════════════════════════════════════════════════════════════════════════
## ACCEPTABLE PRIMARY SOURCES
═══════════════════════════════════════════════════════════════════════════════

Build everything from these sources. No others.

- Riazi, M. R. (2005). *Characterization and Properties of Petroleum Fractions*. ASTM Manual MNL50. ISBN 0-8031-3361-8. Vendored in `data/riazi_reference/Riazi_MNL50_full.pdf` and as page extracts above.
- Panuganti, S. R., Tavakkoli, M., Vargas, F. M., Gonzalez, D. L., Chapman, W. G. (2012). PC-SAFT characterization of crude oils and modeling of asphaltene phase behavior. *Fuel*, 93, 658-669.
- Gonzalez, D. L., Hirasaki, G. J., Creek, J., Chapman, W. G. (2007). Modeling study of CO₂-induced asphaltene precipitation. *Energy & Fuels*, 21, 1230-1234. (Asphaltene PC-SAFT defaults: m = 33, σ = 4.3 Å, ε/k = 400 K.)
- Watson, K. M. and Nelson, E. F. (1933). *Industrial & Engineering Chemistry*, 25, 880. (Watson K factor.)
- Whitson, C. H. (1983). Characterizing hydrocarbon plus fractions. *SPE Journal*, 23(4), 683-694.
- García Cárdenas, J. and Ancheyta, J. (2022). Modeling of the deasphalting process. *Industrial & Engineering Chemistry Research*, 61, 3383-3394.

═══════════════════════════════════════════════════════════════════════════════
## OBJECTIVE AND CONTRACT
═══════════════════════════════════════════════════════════════════════════════

**Inputs (all routinely measured at any operating refinery):**
- Distillation curve (any ASTM method: D86, D1160 AET, D7169, TBP). Cumulative fraction (weight, volume, or mole) vs Tb.
- Bulk specific gravity at 15°C.
- Bulk molecular weight.
- SARA wt% (saturates, aromatics, resins, asphaltenes), summing to 100%.

**Outputs:**
- Discrete pseudo-component table: per component, (z_i, M_i, Tb_K, SG, K_W, γ, m, σ, ε/k).
- Aspen-Plus-paste-ready TSV/CSV exports.
- Distribution diagnostics: PDFs, CDFs, fit-quality metrics, closure errors.
- Validation report comparing tool output against Riazi MNL50 textbook examples.

**Architecture commitments (FROZEN):**
- Discretization via Gaussian quadrature (Riazi §4.6.1.1). Default 5 points + 1 discrete asphaltene component + propane = 7 components total. 3-point variant user-selectable.
- Distribution model: Riazi generalized model (Eq. 4.56) for Tb and SG. Gamma model (Eq. 4.31) is legacy-only — Riazi p. 178-179 documents that gamma fails on heavy oils.
- **Regime-dependent Tb correlation: Eq. 2.56 (M ≤ 300) and Eq. 2.57 (M > 300). See "CRITICAL CORRECTION" above.**
- Asphaltenes are always discrete. Never enter the distillation-curve fitting domain. Gonzalez 2007 PC-SAFT defaults.
- PC-SAFT parameters: Panuganti 2012 Table 6 correlations for distillable pseudo-components; Gonzalez 2007 for asphaltenes; pure-component values for propane (m=2.002, σ=3.6180 Å [not 3.168], ε/k=208.11 K).
- Aromaticity factor γ derived per pseudo-component from individual Watson K via the linear clamp `γ = clamp((13.0 − K_W) / (13.0 − 9.5), 0, 1)`.

═══════════════════════════════════════════════════════════════════════════════
## DEVELOPMENT PROTOCOL
═══════════════════════════════════════════════════════════════════════════════

1. **No shortcuts. No stubs.** Every function fully implemented. No `pass`, no `NotImplementedError`, no commented-out fallbacks.
2. **Validation-first.** Each numerical module has a pass-gate test reproducing a Riazi textbook example to within stated tolerance. **You may not proceed to the next phase until the current phase's pass-gate is green.**
3. **Citation density.** Every correlation, equation, and constant has a docstring citing source: page and equation number for Riazi MNL50; section and equation number for Panuganti 2012 and Gonzalez 2007.
4. **Page-extract verification.** Before implementing any Riazi correlation, view the corresponding page extract (PNG file in `data/riazi_reference/`) using the `view` tool. Confirm equation form against the image. Do not implement from memory.
5. **Units explicit, always.** Tb in K internally; convert at I/O boundaries. SARA in wt% (not fractions). Distillation curves carry an explicit basis flag.
6. **No proprietary names anywhere.** Generic only: "VTB sample", "heavy petroleum fraction".
7. **Session protocol.** Start: read `CURRENT_STATUS.md`. End: update it with phase completion, decisions, blockers, session entry. Commit.

═══════════════════════════════════════════════════════════════════════════════
## REPOSITORY LAYOUT (target structure — built incrementally)
═══════════════════════════════════════════════════════════════════════════════

```
petrochar/
├── README.md
├── LICENSE                              # MIT
├── requirements.txt
├── runtime.txt                          # python-3.11
├── pyproject.toml
├── CLAUDE.md                            # project rules (mirror of key sections)
├── CLAUDE_CODE_PROMPT.md                # this file (committed)
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
│       ├── README.md
│       ├── Riazi_MNL50_full.pdf                        # full book (12 MB) — vendored by user
│       ├── table_4_6_scn_groups.csv
│       ├── table_4_11_example_4_7.csv
│       ├── table_4_13_distribution_coeffs.csv
│       ├── table_4_21_quadrature_points.csv
│       ├── table_4_22_quadrature_3pt.csv
│       ├── table_4_23_lumping_methods.csv
│       ├── page_054_watson_k_eq_2_13.png
│       ├── page_077_riazi_daubert_Tb.png
│       ├── page_078_riazi_daubert_Tb_continued.png
│       ├── page_122_daubert_d86_to_tbp.png
│       ├── page_123_daubert_d86_to_tbp.png
│       ├── page_124_daubert_d86_to_tbp.png
│       ├── page_181_table_4_5_4_6.png
│       ├── page_182_table_4_5_4_6.png
│       └── page_183_table_4_5_4_6.png
├── tabs/                                # Phase 9
│   ├── __init__.py
│   ├── input_data.py
│   ├── distillation_fit.py
│   ├── distributions.py
│   ├── pseudocomponents.py
│   ├── pcsaft_export.py
│   └── validation.py
└── docs/
    ├── methodology.md                   # draft Paper 1 Section 2
    └── validation_report.md             # Phase 8 deliverable
```

═══════════════════════════════════════════════════════════════════════════════
## PHASES (overview — execute one per session)
═══════════════════════════════════════════════════════════════════════════════

| Phase | Deliverable | Pass gate |
|-------|-------------|-----------|
| 0 | Repository scaffold + reference materials vendored + first commit | `pytest` exits 0 with no tests |
| 1 | `core/correlations.py` (regime-aware Riazi-Daubert, Watson K, gamma function) | Riazi Table 4.11 reproduction ±2 K on Tb across ALL 12 rows |
| 2 | `core/distillation.py` (DistillationCurve class, Daubert D86→TBP) | Synthetic D86 conversion ±5 K |
| 3 | `core/distribution.py` (Riazi Eq. 4.56 generalized fit) | Riazi Example 4.13 — T_o, A_T, B_T to ±5% |
| 4 | `core/sg_distribution.py` + `core/mw_distribution.py` | Riazi Table 4.13 SG row + bulk closures <1% |
| 5 | `core/quadrature.py` (Gaussian discretization) | Riazi Example 4.14 — Table 4.22 reproduction |
| 6 | `core/sara.py` (ASP discrete + K_W bin closure check) | Synthetic three-component fluid recovery ±3 wt% |
| 7 | `core/watson_k.py` + `core/pcsaft_params.py` | Panuganti Table 6 reproduction ±2% |
| 8 | End-to-end pipeline + validation report | All metrics within tolerances on synthetic VTB |
| 9 | Streamlit UI (6 tabs) | App boots, validation tab passes |
| 10 | Final docs + pyproject packaging | `pip install -e .` works |

Detailed pass-gate criteria for each phase are in the corresponding test files (which are written as part of that phase's deliverable). Future-phase content is **out of scope** for the current phase.

═══════════════════════════════════════════════════════════════════════════════
## PHASE 0 — DETAILED INSTRUCTIONS (this is what gets executed when user says "execute Phase 0")
═══════════════════════════════════════════════════════════════════════════════

**Phase 0 is intentionally minimal. Do not expand it. Do not write any code in `core/`, `tests/`, or `tabs/`. Do not write correlations. Do not write tests. Just the scaffold.**

### Step 1 — Confirm working directory

The user has placed you in `~/projects/petrochar/` (or equivalent). Confirm:
```bash
pwd        # should end in /petrochar
ls -la
```

You should see, at minimum:
- `CLAUDE_CODE_PROMPT.md`, `CLAUDE.md`, `CURRENT_STATUS.md`, `README.md` at root
- `data/riazi_reference/` containing 6 CSVs, 9 PNGs, the full PDF, and a README

If `pwd` does not end in `/petrochar`, STOP and ask the user.

### Step 2 — Create directory tree

```bash
mkdir -p core tests tabs docs
touch core/__init__.py tests/__init__.py tabs/__init__.py
```

Empty `__init__.py` files only. Do NOT create `core/correlations.py` or any other module file.

### Step 3 — Write `LICENSE` (MIT, copyright "Piyush Paradkar 2026")

### Step 4 — Write `requirements.txt`
```
numpy>=1.24
scipy>=1.10
pandas>=2.0
matplotlib>=3.7
pytest>=7.4
streamlit>=1.30
openpyxl>=3.1
```

### Step 5 — Write `runtime.txt`: `python-3.11`

### Step 6 — Write `pyproject.toml`
Project name `petrochar`, version `0.1.0`, description per project identity, author Piyush Paradkar, MIT license, dependencies mirror `requirements.txt`. Standard PEP-621 format.

### Step 7 — Write `.gitignore`
Python standard (`__pycache__/`, `*.pyc`, `.pytest_cache/`, `.venv/`, `.streamlit/`, `*.egg-info/`, `dist/`, `build/`). Plus: ignore `*.xlsx` and `*.csv` at repo root, but **not** inside `data/`. Add an exception line `!data/riazi_reference/*.csv`. Also `!data/riazi_reference/*.png` and `!data/riazi_reference/*.pdf`.

### Step 8 — Verify `CLAUDE.md` and `CURRENT_STATUS.md` are present
Both should already be in place from user setup. If not, create from templates at the bottom of this prompt.

### Step 9 — Verify reference materials exist
```bash
ls data/riazi_reference/
```
Expected files:
- `README.md`
- `Riazi_MNL50_full.pdf` (full book, ~12 MB)
- 6 CSVs: `table_4_6_scn_groups.csv`, `table_4_11_example_4_7.csv`, `table_4_13_distribution_coeffs.csv`, `table_4_21_quadrature_points.csv`, `table_4_22_quadrature_3pt.csv`, `table_4_23_lumping_methods.csv`
- 9 PNGs: page extracts as listed in "REFERENCE MATERIALS" section above

If any item is missing, log in `CURRENT_STATUS.md` as a Phase 0 blocker, but proceed with the rest of the scaffold. **Do not invent CSV content or generate PNG files.**

### Step 10 — Initialize git and commit
```bash
git init
git add -A
git commit -m "Phase 0: repository scaffold and Riazi reference materials vendored"
```

### Step 11 — Verify environment
```bash
pip install -r requirements.txt    # if not already in a venv with deps
python -c "import numpy, scipy, pandas, matplotlib, pytest, streamlit, openpyxl; print('imports OK')"
pytest    # should exit 0 with no tests
```

If imports fail, log in `CURRENT_STATUS.md`. If `pytest` returns non-zero with no tests, investigate.

### Step 12 — Update `CURRENT_STATUS.md`
- Mark Phase 0 complete with today's date.
- Set "Current Phase" to "Phase 1 — Shared correlations".
- Append session log entry: `YYYY-MM-DD | Phase 0 complete: scaffold + Riazi materials vendored + git initialised | next: Phase 1`
- Note any blockers.

### Step 13 — Self-check before final commit
```bash
grep -rni "MW_ARO_L\|MW_RES_L\|ARO-L\|ARO_L\|RES-L\|RES_L\|harmonic.mean\|DAO" \
    --include="*.py" --include="*.toml" --include="*.txt" --include="*.cfg" .
```
Expected: zero hits. If any non-md file has a forbidden term, that is a contamination bug — stop, investigate, fix, retry.

### Step 14 — Final commit
```bash
git add -A
git commit -m "Phase 0: complete — status updated, env verified"
```

═══════════════════════════════════════════════════════════════════════════════
## ▼▼▼  END OF PHASE 0. STOP HERE.  ▼▼▼
═══════════════════════════════════════════════════════════════════════════════

**Do not start Phase 1.** Do not write `core/correlations.py`. Do not write any test. Do not write any tab.

If you have time, energy, or capability remaining: stop anyway. The user will start a new session for Phase 1. That session will read `CURRENT_STATUS.md` and proceed with focused context.

Print the contents of `CURRENT_STATUS.md` to the chat as your final action so the user can verify.

═══════════════════════════════════════════════════════════════════════════════
## CURRENT_STATUS.md TEMPLATE (use only if the file is missing)
═══════════════════════════════════════════════════════════════════════════════

```markdown
# petrochar — Current Status

## Project
Standalone Python tool for heavy petroleum fraction characterization (Paper 1).
Target journal: Energy & Fuels or Fuel.

## Current Phase
Phase 0 — Repository scaffold

## Phases Completed
(none yet)

## Decisions Made (frozen)
1. Distribution model: Riazi generalized (Eq. 4.56). Gamma is legacy-only — fails on heavy oils per Riazi p. 178-179.
2. Default quadrature: 5-point + 1 ASP discrete + propane = 7 components.
3. Asphaltenes always discrete: Gonzalez 2007 defaults (m=33, σ=4.3, ε/k=400). Tb=800°C and SG=1.15 are conventions.
4. Aromaticity γ via per-component Watson K, linear clamp 9.5-13.0.
5. SARA: ASP wt% is hard tail constraint; SAT/ARO/RES are K_W-bin closure checks only (NOT tuning targets).
6. No DAO data anywhere in this repository.
7. No reading of any external project (especially `~/projects/pda_pcsaft_tool/`).
8. Validation-first: each phase requires reproducing a Riazi textbook example before proceeding.
9. Reference data: 6 CSVs + 9 PNG page extracts + full Riazi PDF in `data/riazi_reference/`.
10. Tb correlation is REGIME-DEPENDENT: Eq. 2.56 for M ≤ 300, Eq. 2.57 for M > 300. SG correlation similarly two-regime.

## Known Issues / Blockers
(none yet)

## Session Log
(empty — first session pending)
```

═══════════════════════════════════════════════════════════════════════════════
## SUGGESTED INSTRUCTION LINE FOR USER TO PASTE INTO CLAUDE CODE
═══════════════════════════════════════════════════════════════════════════════

> Read `CLAUDE_CODE_PROMPT.md` in full before doing anything. Then execute Phase 0 only — exactly the steps in the section labelled "PHASE 0 — DETAILED INSTRUCTIONS". Stop at the explicit STOP marker. Do not start Phase 1. Print `CURRENT_STATUS.md` at the end so I can review.

═══════════════════════════════════════════════════════════════════════════════
## FINAL NOTES
═══════════════════════════════════════════════════════════════════════════════

- If any instruction in this file conflicts with a docstring or in-code comment: this file wins.
- If user instructions in chat conflict with this file: follow the user, but flag the conflict and propose an update.
- This is a slow, careful build. Do not pattern-match against your training data on similar projects. petrochar's specific architecture is uncommon, and the Riazi textbook examples + page extracts are the only correct check.
