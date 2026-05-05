# petrochar — Current Status

## Project
Standalone Python tool for heavy petroleum fraction characterization. Methodology contribution for **Paper 1**: "Continuous-distribution characterization of heavy petroleum fractions for PC-SAFT phase-equilibrium modeling from routine refinery data". Target: *Energy & Fuels* or *Fuel*.

This is a separate repository from any other PC-SAFT project. Do not import, read, or reference code from outside this directory. See `CLAUDE.md` for forbidden paths.

## Current Phase
Phase 0 — Repository scaffold (not yet started)

## Phases Completed
(none yet)

## Decisions Made (frozen)

1. **Distribution model:** Riazi generalized (Eq. 4.56) is default. Gamma is legacy-only — fails on heavy oils per Riazi p. 178-179.
2. **Quadrature:** 5-point default + 1 ASP discrete + propane = 7 components total. 3-point selectable.
3. **Asphaltenes:** always discrete, Gonzalez 2007 defaults (m=33, σ=4.3 Å, ε/k=400 K). Tb=800°C and SG=1.15 are numerical conventions.
4. **Aromaticity γ:** per-component Watson K, linear clamp `γ = clamp((13.0 − K_W) / (13.0 − 9.5), 0, 1)`.
5. **SARA:** ASP wt% used as hard tail constraint. SAT/ARO/RES used as K_W-bin closure check only — never tuned to.
6. **No DAO data anywhere** in this repository.
7. **No reading external projects** (specifically `~/projects/pda_pcsaft_tool/` or any file outside this repo).
8. **Validation-first:** each phase requires reproducing a Riazi textbook example before proceeding.
9. **Repository name:** `petrochar`.
10. **Tb correlation is REGIME-DEPENDENT:** Eq. 2.56 for M ≤ 300, Eq. 2.57 for M > 300 (recommended for M > 300). Same for SG: Eq. 2.59 and Eq. 2.60. **Critical correction; verified against PDF page extract.**
11. **Reference materials vendored:** full Riazi MNL50 PDF, 6 textbook tables as CSVs (4.6, 4.11, 4.13, 4.21, 4.22, 4.23), 9 PNG page extracts at 200 DPI for visual verification of equations during implementation.

## Reference Materials Inventory

`data/riazi_reference/` contains:

**Full reference:**
- `Riazi_MNL50_full.pdf` — complete book, 427 pages, ~12 MB

**Page extracts (PNG, 200 DPI):**
- `page_054_watson_k_eq_2_13.png` — Eq. 2.13 (Watson K). Phase 1.
- `page_077_riazi_daubert_Tb.png` — Eqs. 2.56, 2.57 (Tb), 2.59, 2.60 (SG). **Phase 1 critical reference.**
- `page_078_riazi_daubert_Tb_continued.png` — continuation. Phase 1.
- `page_122_daubert_d86_to_tbp.png` — Section 3.2.2.2, Eqs. 3.20-3.22, Table 3.7. Phase 2.
- `page_123_daubert_d86_to_tbp.png` — continuation.
- `page_124_daubert_d86_to_tbp.png` — continuation, Eqs. 3.23-3.25 (SD-to-TBP).
- `page_181_table_4_5_4_6.png` — Table 4.5 SCN coefficients + Table 4.6 SCN data. Phase 4 fallback.
- `page_182_table_4_5_4_6.png` — continuation.
- `page_183_table_4_5_4_6.png` — continuation.

**CSVs (textbook tables for validation tests):**
- `table_4_6_scn_groups.csv` — SCN group properties C6-C50, fallback resource.
- `table_4_11_example_4_7.csv` — North Sea gas condensate C7+ data. Phases 1, 3, 4 pass-gate.
- `table_4_13_distribution_coeffs.csv` — Distribution coefficients (8 rows: 4 three-param + 4 two-param). Phases 3, 4.
- `table_4_21_quadrature_points.csv` — Gaussian quadrature roots/weights (3-point and 5-point). Phase 5.
- `table_4_22_quadrature_3pt.csv` — Expected output for Example 4.14. Phase 5 pass-gate.
- `table_4_23_lumping_methods.csv` — Five-component lumping, Example 4.15. Supplementary.

## Known Issues / Blockers

(none — all reference materials vendored at project setup)

## Session Log

(empty — first session pending)
