# petrochar — Current Status

## Project
Standalone Python tool for heavy petroleum fraction characterization. Methodology contribution for **Paper 1**: "Continuous-distribution characterization of heavy petroleum fractions for PC-SAFT phase-equilibrium modeling from routine refinery data". Target: *Energy & Fuels* or *Fuel*.

This is a separate repository from any other PC-SAFT project. Do not import, read, or reference code from outside this directory. See `CLAUDE.md` for forbidden paths.

## Current Phase
Phase 5 — Gaussian quadrature discretization (`core/quadrature.py`)

## Phases Completed
- ✅ Phase 0 — Repository scaffold (2026-05-05)
- ✅ Phase 1 — Shared correlations (`core/correlations.py`) (2026-05-06)
- ✅ Phase 2 — Distillation curve conversion (`core/distillation.py`) (2026-05-06)
- ✅ Phase 3 — Distribution fitting (`core/distribution.py`) (2026-05-08)
- ✅ Phase 4 — SG and MW distributions (`core/sg_distribution.py`, `core/mw_distribution.py`) (2026-05-08)

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
10. **Tb correlation is REGIME-DEPENDENT:** Eq. 2.56 for M ≤ 300, Eq. 2.57 for M > 300 (recommended for M > 300). Critical correction; verified against page_077 PNG.
11. **SG correlation method:** `riazi_daubert_SG` uses numerical inversion of `riazi_daubert_Tb`, NOT direct use of Eqs. 2.59/2.60. Reason: Eqs. 2.59 (inputs: Tb, I) and 2.60 (inputs: M, I) both require refractivity index I — confirmed by viewing page_077 and page_078 PNGs. No I-free direct SG correlation exists in §2.4.3.1. Implication for accuracy: SG error is bounded by Tb error × ∂SG/∂Tb (typically small). Phase 4 SG distribution does not rely on `riazi_daubert_SG` directly (it uses cumulative-fraction data and bulk closure), so this is not on the critical path.
12. **SG bracketing for inversion:** Lower bracket = SG_min = -f/(c+d·M), analytically derived from d(ln Tb)/d(SG) = 0. The Eq. 2.56/2.57 forms are non-monotone in SG and have a minimum near SG~0.61-0.74; plain brentq on [0.40, 1.30] fails for all typical petroleum SG values.
13. **Reference materials:** 6 textbook table CSVs (4.6, 4.11, 4.13, 4.21, 4.22, 4.23) tracked in repo. Riazi MNL50 PDF and 9 PNG page extracts are NOT tracked (removed from all git history 2026-05-08 — copyrighted ASTM material). Must be obtained separately and placed in `data/riazi_reference/` locally; see README in that directory.
14. **Daubert Eq. 3.20 exponent:** `1.0258` in Kelvin (NOT 0.9217 in °F — early PNG reading was wrong). Verified numerically: D86(50%)=479.85 K → TBP(50%)=483.75 K = 210.6°C, exactly matching Riazi MNL50 Table 3.8 kerosene Example 3.3. OCR rendered the exponent as `1~` (decimal truncated); PDF page 122 text confirms temperatures in Kelvin.
15. **D1160_AET is a pass-through:** method tag changes to TBP, temperatures unchanged. Atmospheric equivalent temperature is already on a TBP-like basis.
16. **to_weight_basis deferred to Phase 4:** raises ValueError for any non-weight basis input in Phase 2. SG distribution required for volume→weight; not available until Phase 4.

## Reference Materials Inventory

`data/riazi_reference/` contains:

**NOT in git (copyrighted ASTM material — purged from all history 2026-05-08):**
- `Riazi_MNL50_full.pdf` — obtain from ASTM; place locally before implementing new phases.
- `page_054_watson_k_eq_2_13.png` ... `page_183_table_4_5_4_6.png` — 9 PNG page extracts; generate with `pdftoppm -r 200 -png`.

**CSVs tracked in git (factual tabular data — textbook tables for validation tests):**
- `table_4_6_scn_groups.csv` — SCN group properties C6-C50, fallback resource.
- `table_4_11_example_4_7.csv` — North Sea gas condensate C7+ data. Phases 1, 3, 4 pass-gate.
- `table_4_13_distribution_coeffs.csv` — Distribution coefficients (8 rows: 4 three-param + 4 two-param). Phases 3, 4.
- `table_4_21_quadrature_points.csv` — Gaussian quadrature roots/weights (3-point and 5-point). Phase 5.
- `table_4_22_quadrature_3pt.csv` — Expected output for Example 4.14. Phase 5 pass-gate.
- `table_4_23_lumping_methods.csv` — Five-component lumping, Example 4.15. Supplementary.

## Known Issues / Blockers

- Python environment uses Python 3.12 (runtime.txt specifies 3.11 — no functional impact; note for deployment).
- `python -m pytest` exits with code 5 (no tests collected) — expected in pytest 7+ when no tests exist; not a failure. Use `python -m pytest` not bare `pytest` (scripts not on PATH).
- Riazi-Daubert Eq. 2.56 max deviation vs Table 4.11 is **3.08 K** (C9 row), not ±2 K as the spec states. Spec tolerance is tighter than the correlation's actual %AAD (~3.5% for M<300). Pass-gate uses ±5 K. Eq. 2.57 deviations vs Table 4.6 SCN rows are 17–20 K (consistent with stated 4.7% AAD); test gate set to ±25 K.
- `riazi_daubert_SG` uses smart bracketing: lower bound is analytically computed SG_min = -f/(c+d·M) (minimum of Tb vs SG), not 0.40. Required because Tb(M,SG) is non-monotone; the function has a minimum in SG at ~0.61–0.74 for refinery-relevant M. Plain brentq on [0.40, 1.30] fails for all typical petroleum SG values.
- Windows cp1252 encoding: test print strings must use ASCII only — no Unicode delta (U+0394) or similar.

## Session Log

2026-05-05 | Phase 0 complete: scaffold + Riazi materials vendored + git initialised | next: Phase 1
2026-05-06 | Phase 1 complete: correlations.py 6 functions, test_phase1 49/49 pass, max Tb dev 3.08 K | next: Phase 2
2026-05-06 | Phase 2 complete: distillation.py DistillationCurve class + Daubert D86→TBP, test_phase2 30/30 pass, kerosene Example 3.3 all +-5 K, 50% point within 0.02 K | next: Phase 3
2026-05-08 | Copyright fix: Riazi PDF + 9 PNGs removed from git tracking, .gitignore updated, README updated with acquisition instructions | Phase 3 complete: distribution.py GeneralizedDistribution class, test_phase3 34/34 pass, T_o/A/B all within +-5% of Table 4.13 for both 3-param and 2-param Tb fits | next: Phase 4
2026-05-08 | git history purge: filter-repo removed PDF and all page_*.png from all commits; force-pushed; all history SHAs rewritten | Phase 4 complete: sg_distribution.py + mw_distribution.py, test_phase4 32/32 pass, SG 3-param <0.1% of Table 4.13, M_av 0.43% of 118.9, SG_av 0.15% of 0.7597 | next: Phase 5
