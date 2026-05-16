# petrochar — Current Status

## Project
Standalone Python tool for heavy petroleum fraction characterization. Methodology contribution for **Paper 1**: "Continuous-distribution characterization of heavy petroleum fractions for PC-SAFT phase-equilibrium modeling from routine refinery data". Target: *Energy & Fuels* or *Fuel*.

This is a separate repository from any other PC-SAFT project. Do not import, read, or reference code from outside this directory. See `CLAUDE.md` for forbidden paths.

## Current Phase
COMPLETE — all phases 0-10 done. Next: Paper 1 drafting.

## Phases Completed
- ✅ Phase 0 — Repository scaffold (2026-05-05)
- ✅ Phase 1 — Shared correlations (`core/correlations.py`) (2026-05-06)
- ✅ Phase 2 — Distillation curve conversion (`core/distillation.py`) (2026-05-06)
- ✅ Phase 3 — Distribution fitting (`core/distribution.py`) (2026-05-08)
- ✅ Phase 4 — SG and MW distributions (`core/sg_distribution.py`, `core/mw_distribution.py`) (2026-05-08)
- ✅ Phase 5 — Gaussian quadrature discretization (`core/quadrature.py`) (2026-05-08)
- ✅ Phase 6 — SARA closure check and asphaltene assembly (`core/sara.py`) (2026-05-08)
- ✅ Phase 7 — Watson K → γ + PC-SAFT parameters (`core/watson_k.py`, `core/pcsaft_params.py`) (2026-05-15)
- ✅ Phase 8 — End-to-end pipeline integration test + validation report (reworked 2026-05-16)
- ✅ Phase 9 — Streamlit UI: 6 tabs, app.py + tabs/*.py (2026-05-16)
- ✅ Phase 10 — Final docs, packaging, D7169 fix (2026-05-16)

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
17. **`from_params` is public API:** GeneralizedDistribution.from_params() is documented as a public classmethod for instantiating from external parameters (published tables, external sources) without re-fitting. Not just a test utility.
18. **Pseudocomponent.z is always a true mole fraction in the full mixture:** After sara.append_asphaltene() is called, ALL z values (distillable + ASP) are true mole fractions summing to 1. append_asphaltene performs the basis conversion internally: z_i_true = z_i * n_dist/n_tot; z_asp = n_asp/n_tot where n_dist = (1-f_asp)/M_dist_av, n_asp = f_asp/M_asp. Phase 7 can use z directly without any further conversion.
19. **ASP Tb_K = 1073.15 K (800°C) is a numerical convention:** Asphaltenes do not boil. Value chosen to sort above all distillable pseudo-components. Asphaltene component identified in kw_bin_check by Tb_K > 1000 K threshold.
20. **K_W bin thresholds are parameterizable conventions:** Default SAT >= 12.0, 11.0 <= ARO < 12.0, RES < 11.0 (Riazi p. 75). User can override via kw_sat/kw_aro parameters in kw_bin_check.
21. **γ-from-Watson-K is a deliberate architectural deviation from Panuganti 2012 — must be disclosed in Paper 1 §2:** Panuganti defines γ as a free parameter fitted simultaneously to crude density and bubble-point pressure (AOP data) for each crude (Panuganti 2012, p. 4). That fitting procedure allows γ to absorb crude-specific aromaticity information that Watson K cannot capture — in particular, two crudes with identical boiling-point distributions but different aromatic content will have different optimal γ but identical K_W. petrochar abandons that fitting and instead uses a deterministic per-component mapping: γ = clamp((13.0 − K_W) / (13.0 − 9.5), 0, 1), derived from Watson K (Watson & Nelson 1933). This eliminates the need for AOP data entirely — the characterization runs from routine refinery assay data alone — but at the cost of losing crude-level density-matching capability: the PC-SAFT parameters will not reproduce measured crude densities as accurately as Panuganti's AOP-fitted values, and the model cannot be calibrated to a specific crude's phase behaviour without external AOP data. This trade-off is the defining methodological choice of petrochar and must be stated plainly in Paper 1 §2. The single γ-interpolated A+R correlation (Panuganti Table 6, rows 4–9) is used for ALL distillable components; the separate Saturates form (rows 1–3) is retained only in `panuganti_saturate_params` for reference and independent testing.
22. **Panuganti reference data vendored:** `data/panuganti_2012/` contains four CSVs (table_5_light_component_pcsaft_params.csv, table_6_correlations.csv, tables_10_11_12_crude_pcsaft.csv, README.md). The PDF is gitignored (copyrighted). CSVs are tracked via `!data/panuganti_2012/*.csv` exception in .gitignore.
23. **Propane σ = 3.6180 Å (NOT 3.168):** The value 3.168 appears in some Aspen Plus built-in databases as a transcription error (14% wrong). The correct value from Gross & Sadowski 2001, reproduced in Panuganti Table 5, is 3.618 Å. This is explicitly guarded by `test_sigma_is_not_aspen_typo` in the Phase 7 test suite.
24. **riazi_daubert_M non-monotone bracket fix (Phase 8):** Eq. 2.57 Tb(M) peaks at M_peak = 0.5369/(7.5152e-4·SG - 1.6514e-4). Bracket now uses ascending branch [300.01, M_peak×0.9999] only. Regime-gap fallback returns M=300 with UserWarning when target Tb falls in the discontinuity between Eq.2.56(M=300) and Eq.2.57(M=300+). Both behaviours tested in test_phase8_pipeline.py.
25. **Self-consistent (M, K_W) → (Tb, SG) solve (Phase 8 rework, 2026-05-16):** Pseudo-component Tb_i is derived from riazi_daubert_Tb(M_i, SG_i) directly, NOT by evaluating the Tb distribution at quadrature nodes. Under constant Watson K: SG_i = (1.8·Tb_i)^(1/3)/K_W_bulk. The coupled system is solved per node via brentq on [300, 990] K. The Tb distribution is diagnostic only — it characterises the feed but does not generate pseudo-component properties. Result: all 5 distillable GL nodes have Tb ∈ [629, 870] K (< 1000 K); the Tb>1000 K reclassification path is eliminated at source.
26. **is_asphaltene flag replaces Tb>1000 K threshold everywhere (Phase 8 rework, 2026-05-16):** Asphaltene identity is carried by Pseudocomponent.is_asphaltene=True, set exclusively by sara.append_asphaltene(). No code in core/ uses Tb>1000 K to identify asphaltenes. A hard ValueError is raised if a non-asphaltene component has Tb>1000 K. Decision 19 (ASP Tb=1073.15 K as sorting convention) is retained but the identification logic uses the flag, not the threshold.
27. **M_av pass-gate compares GL result vs distribution analytic mean (Phase 8 rework, 2026-05-16):** The gate is GL_M_av vs m_dist.average() ≤ 0.5% (quadrature accuracy criterion). M_DIST_TARGET=563.6 g/mol is reported as a diagnostic but is not a gate — the synthetic test feed is internally inconsistent (D1160 endpoint at xc=0.95 implies M~527; bulk MW=700 requires ~4200 g/mol in the 5% tail, which is unphysical). Verified: GL_M_av=371.66, analytic mean=369.95 → 0.46% deviation (PASS).
28. **K_W-bin test restricted to ASP placement and closure sum (Phase 8 rework, 2026-05-16):** Under constant Watson K all distillable components share K_W_bulk and fall in the same bin — SAT/ARO/RES wt% tests on the constant-K_W synthetic feed have no discriminating power. Tests check: (a) ASP in ASP class; (b) sum = 100 wt%; (c) no NaN deltas. SAT/ARO/RES bin deviation tests removed (were xfail).

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
- **Phase 8 informational notes (not blockers, not xfail):**
  - Tb distribution fit: 3-param free-B fit gives P_o=424.3 K (onset below IBP=553.15 K). This is expected for the free-B form. Invariant enforced: 0 < P_o < Tb_min.
  - Constant Watson K assigns K_W_bulk to ALL distillable components (K_W=11.331). This means all fall in the ARO K_W bin for this feed — physically reasonable for a VTB fraction. SAT/RES bins are empty by construction, not by error.
  - Regime gap at Tb=663 K, SG=0.936 (xc=0.30) returns M=300 as boundary estimate with UserWarning. Reported as diagnostic in pipeline output.
  - Synthetic test feed inconsistency: D1160 endpoint at xc=0.95 (Tb~663 K → M~527) combined with bulk MW=700 implies ~4200 g/mol in the 5% tail. Unphysical. M_DIST_TARGET=563.6 is a diagnostic, not a gate (Decision 27).

## Phase 9 deferred items (scope-creep rule enforced)

Classified into three subsections. "Paper 1 disclosures" must appear in the paper
before submission. "Pre-submission fixes" are bugs addressed in Phase 10.
"Future enhancements" do not affect publishability.

### Paper 1 disclosures required

- **Constant Watson K only in current UI:** The generalized SG distribution
  (Riazi Eq. 4.56 fitted to volume-basis SG data) is architecturally supported by
  `core/sg_distribution.py` but is not exposed in the current pipeline because it
  requires per-cut SG measurements not routinely available for heavy residue streams.
  The current implementation uses the constant Watson K method exclusively.
  Paper 1 §2 disclosure: state this as a deliberate scope decision with stated reason,
  not a gap. Under constant Watson K, all distillable components share K_W_bulk —
  SAT/ARO/RES K_W-bin classification is degenerate by construction; users who need
  finer SARA classification should supply per-cut SG data and use generalized mode
  (future version).
- **Weight-basis distillation input only:** `to_weight_basis()` raises ValueError for
  volume and mole basis inputs (Phase 2 design — SG distribution needed for conversion).
  The `DistillationCurve` class supports the basis tag architecturally, but the
  volume→weight conversion path is not wired through the current pipeline.
  Paper 1 §2 disclosure: "The current pipeline requires weight-basis distillation data.
  Volume-basis input is supported by the underlying DistillationCurve class but requires
  density-distribution input not exposed in the current pipeline."

### Pre-submission fixes (addressed in Phase 10)

- **D7169 option raises ValueError silently:** Tab 1 shows D7169 as a method option
  but triggers an unhandled exception on pipeline run (Phase 2 D7169→TBP conversion
  not implemented). Fix: keep D7169 visible (signals architectural intent) but show
  st.warning() and disable the Run button when selected. Fixed in Phase 10.

### Future enhancements (no Paper 1 impact)

- **PDF/SVG figure download buttons:** Tabs 2-3 render matplotlib figures as PNG via
  st.image(). For Paper 1 figure generation, write a standalone script calling
  core/ functions directly with matplotlib.pyplot.savefig(). Add UI download button
  only after knowing which figures make the paper.
- **Component name customization in Tab 5:** PC1, PC2, ..., ASP, C3 are
  auto-generated. User-editable names for Aspen Plus compatibility is a UX improvement.
- **k_ij calibration guidance:** Default propane k_ij=0.010 is a literature starting
  point. A reference table from Panuganti 2012 Tables 11-12 fitted values for Crudes
  A/B/C would help users calibrate for specific feeds. Add to docs, not to UI.

## Session Log

2026-05-05 | Phase 0 complete: scaffold + Riazi materials vendored + git initialised | next: Phase 1
2026-05-06 | Phase 1 complete: correlations.py 6 functions, test_phase1 49/49 pass, max Tb dev 3.08 K | next: Phase 2
2026-05-06 | Phase 2 complete: distillation.py DistillationCurve class + Daubert D86→TBP, test_phase2 30/30 pass, kerosene Example 3.3 all +-5 K, 50% point within 0.02 K | next: Phase 3
2026-05-08 | Copyright fix: Riazi PDF + 9 PNGs removed from git tracking, .gitignore updated, README updated with acquisition instructions | Phase 3 complete: distribution.py GeneralizedDistribution class, test_phase3 34/34 pass, T_o/A/B all within +-5% of Table 4.13 for both 3-param and 2-param Tb fits | next: Phase 4
2026-05-08 | git history purge: filter-repo removed PDF and all page_*.png from all commits; force-pushed; all history SHAs rewritten | Phase 4 complete: sg_distribution.py + mw_distribution.py, test_phase4 32/32 pass, SG 3-param <0.1% of Table 4.13, M_av 0.43% of 118.9, SG_av 0.15% of 0.7597 | next: Phase 5
2026-05-08 | Phase 5 complete: quadrature.py (quadrature_points, Pseudocomponent, discretize_generalized), distribution.py +from_params classmethod, test_phase5 37/37 pass, 182/182 total pass, 3-pt M_i all within 1% of Table 4.22, M_av 0.38% of 118.9 | next: Phase 6
2026-05-08 | Phase 6 complete: sara.py (validate_sara, append_asphaltene, kw_bin_check), test_phase6 57/57 pass, 239/239 total pass, K_W binning recovers SAT/ARO/RES/ASP wt% to <0.01 wt% on synthetic input, no retuning logic | next: Phase 7
2026-05-08 | Phase 6 z-convention fix: append_asphaltene now converts all z to true full-mixture mole fractions (uniform basis), kw_bin_check updated to use z_i*M_i/M_mix uniformly; test_phase6 60/60 pass, 242/242 total pass | next: Phase 7
2026-05-15 | Phase 7 complete: watson_k.py (compute_K_W_per_pseudocomponent), pcsaft_params.py (panuganti_saturate_params, panuganti_aromatic_resin_params, panuganti_distillable_params, gonzalez_asphaltene_params, propane_params, generate_pcsaft_table), quadrature.py extended with K_W/gamma fields; test_phase7 71/71 pass, 314/314 total pass; A+R params within 1%/0.5%/2% of Panuganti Tables 10-12 for Crudes A/B/C | next: Phase 8
2026-05-15 | Phase 8 complete: riazi_daubert_M non-monotone bracket fix + regime-gap fallback; tests/test_phase8_pipeline.py 40/40 pass + 6 xfail (documented limitations); docs/validation_report.md written (draft for Paper 1 Section 3); total 354 passed, 6 xfailed | next: Phase 9
2026-05-16 | Phase 8 rework: self-consistent (M,K_W)→(Tb,SG) solve (Decision 25); is_asphaltene flag replaces Tb>1000K threshold everywhere (Decision 26); M_av gate vs analytic mean not M_DIST_TARGET (Decision 27); K_W-bin tests restricted to ASP+closure (Decision 28); all 6 xfail removed; validation_report.md updated; total 372 passed, 0 xfailed | next: Phase 9
2026-05-16 | Phase 9 complete: app.py + 6 tab files; boot health=ok; 372 passed (unchanged); forbidden-terms 0 hits | next: Phase 10
2026-05-16 | Phase 10 complete: README.md, docs/methodology.md, CHANGELOG.md, pyproject.toml PyPI-ready, petrochar/_cli.py, D7169 fix; pip install -e . OK; 372 passed; 0 forbidden-term hits | IMPLEMENTATION COMPLETE
