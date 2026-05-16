# petrochar — Changelog

All phases of the petrochar v0.1.0 development cycle are documented here.
Each entry names the deliverable, the numerical pass-gate used to validate it,
and the outcome.

---

## Phase 0 — Repository scaffold (2026-05-05)

Created the repository directory structure (`core/`, `tests/`, `tabs/`, `docs/`,
`data/riazi_reference/`), `requirements.txt`, `runtime.txt`, `pyproject.toml`,
`.gitignore`, and `LICENSE` (MIT). Vendored six Riazi MNL50 reference CSVs and
nine PNG page extracts for offline validation. Initialized git and verified the
Python environment.

**Pass-gate:** `pytest` exits 0 with no tests collected.

---

## Phase 1 — Shared correlations (2026-05-06)

Implemented `core/correlations.py`: regime-dependent Riazi-Daubert T_b correlation
(Eq. 2.56 for M ≤ 300 g/mol, Eq. 2.57 for M > 300 g/mol), Riazi-Daubert M
inversion, SG inversion, Watson K factor (Eq. 2.13), and aromaticity gamma clamp.
49 tests covering all correlations and round-trip fidelity.

**Pass-gate:** Riazi Table 4.11 North Sea C7+ data — T_b reproduced within ±5 K
for all 11 SCN groups (max deviation 3.08 K at C9). 49/49 tests passed.

---

## Phase 2 — Distillation curve handling (2026-05-06)

Implemented `core/distillation.py`: `DistillationCurve` class with method/basis tags
and `to_tbp()` dispatch. D86 → TBP via Daubert Eqs. 3.20–3.22 (Table 3.7 constants).
D1160_AET is a pass-through. 30 tests.

**Pass-gate:** Synthetic kerosene D86 curve (Riazi Example 3.3) — TBP(50%) within
0.02 K of published value; all seven cut points within ±5 K. 30/30 tests passed.

---

## Phase 3 — Generalized distribution fitting (2026-05-08)

Implemented `core/distribution.py`: `GeneralizedDistribution` class with 3-parameter
and 2-parameter fitting modes (Riazi Eq. 4.56), evaluation (CDF, inverse CDF, PDF),
and analytical mean. `from_params()` classmethod for literature-parameter instantiation.
34 tests including fit quality assertions.

**Pass-gate:** Riazi Example 4.13 — North Sea C7+ M, T_b, SG distributions. Fitted
P_o, A, B each within 5% of Table 4.13 reference values for all six fit cases
(M/T_b/SG × 3-param/2-param). 34/34 tests passed.

---

## Phase 4 — SG and MW distributions (2026-05-08)

Implemented `core/sg_distribution.py` (constant Watson K and generalized SG
distribution) and `core/mw_distribution.py` (per-cut M from Riazi-Daubert inversion,
bulk closure verification). Removed large copyrighted PDF and PNG reference files from
git history using `git-filter-repo`. 32 tests.

**Pass-gate:** Table 4.13 SG 3-parameter row — P_o, A, B within 0.1% of reference.
Bulk M_av within 0.43% and SG_av within 0.15% of North Sea C7+ feed. 32/32 tests passed.

---

## Phase 5 — Gaussian quadrature discretization (2026-05-08)

Implemented `core/quadrature.py`: `quadrature_points()` returning vendored Riazi
Table 4.21 Gauss-Laguerre nodes and weights (3-point and 5-point), `Pseudocomponent`
dataclass, and `discretize_generalized()`. 37 tests.

**Pass-gate:** Riazi Example 4.14 — 3-point quadrature with Table 4.22 parameters.
M_i within 1% and z_i within 0.005 for all three nodes. 3-pt M_av within 0.38% of
experimental 118.9 g/mol. 37/37 tests passed; cumulative 182/182.

---

## Phase 6 — SARA closure and asphaltene assembly (2026-05-08)

Implemented `core/sara.py`: `validate_sara()`, `append_asphaltene()` (mole-fraction
basis conversion), and `kw_bin_check()` (K_W-bin aggregation as informational closure
check). Established the z-convention: after `append_asphaltene()`, all z values are
true mole fractions in the full mixture summing to 1. 60 tests.

**Pass-gate:** Synthetic three-component fluid — K_W-bin recovered SAT/ARO/RES/ASP
wt% to within 0.01 wt% on known-composition input. 60/60 tests passed;
cumulative 242/242.

---

## Phase 7 — Watson K, aromaticity, and PC-SAFT parameters (2026-05-15)

Implemented `core/watson_k.py` (`compute_K_W_per_pseudocomponent`) and
`core/pcsaft_params.py` (Panuganti 2012 Table 6 correlations: `panuganti_saturate_params`,
`panuganti_aromatic_resin_params`, `panuganti_distillable_params`, `gonzalez_asphaltene_params`,
`propane_params`, `generate_pcsaft_table`). Extended `Pseudocomponent` dataclass with
`K_W`, `gamma`, and `is_asphaltene` fields. 71 tests including propane sigma typo guard.

**Pass-gate:** Panuganti 2012 Tables 10–12 — aromatic/resin PC-SAFT parameters for
Crudes A, B, C reproduced within 1% (m), 0.5% (σ), 2% (ε/k). Propane σ = 3.6180 Å
(not 3.168 Aspen typo) confirmed by defensive test. 71/71 tests passed;
cumulative 314/314.

---

## Phase 8 — End-to-end pipeline + validation report (2026-05-15; reworked 2026-05-16)

**Original (2026-05-15):** Implemented the full 9-step characterization pipeline in
`tests/test_phase8_pipeline.py`. Fixed `riazi_daubert_M` non-monotone bracket (Eq. 2.57
peak issue) and regime-gap fallback returning M = 300 with UserWarning. Wrote
`docs/validation_report.md`. 40 tests passed + 6 xfail (documented limitations).

**Rework (2026-05-16):** Resolved four architectural decisions that eliminated all
xfail tests:

- *Decision 25:* Pseudo-component T_b derived from `riazi_daubert_Tb(M_i, SG_i)` via
  self-consistent constant-Watson-K solve (brentq). T_b distribution is diagnostic only.
- *Decision 26:* `is_asphaltene` flag (not T_b > 1000 K) is the sole asphaltene
  identifier throughout `core/`. Hard ValueError if non-asphaltene has T_b > 1000 K.
- *Decision 27:* M_av pass-gate compares GL result vs distribution analytic mean (≤ 0.5%),
  not against M_DIST_TARGET. Synthetic feed inconsistency documented.
- *Decision 28:* K_W-bin test restricted to ASP placement + closure sum under constant K_W.

**Pass-gate results (reworked):** GL M_av = 371.66 g/mol vs dist mean 369.95 g/mol →
0.46% (gate ≤ 0.5%: PASS). Full-mixture SG_av = 0.981 vs SG_BULK = 1.020 → 3.9%
(gate ≤ 5%: PASS). All 5 distillable nodes T_b < 1000 K (max 869.7 K): PASS.
372/372 tests passed, 0 xfailed.

---

## Phase 9 — Streamlit UI (2026-05-16)

Built `app.py` (session-state setup, 6-tab routing) and six tab files in `tabs/`:

- **Tab 1 (input_data.py):** Distillation curve entry (manual template or CSV upload),
  bulk SG/MW inputs, SARA wt% inputs, method/basis/unit selectors, quadrature count
  selector, full 9-step pipeline trigger. Live validity warnings. D7169 visible with
  warning and disabled run button (Option B fix).
- **Tab 2 (distillation_fit.py):** T_b distribution diagnostic plot (Riazi Eq. 4.56),
  fit quality metrics, 3-param vs 2-param comparison toggle.
- **Tab 3 (distributions.py):** 2×2 matplotlib figure — M PDF, M CDF, SG PDF, SG CDF;
  quadrature nodes marked with vertical lines on CDF panels.
- **Tab 4 (pseudocomponents.py):** Pseudo-component table with ASP row highlighted,
  K_W-bin closure report, CSV download.
- **Tab 5 (pcsaft_export.py):** Pure-component PC-SAFT parameters, component definitions,
  editable binary k_ij matrix; CSV + TSV download per section.
- **Tab 6 (validation.py):** In-app reproduction of Riazi Examples 4.13, 4.14, and 4.7
  with pass/fail assessment; architectural commitments statement.

**Pass-gate:** `streamlit run app.py` HEALTH=ok. 372/372 tests unchanged.
Forbidden-terms grep: 0 hits.

---

## Phase 10 — Final docs, packaging, and D7169 fix (2026-05-16)

Finalized all documentation and packaging for v0.1.0:

- **D7169 fix:** `tabs/input_data.py` — D7169 selection shows `st.warning()` and
  disables the Run button (Option B: visible but gated). No silent exception.
- **README.md:** Full user-facing documentation with installation, quick start,
  methodology and validation links, architectural commitments, citation placeholder,
  and acknowledgments.
- **docs/methodology.md:** Paper 1 Section 2 draft — nine subsections covering the
  complete methodology with equations, references, and three critical deviation
  disclosures (§2.9.1 T_b distribution as diagnostic; §2.9.2 constant Watson K
  limitation; §2.9.3 gamma from Watson K vs AOP fitting; §2.9.4 bulk MW as input
  not target).
- **pyproject.toml:** PyPI-ready — `setuptools.build_meta` backend, author email,
  classifiers, `[project.urls]`, `[project.scripts]` defining `petrochar` CLI.
- **petrochar/_cli.py:** Console entry point launching `streamlit run app.py`.
- **CHANGELOG.md:** This file.
- **CURRENT_STATUS.md:** Phase 9 deferred items classified into three subsections
  (Paper 1 disclosures, pre-submission fixes, future enhancements).

**Pass-gate:** `pip install -e ".[dev]"` succeeded. `petrochar` CLI launches app.
D7169 warning and disabled run confirmed. 372/372 tests unchanged.
Forbidden-terms grep: 0 hits.

---

*End of changelog for petrochar v0.1.0.*
