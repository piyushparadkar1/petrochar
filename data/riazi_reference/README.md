# Riazi MNL50 Reference Materials

This folder contains reference materials transcribed and rasterized from:

> Riazi, M. R. (2005). *Characterization and Properties of Petroleum Fractions*. ASTM Manual MNL50. ASTM International, West Conshohocken, PA. ISBN 0-8031-3361-8.

These materials enable petrochar's validation-first development protocol: every numerical correlation must be verified against the textbook before being used in production code.

## Full reference

- **`Riazi_MNL50_full.pdf`** — complete 427-page book. Use only when a page outside the vendored extracts is needed; in that case, rasterize the relevant page yourself with `pdftoppm` and add it to this folder.

## Page extracts (PNG, 200 DPI)

These are the equations and tables most frequently referenced during petrochar implementation:

| File | Section / Content | Phase |
|------|-------------------|-------|
| `page_054_watson_k_eq_2_13.png` | §2.1.15 Watson K factor (Eq. 2.13) | 1 |
| `page_077_riazi_daubert_Tb.png` | §2.4.2.1 Eqs. 2.56, 2.57 (Tb) and §2.4.3.1 Eqs. 2.59, 2.60 (SG) | 1 |
| `page_078_riazi_daubert_Tb_continued.png` | continuation of above | 1 |
| `page_122_daubert_d86_to_tbp.png` | §3.2.2.2 Daubert D86→TBP, Eqs. 3.20-3.22, Table 3.7 | 2 |
| `page_123_daubert_d86_to_tbp.png` | continuation | 2 |
| `page_124_daubert_d86_to_tbp.png` | continuation, Eqs. 3.23-3.25 (SD→TBP) | 2 |
| `page_181_table_4_5_4_6.png` | Table 4.5 SCN coefficients + Table 4.6 SCN data | 4 |
| `page_182_table_4_5_4_6.png` | continuation | 4 |
| `page_183_table_4_5_4_6.png` | continuation | 4 |

**Workflow rule:** before implementing any correlation in `core/correlations.py` or related modules, view the corresponding PNG with the `view` tool and verify the equation form. Do not implement equations from memory.

## CSVs (for automated test fixtures)

Each CSV begins with comment rows (lines starting with `#`) documenting source citation and bulk-property context.

| File | Content | Used in |
|------|---------|---------|
| `table_4_6_scn_groups.csv` | SCN group properties C6-C50 (M, Tb, SG, n_20, d_20, Tc, Pc, dc, Zc, ω, σ, δ) | Phase 4 fallback |
| `table_4_11_example_4_7.csv` | North Sea gas condensate C7+, 12 rows. Bulk M_7+ = 118.9, SG_7+ = 0.7597 | Phase 1, 3, 4 pass-gates |
| `table_4_13_distribution_coeffs.csv` | Distribution coefficients (P_o, A, B) for M, Tb, SG. 8 rows (4 three-param + 4 two-param) | Phase 3, 4 pass-gates |
| `table_4_21_quadrature_points.csv` | Gaussian quadrature roots/weights for N=3 and N=5 | Phase 5 (used in implementation, not just testing) |
| `table_4_22_quadrature_3pt.csv` | Expected 3-point quadrature output for Example 4.14 (z_i, M_i) | Phase 5 pass-gate |
| `table_4_23_lumping_methods.csv` | 5-component lumping by two methods, Example 4.15 | Supplementary |

Read CSVs with `pandas.read_csv(path, comment='#')`.

## Tolerances for pass-gate tests

Following Riazi MNL50 conventions:
- Tb: ±2 K typical, ±5 K worst case
- M: ±1.5% typical
- SG: ±0.005 typical
- Distribution coefficients (P_o, A, B): ±5% on each
- Quadrature output (z_i, M_i): ±0.005 on z_i, ±1% on M_i

## Critical implementation note

Riazi MNL50 §2.4.2.1 provides **two** Riazi-Daubert correlations for Tb:
- **Eq. 2.56** valid for M = 70-300
- **Eq. 2.57** valid for M = 300-700

VTB pseudo-components span both ranges. The `riazi_daubert_Tb` function in petrochar must select correlation based on input M. Same for the inverse (`riazi_daubert_M`) and for SG (Eqs. 2.59 vs 2.60). See `page_077_riazi_daubert_Tb.png` for both equation forms.

## Fair use

These materials are reproduced strictly for software validation purposes. The full PDF is vendored for offline reference during development. If petrochar is published, the reference materials remain in the repository under standard academic fair use; users of petrochar should cite Riazi MNL50 in any work derived from this tool.
