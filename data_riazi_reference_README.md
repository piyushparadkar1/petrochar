# Riazi MNL50 Reference Data

This folder contains numerical data tables transcribed from:

> Riazi, M. R. (2005). *Characterization and Properties of Petroleum Fractions*. ASTM Manual MNL50. ASTM International, West Conshohocken, PA. ISBN 0-8031-3361-8.

These tables are reproduced here strictly for **software validation purposes** — to enable automated regression tests that confirm petrochar reproduces the textbook's worked examples to within published tolerances.

## Files

- `table_4_6_scn_groups.csv` — SCN group properties C5-C45 (M, Tb, SG, n_20). From Riazi MNL50 Table 4.6.
- `table_4_11_example_4_7.csv` — North Sea gas condensate C7+ fraction characterization data. From Riazi MNL50 Table 4.11 (Example 4.7).
- `table_4_13_distribution_coeffs.csv` — Distribution coefficients for Example 4.7 system (P_o, A, B for M, Tb, SG with various basis options). From Riazi MNL50 Table 4.13.
- `table_4_22_quadrature_3pt.csv` — Three-point Gaussian quadrature output for Example 4.14 (y_i, w_i=z_i, M_i). From Riazi MNL50 Table 4.22.
- `table_4_23_lumping_methods.csv` — Five-component lumping for Example 4.15 (mole fraction, weight fraction, M_i, SG_i for two methods). From Riazi MNL50 Table 4.23.

## CSV format conventions

Each CSV begins with a comment row starting `# SOURCE:` citing the table.

Column names include explicit units in parentheses: e.g., `M_g_per_mol`, `Tb_K`, `SG_dimensionless`.

No proprietary data, no plant data, no client-specific information appears in any file in this folder. All values are textbook reproductions for validation testing.

## Usage in tests

The `tests/test_phase*` modules load these CSVs to obtain reference inputs and expected outputs for pass-gate validation. Test tolerances follow the conventions stated in Riazi MNL50:
- Tb: ±2 K typical, ±5 K worst case
- M: ±1.5% typical
- SG: ±0.005 typical
- Distribution coefficients: ±5% on (P_o, A, B)
