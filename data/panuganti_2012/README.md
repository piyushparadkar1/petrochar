# Panuganti 2012 Reference Materials

This folder contains reference materials transcribed from:

> Panuganti, S. R., Vargas, F. M., Gonzalez, D. L., Kurup, A. S., Chapman, W. G. (2012). PC-SAFT characterization of crude oils and modeling of asphaltene phase behavior. *Fuel*, 93, 658-669. DOI: 10.1016/j.fuel.2011.09.028.

## Files

**`Panuganti_2012_Fuel_full.pdf`** — full paper, vendored locally for offline reference during Phase 7 implementation. **Do NOT commit to git.** Add to `.gitignore` for the same fair-use reasons that apply to the Riazi MNL50 PDF.

**`table_5_light_component_pcsaft_params.csv`** — PC-SAFT parameters for N₂, CO₂, H₂S, C₁, C₂, C₃. Reference values from Gross & Sadowski 2001 reproduced in Panuganti Table 5.

**`table_6_correlations.csv`** — Saturates and Aromatics+Resins correlation formulas. Saturates: single correlation set; Aro+Resins: γ-interpolated between benzene-derivative endpoint (γ=0) and PNA endpoint (γ=1).

**`tables_10_11_12_crude_pcsaft.csv`** — PC-SAFT characterized parameters for Crudes A, B, C from Panuganti's published results. Used as Phase 7 numerical pass-gate.

## Phase 7 pass-gate methodology

For each of Crudes A, B, C:
1. Read (MW, γ) for Saturates and Aromatics+Resins from `tables_10_11_12_crude_pcsaft.csv`.
2. Plug into the petrochar implementation of Panuganti Table 6 correlations.
3. Compare computed (m, σ, ε/k) against the published values.
4. Tolerance: ±1% on m, ±0.5% on σ, ±2% on ε/k. (The paper reports values to 3-4 significant figures; tighter tolerance would test publication rounding, not implementation correctness.)

## Note on γ source

Panuganti 2012 defines γ as a free parameter adjusted to match crude density and bubble pressure. petrochar does NOT do this — it derives γ deterministically per pseudo-component from Watson K (Watson 1933) via the linear clamp `γ = clamp((13.0 − K_W) / (13.0 − 9.5), 0, 1)`. This is a deliberate methodological substitution to preserve petrochar's feed-side calibration commitment.

For Phase 7's pass-gate, γ values are taken directly from Panuganti's published Tables 10/11/12 to isolate testing of the correlation formulas from testing of γ derivation. Tests of γ-from-Watson-K were already validated in Phase 1.

## Fair use

These materials are reproduced strictly for software validation purposes. The full PDF is vendored for offline reference during development. petrochar's users must cite Panuganti 2012 directly in any work derived from this tool.
