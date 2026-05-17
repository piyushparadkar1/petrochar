"""
Tab 6 — Validation.

In-app reproduction of Riazi MNL50 textbook examples.  Lets reviewers and
end users verify that petrochar reproduces published numerical results without
leaving the application.

Three examples:
  1. Riazi Example 4.13 — Generalized distribution fit (Eq. 4.56) on North Sea
     C7+ M, Tb, SG data.  Compares fitted (P_o, A, B) against Table 4.13.
  2. Riazi Example 4.14 — 3-point Gaussian quadrature discretization.  Uses
     Table 4.22 parameters directly (from_params) and compares M_i, z_i.
  3. Riazi Example 4.7 (Table 4.11) — North Sea gas condensate C7+ data table;
     reproduces the cumulative fraction data as reported.

Reference data loaded from data/riazi_reference/*.csv at runtime.

References
----------
Riazi MNL50 §4.5.4.1 (Eq. 4.56); Table 4.13 (p. 173); Table 4.22 (p. 186);
Table 4.11 (p. 171).
"""

from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd
import streamlit as st

from core.distribution import GeneralizedDistribution
from core.quadrature import discretize_generalized

# ── Data paths ────────────────────────────────────────────────────────────────
_DATA_DIR = pathlib.Path(__file__).parent.parent / "data" / "riazi_reference"
_TABLE_4_11 = _DATA_DIR / "table_4_11_example_4_7.csv"
_TABLE_4_13 = _DATA_DIR / "table_4_13_distribution_coeffs.csv"
_TABLE_4_22 = _DATA_DIR / "table_4_22_quadrature_3pt.csv"

# ── Table 4.22 parameters (used in Example 4.14) ─────────────────────────────
_M_O_4_22  = 90.0
_A_M_4_22  = 0.3324
_B_M_4_22  = 1.096
_M_AV_EXP  = 118.9   # experimental bulk MW (Table 4.11)


def render(ss) -> None:
    """Render Tab 6 — Validation.

    Parameters
    ----------
    ss : streamlit.session_state
    """
    st.header("Validation")
    st.write(
        "In-app reproduction of Riazi MNL50 textbook examples.  "
        "Run any example to verify that petrochar matches published values."
    )

    # ── Architectural commitments statement ───────────────────────────────────
    st.subheader("Architectural commitments")
    st.markdown(
        """
        petrochar operates under the following non-negotiable architectural rules,
        which govern every characterization:

        - **No solvent-deasphalting data required.**  No input from
          solvent-deasphalting operations is used or needed.
          Inputs are routine refinery assay data only (distillation
          curve, bulk SG, bulk MW, SARA wt%, recovery fraction).
        - **No tuning against process observables.**  Pseudo-component properties
          (Tb, SG, PC-SAFT parameters) are determined entirely by published
          correlations.  There are no adjustable parameters fitted to AOP or
          solubility data.
        - **Asphaltenes are always discrete.**  The asphaltene component is
          appended as a single discrete pseudo-component using Gonzalez 2007
          nanoaggregate defaults (m=33, sigma=4.3 A, eps/k=400 K).  It does not
          enter the distillation-curve fitting domain.
        - **Watson K aromaticity mapping (Decision 21).**  The aromaticity factor
          gamma is derived deterministically from per-component Watson K via a
          linear clamp, not fitted to AOP data.  This is a deliberate deviation
          from Panuganti 2012 and must be disclosed in Paper 1 §2.
        - **Self-consistent Tb/SG assignment (Decision 25).**  Pseudo-component
          Tb_i is derived from riazi_daubert_Tb(M_i, SG_i) via root finding, not
          by evaluating the Tb distribution at quadrature nodes.
        - **Recovery-aware quadrature with heavy-resin lump (Decisions 30, 31,
          Phase 11).**  When the distillation curve covers only
          `recovery_fraction` < 1.0 of the total feed mass, the M distribution
          is fitted on a scaled xc basis (xc_scaled = xc_raw / recovery_fraction)
          so it characterizes the distillable subfraction only.  The unmeasured
          tail mass (1 - recovery_fraction - f_asp) is represented by a single
          discrete heavy-resin lump whose MW and SG are fixed by mass-balance
          closure on bulk MW and bulk SG.  Heavy-resin Tb is derived from
          constant Watson K (Tb_hr = (K_W_bulk · SG_hr)^3 / 1.8).  Bulk MW and
          bulk SG become hard closure constraints — the previous "MW is a
          diagnostic" stance from Decision 27 is superseded.

        For complete details see `docs/validation_report.md`.
        """
    )
    st.divider()

    # ── Example selector ──────────────────────────────────────────────────────
    example = st.selectbox(
        "Select Riazi textbook example to reproduce",
        options=[
            "Riazi Example 4.13 — Distribution fit (Table 4.13)",
            "Riazi Example 4.14 — 3-point quadrature (Table 4.22)",
            "Riazi Example 4.7 (Table 4.11) — North Sea C7+ data",
        ],
        index=0,
    )

    if st.button("Run validation", type="primary"):
        if "Example 4.13" in example:
            _run_example_4_13()
        elif "Example 4.14" in example:
            _run_example_4_14()
        else:
            _run_example_4_7()


# ── Example 4.13 ──────────────────────────────────────────────────────────────

def _run_example_4_13() -> None:
    """Reproduce Riazi MNL50 Example 4.13.

    Fits Riazi Eq. 4.56 to North Sea C7+ M, Tb, SG data (Table 4.11).
    Compares fitted parameters against published Table 4.13 values.

    Pass criteria (from test_phase3_riazi_example_4_13.py):
      - P_o, A, B each within 5% of Table 4.13 reference values.
    """
    st.subheader("Example 4.13 — Riazi Eq. 4.56 distribution fit")
    st.caption("Source: Riazi MNL50 §4.5.4.1, Table 4.13 (p. 173).")

    # Load reference data
    try:
        df11 = pd.read_csv(_TABLE_4_11, comment="#")
        df13 = pd.read_csv(_TABLE_4_13, comment="#")
    except FileNotFoundError as exc:
        st.error(f"Reference CSV not found: {exc}.  "
                 "Place Riazi reference CSVs in data/riazi_reference/.")
        return

    # Prepare fit inputs (strip rows with NaN Tb before Tb fit)
    df_valid = df11.dropna(subset=["Tb_K"]).copy()

    xc_mole   = df11["x_cm"].values.astype(float)    # mole basis CDF
    M_data    = df11["M_g_per_mol"].values.astype(float)
    xc_weight = df_valid["x_cw"].values.astype(float)
    Tb_data   = df_valid["Tb_K"].values.astype(float)
    xc_vol    = df11["x_cv"].values.astype(float)
    SG_data   = df11["SG_dimensionless"].values.astype(float)

    # Fit three distributions (3-param)
    m_fit  = GeneralizedDistribution().fit(xc_mole,   M_data,  mode="3param")
    tb_fit = GeneralizedDistribution().fit(xc_weight, Tb_data, mode="3param")
    sg_fit = GeneralizedDistribution().fit(xc_vol,    SG_data, mode="3param")

    # 2-param fits for reference
    m_fit2  = GeneralizedDistribution().fit(xc_mole,   M_data,  mode="2param", B_fixed=1.0)
    tb_fit2 = GeneralizedDistribution().fit(xc_weight, Tb_data, mode="2param", B_fixed=1.5)
    sg_fit2 = GeneralizedDistribution().fit(xc_vol,    SG_data, mode="2param", B_fixed=3.0)

    # Reference rows from Table 4.13
    ref_M_3p  = df13[df13["Property"] == "M"][df13["Param_count"] == 3].iloc[0]
    ref_Tb_3p = df13[df13["Property"] == "Tb_K"][df13["Param_count"] == 3].iloc[0]
    ref_SG_3p = df13[df13["Property"] == "SG"][df13["xc_basis"] == "volume"][df13["Param_count"] == 3].iloc[0]
    ref_M_2p  = df13[df13["Property"] == "M"][df13["Param_count"] == 2].iloc[0]
    ref_Tb_2p = df13[df13["Property"] == "Tb_K"][df13["Param_count"] == 2].iloc[0]
    ref_SG_2p = df13[df13["Property"] == "SG"][df13["xc_basis"] == "volume"][df13["Param_count"] == 2].iloc[0]

    # Build comparison table
    rows = []
    for label, fit, ref in [
        ("M 3-param",   m_fit,  ref_M_3p),
        ("M 2-param",   m_fit2, ref_M_2p),
        ("Tb 3-param",  tb_fit, ref_Tb_3p),
        ("Tb 2-param",  tb_fit2, ref_Tb_2p),
        ("SG 3-param",  sg_fit, ref_SG_3p),
        ("SG 2-param",  sg_fit2, ref_SG_2p),
    ]:
        fq = fit.fit_quality
        rows.append({
            "Case":           label,
            "P_o calc":       round(fit.P_o, 4),
            "P_o ref":        round(float(ref["P_o"]), 4),
            "P_o dev %":      f"{abs(fit.P_o - ref['P_o']) / ref['P_o'] * 100:.2f}",
            "A calc":         round(fit.A, 4),
            "A ref":          round(float(ref["A"]), 4),
            "A dev %":        f"{abs(fit.A - ref['A']) / ref['A'] * 100:.2f}",
            "B calc":         round(fit.B, 4),
            "B ref":          round(float(ref["B"]), 4),
            "B dev %":        f"{abs(fit.B - ref['B']) / ref['B'] * 100:.2f}",
            "RMS calc":       round(fq["RMS"], 4),
            "RMS ref":        round(float(ref["RMS"]), 3),
            "Pass (5%)?":     _all_pass_5pct(fit, ref),
        })

    result_df = pd.DataFrame(rows)
    st.dataframe(result_df, use_container_width=True, hide_index=True)

    n_pass = int(result_df["Pass (5%)?"].str.contains("PASS").sum())
    n_fail = len(result_df) - n_pass
    if n_fail == 0:
        st.success(f"All {len(result_df)} cases: P_o, A, B each within 5% of Table 4.13. PASS")
    else:
        st.warning(
            f"{n_pass}/{len(result_df)} cases pass the 5% tolerance gate. "
            f"{n_fail} case(s) outside tolerance."
        )


def _all_pass_5pct(fit, ref) -> str:
    """Return 'PASS' if P_o, A, B each within 5% of reference."""
    tol = 0.05
    po_ok = abs(fit.P_o - ref["P_o"]) / ref["P_o"] <= tol
    a_ok  = abs(fit.A   - ref["A"])   / ref["A"]   <= tol
    b_ok  = abs(fit.B   - ref["B"])   / ref["B"]   <= tol
    if po_ok and a_ok and b_ok:
        return "PASS"
    parts = []
    if not po_ok:
        parts.append("P_o")
    if not a_ok:
        parts.append("A")
    if not b_ok:
        parts.append("B")
    return f"FAIL ({', '.join(parts)})"


# ── Example 4.14 ──────────────────────────────────────────────────────────────

def _run_example_4_14() -> None:
    """Reproduce Riazi MNL50 Example 4.14.

    Uses Table 4.22 parameters (M_o=90, A_M=0.3324, B_M=1.096) to run
    3-point Gauss-Laguerre discretization and compares M_i, z_i against
    published Table 4.22 values.

    Pass criteria (from test_phase5_riazi_example_4_14.py):
      - M_i within 1% of Table 4.22 values (103.6, 154.6, 252.2 g/mol).
      - z_i within 0.005 of Table 4.22 values (0.711, 0.279, 0.010).
      - 3-pt M_av within 0.5% of experimental 118.9 g/mol.
    """
    st.subheader("Example 4.14 — 3-point Gauss-Laguerre quadrature")
    st.caption(
        "Source: Riazi MNL50 §4.6.1.1, Table 4.22 (p. 186).  "
        "Parameters: M_o=90, A_M=0.3324, B_M=1.096."
    )

    try:
        df_ref = pd.read_csv(_TABLE_4_22, comment="#")
    except FileNotFoundError as exc:
        st.error(f"Reference CSV not found: {exc}")
        return

    # Instantiate distribution from published Table 4.22 parameters
    dist = GeneralizedDistribution.from_params(_M_O_4_22, _A_M_4_22, _B_M_4_22)
    comps = discretize_generalized(3, dist)

    M_ref = [103.6, 154.6, 252.2]
    z_ref = [0.711, 0.279, 0.010]

    M_av_calc = float(sum(c.z * c.M for c in comps))
    M_av_dev  = abs(M_av_calc - _M_AV_EXP) / _M_AV_EXP * 100.0

    rows = []
    for i, (c, M_r, z_r) in enumerate(zip(comps, M_ref, z_ref)):
        M_dev = abs(c.M - M_r) / M_r * 100.0
        z_dev = abs(c.z - z_r)
        rows.append({
            "Node":       i + 1,
            "M calc (g/mol)":  round(c.M, 2),
            "M ref (g/mol)":   M_r,
            "M dev %":         f"{M_dev:.3f}",
            "M pass (1%)?":    "PASS" if M_dev <= 1.0 else "FAIL",
            "z calc":          round(c.z, 5),
            "z ref":           z_r,
            "|z_calc - z_ref|": f"{z_dev:.5f}",
            "z pass (0.005)?":  "PASS" if z_dev <= 0.005 else "FAIL",
        })

    result_df = pd.DataFrame(rows)
    st.dataframe(result_df, use_container_width=True, hide_index=True)

    # M_av closure
    m_av_col, _ = st.columns([1, 2])
    with m_av_col:
        m_av_pass = M_av_dev <= 0.5
        st.metric(
            label="3-pt M_av vs experimental 118.9 g/mol",
            value=f"{M_av_calc:.2f} g/mol",
            delta=f"{'PASS' if m_av_pass else 'FAIL'} ({M_av_dev:.3f}%)",
            delta_color="normal" if m_av_pass else "inverse",
        )

    # Overall pass/fail
    n_pass_M = sum(1 for r in rows if r["M pass (1%)?"] == "PASS")
    n_pass_z = sum(1 for r in rows if r["z pass (0.005)?"] == "PASS")
    all_pass = (n_pass_M == 3) and (n_pass_z == 3) and m_av_pass

    if all_pass:
        st.success(
            "All 3 nodes: M_i within 1%, z_i within 0.005.  "
            f"M_av = {M_av_calc:.2f} vs 118.9 g/mol ({M_av_dev:.3f}% < 0.5%).  "
            "PASS"
        )
    else:
        st.error(
            f"M pass: {n_pass_M}/3 | z pass: {n_pass_z}/3 | "
            f"M_av dev: {M_av_dev:.3f}%.  One or more criteria FAILED."
        )


# ── Example 4.7 (Table 4.11) ─────────────────────────────────────────────────

def _run_example_4_7() -> None:
    """Display Riazi Example 4.7 — North Sea gas condensate C7+ data.

    Loads Table 4.11 data, shows the raw table as published, and reproduces
    computed Tb/M/SG distributions fitted to it (cross-reference with
    Example 4.13 results above).
    """
    st.subheader("Example 4.7 — North Sea gas condensate C7+ data (Table 4.11)")
    st.caption(
        "Source: Riazi MNL50 Table 4.11 (p. 171).  "
        "12-fraction C7+ North Sea gas condensate."
    )

    try:
        df = pd.read_csv(_TABLE_4_11, comment="#")
    except FileNotFoundError as exc:
        st.error(f"Reference CSV not found: {exc}")
        return

    # Display raw table with nicer column names
    display_df = df.rename(columns={
        "Fraction_No": "Frac.",
        "Carbon_No_label": "SCN",
        "x_w": "x_w",
        "M_g_per_mol": "M (g/mol)",
        "Tb_K": "Tb (K)",
        "SG_dimensionless": "SG",
        "x_m": "x_m",
        "x_v": "x_v",
        "x_cm": "x_cm",
        "x_cw": "x_cw",
        "x_cv": "x_cv",
    })
    st.dataframe(display_df, use_container_width=True, hide_index=True)

    st.caption(
        "x_w: weight fraction; x_m: mole fraction; x_v: volume fraction.  "
        "x_cm/x_cw/x_cv: cumulative mole/weight/volume fractions."
    )

    # Reproduce bulk properties from the table
    x_w  = df["x_w"].values.astype(float)
    M    = df["M_g_per_mol"].values.astype(float)
    SG   = df["SG_dimensionless"].values.astype(float)

    # Bulk MW: number-average (1 / sum(x_w/M))
    M_av  = 1.0 / float(np.sum(x_w / M))
    # Bulk SG: volume-additive (1 / sum(x_w/SG))
    SG_av = 1.0 / float(np.sum(x_w / SG))

    st.subheader("Bulk closure check")
    col1, col2 = st.columns(2)
    with col1:
        # Riazi Table 4.11 gives M_C7+ = 119 g/mol
        ref_M = 119.0
        M_dev = abs(M_av - ref_M) / ref_M * 100.0
        st.metric(
            "M_av (number-average)",
            f"{M_av:.2f} g/mol",
            delta=f"ref 119 g/mol | dev {M_dev:.2f}%",
            delta_color="normal" if M_dev < 1.0 else "off",
        )
    with col2:
        ref_SG = 0.7597
        SG_dev = abs(SG_av - ref_SG) / ref_SG * 100.0
        st.metric(
            "SG_av (volume-additive)",
            f"{SG_av:.4f}",
            delta=f"ref 0.7597 | dev {SG_dev:.2f}%",
            delta_color="normal" if SG_dev < 1.0 else "off",
        )

    if M_dev < 1.0 and SG_dev < 1.0:
        st.success(
            "Bulk M_av and SG_av reproduce Table 4.11 within 1%.  PASS"
        )
    else:
        st.warning(
            "One or more bulk properties deviate > 1% from Table 4.11 reference."
        )
