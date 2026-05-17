"""
Tab 1 — Input Data.

Accepts all user inputs and triggers the full characterization pipeline on
"Run pipeline".  Stores all intermediates in session_state for Tabs 2-6.

Session state written:
    session_state['inputs']          — validated user inputs dict
    session_state['pipeline_result'] — complete pipeline output dict
    session_state['pipeline_ready']  — True after successful run

Pipeline (mirrors test_phase8_pipeline.py fixture, 9 steps):
    1.  DistillationCurve              (core/distillation.py)
    2.  Tb distribution fit (3-param)  (core/distribution.py) — DIAGNOSTIC ONLY
    3.  SG per cut (constant Watson K) (core/sg_distribution.py)
    4.  M per cut + M distribution fit (core/mw_distribution.py + distribution.py)
    5.  5-pt GL discretization + self-consistent (Tb, SG) solve
                                       (core/quadrature.py + correlations.py)
    6.  Append discrete asphaltene     (core/sara.py)
    7.  K_W bin closure check          (core/sara.py)
    8.  Per-component Watson K + gamma (core/watson_k.py)
    9.  PC-SAFT parameter table        (core/pcsaft_params.py)

References
----------
Decision 25: Tb_i derived from riazi_daubert_Tb(M_i, SG_i) via constant K_W.
Decision 26: is_asphaltene flag (not Tb > 1000 K) identifies the ASP component.
Decision 27: M_av gate is GL vs analytic mean; M_DIST_TARGET is diagnostic.
Architecture: CLAUDE.md §ARCHITECTURE COMMITMENTS.
"""

from __future__ import annotations

import io
import time
import warnings

import numpy as np
import pandas as pd
import streamlit as st
from scipy.optimize import brentq

from core.correlations import riazi_daubert_Tb
from core.distillation import DistillationCurve
from core.distribution import GeneralizedDistribution
from core.mw_distribution import compute_M_array
from core.pcsaft_params import generate_pcsaft_table, propane_params
from core.quadrature import Pseudocomponent, discretize_generalized
from core.sara import (
    append_asphaltene,
    append_heavy_resin_and_asphaltene,
    kw_bin_check,
    validate_sara,
)
from core.sg_distribution import sg_from_watson_k
from core.watson_k import compute_K_W_per_pseudocomponent

# ── Default 11-point VTB-like template (D1160 AET, weight basis) ─────────────
# Used as the pre-populated data_editor starting point.
_TEMPLATE_PCT  = [5,  10, 20,  30,  40,  50,  60,  70,  80,  90,  95]
_TEMPLATE_T_C  = [305, 325, 360, 390, 415, 440, 465, 488, 510, 530, 540]

# ── Self-consistent (M, K_W) → (Tb, SG) root-finding ────────────────────────
# Decision 25 (Phase 8 rework): Tb_i is derived from riazi_daubert_Tb(M_i, SG_i),
# not from evaluation of the Tb distribution at a quadrature node.

def _solve_tb_from_M_constant_kw(M_i: float, K_W_bulk: float) -> tuple[float, float]:
    """Solve the coupled (Tb_i, SG_i) pair for a pseudo-component.

    Under constant Watson K:
        SG_i = (1.8 * Tb_i)^(1/3) / K_W_bulk
        Tb_i = riazi_daubert_Tb(M_i, SG_i)

    Solved by 1-D root finding on the ascending branch of Eq. 2.57 only.
    For M_i exceeding Eq. 2.57's M_peak at SG_at_Tb_hi (the constant-K_W SG
    at Tb = 990 K), the residual has the same sign at both bracket ends and
    brentq could lock onto an unphysical descending-branch root.  Such above-
    peak components are capped at Tb_hi = 990 K with corresponding SG (Phase
    11 above-peak fallback, consistent with riazi_daubert_M's pattern).

    Parameters
    ----------
    M_i : float        Pseudo-component molecular weight, g/mol.
    K_W_bulk : float   Bulk Watson K factor.

    Returns
    -------
    Tb_i, SG_i : float, float
    """
    Tb_lo, Tb_hi = 300.0, 990.0
    SG_at_Tb_hi = (1.8 * Tb_hi) ** (1.0 / 3.0) / K_W_bulk
    denom_at_hi = 7.5152e-4 * SG_at_Tb_hi - 1.6514e-4
    if denom_at_hi > 0:
        M_peak_at_Tb_hi = 0.5369 / denom_at_hi
    else:
        M_peak_at_Tb_hi = float("inf")

    # Above-peak fallback: M_i exceeds Eq. 2.57's ascending-branch upper bound.
    if M_i >= M_peak_at_Tb_hi * 0.9999:
        warnings.warn(
            f"_solve_tb_from_M_constant_kw: M_i={M_i:.1f} g/mol exceeds "
            f"Eq. 2.57 M_peak ({M_peak_at_Tb_hi:.1f} g/mol) at SG_at_Tb_hi="
            f"{SG_at_Tb_hi:.4f}.  Capping Tb at {Tb_hi:.0f} K (above-peak "
            f"fallback).",
            UserWarning, stacklevel=2,
        )
        return float(Tb_hi), float(SG_at_Tb_hi)

    def _residual(Tb_try: float) -> float:
        SG_try = (1.8 * Tb_try) ** (1.0 / 3.0) / K_W_bulk
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            return riazi_daubert_Tb(M_i, SG_try) - Tb_try

    try:
        Tb_i = float(brentq(_residual, Tb_lo, Tb_hi, xtol=0.01))
    except ValueError:
        warnings.warn(
            f"_solve_tb_from_M_constant_kw: brentq bracket failed for "
            f"M_i={M_i:.1f} g/mol; capping Tb at {Tb_hi:.0f} K.",
            UserWarning, stacklevel=2,
        )
        Tb_i = Tb_hi
    SG_i = float((1.8 * Tb_i) ** (1.0 / 3.0) / K_W_bulk)
    return Tb_i, SG_i


# ── Full pipeline execution ───────────────────────────────────────────────────

def _run_pipeline(
    pct: np.ndarray,
    temps_K: np.ndarray,
    method: str,
    basis: str,
    sg_bulk: float,
    mw_bulk: float,
    sara: dict,
    n_quad: int,
    recovery_fraction: float,
) -> dict:
    """Execute the full 9-step characterization pipeline.

    When recovery_fraction < 1.0, the M distribution is fitted on a SCALED xc
    basis (xc_scaled = xc_raw / recovery_fraction) so it characterizes the
    distillable subfraction only.  The unmeasured tail mass
    (1 - recovery_fraction - f_asp) becomes a single discrete heavy-resin lump
    whose MW, SG, and Tb are fixed by mass-balance closure on bulk MW, bulk SG,
    and constant Watson K (Decisions 30, 31, Phase 11).

    Parameters
    ----------
    pct : ndarray       Cumulative percent distilled (incl. any IBP 0% point).
    temps_K : ndarray   Temperatures at each cut point, K.
    method : str        Distillation method tag.
    basis : str         Distillation basis tag.
    sg_bulk : float     Bulk specific gravity at 15 degC.
    mw_bulk : float     Bulk molecular weight, g/mol.
    sara : dict         Keys SAT, ARO, RES, ASP; values in wt%.
    n_quad : int        Number of Gauss-Laguerre quadrature points (3 or 5).
    recovery_fraction : float
        Mass fraction of total feed covered by the distillation curve.  Must
        be in (0, 1].  recovery_fraction = 1.0 dispatches to the legacy Phase
        8 path (no heavy-resin lump).

    Returns
    -------
    dict  — all intermediate and final objects needed by Tabs 2-6, including
    `has_heavy_resin`, `M_hr`, `SG_hr`, `Tb_hr`, `f_hr` keys when partial
    recovery is detected.
    """
    t0 = time.perf_counter()

    # ── Step 1: DistillationCurve ─────────────────────────────────────────────
    dc  = DistillationCurve(pct, temps_K, method=method, basis=basis)
    tbp = dc.to_tbp()

    # Interior points for distribution fitting (exclude 0% IBP if present).
    mask        = tbp.pct > 0.0
    xc_int_raw  = tbp.pct[mask] / 100.0
    tb_int_raw  = tbp.temps_K[mask]

    # Recovery-aware xc scaling (Phase 11).  When recovery_fraction < 1.0, the
    # measured xc points span [0, recovery_fraction] of the total feed; scale
    # them to [0, 1] of the distillable subfraction so GL quadrature samples
    # within the distillable range.  The boundary point at xc_scaled ≈ 1 is
    # excluded to avoid singular fit behaviour.
    if recovery_fraction < 1.0 - 1e-9:
        xc_scaled  = xc_int_raw / recovery_fraction
        keep       = xc_scaled < 1.0 - 1e-9
        xc_int     = xc_scaled[keep]
        tb_int     = tb_int_raw[keep]
    else:
        xc_int     = xc_int_raw
        tb_int     = tb_int_raw

    # ── Step 2: Tb distribution (3-param, DIAGNOSTIC ONLY) ───────────────────
    tb_dist = GeneralizedDistribution().fit(xc_int, tb_int, mode="3param")
    Tb_pred = np.array([tb_dist.P(x) for x in xc_int])
    rms_tb  = float(np.sqrt(np.mean((Tb_pred - tb_int) ** 2)))
    tb_dist_2param = GeneralizedDistribution().fit(xc_int, tb_int,
                                                    mode="2param", B_fixed=1.5)
    Tb_pred_2p = np.array([tb_dist_2param.P(x) for x in xc_int])
    rms_tb_2p  = float(np.sqrt(np.mean((Tb_pred_2p - tb_int) ** 2)))

    # ── Step 3: SG per cut (constant Watson K) ────────────────────────────────
    SG_cuts, K_W_bulk = sg_from_watson_k(tb_int, mw_bulk, sg_bulk)

    # ── Step 4: M per cut + M distribution fit ────────────────────────────────
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        M_cuts = compute_M_array(tb_int, SG_cuts)
    m_dist = GeneralizedDistribution().fit(xc_int, M_cuts, mode="3param")

    # ── Step 5: Quadrature discretization + self-consistent (Tb, SG) ─────────
    comps_raw  = discretize_generalized(n_quad, m_dist)
    comps_dist = []
    for c in comps_raw:
        Tb_i, SG_i = _solve_tb_from_M_constant_kw(c.M, K_W_bulk)
        comps_dist.append(Pseudocomponent(
            z=c.z, M=c.M, Tb_K=Tb_i, SG=SG_i,
            xc_lower=c.xc_lower, xc_upper=c.xc_upper,
        ))

    # ── Step 6: Recovery-aware assembly (HR lump + asphaltene) ────────────────
    validate_sara(sara["SAT"], sara["ARO"], sara["RES"], sara["ASP"])
    asm = append_heavy_resin_and_asphaltene(
        comps_dist,
        recovery_fraction=recovery_fraction,
        asp_wt_pct=sara["ASP"],
        M_bulk=mw_bulk,
        SG_bulk=sg_bulk,
        K_W_bulk=K_W_bulk,
    )
    comps6 = asm["components"]

    # ── Step 7: K_W bin closure check ─────────────────────────────────────────
    kw_result = kw_bin_check(comps6, sara, flag_tol=5.0)

    # ── Step 8: Per-component Watson K and gamma ──────────────────────────────
    comps8 = compute_K_W_per_pseudocomponent(comps6)

    # ── Step 9: PC-SAFT parameter table ───────────────────────────────────────
    df = generate_pcsaft_table(comps8)
    m_c3, s_c3, e_c3 = propane_params()
    prop_row = pd.DataFrame([{
        "component_type": "propane",
        "z": float("nan"), "M": 44.096, "Tb_K": 231.11, "SG": 0.507,
        "K_W": float("nan"), "gamma": float("nan"),
        "m": m_c3, "sigma_A": s_c3, "eps_over_k_K": e_c3,
    }])
    df_full = pd.concat([df, prop_row], ignore_index=True)

    # ── Bulk closure metrics ───────────────────────────────────────────────────
    # Distillable nodes only (exclude HR and ASP) for the GL-vs-analytic gate.
    gl_dist_only = [c for c in comps8
                    if not c.is_asphaltene and not c.is_heavy_resin]
    z_d         = np.array([c.z for c in gl_dist_only])
    M_d         = np.array([c.M for c in gl_dist_only])
    GL_M_av     = float(np.dot(z_d, M_d) / z_d.sum())
    dist_M_mean = m_dist.average()

    all_comps = comps8
    M_mix     = float(sum(c.z * c.M for c in all_comps))
    xw_all    = np.array([c.z * c.M / M_mix for c in all_comps])
    sg_all    = np.array([c.SG for c in all_comps])
    SG_av_full = float(1.0 / np.dot(xw_all, 1.0 / sg_all))

    M_dist_target = (mw_bulk - sara["ASP"] / 100.0 * 1700.0) / (1.0 - sara["ASP"] / 100.0)

    elapsed = time.perf_counter() - t0

    return dict(
        # Raw inputs preserved for display
        pct=pct, temps_K=temps_K, method=method, basis=basis,
        sg_bulk=sg_bulk, mw_bulk=mw_bulk, sara=sara, n_quad=n_quad,
        recovery_fraction=recovery_fraction,
        # Step outputs
        dc=dc, tbp=tbp,
        xc_int=xc_int, tb_int=tb_int,
        xc_int_raw=xc_int_raw, tb_int_raw=tb_int_raw,
        tb_dist=tb_dist, rms_tb=rms_tb,
        tb_dist_2param=tb_dist_2param, rms_tb_2p=rms_tb_2p,
        SG_cuts=SG_cuts, K_W_bulk=K_W_bulk,
        M_cuts=M_cuts, m_dist=m_dist,
        comps_dist=comps_dist, comps6=comps6,
        kw_result=kw_result, comps8=comps8,
        df=df, df_full=df_full,
        GL_M_av=GL_M_av, dist_M_mean=dist_M_mean,
        SG_av_full=SG_av_full,
        M_dist_target=M_dist_target,
        # Phase 11 heavy-resin diagnostics
        has_heavy_resin=asm["has_heavy_resin"],
        f_hr=asm["f_hr"],
        M_hr=asm["M_hr"],
        SG_hr=asm["SG_hr"],
        Tb_hr=asm["Tb_hr"],
        M_dist_av=asm["M_dist_av"],
        elapsed=elapsed,
    )


# ── Tab render ────────────────────────────────────────────────────────────────

def render(ss) -> None:
    """Render Tab 1 — Input Data.

    Parameters
    ----------
    ss : streamlit.session_state
        Session state shared across all tabs.
    """
    st.header("Input Data")
    st.write(
        "Enter the distillation curve and bulk properties for your heavy "
        "petroleum fraction, then click **Run pipeline** to characterize it."
    )

    # ── Distillation curve input ───────────────────────────────────────────────
    st.subheader("Distillation curve")

    col_method, col_basis, col_unit = st.columns(3)
    with col_method:
        method = st.selectbox(
            "Distillation method",
            options=["D1160_AET", "TBP", "D86", "D7169"],
            index=0,
            help=(
                "D1160_AET and TBP are pass-throughs (temperatures unchanged). "
                "D86 uses Daubert Eqs. 3.20-3.22 conversion. "
                "D7169 (simulated distillation) is architecturally supported "
                "but the conversion path is not yet implemented."
            ),
        )
    with col_basis:
        basis = st.selectbox(
            "Basis",
            options=["weight", "volume", "mole"],
            index=0,
            help=(
                "Only 'weight' basis is fully supported. The DistillationCurve "
                "class carries the basis tag for volume and mole, but the "
                "density-distribution conversion is not wired through the "
                "current pipeline."
            ),
        )
    with col_unit:
        temp_unit = st.selectbox(
            "Temperature unit",
            options=["degC", "K", "degF"],
            index=0,
        )

    # D7169 warning — Option B: keep visible, warn, disable run button.
    # Signals architectural intent without crashing the pipeline.
    _d7169_selected = (method == "D7169")
    if _d7169_selected:
        st.warning(
            "**D7169 (simulated distillation) input is not yet implemented in "
            "the pipeline.** The `DistillationCurve` class supports D7169 "
            "architecturally, but the Phase 2 D7169 → TBP conversion path is "
            "not wired through. Use **D1160_AET**, **D86**, or **TBP** for now."
        )

    if basis != "weight":
        st.warning(
            "Non-weight basis: the DistillationCurve class carries the basis tag, "
            "but the volume/mole → weight conversion requires a density distribution "
            "not exposed in the current pipeline. Only **weight** basis is fully "
            "operational."
        )

    input_mode = st.radio(
        "Data entry",
        options=["Manual entry (template)", "Upload CSV / Excel"],
        horizontal=True,
        label_visibility="collapsed",
    )

    if input_mode == "Upload CSV / Excel":
        uploaded = st.file_uploader(
            "CSV or Excel (.xlsx) with columns: cumulative_fraction (0–100) and temperature",
            type=["csv", "xlsx"],
            help=(
                "Header row required.  Column names must be exactly "
                "'cumulative_fraction' (percent distilled) and 'temperature' "
                "(in the unit selected above).  For Excel, the first sheet is read."
            ),
        )
        if uploaded is not None:
            try:
                fname = uploaded.name.lower()
                if fname.endswith(".xlsx"):
                    df_upload = pd.read_excel(uploaded, sheet_name=0)
                else:
                    df_upload = pd.read_csv(uploaded)

                if "cumulative_fraction" not in df_upload.columns or "temperature" not in df_upload.columns:
                    st.error(
                        "File must have columns named **'cumulative_fraction'** and "
                        "**'temperature'** (exact match, case-sensitive).  "
                        f"Columns found: {list(df_upload.columns)}.  "
                        "Rename your columns and re-upload."
                    )
                    pct_raw  = np.array(_TEMPLATE_PCT, dtype=float)
                    temp_raw = np.array(_TEMPLATE_T_C, dtype=float)
                    has_data = False
                else:
                    pct_raw  = df_upload["cumulative_fraction"].values.astype(float)
                    temp_raw = df_upload["temperature"].values.astype(float)
                    has_data = True
                    st.success(f"Loaded {len(pct_raw)} points from {uploaded.name}.")
            except Exception as exc:
                st.error(f"File read error: {exc}")
                pct_raw  = np.array(_TEMPLATE_PCT, dtype=float)
                temp_raw = np.array(_TEMPLATE_T_C, dtype=float)
                has_data = False
        else:
            pct_raw  = np.array(_TEMPLATE_PCT, dtype=float)
            temp_raw = np.array(_TEMPLATE_T_C, dtype=float)
            has_data = False
    else:
        # Manual data_editor with 11-point template
        df_template = pd.DataFrame({
            "cumulative_fraction": _TEMPLATE_PCT,
            "temperature":         _TEMPLATE_T_C,
        })
        edited_df = st.data_editor(
            df_template,
            num_rows="dynamic",
            column_config={
                "cumulative_fraction": st.column_config.NumberColumn(
                    "Cumulative fraction (%)", min_value=0.0, max_value=100.0
                ),
                "temperature": st.column_config.NumberColumn(
                    f"Temperature ({temp_unit})"
                ),
            },
            key="dist_editor",
            hide_index=True,
        )
        pct_raw  = edited_df["cumulative_fraction"].values.astype(float)
        temp_raw = edited_df["temperature"].values.astype(float)
        has_data = True

    # Convert temperature to K
    if temp_unit == "degC":
        temps_K_input = temp_raw + 273.15
    elif temp_unit == "degF":
        temps_K_input = (temp_raw - 32.0) * 5.0 / 9.0 + 273.15
    else:
        temps_K_input = temp_raw.copy()

    # Live input validation
    _warnings = []
    if len(pct_raw) < 3:
        _warnings.append("At least 3 data points are required for distribution fitting.")
    if np.any(np.diff(temps_K_input) < 0):
        _warnings.append("Temperatures are not monotonically increasing — check input.")
    if np.any(pct_raw < 0) or np.any(pct_raw > 100):
        _warnings.append("Cumulative fraction values must be in [0, 100]%.")

    # Maximum measured cumulative fraction (used as default for recovery_fraction).
    if has_data and len(pct_raw) > 0:
        _max_pct = float(np.max(pct_raw))
    else:
        _max_pct = 100.0
    _max_recovery_default = max(0.05, min(_max_pct / 100.0, 1.0))

    # ── Bulk properties ────────────────────────────────────────────────────────
    st.subheader("Bulk properties")

    col_sg, col_mw = st.columns(2)
    with col_sg:
        sg_bulk = st.number_input(
            "Bulk SG at 15 degC",
            min_value=0.6,
            max_value=1.5,
            value=1.020,
            step=0.001,
            format="%.3f",
            help="Specific gravity of the full feed including asphaltenes.",
        )
    with col_mw:
        mw_bulk = st.number_input(
            "Bulk MW (g/mol)",
            min_value=50.0,
            max_value=2000.0,
            value=700.0,
            step=1.0,
            format="%.1f",
            help="Molecular weight of the full feed including asphaltenes.",
        )

    # Live validity warnings
    if sg_bulk < 0.7 or sg_bulk > 1.3:
        _warnings.append(
            f"Bulk SG = {sg_bulk:.3f} is outside typical VTB range [0.7, 1.3]."
        )
    if mw_bulk < 100.0 or mw_bulk > 1500.0:
        _warnings.append(
            f"Bulk MW = {mw_bulk:.0f} g/mol is outside typical range [100, 1500]."
        )

    # ── SARA wt% ──────────────────────────────────────────────────────────────
    st.subheader("SARA wt%")

    col_sat, col_aro, col_res, col_asp = st.columns(4)
    with col_sat:
        sara_sat = st.number_input("Saturates", min_value=0.0, max_value=100.0,
                                   value=12.0, step=0.1, format="%.1f")
    with col_aro:
        sara_aro = st.number_input("Aromatics", min_value=0.0, max_value=100.0,
                                   value=38.0, step=0.1, format="%.1f")
    with col_res:
        sara_res = st.number_input("Resins", min_value=0.0, max_value=100.0,
                                   value=38.0, step=0.1, format="%.1f")
    with col_asp:
        sara_asp = st.number_input("Asphaltenes", min_value=0.0, max_value=100.0,
                                   value=12.0, step=0.1, format="%.1f")

    sara_sum = sara_sat + sara_aro + sara_res + sara_asp
    sara_diff = abs(sara_sum - 100.0)
    if sara_diff > 0.5:
        _warnings.append(
            f"SARA components sum to {sara_sum:.2f} wt% (expected 100.0 +/- 0.5 wt%)."
        )
    elif sara_diff > 0.01:
        st.info(f"SARA sum = {sara_sum:.2f} wt% (within 0.5 wt% tolerance).")

    # ── Options ───────────────────────────────────────────────────────────────
    st.subheader("Options")

    col_sgm, col_nq = st.columns(2)
    with col_sgm:
        sg_method = st.selectbox(
            "SG distribution method",
            options=["Constant Watson K", "Generalized SG distribution"],
            index=0,
            help=(
                "Constant Watson K assigns the same K_W to all pseudo-components "
                "— recommended for VTB fractions (Whitson 1983). Generalized SG "
                "distribution fits Riazi Eq. 4.56 to measured SG data."
            ),
        )
        if sg_method != "Constant Watson K":
            st.info(
                "Generalized SG distribution requires measured (xc, SG) data. "
                "Only constant Watson K is wired in this version."
            )
    with col_nq:
        n_quad = st.selectbox(
            "Quadrature points",
            options=[5, 3],
            index=0,
            help="5-point default. 3-point selectable for faster runs.",
        )

    # Recovery-aware quadrature (Phase 11, Option D).
    # Default = maximum measured cumulative fraction.  Set 1.0 only if the
    # distillation curve reaches the end of the distillable subfraction
    # (i.e. asphaltenes plus heavy-resin tail are negligible).
    st.markdown("**Distillation recovery (Phase 11)**")
    recovery_fraction = st.number_input(
        "Recovery fraction (mass fraction of feed covered by the distillation curve)",
        min_value=0.05,
        max_value=1.00,
        value=float(_max_recovery_default),
        step=0.001,
        format="%.3f",
        help=(
            "Mass fraction of the total feed mass covered by the distillation "
            "curve.  When < 1.0, an unmeasured tail (1 - recovery - ASP wt%) "
            "is represented by a discrete heavy-resin lump whose MW and SG are "
            "fixed by mass-balance closure on bulk MW and bulk SG (Phase 11 "
            "Option D, Decisions 30, 31).  Defaults to the maximum cumulative "
            "fraction in your data (here: "
            f"{_max_recovery_default*100:.1f}%).  Set 1.0 only if the curve "
            "endpoint already reaches the start of the asphaltene fraction."
        ),
    )

    # Live validation for recovery_fraction
    f_asp_input = (sara_asp) / 100.0
    if recovery_fraction + f_asp_input > 1.0 + 1e-6:
        _warnings.append(
            f"recovery_fraction ({recovery_fraction:.3f}) + ASP wt% "
            f"({sara_asp:.1f}%) exceeds 100%.  Reduce one to avoid overlap "
            f"between distillation data and asphaltenes."
        )
    elif recovery_fraction < _max_recovery_default - 1e-6:
        _warnings.append(
            f"recovery_fraction ({recovery_fraction:.3f}) is less than the "
            f"maximum measured cumulative fraction "
            f"({_max_recovery_default:.3f}).  Some distillation data will be "
            f"discarded — increase recovery_fraction to keep all measured "
            f"cuts in the fit."
        )
    elif recovery_fraction < 1.0 - 1e-9 and (1.0 - recovery_fraction - f_asp_input) < 0.01:
        st.info(
            f"Heavy-resin mass fraction is small "
            f"(f_hr ≈ {(1.0 - recovery_fraction - f_asp_input)*100:.2f} wt%); "
            f"the HR lump will contribute little to the mixture."
        )

    # ── Warnings display ─────────────────────────────────────────────────────
    for w in _warnings:
        st.warning(w)

    # ── Run pipeline button ───────────────────────────────────────────────────
    st.divider()

    run_disabled = bool(_warnings) or _d7169_selected
    if st.button(
        "Run pipeline",
        type="primary",
        disabled=run_disabled,
        help=(
            "D7169 conversion is not yet implemented — select a different method."
            if _d7169_selected
            else ("Blocked until all input warnings are resolved." if _warnings else "")
        ),
    ):
        sara_dict = {
            "SAT": float(sara_sat),
            "ARO": float(sara_aro),
            "RES": float(sara_res),
            "ASP": float(sara_asp),
        }

        # Ensure pct starts from >0 for interior-point distribution fitting.
        # IBP 0% is kept if present for plotting but stripped in the pipeline.
        with st.spinner("Running pipeline..."):
            try:
                result = _run_pipeline(
                    pct=pct_raw,
                    temps_K=temps_K_input,
                    method=method,
                    basis=basis,
                    sg_bulk=float(sg_bulk),
                    mw_bulk=float(mw_bulk),
                    sara=sara_dict,
                    n_quad=int(n_quad),
                    recovery_fraction=float(recovery_fraction),
                )
                ss["inputs"] = dict(
                    pct=pct_raw, temps_K=temps_K_input,
                    method=method, basis=basis, temp_unit=temp_unit,
                    sg_bulk=float(sg_bulk), mw_bulk=float(mw_bulk),
                    sara=sara_dict, n_quad=int(n_quad),
                    recovery_fraction=float(recovery_fraction),
                )
                ss["pipeline_result"] = result
                ss["pipeline_ready"]  = True

                elapsed = result["elapsed"]
                st.success(
                    f"Pipeline complete in {elapsed:.2f} s. "
                    f"Navigate to the other tabs to explore the results."
                )

            except Exception as exc:
                ss["pipeline_ready"] = False
                st.error(f"Pipeline error: {exc}")
                st.exception(exc)

    # ── Last-run summary (displayed when pipeline_ready) ──────────────────────
    if ss.get("pipeline_ready"):
        r = ss["pipeline_result"]
        st.divider()
        st.subheader("Last run — summary")

        c1, c2, c3, c4 = st.columns(4)
        with c1:
            st.metric("Quadrature nodes", r["n_quad"])
            st.metric("K_W (bulk)", f"{r['K_W_bulk']:.3f}")
        with c2:
            m_dev = abs(r["GL_M_av"] - r["dist_M_mean"]) / r["dist_M_mean"] * 100.0
            st.metric("GL M_av vs dist mean", f"{m_dev:.2f}%",
                      delta="PASS" if m_dev <= 0.5 else "CHECK",
                      delta_color="normal" if m_dev <= 0.5 else "inverse")
        with c3:
            sg_dev = abs(r["SG_av_full"] - r["sg_bulk"]) / r["sg_bulk"] * 100.0
            st.metric("SG_av vs SG_bulk", f"{sg_dev:.2f}%",
                      delta="PASS" if sg_dev < 5.0 else "CHECK",
                      delta_color="normal" if sg_dev < 5.0 else "inverse")
        with c4:
            flagged = r["kw_result"]["flagged"]
            st.metric("K_W-bin check", "FLAGGED" if flagged else "OK",
                      delta_color="inverse" if flagged else "normal")

        # ── Phase 11 heavy-resin diagnostics ──────────────────────────────────
        if r.get("has_heavy_resin"):
            st.markdown("**Heavy-resin lump (Phase 11 closure)**")
            h1, h2, h3, h4 = st.columns(4)
            with h1:
                st.metric("f_hr (wt%)", f"{r['f_hr']*100:.2f}%")
            with h2:
                st.metric("M_hr (g/mol)", f"{r['M_hr']:.1f}")
            with h3:
                st.metric("SG_hr", f"{r['SG_hr']:.4f}")
            with h4:
                st.metric("Tb_hr (K)", f"{r['Tb_hr']:.1f}")
            st.caption(
                f"Heavy-resin lump created (recovery = "
                f"{r['recovery_fraction']:.3f}).  MW and SG are set by mass-"
                f"balance closure on bulk MW = {r['mw_bulk']:.1f} g/mol and "
                f"bulk SG = {r['sg_bulk']:.4f}.  Tb_hr derived from constant "
                f"Watson K (Decision 31).  M_dist_av "
                f"(distillable subfraction) = {r['M_dist_av']:.1f} g/mol."
            )
        elif r.get("recovery_fraction", 1.0) < 1.0 - 1e-9:
            st.info(
                f"No heavy-resin lump created (f_hr ≈ 0, recovery + ASP ≈ 100 wt%); "
                f"pipeline dispatched to the Phase 8 path."
            )

        if r["kw_result"]["flagged"]:
            with st.expander("K_W-bin closure flags"):
                for msg in r["kw_result"]["flags"]:
                    st.warning(msg)
