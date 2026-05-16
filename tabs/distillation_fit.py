"""
Tab 2 — Distillation Fit.

Visualises the Tb distribution fit (Riazi Eq. 4.56) against measured
distillation data.  This is a DIAGNOSTIC tab only: the Tb distribution is
not on the critical path for pseudo-component property assignment (Decision 25,
Phase 8 rework).  Pseudo-component Tb_i is derived from riazi_daubert_Tb(M_i,
SG_i), not from evaluating this distribution at quadrature nodes.

Displays:
  - matplotlib plot: measured points vs fitted curve; extrapolation as dashed.
  - Fit quality metrics: RMS (K), AAD (%), R^2.
  - Fitted parameters: P_o, A, B.
  - Toggle to compare 3-param (free B) vs 2-param (B_T=1.5 fixed) fits.
  - Streamlit info box stating the diagnostic-only status.

References
----------
Riazi MNL50 Eq. 4.56; §4.5.4.1.
Decision 25 (Phase 8 rework): Tb_i from riazi_daubert_Tb(M_i, SG_i).
"""

from __future__ import annotations

import io

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st


def render(ss) -> None:
    """Render Tab 2 — Distillation Fit.

    Parameters
    ----------
    ss : streamlit.session_state
    """
    st.header("Distillation Fit")

    if not ss.get("pipeline_ready"):
        st.info("Run the pipeline in **Tab 1 — Input Data** first.")
        return

    r = ss["pipeline_result"]

    # ── Diagnostic-only notice ─────────────────────────────────────────────────
    st.info(
        "**Diagnostic only.** The Tb distribution (Riazi Eq. 4.56) fitted here "
        "characterises the feed's distillation curve but is **not** used to assign "
        "pseudo-component boiling points.  Pseudo-component Tb_i is derived from "
        "riazi_daubert_Tb(M_i, SG_i) via the self-consistent constant-Watson-K "
        "solve (Decision 25, Phase 8 rework).  Fit quality here does not affect "
        "the PC-SAFT parameters."
    )

    # ── Fit mode toggle ────────────────────────────────────────────────────────
    compare = st.checkbox("Compare 3-parameter vs 2-parameter (B_T = 1.5 fixed)", value=False)

    xc_int  = r["xc_int"]
    tb_int  = r["tb_int"]
    tb_3p   = r["tb_dist"]
    tb_2p   = r["tb_dist_2param"]
    rms_3p  = r["rms_tb"]
    rms_2p  = r["rms_tb_2p"]

    # Smooth curve range: from P_o slightly above to extrapolation past data.
    xc_smooth = np.linspace(0.001, 0.999, 300)

    # ── Fit quality table ──────────────────────────────────────────────────────
    fq_3p = tb_3p.fit_quality
    fq_2p = tb_2p.fit_quality

    st.subheader("Fitted parameters")
    if compare:
        col1, col2 = st.columns(2)
        _show_params(col1, "3-parameter (free B)", tb_3p, fq_3p)
        _show_params(col2, "2-parameter (B_T = 1.5)", tb_2p, fq_2p)
    else:
        col_p, col_q = st.columns(2)
        _show_params(col_p, "3-parameter (free B)", tb_3p, fq_3p)
        with col_q:
            st.write("")   # spacer

    # ── Plot ──────────────────────────────────────────────────────────────────
    st.subheader("Tb distribution fit")

    if compare:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=False)
        _plot_tb(axes[0], xc_int, tb_int, tb_3p, xc_smooth,
                 "3-parameter (free B)", rms_3p)
        _plot_tb(axes[1], xc_int, tb_int, tb_2p, xc_smooth,
                 "2-parameter (B_T = 1.5)", rms_2p)
        fig.tight_layout()
    else:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        _plot_tb(ax, xc_int, tb_int, tb_3p, xc_smooth,
                 "3-parameter (free B)", rms_3p)
        fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110)
    plt.close(fig)
    buf.seek(0)
    st.image(buf, use_container_width=True)

    # ── Method note ───────────────────────────────────────────────────────────
    with st.expander("About this fit (Riazi Eq. 4.56)"):
        st.markdown(
            r"""
            **Riazi generalized distribution (Eq. 4.56):**

            $$P(x_c) = P_o \left(1 + \left[\frac{A}{B} \ln\frac{1}{1-x_c}\right]^{1/B}\right)$$

            where $x_c$ is cumulative fraction (weight basis), $P_o$ is the onset
            parameter (property at $x_c \to 0$), $A$ is a scale parameter, and $B$ is
            a shape parameter.

            **3-parameter mode:** $P_o$, $A$, $B$ all free.
            **2-parameter mode:** $B$ fixed at 1.5 (Riazi p. 174 recommendation for Tb).

            The 2-parameter fit tends to overfit the middle range at the expense
            of the tails.  The 3-parameter free-$B$ fit usually gives a lower RMS
            but may extrapolate steeply beyond the measured range.

            **Why B_T = 1.5 is problematic for VTB:** heavy vacuum residues have
            steep Tb gradients above $x_c = 0.8$.  The fixed $B = 1.5$ forces a
            shape that underfits this region (RMS >> 5 K).  Use the 3-parameter
            mode for VTB feeds.
            """
        )


# ── Helpers ───────────────────────────────────────────────────────────────────

def _show_params(container, label: str, dist, fq: dict) -> None:
    with container:
        st.write(f"**{label}**")
        params_df = _params_table(dist, fq)
        st.table(params_df)


def _params_table(dist, fq: dict) -> "pd.DataFrame":
    import pandas as pd
    rows = [
        {"Parameter": "P_o (K)",   "Value": f"{dist.P_o:.2f}"},
        {"Parameter": "A",         "Value": f"{dist.A:.4f}"},
        {"Parameter": "B",         "Value": f"{dist.B:.4f}"},
        {"Parameter": "RMS (K)",   "Value": f"{fq['RMS']:.2f}"},
        {"Parameter": "AAD (%)",   "Value": f"{fq['pct_AAD']:.3f}"},
        {"Parameter": "R^2",       "Value": f"{fq['R_squared']:.6f}"},
    ]
    return pd.DataFrame(rows).set_index("Parameter")


def _plot_tb(ax, xc_int, tb_int, dist, xc_smooth, label: str, rms: float) -> None:
    """Draw measured points + fitted curve on ax."""
    # Measured data points
    ax.plot(
        xc_int * 100.0,
        tb_int - 273.15,
        "o", color="#1f77b4", markersize=6, label="Measured", zorder=5,
    )

    # Fitted curve — solid within data range, dashed extrapolation
    xc_min = xc_int.min()
    xc_max = xc_int.max()
    inside  = xc_smooth[(xc_smooth >= xc_min) & (xc_smooth <= xc_max)]
    outside_lo = xc_smooth[xc_smooth < xc_min]
    outside_hi = xc_smooth[xc_smooth > xc_max]

    if len(inside) > 0:
        tb_fit_in = np.array([dist.P(x) for x in inside]) - 273.15
        ax.plot(inside * 100.0, tb_fit_in, "-", color="#d62728",
                linewidth=1.8, label=f"Fitted (RMS={rms:.1f} K)")

    for xc_extra in (outside_lo, outside_hi):
        if len(xc_extra) > 0:
            tb_extra = np.array([dist.P(x) for x in xc_extra]) - 273.15
            ax.plot(xc_extra * 100.0, tb_extra, "--", color="#d62728",
                    linewidth=1.2, alpha=0.6)

    ax.set_xlabel("Cumulative fraction (%)", fontsize=10)
    ax.set_ylabel("Boiling temperature (degC)", fontsize=10)
    ax.set_title(label, fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
