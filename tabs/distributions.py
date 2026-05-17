"""
Tab 3 — Distributions.

Visualises the M and SG distributions (PDF and CDF) derived from the
pipeline.  Quadrature node positions are marked on each CDF plot.

Layout: 2x2 matplotlib subplots (PDF of M, CDF of M, PDF of SG, CDF of SG).
Node positions shown as vertical lines with labels on CDF panels.

For the constant Watson K SG method: displays K_W_bulk prominently and notes
that all SG values are derived from this constant.

References
----------
Riazi MNL50 Eq. 4.56 (§4.5.4.1) — generalized distribution.
Riazi MNL50 Eq. 2.13 — Watson K factor.
Decision 25: Tb_i derived from riazi_daubert_Tb(M_i, SG_i).
"""

from __future__ import annotations

import io

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import streamlit as st


def render(ss) -> None:
    """Render Tab 3 — Distributions.

    Parameters
    ----------
    ss : streamlit.session_state
    """
    st.header("Distributions")

    if not ss.get("pipeline_ready"):
        st.info("Run the pipeline in **Tab 1 — Input Data** first.")
        return

    r   = ss["pipeline_result"]
    inp = ss["inputs"]

    # ── SG method note ────────────────────────────────────────────────────────
    st.info(
        f"**Constant Watson K** method.  "
        f"Bulk Watson K factor: K_W = **{r['K_W_bulk']:.3f}**  "
        f"(from bulk SG = {inp['sg_bulk']:.3f}, bulk MW = {inp['mw_bulk']:.1f} g/mol).  "
        "All pseudo-component SG values are derived algebraically from their "
        "individual Tb via SG_i = (1.8 × Tb_i)^(1/3) / K_W."
    )

    # ── Phase 11 recovery note ────────────────────────────────────────────────
    rec_frac = r.get("recovery_fraction", 1.0)
    if rec_frac < 1.0 - 1e-9:
        st.warning(
            f"**Recovery-aware fit (Phase 11):** the distillation curve covers "
            f"only **{rec_frac*100:.1f}%** of the total feed mass.  The M and "
            f"SG distributions are fitted on a **scaled** xc basis "
            f"(xc_scaled = xc_raw / {rec_frac:.3f}), so the GL quadrature "
            f"samples the distillable subfraction only.  The unmeasured tail "
            f"is represented by the discrete heavy-resin lump (see Tab 4)."
        )

    # ── Extract data ──────────────────────────────────────────────────────────
    m_dist     = r["m_dist"]
    xc_int     = r["xc_int"]
    M_cuts     = r["M_cuts"]
    SG_cuts    = r["SG_cuts"]
    comps_dist = r["comps_dist"]
    K_W_bulk   = r["K_W_bulk"]

    # M range for plot: slightly below fitted min, up to 1.3x max node
    M_nodes = np.array([c.M for c in comps_dist])
    SG_nodes = np.array([c.SG for c in comps_dist])
    z_nodes  = np.array([c.z for c in comps_dist])

    M_lo = max(m_dist.P_o * 1.01, M_nodes.min() * 0.85)
    M_hi = M_nodes.max() * 1.15
    M_range = np.linspace(M_lo, M_hi, 400)

    # CDF for M
    M_cdf  = np.array([m_dist.x_c(m) for m in M_range])
    # PDF for M
    M_pdf  = np.array([m_dist.pdf(m) for m in M_range])

    # SG range: derive from Tb range via Watson K inversion
    # SG_i = (1.8*Tb_i)^(1/3)/K_W; Tb range from M range via distribution
    SG_lo = SG_nodes.min() * 0.93
    SG_hi = SG_nodes.max() * 1.07
    SG_range = np.linspace(SG_lo, SG_hi, 400)

    # For SG distribution: build a simple analytical SG CDF from M distribution.
    # Under constant K_W: Tb = riazi_daubert_Tb(M, SG(M)), and
    # SG_i = (1.8*Tb_i)^(1/3)/K_W.  We obtain the SG CDF numerically by
    # mapping each measured (xc, SG_cut) pair.
    xc_sg = xc_int.copy()
    sg_fit = SG_cuts.copy()

    # ── Figure ────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    ax_mpdf, ax_mcdf, ax_sgpdf, ax_sgcdf = (
        axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]
    )

    node_colors = plt.cm.tab10(np.linspace(0, 0.5, len(comps_dist)))

    # ── M PDF ─────────────────────────────────────────────────────────────────
    ax_mpdf.plot(M_range, M_pdf, "-", color="#1f77b4", linewidth=1.8)
    for i, (m_node, sg_node, z) in enumerate(zip(M_nodes, SG_nodes, z_nodes)):
        ax_mpdf.axvline(m_node, color=node_colors[i], linestyle="--",
                        linewidth=1.2, alpha=0.7)
    ax_mpdf.set_xlabel("M (g/mol)", fontsize=10)
    ax_mpdf.set_ylabel("f(M)", fontsize=10)
    ax_mpdf.set_title("PDF — Molecular Weight", fontsize=10)
    ax_mpdf.grid(True, alpha=0.3)

    # ── M CDF ─────────────────────────────────────────────────────────────────
    ax_mcdf.plot(M_range, M_cdf, "-", color="#1f77b4", linewidth=1.8,
                 label="Fitted CDF")
    # Measured data points on CDF
    ax_mcdf.plot(M_cuts, xc_int, "s", color="#2ca02c", markersize=5,
                 label="Measured cuts", zorder=4)
    for i, (m_node, z) in enumerate(zip(M_nodes, z_nodes)):
        xc_node = float(m_dist.x_c(m_node))
        ax_mcdf.axvline(m_node, color=node_colors[i], linestyle="--",
                        linewidth=1.2, alpha=0.9,
                        label=f"Node {i+1} (z={z:.3f})")
        ax_mcdf.plot(m_node, xc_node, "^", color=node_colors[i],
                     markersize=8, zorder=5)
    # Mark recovery limit as a horizontal dashed line (Phase 11).
    if rec_frac < 1.0 - 1e-9:
        ax_mcdf.axhline(
            1.0, color="#888888", linestyle=":", linewidth=1.0,
            label="distillable subfraction = 1.0",
        )
    ax_mcdf.set_xlabel("M (g/mol)", fontsize=10)
    ax_mcdf.set_ylabel("Cumulative fraction (within distillable subfraction)"
                        if rec_frac < 1.0 - 1e-9 else "Cumulative fraction",
                        fontsize=10)
    ax_mcdf.set_title("CDF — Molecular Weight", fontsize=10)
    ax_mcdf.legend(fontsize=7, loc="upper left")
    ax_mcdf.grid(True, alpha=0.3)

    # ── SG PDF (numeric from measured cuts via kernel density) ────────────────
    # Use linear interpolation of measured SG vs xc to estimate a PDF shape.
    # dxc/dSG ≈ 1/slope between adjacent measured SG cuts.
    if len(sg_fit) > 2:
        dxc_dSG = np.gradient(xc_sg, sg_fit)
        ax_sgpdf.plot(sg_fit, dxc_dSG, "o-", color="#ff7f0e",
                      markersize=5, linewidth=1.5,
                      label="dx_c/d(SG) approx.")
    else:
        ax_sgpdf.text(0.5, 0.5, "Insufficient data", ha="center",
                      va="center", transform=ax_sgpdf.transAxes)
    for i, (sg_node, z) in enumerate(zip(SG_nodes, z_nodes)):
        ax_sgpdf.axvline(sg_node, color=node_colors[i], linestyle="--",
                         linewidth=1.2, alpha=0.7)
    ax_sgpdf.set_xlabel("SG", fontsize=10)
    ax_sgpdf.set_ylabel("dx_c/d(SG)", fontsize=10)
    ax_sgpdf.set_title("PDF — Specific Gravity (constant K_W)", fontsize=10)
    ax_sgpdf.grid(True, alpha=0.3)

    # ── SG CDF ────────────────────────────────────────────────────────────────
    ax_sgcdf.plot(sg_fit, xc_sg, "o-", color="#ff7f0e",
                  markersize=5, linewidth=1.5, label="Measured cuts")
    for i, (sg_node, z) in enumerate(zip(SG_nodes, z_nodes)):
        ax_sgcdf.axvline(sg_node, color=node_colors[i], linestyle="--",
                         linewidth=1.2, alpha=0.9,
                         label=f"Node {i+1} (z={z:.3f})")
    ax_sgcdf.set_xlabel("SG", fontsize=10)
    ax_sgcdf.set_ylabel("Cumulative fraction", fontsize=10)
    ax_sgcdf.set_title("CDF — Specific Gravity (constant K_W)", fontsize=10)
    ax_sgcdf.legend(fontsize=7, loc="upper left")
    ax_sgcdf.grid(True, alpha=0.3)

    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110)
    plt.close(fig)
    buf.seek(0)
    st.image(buf, use_container_width=True)

    # ── Analytic mean summary ─────────────────────────────────────────────────
    st.subheader("Distribution summary")
    import pandas as pd
    params_data = {
        "Parameter": ["P_o (g/mol)", "A", "B", "Analytic mean (g/mol)",
                       "GL M_av (g/mol)", "M_av vs mean (%)"],
        "Value": [
            f"{m_dist.P_o:.2f}",
            f"{m_dist.A:.4f}",
            f"{m_dist.B:.4f}",
            f"{r['dist_M_mean']:.2f}",
            f"{r['GL_M_av']:.2f}",
            f"{abs(r['GL_M_av'] - r['dist_M_mean']) / r['dist_M_mean'] * 100.0:.3f}%",
        ],
    }
    st.table(pd.DataFrame(params_data).set_index("Parameter"))

    # ── SG note ───────────────────────────────────────────────────────────────
    st.caption(
        "SG values are derived from constant Watson K: "
        f"SG_i = (1.8 × Tb_i)^(1/3) / {K_W_bulk:.3f}.  "
        "No separate SG distribution is fitted for this feed."
    )
