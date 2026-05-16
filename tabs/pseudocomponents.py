"""
Tab 4 — Pseudocomponents.

Displays the discrete pseudo-component table and K_W-bin closure report.

Table columns (per component):
    index, z, MW (g/mol), Tb (degC), Tb (K), SG, K_W, gamma,
    m, sigma (A), eps/k (K), is_asphaltene.

The asphaltene row is visually marked.

Closure report: ASP placement, sum check, per-class wt% (informational).

Download: CSV of the full table.

References
----------
Riazi MNL50 p. 75 (K_W bin thresholds).
Gonzalez 2007 (asphaltene defaults: M=1700, Tb=800 degC, SG=1.15).
Decision 28: K_W-bin test restricted to ASP+closure under constant K_W.
"""

from __future__ import annotations

import io

import numpy as np
import pandas as pd
import streamlit as st


def render(ss) -> None:
    """Render Tab 4 — Pseudocomponents.

    Parameters
    ----------
    ss : streamlit.session_state
    """
    st.header("Pseudocomponents")

    if not ss.get("pipeline_ready"):
        st.info("Run the pipeline in **Tab 1 — Input Data** first.")
        return

    r = ss["pipeline_result"]

    comps8    = r["comps8"]
    df_pcsaft = r["df"]

    # ── Build display table ───────────────────────────────────────────────────
    rows = []
    for i, c in enumerate(comps8):
        pc_row = df_pcsaft.iloc[i]
        rows.append({
            "Index":       i + 1,
            "Type":        "ASP" if c.is_asphaltene else f"D{i+1}",
            "z":           round(c.z,  6),
            "MW (g/mol)":  round(c.M,  2),
            "Tb (degC)":   round(c.Tb_K - 273.15, 1),
            "Tb (K)":      round(c.Tb_K, 1),
            "SG":          round(c.SG, 4),
            "K_W":         round(c.K_W, 3) if c.K_W is not None else float("nan"),
            "gamma":       round(c.gamma, 4) if c.gamma is not None else float("nan"),
            "m":           round(float(pc_row["m"]), 4),
            "sigma (A)":   round(float(pc_row["sigma_A"]), 4),
            "eps/k (K)":   round(float(pc_row["eps_over_k_K"]), 2),
            "is_asphaltene": c.is_asphaltene,
        })
    df_display = pd.DataFrame(rows)

    # ── Highlight asphaltene rows ─────────────────────────────────────────────
    def _highlight_asp(row):
        if row["is_asphaltene"]:
            return ["background-color: #fff3cd"] * len(row)
        return [""] * len(row)

    styled = df_display.style.apply(_highlight_asp, axis=1).format(
        {
            "z":           "{:.6f}",
            "MW (g/mol)":  "{:.2f}",
            "Tb (degC)":   "{:.1f}",
            "Tb (K)":      "{:.1f}",
            "SG":          "{:.4f}",
            "K_W":         "{:.3f}",
            "gamma":       "{:.4f}",
            "m":           "{:.4f}",
            "sigma (A)":   "{:.4f}",
            "eps/k (K)":   "{:.2f}",
        },
        na_rep="—",
    )

    st.subheader(f"Pseudo-component table ({len(comps8)} components)")
    st.caption("Asphaltene row highlighted in yellow.")
    st.dataframe(styled, use_container_width=True, hide_index=True)

    # ── Download CSV ──────────────────────────────────────────────────────────
    csv_bytes = df_display.to_csv(index=False).encode("utf-8")
    st.download_button(
        label="Download CSV",
        data=csv_bytes,
        file_name="petrochar_pseudocomponents.csv",
        mime="text/csv",
    )

    # ── K_W-bin closure report ────────────────────────────────────────────────
    st.subheader("K_W-bin closure report")

    kw_result = r["kw_result"]
    kw_calc   = kw_result["kw_calc"]
    sara_in   = kw_result["sara_input"]
    delta     = kw_result["delta_wt_pct"]

    closure_df = pd.DataFrame({
        "SARA class":         ["SAT", "ARO", "RES", "ASP"],
        "K_W-binned (wt%)":  [f"{kw_calc[k]:.2f}" for k in ("SAT", "ARO", "RES", "ASP")],
        "Measured (wt%)":    [f"{sara_in[k]:.2f}"  for k in ("SAT", "ARO", "RES", "ASP")],
        "Delta (wt%)":       [f"{delta[k]:+.2f}"   for k in ("SAT", "ARO", "RES", "ASP")],
    })
    st.table(closure_df.set_index("SARA class"))

    total_binned = sum(kw_calc.values())
    st.caption(
        f"Sum of K_W-binned wt% = {total_binned:.3f} wt%  "
        f"(should be 100.0 +/- 0.01)."
    )

    if kw_result["flagged"]:
        st.warning(
            "K_W-bin check flagged large deviations.  "
            "Under constant Watson K, all distillable components share K_W_bulk "
            "and fall in one bin — SAT/ARO/RES closure is degenerate by construction "
            "(Decision 28, Phase 8 rework).  This is informational, not an error."
        )
        with st.expander("Flag details"):
            for msg in kw_result["flags"]:
                st.text(msg)
    else:
        st.success("K_W-bin check: no large deviations flagged.")

    # ── Watson K summary ──────────────────────────────────────────────────────
    with st.expander("Watson K and aromaticity summary"):
        kw_vals = [c.K_W for c in comps8 if not c.is_asphaltene]
        g_vals  = [c.gamma for c in comps8 if not c.is_asphaltene]
        if kw_vals:
            st.markdown(
                f"All distillable components share **K_W = {kw_vals[0]:.3f}** "
                f"(constant Watson K method).  "
                f"Corresponding aromaticity: **gamma = {g_vals[0]:.4f}**  "
                f"(linear clamp: gamma = clamp((13.0 - K_W) / (13.0 - 9.5), 0, 1)).  "
                f"ASP K_W = {comps8[-1].K_W:.3f}, ASP gamma = {comps8[-1].gamma:.4f}."
            )
