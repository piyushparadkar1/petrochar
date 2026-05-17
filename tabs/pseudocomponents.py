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
    dist_counter = 0
    for i, c in enumerate(comps8):
        pc_row = df_pcsaft.iloc[i]
        # Type tag: ASP for asphaltene, HR for heavy-resin lump, D{n} for
        # distillable quadrature nodes (Phase 11 type taxonomy).
        if c.is_asphaltene:
            type_tag = "ASP"
        elif c.is_heavy_resin:
            type_tag = "HR"
        else:
            dist_counter += 1
            type_tag = f"D{dist_counter}"
        rows.append({
            "Index":       i + 1,
            "Type":        type_tag,
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
            "is_asphaltene":  c.is_asphaltene,
            "is_heavy_resin": c.is_heavy_resin,
        })
    df_display = pd.DataFrame(rows)

    # ── Highlight asphaltene rows (yellow) and heavy-resin row (light blue) ──
    def _highlight_row(row):
        if row["is_asphaltene"]:
            return ["background-color: #fff3cd"] * len(row)
        if row["is_heavy_resin"]:
            return ["background-color: #d6ecf3"] * len(row)
        return [""] * len(row)

    styled = df_display.style.apply(_highlight_row, axis=1).format(
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
    has_hr = any(c.is_heavy_resin for c in comps8)
    if has_hr:
        st.caption(
            "Asphaltene row highlighted in **yellow**; heavy-resin lump (Phase 11) "
            "highlighted in **light blue**."
        )
    else:
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
        # Distillable components only (exclude both ASP and HR).
        dist_only = [c for c in comps8
                     if not c.is_asphaltene and not c.is_heavy_resin]
        asp_only  = [c for c in comps8 if c.is_asphaltene]
        hr_only   = [c for c in comps8 if c.is_heavy_resin]
        if dist_only:
            kw_d = dist_only[0].K_W
            g_d  = dist_only[0].gamma
            msg = (
                f"All distillable components share **K_W = {kw_d:.3f}** "
                f"(constant Watson K method).  "
                f"Corresponding aromaticity: **gamma = {g_d:.4f}**  "
                f"(linear clamp: gamma = clamp((13.0 - K_W) / (13.0 - 9.5), 0, 1)).  "
            )
            if hr_only:
                msg += (
                    f"\n\n**Heavy-resin lump:** K_W = {hr_only[0].K_W:.3f}, "
                    f"gamma = {hr_only[0].gamma:.4f}.  "
                    f"(HR shares K_W_bulk by construction under constant Watson K, "
                    f"Decision 31.)"
                )
            if asp_only:
                msg += (
                    f"\n\n**Asphaltene:** K_W = {asp_only[0].K_W:.3f}, "
                    f"gamma = {asp_only[0].gamma:.4f}.  "
                    f"(Numerical conventions Tb = 1073.15 K, SG = 1.15.  ASP "
                    f"K_W is recorded only — PC-SAFT uses Gonzalez 2007 defaults.)"
                )
            st.markdown(msg)
