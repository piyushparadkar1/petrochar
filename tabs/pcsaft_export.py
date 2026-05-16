"""
Tab 5 — PC-SAFT Export.

Aspen-Plus-paste-ready PC-SAFT parameter exports.

Three sub-sections:
  1. Pure-component PC-SAFT parameters (m, sigma, eps/k) for all 7 components
     (5 distillable + 1 ASP + propane).
  2. Component definitions (name, MW, Tb, SG).
  3. Binary interaction parameters — default 0.010 for all propane-pseudo-component
     pairs, 0 elsewhere; user-editable via data_editor.

Each sub-section has "Download CSV" and "Copy as TSV" (text area) buttons.

Brief Aspen Plus usage instructions included.

References
----------
Panuganti et al. (2012) Fuel 93, Table 6 — distillable correlations.
Gonzalez et al. (2007) Energy & Fuels 21 — asphaltene defaults.
Gross & Sadowski (2001) — propane sigma = 3.6180 A (not 3.168).
"""

from __future__ import annotations

import io

import numpy as np
import pandas as pd
import streamlit as st


# ── Component name formatter ──────────────────────────────────────────────────

def _comp_name(i: int, comp_type: str, M: float) -> str:
    """Generate a short component name for export."""
    if comp_type == "asphaltene":
        return "ASP"
    if comp_type == "propane":
        return "C3"
    return f"PC{i+1}"   # PC1, PC2, ... for distillable pseudo-components


def render(ss) -> None:
    """Render Tab 5 — PC-SAFT Export.

    Parameters
    ----------
    ss : streamlit.session_state
    """
    st.header("PC-SAFT Export")

    if not ss.get("pipeline_ready"):
        st.info("Run the pipeline in **Tab 1 — Input Data** first.")
        return

    r = ss["pipeline_result"]
    df_full = r["df_full"].copy()

    # ── Assign component names ────────────────────────────────────────────────
    names = []
    asp_count = 0
    for i, row in df_full.iterrows():
        ct = row["component_type"]
        if ct == "asphaltene":
            names.append("ASP")
            asp_count += 1
        elif ct == "propane":
            names.append("C3")
        else:
            names.append(f"PC{i+1}")
    df_full.insert(0, "component", names)

    n_comps = len(df_full)
    n_pc    = n_comps - 2   # distillable pseudo-components only

    # ── Aspen Plus usage instructions ─────────────────────────────────────────
    with st.expander("How to use these exports in Aspen Plus", expanded=False):
        st.markdown(
            """
            **Workflow for Aspen Plus PC-SAFT (POLYPCSF / PCSF) property method:**

            1. **Add components** — In Aspen Plus, go to *Components > Specifications*.
               Add one user-defined component for each row in the Component Definitions
               table (use the component name as the Aspen ID).  For propane (C3), use
               the built-in Aspen component if available.

            2. **Set PC-SAFT parameters** — Navigate to
               *Properties > Methods > Parameters > Pure Components*.
               For each component, enter the values from the Pure-Component Parameters
               table:
               - **PCSEG** = m (number of segments)
               - **PCDIA** = sigma (segment diameter, Angstrom)
               - **PCEDU** = eps/k (dispersion energy, K)

            3. **Set binary interaction parameters** — In
               *Properties > Methods > Parameters > Binary Interaction*.
               Enter the k_ij values from the Binary Interaction Parameters table.
               The default propane k_ij = 0.010 is a literature starting point;
               adjust if measured VLE data are available.

            4. **Download and paste** — Use the TSV format (tab-separated) to paste
               directly into Aspen spreadsheets.  The CSV format is for Excel import.

            > **Note on propane sigma:** The correct value is 3.6180 A.  The value
            > 3.168 A found in some Aspen built-in databases is a transcription error
            > (Gross & Sadowski 2001; Panuganti 2012 Table 5).
            """
        )

    # ── Section 1: Pure-component PC-SAFT parameters ─────────────────────────
    st.subheader("1 — Pure-component PC-SAFT parameters")

    pcsaft_table = df_full[["component", "M", "m", "sigma_A", "eps_over_k_K"]].copy()
    pcsaft_table.columns = ["Component", "MW (g/mol)", "m", "sigma (A)", "eps/k (K)"]

    st.dataframe(pcsaft_table, use_container_width=True, hide_index=True)
    _download_and_tsv(pcsaft_table, "pcsaft_pure_component_params",
                      "Section 1 TSV (copy and paste into Aspen Plus)")

    # ── Section 2: Component definitions ──────────────────────────────────────
    st.subheader("2 — Component definitions")

    comp_defs = df_full[["component", "M", "Tb_K", "SG"]].copy()
    comp_defs["Tb (degC)"] = (comp_defs["Tb_K"] - 273.15).round(1)
    comp_defs = comp_defs[["component", "M", "Tb (degC)", "Tb_K", "SG"]]
    comp_defs.columns = ["Component", "MW (g/mol)", "Tb (degC)", "Tb (K)", "SG"]

    # Replace NaN in Tb for propane with physical value (231.11 K = -42.04 degC)
    st.dataframe(comp_defs, use_container_width=True, hide_index=True)
    _download_and_tsv(comp_defs, "petrochar_component_definitions",
                      "Section 2 TSV")

    # ── Section 3: Binary interaction parameters ───────────────────────────────
    st.subheader("3 — Binary interaction parameters")

    st.caption(
        "Default: k_ij = 0.010 for all propane (C3) — pseudo-component pairs; "
        "0.000 for all other pairs.  Edit cells directly if needed."
    )

    comp_names = list(df_full["component"])
    n = len(comp_names)

    # Build default k_ij matrix
    kij_default = np.zeros((n, n), dtype=float)
    c3_idx = comp_names.index("C3") if "C3" in comp_names else None
    if c3_idx is not None:
        for j in range(n):
            if j != c3_idx:
                kij_default[c3_idx, j] = 0.010
                kij_default[j, c3_idx] = 0.010

    kij_df = pd.DataFrame(kij_default, index=comp_names, columns=comp_names)

    edited_kij = st.data_editor(
        kij_df,
        key="kij_editor",
        use_container_width=True,
        column_config={
            c: st.column_config.NumberColumn(c, format="%.4f")
            for c in comp_names
        },
    )

    # Reset index to make the component column explicit for export
    edited_kij_export = edited_kij.copy()
    edited_kij_export.insert(0, "Component", comp_names)
    _download_and_tsv(edited_kij_export, "petrochar_kij_matrix",
                      "Section 3 TSV (k_ij matrix)")


# ── Helpers ───────────────────────────────────────────────────────────────────

def _df_to_tsv(df: pd.DataFrame) -> str:
    """Convert DataFrame to tab-separated string."""
    return df.to_csv(sep="\t", index=False)


def _download_and_tsv(df: pd.DataFrame, stem: str, tsv_label: str) -> None:
    """Render Download CSV button and TSV text area."""
    col_dl, col_tsv = st.columns([1, 3])
    with col_dl:
        csv_bytes = df.to_csv(index=False).encode("utf-8")
        st.download_button(
            label="Download CSV",
            data=csv_bytes,
            file_name=f"{stem}.csv",
            mime="text/csv",
            key=f"dl_{stem}",
        )
    with col_tsv:
        with st.expander(tsv_label, expanded=False):
            tsv_text = _df_to_tsv(df)
            st.text_area(
                "Tab-separated (select all, copy, paste into Aspen)",
                value=tsv_text,
                height=150,
                key=f"tsv_{stem}",
                label_visibility="collapsed",
            )
