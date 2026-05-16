"""
petrochar — Streamlit entry point.

Responsibilities:
  - Page configuration (title, layout).
  - Session-state key initialisation.
  - Tab routing: instantiate the six tabs and call render() on the active one.
  - No scientific computation lives here.  All computation is in core/.

Usage:
    streamlit run app.py

References
----------
Methodology: see docs/validation_report.md and docs/methodology.md.
Phases 1-8:  all scientific logic in core/*.py, fully tested in tests/*.py.
Phase 9:     this file + tabs/*.py (UI wrapper only).
"""

import streamlit as st

# ── Page configuration (must be the first Streamlit call) ────────────────────
st.set_page_config(
    page_title="petrochar",
    page_icon="🛢",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# ── Session-state initialisation ──────────────────────────────────────────────
# One source of truth for all pipeline data.
# Tab 1 populates these; Tabs 2-6 read them.
#
# Keys:
#   inputs          dict — raw user inputs (distillation curve, SG, MW, SARA, options)
#   pipeline_result dict — complete pipeline output produced by Tab 1 on "Run pipeline"
#   pipeline_ready  bool — True once pipeline has run successfully

_DEFAULTS = {
    "inputs":          None,
    "pipeline_result": None,
    "pipeline_ready":  False,
}

for _key, _val in _DEFAULTS.items():
    if _key not in st.session_state:
        st.session_state[_key] = _val

# ── Tab imports ───────────────────────────────────────────────────────────────
from tabs import input_data         # noqa: E402
from tabs import distillation_fit   # noqa: E402
from tabs import distributions      # noqa: E402
from tabs import pseudocomponents   # noqa: E402
from tabs import pcsaft_export      # noqa: E402
from tabs import validation         # noqa: E402

# ── Header ────────────────────────────────────────────────────────────────────
st.title("petrochar")
st.caption(
    "Continuous-distribution characterization of heavy petroleum fractions "
    "for PC-SAFT phase-equilibrium modeling."
)

# ── Tab routing ───────────────────────────────────────────────────────────────
_tab_labels = [
    "1 — Input Data",
    "2 — Distillation Fit",
    "3 — Distributions",
    "4 — Pseudocomponents",
    "5 — PC-SAFT Export",
    "6 — Validation",
]

(
    tab_input,
    tab_distfit,
    tab_dists,
    tab_pseudo,
    tab_pcsaft,
    tab_valid,
) = st.tabs(_tab_labels)

with tab_input:
    input_data.render(st.session_state)

with tab_distfit:
    distillation_fit.render(st.session_state)

with tab_dists:
    distributions.render(st.session_state)

with tab_pseudo:
    pseudocomponents.render(st.session_state)

with tab_pcsaft:
    pcsaft_export.render(st.session_state)

with tab_valid:
    validation.render(st.session_state)
