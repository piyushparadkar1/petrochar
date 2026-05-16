"""
petrochar tab modules — one file per Streamlit tab.

Each module exposes a single render() function that accepts the Streamlit
session state and renders that tab.  Tabs read from session state; only
Tab 1 (input_data) writes the initial inputs and triggers the pipeline.

Tab layout (exactly six — scope-creep rule enforced):
    1. input_data       — user inputs + pipeline trigger
    2. distillation_fit — Tb distribution diagnostic (informational)
    3. distributions    — M and SG distribution PDFs/CDFs + quadrature nodes
    4. pseudocomponents — discrete pseudo-component table + K_W-bin closure
    5. pcsaft_export    — Aspen-Plus-ready PC-SAFT parameter export
    6. validation       — in-app Riazi textbook example reproduction
"""
