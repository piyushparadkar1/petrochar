# CLAUDE.md — petrochar project rules

## Read this at the start of every Claude Code session.

This file is the project-specific equivalent of `~/.claude/CLAUDE.md`. It establishes rules and design decisions that apply to all sessions on this repository.

If you want full context, also read `CLAUDE_CODE_PROMPT.md` (the original Session 0 brief) and `CURRENT_STATUS.md` (the running session log).

---

## STRICT BOUNDARIES

1. **Do NOT read or reference any file outside this repository.** There is an older project at `~/projects/pda_pcsaft_tool/` (or similar). You are forbidden from reading it. It contains a different architecture (a 7-component "L/H" framework) that this project deliberately moves away from. Reading it would corrupt your design decisions.

2. **Build everything from primary sources only.** Acceptable references:
   - Riazi, M. R. (2005). *Characterization and Properties of Petroleum Fractions*. ASTM Manual MNL50.
   - Panuganti, S. R. et al. (2012). *Fuel*, 93, 658-669.
   - Gonzalez, D. L. et al. (2007). *Energy & Fuels*, 21, 1230-1234.
   - Watson, K. M. and Nelson, E. F. (1933). *Ind. Eng. Chem.*, 25, 880.
   - Whitson, C. H. (1983). *SPE Journal*, 23(4), 683-694.
   - García Cárdenas, J. and Ancheyta, J. (2022). *Ind. Eng. Chem. Res.*, 61, 3383-3394.

3. **No DAO data. No process-side data. Ever.** This tool's contract is feed-side characterization only.

4. **No phase-equilibrium computation.** No flashes, no LLE solvers, no chemical-potential equality root-finding.

5. **No tuning against process observables.** All pseudo-component properties come from feed measurements and published correlations.

---

## DEVELOPMENT PROTOCOL

1. **No shortcuts. No stubs.** Every function fully implemented. No `pass`, no `NotImplementedError`, no commented-out fallbacks.

2. **Validation-first.** Each phase requires reproducing a Riazi textbook example to within stated tolerance before proceeding. **Pass-gates are non-negotiable.**

3. **Citation density.** Every correlation, equation, and constant has a docstring citing source: Riazi MNL50 page/equation, Panuganti 2012 section, Gonzalez 2007 section.

4. **Units explicit, always.** Every function signature documents input and output units. Tb in K internally; convert at I/O. SARA in wt%. Distillation curves carry an explicit basis flag.

5. **No proprietary names anywhere.** No client refinery, no plant numbers, no internal codes. Use generic terms: "VTB sample", "heavy petroleum fraction".

6. **Session protocol.**
   - Start: read `CURRENT_STATUS.md`. Read only files needed for current phase.
   - End: update `CURRENT_STATUS.md` with phase completion, decisions, blockers, session entry. Commit.

---

## ARCHITECTURE COMMITMENTS

- **Distribution model:** Riazi generalized model (Eq. 4.56), not gamma. Gamma fails on heavy oils per Riazi p. 178-179. Gamma kept as legacy option only.
- **Quadrature:** Riazi §4.6.1.1, default 5-point. 3-point selectable. Hard-coded points/weights from Table 4.21.
- **Asphaltenes:** Always a single discrete component. Gonzalez 2007 defaults: m=33, σ=4.3 Å, ε/k=400 K. Tb=800°C and SG=1.15 are numerical conventions, documented as such.
- **Aromaticity factor γ:** computed per-pseudo-component from per-component Watson K factor via linear clamp `γ = clamp((13.0 − K_W) / (13.0 − 9.5), 0, 1)`. Cite Watson 1933, Riazi MNL50 §2.1.17, Panuganti 2012 §2.3.
- **PC-SAFT parameters:** Panuganti 2012 Table 6 for distillable pseudo-components, Gonzalez 2007 for asphaltenes, fixed pure-component values for propane (m=2.002, **σ=3.6180 Å** [not 3.168, beware of Aspen typo], ε/k=208.11 K).
- **SARA closure:** ASP wt% is a hard tail constraint on the distillation distribution. SAT/ARO/RES are K_W-bin closure checks only (not constraints; not tuned to).

---

## VERIFICATION BEFORE COMMIT

- Each phase: run all phase pass-gate tests. Do not commit if any test fails.
- After any change to `core/`: run full test suite via `pytest`.
- After any change to `app.py` or `tabs/`: `python -m py_compile app.py` plus a `streamlit run app.py` boot test.
- Search for proprietary terms with `grep -ni` before any commit: `HPCL`, `KBR`, `Plant 41`, `Mumbai`, plus any other client-specific terms surfaced during dev. Must return zero hits in committed files.
- Confirm no `import` statements reference paths outside this repo.

---

## WHEN INSTRUCTIONS CONFLICT

- If `CLAUDE.md` (this file) conflicts with a docstring or in-code comment: this file wins.
- If `CLAUDE_CODE_PROMPT.md` conflicts with this file: `CLAUDE_CODE_PROMPT.md` wins (it is the original brief; this file is a derived summary).
- If user instructions in chat conflict with this file: follow the user, but flag the conflict explicitly and propose an update to this file or `CLAUDE_CODE_PROMPT.md`.
