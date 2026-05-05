# CLAUDE.md — petrochar project rules

> **Read at the start of every Claude Code session.** Then read `CURRENT_STATUS.md` and the source brief `CLAUDE_CODE_PROMPT.md` if needed.

═══════════════════════════════════════════════════════════════════════════════
## FORBIDDEN PATHS — ABSOLUTE
═══════════════════════════════════════════════════════════════════════════════

**You may not read, reference, import from, or pattern-match against any code in:**

- `~/projects/pda_pcsaft_tool/`
- `../pda_pcsaft_tool/`
- Any folder named `pda*`, `pcsaft_pda*`
- Any file outside the petrochar repository root

That older project uses a **different architecture** (a 7-component "L/H" framework with MW tuning sliders) that this project deliberately abandons. Reading it corrupts design decisions through pattern-matching on wrong abstractions.

═══════════════════════════════════════════════════════════════════════════════
## FORBIDDEN CONCEPTS — DO NOT EXIST IN THIS PROJECT
═══════════════════════════════════════════════════════════════════════════════

If any of these appear in your output, that is contamination from training data. Delete and rewrite:

- `MW_ARO_L`, `MW_RES_L`, or any "tuning slider" variable
- "ARO-L / ARO-H" or "RES-L / RES-H" sub-fraction labels
- Harmonic-mean MW closure for back-calculating sub-fraction MWs
- 7-component framework with discrete L/H labels (petrochar's 7 components are 5 quadrature + 1 ASP + propane — different concept)
- DAO yield as a calibration target
- DAO SARA as a model input
- Any tuning of pseudo-component properties against process observables

**Self-check before every commit:**
```bash
grep -rni "MW_ARO_L\|MW_RES_L\|ARO-L\|ARO_L\|RES-L\|RES_L\|harmonic.mean\|DAO" \
    --include="*.py" --include="*.toml" --include="*.txt" --include="*.cfg" .
```
Zero hits required in committed code (excluding `.md` files which legitimately list these as forbidden).

═══════════════════════════════════════════════════════════════════════════════
## REFERENCE MATERIALS — IN THIS REPO
═══════════════════════════════════════════════════════════════════════════════

`data/riazi_reference/` contains:
- The full Riazi MNL50 PDF (`Riazi_MNL50_full.pdf`)
- 6 CSVs of textbook tables for validation tests
- 9 PNG page extracts at 200 DPI for visual verification

**Before implementing any Riazi correlation, view the corresponding page extract.** Do not implement equations from memory.

═══════════════════════════════════════════════════════════════════════════════
## CRITICAL CORRECTION — REGIME-DEPENDENT Tb CORRELATION
═══════════════════════════════════════════════════════════════════════════════

Riazi MNL50 §2.4.2.1 has TWO Riazi-Daubert Tb correlations:
- **Eq. 2.56** valid for M = 70-300
- **Eq. 2.57** valid for M = 300-700 (recommended for M > 300)

VTB pseudo-components span both ranges. Your `riazi_daubert_Tb` function must select correlation by input M. The same applies to `riazi_daubert_M` (inverse) and to SG correlations (Eqs. 2.59, 2.60). View `data/riazi_reference/page_077_riazi_daubert_Tb.png` before implementing.

═══════════════════════════════════════════════════════════════════════════════
## ACCEPTABLE PRIMARY SOURCES
═══════════════════════════════════════════════════════════════════════════════

- Riazi MNL50 (2005), especially Chapter 4 §4.5-4.6 and Chapter 2 §2.4
- Panuganti et al. (2012), *Fuel*, 93, 658-669
- Gonzalez et al. (2007), *Energy & Fuels*, 21, 1230-1234
- Watson and Nelson (1933), *Ind. Eng. Chem.*, 25, 880
- Whitson (1983), *SPE Journal*, 23(4), 683-694
- García Cárdenas and Ancheyta (2022), *Ind. Eng. Chem. Res.*, 61, 3383-3394

═══════════════════════════════════════════════════════════════════════════════
## DEVELOPMENT PROTOCOL
═══════════════════════════════════════════════════════════════════════════════

1. **No shortcuts. No stubs.** Every function fully implemented.
2. **Validation-first.** Each phase reproduces a Riazi textbook example to within stated tolerance before proceeding. Pass-gates non-negotiable.
3. **Citation density.** Every correlation, equation, constant cites Riazi MNL50 page/equation, Panuganti 2012 section, or Gonzalez 2007 section.
4. **Page-extract verification.** Before implementing any Riazi correlation, view the corresponding PNG in `data/riazi_reference/` with the `view` tool.
5. **Units explicit, always.** Tb in K internally; convert at I/O. SARA in wt%. Distillation curves carry an explicit basis flag.
6. **No proprietary names anywhere.** Generic only.
7. **Session protocol.** Start: read `CURRENT_STATUS.md`. End: update it, commit.

═══════════════════════════════════════════════════════════════════════════════
## ARCHITECTURE COMMITMENTS (FROZEN)
═══════════════════════════════════════════════════════════════════════════════

- **Distribution model:** Riazi generalized (Eq. 4.56). Gamma legacy-only.
- **Quadrature:** Riazi §4.6.1.1, default 5-point. Hard-coded points/weights from Table 4.21 (vendored CSV).
- **Tb correlation:** Regime-dependent — Eq. 2.56 for M ≤ 300, Eq. 2.57 for M > 300. **Critical correction; do not collapse to single equation.**
- **Asphaltenes:** always single discrete component. Gonzalez 2007 defaults: m=33, σ=4.3 Å, ε/k=400 K. Tb=800°C and SG=1.15 are numerical conventions.
- **Aromaticity γ:** per-pseudo-component Watson K, linear clamp `γ = clamp((13.0 − K_W) / (13.0 − 9.5), 0, 1)`.
- **PC-SAFT parameters:** Panuganti 2012 Table 6 for distillable; Gonzalez 2007 for asphaltenes; propane pure-component (m=2.002, **σ=3.6180 Å** [not 3.168 — beware Aspen typo], ε/k=208.11 K).
- **SARA closure:** ASP wt% is hard tail constraint. SAT/ARO/RES are K_W-bin closure checks only — never tuned to.

═══════════════════════════════════════════════════════════════════════════════
## EXECUTION SCOPE — ONE PHASE PER SESSION
═══════════════════════════════════════════════════════════════════════════════

User gives one instruction per session, naming the phase. Execute that phase only. Stop at phase's STOP marker. Do not start the next phase. User is not monitoring; you are responsible for self-policing.

═══════════════════════════════════════════════════════════════════════════════
## VERIFICATION BEFORE COMMIT
═══════════════════════════════════════════════════════════════════════════════

- Run all phase pass-gate tests. Do not commit if any test fails.
- After any change to `core/`: run full test suite via `pytest`.
- After any change to `app.py` or `tabs/`: `python -m py_compile app.py` and a `streamlit run app.py` boot test.
- Run forbidden-terms grep above. Zero non-md hits.
- Confirm no `import` references paths outside the repo.

═══════════════════════════════════════════════════════════════════════════════
## CONFLICT RESOLUTION
═══════════════════════════════════════════════════════════════════════════════

- `CLAUDE_CODE_PROMPT.md` (original brief) supersedes this file.
- This file supersedes any docstring or in-code comment.
- User instructions in chat supersede everything, but flag conflicts and propose updates.
