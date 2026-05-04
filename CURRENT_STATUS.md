# petrochar — Current Status

## Project
Standalone Python tool for heavy petroleum fraction characterization. Methodology contribution for **Paper 1**: "Continuous-distribution characterization of heavy petroleum fractions for PC-SAFT phase-equilibrium modeling from routine refinery data". Target journals: *Energy & Fuels* or *Fuel*.

This is a separate repository from any other PC-SAFT project. Do not import, read, or reference code from outside this directory.

## Current Phase
Phase 0 — Repository scaffold (not yet started)

## Phases Completed
(none yet)

## Decisions Made (frozen)

1. **Distribution model:** Riazi generalized (Eq. 4.56) is the default. Gamma kept as legacy option only — fails on heavy oils per Riazi p. 178-179.
2. **Quadrature:** 5-point default + 1 ASP discrete + propane = 7 components total. 3-point selectable.
3. **Asphaltenes:** always discrete, Gonzalez 2007 defaults (m=33, σ=4.3 Å, ε/k=400 K). Tb=800°C and SG=1.15 are numerical conventions.
4. **Aromaticity γ:** per-component Watson K factor via linear clamp `γ = clamp((13.0 − K_W) / (13.0 − 9.5), 0, 1)`.
5. **SARA:** ASP wt% used as hard tail constraint. SAT/ARO/RES used as K_W-bin closure check only — never tuned to.
6. **No DAO data anywhere** in this repository.
7. **No reading external projects** (specifically `~/projects/pda_pcsaft_tool/` or any file outside this repo).
8. **Validation-first:** each phase requires reproducing a Riazi textbook example before proceeding.
9. **Repository name:** `petrochar`.
10. **Reference data format:** CSV files in `data/riazi_reference/` with citation header rows.

## Known Issues / Blockers

- **Riazi reference CSVs need to be vendored.** Five tables (4.6, 4.11, 4.13, 4.22, 4.23) need to be transcribed from the textbook into CSV. User to provide table values or page scans before Phase 1 can fully pass.

## Session Log

(empty — first session pending)
