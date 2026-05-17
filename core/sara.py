"""
SARA closure check and discrete component assembly.

Responsibilities
----------------
1. Append the asphaltene pseudo-component (discrete, literature-defined).
2. Append the heavy-resin lump + asphaltene (recovery-aware, Phase 11).
3. K_W-bin aggregation of distillable pseudo-components for SARA closure check.

ARCHITECTURAL COMMITMENT (frozen):
    SAT/ARO/RES wt% from SARA are closure checks ONLY.  They are never used
    to adjust pseudo-component properties.  If K_W-binned wt% deviate from
    measured SARA, the deviation is flagged for data inspection.
    Distribution parameters are NEVER retuned to close the SARA balance.

    ASP wt% is the only SARA input that directly affects pseudo-component
    construction: it sets the mass fraction of the discrete asphaltene component.

Mole fraction convention
------------------------
After append_asphaltene() OR append_heavy_resin_and_asphaltene() is called,
ALL Pseudocomponent.z values are true mole fractions in the full mixture
(distillable [+ heavy-resin lump] + asphaltene), summing to 1.

Phase 8 path (append_asphaltene, recovery_fraction = 1.0):
    f_asp    = asp_wt_pct / 100
    M_dist_av = sum_i(z_i * M_i)          (mole-avg MW of distillable)
    n_dist   = (1 - f_asp) / M_dist_av
    n_asp    = f_asp / M_asp
    n_tot    = n_dist + n_asp
    z_i_true = z_i * n_dist / n_tot
    z_asp    = n_asp / n_tot

Phase 11 path (append_heavy_resin_and_asphaltene, recovery_fraction < 1.0):
    f_dist   = recovery_fraction           (mass fraction covered by D1160)
    f_asp    = asp_wt_pct / 100
    f_hr     = 1 - f_dist - f_asp          (unmeasured tail, >= 0)

    Bulk MW closure (number-average):
        M_hr = f_hr / (1/M_bulk - f_dist/M_dist_av - f_asp/M_asp)

    Bulk SG closure (volume-additive):
        SG_hr = f_hr / (1/SG_bulk - f_dist/SG_dist_vavg - f_asp/SG_asp)
        SG_dist_vavg = 1 / sum_i_dist(xw_i_dist / SG_i)
        xw_i_dist    = z_i * M_i / M_dist_av

    Heavy-resin Tb under constant Watson K (Decision 31):
        Tb_hr = (K_W_bulk * SG_hr)^3 / 1.8
        Capped at 1100 K with UserWarning if exceeded.

    Moles per kg total feed:
        n_dist   = f_dist / M_dist_av
        n_hr     = f_hr   / M_hr
        n_asp    = f_asp  / M_asp
        n_total  = n_dist + n_hr + n_asp

    Full-mixture mole fractions:
        z_i_full  = z_i * n_dist / n_total    (distillable)
        z_hr_full = n_hr / n_total            (heavy-resin lump)
        z_asp_full= n_asp / n_total           (asphaltene)

Decision 30 (Phase 11): The M distribution characterises the distillable
    subfraction only (mass = recovery_fraction). The unmeasured tail mass
    (1 - recovery_fraction - f_asp) is represented by a single discrete
    heavy-resin lump whose M and SG are fixed by mass-balance closure.
    When recovery_fraction = 1.0, the pipeline reduces exactly to Phase 8.

Decision 31 (Phase 11): Heavy-resin lump Watson K equals K_W_bulk by
    construction. Tb_hr = (K_W_bulk * SG_hr)^3 / 1.8. Cap at 1100 K.

The xc_lower/xc_upper fields remain on a weight basis (consistent with
Phase 5's cumulative weight-fraction convention).

References
----------
Gonzalez et al. (2007)  : nanoaggregate defaults M=1700 g/mol, SG=1.15.
                          Energy & Fuels, 21, 1230-1234.
Riazi (2005)            : K_W bin thresholds, p. 75.
"""

from __future__ import annotations

import math
import warnings

import numpy as np

from core.correlations import watson_k
from core.quadrature import Pseudocomponent


# ── K_W bin thresholds (Riazi p. 75 conventions — user-overridable) ──────────

KW_SAT_DEFAULT = 12.0   # K_W >= KW_SAT_DEFAULT → SAT
KW_ARO_DEFAULT = 11.0   # KW_ARO_DEFAULT <= K_W < KW_SAT_DEFAULT → ARO
                        # K_W < KW_ARO_DEFAULT → RES

# ── Asphaltene conventions (Gonzalez 2007) ────────────────────────────────────

ASP_M_DEFAULT    = 1700.0    # g/mol; nanoaggregate (Gonzalez 2007)
ASP_TB_K_DEFAULT = 1073.15   # K; 800 deg C — numerical convention only.
                              # Asphaltenes do not have a physical boiling
                              # point; this value ensures ASP sorts above all
                              # distillable pseudo-components (all Tb < 900 K
                              # for normal petroleum fractions).
ASP_SG_DEFAULT   = 1.15      # at 15.5 degC (Gonzalez 2007)


# ── Validation ────────────────────────────────────────────────────────────────

def validate_sara(
    sat: float,
    aro: float,
    res: float,
    asp: float,
    tol: float = 0.5,
) -> None:
    """Validate SARA wt% inputs.

    Parameters
    ----------
    sat, aro, res, asp : float
        Saturates, aromatics, resins, asphaltenes in wt%.
    tol : float
        Allowable deviation from 100 wt% sum (default 0.5 wt%).

    Raises
    ------
    ValueError
        If any component is negative, exceeds 100, or the sum deviates from
        100 by more than tol.

    References
    ----------
    Riazi MNL50 §4.2.3 (SARA analysis definition, p. 73-75).
    """
    labels = ('SAT', 'ARO', 'RES', 'ASP')
    values = (sat, aro, res, asp)
    for label, val in zip(labels, values):
        if val < 0.0:
            raise ValueError(
                f"SARA wt%: {label} = {val:.4f} is negative."
            )
        if val > 100.0:
            raise ValueError(
                f"SARA wt%: {label} = {val:.4f} exceeds 100 wt%."
            )
    total = sat + aro + res + asp
    if abs(total - 100.0) > tol:
        raise ValueError(
            f"SARA wt% do not sum to 100: "
            f"SAT={sat:.3f} + ARO={aro:.3f} + RES={res:.3f} + ASP={asp:.3f} = {total:.3f}."
        )


# ── Asphaltene component ──────────────────────────────────────────────────────

def append_asphaltene(
    components: list[Pseudocomponent],
    asp_wt_pct: float,
    M_asp: float = ASP_M_DEFAULT,
    Tb_K_asp: float = ASP_TB_K_DEFAULT,
    SG_asp: float = ASP_SG_DEFAULT,
) -> list[Pseudocomponent]:
    """Append a discrete asphaltene pseudo-component to the component list.

    The asphaltene is not derived from the distillation curve.  It is a
    single literature-defined component with Gonzalez 2007 nanoaggregate
    conventions.  The appended component occupies the heaviest tail:
    x_cw_lower = 1 - asp_wt_pct/100, x_cw_upper = 1.0.

    All z values in the returned list are TRUE MOLE FRACTIONS in the full
    mixture (distillable + asphaltene), summing to 1.  The input distillable
    components carry z_i = mole fraction within distillable subfraction; this
    function converts them to full-mixture mole fractions before returning.
    See module docstring for the conversion formula.

    Parameters
    ----------
    components : list[Pseudocomponent]
        Distillable pseudo-components from discretize_generalized.
        Each component's z must be a mole fraction within the distillable
        subfraction (z values should sum to 1).
        Not modified in place; a new list is returned.
    asp_wt_pct : float
        Asphaltene content in wt%.  Must be in the open interval (0, 100).
    M_asp : float
        Asphaltene MW in g/mol.  Default 1700 (Gonzalez 2007 nanoaggregate).
    Tb_K_asp : float
        Boiling-point convention in K.  Default 1073.15 K (800 deg C).
        Not a physical boiling point — see ASP_TB_K_DEFAULT note above.
    SG_asp : float
        Specific gravity at 15.5 degC.  Default 1.15 (Gonzalez 2007).

    Returns
    -------
    list[Pseudocomponent]
        New list of N+1 components with z values as true mole fractions in
        the full mixture.  Original input list is not modified.

    References
    ----------
    Gonzalez et al. (2007) Energy & Fuels, 21, 1230-1234.
    """
    if asp_wt_pct <= 0.0 or asp_wt_pct >= 100.0:
        raise ValueError(
            f"asp_wt_pct must be in the open interval (0, 100); got {asp_wt_pct}."
        )

    f_asp = asp_wt_pct / 100.0

    # Mole-average MW of the distillable subfraction
    # (z_i sum to 1 within distillable, so sum(z_i * M_i) = M_dist_av)
    M_dist_av = float(sum(c.z * c.M for c in components))

    # Moles of each subfraction per kg of total feed
    n_dist = (1.0 - f_asp) / M_dist_av
    n_asp  = f_asp          / M_asp
    n_tot  = n_dist + n_asp

    scale_dist = n_dist / n_tot    # multiplier converting distillable z to full-mixture z
    z_asp      = n_asp  / n_tot

    new_components = [
        Pseudocomponent(
            z=c.z * scale_dist,
            M=c.M,
            Tb_K=c.Tb_K,
            SG=c.SG,
            xc_lower=c.xc_lower,
            xc_upper=c.xc_upper,
            is_asphaltene=c.is_asphaltene,    # preserve any pre-existing flag
            is_heavy_resin=c.is_heavy_resin,  # preserve heavy-resin flag
        )
        for c in components
    ]
    asp_pc = Pseudocomponent(
        z=z_asp,
        M=M_asp,
        Tb_K=Tb_K_asp,
        SG=SG_asp,
        xc_lower=1.0 - f_asp,
        xc_upper=1.0,
        is_asphaltene=True,                  # explicit flag; never set by quadrature
    )
    return new_components + [asp_pc]


# ── Recovery-aware assembly (Phase 11) ───────────────────────────────────────

# Tb cap for heavy-resin lump (Decision 31).
_HR_TB_MAX = 1100.0  # K


def append_heavy_resin_and_asphaltene(
    components: list[Pseudocomponent],
    recovery_fraction: float,
    asp_wt_pct: float,
    M_bulk: float,
    SG_bulk: float,
    K_W_bulk: float,
    M_asp: float = ASP_M_DEFAULT,
    Tb_K_asp: float = ASP_TB_K_DEFAULT,
    SG_asp: float = ASP_SG_DEFAULT,
) -> dict:
    """Assemble full-mixture pseudo-components with optional heavy-resin lump.

    When recovery_fraction < 1.0, the distillation curve covers only a
    fraction of the total feed mass.  The unmeasured tail
    (mass = 1 - recovery_fraction - f_asp) is represented by a single
    discrete heavy-resin lump (Decision 30, Phase 11).  Its MW and SG are
    fixed by mass-balance closure on the user-supplied bulk MW and bulk SG.
    When recovery_fraction = 1.0 (or f_hr ≤ 0), dispatches to the existing
    append_asphaltene function — identical to the Phase 8 path.

    Parameters
    ----------
    components : list[Pseudocomponent]
        Distillable pseudo-components from discretize_generalized, with Tb_K
        and SG set (Phase 5 self-consistent solve).  z values must be mole
        fractions within the distillable subfraction (summing to 1).
    recovery_fraction : float
        Mass fraction of total feed covered by the D1160/TBP curve.
        Must be in the open-right interval (0, 1].  Use 1.0 for a complete
        curve (Phase 8 pass-through path).
    asp_wt_pct : float
        Asphaltene content in wt%.
    M_bulk : float
        Bulk molecular weight of the full feed (g/mol).  Becomes a hard
        closure constraint on M_hr when f_hr > 0 (Decision 27 supersession).
    SG_bulk : float
        Bulk specific gravity of the full feed at 15.5 degC.  Becomes a hard
        closure constraint on SG_hr when f_hr > 0.
    K_W_bulk : float
        Bulk Watson K factor.  Used to derive Tb_hr under constant Watson K
        assumption (Decision 31): Tb_hr = (K_W_bulk * SG_hr)^3 / 1.8.
    M_asp : float   Default 1700 g/mol (Gonzalez 2007 nanoaggregate).
    Tb_K_asp : float  Default 1073.15 K (800 deg C convention).
    SG_asp : float    Default 1.15 (Gonzalez 2007).

    Returns
    -------
    dict with keys:
        'components' : list[Pseudocomponent]
            Full component list with true full-mixture mole fractions.
            Order: [distillable_1..N, heavy_resin?, asphaltene].
        'has_heavy_resin' : bool
            True when a heavy-resin lump was created.
        'f_hr' : float
            Heavy-resin mass fraction (0.0 when no lump).
        'M_hr' : float or None
            Heavy-resin lump MW (g/mol), None when no lump.
        'SG_hr' : float or None
            Heavy-resin lump SG, None when no lump.
        'Tb_hr' : float or None
            Heavy-resin lump Tb (K), None when no lump.
        'M_dist_av' : float
            Mole-average MW of distillable subfraction.

    Raises
    ------
    ValueError
        If recovery_fraction + f_asp > 1 (overlap error).
        If M_hr or SG_hr fall outside physically plausible bounds.
        If recovery_fraction is outside (0, 1].

    References
    ----------
    Decision 30/31, Phase 11 (petrochar CURRENT_STATUS.md).
    """
    if not (0.0 < recovery_fraction <= 1.0):
        raise ValueError(
            f"recovery_fraction must be in the open-right interval (0, 1]; "
            f"got {recovery_fraction:.4f}."
        )
    if asp_wt_pct <= 0.0 or asp_wt_pct >= 100.0:
        raise ValueError(
            f"asp_wt_pct must be in the open interval (0, 100); got {asp_wt_pct}."
        )

    # recovery_fraction = 1.0 is the Phase 8 complete-recovery bypass path.
    # The heavy-resin formula f_hr = 1 - f_dist - f_asp is only valid for
    # partial recovery (f_dist < 1 - f_asp); always dispatch to
    # append_asphaltene when the user asserts complete recovery.
    if recovery_fraction >= 1.0 - 1e-9:
        result_comps = append_asphaltene(
            components, asp_wt_pct, M_asp, Tb_K_asp, SG_asp
        )
        M_dist_av = float(sum(c.z * c.M for c in components))
        return {
            'components':      result_comps,
            'has_heavy_resin': False,
            'f_hr':            0.0,
            'M_hr':            None,
            'SG_hr':           None,
            'Tb_hr':           None,
            'M_dist_av':       M_dist_av,
        }

    f_asp  = asp_wt_pct / 100.0
    f_dist = recovery_fraction
    f_hr   = 1.0 - f_dist - f_asp

    # Overlap error: distillation data and ASP fraction together exceed 100%.
    if f_hr < -1e-9:
        raise ValueError(
            f"recovery_fraction ({recovery_fraction:.4f}) + ASP weight fraction "
            f"({f_asp:.4f}) = {recovery_fraction + f_asp:.4f} > 1.0.  "
            f"The distillation data overlaps with the SARA asphaltene class.  "
            f"Reduce recovery_fraction or reduce ASP wt% so they sum to <= 1.0."
        )

    # f_hr ≈ 0: no unmeasured tail — dispatch to Phase 8 path.
    if f_hr < 1e-9:
        result_comps = append_asphaltene(
            components, asp_wt_pct, M_asp, Tb_K_asp, SG_asp
        )
        M_dist_av = float(sum(c.z * c.M for c in components))
        return {
            'components':      result_comps,
            'has_heavy_resin': False,
            'f_hr':            0.0,
            'M_hr':            None,
            'SG_hr':           None,
            'Tb_hr':           None,
            'M_dist_av':       M_dist_av,
        }

    # ── f_hr > 0: compute heavy-resin lump by mass-balance closure ────────────

    # Mole-average MW of the distillable subfraction (z_i sum to 1 within it)
    M_dist_av = float(sum(c.z * c.M for c in components))

    # Volume-additive average SG of the distillable subfraction
    xw_dist = [c.z * c.M / M_dist_av for c in components]
    for i, c in enumerate(components):
        if math.isnan(c.SG):
            raise ValueError(
                f"append_heavy_resin_and_asphaltene: distillable component {i} "
                f"has SG=nan.  Run the self-consistent (Tb, SG) solve (Phase 5) "
                f"before assembling components."
            )
    SG_dist_vavg = 1.0 / sum(xw / c.SG for xw, c in zip(xw_dist, components))

    # ── Bulk MW closure (number-average) ──────────────────────────────────────
    # 1/M_bulk = f_dist/M_dist_av + f_hr/M_hr + f_asp/M_asp
    denom_M = 1.0 / M_bulk - f_dist / M_dist_av - f_asp / M_asp
    if denom_M <= 0.0:
        raise ValueError(
            f"append_heavy_resin_and_asphaltene: bulk MW closure gives "
            f"M_hr <= 0 (denominator = {denom_M:.6f}).  "
            f"Inputs are internally inconsistent.  Checks: "
            f"(1) Is bulk MW = {M_bulk:.1f} g/mol reliable?  "
            f"(2) Is recovery_fraction = {recovery_fraction:.4f} correct?  "
            f"(3) Does the distillation curve reach far enough?  "
            f"Current distillable M_dist_av = {M_dist_av:.1f} g/mol."
        )
    M_hr = f_hr / denom_M

    if not (50.0 < M_hr < 5000.0):
        raise ValueError(
            f"append_heavy_resin_and_asphaltene: bulk MW closure gives "
            f"M_hr = {M_hr:.1f} g/mol outside physically plausible range "
            f"(50, 5000) g/mol.  Inputs are inconsistent.  "
            f"Bulk MW = {M_bulk:.1f}, recovery = {recovery_fraction:.4f}, "
            f"ASP = {asp_wt_pct:.1f} wt%, M_dist_av = {M_dist_av:.1f} g/mol."
        )

    # ── Bulk SG closure (volume-additive) ─────────────────────────────────────
    # 1/SG_bulk = f_dist/SG_dist_vavg + f_hr/SG_hr + f_asp/SG_asp
    denom_SG = 1.0 / SG_bulk - f_dist / SG_dist_vavg - f_asp / SG_asp
    if denom_SG <= 0.0:
        raise ValueError(
            f"append_heavy_resin_and_asphaltene: bulk SG closure gives "
            f"SG_hr <= 0 (denominator = {denom_SG:.6f}).  "
            f"Inputs are inconsistent.  Checks: "
            f"(1) Is bulk SG = {SG_bulk:.4f} reliable?  "
            f"(2) Is the distillable SG distribution correct?  "
            f"Current SG_dist_vavg = {SG_dist_vavg:.4f}."
        )
    SG_hr = f_hr / denom_SG

    if not (0.6 < SG_hr < 1.30):
        raise ValueError(
            f"append_heavy_resin_and_asphaltene: bulk SG closure gives "
            f"SG_hr = {SG_hr:.4f} outside physically plausible range "
            f"(0.6, 1.30).  Inputs are inconsistent.  "
            f"Bulk SG = {SG_bulk:.4f}, recovery = {recovery_fraction:.4f}, "
            f"SG_dist_vavg = {SG_dist_vavg:.4f}, SG_asp = {SG_asp:.4f}."
        )

    # ── Tb_hr from constant Watson K (Decision 31) ────────────────────────────
    # Under SG_i = (1.8 * Tb_i)^(1/3) / K_W_bulk:
    #   Tb_hr = (K_W_bulk * SG_hr)^3 / 1.8
    Tb_hr = (K_W_bulk * SG_hr) ** 3.0 / 1.8
    if Tb_hr > _HR_TB_MAX:
        warnings.warn(
            f"append_heavy_resin_and_asphaltene: Tb_hr = {Tb_hr:.1f} K exceeds "
            f"cap of {_HR_TB_MAX:.0f} K; capping.  "
            f"K_W_bulk = {K_W_bulk:.3f}, SG_hr = {SG_hr:.4f}.  "
            f"This cap is applied for numerical stability; the heavy-resin "
            f"lump does not have a physical boiling point.",
            UserWarning,
            stacklevel=2,
        )
        Tb_hr = _HR_TB_MAX

    # ── Assemble full-mixture mole fractions ──────────────────────────────────
    n_dist  = f_dist / M_dist_av
    n_hr    = f_hr   / M_hr
    n_asp   = f_asp  / M_asp
    n_total = n_dist + n_hr + n_asp

    z_dist_scale = n_dist / n_total
    z_hr_full    = n_hr   / n_total
    z_asp_full   = n_asp  / n_total

    new_components = [
        Pseudocomponent(
            z=c.z * z_dist_scale,
            M=c.M,
            Tb_K=c.Tb_K,
            SG=c.SG,
            xc_lower=c.xc_lower,
            xc_upper=c.xc_upper,
            is_asphaltene=False,
            is_heavy_resin=False,
        )
        for c in components
    ]

    # Heavy-resin lump: xc bounds span [1 - f_asp - f_hr, 1 - f_asp]
    hr_pc = Pseudocomponent(
        z=z_hr_full,
        M=M_hr,
        Tb_K=Tb_hr,
        SG=SG_hr,
        xc_lower=1.0 - f_asp - f_hr,
        xc_upper=1.0 - f_asp,
        is_asphaltene=False,
        is_heavy_resin=True,
    )

    # Asphaltene: xc bounds span [1 - f_asp, 1]
    asp_pc = Pseudocomponent(
        z=z_asp_full,
        M=M_asp,
        Tb_K=Tb_K_asp,
        SG=SG_asp,
        xc_lower=1.0 - f_asp,
        xc_upper=1.0,
        is_asphaltene=True,
        is_heavy_resin=False,
    )

    result_comps = new_components + [hr_pc, asp_pc]

    # Sanity: z sum should equal 1 within numerical precision.  Tolerance 1e-6
    # accommodates floating-point accumulation across 5 GL nodes + HR + ASP
    # after the n_dist / n_total / scale_dist arithmetic above (typically
    # ~1e-7 per node, summing to a few × 1e-7).
    z_sum = sum(c.z for c in result_comps)
    if abs(z_sum - 1.0) > 1e-6:
        raise ValueError(
            f"append_heavy_resin_and_asphaltene: mole fraction sum = {z_sum:.12f} "
            f"deviates from 1.0 by {abs(z_sum - 1.0):.2e}.  Internal error."
        )

    return {
        'components':      result_comps,
        'has_heavy_resin': True,
        'f_hr':            float(f_hr),
        'M_hr':            float(M_hr),
        'SG_hr':           float(SG_hr),
        'Tb_hr':           float(Tb_hr),
        'M_dist_av':       float(M_dist_av),
    }


# ── K_W bin closure check ──────────────────────────────────────────────────────

def kw_bin_check(
    components: list[Pseudocomponent],
    sara_wt_pct: dict,
    kw_sat: float = KW_SAT_DEFAULT,
    kw_aro: float = KW_ARO_DEFAULT,
    flag_tol: float = 3.0,
) -> dict:
    """K_W-bin aggregation of distillable pseudo-components for SARA closure.

    Computes the Watson K factor for each distillable pseudo-component and
    aggregates weight fractions into SAT/ARO/RES bins.  Compares aggregated
    wt% against user-supplied SARA and flags classes where the deviation
    exceeds flag_tol.

    IMPORTANT: this is a closure check, not a retuning step.  Deviations
    indicate inconsistency in the input data and should prompt data inspection.
    Distribution parameters must not be adjusted to close the SARA balance.

    Asphaltene identification
    -------------------------
    Components are classified as asphaltene via the `is_asphaltene` flag,
    which is set exclusively by sara.append_asphaltene().  The Tb > 1000 K
    heuristic is NOT used; using it would silently reclassify quadrature
    nodes that extrapolate beyond the distillation data range.
    All z values are true mole fractions in the full mixture (as produced by
    append_asphaltene).  Weight fractions are derived uniformly as
    wt_i = z_i * M_i / M_mix, where M_mix = sum_j(z_j * M_j) over all
    components.  Distillable components must have non-nan Tb_K and SG.

    Parameters
    ----------
    components : list[Pseudocomponent]
        Pseudo-component list, typically produced by append_asphaltene().
        All z values must be true mole fractions in the full mixture (summing
        to 1), as returned by append_asphaltene().
        Distillable components (is_asphaltene=False) must have Tb_K and SG set.
        Asphaltene components (is_asphaltene=True) are assigned to the ASP class.
    sara_wt_pct : dict
        Keys: 'SAT', 'ARO', 'RES', 'ASP'; values: measured wt%.
    kw_sat : float
        SAT threshold: K_W >= kw_sat → SAT.  Default 12.0 (Riazi p. 75).
    kw_aro : float
        ARO threshold: kw_aro <= K_W < kw_sat → ARO; K_W < kw_aro → RES.
        Default 11.0 (Riazi p. 75).
    flag_tol : float
        Flag any class where |calc - meas| > flag_tol wt%.  Default 3.0 wt%.

    Returns
    -------
    dict with keys:
        'kw_calc'      : dict SAT/ARO/RES/ASP — K_W-binned wt% in total feed.
        'sara_input'   : dict SAT/ARO/RES/ASP — user-supplied wt%.
        'delta_wt_pct' : dict SAT/ARO/RES/ASP — calc minus meas, wt%.
        'flagged'      : bool — True if any |delta| > flag_tol.
        'flag_tol'     : float — threshold used.
        'flags'        : list[str] — message for each flagged class.

    Raises
    ------
    ValueError
        If sara_wt_pct is missing required keys, or if a distillable
        component has nan Tb_K or nan SG (Tb_K and SG must be assigned
        before calling this function — Phase 7 responsibility).

    References
    ----------
    Riazi MNL50 p. 75 (K_W bin thresholds as petroleum classification
    conventions; thresholds are adjustable via kw_sat/kw_aro parameters).
    """
    required = {'SAT', 'ARO', 'RES', 'ASP'}
    missing = required - set(sara_wt_pct.keys())
    if missing:
        raise ValueError(f"sara_wt_pct is missing required keys: {missing}.")

    # Separate ASP, HR, and distillable components.
    # Asphaltene identity set exclusively by append_asphaltene() /
    # append_heavy_resin_and_asphaltene() via is_asphaltene=True.
    # Heavy-resin identity set via is_heavy_resin=True (Phase 11).
    asp_comps  = [c for c in components if c.is_asphaltene]
    hr_comps   = [c for c in components if c.is_heavy_resin]
    distillable = [c for c in components
                   if not c.is_asphaltene and not c.is_heavy_resin]

    # Validate non-asphaltene components have Tb_K and SG populated
    for i, c in enumerate(distillable + hr_comps):
        if math.isnan(c.Tb_K) or math.isnan(c.SG):
            raise ValueError(
                f"kw_bin_check: component {i} has nan Tb_K or SG.  "
                f"Assign Tb_K and SG values (Phase 7) before calling kw_bin_check."
            )

    # Mixture MW — z values are true mole fractions in the full mixture
    M_mix = float(sum(c.z * c.M for c in components))

    # Weight fraction of each component: wt_i = z_i * M_i / M_mix
    def _wt(c: Pseudocomponent) -> float:
        return c.z * c.M / M_mix

    # For K_W-bin aggregation, HR components are treated as part of the
    # binnable mass (they have Tb_K and SG set, so watson_k works).
    # Their K_W = K_W_bulk by construction → they fall in whatever bin
    # K_W_bulk occupies.  HR mass is also reported separately in 'HR' key.
    binnable = distillable + hr_comps  # all non-ASP components

    if len(binnable) == 0:
        kw_calc = {
            'SAT': 0.0,
            'ARO': 0.0,
            'RES': 0.0,
            'ASP': float(sum(_wt(c) for c in asp_comps)) * 100.0,
            'HR':  0.0,
        }
    else:
        kw_arr = np.array([watson_k(c.Tb_K, c.SG) for c in binnable])
        wt_arr = np.array([_wt(c) for c in binnable])

        # Aggregate total-feed wt% per K_W bin (includes HR contribution)
        kw_calc = {
            'SAT': float(np.sum(wt_arr[kw_arr >= kw_sat]))                        * 100.0,
            'ARO': float(np.sum(wt_arr[(kw_arr >= kw_aro) & (kw_arr < kw_sat)]))  * 100.0,
            'RES': float(np.sum(wt_arr[kw_arr < kw_aro]))                          * 100.0,
            'ASP': float(sum(_wt(c) for c in asp_comps))                           * 100.0,
            'HR':  float(sum(_wt(c) for c in hr_comps))                            * 100.0,
        }

    sara_in = {k: float(sara_wt_pct[k]) for k in ('SAT', 'ARO', 'RES', 'ASP')}
    # SARA-closure delta uses only SAT/ARO/RES/ASP (HR is diagnostic only).
    delta   = {k: kw_calc[k] - sara_in[k] for k in ('SAT', 'ARO', 'RES', 'ASP')}

    flags = []
    for k in ('SAT', 'ARO', 'RES', 'ASP'):
        if abs(delta[k]) > flag_tol:
            flags.append(
                f"{k}: K_W-binned = {kw_calc[k]:.2f} wt%, "
                f"measured = {sara_in[k]:.2f} wt%, "
                f"delta = {delta[k]:+.2f} wt% (exceeds +/-{flag_tol:.1f} wt%)"
            )

    return {
        'kw_calc':      kw_calc,
        'sara_input':   sara_in,
        'delta_wt_pct': delta,
        'flagged':      bool(flags),
        'flag_tol':     float(flag_tol),
        'flags':        flags,
    }
