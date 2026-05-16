"""
SARA closure check and discrete asphaltene component assembly.

Two responsibilities:
1. Append the asphaltene pseudo-component (discrete, literature-defined).
2. K_W-bin aggregation of distillable pseudo-components for SARA closure check.

ARCHITECTURAL COMMITMENT (frozen):
    SAT/ARO/RES wt% from SARA are closure checks ONLY.  They are never used
    to adjust pseudo-component properties.  If K_W-binned wt% deviate from
    measured SARA, the deviation is flagged for data inspection.
    Distribution parameters are NEVER retuned to close the SARA balance.

ASP wt% is the only SARA input that directly affects pseudo-component
construction: it sets the mass fraction of the discrete asphaltene component.

Mole fraction convention
------------------------
After append_asphaltene() is called, ALL Pseudocomponent.z values are true
mole fractions in the full mixture (distillable + asphaltene), summing to 1.

append_asphaltene() performs the basis conversion internally.  The input
distillable components carry z_i = Gauss-Laguerre weight (mole fraction
within distillable subfraction).  Given ASP weight fraction f_asp:

    M_dist_av = sum_i(z_i * M_i)          (mole-avg MW of distillable)
    n_dist = (1 - f_asp) / M_dist_av      (moles of distillable per kg feed)
    n_asp  = f_asp / M_asp                (moles of ASP per kg feed)
    n_tot  = n_dist + n_asp

    z_i_true = z_i * n_dist / n_tot       (rescaled distillable)
    z_asp    = n_asp / n_tot              (ASP mole fraction)

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
            is_asphaltene=c.is_asphaltene,   # preserve any pre-existing flag
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

    # Separate ASP from distillable using the explicit is_asphaltene flag.
    # The Tb > 1000 K heuristic is NOT used here; asphaltene identity is set
    # exclusively by sara.append_asphaltene() via is_asphaltene=True.
    asp_comps   = [c for c in components if c.is_asphaltene]
    distillable = [c for c in components if not c.is_asphaltene]

    # Validate distillable components have Tb_K and SG populated
    for i, c in enumerate(distillable):
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

    if len(distillable) == 0:
        kw_calc = {
            'SAT': 0.0,
            'ARO': 0.0,
            'RES': 0.0,
            'ASP': float(sum(_wt(c) for c in asp_comps)) * 100.0,
        }
    else:
        # Watson K for each distillable component
        kw_arr = np.array([watson_k(c.Tb_K, c.SG) for c in distillable])
        wt_arr = np.array([_wt(c) for c in distillable])

        # Aggregate total-feed wt% per K_W bin
        kw_calc = {
            'SAT': float(np.sum(wt_arr[kw_arr >= kw_sat]))                       * 100.0,
            'ARO': float(np.sum(wt_arr[(kw_arr >= kw_aro) & (kw_arr < kw_sat)])) * 100.0,
            'RES': float(np.sum(wt_arr[kw_arr < kw_aro]))                         * 100.0,
            'ASP': float(sum(_wt(c) for c in asp_comps))                          * 100.0,
        }

    sara_in = {k: float(sara_wt_pct[k]) for k in ('SAT', 'ARO', 'RES', 'ASP')}
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
