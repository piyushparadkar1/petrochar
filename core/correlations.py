"""
Core physical property correlations for petroleum fraction characterization.

All functions operate in SI-adjacent internal units:
  - Tb in Kelvin
  - SG dimensionless at 15.5 °C (60 °F)
  - M in g/mol

Unit conversion (to Rankine) happens only inside watson_k; all other
functions receive and return values in the units stated above.

References
----------
Riazi (2005)    : Riazi, M.R. Characterization and Properties of Petroleum
                  Fractions. ASTM Manual MNL50. ISBN 0-8031-3361-8.
Whitson (1983)  : Whitson, C.H. Characterizing hydrocarbon plus fractions.
                  SPE Journal 23(4), 683-694.
Watson (1933)   : Watson, K.M. and Nelson, E.F. Ind. Eng. Chem. 25, 880.
"""

import warnings
import math

import numpy as np
from scipy import optimize, special


# ── Riazi-Daubert Tb correlations ────────────────────────────────────────────

def riazi_daubert_Tb(M: float, SG: float) -> float:
    """Compute normal boiling point from molecular weight and specific gravity.

    Two-regime correlation per Riazi MNL50 §2.4.2.1:
        Eq. 2.56  valid for M in [70, 300]   (light/medium hydrocarbons)
        Eq. 2.57  valid for M in [300, 700]  (heavy hydrocarbons)

    For M > 700 Eq. 2.57 is extrapolated (no better correlation in MNL50);
    a UserWarning is raised.  For M < 70 a ValueError is raised.

    Parameters
    ----------
    M : float
        Molecular weight, g/mol.
    SG : float
        Specific gravity at 15.5 °C, dimensionless.

    Returns
    -------
    float
        Normal boiling point, K.

    References
    ----------
    Riazi MNL50 Eqs. 2.56 and 2.57, page 58 (PDF page 77).
    """
    if M < 70:
        raise ValueError(
            f"M = {M:.1f} g/mol is below the Riazi-Daubert validity range "
            f"(M >= 70 required). Outside refinery-relevant range."
        )
    if M <= 300:
        # Eq. 2.56 — valid M = 70-300
        Tb = (
            3.76587
            * np.exp(3.7741e-3 * M + 2.98404 * SG - 4.25288e-3 * M * SG)
            * M ** 0.40167
            * SG ** (-1.58262)
        )
    else:
        # Eq. 2.57 — valid M = 300-700, recommended for M > 300
        if M > 700:
            warnings.warn(
                f"M = {M:.1f} g/mol exceeds Riazi-Daubert Eq. 2.57 validated "
                f"range (300-700); result is an extrapolation.",
                UserWarning,
                stacklevel=2,
            )
        Tb = (
            9.3369
            * np.exp(1.6514e-4 * M + 1.4103 * SG - 7.5152e-4 * M * SG)
            * M ** 0.5369
            * SG ** (-0.7276)
        )
    return float(Tb)


def riazi_daubert_M(Tb_K: float, SG: float) -> float:
    """Compute molecular weight from normal boiling point and specific gravity.

    Numerically inverts riazi_daubert_Tb via Brent's method.  The root is
    sought separately in each Eq. 2.56 regime (M in [70, 300]) and Eq. 2.57
    regime (M in [300, 2000]).  The candidate whose forward-Tb matches the
    input Tb within tolerance is returned; if both brackets find roots the one
    with the smaller residual is preferred.

    Parameters
    ----------
    Tb_K : float
        Normal boiling point, K.
    SG : float
        Specific gravity at 15.5 °C, dimensionless.

    Returns
    -------
    float
        Molecular weight, g/mol.

    References
    ----------
    Inverts Riazi MNL50 Eqs. 2.56 and 2.57, page 58 (PDF page 77).
    """
    def residual(M_val: float) -> float:
        return riazi_daubert_Tb(M_val, SG) - Tb_K

    best_M = None
    best_res = np.inf

    # Regime 1 — Eq. 2.56, M in [70, 300]
    try:
        f_lo = residual(70.0)
        f_hi = residual(300.0)
        if f_lo * f_hi <= 0:
            M_56 = optimize.brentq(residual, 70.0, 300.0, xtol=0.01, rtol=1e-6)
            res_56 = abs(residual(M_56))
            if res_56 < best_res:
                best_M, best_res = M_56, res_56
    except Exception:
        pass

    # Regime 2 — Eq. 2.57, M in (300, 2000]
    try:
        f_lo = residual(300.01)
        f_hi = residual(2000.0)
        if f_lo * f_hi <= 0:
            M_57 = optimize.brentq(residual, 300.01, 2000.0, xtol=0.01, rtol=1e-6)
            res_57 = abs(residual(M_57))
            if res_57 < best_res:
                best_M, best_res = M_57, res_57
    except Exception:
        pass

    if best_M is None:
        raise ValueError(
            f"riazi_daubert_M: no root found for Tb_K={Tb_K:.2f} K, "
            f"SG={SG:.4f}. Check that inputs are in the valid range."
        )
    return float(best_M)


# Note: Riazi-Daubert §2.4.3.1 Eqs. 2.59 (Tb, I) and 2.60 (M, I) for direct SG
# computation require refractivity index I, which is not a standard refinery
# input. We numerically invert riazi_daubert_Tb instead. Verified 2026-05-06
# against PDF pages 077-078; no I-free SG correlation exists in §2.4.3.1.
def riazi_daubert_SG(M: float, Tb_K: float) -> float:
    """Compute specific gravity from molecular weight and normal boiling point.

    Numerically inverts riazi_daubert_Tb(M, SG) = Tb_K for SG using Brent's
    method.  This is the correct approach for petrochar because the explicit
    Riazi MNL50 inverse correlations (Eqs. 2.59/2.60) require the refractivity
    index I, which is not a standard refinery measurement.

    The Riazi-Daubert form θ = a·exp(…)·M^e·SG^f has a minimum in Tb vs SG
    at SG_min = −f/(c + d·M).  Brent's method is applied on [SG_min, 1.30]
    (the monotone-increasing branch), which covers all refinery-relevant SG
    values (SG_min ≈ 0.61–0.74 for typical M).

    Parameters
    ----------
    M : float
        Molecular weight, g/mol.
    Tb_K : float
        Normal boiling point, K.

    Returns
    -------
    float
        Specific gravity at 15.5 °C, dimensionless.

    References
    ----------
    Inverts Riazi MNL50 Eqs. 2.56 and 2.57, page 58 (PDF page 77).
    Explicit Eqs. 2.59 (SG from Tb, I) and 2.60 (SG from M, I) require I.
    """
    def residual(SG_val: float) -> float:
        return riazi_daubert_Tb(M, SG_val) - Tb_K

    # Analytically-computed SG where d(Tb)/d(SG) = 0 (minimum of Tb vs SG).
    # Form: ln Tb = … + c·SG + d·M·SG + f·ln(SG)  →  d(ln Tb)/d(SG) = c + d·M + f/SG = 0
    if M <= 300:
        c, d, f = 2.98404, -4.25288e-3, -1.58262   # Eq. 2.56 coefficients
    else:
        c, d, f = 1.4103, -7.5152e-4, -0.7276      # Eq. 2.57 coefficients
    denom = c + d * M
    SG_min = (-f / denom) if denom > 0 else 0.40
    SG_lo = max(SG_min, 0.40)   # lower bracket on monotone-increasing branch
    SG_hi = 1.30

    try:
        SG_out = optimize.brentq(residual, SG_lo, SG_hi, xtol=1e-6, rtol=1e-8)
    except ValueError as exc:
        raise ValueError(
            f"riazi_daubert_SG: no solution in SG=[{SG_lo:.4f}, {SG_hi}] for "
            f"M={M:.1f} g/mol, Tb_K={Tb_K:.2f} K. "
            f"Check that inputs are in the Riazi-Daubert valid range."
        ) from exc
    return float(SG_out)


# ── Watson K factor ───────────────────────────────────────────────────────────

def watson_k(Tb_K: float, SG: float) -> float:
    """Compute the Watson (UOP) characterisation factor K_W.

    Defined as K_W = Tb_R^(1/3) / SG where Tb_R is the normal boiling point
    in Rankine.  The factor 1.8 converts Kelvin to Rankine.

    Parameters
    ----------
    Tb_K : float
        Normal boiling point, K.
    SG : float
        Specific gravity at 15.5 °C, dimensionless.

    Returns
    -------
    float
        Watson characterisation factor, dimensionless.

    References
    ----------
    Watson & Nelson (1933) Ind. Eng. Chem. 25, 880.
    Riazi MNL50 Eq. 2.13, page 35 (PDF page 54).
    """
    Tb_R = 1.8 * Tb_K  # convert K → Rankine; 1.8 is exact by definition
    return float(Tb_R ** (1.0 / 3.0) / SG)


# ── Aromaticity factor ────────────────────────────────────────────────────────

def aromaticity_gamma(K_W: float) -> float:
    """Compute the aromaticity factor γ from the Watson K factor.

    Linear clamp between the paraffinic anchor (K_W = 13.0 → γ = 0) and the
    aromatic anchor (K_W = 9.5 → γ = 1).  Values outside [9.5, 13.0] are
    clamped to [0, 1].

    Parameters
    ----------
    K_W : float
        Watson characterisation factor, dimensionless.

    Returns
    -------
    float
        Aromaticity factor γ in [0, 1]; 0 = paraffinic, 1 = aromatic.

    References
    ----------
    Architecture commitment in CLAUDE_CODE_PROMPT.md §ARCHITECTURE COMMITMENTS.
    γ = clamp((13.0 − K_W) / (13.0 − 9.5), 0, 1)
    """
    gamma = (13.0 - K_W) / (13.0 - 9.5)
    return float(np.clip(gamma, 0.0, 1.0))


# ── Gamma function ────────────────────────────────────────────────────────────

def gamma_function(x: float) -> float:
    """Evaluate the gamma function Γ(x).

    Uses scipy.special.gamma, which implements the recurrence relation
    Γ(x + 1) = x · Γ(x) to extend the polynomial approximation of Eq. 4.44
    (Whitson 1983) to all positive real x.  The polynomial approximation alone
    is valid only for x in (0, 1]; scipy reproduces it to machine precision.

    Parameters
    ----------
    x : float
        Argument.  Must satisfy x > 0; x ≤ 0 raises ValueError because
        negative-integer poles and x = 0 are outside petrochar's domain.

    Returns
    -------
    float
        Γ(x).

    References
    ----------
    Whitson (1983) SPE Journal 23(4), 683-694, Eq. 4.44.
    Riazi MNL50 Eq. 4.44 (equivalent polynomial approximation for x ∈ (0,1]).
    """
    if x <= 0:
        raise ValueError(
            f"gamma_function: x = {x} is not in the valid domain (x > 0 required). "
            f"Negative and zero arguments are outside petrochar's use range."
        )
    return float(special.gamma(x))
