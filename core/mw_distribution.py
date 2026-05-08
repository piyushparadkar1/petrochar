"""
Molecular weight computation for pseudo-components and bulk closure verification.

Per-component M_i is derived from (Tb_i, SG_i) via the Riazi-Daubert inverse
correlation (riazi_daubert_M from Phase 1).  Bulk-average properties (M_av,
SG_av) are computed from weight fractions using mixture-rule formulas
consistent with Riazi MNL50 §4.5.4.4 (Eqs. 4.74-4.76):

    M_av  = 1 / sum(x_w_i / M_i)          ...Eq. 4.74 (number-average MW)
    SG_av = 1 / sum(x_w_i / SG_i)         ...Eq. 4.76 (density-additive mixing)

These bulk averages serve as closure checks — not constraints — since M_i and
SG_i are derived independently from Tb_i.

References
----------
Riazi (2005) : §4.5.4.4 Eqs. 4.74-4.76 (page 174-175).
"""

import warnings

import numpy as np

from core.correlations import riazi_daubert_M


def compute_M_array(Tb_K, SG) -> np.ndarray:
    """Compute molecular weight for each pseudo-component from Tb and SG.

    Calls riazi_daubert_M (which inverts Riazi-Daubert Eq. 2.56/2.57) for
    each (Tb_i, SG_i) pair.  Operates on 1-D arrays via a Python loop;
    suitable for the typical 5-50 pseudo-component count.

    Parameters
    ----------
    Tb_K : array-like
        Normal boiling points in K.  Must not contain NaN; exclude any
        plus-fraction row that has no measured Tb before calling.
    SG : array-like
        Specific gravity at 15.5 °C, corresponding to each Tb.  Must be
        the same length as Tb_K.

    Returns
    -------
    M : ndarray
        Molecular weights in g/mol, same length as Tb_K.

    References
    ----------
    Riazi MNL50 Eqs. 2.56, 2.57 (inverted numerically by riazi_daubert_M).
    """
    Tb_K = np.asarray(Tb_K, dtype=float)
    SG   = np.asarray(SG,   dtype=float)

    if Tb_K.ndim != 1 or Tb_K.shape != SG.shape:
        raise ValueError(
            f"Tb_K and SG must be 1-D arrays of equal length; "
            f"got shapes {Tb_K.shape} and {SG.shape}."
        )
    if np.any(np.isnan(Tb_K)):
        raise ValueError(
            "Tb_K contains NaN.  Exclude plus-fraction rows with no measured "
            "Tb before calling compute_M_array; use the tabulated M for those rows."
        )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)   # suppress extrapolation warns
        M = np.array([riazi_daubert_M(float(t), float(s)) for t, s in zip(Tb_K, SG)])
    return M


def bulk_M_av(x_w, M) -> float:
    """Number-average molecular weight from weight fractions and M_i.

    M_av = 1 / sum(x_w_i / M_i)

    This is the standard formula for the number-average MW of a mixture
    when weight fractions x_w_i and component molecular weights M_i are
    known.  It is equivalent to total mass divided by total moles.

    Parameters
    ----------
    x_w : array-like
        Weight fractions.  Need not sum exactly to 1 (minor rounding in
        tabulated data is tolerated); a warning is issued if |sum-1| > 0.01.
    M : array-like
        Molecular weights in g/mol, same length as x_w.

    Returns
    -------
    float
        Bulk number-average molecular weight, g/mol.

    References
    ----------
    Riazi MNL50 Eq. 4.74 (generalized distribution; reduces to this discrete
    form for a finite set of pseudo-components).
    """
    x_w = np.asarray(x_w, dtype=float)
    M   = np.asarray(M,   dtype=float)

    if x_w.shape != M.shape:
        raise ValueError(
            f"x_w and M must have the same shape; got {x_w.shape} vs {M.shape}."
        )
    if np.any(M <= 0):
        raise ValueError("All M values must be positive.")

    total = float(x_w.sum())
    if abs(total - 1.0) > 0.01:
        warnings.warn(
            f"bulk_M_av: weight fractions sum to {total:.4f}, not 1.0; "
            f"results may reflect rounding in tabulated data.",
            UserWarning,
            stacklevel=2,
        )

    return float(1.0 / np.sum(x_w / M))


def bulk_SG_av(x_w, SG) -> float:
    """Density-additive (weight-fraction-harmonic) bulk specific gravity.

    SG_av = 1 / sum(x_w_i / SG_i)

    Derived from ideal volume additivity: V_total = sum(m_i / rho_i).
    Equivalent to total mass divided by total volume.  This is the formula
    to use when checking against a measured bulk SG because laboratory
    density measurements are volume-based.

    Parameters
    ----------
    x_w : array-like
        Weight fractions.
    SG : array-like
        Specific gravity at 15.5 °C, same length as x_w.

    Returns
    -------
    float
        Bulk specific gravity (dimensionless).

    References
    ----------
    Riazi MNL50 Eq. 4.76 (density-additive mixing rule in the discrete-
    component limit of the generalized distribution framework).
    """
    x_w = np.asarray(x_w, dtype=float)
    SG  = np.asarray(SG,  dtype=float)

    if x_w.shape != SG.shape:
        raise ValueError(
            f"x_w and SG must have the same shape; got {x_w.shape} vs {SG.shape}."
        )
    if np.any(SG <= 0):
        raise ValueError("All SG values must be positive.")

    return float(1.0 / np.sum(x_w / SG))
