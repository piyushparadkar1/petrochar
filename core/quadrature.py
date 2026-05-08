"""
Gaussian quadrature discretization of a generalized distribution into
N pseudo-components (Riazi MNL50 §4.6.1.1, Eqs. 4.83-4.91).

The Gauss-Laguerre quadrature rule approximates integrals of the form

    int_0^inf  f(t) * e^{-t}  dt  ~  sum_i  w_i * f(y_i)

Applied to the generalized distribution with the substitution
t = (B/A) * P*^B, where P* = (P - P_o) / P_o, the property at each
quadrature node is the analytically correct sub-fraction average (Eq. 4.84):

    P_i = P_o * (1 + (A/B * y_i)^{1/B})

The Gauss-Laguerre weights w_i are the mole fractions z_i (Eq. 4.83).
Cumulative weight-fraction bounds follow from z_i and M_i:

    x_wt_i  = z_i * M_i / sum_j(z_j * M_j)
    xc_lower_i = sum_{j<i} x_wt_j
    xc_upper_i = xc_lower_i + x_wt_i

References
----------
Riazi (2005) : §4.6.1.1, Table 4.21 (quadrature roots/weights),
               Table 4.22 (3-point validation, Example 4.14).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from core.distribution import GeneralizedDistribution


# ── Quadrature tables (Riazi Table 4.21) ─────────────────────────────────────

_ROOTS_3   = np.array([0.41577,  2.29428,  6.28995])
_WEIGHTS_3 = np.array([0.711093, 0.278518, 0.0103893])

_ROOTS_5   = np.array([0.26356,   1.41340,    3.59643,     7.08581,      12.64080])
_WEIGHTS_5 = np.array([0.521756,  0.398667,   0.0759424,   0.00361176,   0.0000233700])

_TABLES: dict[int, tuple[np.ndarray, np.ndarray]] = {
    3: (_ROOTS_3,   _WEIGHTS_3),
    5: (_ROOTS_5,   _WEIGHTS_5),
}


def quadrature_points(N: int) -> tuple[np.ndarray, np.ndarray]:
    """Return Gauss-Laguerre roots y_i and weights w_i for the N-point rule.

    Parameters
    ----------
    N : int
        Number of quadrature points.  Supported values: 3, 5.

    Returns
    -------
    y : ndarray, shape (N,)
        Gauss-Laguerre roots (dimensionless), in ascending order.
    w : ndarray, shape (N,)
        Gauss-Laguerre weights; sum to 1.

    References
    ----------
    Riazi MNL50 Table 4.21 (page 184).
    """
    if N not in _TABLES:
        raise ValueError(
            f"N must be 3 or 5; got {N}.  "
            f"Only Riazi Table 4.21 values are supported."
        )
    y, w = _TABLES[N]
    return y.copy(), w.copy()


# ── Pseudo-component record ───────────────────────────────────────────────────

@dataclass
class Pseudocomponent:
    """One pseudo-component produced by Gaussian quadrature discretization.

    Fields
    ------
    z        : mole fraction (= Gauss-Laguerre weight w_i).
    M        : molecular weight, g/mol (sub-fraction average via Eq. 4.84).
    Tb_K     : normal boiling point, K.  nan if not yet assigned.
    SG       : specific gravity at 15.5 degC.  nan if not yet assigned.
    xc_lower : lower bound of cumulative weight-fraction interval [0, 1).
    xc_upper : upper bound of cumulative weight-fraction interval (0, 1].
    """
    z:        float
    M:        float
    Tb_K:     float = float('nan')
    SG:       float = float('nan')
    xc_lower: float = 0.0
    xc_upper: float = 0.0


# ── Discretization ────────────────────────────────────────────────────────────

def discretize_generalized(
    N: int,
    distribution: GeneralizedDistribution,
) -> list[Pseudocomponent]:
    """Discretize a GeneralizedDistribution into N pseudo-components.

    Uses the N-point Gauss-Laguerre rule (Riazi §4.6.1.1) to replace the
    continuous M distribution with a finite set of fractions, each carrying
    a mole fraction z_i and a mean molecular weight M_i.

    Mole fractions equal the Gauss-Laguerre weights (Eq. 4.83):

        z_i = w_i

    Property at each node is the analytically correct sub-fraction average
    in the Gauss-Laguerre framework, NOT a sample value at that point
    (Eq. 4.84):

        M_i = P_o * (1 + (A/B * y_i)^{1/B})

    Cumulative weight-fraction bounds are ordered lightest to heaviest
    (y_i are already ascending in Table 4.21).

    Parameters
    ----------
    N : int
        Number of quadrature points.  Supported: 3, 5.
    distribution : GeneralizedDistribution
        Fitted (or manually parameterized via from_params) distribution for
        molecular weight.  P_o, A, B must be set.

    Returns
    -------
    list[Pseudocomponent]
        N records, ordered lightest to heaviest.  Tb_K and SG fields are nan;
        assign separately from Tb and SG distributions in later workflow steps.

    References
    ----------
    Riazi MNL50 §4.6.1.1, Eqs. 4.83-4.84; Table 4.22 (3-pt Example 4.14).
    """
    y, w = quadrature_points(N)
    P_o, A, B = distribution.P_o, distribution.A, distribution.B

    # Eq. 4.84: sub-fraction average property at each Gauss-Laguerre node
    M_i = P_o * (1.0 + (A / B * y) ** (1.0 / B))

    # Mole fractions = Gauss-Laguerre weights (Eq. 4.83)
    z_i = w.copy()

    # Weight fractions for xc bounds
    zM    = z_i * M_i
    x_wt  = zM / zM.sum()

    xc_lower = np.concatenate([[0.0], np.cumsum(x_wt)[:-1]])
    xc_upper = np.cumsum(x_wt)

    return [
        Pseudocomponent(
            z=float(z_i[i]),
            M=float(M_i[i]),
            xc_lower=float(xc_lower[i]),
            xc_upper=float(xc_upper[i]),
        )
        for i in range(N)
    ]
