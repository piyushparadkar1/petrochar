"""
Specific gravity distribution for petroleum fraction pseudo-components.

Two methods are supported, user-toggleable:

1. **Constant Watson K** (default, recommended for VTB-scale fractions)
   Watson K factor K_W is computed from bulk M and bulk SG via Eq. 2.13.
   All pseudo-components are assumed to share the same K_W.  Per-component
   SG is then recovered algebraically from each component's Tb:

       SG_i = (1.8 * Tb_i)^(1/3) / K_W        ...inversion of Eq. 2.13

   This is a one-liner rearrangement of Eq. 2.13; no root-finding required.

2. **Generalized distribution fit**
   Fits Riazi Eq. 4.56 to (x_cv, SG) data using GeneralizedDistribution
   from core/distribution.py.  Cumulative volume fraction x_cv is the
   recommended basis for SG per Riazi MNL50 Table 4.13; weight basis is
   also supported by passing the appropriate xc array.

   Default two-parameter fit uses B_SG = 3 per Riazi (2005) p. 174.
   Three-parameter fit is available for pass-gate validation.

References
----------
Riazi (2005)    : §4.5.4.1 Eq. 4.56, Table 4.13; §2.4.2.1 Eq. 2.13.
Whitson (1983)  : SPE Journal 23(4), 683-694 (constant K_W method, p. 178).
"""

import numpy as np

from core.correlations import riazi_daubert_Tb, watson_k
from core.distribution import GeneralizedDistribution


# ── Method 1: Constant Watson K ───────────────────────────────────────────────

def sg_from_watson_k(
    Tb_K,
    M_bulk: float,
    SG_bulk: float,
) -> tuple[np.ndarray, float]:
    """Compute per-component SG assuming a constant Watson K factor.

    The bulk Watson K factor K_W is derived from bulk M and SG (via a
    riazi_daubert_Tb forward call to get a representative bulk Tb, then
    Eq. 2.13).  Each component's SG follows algebraically from its own Tb:

        SG_i = (1.8 * Tb_i [K])^(1/3) / K_W

    This is the direct algebraic inversion of K_W = (1.8*Tb)^(1/3) / SG
    and requires no root-finding.

    Parameters
    ----------
    Tb_K : array-like
        Normal boiling points in K for each pseudo-component.
    M_bulk : float
        Bulk molecular weight of the C7+ fraction, g/mol.
    SG_bulk : float
        Bulk specific gravity of the C7+ fraction at 15.5 °C.

    Returns
    -------
    SG : ndarray
        Per-component specific gravity array, same length as Tb_K.
    K_W : float
        Watson K factor used (derived from bulk properties).

    References
    ----------
    Riazi MNL50 Eq. 2.13 (Watson K); §4.5.4.4 constant-K_W SG method;
    Whitson (1983) SPE Journal 23(4), 683-694 (p. 178 in Riazi MNL50).
    """
    Tb_K = np.asarray(Tb_K, dtype=float)

    if M_bulk <= 0 or SG_bulk <= 0:
        raise ValueError(
            f"M_bulk and SG_bulk must be positive; "
            f"got M_bulk={M_bulk}, SG_bulk={SG_bulk}."
        )
    if np.any(Tb_K <= 0):
        raise ValueError("All Tb_K values must be positive.")

    Tb_bulk = riazi_daubert_Tb(M_bulk, SG_bulk)
    K_W_val = watson_k(Tb_bulk, SG_bulk)
    SG = (1.8 * Tb_K) ** (1.0 / 3.0) / K_W_val
    return SG, float(K_W_val)


# ── Method 2: Generalized distribution fit ────────────────────────────────────

class SGDistributionFit:
    """Specific gravity generalized distribution (Riazi Eq. 4.56).

    Thin wrapper around GeneralizedDistribution that enforces the Riazi
    recommendation for SG: two-parameter fit with B_SG = 3 as default
    (Riazi 2005, p. 174).  Three-parameter fit is available via mode='3param'.

    Cumulative volume fraction x_cv is the canonical basis for SG per
    Riazi MNL50 Table 4.13.  Cumulative weight fraction x_cw is also
    accepted (pass xc_weight as the xc argument).

    Attributes (populated after fit)
    ----------------------------------
    SG_o : float   Onset SG parameter (= GeneralizedDistribution.P_o).
    A    : float   Scale parameter.
    B    : float   Shape parameter.
    """

    def __init__(self) -> None:
        self._dist: GeneralizedDistribution | None = None
        self.SG_o: float | None = None
        self.A:    float | None = None
        self.B:    float | None = None

    def fit(
        self,
        xc,
        SG_data,
        mode: str = '2param',
        B_fixed: float = 3.0,
    ) -> 'SGDistributionFit':
        """Fit the SG generalized distribution to (x_c, SG) data.

        Parameters
        ----------
        xc : array-like
            Cumulative fraction values in (0, 1).  Use x_cv (volume basis)
            to reproduce Riazi Table 4.13 volume rows; x_cw for weight rows.
        SG_data : array-like
            Specific gravity at each cumulative fraction.
        mode : str
            '2param' — fix B = B_fixed (default, recommended).
            '3param' — fit B freely along with SG_o and A.
        B_fixed : float
            B value for two-parameter fit.  Default 3.0 (Riazi p. 174).

        Returns
        -------
        SGDistributionFit
            self.

        References
        ----------
        Riazi MNL50 §4.5.4.1 Eq. 4.56; Table 4.13 SG rows.
        Default B_SG = 3 per Riazi (2005) p. 174.
        """
        self._dist = GeneralizedDistribution().fit(
            xc, SG_data, mode=mode, B_fixed=B_fixed
        )
        self.SG_o = self._dist.P_o
        self.A    = self._dist.A
        self.B    = self._dist.B
        return self

    def SG(self, xc):
        """SG at cumulative fraction xc (forward CDF evaluation).

        Parameters
        ----------
        xc : float or array-like
            Cumulative fraction in (0, 1).

        Returns
        -------
        float or ndarray

        References
        ----------
        Riazi MNL50 Eq. 4.56.
        """
        self._require_fit()
        return self._dist.P(xc)

    def xc(self, SG_val):
        """Cumulative fraction at SG value (inverse CDF).

        Parameters
        ----------
        SG_val : float or array-like

        Returns
        -------
        float or ndarray
        """
        self._require_fit()
        return self._dist.x_c(SG_val)

    def pdf(self, SG_val):
        """Probability density f(SG) = dx_c / d(SG).

        Parameters
        ----------
        SG_val : float or array-like

        Returns
        -------
        float or ndarray
        """
        self._require_fit()
        return self._dist.pdf(SG_val)

    def average_SG(self) -> float:
        """Analytical mean SG from the fitted distribution.

        <SG> = SG_o * (1 + (A/B)^(1/B) * Gamma(1 + 1/B))

        For volume-basis fits this is the volume-weighted average SG.
        The density-additive bulk SG (1/sum(x_w/SG_i)) is computed
        separately by bulk_SG_av() in core/mw_distribution.py.

        Returns
        -------
        float
        """
        self._require_fit()
        return self._dist.average()

    @property
    def params(self) -> dict:
        """Fitted parameters: {'P_o' (= SG_o), 'A', 'B', 'mode'}."""
        self._require_fit()
        return self._dist.params

    @property
    def fit_quality(self) -> dict:
        """Fit quality metrics (RMS, pct_AAD, R_squared) in SG units."""
        self._require_fit()
        return self._dist.fit_quality

    def _require_fit(self) -> None:
        if self._dist is None:
            raise RuntimeError(
                "SGDistributionFit has not been fitted. Call fit() first."
            )
