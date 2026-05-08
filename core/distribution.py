"""
Riazi generalized distribution model (Eq. 4.56) for petroleum fraction properties.

The model relates a property P (Tb in K, M in g/mol, or SG) to cumulative
fraction x_c via:

    P = P_o * (1 + P*)
    P* = [(A/B) * ln(1/(1-x_c))]^(1/B)      ...Eq. 4.56

Fitting uses the linearized form (Eq. 4.57):

    Y = B * X + c
    Y = ln(-ln(1 - x_c))
    X = ln(P*) = ln((P - P_o) / P_o)
    c = ln(B/A)  →  A = B * exp(-c)

scipy.optimize.least_squares minimizes sum(r_i^2) where r_i = Y_i - B*X_i - c.
P_o is optimized jointly with B and c; all three parameters (P_o, B, c) are free
in the three-parameter fit.  For the two-parameter fit, B is fixed externally and
only P_o and c are free.

References
----------
Riazi (2005)   : §4.5.4.1, Eqs. 4.56–4.57, Table 4.13 (page 172–175).
Whitson (1983) : SPE Journal 23(4), 683-694 (gamma model; Eq. 4.56 is the
                 generalized successor used here per Riazi p. 178-179).
"""

import numpy as np
from scipy import optimize, special


class GeneralizedDistribution:
    """Riazi generalized distribution (Eq. 4.56) for one petroleum property.

    Fit a CDF P(x_c) to discrete (x_c, P) data and expose evaluation,
    inversion, PDF, and analytical average.

    Attributes (populated after fit)
    ---------------------------------
    P_o : float  Onset parameter (property value at x_c → 0).
    A   : float  Scale parameter.
    B   : float  Shape parameter.
    """

    def __init__(self) -> None:
        self.P_o: float | None = None
        self.A:   float | None = None
        self.B:   float | None = None
        self._mode:   str | None = None
        self._xc_fit: np.ndarray | None = None
        self._P_fit:  np.ndarray | None = None

    # ── Fitting ───────────────────────────────────────────────────────────────

    def fit(
        self,
        xc,
        P_data,
        mode: str = '3param',
        B_fixed: float = 1.5,
    ) -> 'GeneralizedDistribution':
        """Fit the generalized distribution to discrete (x_c, P) data.

        Parameters
        ----------
        xc : array-like
            Cumulative fraction values in the open interval (0, 1).
            For weight-basis fits, pass cumulative weight fractions directly.
        P_data : array-like
            Property values (Tb in K, M in g/mol, or SG) at each x_c.
        mode : str
            '3param' — fit P_o, A, B freely (default).
            '2param' — fix B = B_fixed and fit P_o, A only.
        B_fixed : float
            B value for two-parameter fit.  Riazi (2005) p. 174 recommends
            B = 1.5 for Tb, B = 1 for M, B = 3 for SG.

        Returns
        -------
        GeneralizedDistribution
            self (allows chaining).

        References
        ----------
        Riazi MNL50 §4.5.4.1, Eqs. 4.56–4.57 (PDF pages ~191–194).
        Default B values per Riazi p. 174.
        """
        xc     = np.asarray(xc,     dtype=float)
        P_data = np.asarray(P_data, dtype=float)

        if xc.ndim != 1 or xc.shape != P_data.shape:
            raise ValueError(
                f"xc and P_data must be 1-D arrays of equal length; "
                f"got shapes {xc.shape} and {P_data.shape}."
            )
        if xc.size < 2:
            raise ValueError("At least 2 data points are required for fitting.")
        if np.any(xc <= 0.0) or np.any(xc >= 1.0):
            bad_lo = int(np.sum(xc <= 0.0))
            bad_hi = int(np.sum(xc >= 1.0))
            raise ValueError(
                f"xc values must be in the open interval (0, 1); "
                f"found {bad_lo} value(s) <= 0 and {bad_hi} value(s) >= 1."
            )
        if mode not in ('3param', '2param'):
            raise ValueError(
                f"mode must be '3param' or '2param'; got {mode!r}."
            )

        self._mode   = mode
        self._xc_fit = xc.copy()
        self._P_fit  = P_data.copy()

        P_min = float(P_data.min())
        Y = np.log(-np.log(1.0 - xc))      # linearized x_c side; shape (n,)

        if mode == '3param':
            P_o, A, B = self._fit_3param(xc, P_data, Y, P_min)
        else:
            P_o, A, B = self._fit_2param(xc, P_data, Y, P_min, B_fixed)

        self.P_o, self.A, self.B = P_o, A, B
        return self

    # ── Internal fit helpers ──────────────────────────────────────────────────

    def _fit_3param(
        self,
        xc: np.ndarray,
        P_data: np.ndarray,
        Y: np.ndarray,
        P_min: float,
    ) -> tuple[float, float, float]:
        """Three-parameter fit: P_o, A, B free.

        Minimizes sum of squared linearized residuals r_i = Y_i - B*X_i - c
        where X_i = ln((P_i - P_o) / P_o).  Returns (P_o, A, B).
        """
        def residuals(params):
            P_o, B, c = params
            if P_o <= 0.0 or P_o >= P_min or B <= 0.0:
                return np.full(len(xc), 1e6)
            P_star = (P_data - P_o) / P_o
            if np.any(P_star <= 0.0):
                return np.full(len(xc), 1e6)
            X = np.log(P_star)
            return Y - (B * X + c)

        # Smart initial guess: OLS regression at P_o = 0.95 * P_min
        P_o0 = 0.95 * P_min
        X0 = np.log((P_data - P_o0) / P_o0)
        X0_mean, Y_mean = X0.mean(), Y.mean()
        B0 = float(np.sum((X0 - X0_mean) * (Y - Y_mean)) / np.sum((X0 - X0_mean) ** 2))
        c0 = float(Y_mean - B0 * X0_mean)

        lb = [1e-3,         0.05, -50.0]
        ub = [P_min - 1e-6, 20.0,  50.0]

        result = optimize.least_squares(
            residuals,
            [P_o0, B0, c0],
            bounds=(lb, ub),
            method='trf',
            ftol=1e-12, xtol=1e-12, gtol=1e-12,
            max_nfev=20000,
        )

        P_o, B, c = result.x
        A = B * np.exp(-c)
        return float(P_o), float(A), float(B)

    def _fit_2param(
        self,
        xc: np.ndarray,
        P_data: np.ndarray,
        Y: np.ndarray,
        P_min: float,
        B_fixed: float,
    ) -> tuple[float, float, float]:
        """Two-parameter fit: B fixed at B_fixed; fit P_o and A.

        Minimizes sum of squared linearized residuals r_i = Y_i - B_fixed*X_i - c.
        Returns (P_o, A, B_fixed).
        """
        def residuals(params):
            P_o, c = params
            if P_o <= 0.0 or P_o >= P_min:
                return np.full(len(xc), 1e6)
            P_star = (P_data - P_o) / P_o
            if np.any(P_star <= 0.0):
                return np.full(len(xc), 1e6)
            X = np.log(P_star)
            return Y - (B_fixed * X + c)

        P_o0 = 0.95 * P_min
        X0 = np.log((P_data - P_o0) / P_o0)
        c0 = float(np.mean(Y - B_fixed * X0))

        lb = [1e-3,         -50.0]
        ub = [P_min - 1e-6,  50.0]

        result = optimize.least_squares(
            residuals,
            [P_o0, c0],
            bounds=(lb, ub),
            method='trf',
            ftol=1e-12, xtol=1e-12, gtol=1e-12,
            max_nfev=20000,
        )

        P_o, c = result.x
        A = B_fixed * np.exp(-c)
        return float(P_o), float(A), float(B_fixed)

    # ── Evaluation methods ────────────────────────────────────────────────────

    def P(self, xc):
        """Property value at cumulative fraction x_c (forward CDF).

        P(x_c) = P_o * (1 + [(A/B) * (-ln(1 - x_c))]^(1/B))

        Parameters
        ----------
        xc : float or array-like
            Cumulative fraction in (0, 1).

        Returns
        -------
        float or ndarray
            Property value.

        References
        ----------
        Riazi MNL50 Eq. 4.56.
        """
        self._require_fit()
        scalar = np.ndim(xc) == 0
        xc = np.asarray(xc, dtype=float)
        P_star = ((self.A / self.B) * (-np.log(1.0 - xc))) ** (1.0 / self.B)
        result = self.P_o * (1.0 + P_star)
        return float(result) if scalar else result

    def x_c(self, P):
        """Cumulative fraction at property value P (inverse CDF).

        x_c(P) = 1 - exp(-(B/A) * [(P - P_o)/P_o]^B)

        Parameters
        ----------
        P : float or array-like
            Property value. Must satisfy P > P_o.

        Returns
        -------
        float or ndarray
            Cumulative fraction in (0, 1).

        References
        ----------
        Riazi MNL50 Eq. 4.56 inverted analytically.
        """
        self._require_fit()
        scalar = np.ndim(P) == 0
        P = np.asarray(P, dtype=float)
        P_star = (P - self.P_o) / self.P_o
        result = 1.0 - np.exp(-(self.B / self.A) * P_star ** self.B)
        return float(result) if scalar else result

    def pdf(self, P):
        """Probability density function f(P) = dx_c/dP in property space.

        f(P) = (B^2 / (A * P_o)) * P*^(B-1) * exp(-(B/A) * P*^B)

        where P* = (P - P_o) / P_o.

        Parameters
        ----------
        P : float or array-like
            Property value. Must satisfy P > P_o.

        Returns
        -------
        float or ndarray
            PDF value (units: 1/[P units]).

        References
        ----------
        Riazi MNL50 Eq. 4.56 differentiated with respect to P.
        """
        self._require_fit()
        scalar = np.ndim(P) == 0
        P = np.asarray(P, dtype=float)
        P_star = (P - self.P_o) / self.P_o
        result = (
            (self.B ** 2 / (self.A * self.P_o))
            * P_star ** (self.B - 1.0)
            * np.exp(-(self.B / self.A) * P_star ** self.B)
        )
        return float(result) if scalar else result

    def average(self) -> float:
        """Analytical mean property value.

        <P> = P_o * (1 + (A/B)^(1/B) * Gamma(1 + 1/B))

        Derived by substituting u = -ln(1-x_c) in the integral
        <P> = int_0^1 P(x_c) dx_c.

        Returns
        -------
        float
            Mean of P over the fitted distribution.

        References
        ----------
        Riazi MNL50 §4.5.4.1; result from standard gamma-integral identity.
        """
        self._require_fit()
        return float(
            self.P_o * (1.0 + (self.A / self.B) ** (1.0 / self.B) * special.gamma(1.0 + 1.0 / self.B))
        )

    @property
    def params(self) -> dict:
        """Fitted parameters: {'P_o', 'A', 'B', 'mode'}."""
        self._require_fit()
        return {'P_o': self.P_o, 'A': self.A, 'B': self.B, 'mode': self._mode}

    @property
    def fit_quality(self) -> dict:
        """Fit quality metrics evaluated on the original (un-transformed) P scale.

        Returns
        -------
        dict
            'RMS'       : root-mean-square error in P units.
            'pct_AAD'   : mean absolute percent deviation (%).
            'R_squared' : coefficient of determination (R²) in P space.

        References
        ----------
        Riazi MNL50 Table 4.13 header for metric conventions.
        """
        self._require_fit()
        P_calc    = self.P(self._xc_fit)
        residuals = P_calc - self._P_fit
        rms       = float(np.sqrt(np.mean(residuals ** 2)))
        pct_aad   = float(np.mean(np.abs(residuals / self._P_fit)) * 100.0)
        ss_res    = float(np.sum(residuals ** 2))
        ss_tot    = float(np.sum((self._P_fit - self._P_fit.mean()) ** 2))
        r2        = float(1.0 - ss_res / ss_tot)
        return {'RMS': rms, 'pct_AAD': pct_aad, 'R_squared': r2}

    # ── Utility ───────────────────────────────────────────────────────────────

    def _require_fit(self) -> None:
        if self.P_o is None:
            raise RuntimeError(
                "GeneralizedDistribution has not been fitted. Call fit() first."
            )
