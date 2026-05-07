"""
Distillation curve container and D86→TBP inter-conversion.

Internal units: temperatures in Kelvin, cumulative percent in [0, 100].

Unit convention for Daubert's method:
  - Eq. 3.20 operates in Kelvin (both input and output).
  - Eq. 3.21 operates on temperature *differences* (K ≡ °C difference — same
    numerical value).
  - Table 3.7 max-delta values are in °C; applied as K caps numerically.

References
----------
Riazi (2005)    : §3.2.2.2 Daubert D86→TBP, Eqs. 3.20-3.22, Table 3.7.
                  Pages 103-104 (PDF pages 122-123).
Daubert (1994)  : API Technical Data Book, 6th ed., API, Washington DC.
                  Primary source for Table 3.7 constants; tabulated in Riazi.
"""

import numpy as np


# ── Daubert Eq. 3.20 constants ───────────────────────────────────────────────
# TBP(50%, K) = 255.4 + 0.8851 × [D86(50%, K) − 255.4]^{1.0258}
# Source: Riazi MNL50 §3.2.2.2.1, Eq. 3.20 (page 103, PDF page 122).
# Exponent 1.0258 is verified numerically:
#   D86(50%) = 479.85 K (206.7°C) → TBP(50%) = 483.75 K (210.6°C)
#   matches Table 3.8 kerosene result in Example 3.3.
_EQ320_SHIFT = 255.4    # K
_EQ320_COEFF = 0.8851
_EQ320_EXP   = 1.0258

# ── Table 3.7 — Daubert Eq. 3.21 constants (Y_i = A_i × X_i^B_i) ───────────
# Rows indexed 0-5 correspond to cut-point intervals:
#   row 0: 100-90%,  row 1: 90-70%,  row 2: 70-50%
#   row 3: 50-30%,   row 4: 30-10%,  row 5: 10-0%
# Each tuple: (A, B, max_delta) where max_delta is the cap on X_i in °C (≡ K).
# Source: Riazi MNL50 Table 3.7, page 103 (PDF page 122).
# Verified against Example 3.3 kerosene: Y4=21.6°C, Y3=22.2°C, Y2=25.3°C.
_TABLE_3_7 = (
    (0.1403, 1.6606, 55.0),   # row 0: 100-90%
    (2.6339, 0.7550, 85.0),   # row 1: 90-70%
    (2.2744, 0.8200, 140.0),  # row 2: 70-50%
    (2.6956, 0.8008, 140.0),  # row 3: 50-30%
    (4.1481, 0.7164, 55.0),   # row 4: 30-10%
    (5.8589, 0.6024, 55.0),   # row 5: 10-0%  (max estimated; last row OCR unclear)
)

# Standard cumulative percent cut points required for Daubert's method.
_DAUBERT_PCT = (0.0, 10.0, 30.0, 50.0, 70.0, 90.0, 100.0)


class DistillationCurve:
    """Petroleum distillation curve with explicit method and basis tags.

    Stores cumulative percent distilled vs boiling temperature.  All
    temperatures are in Kelvin internally.

    Parameters
    ----------
    pct : array-like
        Cumulative percent distilled, values in [0, 100].
    temps_K : array-like
        Boiling temperatures at each percent point, K.
    method : str
        Distillation method: 'D86', 'D1160_AET', 'D7169', or 'TBP'.
    basis : str
        Distillation basis: 'volume', 'weight', or 'mole'.
    """

    _VALID_METHODS = frozenset({'D86', 'D1160_AET', 'D7169', 'TBP'})
    _VALID_BASES   = frozenset({'volume', 'weight', 'mole'})

    def __init__(self, pct, temps_K, method: str, basis: str) -> None:
        pct     = np.asarray(pct,     dtype=float)
        temps_K = np.asarray(temps_K, dtype=float)

        if pct.ndim != 1 or pct.shape != temps_K.shape:
            raise ValueError(
                "pct and temps_K must be 1-D arrays of equal length; "
                f"got shapes {pct.shape} and {temps_K.shape}."
            )
        if pct.size == 0:
            raise ValueError("pct and temps_K must not be empty.")
        if pct.min() < 0.0 or pct.max() > 100.0:
            raise ValueError(
                f"pct values must be in [0, 100]; "
                f"got min={pct.min():.2f}, max={pct.max():.2f}."
            )
        if method not in self._VALID_METHODS:
            raise ValueError(
                f"method must be one of {sorted(self._VALID_METHODS)}; "
                f"got {method!r}."
            )
        if basis not in self._VALID_BASES:
            raise ValueError(
                f"basis must be one of {sorted(self._VALID_BASES)}; "
                f"got {basis!r}."
            )

        self.pct     = pct
        self.temps_K = temps_K
        self.method  = method
        self.basis   = basis

    # ── Public conversion methods ─────────────────────────────────────────────

    def to_tbp(self) -> 'DistillationCurve':
        """Return the curve expressed as a TBP distillation.

        Dispatch table:
          - TBP       → copy unchanged (already TBP).
          - D1160_AET → pass-through; temperatures already at atmospheric
                        equivalent boiling point, equivalent to TBP for
                        characterization purposes.
          - D86       → Daubert Eqs. 3.20-3.22 (requires 7 standard cut points).
          - D7169     → not implemented in Phase 2; raises ValueError.

        Returns
        -------
        DistillationCurve
            New instance with method='TBP' and the same basis.
        """
        if self.method == 'TBP':
            return DistillationCurve(
                self.pct.copy(), self.temps_K.copy(), 'TBP', self.basis
            )
        if self.method == 'D1160_AET':
            # D1160 at atmospheric equivalent temperature maps directly to TBP.
            # No numerical conversion required.
            return DistillationCurve(
                self.pct.copy(), self.temps_K.copy(), 'TBP', self.basis
            )
        if self.method == 'D86':
            return self._d86_to_tbp_daubert()
        raise ValueError(
            f"to_tbp(): D7169 → TBP conversion is not implemented in Phase 2. "
            f"Got method={self.method!r}."
        )

    def to_weight_basis(self) -> 'DistillationCurve':
        """Convert to weight basis.

        If the curve is already on a weight basis, returns a copy unchanged.
        Volume and mole bases require per-cut SG (and MW for mole), which
        become available only in Phase 4 (SG distribution).

        Raises
        ------
        ValueError
            Always raised when basis is not 'weight', because the SG and MW
            distributions needed for the conversion are not yet available.
            This is the correct Phase-2 behaviour; the implementation is
            completed in Phase 4.
        """
        if self.basis == 'weight':
            return DistillationCurve(
                self.pct.copy(), self.temps_K.copy(), self.method, 'weight'
            )
        raise ValueError(
            f"to_weight_basis: '{self.basis}' → 'weight' conversion requires "
            f"the SG distribution (volume basis) or MW distribution (mole basis), "
            f"neither of which is available until Phase 4. "
            f"Input basis: {self.basis!r}."
        )

    # ── Internal conversion ───────────────────────────────────────────────────

    def _d86_to_tbp_daubert(self) -> 'DistillationCurve':
        """Daubert D86→TBP via Riazi MNL50 §3.2.2.2.1, Eqs. 3.20-3.22.

        Algorithm
        ---------
        1. Eq. 3.20 — anchor: TBP(50%, K) from D86(50%, K).
        2. Eq. 3.21 — gaps: Y_i = A_i × X_i^B_i where X_i is the D86 gap
           at interval i (K ≡ °C).  X_i is capped at max_delta (Table 3.7).
        3. Eq. 3.22 — summation: TBP at all other cut points by adding/
           subtracting the Y_i values around the 50% anchor.

        Parameters (via self)
        ---------------------
        self.pct must include exactly 0, 10, 30, 50, 70, 90, 100%.
        D86 temperatures must be monotonically increasing.
        """
        # ── Validate that all 7 standard cut points are present ──────────────
        required = set(_DAUBERT_PCT)
        present  = set(self.pct.tolist())
        missing  = required - present
        if missing:
            raise ValueError(
                f"D86→TBP via Daubert requires cut points at "
                f"{sorted(required)}%; missing: {sorted(missing)}%."
            )

        # ── Extract D86 temperatures at the 7 cut points ─────────────────────
        def _T(p: float) -> float:
            idx = int(np.where(self.pct == p)[0][0])
            return float(self.temps_K[idx])

        D = {p: _T(p) for p in _DAUBERT_PCT}

        # ── Eq. 3.20: TBP(50%) in K ──────────────────────────────────────────
        tbp_50 = (
            _EQ320_SHIFT
            + _EQ320_COEFF * (D[50.0] - _EQ320_SHIFT) ** _EQ320_EXP
        )

        # ── Eq. 3.21: temperature gaps (K / °C) ──────────────────────────────
        # X values: D86 temperature differences at each interval.
        # Intervals ordered top to bottom: 100-90, 90-70, 70-50, 50-30, 30-10, 10-0%.
        X = [
            D[100.0] - D[90.0],   # X1 for row 0 (100-90%)
            D[90.0]  - D[70.0],   # X2 for row 1 (90-70%)
            D[70.0]  - D[50.0],   # X3 for row 2 (70-50%)
            D[50.0]  - D[30.0],   # X4 for row 3 (50-30%)
            D[30.0]  - D[10.0],   # X5 for row 4 (30-10%)
            D[10.0]  - D[0.0],    # X6 for row 5 (10-0%)
        ]

        Y = []
        for i, (xi, (A, B, max_delta)) in enumerate(zip(X, _TABLE_3_7)):
            if xi < 0.0:
                raise ValueError(
                    f"D86 temperatures are not monotonically increasing: "
                    f"interval {i} (row {i+1} of Table 3.7) has gap "
                    f"X={xi:.3f} K < 0."
                )
            xi_eff = min(xi, max_delta)
            Y.append(float(A * xi_eff ** B))

        Y1, Y2, Y3, Y4, Y5, Y6 = Y[0], Y[1], Y[2], Y[3], Y[4], Y[5]

        # ── Eq. 3.22: TBP at all seven cut points ────────────────────────────
        tbp_map = {
            0.0:   tbp_50 - Y4 - Y5 - Y6,
            10.0:  tbp_50 - Y4 - Y5,
            30.0:  tbp_50 - Y4,
            50.0:  tbp_50,
            70.0:  tbp_50 + Y3,
            90.0:  tbp_50 + Y3 + Y2,
            100.0: tbp_50 + Y3 + Y2 + Y1,
        }

        tbp_pct  = np.array(sorted(tbp_map.keys()))
        tbp_temp = np.array([tbp_map[p] for p in tbp_pct])

        return DistillationCurve(tbp_pct, tbp_temp, 'TBP', self.basis)
