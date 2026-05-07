"""
Phase 2 pass-gate tests for core/distillation.py.

Pass-gate: Daubert D86->TBP conversion reproduces Example 3.3 kerosene
Table 3.8 results within +-5 K at all six tabulated cut points.

Additional coverage:
  - D1160_AET pass-through (no-op)
  - TBP.to_tbp() identity copy
  - Input validation (method, basis, cut points, monotonicity)
  - to_weight_basis error when basis != 'weight'
  - Synthetic D86 curve round-trip check

References
----------
Riazi MNL50 Table 3.8 (page 105): kerosene Example 3.3 D86 -> TBP via
Daubert Eqs. 3.20-3.22.  D86 data and experimental TBP are input;
Daubert-predicted TBP values are the check targets.
"""

import numpy as np
import pytest

from core.distillation import DistillationCurve


# ── Kerosene Example 3.3 data (Table 3.8, Riazi MNL50 page 105) ──────────────
# ASTM D86 input (experimental), in Celsius then converted to K.
# Daubert-predicted TBP values (Eqs. 3.20-3.22) from Table 3.8, col "TBP calc".
# Experimental TBP (not used in pass-gate, included for reference).
_KEROSENE_PCT   = (0.0,   10.0,  30.0,  50.0,  70.0,  90.0)
_KEROSENE_D86_C = (165.6, 176.7, 193.3, 206.7, 222.8, 242.8)
_KEROSENE_TBP_CALC_C = (133.1, 158.1, 189.2, 210.6, 232.9, 258.1)  # Daubert predicted
_KEROSENE_TBP_EXP_C  = (146.1, 160.6, 188.3, 209.4, 230.6, 255.0)  # experimental

# Table 3.8 includes 100% point in D86 but not in the summary columns.
# For the pass-gate, use a representative 100% D86 value (estimated from kerosene data).
# Riazi MNL50 Example 3.3 does not list D86(100%); we synthesize it consistently.
# Synthesis: extend linearly from 90% using the 70-90 slope.
_D86_90 = _KEROSENE_D86_C[5]  # 242.8 C
_D86_70 = _KEROSENE_D86_C[4]  # 222.8 C
_D86_100_C = _D86_90 + (_D86_90 - _D86_70)  # = 262.8 C (linear extension)

_K = 273.15  # offset


def _make_kerosene_d86() -> DistillationCurve:
    """Build the kerosene D86 input curve from Table 3.8 data."""
    pct_full = list(_KEROSENE_PCT) + [100.0]
    d86_C    = list(_KEROSENE_D86_C) + [_D86_100_C]
    temps_K  = [T + _K for T in d86_C]
    return DistillationCurve(pct_full, temps_K, method='D86', basis='volume')


# ── Pass-gate: Example 3.3 kerosene (Table 3.8) ──────────────────────────────

class TestDaubertKeroseneExample33:
    """Reproduce Table 3.8 kerosene results within +-5 K."""

    @pytest.fixture(scope='class')
    def tbp_kerosene(self):
        return _make_kerosene_d86().to_tbp()

    def test_method_is_tbp(self, tbp_kerosene):
        assert tbp_kerosene.method == 'TBP'

    def test_basis_unchanged(self, tbp_kerosene):
        assert tbp_kerosene.basis == 'volume'

    @pytest.mark.parametrize('pct,tbp_calc_C', zip(
        _KEROSENE_PCT, _KEROSENE_TBP_CALC_C
    ))
    def test_tbp_within_5K(self, tbp_kerosene, pct, tbp_calc_C):
        """Each cut point: |computed - Table 3.8 Daubert value| <= 5 K."""
        idx = int(np.where(tbp_kerosene.pct == pct)[0][0])
        tbp_K = float(tbp_kerosene.temps_K[idx])
        ref_K = tbp_calc_C + _K
        delta = abs(tbp_K - ref_K)
        assert delta <= 5.0, (
            f"pct={pct}%: TBP_calc={tbp_K - _K:.2f} C, "
            f"Table_3.8={tbp_calc_C} C, |delta|={delta:.2f} K"
        )

    def test_50pct_within_1K(self, tbp_kerosene):
        """50% point (direct Eq. 3.20) should agree with Table 3.8 within 1 K."""
        idx = int(np.where(tbp_kerosene.pct == 50.0)[0][0])
        tbp_K = float(tbp_kerosene.temps_K[idx])
        ref_K = 210.6 + _K
        delta = abs(tbp_K - ref_K)
        assert delta <= 1.0, (
            f"50% TBP={tbp_K - _K:.2f} C, ref=210.6 C, |delta|={delta:.2f} K"
        )

    def test_max_deviation_report(self, tbp_kerosene, capsys):
        """Print deviations for all tabulated cut points (informational)."""
        print("\nExample 3.3 kerosene TBP deviations (vs Table 3.8 Daubert):")
        print(f"{'pct%':>6} {'D86 C':>8} {'TBP_ref C':>10} {'TBP_calc C':>11} {'delta K':>8}")
        max_delta = 0.0
        for pct, d86_C, tbp_ref_C in zip(
            _KEROSENE_PCT, _KEROSENE_D86_C, _KEROSENE_TBP_CALC_C
        ):
            idx = int(np.where(tbp_kerosene.pct == pct)[0][0])
            tbp_calc_C = float(tbp_kerosene.temps_K[idx]) - _K
            delta = tbp_calc_C - tbp_ref_C
            max_delta = max(max_delta, abs(delta))
            print(f"{pct:>6.0f} {d86_C:>8.1f} {tbp_ref_C:>10.1f} {tbp_calc_C:>11.2f} {delta:>+8.2f}")
        print(f"\nMax |delta| = {max_delta:.2f} K  (pass gate: <= 5 K)")
        assert max_delta <= 5.0


# ── D1160_AET pass-through ────────────────────────────────────────────────────

class TestD1160AETPassThrough:
    """D1160_AET.to_tbp() returns identical temperatures, method='TBP'."""

    def test_method_becomes_tbp(self):
        pct = [0.0, 50.0, 100.0]
        T   = [500.0, 600.0, 700.0]
        c = DistillationCurve(pct, T, method='D1160_AET', basis='weight')
        tbp = c.to_tbp()
        assert tbp.method == 'TBP'

    def test_temperatures_unchanged(self):
        pct = [0.0, 10.0, 30.0, 50.0, 70.0, 90.0, 100.0]
        T   = [450.0, 480.0, 530.0, 580.0, 630.0, 700.0, 780.0]
        c = DistillationCurve(pct, T, method='D1160_AET', basis='volume')
        tbp = c.to_tbp()
        np.testing.assert_array_equal(tbp.temps_K, c.temps_K)

    def test_basis_preserved(self):
        pct = [0.0, 50.0, 100.0]
        T   = [500.0, 600.0, 700.0]
        for basis in ('volume', 'weight', 'mole'):
            c = DistillationCurve(pct, T, method='D1160_AET', basis=basis)
            assert c.to_tbp().basis == basis

    def test_returns_new_instance(self):
        pct = [0.0, 50.0, 100.0]
        T   = [500.0, 600.0, 700.0]
        c   = DistillationCurve(pct, T, method='D1160_AET', basis='weight')
        tbp = c.to_tbp()
        assert tbp is not c


# ── TBP identity ──────────────────────────────────────────────────────────────

class TestTBPIdentity:
    """TBP.to_tbp() returns a copy with unchanged temperatures."""

    def test_method_unchanged(self):
        pct = [0.0, 50.0, 100.0]
        T   = [500.0, 600.0, 700.0]
        c   = DistillationCurve(pct, T, method='TBP', basis='weight')
        tbp = c.to_tbp()
        assert tbp.method == 'TBP'

    def test_temperatures_unchanged(self):
        pct = [0.0, 10.0, 30.0, 50.0, 70.0, 90.0, 100.0]
        T   = [450.0, 480.0, 530.0, 580.0, 630.0, 700.0, 780.0]
        c   = DistillationCurve(pct, T, method='TBP', basis='volume')
        tbp = c.to_tbp()
        np.testing.assert_array_equal(tbp.temps_K, c.temps_K)

    def test_returns_new_instance(self):
        pct = [0.0, 50.0, 100.0]
        T   = [500.0, 600.0, 700.0]
        c   = DistillationCurve(pct, T, method='TBP', basis='weight')
        assert c.to_tbp() is not c


# ── Input validation ──────────────────────────────────────────────────────────

class TestInputValidation:
    def test_invalid_method(self):
        with pytest.raises(ValueError, match="method"):
            DistillationCurve([0.0, 50.0], [300.0, 500.0], method='ASTM', basis='volume')

    def test_invalid_basis(self):
        with pytest.raises(ValueError, match="basis"):
            DistillationCurve([0.0, 50.0], [300.0, 500.0], method='D86', basis='mass')

    def test_shape_mismatch(self):
        with pytest.raises(ValueError):
            DistillationCurve([0.0, 50.0], [300.0, 400.0, 500.0], method='D86', basis='volume')

    def test_pct_out_of_range(self):
        with pytest.raises(ValueError, match="\\[0, 100\\]"):
            DistillationCurve([-5.0, 50.0], [300.0, 500.0], method='D86', basis='volume')

    def test_empty_arrays(self):
        with pytest.raises(ValueError):
            DistillationCurve([], [], method='D86', basis='volume')

    def test_d7169_to_tbp_raises(self):
        pct = [0.0, 10.0, 30.0, 50.0, 70.0, 90.0, 100.0]
        T   = [400.0, 430.0, 470.0, 520.0, 570.0, 640.0, 730.0]
        c   = DistillationCurve(pct, T, method='D7169', basis='weight')
        with pytest.raises(ValueError, match="D7169"):
            c.to_tbp()


# ── D86 conversion error cases ────────────────────────────────────────────────

class TestD86ConversionErrors:
    def test_missing_cut_points(self):
        """Omitting a required cut point raises ValueError."""
        pct = [0.0, 10.0, 50.0, 70.0, 90.0, 100.0]  # missing 30%
        T   = [400.0, 430.0, 520.0, 570.0, 640.0, 730.0]
        c   = DistillationCurve(pct, T, method='D86', basis='volume')
        with pytest.raises(ValueError, match="30"):
            c.to_tbp()

    def test_non_monotonic_d86(self):
        """Non-monotonic D86 (decreasing gap) raises ValueError."""
        pct = [0.0, 10.0, 30.0, 50.0, 70.0, 90.0, 100.0]
        # Make D86(70%) < D86(50%) to create a negative gap
        T   = [400.0, 430.0, 470.0, 520.0, 510.0, 560.0, 640.0]
        c   = DistillationCurve(pct, T, method='D86', basis='volume')
        with pytest.raises(ValueError, match="monoton"):
            c.to_tbp()


# ── to_weight_basis errors ────────────────────────────────────────────────────

class TestToWeightBasisErrors:
    def test_volume_raises(self):
        pct = [0.0, 50.0, 100.0]
        T   = [400.0, 500.0, 600.0]
        c   = DistillationCurve(pct, T, method='TBP', basis='volume')
        with pytest.raises(ValueError, match="weight"):
            c.to_weight_basis()

    def test_mole_raises(self):
        pct = [0.0, 50.0, 100.0]
        T   = [400.0, 500.0, 600.0]
        c   = DistillationCurve(pct, T, method='TBP', basis='mole')
        with pytest.raises(ValueError, match="weight"):
            c.to_weight_basis()

    def test_weight_no_op(self):
        """Weight basis returns a copy without error."""
        pct = [0.0, 50.0, 100.0]
        T   = [400.0, 500.0, 600.0]
        c   = DistillationCurve(pct, T, method='TBP', basis='weight')
        result = c.to_weight_basis()
        assert result.basis == 'weight'
        np.testing.assert_array_equal(result.temps_K, c.temps_K)


# ── Synthetic D86 curve pass-gate (independent of kerosene data) ──────────────

class TestSyntheticD86PassGate:
    """Pass-gate using a hand-computed reference TBP curve.

    Reference TBP is computed independently using the Daubert formulae with
    a synthetic D86 input chosen to give representative petroleum fractions.

    Tolerance: +-5 K per cut point (spec criterion for Phase 2).
    """

    # Synthetic D86 input (typical medium gas oil, volume basis), in K.
    _D86_K = {
        0.0:   273.15 + 180.0,   # 453.15 K
        10.0:  273.15 + 210.0,   # 483.15 K
        30.0:  273.15 + 245.0,   # 518.15 K
        50.0:  273.15 + 275.0,   # 548.15 K
        70.0:  273.15 + 305.0,   # 578.15 K
        90.0:  273.15 + 345.0,   # 618.15 K
        100.0: 273.15 + 380.0,   # 653.15 K
    }

    @classmethod
    def _hand_compute_tbp(cls) -> dict:
        """Hand-compute expected TBP using Daubert formulae."""
        D = cls._D86_K

        # Eq. 3.20
        from core.distillation import _EQ320_SHIFT, _EQ320_COEFF, _EQ320_EXP
        tbp_50 = (
            _EQ320_SHIFT
            + _EQ320_COEFF * (D[50.0] - _EQ320_SHIFT) ** _EQ320_EXP
        )

        # Eq. 3.21
        from core.distillation import _TABLE_3_7
        X = [
            D[100.0] - D[90.0],
            D[90.0]  - D[70.0],
            D[70.0]  - D[50.0],
            D[50.0]  - D[30.0],
            D[30.0]  - D[10.0],
            D[10.0]  - D[0.0],
        ]
        Y = [A * min(xi, md) ** B for xi, (A, B, md) in zip(X, _TABLE_3_7)]
        Y1, Y2, Y3, Y4, Y5, Y6 = Y

        # Eq. 3.22
        return {
            0.0:   tbp_50 - Y4 - Y5 - Y6,
            10.0:  tbp_50 - Y4 - Y5,
            30.0:  tbp_50 - Y4,
            50.0:  tbp_50,
            70.0:  tbp_50 + Y3,
            90.0:  tbp_50 + Y3 + Y2,
            100.0: tbp_50 + Y3 + Y2 + Y1,
        }

    def test_synthetic_within_5K(self):
        """Computed TBP matches hand-computed reference within +-5 K."""
        D = self._D86_K
        pct    = sorted(D.keys())
        temps  = [D[p] for p in pct]
        curve  = DistillationCurve(pct, temps, method='D86', basis='volume')
        tbp    = curve.to_tbp()
        ref    = self._hand_compute_tbp()

        for p in pct:
            idx     = int(np.where(tbp.pct == p)[0][0])
            tbp_K   = float(tbp.temps_K[idx])
            ref_K   = ref[p]
            delta   = abs(tbp_K - ref_K)
            assert delta < 1e-9, (
                f"pct={p}%: computed={tbp_K:.4f} K, expected={ref_K:.4f} K, "
                f"|delta|={delta:.2e} K (should be exact)"
            )

    def test_tbp_monotonically_increasing(self):
        """TBP temperatures must increase with percent distilled."""
        D = self._D86_K
        pct   = sorted(D.keys())
        temps = [D[p] for p in pct]
        tbp   = DistillationCurve(pct, temps, 'D86', 'volume').to_tbp()
        diffs = np.diff(tbp.temps_K)
        assert np.all(diffs > 0), (
            f"TBP is not monotone: diffs={diffs.tolist()}"
        )
