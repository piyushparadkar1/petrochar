"""
Phase 3 pass-gate: reproduce Riazi Table 4.13 distribution coefficients
from Table 4.11 North Sea gas condensate C7+ data.

Pass-gate criteria (Riazi MNL50 §4.5.4.1, Table 4.13):
  Three-parameter Tb fit: T_o=350 K, A_T=0.1679, B_T=1.2586
  Two-parameter Tb fit (B_T=1.5 fixed): T_o=340 K, A_T=0.1875
  Each parameter within ±5% of table value.

Data source:
  Table 4.11 columns Tb_K (property) and x_cw (cumulative weight fraction).
  Rows with NaN Tb are excluded (C18+, row 12).

References
----------
Riazi MNL50 Table 4.11 (page 171) — input data.
Riazi MNL50 Table 4.13 (page 173) — expected coefficient output.
Riazi MNL50 §4.5.4.1  (pages 172–175) — model equations.
"""

import pathlib

import numpy as np
import pandas as pd
import pytest

from core.distribution import GeneralizedDistribution

DATA_DIR = pathlib.Path(__file__).parent.parent / "data" / "riazi_reference"


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def table_411():
    df = pd.read_csv(DATA_DIR / "table_4_11_example_4_7.csv", comment="#")
    return df.dropna(subset=["Tb_K"])


@pytest.fixture(scope="module")
def tb_data(table_411):
    xc     = table_411["x_cw"].values.astype(float)
    Tb     = table_411["Tb_K"].values.astype(float)
    return xc, Tb


@pytest.fixture(scope="module")
def dist_3param(tb_data):
    xc, Tb = tb_data
    return GeneralizedDistribution().fit(xc, Tb, mode='3param')


@pytest.fixture(scope="module")
def dist_2param(tb_data):
    xc, Tb = tb_data
    return GeneralizedDistribution().fit(xc, Tb, mode='2param', B_fixed=1.5)


# ── Data integrity ────────────────────────────────────────────────────────────

class TestDataIntegrity:
    def test_row_count(self, table_411):
        assert len(table_411) == 11, (
            f"Expected 11 non-NaN Tb rows from Table 4.11; got {len(table_411)}"
        )

    def test_xc_range(self, tb_data):
        xc, _ = tb_data
        assert xc.min() > 0.0
        assert xc.max() < 1.0

    def test_Tb_increasing(self, tb_data):
        xc, Tb = tb_data
        assert np.all(np.diff(Tb) > 0), "Tb must be monotonically increasing with x_cw"

    def test_first_row_c7(self, table_411):
        row = table_411[table_411["Carbon_No_label"] == "C7"].iloc[0]
        assert abs(row["Tb_K"] - 365) < 0.5
        assert abs(row["x_cw"] - 0.130) < 0.001


# ── Three-parameter fit pass-gate ─────────────────────────────────────────────

class TestThreeParamFit:
    """Pass-gate: reproduce Table 4.13 three-parameter row (P_o=350, A=0.1679, B=1.2586)."""

    _T_o_ref = 350.0
    _A_ref   = 0.1679
    _B_ref   = 1.2586
    _TOL     = 0.05     # 5% tolerance per Riazi pass-gate spec

    def test_T_o_within_5pct(self, dist_3param):
        T_o = dist_3param.P_o
        err = abs(T_o - self._T_o_ref) / self._T_o_ref
        assert err <= self._TOL, (
            f"T_o={T_o:.2f} K vs reference {self._T_o_ref} K; "
            f"relative error {err*100:.2f}% > 5% tolerance"
        )

    def test_A_within_5pct(self, dist_3param):
        A = dist_3param.A
        err = abs(A - self._A_ref) / self._A_ref
        assert err <= self._TOL, (
            f"A={A:.5f} vs reference {self._A_ref}; "
            f"relative error {err*100:.2f}% > 5% tolerance"
        )

    def test_B_within_5pct(self, dist_3param):
        B = dist_3param.B
        err = abs(B - self._B_ref) / self._B_ref
        assert err <= self._TOL, (
            f"B={B:.5f} vs reference {self._B_ref}; "
            f"relative error {err*100:.2f}% > 5% tolerance"
        )

    def test_mode_is_3param(self, dist_3param):
        assert dist_3param.params['mode'] == '3param'

    def test_fit_quality_R2_above_0p99(self, dist_3param):
        fq = dist_3param.fit_quality
        assert fq['R_squared'] >= 0.99, (
            f"R^2 = {fq['R_squared']:.5f} < 0.99 for three-param Tb fit"
        )

    def test_fit_quality_pct_AAD_below_2(self, dist_3param):
        fq = dist_3param.fit_quality
        assert fq['pct_AAD'] < 2.0, (
            f"%AAD = {fq['pct_AAD']:.3f}% > 2.0% threshold "
            f"(Riazi Table 4.13 shows 0.62%)"
        )

    def test_print_params(self, dist_3param, capsys):
        p  = dist_3param.params
        fq = dist_3param.fit_quality
        print("\nThree-parameter Tb fit (Table 4.13 reference: T_o=350, A=0.1679, B=1.2586):")
        print(f"  T_o   = {p['P_o']:.3f} K  (ref 350 K)")
        print(f"  A     = {p['A']:.5f}    (ref 0.1679)")
        print(f"  B     = {p['B']:.5f}    (ref 1.2586)")
        print(f"  RMS   = {fq['RMS']:.4f} K  (ref 3.794 K)")
        print(f"  %AAD  = {fq['pct_AAD']:.4f}%  (ref 0.62%)")
        print(f"  R^2   = {fq['R_squared']:.6f}  (ref 0.998)")


# ── Two-parameter fit pass-gate ───────────────────────────────────────────────

class TestTwoParamFit:
    """Pass-gate: reproduce Table 4.13 two-parameter row (P_o=340, A=0.1875, B=1.5 fixed)."""

    _T_o_ref   = 340.0
    _A_ref     = 0.1875
    _B_expected = 1.5
    _TOL       = 0.05

    def test_B_is_fixed_at_1p5(self, dist_2param):
        assert abs(dist_2param.B - self._B_expected) < 1e-10, (
            f"B should be fixed at 1.5 for two-param fit; got {dist_2param.B}"
        )

    def test_T_o_within_5pct(self, dist_2param):
        T_o = dist_2param.P_o
        err = abs(T_o - self._T_o_ref) / self._T_o_ref
        assert err <= self._TOL, (
            f"T_o={T_o:.2f} K vs reference {self._T_o_ref} K; "
            f"relative error {err*100:.2f}% > 5% tolerance"
        )

    def test_A_within_5pct(self, dist_2param):
        A = dist_2param.A
        err = abs(A - self._A_ref) / self._A_ref
        assert err <= self._TOL, (
            f"A={A:.5f} vs reference {self._A_ref}; "
            f"relative error {err*100:.2f}% > 5% tolerance"
        )

    def test_mode_is_2param(self, dist_2param):
        assert dist_2param.params['mode'] == '2param'

    def test_fit_quality_R2_above_0p98(self, dist_2param):
        fq = dist_2param.fit_quality
        assert fq['R_squared'] >= 0.98, (
            f"R^2 = {fq['R_squared']:.5f} < 0.98 for two-param Tb fit"
        )

    def test_print_params(self, dist_2param, capsys):
        p  = dist_2param.params
        fq = dist_2param.fit_quality
        print("\nTwo-parameter Tb fit (Table 4.13 reference: T_o=340, A=0.1875, B=1.5 fixed):")
        print(f"  T_o   = {p['P_o']:.3f} K  (ref 340 K)")
        print(f"  A     = {p['A']:.5f}    (ref 0.1875)")
        print(f"  B     = {p['B']:.5f}    (fixed 1.5)")
        print(f"  RMS   = {fq['RMS']:.4f} K  (ref 5.834 K)")
        print(f"  %AAD  = {fq['pct_AAD']:.4f}%  (ref 1.15%)")
        print(f"  R^2   = {fq['R_squared']:.6f}  (ref 0.993)")


# ── Evaluation methods ────────────────────────────────────────────────────────

class TestEvaluationMethods:
    def test_P_at_zero_cumulative_approaches_P_o(self, dist_3param):
        """P(x_c) → P_o as x_c → 0."""
        P_small = dist_3param.P(1e-9)
        assert abs(P_small - dist_3param.P_o) < 0.1, (
            f"P(1e-9) = {P_small:.3f} K; P_o = {dist_3param.P_o:.3f} K"
        )

    def test_P_increases_with_xc(self, dist_3param):
        xc_vals = np.linspace(0.05, 0.95, 20)
        P_vals = dist_3param.P(xc_vals)
        assert np.all(np.diff(P_vals) > 0), "P(x_c) must be monotonically increasing"

    def test_xc_roundtrip(self, dist_3param):
        """x_c(P(x_c)) ≈ x_c for a range of cumulative fractions."""
        xc_vals = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        P_vals  = dist_3param.P(xc_vals)
        xc_back = dist_3param.x_c(P_vals)
        np.testing.assert_allclose(xc_back, xc_vals, rtol=1e-8,
                                   err_msg="x_c(P(x_c)) round-trip failed")

    def test_pdf_positive(self, dist_3param):
        P_vals = np.linspace(dist_3param.P_o + 1.0, 700.0, 50)
        f_vals = dist_3param.pdf(P_vals)
        assert np.all(f_vals > 0), "PDF must be positive for P > P_o"

    def test_pdf_integrates_to_unity(self, dist_3param):
        """Numerical integral of f(P) from P_o to large P ≈ 1."""
        from scipy.integrate import quad
        P_lo = dist_3param.P_o + 0.01
        P_hi = 1500.0
        integral, _ = quad(dist_3param.pdf, P_lo, P_hi)
        assert abs(integral - 1.0) < 0.01, (
            f"PDF integral over [{P_lo:.1f}, {P_hi}] = {integral:.5f}; expected ~1"
        )

    def test_average_above_P_o(self, dist_3param):
        avg = dist_3param.average()
        assert avg > dist_3param.P_o, (
            f"average {avg:.2f} K must exceed P_o = {dist_3param.P_o:.2f} K"
        )

    def test_average_plausible_range(self, dist_3param):
        """For this gas condensate, mean Tb should be in the 400-500 K range."""
        avg = dist_3param.average()
        assert 380.0 <= avg <= 550.0, (
            f"average Tb = {avg:.2f} K; expected in [380, 550] K for this condensate"
        )

    def test_scalar_return_types(self, dist_3param):
        assert isinstance(dist_3param.P(0.5),       float)
        assert isinstance(dist_3param.x_c(450.0),   float)
        assert isinstance(dist_3param.pdf(450.0),   float)
        assert isinstance(dist_3param.average(),    float)


# ── Input validation ──────────────────────────────────────────────────────────

class TestInputValidation:
    def test_raises_if_not_fitted(self):
        dist = GeneralizedDistribution()
        with pytest.raises(RuntimeError, match="not been fitted"):
            dist.P(0.5)

    def test_raises_if_not_fitted_xc(self):
        dist = GeneralizedDistribution()
        with pytest.raises(RuntimeError, match="not been fitted"):
            dist.x_c(400.0)

    def test_raises_if_not_fitted_pdf(self):
        dist = GeneralizedDistribution()
        with pytest.raises(RuntimeError, match="not been fitted"):
            dist.pdf(400.0)

    def test_raises_if_not_fitted_average(self):
        dist = GeneralizedDistribution()
        with pytest.raises(RuntimeError, match="not been fitted"):
            dist.average()

    def test_raises_for_xc_zero(self, tb_data):
        xc, Tb = tb_data
        bad_xc = np.concatenate([[0.0], xc[1:]])
        with pytest.raises(ValueError, match="open interval"):
            GeneralizedDistribution().fit(bad_xc, Tb)

    def test_raises_for_xc_one(self, tb_data):
        xc, Tb = tb_data
        bad_xc = np.concatenate([xc[:-1], [1.0]])
        with pytest.raises(ValueError, match="open interval"):
            GeneralizedDistribution().fit(bad_xc, Tb)

    def test_raises_for_mismatched_lengths(self, tb_data):
        xc, Tb = tb_data
        with pytest.raises(ValueError, match="equal length"):
            GeneralizedDistribution().fit(xc[:-1], Tb)

    def test_raises_for_invalid_mode(self, tb_data):
        xc, Tb = tb_data
        with pytest.raises(ValueError, match="mode"):
            GeneralizedDistribution().fit(xc, Tb, mode='bad')

    def test_raises_for_single_point(self, tb_data):
        xc, Tb = tb_data
        with pytest.raises(ValueError, match="2 data points"):
            GeneralizedDistribution().fit(xc[:1], Tb[:1])
