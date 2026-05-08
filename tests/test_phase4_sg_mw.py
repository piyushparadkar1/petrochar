"""
Phase 4 pass-gate: SG and MW distributions from Riazi Table 4.11 data.

Pass-gate criteria:
  1. SGDistributionFit (3-param, volume basis): SG_o=0.705, A=0.0232, B=1.811
     within ±5% of Riazi Table 4.13 values.
  2. Bulk M_av (riazi_daubert_M for C7-C17, tabulated M for C18+) within 1%
     of Riazi Example 4.7 measured value M_7+ = 118.9 g/mol.
  3. Bulk SG_av (density-additive, tabulated SG) within 1% of SG_7+ = 0.7597.

Additional coverage:
  - sg_from_watson_k: per-component SG within ±0.010 of tabulated SG.
  - Watson K SG bulk closure within 1% of 0.7597.
  - compute_M_array: M values computed without NaN inputs.
  - Input validation for all public functions.

References
----------
Riazi MNL50 Table 4.11 (page 171) — input data.
Riazi MNL50 Table 4.13 (page 173) — expected SG distribution coefficients.
Riazi MNL50 §4.5.4.4 Eqs. 4.74-4.76 — bulk property closure formulas.
"""

import pathlib
import warnings

import numpy as np
import pandas as pd
import pytest

from core.sg_distribution import SGDistributionFit, sg_from_watson_k
from core.mw_distribution import bulk_M_av, bulk_SG_av, compute_M_array

DATA_DIR = pathlib.Path(__file__).parent.parent / "data" / "riazi_reference"

_M_BULK  = 118.9    # Riazi Example 4.7, g/mol
_SG_BULK = 0.7597   # Riazi Example 4.7


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def table_411():
    return pd.read_csv(DATA_DIR / "table_4_11_example_4_7.csv", comment="#")


@pytest.fixture(scope="module")
def table_411_with_tb(table_411):
    return table_411.dropna(subset=["Tb_K"])


@pytest.fixture(scope="module")
def sg_fit_3param(table_411):
    xc  = table_411["x_cv"].values.astype(float)
    SG  = table_411["SG_dimensionless"].values.astype(float)
    return SGDistributionFit().fit(xc, SG, mode="3param")


@pytest.fixture(scope="module")
def sg_fit_2param(table_411):
    xc  = table_411["x_cv"].values.astype(float)
    SG  = table_411["SG_dimensionless"].values.astype(float)
    return SGDistributionFit().fit(xc, SG, mode="2param", B_fixed=3.0)


@pytest.fixture(scope="module")
def watson_k_SG(table_411_with_tb):
    Tb  = table_411_with_tb["Tb_K"].values.astype(float)
    SG, K_W = sg_from_watson_k(Tb, _M_BULK, _SG_BULK)
    return SG, K_W


# ── SG generalized 3-param pass-gate ─────────────────────────────────────────

class TestSGGeneralizedFitPassGate:
    """Pass-gate: reproduce Table 4.13 SG row (volume basis, 3-param)."""

    _SG_o_ref = 0.705
    _A_ref    = 0.0232
    _B_ref    = 1.811
    _TOL      = 0.05

    def test_SG_o_within_5pct(self, sg_fit_3param):
        err = abs(sg_fit_3param.SG_o - self._SG_o_ref) / self._SG_o_ref
        assert err <= self._TOL, (
            f"SG_o={sg_fit_3param.SG_o:.5f} vs ref {self._SG_o_ref}; "
            f"error {err*100:.2f}% > 5%"
        )

    def test_A_within_5pct(self, sg_fit_3param):
        err = abs(sg_fit_3param.A - self._A_ref) / self._A_ref
        assert err <= self._TOL, (
            f"A={sg_fit_3param.A:.6f} vs ref {self._A_ref}; "
            f"error {err*100:.2f}% > 5%"
        )

    def test_B_within_5pct(self, sg_fit_3param):
        err = abs(sg_fit_3param.B - self._B_ref) / self._B_ref
        assert err <= self._TOL, (
            f"B={sg_fit_3param.B:.5f} vs ref {self._B_ref}; "
            f"error {err*100:.2f}% > 5%"
        )

    def test_R2_above_0p98(self, sg_fit_3param):
        assert sg_fit_3param.fit_quality["R_squared"] >= 0.98

    def test_print_params(self, sg_fit_3param, capsys):
        p  = sg_fit_3param.params
        fq = sg_fit_3param.fit_quality
        print("\nSG 3-param fit (Table 4.13 ref: SG_o=0.705, A=0.0232, B=1.811):")
        print(f"  SG_o = {p['P_o']:.5f}  (ref 0.705)   err {abs(p['P_o']-0.705)/0.705*100:.2f}%")
        print(f"  A    = {p['A']:.6f}  (ref 0.0232)  err {abs(p['A']-0.0232)/0.0232*100:.2f}%")
        print(f"  B    = {p['B']:.5f}  (ref 1.811)   err {abs(p['B']-1.811)/1.811*100:.2f}%")
        print(f"  RMS  = {fq['RMS']:.5f}   %AAD = {fq['pct_AAD']:.4f}%   R2 = {fq['R_squared']:.6f}")


# ── SG generalized 2-param pass-gate ─────────────────────────────────────────

class TestSGGeneralizedTwoParam:
    """Two-parameter SG fit (B_SG=3 fixed per Riazi p.174).

    Table 4.13 reference (volume basis, 2-param): SG_o=0.665, A=0.0132, B=3.
    """

    _SG_o_ref = 0.665
    _A_ref    = 0.0132
    _TOL      = 0.05

    def test_B_fixed_at_3(self, sg_fit_2param):
        assert abs(sg_fit_2param.B - 3.0) < 1e-10

    def test_SG_o_within_5pct(self, sg_fit_2param):
        err = abs(sg_fit_2param.SG_o - self._SG_o_ref) / self._SG_o_ref
        assert err <= self._TOL, (
            f"SG_o={sg_fit_2param.SG_o:.5f} vs ref {self._SG_o_ref}; "
            f"error {err*100:.2f}% > 5%"
        )

    def test_A_within_5pct(self, sg_fit_2param):
        err = abs(sg_fit_2param.A - self._A_ref) / self._A_ref
        assert err <= self._TOL, (
            f"A={sg_fit_2param.A:.6f} vs ref {self._A_ref}; "
            f"error {err*100:.2f}% > 5%"
        )


# ── Constant Watson K method ──────────────────────────────────────────────────

class TestWatsonKSGMethod:
    """sg_from_watson_k: per-component SG close to Table 4.11 tabulated values."""

    def test_K_W_plausible_range(self, watson_k_SG):
        _, K_W = watson_k_SG
        assert 10.0 < K_W < 14.0, f"K_W = {K_W:.4f} outside expected [10, 14]"

    def test_SG_array_increasing(self, watson_k_SG, table_411_with_tb):
        SG, _ = watson_k_SG
        assert np.all(np.diff(SG) > 0), "Watson K SG must increase with Tb"

    def test_SG_close_to_tabulated(self, watson_k_SG, table_411_with_tb):
        """Per-component SG within ±0.010 of Table 4.11 tabulated values."""
        SG_computed, _ = watson_k_SG
        SG_tab = table_411_with_tb["SG_dimensionless"].values.astype(float)
        labels = table_411_with_tb["Carbon_No_label"].tolist()
        tol = 0.010
        for lbl, sg_c, sg_t in zip(labels, SG_computed, SG_tab):
            assert abs(sg_c - sg_t) <= tol, (
                f"{lbl}: SG_watson_k={sg_c:.4f}, SG_tab={sg_t:.4f}, "
                f"|diff|={abs(sg_c-sg_t):.4f} > {tol}"
            )

    def test_SG_array_length(self, watson_k_SG, table_411_with_tb):
        SG, _ = watson_k_SG
        assert len(SG) == len(table_411_with_tb)

    def test_bulk_SG_watson_k_within_1pct(self, watson_k_SG, table_411):
        """Bulk SG_av from Watson K SG (C7-C17) + tabulated C18+ within 1% of 0.7597."""
        SG_11, _ = watson_k_SG
        xw_11  = table_411.dropna(subset=["Tb_K"])["x_w"].values.astype(float)
        xw_12  = float(table_411.iloc[11]["x_w"])
        SG_12  = float(table_411.iloc[11]["SG_dimensionless"])
        SG_all = np.append(SG_11, SG_12)
        xw_all = np.append(xw_11, xw_12)
        SG_av  = bulk_SG_av(xw_all, SG_all)
        err = abs(SG_av - _SG_BULK) / _SG_BULK
        assert err < 0.01, (
            f"Bulk SG_av (Watson K) = {SG_av:.5f}, target {_SG_BULK}; "
            f"error {err*100:.2f}% > 1%"
        )


# ── Bulk closure pass-gate ────────────────────────────────────────────────────

class TestBulkClosurePassGate:
    """Pass-gate: M_av and SG_av within 1% of Riazi Example 4.7 bulk values."""

    def test_bulk_SG_av_tabulated_within_1pct(self, table_411):
        """Eq. 4.76: SG_av = 1/sum(x_w/SG) from tabulated data within 1% of 0.7597."""
        x_w = table_411["x_w"].values.astype(float)
        SG  = table_411["SG_dimensionless"].values.astype(float)
        SG_av = bulk_SG_av(x_w, SG)
        err = abs(SG_av - _SG_BULK) / _SG_BULK
        assert err < 0.01, (
            f"Bulk SG_av = {SG_av:.5f}, target {_SG_BULK}; error {err*100:.2f}% > 1%"
        )

    def test_bulk_M_av_computed_within_1pct(self, table_411, table_411_with_tb):
        """Eq. 4.74: M_av = 1/sum(x_w/M) within 1% of 118.9 g/mol.

        M_i computed via riazi_daubert_M for C7-C17 (11 rows with Tb);
        tabulated M=264 used for C18+ (no Tb available).
        """
        Tb11 = table_411_with_tb["Tb_K"].values.astype(float)
        SG11 = table_411_with_tb["SG_dimensionless"].values.astype(float)
        xw11 = table_411_with_tb["x_w"].values.astype(float)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            M11 = compute_M_array(Tb11, SG11)

        xw12 = float(table_411.iloc[11]["x_w"])
        M12  = float(table_411.iloc[11]["M_g_per_mol"])

        x_w_all = np.append(xw11, xw12)
        M_all   = np.append(M11, M12)
        M_av    = bulk_M_av(x_w_all, M_all)

        err = abs(M_av - _M_BULK) / _M_BULK
        assert err < 0.01, (
            f"Bulk M_av = {M_av:.3f} g/mol, target {_M_BULK}; error {err*100:.2f}% > 1%"
        )

    def test_print_bulk_values(self, table_411, table_411_with_tb, capsys):
        x_w  = table_411["x_w"].values.astype(float)
        SG   = table_411["SG_dimensionless"].values.astype(float)
        SG_av = bulk_SG_av(x_w, SG)

        Tb11 = table_411_with_tb["Tb_K"].values.astype(float)
        SG11 = table_411_with_tb["SG_dimensionless"].values.astype(float)
        xw11 = table_411_with_tb["x_w"].values.astype(float)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            M11 = compute_M_array(Tb11, SG11)
        xw12 = float(table_411.iloc[11]["x_w"])
        M12  = float(table_411.iloc[11]["M_g_per_mol"])
        M_av = bulk_M_av(np.append(xw11, xw12), np.append(M11, M12))

        print("\nBulk closure (Eqs. 4.74-4.76 vs Example 4.7 measured values):")
        print(f"  M_av  = {M_av:.3f} g/mol  (ref 118.9)   err {abs(M_av-118.9)/118.9*100:.2f}%")
        print(f"  SG_av = {SG_av:.5f}         (ref 0.7597)  err {abs(SG_av-0.7597)/0.7597*100:.2f}%")


# ── compute_M_array ───────────────────────────────────────────────────────────

class TestComputeMArray:
    def test_output_length(self, table_411_with_tb):
        Tb  = table_411_with_tb["Tb_K"].values.astype(float)
        SG  = table_411_with_tb["SG_dimensionless"].values.astype(float)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            M = compute_M_array(Tb, SG)
        assert len(M) == len(Tb)

    def test_M_increasing_with_Tb(self, table_411_with_tb):
        Tb  = table_411_with_tb["Tb_K"].values.astype(float)
        SG  = table_411_with_tb["SG_dimensionless"].values.astype(float)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            M = compute_M_array(Tb, SG)
        assert np.all(np.diff(M) > 0), "M must increase with Tb for this condensate"

    def test_M_range_plausible(self, table_411_with_tb):
        Tb  = table_411_with_tb["Tb_K"].values.astype(float)
        SG  = table_411_with_tb["SG_dimensionless"].values.astype(float)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            M = compute_M_array(Tb, SG)
        assert M.min() > 70, f"Minimum M = {M.min():.1f} below valid range"
        assert M.max() < 500, f"Maximum M = {M.max():.1f} unexpectedly high"

    def test_raises_on_nan_Tb(self, table_411_with_tb):
        Tb  = table_411_with_tb["Tb_K"].values.astype(float)
        SG  = table_411_with_tb["SG_dimensionless"].values.astype(float)
        Tb_nan = np.append(Tb, np.nan)
        SG_nan = np.append(SG, 0.86)
        with pytest.raises(ValueError, match="NaN"):
            compute_M_array(Tb_nan, SG_nan)


# ── Evaluation methods of SGDistributionFit ───────────────────────────────────

class TestSGDistributionEvaluation:
    def test_SG_increases_with_xc(self, sg_fit_3param):
        xc = np.linspace(0.05, 0.95, 20)
        sg = sg_fit_3param.SG(xc)
        assert np.all(np.diff(sg) > 0)

    def test_xc_roundtrip(self, sg_fit_3param):
        xc_in  = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
        SG_val = sg_fit_3param.SG(xc_in)
        xc_out = sg_fit_3param.xc(SG_val)
        np.testing.assert_allclose(xc_out, xc_in, rtol=1e-8)

    def test_average_SG_above_SG_o(self, sg_fit_3param):
        avg = sg_fit_3param.average_SG()
        assert avg > sg_fit_3param.SG_o

    def test_average_SG_plausible(self, sg_fit_3param):
        avg = sg_fit_3param.average_SG()
        assert 0.75 < avg < 0.85, f"average_SG = {avg:.4f} outside expected [0.75, 0.85]"

    def test_scalar_returns(self, sg_fit_3param):
        assert isinstance(sg_fit_3param.SG(0.5), float)
        assert isinstance(sg_fit_3param.xc(0.80), float)
        assert isinstance(sg_fit_3param.average_SG(), float)


# ── Input validation ──────────────────────────────────────────────────────────

class TestInputValidation:
    def test_sg_fit_raises_if_not_fitted(self):
        dist = SGDistributionFit()
        with pytest.raises(RuntimeError, match="not been fitted"):
            dist.SG(0.5)

    def test_sg_from_watson_k_raises_negative_M(self, table_411_with_tb):
        Tb = table_411_with_tb["Tb_K"].values.astype(float)
        with pytest.raises(Exception):
            sg_from_watson_k(Tb, M_bulk=-1.0, SG_bulk=0.76)

    def test_bulk_M_av_raises_mismatched_shapes(self):
        with pytest.raises(ValueError, match="same shape"):
            bulk_M_av([0.5, 0.5], [100.0, 200.0, 300.0])

    def test_bulk_SG_av_raises_mismatched_shapes(self):
        with pytest.raises(ValueError, match="same shape"):
            bulk_SG_av([0.5, 0.5], [0.8])

    def test_bulk_M_av_raises_nonpositive_M(self):
        with pytest.raises(ValueError, match="positive"):
            bulk_M_av([1.0], [-100.0])

    def test_bulk_SG_av_raises_nonpositive_SG(self):
        with pytest.raises(ValueError, match="positive"):
            bulk_SG_av([1.0], [0.0])

    def test_compute_M_raises_shape_mismatch(self):
        with pytest.raises(ValueError, match="equal length"):
            compute_M_array([365.0, 390.0], [0.727])
