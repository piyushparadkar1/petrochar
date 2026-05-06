"""
Phase 1 pass-gate tests for core/correlations.py.

Pass-gate: riazi_daubert_Tb reproduces Riazi Table 4.11 within ±5 K for all
non-NaN rows (11 of 12; C18+ has no tabulated Tb).

Additional coverage:
  - Eq. 2.57 regime (M > 300) via Table 4.6 C22–C30 rows
  - riazi_daubert_M round-trip
  - watson_k spot-checks
  - aromaticity_gamma boundary values
  - gamma_function basic values and domain error

References
----------
Riazi MNL50 Table 4.11 (page 171): North Sea gas condensate C7+ example
Riazi MNL50 Table 4.6  (page 162): SCN group properties C6–C50
"""

import math
import pathlib
import warnings

import numpy as np
import pandas as pd
import pytest

from core.correlations import (
    aromaticity_gamma,
    gamma_function,
    riazi_daubert_M,
    riazi_daubert_SG,
    riazi_daubert_Tb,
    watson_k,
)

DATA_DIR = pathlib.Path(__file__).parent.parent / "data" / "riazi_reference"


# ── Table 4.11 — Eq. 2.56 regime pass-gate ────────────────────────────────────

class TestRiaziDaubertTbTable411:
    """Eq. 2.56 pass-gate: all 11 non-NaN rows from Riazi Table 4.11."""

    @pytest.fixture(scope="class")
    def table_411(self):
        df = pd.read_csv(DATA_DIR / "table_4_11_example_4_7.csv", comment="#")
        return df.dropna(subset=["Tb_K"])

    def test_row_count(self, table_411):
        assert len(table_411) == 11, (
            f"Expected 11 non-NaN Tb rows; got {len(table_411)}"
        )

    def test_all_rows_eq256_regime(self, table_411):
        """All Table 4.11 rows have M ≤ 300 so Eq. 2.56 is used."""
        assert (table_411["M_g_per_mol"] <= 300).all()

    @pytest.mark.parametrize("label,M,Tb_ref,SG", [
        ("C7",  95,  365, 0.727),
        ("C8",  107, 390, 0.749),
        ("C9",  121, 416, 0.768),
        ("C10", 136, 440, 0.782),
        ("C11", 149, 461, 0.793),
        ("C12", 163, 482, 0.804),
        ("C13", 176, 500, 0.815),
        ("C14", 191, 520, 0.826),
        ("C15", 207, 539, 0.836),
        ("C16", 221, 556, 0.843),
        ("C17", 237, 573, 0.851),
    ])
    def test_tb_within_5K(self, label, M, Tb_ref, SG):
        """Each Table 4.11 row: |Tb_computed − Tb_tabulated| ≤ 5 K."""
        Tb_calc = riazi_daubert_Tb(M, SG)
        delta = abs(Tb_calc - Tb_ref)
        assert delta <= 5.0, (
            f"{label}: Tb_calc={Tb_calc:.2f} K, Tb_ref={Tb_ref} K, "
            f"|Δ|={delta:.2f} K > 5.0 K tolerance"
        )

    def test_max_deviation_report(self, table_411, capsys):
        """Print actual deviations for all rows (informational)."""
        rows = []
        for _, row in table_411.iterrows():
            M, Tb_ref, SG = row["M_g_per_mol"], row["Tb_K"], row["SG_dimensionless"]
            Tb_calc = riazi_daubert_Tb(M, SG)
            delta = Tb_calc - Tb_ref
            rows.append((row["Carbon_No_label"], M, Tb_ref, Tb_calc, delta))

        print("\nTable 4.11 Eq. 2.56 deviations:")
        print(f"{'Label':<8} {'M':>6} {'Tb_ref':>8} {'Tb_calc':>8} {'delta':>8}")
        for label, M, Tb_ref, Tb_calc, delta in rows:
            print(f"{label:<8} {M:>6.0f} {Tb_ref:>8.1f} {Tb_calc:>8.2f} {delta:>+8.2f}")

        max_abs = max(abs(d) for *_, d in rows)
        print(f"\nMax |delta| = {max_abs:.2f} K  (spec target <= 2 K, test gate <= 5 K)")
        assert max_abs <= 5.0


# ── Table 4.6 — Eq. 2.57 regime ──────────────────────────────────────────────

class TestRiaziDaubertTbEq257:
    """Eq. 2.57 regime (M > 300) via Table 4.6 SCN rows C22–C30."""

    @pytest.mark.parametrize("cn,M,Tb_ref,SG", [
        (22, 303, 637, 0.876),
        (24, 331, 660, 0.885),
        (26, 359, 681, 0.892),
        (28, 387, 701, 0.899),
        (30, 415, 720, 0.905),
    ])
    def test_tb_eq257_within_25K(self, cn, M, Tb_ref, SG):
        """Eq. 2.57 rows: |Δ| ≤ 25 K (Riazi states ~4.7% AAD; at Tb≈650 K that is ~30 K)."""
        Tb_calc = riazi_daubert_Tb(M, SG)
        delta = abs(Tb_calc - Tb_ref)
        assert delta <= 25.0, (
            f"C{cn}: Tb_calc={Tb_calc:.2f} K, Tb_ref={Tb_ref} K, "
            f"|Δ|={delta:.2f} K > 25 K tolerance"
        )

    def test_eq257_no_warning_for_valid_range(self):
        """No UserWarning for M = 350 (within Eq. 2.57 validated range 300–700)."""
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            riazi_daubert_Tb(350.0, 0.88)

    def test_eq257_warns_above_700(self):
        """UserWarning raised for M > 700 (extrapolation)."""
        with pytest.warns(UserWarning, match="extrapolation"):
            riazi_daubert_Tb(750.0, 0.95)


# ── Boundary / validation tests ───────────────────────────────────────────────

class TestRiaziDaubertTbBoundaries:
    def test_raises_below_70(self):
        with pytest.raises(ValueError, match="70"):
            riazi_daubert_Tb(50.0, 0.70)

    def test_boundary_at_70(self):
        Tb = riazi_daubert_Tb(70.0, 0.70)
        assert 200 < Tb < 450

    def test_boundary_at_300_continuity(self):
        """Eq. 2.56 and 2.57 should not have a catastrophic discontinuity at M=300."""
        Tb_56 = riazi_daubert_Tb(299.9, 0.87)
        Tb_57 = riazi_daubert_Tb(300.1, 0.87)
        assert abs(Tb_56 - Tb_57) < 30.0, (
            f"Discontinuity at M=300: Tb(299.9)={Tb_56:.2f}, Tb(300.1)={Tb_57:.2f}"
        )


# ── riazi_daubert_M round-trip ────────────────────────────────────────────────

class TestRiaziDaubertMRoundTrip:
    @pytest.mark.parametrize("M_in,SG", [
        (95,  0.727),
        (150, 0.80),
        (250, 0.86),
        (310, 0.875),
        (400, 0.90),
    ])
    def test_round_trip(self, M_in, SG):
        """riazi_daubert_M(riazi_daubert_Tb(M, SG), SG) ≈ M within 1 g/mol."""
        Tb = riazi_daubert_Tb(M_in, SG)
        M_out = riazi_daubert_M(Tb, SG)
        assert abs(M_out - M_in) < 1.0, (
            f"Round-trip failed: M_in={M_in}, Tb={Tb:.2f}, M_out={M_out:.2f}"
        )


# ── riazi_daubert_SG round-trip ───────────────────────────────────────────────

class TestRiaziDaubertSGRoundTrip:
    @pytest.mark.parametrize("M,SG_in", [
        (95,  0.727),
        (200, 0.855),
        (310, 0.875),
    ])
    def test_round_trip(self, M, SG_in):
        """riazi_daubert_SG(M, riazi_daubert_Tb(M, SG)) ≈ SG within 0.001."""
        Tb = riazi_daubert_Tb(M, SG_in)
        SG_out = riazi_daubert_SG(M, Tb)
        assert abs(SG_out - SG_in) < 0.001, (
            f"Round-trip failed: M={M}, SG_in={SG_in}, SG_out={SG_out:.5f}"
        )


# ── watson_k ──────────────────────────────────────────────────────────────────

class TestWatsonK:
    def test_c7_spot_check(self):
        """C7 SCN: Tb=365 K, SG=0.727. K_W should be in paraffinic-aromatic range."""
        Kw = watson_k(365.0, 0.727)
        assert 10.0 < Kw < 14.0, f"K_W={Kw:.3f} out of expected range (10–14)"

    def test_rankine_conversion(self):
        """K_W = (1.8 * Tb_K)^(1/3) / SG — verify against manual calculation."""
        Tb_K, SG = 500.0, 0.85
        expected = (1.8 * 500.0) ** (1.0 / 3.0) / 0.85
        assert abs(watson_k(Tb_K, SG) - expected) < 1e-10

    def test_higher_kw_for_paraffin_character(self):
        """Lower SG (more paraffinic) at same Tb should give higher K_W."""
        Kw_low_sg  = watson_k(500.0, 0.80)
        Kw_high_sg = watson_k(500.0, 0.90)
        assert Kw_low_sg > Kw_high_sg


# ── aromaticity_gamma ─────────────────────────────────────────────────────────

class TestAromaticityGamma:
    def test_paraffinic_anchor(self):
        assert aromaticity_gamma(13.0) == pytest.approx(0.0)

    def test_aromatic_anchor(self):
        assert aromaticity_gamma(9.5) == pytest.approx(1.0)

    def test_midpoint(self):
        assert aromaticity_gamma(11.25) == pytest.approx(0.5)

    def test_clamp_below_lower_kw(self):
        """K_W < 9.5 → γ clamped to 1.0."""
        assert aromaticity_gamma(8.0) == pytest.approx(1.0)

    def test_clamp_above_upper_kw(self):
        """K_W > 13.0 → γ clamped to 0.0."""
        assert aromaticity_gamma(15.0) == pytest.approx(0.0)

    def test_monotone_decreasing(self):
        """γ decreases as K_W increases from 9.5 to 13.0."""
        kw_vals = [9.5, 10.5, 11.25, 12.0, 13.0]
        gamma_vals = [aromaticity_gamma(k) for k in kw_vals]
        for i in range(len(gamma_vals) - 1):
            assert gamma_vals[i] >= gamma_vals[i + 1]


# ── gamma_function ────────────────────────────────────────────────────────────

class TestGammaFunction:
    def test_gamma_1(self):
        assert gamma_function(1.0) == pytest.approx(1.0, rel=1e-10)

    def test_gamma_2(self):
        assert gamma_function(2.0) == pytest.approx(1.0, rel=1e-10)

    def test_gamma_3(self):
        assert gamma_function(3.0) == pytest.approx(2.0, rel=1e-10)

    def test_gamma_half(self):
        assert gamma_function(0.5) == pytest.approx(math.sqrt(math.pi), rel=1e-10)

    def test_gamma_recurrence(self):
        """Γ(x+1) = x · Γ(x) for x = 2.7."""
        x = 2.7
        assert gamma_function(x + 1) == pytest.approx(x * gamma_function(x), rel=1e-10)

    def test_raises_for_zero(self):
        with pytest.raises(ValueError):
            gamma_function(0.0)

    def test_raises_for_negative(self):
        with pytest.raises(ValueError):
            gamma_function(-1.5)

    def test_positive_for_all_valid(self):
        for x in [0.01, 0.5, 1.0, 2.0, 5.0, 10.0]:
            assert gamma_function(x) > 0.0
