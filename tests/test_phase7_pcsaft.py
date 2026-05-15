"""
Phase 7 pass-gate: PC-SAFT parameter assignment (Panuganti 2012 Table 6).

Numerical pass criteria (per user Phase 7 instruction):
  m        : within ±1%  of Panuganti Tables 10-12 published values
  sigma    : within ±0.5% of published values
  eps/k    : within ±2%  of published values

Parametrized over Saturates and Aromatics+Resins rows from Crudes A, B, C
loaded from data/panuganti_2012/tables_10_11_12_crude_pcsaft.csv.

Additional tests:
  - propane sigma is 3.6180 Å (explicit guard against the 3.168 Aspen typo)
  - propane/asphaltene params return exact fixed values
  - panuganti_distillable_params routes to A+R form (not Saturates form)
  - monotonicity: m increases with MW at fixed gamma; sigma increases with MW
  - gamma endpoint exactness: gamma=0 → benzene-derivative form, gamma=1 → PNA form
  - compute_K_W_per_pseudocomponent assigns K_W and gamma; immutable input
  - generate_pcsaft_table returns correct columns and dispatches correctly

References
----------
Panuganti et al. (2012) Fuel 93, 658-669, Tables 6, 10, 11, 12.
Gonzalez et al. (2007) Energy & Fuels 21, 1230-1234.
Gross J. and Sadowski G. (2001) Ind. Eng. Chem. Res. 40, 1244.
"""

from __future__ import annotations

import math
import pathlib

import numpy as np
import pandas as pd
import pytest

from core.pcsaft_params import (
    generate_pcsaft_table,
    gonzalez_asphaltene_params,
    panuganti_aromatic_resin_params,
    panuganti_distillable_params,
    panuganti_saturate_params,
    propane_params,
)
from core.quadrature import Pseudocomponent
from core.watson_k import compute_K_W_per_pseudocomponent


# ── Load Panuganti reference data ─────────────────────────────────────────────

_DATA_DIR = pathlib.Path(__file__).parent.parent / 'data' / 'panuganti_2012'


def _load_crude_csv() -> pd.DataFrame:
    path = _DATA_DIR / 'tables_10_11_12_crude_pcsaft.csv'
    return pd.read_csv(path, comment='#')


_DF = _load_crude_csv()

# Saturates rows: (crude_label, MW, m_ref, sigma_ref, epsk_ref)
_SAT_ROWS = [
    (
        str(row['Crude']),
        float(row['MW_g_per_mol']),
        float(row['m']),
        float(row['sigma_Angstrom']),
        float(row['epsilon_over_k_K']),
    )
    for _, row in _DF[_DF['Component'] == 'Saturates'].iterrows()
]

# Aromatics+Resins rows: (crude_label, MW, gamma, m_ref, sigma_ref, epsk_ref)
_AR_ROWS = [
    (
        str(row['Crude']),
        float(row['MW_g_per_mol']),
        float(row['gamma']),
        float(row['m']),
        float(row['sigma_Angstrom']),
        float(row['epsilon_over_k_K']),
    )
    for _, row in _DF[_DF['Component'] == 'Aromatics_plus_Resins'].iterrows()
]


# ── Helper: build minimal Pseudocomponent with Tb_K and SG populated ─────────

def _make_pc(z: float, M: float, Tb_K: float, SG: float) -> Pseudocomponent:
    """Return a Pseudocomponent with Tb_K and SG set (K_W/gamma still None)."""
    return Pseudocomponent(z=z, M=M, Tb_K=Tb_K, SG=SG)


# ── Propane params ────────────────────────────────────────────────────────────

class TestPropaneParams:
    def test_returns_tuple_of_three(self):
        result = propane_params()
        assert len(result) == 3

    def test_m_exact(self):
        m, _, _ = propane_params()
        assert m == pytest.approx(2.002, abs=0.0)

    def test_sigma_exact_3618_not_3168(self):
        # 3.618 Å is correct (Gross & Sadowski 2001).
        # 3.168 is the Aspen typo — a 14% error.  This test guards against it.
        _, sigma, _ = propane_params()
        assert sigma == pytest.approx(3.6180, abs=1e-9), (
            f"sigma={sigma:.4f}: should be 3.6180, not 3.168 (Aspen typo)"
        )

    def test_sigma_is_not_aspen_typo(self):
        _, sigma, _ = propane_params()
        assert abs(sigma - 3.168) > 0.1, (
            "sigma=3.168 is the Aspen built-in database typo. Correct value is 3.6180."
        )

    def test_epsk_exact(self):
        _, _, epsk = propane_params()
        assert epsk == pytest.approx(208.11, abs=0.0)


# ── Gonzalez 2007 asphaltene defaults ─────────────────────────────────────────

class TestGonzalezAsphalteneParams:
    def test_returns_tuple_of_three(self):
        result = gonzalez_asphaltene_params()
        assert len(result) == 3

    def test_m_exact(self):
        m, _, _ = gonzalez_asphaltene_params()
        assert m == pytest.approx(33.0, abs=0.0)

    def test_sigma_exact(self):
        _, sigma, _ = gonzalez_asphaltene_params()
        assert sigma == pytest.approx(4.3, abs=0.0)

    def test_epsk_exact(self):
        _, _, epsk = gonzalez_asphaltene_params()
        assert epsk == pytest.approx(400.0, abs=0.0)


# ── Saturates correlation (Table 6 rows 1-3) ──────────────────────────────────

class TestPanuganti_SaturateParams:
    @pytest.mark.parametrize(
        "crude,MW,m_ref,sigma_ref,epsk_ref",
        _SAT_ROWS,
        ids=[f"Crude{r[0]}" for r in _SAT_ROWS],
    )
    def test_m_within_1pct(self, crude, MW, m_ref, sigma_ref, epsk_ref):
        m, _, _ = panuganti_saturate_params(MW)
        pct = abs(m - m_ref) / m_ref * 100.0
        assert pct <= 1.0, (
            f"Crude {crude}: m={m:.4f}, ref={m_ref:.3f}, dev={pct:.3f}% > 1%"
        )

    @pytest.mark.parametrize(
        "crude,MW,m_ref,sigma_ref,epsk_ref",
        _SAT_ROWS,
        ids=[f"Crude{r[0]}" for r in _SAT_ROWS],
    )
    def test_sigma_within_05pct(self, crude, MW, m_ref, sigma_ref, epsk_ref):
        _, sigma, _ = panuganti_saturate_params(MW)
        pct = abs(sigma - sigma_ref) / sigma_ref * 100.0
        assert pct <= 0.5, (
            f"Crude {crude}: sigma={sigma:.4f}, ref={sigma_ref:.3f}, dev={pct:.3f}% > 0.5%"
        )

    @pytest.mark.parametrize(
        "crude,MW,m_ref,sigma_ref,epsk_ref",
        _SAT_ROWS,
        ids=[f"Crude{r[0]}" for r in _SAT_ROWS],
    )
    def test_epsk_within_2pct(self, crude, MW, m_ref, sigma_ref, epsk_ref):
        _, _, epsk = panuganti_saturate_params(MW)
        pct = abs(epsk - epsk_ref) / epsk_ref * 100.0
        assert pct <= 2.0, (
            f"Crude {crude}: eps/k={epsk:.2f}, ref={epsk_ref:.2f}, dev={pct:.3f}% > 2%"
        )

    def test_m_increases_with_MW(self):
        # m = 0.0257*MW + 0.8444 is strictly linear increasing in MW.
        m100, _, _ = panuganti_saturate_params(100.0)
        m300, _, _ = panuganti_saturate_params(300.0)
        m700, _, _ = panuganti_saturate_params(700.0)
        assert m100 < m300 < m700

    def test_sigma_positive(self):
        for MW in [100.0, 200.0, 400.0, 700.0]:
            _, sigma, _ = panuganti_saturate_params(MW)
            assert sigma > 0.0, f"MW={MW}: sigma={sigma:.4f} <= 0"

    def test_epsk_positive(self):
        for MW in [100.0, 200.0, 400.0, 700.0]:
            _, _, epsk = panuganti_saturate_params(MW)
            assert epsk > 0.0, f"MW={MW}: eps/k={epsk:.2f} <= 0"

    def test_returns_floats(self):
        result = panuganti_saturate_params(200.0)
        assert all(isinstance(v, float) for v in result)


# ── Aromatics+Resins correlation (Table 6 rows 4-9) ──────────────────────────

class TestPanuganti_AromaticResinParams:
    @pytest.mark.parametrize(
        "crude,MW,gamma,m_ref,sigma_ref,epsk_ref",
        _AR_ROWS,
        ids=[f"Crude{r[0]}" for r in _AR_ROWS],
    )
    def test_m_within_1pct(self, crude, MW, gamma, m_ref, sigma_ref, epsk_ref):
        m, _, _ = panuganti_aromatic_resin_params(MW, gamma)
        pct = abs(m - m_ref) / m_ref * 100.0
        assert pct <= 1.0, (
            f"Crude {crude} (gamma={gamma}): m={m:.4f}, ref={m_ref:.3f}, dev={pct:.3f}% > 1%"
        )

    @pytest.mark.parametrize(
        "crude,MW,gamma,m_ref,sigma_ref,epsk_ref",
        _AR_ROWS,
        ids=[f"Crude{r[0]}" for r in _AR_ROWS],
    )
    def test_sigma_within_05pct(self, crude, MW, gamma, m_ref, sigma_ref, epsk_ref):
        _, sigma, _ = panuganti_aromatic_resin_params(MW, gamma)
        pct = abs(sigma - sigma_ref) / sigma_ref * 100.0
        assert pct <= 0.5, (
            f"Crude {crude} (gamma={gamma}): sigma={sigma:.4f}, ref={sigma_ref:.3f}, "
            f"dev={pct:.3f}% > 0.5%"
        )

    @pytest.mark.parametrize(
        "crude,MW,gamma,m_ref,sigma_ref,epsk_ref",
        _AR_ROWS,
        ids=[f"Crude{r[0]}" for r in _AR_ROWS],
    )
    def test_epsk_within_2pct(self, crude, MW, gamma, m_ref, sigma_ref, epsk_ref):
        _, _, epsk = panuganti_aromatic_resin_params(MW, gamma)
        pct = abs(epsk - epsk_ref) / epsk_ref * 100.0
        assert pct <= 2.0, (
            f"Crude {crude} (gamma={gamma}): eps/k={epsk:.2f}, ref={epsk_ref:.2f}, "
            f"dev={pct:.3f}% > 2%"
        )

    def test_gamma_zero_equals_benzene_derivative_form(self):
        # At gamma=0, all three expressions collapse to the benzene-derivative end-member.
        MW = 250.0
        m, sigma, epsk = panuganti_aromatic_resin_params(MW, 0.0)
        assert m     == pytest.approx(0.0223 * MW + 0.751,          rel=1e-9)
        assert sigma == pytest.approx(4.1377 - 38.1483 / MW,         rel=1e-9)
        assert epsk  == pytest.approx(0.00436 * MW + 283.93,          rel=1e-9)

    def test_gamma_one_equals_pna_form(self):
        # At gamma=1, all three expressions collapse to the PNA end-member.
        MW = 250.0
        m, sigma, epsk = panuganti_aromatic_resin_params(MW, 1.0)
        assert m     == pytest.approx(0.0101 * MW + 1.7296,           rel=1e-9)
        assert sigma == pytest.approx(4.6169 - 93.98 / MW,            rel=1e-9)
        assert epsk  == pytest.approx(508.0 - 234100.0 / MW ** 1.5,   rel=1e-9)

    def test_linear_interpolation_midpoint(self):
        # At gamma=0.5, result must be exact average of the two end-members.
        MW = 300.0
        m0, s0, e0 = panuganti_aromatic_resin_params(MW, 0.0)
        m1, s1, e1 = panuganti_aromatic_resin_params(MW, 1.0)
        mh, sh, eh = panuganti_aromatic_resin_params(MW, 0.5)
        assert mh == pytest.approx(0.5 * m0 + 0.5 * m1, rel=1e-9)
        assert sh == pytest.approx(0.5 * s0 + 0.5 * s1, rel=1e-9)
        assert eh == pytest.approx(0.5 * e0 + 0.5 * e1, rel=1e-9)

    def test_m_increases_with_MW_at_fixed_gamma(self):
        # Both end-members have a positive MW coefficient for m.
        for gamma in [0.0, 0.5, 1.0]:
            m100, _, _ = panuganti_aromatic_resin_params(100.0, gamma)
            m300, _, _ = panuganti_aromatic_resin_params(300.0, gamma)
            m700, _, _ = panuganti_aromatic_resin_params(700.0, gamma)
            assert m100 < m300 < m700, (
                f"gamma={gamma}: m not monotone in MW: {m100:.3f}, {m300:.3f}, {m700:.3f}"
            )

    def test_sigma_increases_with_MW_at_fixed_gamma(self):
        # Both end-members have the form: C1 - C2/MW, which is monotone increasing in MW.
        for gamma in [0.0, 0.5, 1.0]:
            _, s100, _ = panuganti_aromatic_resin_params(100.0, gamma)
            _, s300, _ = panuganti_aromatic_resin_params(300.0, gamma)
            _, s700, _ = panuganti_aromatic_resin_params(700.0, gamma)
            assert s100 < s300 < s700, (
                f"gamma={gamma}: sigma not monotone in MW: {s100:.4f}, {s300:.4f}, {s700:.4f}"
            )

    def test_returns_floats(self):
        result = panuganti_aromatic_resin_params(200.0, 0.3)
        assert all(isinstance(v, float) for v in result)


# ── panuganti_distillable_params routes to A+R (not Saturates) ────────────────

class TestPanuganti_DistillableParams:
    """Verify that panuganti_distillable_params returns the same result as
    panuganti_aromatic_resin_params, not panuganti_saturate_params.

    Design rationale: the gamma-interpolated A+R correlation already encodes
    the full paraffinic-to-PNA continuum; using the Saturates correlation in
    parallel would double-count the paraffinic character captured by gamma.
    (See core/pcsaft_params.py module docstring, Design note section.)
    """

    @pytest.mark.parametrize("MW,gamma", [
        (150.0, 0.0),
        (250.0, 0.1),
        (350.0, 0.5),
        (500.0, 1.0),
    ])
    def test_matches_aromatic_resin_params(self, MW, gamma):
        m_d, s_d, e_d = panuganti_distillable_params(MW, gamma)
        m_ar, s_ar, e_ar = panuganti_aromatic_resin_params(MW, gamma)
        assert m_d == pytest.approx(m_ar, rel=1e-12)
        assert s_d == pytest.approx(s_ar, rel=1e-12)
        assert e_d == pytest.approx(e_ar, rel=1e-12)

    @pytest.mark.parametrize("MW,gamma", [
        (150.0, 0.0),
        (250.0, 0.0),
        (350.0, 0.0),
    ])
    def test_does_not_match_saturate_params(self, MW, gamma):
        # The Saturates form uses slope 0.0257 for m; the A+R form uses 0.0223.
        # At gamma=0 (paraffinic end), m values from the two forms differ by >0.1
        # across the full MW range — a clear structural distinction.
        m_d, _, _ = panuganti_distillable_params(MW, gamma)
        m_s, _, _ = panuganti_saturate_params(MW)
        assert abs(m_d - m_s) > 0.1, (
            f"MW={MW}: distillable m {m_d:.4f} too close to saturate m {m_s:.4f} "
            f"— routing may be wrong."
        )


# ── compute_K_W_per_pseudocomponent ───────────────────────────────────────────

class TestComputeKWPerPseudocomponent:
    """Watson K and gamma assignment for pseudo-components."""

    @pytest.fixture
    def three_components(self):
        """Three distillable pseudo-components with Tb_K and SG set."""
        return [
            _make_pc(z=0.5, M=120.0, Tb_K=450.0, SG=0.82),
            _make_pc(z=0.3, M=200.0, Tb_K=600.0, SG=0.90),
            _make_pc(z=0.2, M=350.0, Tb_K=750.0, SG=1.00),
        ]

    def test_returns_new_list(self, three_components):
        result = compute_K_W_per_pseudocomponent(three_components)
        assert result is not three_components

    def test_same_length(self, three_components):
        result = compute_K_W_per_pseudocomponent(three_components)
        assert len(result) == len(three_components)

    def test_input_list_not_modified(self, three_components):
        original_kw = [c.K_W for c in three_components]
        compute_K_W_per_pseudocomponent(three_components)
        assert [c.K_W for c in three_components] == original_kw

    def test_kw_populated(self, three_components):
        result = compute_K_W_per_pseudocomponent(three_components)
        assert all(c.K_W is not None for c in result)

    def test_gamma_populated(self, three_components):
        result = compute_K_W_per_pseudocomponent(three_components)
        assert all(c.gamma is not None for c in result)

    def test_kw_positive(self, three_components):
        result = compute_K_W_per_pseudocomponent(three_components)
        assert all(c.K_W > 0 for c in result)

    def test_gamma_in_unit_interval(self, three_components):
        result = compute_K_W_per_pseudocomponent(three_components)
        assert all(0.0 <= c.gamma <= 1.0 for c in result), (
            f"gamma out of [0,1]: {[c.gamma for c in result]}"
        )

    def test_other_fields_preserved(self, three_components):
        result = compute_K_W_per_pseudocomponent(three_components)
        for orig, new in zip(three_components, result):
            assert new.z     == orig.z
            assert new.M     == orig.M
            assert new.Tb_K  == orig.Tb_K
            assert new.SG    == orig.SG

    def test_kw_numerical_value(self):
        # K_W = (1.8 * Tb_K)^(1/3) / SG (Watson & Nelson 1933).
        # For Tb_K=500, SG=1.0: K_W = (900)^(1/3)/1.0 = 9.6549...
        # gamma = clamp((13.0-9.6549)/3.5, 0, 1) = clamp(0.9557, 0, 1) = 0.9557
        from core.correlations import watson_k, aromaticity_gamma
        Tb_K, SG = 500.0, 1.0
        expected_kw    = watson_k(Tb_K, SG)
        expected_gamma = aromaticity_gamma(expected_kw)
        pc = _make_pc(z=1.0, M=200.0, Tb_K=Tb_K, SG=SG)
        result = compute_K_W_per_pseudocomponent([pc])
        assert result[0].K_W   == pytest.approx(expected_kw,    rel=1e-9)
        assert result[0].gamma == pytest.approx(expected_gamma,  rel=1e-9)

    def test_raises_on_nan_Tb_K(self):
        bad = _make_pc(z=1.0, M=200.0, Tb_K=float('nan'), SG=0.9)
        with pytest.raises(ValueError, match="nan Tb_K"):
            compute_K_W_per_pseudocomponent([bad])

    def test_raises_on_nan_SG(self):
        bad = _make_pc(z=1.0, M=200.0, Tb_K=500.0, SG=float('nan'))
        with pytest.raises(ValueError, match="nan"):
            compute_K_W_per_pseudocomponent([bad])

    def test_high_kw_gives_low_gamma(self):
        # High K_W (paraffinic, K_W close to 13) → gamma near 0
        pc = _make_pc(z=1.0, M=120.0, Tb_K=390.0, SG=0.72)   # expect K_W ~ 12.9
        result = compute_K_W_per_pseudocomponent([pc])
        assert result[0].gamma < 0.3, (
            f"Expected low gamma for paraffinic component, got {result[0].gamma:.4f}"
        )

    def test_low_kw_gives_high_gamma(self):
        # Low K_W (aromatic, K_W close to 9.5) → gamma near 1
        pc = _make_pc(z=1.0, M=300.0, Tb_K=700.0, SG=1.10)   # expect K_W ~ 9.9
        result = compute_K_W_per_pseudocomponent([pc])
        assert result[0].gamma > 0.7, (
            f"Expected high gamma for aromatic component, got {result[0].gamma:.4f}"
        )

    def test_asphaltene_convention_gets_kw_and_gamma(self):
        # Asphaltene: Tb_K=1073.15, SG=1.15 (numerical convention).
        # K_W and gamma are computed but not used for PC-SAFT; still populated.
        asp = _make_pc(z=0.01, M=1700.0, Tb_K=1073.15, SG=1.15)
        result = compute_K_W_per_pseudocomponent([asp])
        assert result[0].K_W is not None
        assert result[0].gamma is not None
        assert 0.0 <= result[0].gamma <= 1.0


# ── generate_pcsaft_table ─────────────────────────────────────────────────────

class TestGeneratePCSAFTTable:
    """DataFrame output and dispatch logic for generate_pcsaft_table."""

    @pytest.fixture
    def distillable_components(self):
        """Three distillable pseudo-components with K_W and gamma set."""
        raw = [
            _make_pc(z=0.6, M=150.0, Tb_K=500.0, SG=0.85),
            _make_pc(z=0.3, M=250.0, Tb_K=650.0, SG=0.92),
            _make_pc(z=0.1, M=400.0, Tb_K=800.0, SG=1.00),
        ]
        return compute_K_W_per_pseudocomponent(raw)

    @pytest.fixture
    def mixed_components(self, distillable_components):
        """Two distillable + one asphaltene (Tb_K=1073.15)."""
        asp = Pseudocomponent(
            z=0.02, M=1700.0, Tb_K=1073.15, SG=1.15,
            xc_lower=0.98, xc_upper=1.0,
            K_W=10.5, gamma=0.7,   # K_W/gamma set (compute would give them too)
        )
        return distillable_components[:2] + [asp]

    def test_returns_dataframe(self, distillable_components):
        df = generate_pcsaft_table(distillable_components)
        assert isinstance(df, pd.DataFrame)

    def test_row_count_matches_input(self, distillable_components):
        df = generate_pcsaft_table(distillable_components)
        assert len(df) == len(distillable_components)

    def test_expected_columns(self, distillable_components):
        df = generate_pcsaft_table(distillable_components)
        expected = {
            'component_type', 'z', 'M', 'Tb_K', 'SG',
            'K_W', 'gamma', 'm', 'sigma_A', 'eps_over_k_K',
        }
        assert set(df.columns) == expected

    def test_distillable_type_label(self, distillable_components):
        df = generate_pcsaft_table(distillable_components)
        assert (df['component_type'] == 'distillable').all(), (
            f"Expected all 'distillable': {df['component_type'].tolist()}"
        )

    def test_asphaltene_type_label(self, mixed_components):
        df = generate_pcsaft_table(mixed_components)
        asp_rows = df[df['Tb_K'] > 1000.0]
        assert len(asp_rows) == 1
        assert asp_rows.iloc[0]['component_type'] == 'asphaltene'

    def test_asphaltene_gets_gonzalez_defaults(self, mixed_components):
        df = generate_pcsaft_table(mixed_components)
        asp_row = df[df['Tb_K'] > 1000.0].iloc[0]
        m_asp, sigma_asp, epsk_asp = gonzalez_asphaltene_params()
        assert asp_row['m']            == pytest.approx(m_asp,    abs=0.0)
        assert asp_row['sigma_A']      == pytest.approx(sigma_asp, abs=0.0)
        assert asp_row['eps_over_k_K'] == pytest.approx(epsk_asp,  abs=0.0)

    def test_distillable_params_match_direct_call(self, distillable_components):
        df = generate_pcsaft_table(distillable_components)
        for i, c in enumerate(distillable_components):
            m_exp, s_exp, e_exp = panuganti_aromatic_resin_params(c.M, c.gamma)
            row = df.iloc[i]
            assert row['m']            == pytest.approx(m_exp, rel=1e-9)
            assert row['sigma_A']      == pytest.approx(s_exp, rel=1e-9)
            assert row['eps_over_k_K'] == pytest.approx(e_exp, rel=1e-9)

    def test_raises_if_gamma_none(self):
        # Distillable component with gamma=None should raise ValueError.
        bad = Pseudocomponent(z=0.5, M=200.0, Tb_K=500.0, SG=0.88,
                              K_W=11.5, gamma=None)
        with pytest.raises(ValueError, match="gamma=None"):
            generate_pcsaft_table([bad])

    def test_m_column_positive(self, distillable_components):
        df = generate_pcsaft_table(distillable_components)
        assert (df['m'] > 0).all()

    def test_sigma_column_positive(self, distillable_components):
        df = generate_pcsaft_table(distillable_components)
        assert (df['sigma_A'] > 0).all()

    def test_epsk_column_positive(self, distillable_components):
        df = generate_pcsaft_table(distillable_components)
        assert (df['eps_over_k_K'] > 0).all()

    def test_input_field_passthrough(self, distillable_components):
        # z, M, Tb_K, SG, K_W, gamma should pass through unchanged.
        df = generate_pcsaft_table(distillable_components)
        for i, c in enumerate(distillable_components):
            row = df.iloc[i]
            assert row['z']     == pytest.approx(c.z,     rel=1e-12)
            assert row['M']     == pytest.approx(c.M,     rel=1e-12)
            assert row['Tb_K']  == pytest.approx(c.Tb_K,  rel=1e-12)
            assert row['SG']    == pytest.approx(c.SG,    rel=1e-12)
            assert row['K_W']   == pytest.approx(c.K_W,   rel=1e-12)
            assert row['gamma'] == pytest.approx(c.gamma, rel=1e-12)

    def test_row_order_preserved(self, distillable_components):
        df = generate_pcsaft_table(distillable_components)
        for i, c in enumerate(distillable_components):
            assert df.iloc[i]['M'] == pytest.approx(c.M, rel=1e-12)
