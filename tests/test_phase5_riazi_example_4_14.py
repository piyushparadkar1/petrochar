"""
Phase 5 pass-gate: Gaussian quadrature discretization (Riazi Example 4.14).

Reference: Riazi MNL50 Table 4.22 (3-point quadrature on M distribution).
Input parameters (Table 4.22 header): M_o=90, A_M=0.3324, B_M=1.096.
These differ from Table 4.13 fitted values and are used directly via from_params.

Pass criteria
-------------
M_i  : within +-1%   of Table 4.22 values (103.6, 154.6, 252.2 g/mol).
z_i  : within +-0.005 of Table 4.22 values (0.711, 0.279, 0.010).
M_av : within +-0.5%  of experimental 118.9 g/mol (Riazi 3-pt result: 119.4).
5-pt : M_av also within +-0.5% of 118.9 g/mol.
"""

import math

import numpy as np
import pytest

from core.distribution import GeneralizedDistribution
from core.quadrature import Pseudocomponent, discretize_generalized, quadrature_points

# ── Table 4.22 parameters and reference data ──────────────────────────────────

M_O   = 90.0
A_M   = 0.3324
B_M   = 1.096

# (index 1-based, y_i, z_i_ref, M_i_ref)
TABLE_4_22 = [
    (1, 0.416, 0.711, 103.6),
    (2, 2.294, 0.279, 154.6),
    (3, 6.290, 0.010, 252.2),
]

M_AV_EXP = 118.9   # experimental bulk (Table 4.11)


# ── quadrature_points ─────────────────────────────────────────────────────────

class TestQuadraturePoints:
    def test_3pt_roots(self):
        y, _ = quadrature_points(3)
        assert y.shape == (3,)
        np.testing.assert_allclose(y, [0.41577, 2.29428, 6.28995], rtol=1e-5)

    def test_3pt_weights(self):
        _, w = quadrature_points(3)
        assert w.shape == (3,)
        np.testing.assert_allclose(w, [0.711093, 0.278518, 0.0103893], rtol=1e-5)

    def test_5pt_roots(self):
        y, _ = quadrature_points(5)
        assert y.shape == (5,)
        np.testing.assert_allclose(
            y,
            [0.26356, 1.41340, 3.59643, 7.08581, 12.64080],
            rtol=1e-5,
        )

    def test_5pt_weights(self):
        _, w = quadrature_points(5)
        assert w.shape == (5,)
        np.testing.assert_allclose(
            w,
            [0.521756, 0.398667, 0.0759424, 0.00361176, 0.0000233700],
            rtol=1e-5,
        )

    def test_3pt_weights_sum_to_one(self):
        _, w = quadrature_points(3)
        assert abs(w.sum() - 1.0) < 1e-6

    def test_5pt_weights_sum_to_one(self):
        _, w = quadrature_points(5)
        assert abs(w.sum() - 1.0) < 1e-6

    def test_invalid_N_raises(self):
        with pytest.raises(ValueError, match="3 or 5"):
            quadrature_points(4)

    def test_returns_copies(self):
        y1, w1 = quadrature_points(3)
        y2, w2 = quadrature_points(3)
        y1[0] = 999.0
        y3, _ = quadrature_points(3)
        assert y3[0] != 999.0   # original table unchanged


# ── Pseudocomponent dataclass ─────────────────────────────────────────────────

class TestPseudomponentDataclass:
    def test_required_fields(self):
        pc = Pseudocomponent(z=0.5, M=120.0)
        assert pc.z == 0.5
        assert pc.M == 120.0
        assert math.isnan(pc.Tb_K)
        assert math.isnan(pc.SG)
        assert pc.xc_lower == 0.0
        assert pc.xc_upper == 0.0

    def test_all_fields_settable(self):
        pc = Pseudocomponent(z=0.1, M=150.0, Tb_K=520.0, SG=0.85,
                             xc_lower=0.2, xc_upper=0.5)
        assert pc.Tb_K == 520.0
        assert pc.SG == 0.85
        assert pc.xc_lower == 0.2
        assert pc.xc_upper == 0.5

    def test_eight_fields_present(self):
        # Phase 5 original fields plus K_W and gamma added in Phase 7.
        pc = Pseudocomponent(z=0.5, M=100.0)
        expected = {'z', 'M', 'Tb_K', 'SG', 'xc_lower', 'xc_upper', 'K_W', 'gamma'}
        assert set(pc.__dataclass_fields__) == expected

    def test_kw_gamma_default_none(self):
        pc = Pseudocomponent(z=0.5, M=100.0)
        assert pc.K_W is None
        assert pc.gamma is None


# ── 3-point discretization (Riazi Table 4.22) ────────────────────────────────

class TestDiscretize3pt:
    @pytest.fixture(scope='class')
    def components(self):
        dist = GeneralizedDistribution.from_params(M_O, A_M, B_M)
        return discretize_generalized(3, dist)

    def test_returns_three(self, components):
        assert len(components) == 3

    def test_all_pseudocomponent_instances(self, components):
        assert all(isinstance(c, Pseudocomponent) for c in components)

    @pytest.mark.parametrize("idx,_y,z_ref,M_ref", TABLE_4_22)
    def test_z_i_within_tolerance(self, components, idx, _y, z_ref, M_ref):
        z_calc = components[idx - 1].z
        delta = abs(z_calc - z_ref)
        assert delta <= 0.005, (
            f"i={idx}: z={z_calc:.5f}, ref={z_ref:.3f}, delta={delta:.5f} > 0.005"
        )

    @pytest.mark.parametrize("idx,_y,z_ref,M_ref", TABLE_4_22)
    def test_M_i_within_1pct(self, components, idx, _y, z_ref, M_ref):
        M_calc = components[idx - 1].M
        pct = abs(M_calc - M_ref) / M_ref * 100.0
        assert pct <= 1.0, (
            f"i={idx}: M={M_calc:.2f}, ref={M_ref:.1f}, dev={pct:.3f}% > 1%"
        )

    def test_mixture_M_av_within_half_pct(self, components):
        z = np.array([c.z for c in components])
        M = np.array([c.M for c in components])
        M_av = float(np.sum(z * M))   # = sum z_i*M_i (mole-fraction weighted)
        pct = abs(M_av - M_AV_EXP) / M_AV_EXP * 100.0
        assert pct <= 0.5, (
            f"3-pt M_av={M_av:.2f}, exp={M_AV_EXP}, dev={pct:.3f}% > 0.5%"
        )

    def test_xc_bounds_ascending(self, components):
        for i, c in enumerate(components):
            assert c.xc_lower < c.xc_upper, (
                f"i={i}: xc_lower={c.xc_lower:.6f} >= xc_upper={c.xc_upper:.6f}"
            )

    def test_xc_bounds_contiguous(self, components):
        for i in range(len(components) - 1):
            gap = abs(components[i].xc_upper - components[i + 1].xc_lower)
            assert gap < 1e-10, (
                f"gap between component {i} upper and {i+1} lower = {gap:.2e}"
            )

    def test_xc_spans_zero_to_one(self, components):
        assert abs(components[0].xc_lower) < 1e-10
        assert abs(components[-1].xc_upper - 1.0) < 1e-6

    def test_Tb_K_is_nan(self, components):
        assert all(math.isnan(c.Tb_K) for c in components)

    def test_SG_is_nan(self, components):
        assert all(math.isnan(c.SG) for c in components)

    def test_z_sum_to_one(self, components):
        total = sum(c.z for c in components)
        assert abs(total - 1.0) < 1e-6, f"z sum = {total:.8f}"

    def test_M_i_ascending(self, components):
        M_vals = [c.M for c in components]
        assert all(M_vals[i] < M_vals[i + 1] for i in range(len(M_vals) - 1)), (
            f"M values not ascending: {M_vals}"
        )


# ── 5-point discretization ────────────────────────────────────────────────────

class TestDiscretize5pt:
    @pytest.fixture(scope='class')
    def components(self):
        dist = GeneralizedDistribution.from_params(M_O, A_M, B_M)
        return discretize_generalized(5, dist)

    def test_returns_five(self, components):
        assert len(components) == 5

    def test_M_av_within_half_pct(self, components):
        z = np.array([c.z for c in components])
        M = np.array([c.M for c in components])
        M_av = float(np.sum(z * M))
        pct = abs(M_av - M_AV_EXP) / M_AV_EXP * 100.0
        assert pct <= 0.5, (
            f"5-pt M_av={M_av:.2f}, exp={M_AV_EXP}, dev={pct:.3f}% > 0.5%"
        )

    def test_xc_spans_zero_to_one(self, components):
        assert abs(components[0].xc_lower) < 1e-10
        assert abs(components[-1].xc_upper - 1.0) < 1e-6

    def test_z_sum_to_one(self, components):
        total = sum(c.z for c in components)
        assert abs(total - 1.0) < 1e-6

    def test_M_i_ascending(self, components):
        M_vals = [c.M for c in components]
        assert all(M_vals[i] < M_vals[i + 1] for i in range(len(M_vals) - 1))

    def test_Tb_K_is_nan(self, components):
        assert all(math.isnan(c.Tb_K) for c in components)

    def test_SG_is_nan(self, components):
        assert all(math.isnan(c.SG) for c in components)


# ── from_params used in discretization (sanity) ───────────────────────────────

class TestFromParamsIntegration:
    def test_from_params_sets_attributes(self):
        dist = GeneralizedDistribution.from_params(M_O, A_M, B_M)
        assert dist.P_o == M_O
        assert dist.A   == A_M
        assert dist.B   == B_M

    def test_from_params_forward_eval(self):
        dist = GeneralizedDistribution.from_params(M_O, A_M, B_M)
        M_at_half = dist.P(0.5)
        assert M_at_half > M_O      # monotonically increasing CDF

    def test_from_params_no_fit_quality(self):
        dist = GeneralizedDistribution.from_params(M_O, A_M, B_M)
        with pytest.raises((RuntimeError, TypeError, AttributeError)):
            _ = dist.fit_quality
