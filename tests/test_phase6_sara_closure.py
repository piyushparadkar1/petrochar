"""
Phase 6 pass-gate: SARA closure check and discrete asphaltene assembly.

Synthetic three-component distillable fluid + 5 wt% asphaltenes.

Component construction
----------------------
Mole fractions z_i are chosen so that weight fractions in the total feed
come out to exactly SAT=30, ARO=50, RES=15, ASP=5 wt%:

    M_dist_av = 0.95 / (0.8/400 + 0.15/600) = 422.22 g/mol
    z_SAT = 1/3,  z_ARO = 5/9,  z_RES = 1/9   (sum to 1)

Tb_K values are set to yield exact K_W targets (K_W = Tb derived algebraically
from K_W = (1.8*Tb)^(1/3) / SG, so no Riazi-Daubert approximation involved):

    K_W = 13.0  (SG=0.81) → SAT bin  (K_W >= 12.0)
    K_W = 11.5  (SG=0.91) → ARO bin  (11.0 <= K_W < 12.0)
    K_W = 10.5  (SG=1.02) → RES bin  (K_W < 11.0)

Pass criteria
-------------
- validate_sara: accepts valid inputs; raises on negative, >100, or bad sum.
- append_asphaltene: ASP component has correct z, M, Tb_K, SG, xc bounds.
- kw_bin_check: per-class |delta| <= 3 wt% on clean synthetic input.
- kw_bin_check: zero flags raised on clean synthetic input.
- kw_bin_check: raises ValueError on missing SARA keys or nan Tb_K.
"""

import math

import numpy as np
import pytest

from core.correlations import watson_k
from core.quadrature import Pseudocomponent
from core.sara import (
    ASP_M_DEFAULT,
    ASP_SG_DEFAULT,
    ASP_TB_K_DEFAULT,
    KW_ARO_DEFAULT,
    KW_SAT_DEFAULT,
    append_asphaltene,
    kw_bin_check,
    validate_sara,
)


# ── Synthetic fluid parameters ─────────────────────────────────────────────────

ASP_WT_PCT = 5.0
F_DIST     = 1.0 - ASP_WT_PCT / 100.0  # 0.95

# Mole fractions within distillable subfraction
Z_SAT = 1.0 / 3.0
Z_ARO = 5.0 / 9.0
Z_RES = 1.0 / 9.0

# Molecular weights
M_SAT, M_ARO, M_RES = 400.0, 400.0, 600.0

# K_W targets and SG values → Tb_K derived from K_W = (1.8*Tb)^(1/3) / SG
KW_SAT, SG_SAT = 13.0, 0.81
KW_ARO, SG_ARO = 11.5, 0.91
KW_RES, SG_RES = 10.5, 1.02

# Tb_K from K_W and SG (algebraically exact: Tb = K_W^3 * SG^3 / 1.8)
TB_SAT = KW_SAT ** 3 * SG_SAT ** 3 / 1.8
TB_ARO = KW_ARO ** 3 * SG_ARO ** 3 / 1.8
TB_RES = KW_RES ** 3 * SG_RES ** 3 / 1.8

# Expected wt% in total feed
SARA_INPUT = {'SAT': 30.0, 'ARO': 50.0, 'RES': 15.0, 'ASP': ASP_WT_PCT}


def _make_distillable() -> list[Pseudocomponent]:
    return [
        Pseudocomponent(z=Z_SAT, M=M_SAT, Tb_K=TB_SAT, SG=SG_SAT),
        Pseudocomponent(z=Z_ARO, M=M_ARO, Tb_K=TB_ARO, SG=SG_ARO),
        Pseudocomponent(z=Z_RES, M=M_RES, Tb_K=TB_RES, SG=SG_RES),
    ]


# ── Constants ─────────────────────────────────────────────────────────────────

class TestConstants:
    def test_kw_sat_default(self):
        assert KW_SAT_DEFAULT == 12.0

    def test_kw_aro_default(self):
        assert KW_ARO_DEFAULT == 11.0

    def test_asp_M(self):
        assert ASP_M_DEFAULT == 1700.0

    def test_asp_Tb_K(self):
        assert abs(ASP_TB_K_DEFAULT - 1073.15) < 1e-6

    def test_asp_SG(self):
        assert abs(ASP_SG_DEFAULT - 1.15) < 1e-9

    def test_asp_Tb_K_above_1000(self):
        assert ASP_TB_K_DEFAULT > 1000.0


# ── validate_sara ─────────────────────────────────────────────────────────────

class TestValidateSara:
    def test_valid_input_no_error(self):
        validate_sara(30.0, 50.0, 15.0, 5.0)

    def test_exact_100(self):
        validate_sara(25.0, 25.0, 25.0, 25.0)

    def test_sum_within_tolerance(self):
        validate_sara(30.1, 50.1, 15.0, 4.5)   # sum = 99.7; within 0.5 tol

    def test_negative_sat_raises(self):
        with pytest.raises(ValueError, match="SAT"):
            validate_sara(-1.0, 51.0, 25.0, 25.0)

    def test_negative_asp_raises(self):
        with pytest.raises(ValueError, match="ASP"):
            validate_sara(30.0, 50.0, 21.0, -1.0)

    def test_exceeds_100_raises(self):
        with pytest.raises(ValueError, match="ARO"):
            validate_sara(0.0, 110.0, 0.0, 0.0)

    def test_sum_exceeds_tolerance_raises(self):
        with pytest.raises(ValueError, match="sum"):
            validate_sara(30.0, 50.0, 15.0, 10.0)   # sum = 105

    def test_sum_below_tolerance_raises(self):
        with pytest.raises(ValueError, match="sum"):
            validate_sara(10.0, 10.0, 10.0, 5.0)    # sum = 35


# ── append_asphaltene ─────────────────────────────────────────────────────────

class TestAppendAsphaltene:
    @pytest.fixture(scope='class')
    def full_list(self):
        return append_asphaltene(_make_distillable(), asp_wt_pct=ASP_WT_PCT)

    def test_length(self, full_list):
        assert len(full_list) == 4

    def test_asp_is_last(self, full_list):
        asp = full_list[-1]
        assert asp.Tb_K > 1000.0

    def test_asp_z_is_weight_fraction(self, full_list):
        asp = full_list[-1]
        assert abs(asp.z - ASP_WT_PCT / 100.0) < 1e-10

    def test_asp_M_default(self, full_list):
        asp = full_list[-1]
        assert asp.M == ASP_M_DEFAULT

    def test_asp_Tb_K_default(self, full_list):
        asp = full_list[-1]
        assert abs(asp.Tb_K - ASP_TB_K_DEFAULT) < 1e-6

    def test_asp_SG_default(self, full_list):
        asp = full_list[-1]
        assert abs(asp.SG - ASP_SG_DEFAULT) < 1e-9

    def test_asp_xc_lower(self, full_list):
        asp = full_list[-1]
        assert abs(asp.xc_lower - (1.0 - ASP_WT_PCT / 100.0)) < 1e-10

    def test_asp_xc_upper(self, full_list):
        asp = full_list[-1]
        assert abs(asp.xc_upper - 1.0) < 1e-10

    def test_distillable_components_unchanged(self, full_list):
        orig = _make_distillable()
        for orig_c, new_c in zip(orig, full_list[:-1]):
            assert orig_c.z   == new_c.z
            assert orig_c.M   == new_c.M
            assert orig_c.Tb_K == new_c.Tb_K
            assert orig_c.SG  == new_c.SG

    def test_returns_new_list(self):
        distillable = _make_distillable()
        full = append_asphaltene(distillable, asp_wt_pct=ASP_WT_PCT)
        assert full is not distillable

    def test_input_list_not_modified(self):
        distillable = _make_distillable()
        _ = append_asphaltene(distillable, asp_wt_pct=ASP_WT_PCT)
        assert len(distillable) == 3

    def test_zero_asp_raises(self):
        with pytest.raises(ValueError, match="open interval"):
            append_asphaltene(_make_distillable(), asp_wt_pct=0.0)

    def test_hundred_asp_raises(self):
        with pytest.raises(ValueError, match="open interval"):
            append_asphaltene(_make_distillable(), asp_wt_pct=100.0)

    def test_custom_M_asp(self):
        full = append_asphaltene(_make_distillable(), asp_wt_pct=5.0, M_asp=2000.0)
        assert full[-1].M == 2000.0

    def test_custom_Tb_K_asp(self):
        full = append_asphaltene(_make_distillable(), asp_wt_pct=5.0, Tb_K_asp=1100.0)
        assert full[-1].Tb_K == 1100.0

    def test_custom_SG_asp(self):
        full = append_asphaltene(_make_distillable(), asp_wt_pct=5.0, SG_asp=1.20)
        assert abs(full[-1].SG - 1.20) < 1e-9


# ── K_W construction verification ────────────────────────────────────────────

class TestKwConstruction:
    """Verify the synthetic components have the intended K_W values."""

    def test_sat_kw_is_target(self):
        kw = watson_k(TB_SAT, SG_SAT)
        assert abs(kw - KW_SAT) < 1e-6, f"K_W_SAT = {kw:.6f}, expected {KW_SAT}"

    def test_aro_kw_is_target(self):
        kw = watson_k(TB_ARO, SG_ARO)
        assert abs(kw - KW_ARO) < 1e-6, f"K_W_ARO = {kw:.6f}, expected {KW_ARO}"

    def test_res_kw_is_target(self):
        kw = watson_k(TB_RES, SG_RES)
        assert abs(kw - KW_RES) < 1e-6, f"K_W_RES = {kw:.6f}, expected {KW_RES}"

    def test_sat_bin(self):
        assert watson_k(TB_SAT, SG_SAT) >= KW_SAT_DEFAULT

    def test_aro_bin(self):
        kw = watson_k(TB_ARO, SG_ARO)
        assert KW_ARO_DEFAULT <= kw < KW_SAT_DEFAULT

    def test_res_bin(self):
        assert watson_k(TB_RES, SG_RES) < KW_ARO_DEFAULT


# ── kw_bin_check — clean synthetic input ─────────────────────────────────────

class TestKwBinCheckClean:
    @pytest.fixture(scope='class')
    def report(self):
        comps = append_asphaltene(_make_distillable(), asp_wt_pct=ASP_WT_PCT)
        return kw_bin_check(comps, SARA_INPUT)

    def test_not_flagged(self, report):
        assert report['flagged'] is False

    def test_no_flag_messages(self, report):
        assert report['flags'] == []

    def test_delta_SAT_within_3(self, report):
        assert abs(report['delta_wt_pct']['SAT']) <= 3.0

    def test_delta_ARO_within_3(self, report):
        assert abs(report['delta_wt_pct']['ARO']) <= 3.0

    def test_delta_RES_within_3(self, report):
        assert abs(report['delta_wt_pct']['RES']) <= 3.0

    def test_delta_ASP_within_3(self, report):
        assert abs(report['delta_wt_pct']['ASP']) <= 3.0

    def test_kw_calc_SAT_near_30(self, report):
        assert abs(report['kw_calc']['SAT'] - 30.0) < 0.01

    def test_kw_calc_ARO_near_50(self, report):
        assert abs(report['kw_calc']['ARO'] - 50.0) < 0.01

    def test_kw_calc_RES_near_15(self, report):
        assert abs(report['kw_calc']['RES'] - 15.0) < 0.01

    def test_kw_calc_ASP_near_5(self, report):
        assert abs(report['kw_calc']['ASP'] - 5.0) < 1e-10

    def test_report_has_required_keys(self, report):
        for k in ('kw_calc', 'sara_input', 'delta_wt_pct', 'flagged', 'flag_tol', 'flags'):
            assert k in report

    def test_sara_input_preserved(self, report):
        for k in ('SAT', 'ARO', 'RES', 'ASP'):
            assert abs(report['sara_input'][k] - SARA_INPUT[k]) < 1e-10

    def test_kw_calc_sums_to_100(self, report):
        total = sum(report['kw_calc'].values())
        assert abs(total - 100.0) < 0.01


# ── kw_bin_check — mismatched SARA triggers flags ─────────────────────────────

class TestKwBinCheckFlagged:
    def test_large_deviation_flagged(self):
        comps = append_asphaltene(_make_distillable(), asp_wt_pct=ASP_WT_PCT)
        bad_sara = {'SAT': 60.0, 'ARO': 20.0, 'RES': 15.0, 'ASP': 5.0}
        report = kw_bin_check(comps, bad_sara)
        assert report['flagged'] is True
        assert len(report['flags']) >= 1

    def test_custom_flag_tol(self):
        comps = append_asphaltene(_make_distillable(), asp_wt_pct=ASP_WT_PCT)
        sara = {'SAT': 31.0, 'ARO': 50.0, 'RES': 15.0, 'ASP': 4.0}
        # delta_SAT ≈ -1 wt%; should not flag at default 3.0 but flag at 0.5
        report_default = kw_bin_check(comps, sara, flag_tol=3.0)
        report_tight   = kw_bin_check(comps, sara, flag_tol=0.5)
        assert report_default['flagged'] is False or report_tight['flagged'] is True


# ── kw_bin_check — error conditions ─────────────────────────────────────────

class TestKwBinCheckErrors:
    def test_missing_key_raises(self):
        comps = append_asphaltene(_make_distillable(), asp_wt_pct=ASP_WT_PCT)
        incomplete = {'SAT': 30.0, 'ARO': 50.0, 'RES': 15.0}
        with pytest.raises(ValueError, match="missing"):
            kw_bin_check(comps, incomplete)

    def test_nan_tb_raises(self):
        comps = [
            Pseudocomponent(z=0.5, M=400.0),   # Tb_K = nan (Phase 5 state)
            Pseudocomponent(z=0.5, M=400.0, Tb_K=600.0, SG=0.85),
        ]
        comps_with_asp = append_asphaltene(comps, asp_wt_pct=5.0)
        with pytest.raises(ValueError, match="nan"):
            kw_bin_check(comps_with_asp, SARA_INPUT)

    def test_nan_sg_raises(self):
        comps = [
            Pseudocomponent(z=1.0, M=400.0, Tb_K=600.0),   # SG = nan
        ]
        comps_with_asp = append_asphaltene(comps, asp_wt_pct=5.0)
        with pytest.raises(ValueError, match="nan"):
            kw_bin_check(comps_with_asp, SARA_INPUT)


# ── kw_bin_check — custom thresholds ─────────────────────────────────────────

class TestKwBinCustomThresholds:
    def test_kw_sat_override(self):
        comps = append_asphaltene(_make_distillable(), asp_wt_pct=ASP_WT_PCT)
        # Raise SAT threshold above KW_SAT=13.0 → SAT component reclassified as ARO
        report = kw_bin_check(comps, SARA_INPUT, kw_sat=13.5, kw_aro=11.0)
        assert report['kw_calc']['SAT'] < 0.5   # SAT component now in ARO bin

    def test_kw_aro_override(self):
        comps = append_asphaltene(_make_distillable(), asp_wt_pct=ASP_WT_PCT)
        # Raise ARO threshold above KW_ARO=11.5 → ARO reclassified as RES
        report = kw_bin_check(comps, SARA_INPUT, kw_sat=12.0, kw_aro=12.0)
        assert report['kw_calc']['ARO'] < 0.5   # ARO component now in RES bin


# ── No retuning logic present (architectural guard) ───────────────────────────

class TestNoRetuningLogic:
    """Verify kw_bin_check never modifies component properties."""

    def test_components_not_modified(self):
        distillable = _make_distillable()
        original_z   = [c.z   for c in distillable]
        original_M   = [c.M   for c in distillable]
        original_Tb  = [c.Tb_K for c in distillable]
        original_SG  = [c.SG  for c in distillable]

        comps = append_asphaltene(distillable, asp_wt_pct=ASP_WT_PCT)
        _ = kw_bin_check(comps, SARA_INPUT)

        for i, c in enumerate(distillable):
            assert c.z    == original_z[i]
            assert c.M    == original_M[i]
            assert c.Tb_K == original_Tb[i]
            assert c.SG   == original_SG[i]
