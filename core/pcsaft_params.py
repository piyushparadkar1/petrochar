"""
PC-SAFT parameter assignment for petroleum pseudo-components.

Three parameter sources:

1. Distillable pseudo-components (is_asphaltene=False):
   Panuganti 2012 Table 6 gamma-interpolated correlations.
   gamma-interpolation spans gamma=0 (benzene-derivative/paraffinic)
   to gamma=1 (PNA/poly-nuclear-aromatic).

2. Asphaltene discrete component (is_asphaltene=True):
   Gonzalez 2007 nanoaggregate initial defaults — m=33, sigma=4.3 A,
   eps/k=400 K.  These are the INITIAL values before AOP fitting in
   Panuganti's procedure.  petrochar does not fit to AOP data.
   Asphaltenes are identified by the is_asphaltene flag (set exclusively
   by sara.append_asphaltene), NOT by a Tb > 1000 K threshold.

3. Propane (pure component, from Table 5):
   m=2.002, sigma=3.6180 A, eps/k=208.11 K.
   WARNING: sigma=3.618 NOT 3.168 — 3.168 is a known typo in some Aspen
   built-in databases.  The correct value is from Gross and Sadowski 2001,
   reproduced as Panuganti Table 5.

Design note — single gamma-interpolated form for all distillable components
------------------------------------------------------------------
Panuganti 2012 uses two separate correlation sets: one for Saturates and one
for Aromatics+Resins.  petrochar uses only the gamma-interpolated A+R form
(Table 6 rows 4-9) for ALL distillable components, with gamma derived from
Watson K.  The gamma-interpolation already encodes the saturate-to-PNA
continuum: at gamma=0, the A+R form reduces to the benzene-derivative
(paraffinic) end-member; at gamma=1 it reaches the PNA end-member.  Using the
separate Saturates form (rows 1-3) in parallel would double-count the
paraffinic character already captured by gamma.  The Saturate-specific
correlation is retained as `panuganti_saturate_params` for reference and
independent testing.

References
----------
Panuganti et al. (2012) Fuel 93, 658-669.  Table 6 (page 663): saturate and
   aromatic+resin correlations.  Table 5 (page 662): propane params.
Gonzalez et al. (2007) Energy & Fuels 21, 1230-1234 (asphaltene defaults,
   cited as ref [32] in Panuganti 2012).
Gross J. and Sadowski G. (2001) Ind. Eng. Chem. Res. 40, 1244 (propane
   pure-component values, reproduced in Panuganti Table 5).
"""

from __future__ import annotations

import math

import pandas as pd

from core.quadrature import Pseudocomponent


# ── Gonzalez 2007 asphaltene defaults ────────────────────────────────────────

_ASP_M_DEFAULT     = 33.0    # dimensionless
_ASP_SIGMA_DEFAULT =  4.3    # Angstrom
_ASP_EPSK_DEFAULT  = 400.0   # K

# ── Propane (Table 5, Gross & Sadowski 2001 via Panuganti 2012) ───────────────

_C3_M      =   2.002   # dimensionless
_C3_SIGMA  =   3.6180  # Angstrom; NOT 3.168 (Aspen typo)
_C3_EPSK   = 208.11    # K

# ── Guard thresholds for non-asphaltene Tb values ────────────────────────────
# Distillable (is_asphaltene=False, is_heavy_resin=False): hard ceiling 1000 K.
# Tb > 1000 K for a distillable signals distribution extrapolation beyond the
# measured range — characterisation error, raise ValueError.
# Heavy-resin lump (is_heavy_resin=True): ceiling relaxed to 1150 K because
# Tb_hr is derived from constant Watson K closure and is not a physical
# boiling point (Decision 31, Phase 11).

_DISTILLABLE_TB_MAX = 1000.0   # K; hard ceiling for ordinary distillable nodes
_HR_TB_MAX_PCSAFT   = 1150.0   # K; relaxed ceiling for heavy-resin lump only


# ── Saturates (Panuganti Table 6, rows 1-3) ───────────────────────────────────

def panuganti_saturate_params(MW: float) -> tuple[float, float, float]:
    """PC-SAFT parameters for saturate pseudo-components (Panuganti Table 6).

    Single correlation set, no aromaticity parameter.

        m      = 0.0257 * MW + 0.8444
        sigma  = 4.047 - 4.8013 * ln(MW) / MW           [Angstrom]
        eps/k  = exp(5.5769 - 9.523 / MW)               [K]

    Note on sigma: the natural-log form (ln(MW)/MW) is as published in
    Table 6; it is distinct from the linear forms used in the A+R correlation.

    Note on eps/k: the left-hand side in Table 6 is ln(eps/k), so the
    correlation gives eps/k = exp(5.5769 - 9.523 / MW).

    Parameters
    ----------
    MW : float
        Molecular weight, g/mol.  Valid range per Panuganti: ~100-700 g/mol.

    Returns
    -------
    m : float   Number of PC-SAFT segments (dimensionless).
    sigma : float  Segment diameter, Angstrom.
    eps_over_k : float  Dispersion energy parameter, K.

    References
    ----------
    Panuganti et al. (2012) Fuel 93, 658-669, Table 6, page 663.
    """
    m         = 0.0257 * MW + 0.8444
    sigma     = 4.047 - 4.8013 * math.log(MW) / MW
    eps_over_k = math.exp(5.5769 - 9.523 / MW)
    return float(m), float(sigma), float(eps_over_k)


# ── Aromatics + Resins (Panuganti Table 6, rows 4-9) ──────────────────────────

def panuganti_aromatic_resin_params(
    MW: float,
    gamma: float,
) -> tuple[float, float, float]:
    """PC-SAFT parameters for aromatics+resins via gamma-interpolation.

    Linear interpolation between benzene-derivative (gamma=0) and
    PNA/poly-nuclear-aromatic (gamma=1) end-members (Panuganti Table 6):

        m     = (1-g)*(0.0223*MW + 0.751)        + g*(0.0101*MW + 1.7296)
        sigma = (1-g)*(4.1377 - 38.1483/MW)      + g*(4.6169 - 93.98/MW)   [A]
        eps/k = (1-g)*(0.00436*MW + 283.93)       + g*(508 - 234100/MW^1.5) [K]

    At gamma=0 (paraffinic, K_W~13), returns the benzene-derivative form.
    At gamma=1 (aromatic/PNA, K_W~9.5), returns the PNA form.

    Parameters
    ----------
    MW : float
        Molecular weight, g/mol.
    gamma : float
        Aromaticity factor in [0, 1]; derived from Watson K via
        aromaticity_gamma(K_W) (core/correlations.py).

    Returns
    -------
    m, sigma, eps_over_k : float, float, float
        PC-SAFT segment number (dimensionless), segment diameter (Angstrom),
        dispersion energy (K).

    References
    ----------
    Panuganti et al. (2012) Fuel 93, 658-669, Table 6, page 663.
    """
    g = float(gamma)
    m         = (1.0 - g) * (0.0223 * MW + 0.751)       + g * (0.0101 * MW + 1.7296)
    sigma     = (1.0 - g) * (4.1377 - 38.1483 / MW)     + g * (4.6169 - 93.98 / MW)
    eps_over_k = (1.0 - g) * (0.00436 * MW + 283.93)    + g * (508.0 - 234100.0 / MW ** 1.5)
    return float(m), float(sigma), float(eps_over_k)


# ── Top-level distillable dispatcher ──────────────────────────────────────────

def panuganti_distillable_params(
    MW: float,
    gamma: float,
) -> tuple[float, float, float]:
    """PC-SAFT parameters for a distillable pseudo-component.

    Routes to panuganti_aromatic_resin_params(MW, gamma) for all distillable
    components.  The gamma-interpolated A+R correlation spans the full
    chemistry continuum from paraffinic (gamma=0) to PNA (gamma=1), so no
    separate Saturates branch is needed.  The separate Saturates correlation
    (panuganti_saturate_params) is retained for reference testing.

    Parameters
    ----------
    MW : float   Molecular weight, g/mol.
    gamma : float   Aromaticity factor in [0, 1].

    Returns
    -------
    m, sigma, eps_over_k : float, float, float

    References
    ----------
    Panuganti et al. (2012) Fuel 93, 658-669, Table 6.
    """
    return panuganti_aromatic_resin_params(MW, gamma)


# ── Fixed-parameter components ────────────────────────────────────────────────

def gonzalez_asphaltene_params() -> tuple[float, float, float]:
    """PC-SAFT parameters for the discrete asphaltene component.

    Returns the Gonzalez 2007 nanoaggregate initial defaults, which are also
    the starting values in Panuganti 2012's AOP-fitting procedure.  petrochar
    uses these defaults without adjustment because AOP data are outside scope.

        m     = 33      (dimensionless)
        sigma = 4.3     (Angstrom)
        eps/k = 400     (K)

    Returns
    -------
    m : float   = 33.0
    sigma : float   = 4.3 Angstrom
    eps_over_k : float   = 400.0 K

    References
    ----------
    Gonzalez et al. (2007) Energy & Fuels 21, 1230-1234 (nanoaggregate
    initial values, cited as ref [32] in Panuganti 2012, page 4).
    """
    return float(_ASP_M_DEFAULT), float(_ASP_SIGMA_DEFAULT), float(_ASP_EPSK_DEFAULT)


def propane_params() -> tuple[float, float, float]:
    """PC-SAFT parameters for pure propane (C3).

        m     = 2.002      (dimensionless)
        sigma = 3.6180     (Angstrom)   *** NOT 3.168 — Aspen typo ***
        eps/k = 208.11     (K)

    The sigma value 3.168 appears in some Aspen Plus built-in databases as a
    transcription error.  The correct value from Gross and Sadowski 2001 is
    3.618 (also reported in Panuganti 2012 Table 5 as 3.618).

    Returns
    -------
    m : float   = 2.002
    sigma : float   = 3.6180 Angstrom
    eps_over_k : float   = 208.11 K

    References
    ----------
    Gross J. and Sadowski G. (2001) Ind. Eng. Chem. Res. 40, 1244.
    Reproduced in Panuganti et al. (2012) Fuel 93, Table 5, page 662.
    """
    return float(_C3_M), float(_C3_SIGMA), float(_C3_EPSK)


# ── DataFrame generator ───────────────────────────────────────────────────────

def generate_pcsaft_table(
    components: list[Pseudocomponent],
) -> pd.DataFrame:
    """Assign PC-SAFT parameters to all pseudo-components; return DataFrame.

    For each component:
    - If is_asphaltene=True: gonzalez_asphaltene_params().
    - Otherwise (distillable): panuganti_distillable_params(M, gamma).

    Asphaltene identity is determined by the is_asphaltene flag (set by
    sara.append_asphaltene), NOT by a Tb threshold.  If any distillable
    component (is_asphaltene=False) has Tb_K > 1000 K, a ValueError is
    raised — this indicates a characterisation error (distribution
    extrapolation producing unphysical pseudo-component properties).

    Requires K_W and gamma to be populated on each distillable component
    (i.e., compute_K_W_per_pseudocomponent must have been called first).
    Raises ValueError if a distillable component has gamma = None.

    Parameters
    ----------
    components : list[Pseudocomponent]
        Augmented pseudo-component list with K_W and gamma set (Phase 7).

    Returns
    -------
    pd.DataFrame
        Columns: component_type, z, M, Tb_K, SG, K_W, gamma, m, sigma_A,
        eps_over_k_K.  One row per pseudo-component, ordered as supplied.

    References
    ----------
    Panuganti et al. (2012) Fuel 93, Table 6 (distillable correlations).
    Gonzalez et al. (2007) Energy & Fuels 21 (asphaltene defaults).
    """
    rows = []
    for i, c in enumerate(components):
        if c.is_asphaltene:
            m, sigma, epsk = gonzalez_asphaltene_params()
            comp_type = 'asphaltene'

        elif c.is_heavy_resin:
            # Heavy-resin lump uses the same Panuganti A+R gamma-interpolated
            # correlations as distillable nodes.  Its Tb_hr is a Watson-K-closure
            # convention, not a physical boiling point, so the Tb ceiling is
            # relaxed to _HR_TB_MAX_PCSAFT = 1150 K (Decision 31, Phase 11).
            if not math.isnan(c.Tb_K) and c.Tb_K > _HR_TB_MAX_PCSAFT:
                raise ValueError(
                    f"generate_pcsaft_table: heavy-resin lump (component {i}, "
                    f"M={c.M:.1f} g/mol) has Tb_K={c.Tb_K:.1f} K > "
                    f"{_HR_TB_MAX_PCSAFT:.0f} K.  Closure may be unreliable.  "
                    f"Check bulk MW/SG inputs."
                )
            if c.gamma is None:
                raise ValueError(
                    f"generate_pcsaft_table: heavy-resin component {i} "
                    f"(M={c.M:.1f}) has gamma=None.  Call "
                    f"compute_K_W_per_pseudocomponent before generate_pcsaft_table."
                )
            m, sigma, epsk = panuganti_distillable_params(c.M, c.gamma)
            comp_type = 'heavy_resin'

        else:
            # Distillable quadrature node — hard 1000 K ceiling.
            # Tb > 1000 K signals distribution extrapolation beyond the measured
            # range (Decision 26, Phase 8 rework).
            if not math.isnan(c.Tb_K) and c.Tb_K > _DISTILLABLE_TB_MAX:
                raise ValueError(
                    f"generate_pcsaft_table: distillable component {i} "
                    f"(M={c.M:.1f} g/mol) has Tb_K={c.Tb_K:.1f} K > "
                    f"{_DISTILLABLE_TB_MAX:.0f} K.  This indicates the M or "
                    f"Tb distribution extrapolated beyond the measured range. "
                    f"Use 3-parameter distribution fits; asphaltenes must be "
                    f"added exclusively via sara.append_asphaltene()."
                )
            if c.gamma is None:
                raise ValueError(
                    f"generate_pcsaft_table: component {i} (M={c.M:.1f}) "
                    f"has gamma=None.  Call compute_K_W_per_pseudocomponent "
                    f"before generate_pcsaft_table."
                )
            m, sigma, epsk = panuganti_distillable_params(c.M, c.gamma)
            comp_type = 'distillable'

        rows.append({
            'component_type': comp_type,
            'z':              c.z,
            'M':              c.M,
            'Tb_K':           c.Tb_K,
            'SG':             c.SG,
            'K_W':            c.K_W,
            'gamma':          c.gamma,
            'm':              m,
            'sigma_A':        sigma,
            'eps_over_k_K':   epsk,
        })

    return pd.DataFrame(rows)
