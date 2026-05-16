"""
Per-pseudo-component Watson K factor and aromaticity gamma assignment.

For each pseudo-component in the list:
    K_W_i = watson_k(Tb_K_i, SG_i)          ...Phase 1 Eq. 2.13
    gamma_i = aromaticity_gamma(K_W_i)       ...linear clamp [9.5, 13.0]

The aromaticity factor gamma is petrochar's deterministic proxy for the
Panuganti 2012 gamma parameter.  Panuganti defines gamma as a free parameter
fitted to crude density and bubble pressure simultaneously (Panuganti 2012,
page 4).  petrochar substitutes a deterministic per-component mapping from
Watson K (Watson and Nelson 1933) so that gamma requires no AOP data and no
tuning against process observables.

This is an intentional architectural deviation from Panuganti's original
procedure.  See Decision 21 in CURRENT_STATUS.md for the explicit
methodological note required for Paper 1 disclosure.

References
----------
Watson and Nelson (1933) Ind. Eng. Chem. 25, 880 (Watson K).
Panuganti et al. (2012) Fuel 93, 658-669 (gamma parameter definition, p. 4).
Architecture commitment: CLAUDE_CODE_PROMPT.md §ARCHITECTURE COMMITMENTS.
"""

from __future__ import annotations

import math

from core.correlations import aromaticity_gamma, watson_k
from core.quadrature import Pseudocomponent


def compute_K_W_per_pseudocomponent(
    components: list[Pseudocomponent],
) -> list[Pseudocomponent]:
    """Assign Watson K and aromaticity gamma to every pseudo-component.

    Creates a new list of Pseudocomponent objects with K_W and gamma
    populated.  Input list is not modified in place.

    The asphaltene component (Tb_K > 1000 K convention) receives K_W and
    gamma computed from its conventional Tb_K = 1073.15 K and SG = 1.15.
    These values are not used in PC-SAFT parameter assignment (the Gonzalez
    2007 defaults are used instead), but they are recorded for completeness
    and consistency.

    Parameters
    ----------
    components : list[Pseudocomponent]
        Pseudo-component list, typically from sara.append_asphaltene().
        Every component must have Tb_K and SG set (non-nan).  Components
        with nan Tb_K or SG raise ValueError — assign Tb_K and SG before
        calling (e.g., from the Tb distribution quadrature and SG from the
        constant-Watson-K or SG distribution methods of Phase 4).

    Returns
    -------
    list[Pseudocomponent]
        New list of the same length; K_W and gamma populated for every
        component.  All other fields (z, M, Tb_K, SG, xc_lower, xc_upper)
        are preserved unchanged.

    Raises
    ------
    ValueError
        If any component has nan Tb_K or nan SG.

    References
    ----------
    Watson and Nelson (1933) Ind. Eng. Chem. 25, 880 (Eq. 2.13).
    Architecture: gamma = clamp((13.0 - K_W) / (13.0 - 9.5), 0, 1).
    """
    result = []
    for i, c in enumerate(components):
        if math.isnan(c.Tb_K) or math.isnan(c.SG):
            raise ValueError(
                f"compute_K_W_per_pseudocomponent: component {i} "
                f"(M={c.M:.1f} g/mol) has nan Tb_K or SG.  "
                f"Assign Tb_K and SG before calling this function."
            )
        kw    = watson_k(c.Tb_K, c.SG)
        gamma = aromaticity_gamma(kw)
        result.append(
            Pseudocomponent(
                z=c.z,
                M=c.M,
                Tb_K=c.Tb_K,
                SG=c.SG,
                xc_lower=c.xc_lower,
                xc_upper=c.xc_upper,
                K_W=kw,
                gamma=gamma,
                is_asphaltene=c.is_asphaltene,   # preserve asphaltene identity
            )
        )
    return result
