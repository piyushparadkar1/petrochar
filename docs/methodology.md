# petrochar — Methodology

*Draft for Paper 1, Section 2 (Methodology)*  
*Date: 2026-05-16*  
*Code version: Phase 10 (complete)*

---

## 2.1 Architecture overview

petrochar characterizes heavy petroleum fractions into discrete pseudo-components
with PC-SAFT molecular parameters from four routine refinery assay measurements:
a distillation curve (any ASTM method), bulk specific gravity at 15 °C, bulk
molecular weight, and SARA (saturates, aromatics, resins, asphaltenes) weight
fractions.

The characterization pipeline executes nine sequential steps:

1. Distillation curve ingestion and TBP conversion
2. Generalized Riazi distribution fit to the TBP curve (diagnostic)
3. Per-cut SG assignment via constant Watson K
4. Per-cut M via Riazi-Daubert inversion; M distribution fit
5. Gauss-Laguerre quadrature discretization of the M distribution
6. Discrete asphaltene component assembly
7. K_W-bin SARA closure check
8. Per-component Watson K and aromaticity factor
9. PC-SAFT parameter assignment

The outputs are seven components: five distillable pseudo-components from
quadrature, one discrete asphaltene, and one propane solvent component (the
last added for deasphalting or injection simulations). All scientific computation
is in `core/`; the Streamlit UI in `tabs/` and `app.py` is a wrapper with no
embedded correlations.

---

## 2.2 Input data and distillation curve handling

### 2.2.1 Distillation curve

A distillation curve supplies boiling temperature T_b as a function of cumulative
fraction x_c. Four ASTM method tags are accepted: D86, D1160 AET (atmospheric
equivalent temperature), D7169 (simulated distillation, conversion pending), and
TBP. The curve also carries a basis tag (weight, volume, or mole); the current
pipeline operates on weight-basis data. The DistillationCurve class (Phase 2)
enforces monotonicity and range checks at construction time.

**D86 → TBP conversion:** Daubert's method (Riazi MNL50 §3.2.2.2.1, Eqs. 3.20–3.22,
Table 3.7) is applied. Eq. 3.20 anchors TBP(50%) in Kelvin. Eq. 3.21 converts
D86 temperature gaps at six standard intervals (Table 3.7 constants; verified
against Riazi Example 3.3 kerosene). Eq. 3.22 assembles TBP values at all cut
points by summation around the anchor.

**D1160 AET → TBP:** D1160 at atmospheric equivalent temperature is a direct
pass-through; no numerical conversion is required (temperatures are already on
a TBP-equivalent basis for characterization purposes).

### 2.2.2 Interior points for distribution fitting

IBP (0%) data are excluded before fitting. The fitting domain is the open interval
(0, 1) in cumulative fraction, which excludes both endpoints where the Riazi
generalized distribution diverges. Eleven interior points (x_c = 0.05 to 0.95)
are used for the VTB synthetic feed; any positive integer of interior points ≥ 2
is accepted.

---

## 2.3 Generalized distribution fitting (Riazi Eq. 4.56)

### 2.3.1 The model

The Riazi (2005) generalized distribution (MNL50 §4.5.4.1, Eq. 4.56) relates any
petroleum property P to cumulative fraction x_c:

$$P(x_c) = P_o \left(1 + \left[\frac{A}{B}\ln\frac{1}{1-x_c}\right]^{1/B}\right)$$

where P_o is the onset parameter (P at x_c → 0), A is a scale parameter, and B is
a shape parameter. The corresponding PDF, inverse CDF (x_c as a function of P), and
analytical mean are all available in closed form:

$$\langle P \rangle = P_o \left(1 + \left(\frac{A}{B}\right)^{1/B} \Gamma\!\left(1+\frac{1}{B}\right)\right)$$

### 2.3.2 Fitting procedure

Fitting linearizes Eq. 4.56 via the substitution Y = ln(−ln(1−x_c)),
X = ln((P − P_o)/P_o), giving Y = B·X + ln(B/A). Non-linear least squares
(`scipy.optimize.least_squares`, trust-region-reflective method) minimizes
Σ r_i² where r_i = Y_i − B·X_i − c for three-parameter mode (P_o, A, B free) or
two-parameter mode (B fixed at a literature value; Riazi p. 174 recommends
B_T = 1.5 for T_b, B_M = 1 for M, B_SG = 3 for SG).

**Validation:** Three-parameter fits reproduce Riazi Table 4.13 reference values
for the North Sea C7+ data (Example 4.13) with P_o, A, B each within 5% for M,
T_b, and SG distributions. The fitted M distribution for the VTB synthetic feed
gives an analytical mean of 369.95 g/mol.

### 2.3.3 Role of the T_b distribution

The T_b distribution characterizes the feed's distillation curve and is used
diagnostically (RMS, AAD, R² reported in the UI). It is **not** used to assign
pseudo-component boiling points. Pseudo-component T_b values are derived from
the Riazi-Daubert correlation via the self-consistent (M, K_W) solve described
in §2.4.1. This distinction is Decision 25 in the architectural record and is
the subject of §2.9.1.

---

## 2.4 SG and MW per pseudo-component

### 2.4.1 Constant Watson K method (default)

The Watson K factor (Watson and Nelson 1933; Riazi MNL50 Eq. 2.13) couples T_b
and SG for petroleum fractions of similar chemical character:

$$K_W = \frac{(1.8 \, T_b)^{1/3}}{SG}$$

Under the constant Watson K assumption (Whitson 1983), all cuts share the bulk
K_W derived from bulk T_b and bulk SG. Per-cut SG then follows algebraically:
SG_i = (1.8·T_b,i)^(1/3) / K_W.

**Bulk T_b estimation:** The bulk T_b corresponding to the bulk M and SG is
obtained from a single Riazi-Daubert call (Eq. 2.56 for M ≤ 300, Eq. 2.57 for
M > 300). The resulting K_W is applied uniformly to all distillable cuts.

**Consequence for SARA classification:** Under constant Watson K, all distillable
pseudo-components share the same K_W and therefore fall in the same K_W bin.
SAT/ARO/RES K_W-bin classification is degenerate by construction. The K_W-bin
check (§2.8) detects and flags this. Users requiring finer SARA-class resolution
should supply per-cut SG measurements and use the generalized SG distribution
mode (§2.4.2).

### 2.4.2 Generalized SG distribution (alternative)

`core/sg_distribution.py` implements a generalized distribution fit to
(x_cv, SG) data for cases where measured per-cut SG values are available.
This path requires volume-basis cumulative fraction data and measured SG at
each cut point. It is architecturally supported but not exposed in the current
UI pipeline, because per-cut SG measurements are not routinely available for
heavy residue streams (see §2.9.2).

### 2.4.3 Molecular weight per cut: Riazi-Daubert inversion

Per-cut M is obtained by numerically inverting the Riazi-Daubert T_b(M, SG)
correlation. Two equations apply depending on M:

- **Eq. 2.56** (valid M = 70–300 g/mol):
  T_b = 3.76587·exp(3.7741×10⁻³·M + 2.98404·SG − 4.25288×10⁻³·M·SG)·M^0.40167·SG^−1.58262
- **Eq. 2.57** (valid M = 300–700 g/mol, recommended for M > 300):
  T_b = 9.3369·exp(1.6514×10⁻⁴·M + 1.4103·SG − 7.5152×10⁻⁴·M·SG)·M^0.5369·SG^−0.7276

The inversion uses `scipy.optimize.brentq`. Eq. 2.57 is non-monotone in M (peak
at M_peak = 0.5369/(7.5152×10⁻⁴·SG − 1.6514×10⁻⁴)); the search is confined to
the ascending branch [300.01, M_peak·0.9999]. A regime-gap fallback returns M = 300
with a UserWarning when the target T_b falls in the discontinuity between the
two equations (~636–668 K at SG ≈ 0.94).

---

## 2.5 Gauss-Laguerre quadrature discretization

### 2.5.1 Quadrature scheme

The fitted M distribution is discretized by 5-point (default) or 3-point
Gauss-Laguerre quadrature following Riazi (2005) §4.6.1.1. The quadrature
nodes y_i and weights w_i are taken directly from Riazi Table 4.21 (vendored
as `data/riazi_reference/table_4_21_quadrature_points.csv`). The transformation
from the standard Laguerre domain to cumulative fraction and property value is:

$$x_{c,i} = 1 - e^{-y_i}, \qquad M_i = P(x_{c,i}), \qquad z_i = w_i$$

where P(·) is the fitted generalized distribution CDF.

**Validation:** The 3-point quadrature reproduces Riazi Table 4.22 (Example 4.14)
with M_i within 1% and z_i within 0.005 for all three nodes. The 5-point result
gives M_av within 0.5% of the distribution's analytical mean. The in-app
Validation tab reproduces both results live.

### 2.5.2 Self-consistent (T_b, SG) assignment per node

Each quadrature node carries M_i and z_i from the M distribution. The corresponding
T_b,i and SG_i are obtained by solving the coupled system:

$$SG_i = \frac{(1.8 \, T_{b,i})^{1/3}}{K_W}, \qquad T_{b,i} = T_b^{RD}(M_i,\, SG_i)$$

where T_b^{RD} is the Riazi-Daubert correlation (regime-selecting). Substituting
the first equation into the second gives a single implicit equation in T_b,i,
solved by `scipy.optimize.brentq` on [300, 990] K. SG_i then follows algebraically.

This approach — deriving (T_b, SG) from (M, K_W) rather than evaluating the
T_b distribution at the quadrature node — ensures that pseudo-component boiling
points are physically consistent with their molecular weights and lie within the
Riazi-Daubert validity range. The T_b distribution is a diagnostic of the feed's
distillation curve, not a path to pseudo-component property assignment.

For the VTB synthetic feed: all five quadrature nodes produce T_b ∈ [629, 870] K,
well below the 1000 K threshold that would indicate extrapolation beyond the
physical petroleum regime.

---

## 2.6 Asphaltene treatment

The asphaltene pseudo-component is discrete and literature-defined. It is appended
after quadrature discretization via `sara.append_asphaltene()` and carries the
Gonzalez et al. (2007) nanoaggregate default properties:

| Property | Value | Source |
|----------|-------|--------|
| M_asp    | 1700 g/mol   | Gonzalez et al. (2007) |
| T_b,asp  | 1073.15 K (800 °C) | Numerical convention (asphaltenes do not boil) |
| SG_asp   | 1.15         | Gonzalez et al. (2007) |
| m        | 33           | Gonzalez et al. (2007) |
| σ        | 4.3 Å        | Gonzalez et al. (2007) |
| ε/k      | 400 K        | Gonzalez et al. (2007) |

The asphaltene mass fraction is set by the measured SARA ASP wt% input. The
mole fraction z_asp follows from a mole-basis conversion using M_asp and the
distillable sub-fraction MW (see `core/sara.py`). All z values in the final
pseudo-component list are true mole fractions in the full mixture (distillable
+ asphaltene), summing to 1.

**Asphaltene identity flag:** The `Pseudocomponent.is_asphaltene` boolean field
identifies the asphaltene component throughout the pipeline. It is set exclusively
by `sara.append_asphaltene()` via `is_asphaltene=True`. No code in `core/` uses a
T_b > 1000 K threshold to classify asphaltenes; a ValueError is raised if any
non-asphaltene component is found with T_b > 1000 K (which would indicate
distribution extrapolation beyond the physical regime).

---

## 2.7 PC-SAFT parameter assignment

### 2.7.1 Distillable pseudo-components

PC-SAFT parameters (m, σ, ε/k) for distillable pseudo-components use the
gamma-interpolated aromatic/resin correlations from Panuganti et al. (2012)
Table 6 (rows 4–9). The correlations interpolate linearly between two
end-members as a function of the aromaticity factor γ:

$$m     = (1-\gamma)(0.0223\,M + 0.751)       + \gamma(0.0101\,M + 1.7296)$$
$$\sigma = (1-\gamma)(4.1377 - 38.1483/M)      + \gamma(4.6169 - 93.98/M) \quad\text{[Å]}$$
$$\varepsilon/k = (1-\gamma)(0.00436\,M + 283.93) + \gamma(508 - 234100/M^{1.5}) \quad\text{[K]}$$

At γ = 0 (paraffinic character, K_W → 13) the benzene-derivative end-member is
recovered. At γ = 1 (poly-nuclear-aromatic character, K_W → 9.5) the PNA
end-member is recovered.

The separate Saturates correlation from Panuganti Table 6 rows 1–3 is retained
in `core/pcsaft_params.py` for reference and independent testing but is not
used on the main characterization path. The gamma-interpolated form already
encodes the saturate–PNA continuum; using the Saturates form in parallel would
double-count the paraffinic character already captured by γ.

### 2.7.2 Aromaticity factor γ from Watson K

The aromaticity factor γ is derived deterministically per pseudo-component from
the Watson K factor via a linear clamp:

$$\gamma_i = \text{clamp}\!\left(\frac{13.0 - K_{W,i}}{13.0 - 9.5},\; 0,\; 1\right)$$

This mapping assigns γ = 0 to the most paraffinic fractions (K_W ≥ 13) and
γ = 1 to the most aromatic (K_W ≤ 9.5), with linear interpolation in between.
The boundaries 9.5 and 13.0 correspond to the K_W thresholds conventionally
associated with poly-nuclear-aromatic and paraffinic petroleum fractions
respectively (Riazi MNL50 p. 75).

### 2.7.3 Asphaltene component

The discrete asphaltene component receives the Gonzalez et al. (2007) nanoaggregate
defaults (m = 33, σ = 4.3 Å, ε/k = 400 K) unconditionally. These are the initial
PC-SAFT parameters reported before AOP fitting in Panuganti (2012). petrochar uses
these defaults without adjustment because AOP data are outside the tool's scope.

### 2.7.4 Propane (solvent component)

Pure-component propane PC-SAFT parameters from Gross and Sadowski (2001) as
reproduced in Panuganti (2012) Table 5: m = 2.002, σ = 3.6180 Å, ε/k = 208.11 K.

**Caution:** The value σ = 3.168 Å appears in some Aspen Plus built-in component
databases as a transcription error. The correct value is 3.618 Å. This is enforced
by a defensive test (`test_sigma_is_not_aspen_typo` in `tests/test_phase7_pcsaft.py`).

---

## 2.8 SARA closure as informational check

The Watson K bin check (`core/sara.py`, `kw_bin_check`) aggregates distillable
pseudo-components into SAT (K_W ≥ 12), ARO (11 ≤ K_W < 12), and RES (K_W < 11)
bins (Riazi MNL50 p. 75 conventions) and compares the K_W-binned wt% against the
measured SARA input. The comparison is purely informational: deviations flag
potential inconsistency in the input data and should prompt inspection. Distribution
parameters are never adjusted to close the SARA balance.

The asphaltene mass fraction in the K_W-bin check exactly matches the SARA ASP input
because the asphaltene component is assembled from that fraction directly. SAT/ARO/RES
closure depends on the feed's actual boiling-point–Watson-K relationship and will only
reproduce SARA values accurately if the constant Watson K assumption is consistent
with the feed's compositional gradient.

---

## 2.9 Methodological deviations from prior work

### 2.9.1 T_b distribution as diagnostic, not assignment mechanism

Panuganti (2012) and most prior characterization frameworks assign pseudo-component
T_b by evaluating the fitted T_b distribution at the quadrature node positions.
petrochar does not do this.

The reason is architectural: evaluating the T_b distribution at Gauss-Laguerre
nodes (which are positioned in the M distribution's CDF space) conflates two
independent distribution fits. The M distribution places nodes at specific M values;
the T_b distribution places nodes at specific cumulative fraction values. These are
only equivalent if T_b and M co-vary monotonically and if the two fitted distributions
are mutually consistent — a coincidence not generally guaranteed for VTB feeds where
the T_b distribution is fitted on the distillation curve and the M distribution is
derived from it via the Riazi-Daubert inverse.

petrochar instead derives T_b,i from (M_i, SG_i) directly via Riazi-Daubert, which
guarantees that each pseudo-component's T_b is physically consistent with its MW
and lies within the correlation's validity range. The T_b distribution fit is
retained as a diagnostic: it characterizes the shape of the distillation curve and
its fit quality is reported, but it has no role in the numerical values of the
pseudo-component properties.

Consequence: all five quadrature nodes for the VTB synthetic feed have T_b ∈
[629, 870] K. No node exceeds 1000 K. The pathological "high-MW quadrature node
above 1000 K is silently reclassified as asphaltene" failure mode is structurally
impossible in petrochar's pipeline.

### 2.9.2 Constant Watson K as the only operative SG method

`core/sg_distribution.py` implements both the constant Watson K method and the
generalized SG distribution (Riazi Eq. 4.56 fitted to measured SG data). The current
UI pipeline exposes only the constant Watson K method, because per-cut SG measurements
are not routinely available for heavy residue streams in refinery assay reports.

Under the constant Watson K assumption, all distillable pseudo-components share
K_W,bulk and therefore all fall in the same K_W bin. SAT/ARO/RES K_W-bin
classification is degenerate by construction for this method. Users requiring
finer SARA classification should provide per-cut SG data and select the generalized
SG distribution mode (future version).

### 2.9.3 Aromaticity factor γ from Watson K (not from AOP fitting)

Panuganti (2012) defines γ as a free parameter fitted simultaneously to crude
density and asphaltene onset pressure (AOP) for each crude. This fitting allows
γ to capture crude-level aromatic character information that Watson K cannot
distinguish — two crudes with identical distillation curves but different aromatic
contents will have different optimal γ values but identical Watson K sequences.

petrochar uses a deterministic per-component mapping γ = clamp((13.0 − K_W) /
(13.0 − 9.5), 0, 1) derived from Watson K. This eliminates the need for AOP data
entirely — the characterization runs from routine refinery assay data alone — but
at the cost of losing the crude-level density-matching capability that AOP fitting
provides. The PC-SAFT parameters will not reproduce measured crude densities as
accurately as Panuganti's AOP-fitted values, and the model cannot be calibrated to
a specific crude's phase behaviour without external AOP data.

This trade-off is the defining methodological choice of petrochar and must be
stated plainly in Paper 1 §2 so that readers understand what the tool does and
does not do.

### 2.9.4 Bulk MW as an input, not a pass-gate target

The bulk MW input (M_bulk) is a feed characterization measurement. It is used to
derive the bulk Watson K factor but is not used as a closure target for the
pseudo-component distribution. The distribution is fitted to the distillation curve;
the bulk MW serves as a consistency diagnostic (reported as the ratio of GL M_av
to M_DIST_TARGET, where M_DIST_TARGET = (M_bulk − f_asp·M_asp)/(1−f_asp)).

For the synthetic VTB test feed, the test data are internally inconsistent (the
D1160 endpoint at x_c = 0.95 gives M ≈ 527 g/mol, but bulk MW = 700 g/mol implies
M ≈ 4200 g/mol in the unmeasured 5% tail — unphysical). The distillable-fraction
M distribution analytic mean (369.95 g/mol) is correctly derived from the measured
distillation data; the departure from M_DIST_TARGET = 563.6 g/mol is a known
inconsistency in the synthetic test inputs, not a model error.

---

## References

- Gonzalez, D. L., Hirasaki, G. J., Creek, J., Chapman, W. G. (2007). Modeling study
  of CO₂-induced asphaltene precipitation utilizing the PC-SAFT equation of state.
  *Energy & Fuels*, 21, 1230–1234.
- Gross, J. and Sadowski, G. (2001). Perturbed-chain SAFT: an equation of state based
  on a perturbation theory for chain molecules. *Industrial & Engineering Chemistry
  Research*, 40, 1244–1260.
- Panuganti, S. R., Tavakkoli, M., Vargas, F. M., Gonzalez, D. L., Chapman, W. G.
  (2012). PC-SAFT characterization of crude oils and modeling of asphaltene phase
  behavior. *Fuel*, 93, 658–669.
- Riazi, M. R. (2005). *Characterization and Properties of Petroleum Fractions*.
  ASTM Manual MNL50. ISBN 0-8031-3361-8.
- Watson, K. M. and Nelson, E. F. (1933). Improved methods for approximating critical
  and thermal properties of petroleum fractions. *Industrial & Engineering Chemistry*,
  25, 880.
- Whitson, C. H. (1983). Characterizing hydrocarbon plus fractions. *SPE Journal*,
  23(4), 683–694.
