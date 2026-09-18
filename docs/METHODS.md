# Model and interpretation

## Explicit line-shape convention

With natural units for energy-like quantities, the implemented function is

```text
sigma(E) = A * E² * Gamma² / [(E² - m²)² + m² Gamma²].
```

E, m and Gamma are in GeV; A is in nb. The fraction is dimensionless and sigma(m) = A. Because of the E² numerator, m is not exactly the position of the maximum. The width in the denominator is constant. This is the original coursework shape with a directly identifiable normalization, not a full electroweak scattering calculation or an s-dependent-width LEP parameterization.

The original function used `12*pi*E**2*gamma_ee**2 / denominator * 389.4`. Its fraction was dimensionless, yet its multiplier was described as a GeV⁻²-to-nb conversion. In addition, 0.3894 mb corresponds to approximately 389,400 nb, not 389.4 nb. Consequently the original fitted `gamma_ee` cannot be reported as a calibrated physical partial width.

Reparameterizing the same original shape uses `A = 12*pi*389.4*gamma_ee_old**2/Gamma**2`; it preserves the attainable curves without giving the legacy amplitude a physical partial-width interpretation. The new fit does not require a cross-section unit conversion because A is fitted directly in the CSV's stated unit. A physical partial width would require a verified channel, a fully specified normalization and appropriate physics corrections.

For background, see the [PDG resonance formation review](https://pdg.lbl.gov/2024/reviews/rpp2024-rev-cross-section-formulae.pdf) and the [LEP/SLC precision Z-resonance report](https://cds.cern.ch/record/892831). They provide physical context, not a provenance claim for these CSVs.

## Selection and fitting

1. Read the header followed by every three-column CSV row; preserve source file and line number.
2. Reject nonnumeric/nonfinite values and nonpositive energy or uncertainty.
3. Retain only `85 < E < 95` and `0 < sigma < 3`, as in the original submission. These are coursework analysis choices. In particular, a negative background-subtracted measurement is not universally invalid, and 3 nb is not a fundamental physical bound.
4. Fit mass, width and pole cross section by weighted nonlinear least squares with `absolute_sigma=True`. Bounds are 85–95 GeV for mass, 0.01–20 GeV for width, and nonnegative A. The included fit is away from the bounds.
5. In the default mode, discard points with an absolute standardized residual at least 5, based on that one pre-fit. Re-fit once. This is not iterative clipping. Failure of the pre-fit stops execution instead of silently changing the procedure.
6. Save both fits. Two-sigma reference lines in the plot do not cause additional exclusions.

Observations are not automatically deduplicated: repeated measurements can be legitimate. Rows outside validity/windows are listed in the CSV audit but omitted from the fit figure. Residual-cut exclusions remain visible in that figure.

## Uncertainty and lifetime

The parameter covariance is the local least-squares approximation using the supplied absolute measurement errors. Chi-square uses N−3 degrees of freedom. Correlated detector errors, energy uncertainties, backgrounds and radiative corrections are not modeled.

`tau = hbar_GeV_seconds / Gamma`, with `hbar_GeV_seconds = scipy.constants.hbar / (1e9 * scipy.constants.electron_volt)`. A width in GeV is already an energy; no additional c² belongs in this conversion. First-order propagation gives `sigma_tau = tau * sigma_Gamma / Gamma`.

For each mass and width in the contour grid, the model is linear in A. With shape f and weights w = 1/error², its optimum is `A = max(0, sum(w*f*y)/sum(w*f*f))`. This profiles out A instead of holding it fixed. The reference minimum is the fitted chi-square, not the smallest sampled grid value.

The two displayed levels are nominal Gaussian/asymptotic two-parameter contours. Covariance errors and contour coverage condition on the retained sample; neither accounts for the residual-based selection. The pre-cut comparison reveals some sensitivity but does not quantify total systematic uncertainty.
