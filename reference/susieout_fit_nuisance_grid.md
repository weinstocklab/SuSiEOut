# Fit a SuSiE-OUT model averaged over local nuisance states

This post-lock extension gives each locus an explicit spike at zero
dense variance and positive infinitesimal/mismatch alternatives.
Conditional fits are averaged using a covariance-comparable posterior
expected-likelihood score, so a global nuisance estimate need not be
imposed on candidate-complete sparse loci.

## Usage

``` r
susieout_fit_nuisance_grid(
  z,
  n,
  R_analysis,
  wgs_ld = NULL,
  examined_ids = NULL,
  effect_prior = "robust_mixture",
  L = 10L,
  count_prior = susieout_count_prior(),
  infinitesimal_grid = c(0, 0.002, 0.008, 0.032),
  ld_mismatch_grid = c(0, 0.04, 0.12),
  infinitesimal_prior = NULL,
  ld_mismatch_prior = NULL,
  nuisance_prior = NULL,
  trained_prior = NULL,
  ...
)
```

## Arguments

- z, n, R_analysis, wgs_ld, examined_ids, effect_prior, L, count_prior:

  Arguments passed to
  [`susieout_fit()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit.md).

- infinitesimal_grid:

  Grid of total local infinitesimal variances.

- ld_mismatch_grid:

  Grid of independent mismatch variances.

- infinitesimal_prior:

  Optional prior weights over the dense grid.

- ld_mismatch_prior:

  Optional prior weights over the mismatch grid.

- nuisance_prior:

  Optional learned joint grid prior returned by
  [`susieout_train_nuisance_grid()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_train_nuisance_grid.md).

- trained_prior:

  Optional frozen prior; its nuisance values are replaced by each grid
  point without modifying the original object.

- ...:

  Further arguments passed to
  [`susieout_fit()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit.md).

## Value

An `susieout_fit` whose PIPs and locus quantities are averaged over the
nuisance grid. `diagnostics$nuisance_grid` contains conditional scores
and posterior weights.
