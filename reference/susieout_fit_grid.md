# Fit SuSiE-OUT by averaging sparse-count and nuisance states

Every sparse-count state is combined with every nuisance state.
Posterior state weights use the complete variational evidence stored by
the conditional fit, including likelihood normalization and variational
uncertainty.

## Usage

``` r
susieout_fit_grid(
  z,
  n,
  R_analysis,
  wgs_ld = NULL,
  examined_ids = NULL,
  effect_prior = "robust_mixture",
  L = 10L,
  prior_grid = susieout_prior_grid(),
  nuisance_prior = NULL,
  trained_prior = NULL,
  infinitesimal_variance = 0,
  ld_mismatch_variance = 0,
  starts = c("zero", "inside", "mixed"),
  keep_fits = FALSE,
  ...
)
```

## Arguments

- z, n, R_analysis, wgs_ld, examined_ids, effect_prior, L:

  Arguments passed to
  [`susieout_fit()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit.md).

- prior_grid:

  A grid returned by
  [`susieout_prior_grid()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_prior_grid.md).

- nuisance_prior:

  Optional frozen nuisance prior returned by
  [`susieout_train_nuisance_grid()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_train_nuisance_grid.md).

- trained_prior:

  Optional frozen effect/nuisance prior. Its count prior is replaced by
  each sparse-grid state.

- infinitesimal_variance, ld_mismatch_variance:

  Singleton nuisance values used when `nuisance_prior` is `NULL`.

- starts:

  Deterministic initialization strategies. The highest-evidence
  converged fit is retained for each grid state.

- keep_fits:

  Retain conditional fits in diagnostics.

- ...:

  Additional arguments passed to
  [`susieout_fit()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit.md).

## Value

A model-averaged `susieout_fit`.
