# Train and freeze L-invariant SuSiE-OUT hyperparameters

Sparse count rates are estimated once from known training counts.
Optional variational EM updates only effect-mixture weights and nuisance
covariance; it never feeds inferred OUT responsibilities back into the
OUT count prior.

## Usage

``` r
susieout_train_prior(
  simulation_bank,
  effect_prior = "robust_mixture",
  L = 10L,
  em_iterations = 0L,
  em_max_loci = 500L,
  em_update = c("effect_and_nuisance", "nuisance_only"),
  initial_infinitesimal_variance = 0,
  initial_ld_mismatch_variance = 0,
  em_tol = 1e-04,
  em_pip_tol = 0.001,
  seed = 1L,
  verbose = FALSE
)
```

## Arguments

- simulation_bank:

  Training/validation data frame or full-locus list.

- effect_prior:

  Starting effect-size mixture.

- L:

  Component limit used by the fitted training models.

- em_iterations:

  Number of effect/nuisance EM iterations.

- em_max_loci:

  Maximum number of training loci used by EM.

- em_update:

  Update effect weights and nuisance parameters or nuisance parameters
  only.

- initial_infinitesimal_variance:

  Initial total local dense variance.

- initial_ld_mismatch_variance:

  Initial independent mismatch variance.

- em_tol:

  Relative parameter convergence tolerance.

- em_pip_tol:

  Maximum PIP-change convergence tolerance.

- seed:

  Training-locus subsampling seed.

- verbose:

  Report EM progress.

## Value

A frozen `susieout_prior` accepted by
[`susieout_fit()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit.md).
