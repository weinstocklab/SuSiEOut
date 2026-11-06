# Learn global sparse-grid weights from retained-data evidence

This sensitivity analysis updates only global grid weights. A Dirichlet
anchor prevents weakly identified retained-data likelihoods from
replacing the externally specified prior.

## Usage

``` r
susieout_train_prior_grid(
  simulation_bank,
  prior_grid = susieout_prior_grid(),
  nuisance_prior = NULL,
  trained_prior = NULL,
  L = 10L,
  iterations = 50L,
  anchor_strength = length(simulation_bank),
  tol = 1e-06,
  ...
)
```

## Arguments

- simulation_bank:

  Full training loci.

- prior_grid:

  Starting grid and anchor weights.

- nuisance_prior:

  Optional frozen nuisance grid.

- trained_prior:

  Optional frozen effect/nuisance prior.

- L:

  Component limit.

- iterations:

  Maximum mixture-weight updates.

- anchor_strength:

  Dirichlet effective sample size. By default it equals the number of
  training loci.

- tol:

  Weight convergence tolerance.

- ...:

  Additional
  [`susieout_fit()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit.md)
  arguments.

## Value

A learned, frozen `susieout_prior_grid`.
