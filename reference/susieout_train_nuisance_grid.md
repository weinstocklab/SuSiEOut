# Learn a spike-and-grid nuisance prior from training loci

Conditional dense/mismatch likelihood scores are computed once per
training locus. A finite-mixture EM then estimates joint grid weights
without using architecture labels or nuisance truth.

## Usage

``` r
susieout_train_nuisance_grid(
  simulation_bank,
  trained_prior,
  infinitesimal_grid = c(0, 0.002, 0.008, 0.032),
  ld_mismatch_grid = c(0, 0.04, 0.12),
  L = 10L,
  iterations = 50L,
  max_loci = 500L,
  dirichlet_strength = 1,
  tol = 1e-06,
  seed = 1L,
  ...
)
```

## Arguments

- simulation_bank:

  List of full loci labelled `training`.

- trained_prior:

  Frozen sparse/outside prior used in conditional fits.

- infinitesimal_grid, ld_mismatch_grid:

  Nuisance grids.

- L:

  Number of sparse single-effect components.

- iterations:

  Maximum mixture-weight EM iterations.

- max_loci:

  Maximum number of training loci.

- dirichlet_strength:

  Total pseudocount assigned according to the default spike prior.

- tol:

  Maximum grid-weight change for convergence.

- seed:

  Training-locus subsampling seed.

- ...:

  Additional
  [`susieout_fit()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit.md)
  arguments.

## Value

A frozen `susieout_nuisance_prior` accepted by
[`susieout_fit_nuisance_grid()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit_nuisance_grid.md).
