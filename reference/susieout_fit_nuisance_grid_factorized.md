# Fit the joint SuSiE-OUT nuisance model from a low-rank genotype reference

This is an exact computational backend for LD matrices constructed from
a finite reference panel. It uses a single thin singular-value
decomposition of the reference genotypes while fitting the same null,
tested-variant, OUT, infinitesimal, and LD-mismatch states as
[`susieout_fit_nuisance_grid()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit_nuisance_grid.md).

## Usage

``` r
susieout_fit_nuisance_grid_factorized(
  z,
  n,
  examined_genotypes,
  outside_genotypes,
  examined_ids = NULL,
  excluded_ids = NULL,
  effect_prior = "robust_mixture",
  L = 10L,
  count_prior = susieout_count_prior(),
  nuisance_prior,
  trained_prior = NULL,
  ridge = 1e-06,
  max_iter = 200L,
  tol = 1e-05
)
```

## Arguments

- z, n, examined_ids, effect_prior, L, count_prior, trained_prior:

  Arguments analogous to
  [`susieout_fit_nuisance_grid()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit_nuisance_grid.md).

- examined_genotypes:

  Reference dosage matrix for tested variants, with samples in rows.

- outside_genotypes:

  Reference dosage matrix for omitted variants, with the same samples
  and ordering.

- excluded_ids:

  IDs for columns of `outside_genotypes`.

- nuisance_prior:

  Frozen grid prior from
  [`susieout_train_nuisance_grid()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_train_nuisance_grid.md).

- ridge, max_iter, tol:

  Numerical controls.

## Value

An `susieout_fit` averaged over the learned joint nuisance grid.
