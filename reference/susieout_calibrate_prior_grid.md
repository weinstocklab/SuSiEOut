# Calibrate sparse-grid weights with masked-WGS targets

Conditional fits use masked statistics only. Reference PIPs from
unfiltered fits enter only this across-locus hyperparameter calibration,
and callers must keep calibration and evaluation loci disjoint.

## Usage

``` r
susieout_calibrate_prior_grid(
  calibration_bank,
  prior_grid = susieout_prior_grid(),
  dirichlet_strength = NULL,
  iterations = 500L
)
```

## Arguments

- calibration_bank:

  List of loci containing `evidence`, a variants-by-states `pip_in`
  matrix, and `target_pip`.

- prior_grid:

  Starting grid and Dirichlet center.

- dirichlet_strength:

  Effective strength of the grid-weight penalty.

- iterations:

  Optimizer iteration limit.

## Value

A learned `susieout_prior_grid` with calibration diagnostics.
