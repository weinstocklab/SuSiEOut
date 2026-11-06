# Construct a spike-and-grid sparse-effect prior

The default grid contains a null-locus spike and positive count-rate
states spanning candidate-complete through entirely omitted
architectures.

## Usage

``` r
susieout_prior_grid(
  lambda_total = c(0, 0.5, 1, 2, 4, 6),
  out_fraction = c(0, 0.1, 0.3, 0.6, 1),
  weights = NULL,
  dirichlet_strength = 5
)
```

## Arguments

- lambda_total:

  Unique nonnegative total-effect rates.

- out_fraction:

  Unique OUT fractions in `[0, 1]`.

- weights:

  Optional weights for the resulting grid rows.

- dirichlet_strength:

  Regularization strength used when learning weights.

## Value

A frozen `susieout_prior_grid` object.
