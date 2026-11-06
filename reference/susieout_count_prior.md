# Construct an L-invariant sparse-effect count prior

`lambda_total` is the locus-level sparse-effect rate. `out_fraction` is
the expected fraction of sparse effects assigned to variants outside the
tested candidate set. The prior is converted to component probabilities
only when a model is fitted, so its locus-level interpretation does not
change with the computational component limit `L`.

## Usage

``` r
susieout_count_prior(lambda_total = 1, out_fraction = 0.1)
```

## Arguments

- lambda_total:

  Nonnegative sparse-effect rate for the locus.

- out_fraction:

  Fraction of sparse effects assigned to OUT candidates.

## Value

A `susieout_count_prior` object.

## Examples

``` r
susieout_count_prior(lambda_total = 1, out_fraction = 0.10)
#> $lambda_total
#> [1] 1
#> 
#> $lambda_in
#> [1] 0.9
#> 
#> $lambda_out
#> [1] 0.1
#> 
#> $out_fraction
#> [1] 0.1
#> 
#> attr(,"class")
#> [1] "susieout_count_prior"
```
