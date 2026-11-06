# Construct an effect-size mixture prior

The robust mixture is the primary prior and gives equal architecture
mass to neutral, stabilizing, and recent directional-optimum-shift
components.

## Usage

``` r
susieout_effect_prior(
  architecture = c("robust_mixture", "generic", "neutral", "stabilizing", "directional")
)
```

## Arguments

- architecture:

  One of `robust_mixture`, `generic`, `neutral`, `stabilizing`, or
  `directional`.

## Value

A normalized list of normal-mixture weights, means, and variances.

## Examples

``` r
str(susieout_effect_prior("robust_mixture"))
#> List of 5
#>  $ weights     : num [1:9] 0.2167 0.1 0.0167 0.1167 0.15 ...
#>  $ means       : num [1:9] 0 0 0 0 0 0 0 0.012 -0.012
#>  $ variances   : num [1:9] 1.0e-05 1.0e-04 1.0e-03 2.5e-05 4.0e-04 4.0e-03 5.0e-05 5.0e-04 5.0e-04
#>  $ labels      : chr [1:9] "neutral_small" "neutral_medium" "neutral_large" "stabilizing_small" ...
#>  $ architecture: chr "robust_mixture"
#>  - attr(*, "class")= chr "susieout_effect_prior"
```
