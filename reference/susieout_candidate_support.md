# Variant support under SuSiE-OUT

Lists every tested variant and every omitted variant (OUT state) with
its posterior inclusion probability (PIP), sorted from highest to
lowest.

## Usage

``` r
susieout_candidate_support(fit)
```

## Arguments

- fit:

  A `susieout_fit` object from
  [`susieout_fit()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit.md)
  or one of the grid-averaged fitting functions.

## Value

A data frame with columns `variant_id`, `state` (`"inside"` for tested
variants, `"outside"` for omitted variants), `pip`, and
`maximum_component_probability`, the largest probability that any single
effect assigns to the variant.

## Examples

``` r
locus <- susieout_simulate_locus(1, architecture = "entirely_omitted")
fit <- susieout_fit(
  locus$z, locus$n, locus$R_analysis, locus$wgs_ld, locus$examined_ids
)
susieout_candidate_support(fit)[1:5, ]
#>     variant_id   state        pip maximum_component_probability
#> v40        v40  inside 0.04196349                   0.004277781
#> v28        v28  inside 0.03861921                   0.003930730
#> v29        v29  inside 0.02933064                   0.002972515
#> v41        v41 outside 0.02656404                   0.002688715
#> v44        v44 outside 0.02528099                   0.002557339
```
