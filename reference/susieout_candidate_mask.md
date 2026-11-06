# Apply realistic candidate-set masks

Apply realistic candidate-set masks

## Usage

``` r
susieout_candidate_mask(
  variants,
  maf_min = 0.01,
  snp_only = TRUE,
  info_min = 0.8
)
```

## Arguments

- variants:

  Data frame containing `maf`, `variant_type`, `info`,
  `reference_present`, and `panel_polymorphic`.

- maf_min:

  Minimum minor allele frequency for a variant to be tested.

- snp_only:

  If `TRUE`, only variants with `variant_type == "SNP"` are tested.

- info_min:

  Minimum imputation quality (`info`) for a variant to be tested.

## Value

A logical vector, `TRUE` for variants that pass every filter and are
therefore tested.

## Examples

``` r
variants <- data.frame(
  maf = c(0.2, 0.004, 0.3), variant_type = c("SNP", "SNP", "INDEL"),
  info = c(0.95, 0.9, 0.99), reference_present = TRUE,
  panel_polymorphic = TRUE
)
susieout_candidate_mask(variants)
#> [1]  TRUE FALSE FALSE
```
