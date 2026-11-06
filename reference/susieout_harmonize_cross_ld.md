# Harmonize WGS cross-LD rows to analysis alleles

Harmonize WGS cross-LD rows to analysis alleles

## Usage

``` r
susieout_harmonize_cross_ld(
  cross_ld,
  wgs_examined,
  analysis_examined,
  allow_strand = TRUE
)
```

## Arguments

- cross_ld:

  Tested-by-omitted cross-LD matrix.

- wgs_examined:

  Data frame with `id`, `effect_allele`, and `other_allele`.

- analysis_examined:

  Data frame with the same columns in analysis order.

- allow_strand:

  Whether unambiguous strand complements may be aligned.

## Value

The row-aligned, sign-corrected cross-LD matrix.
