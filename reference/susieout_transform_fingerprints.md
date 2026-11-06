# Ridge-stabilized whitening and recoloring of WGS fingerprints

Computes `R_analysis^(1/2) R_wgs_SS^(-1/2) cross_ld` with eigenvalue
regularization, including rank-deficient input matrices. Use it when LD
among tested variants in the sequencing panel differs from the analysis
LD, so that the fingerprints are expressed on the analysis LD scale.

## Usage

``` r
susieout_transform_fingerprints(R_analysis, R_wgs_SS, cross_ld, ridge = 1e-06)
```

## Arguments

- R_analysis:

  LD matrix among tested variants used in the analysis.

- R_wgs_SS:

  LD matrix among the same tested variants in the sequencing reference
  panel.

- cross_ld:

  Reference-panel LD between tested variants (rows) and omitted variants
  (columns).

- ridge:

  Positive eigenvalue floor used for regularization.

## Value

A matrix of transported fingerprints with the same dimensions as
`cross_ld`.
