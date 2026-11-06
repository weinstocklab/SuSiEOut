# Bootstrap WGS cross-LD fingerprints

Bootstrap WGS cross-LD fingerprints

## Usage

``` r
susieout_bootstrap_fingerprints(
  genotypes,
  examined,
  excluded,
  B = 100L,
  seed = 1L
)
```

## Arguments

- genotypes:

  Individual-by-variant dosage matrix.

- examined:

  Integer or character indices of tested variants.

- excluded:

  Integer or character indices of omitted variants.

- B:

  Number of bootstrap panels.

- seed:

  Reproducibility seed.

## Value

A list containing bootstrapped cross-LD matrices and uncertainty.
