# Simulate a deterministic SuSiE-OUT locus

This compact truth-known generator supports smoke tests and workflow
banks. It represents recent growth through a rare-allele-enriched
frequency draw, prospective causal assignment, derived-allele
orientation, finite WGS panels, ancestry mismatch, and
sparse/null/infinitesimal/LD-error regimes. Variants have first-order
autoregressive LD, and a realistic mix of MAF, indel,
imputation-quality, reference-absence and panel-monomorphism filters
decides which variants are tested and which are omitted.

## Usage

``` r
susieout_simulate_locus(
  seed,
  p = 40L,
  q = 10L,
  n = 1e+05,
  wgs_n = 1000L,
  architecture = c("contained", "partially_omitted", "entirely_omitted", "null",
    "infinitesimal", "ld_mismatch"),
  evolution = c("neutral", "stabilizing", "directional"),
  ancestry = c("matched", "mismatched"),
  causal_count = NULL,
  omitted_count = NULL,
  oracle_wgs = FALSE,
  infinitesimal_h2 = NULL,
  ld_mismatch_variance = NULL
)
```

## Arguments

- seed:

  Random seed; the same seed gives the same locus.

- p:

  Number of tested variants.

- q:

  Number of omitted variants.

- n:

  GWAS sample size.

- wgs_n:

  Size of the sequencing reference panel used to estimate LD
  fingerprints.

- architecture:

  Causal architecture: all causal variants tested (`"contained"`), some
  or all omitted (`"partially_omitted"`, `"entirely_omitted"`), no
  signal (`"null"`), diffuse polygenic signal (`"infinitesimal"`), or
  reference LD error (`"ld_mismatch"`).

- evolution:

  Selection regime used to draw effect sizes.

- ancestry:

  Whether the reference panel matches the GWAS ancestry.

- causal_count, omitted_count:

  Optional numbers of causal variants and of omitted causal variants;
  defaults depend on `architecture`.

- oracle_wgs:

  If `TRUE`, use population LD instead of a finite-panel estimate.

- infinitesimal_h2, ld_mismatch_variance:

  Optional polygenic variance and LD-mismatch variance; defaults depend
  on `architecture`.

## Value

A list with the arguments needed by
[`susieout_fit()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit.md)
(`z`, `n`, `R_analysis`, `wgs_ld`, `examined_ids`), the variant table
`variants`, full-locus `oracle` statistics, and the simulated `truth`.

## Examples

``` r
locus <- susieout_simulate_locus(1, architecture = "entirely_omitted")
locus$truth$causal
#> [1] "v47"
locus$truth$any_omitted
#> [1] TRUE
```
