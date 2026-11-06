# Fit the SuSiE-OUT fine-mapping model

Fit the SuSiE-OUT fine-mapping model

## Usage

``` r
susieout_fit(
  z,
  n,
  R_analysis,
  wgs_ld = NULL,
  examined_ids = NULL,
  effect_prior = "robust_mixture",
  L = 10L,
  count_prior = susieout_count_prior(),
  infinitesimal_variance = 0,
  ld_mismatch_variance = 0,
  transform_wgs = TRUE,
  ridge = 1e-06,
  max_iter = 200L,
  tol = 1e-05,
  initialization = c("zero", "inside", "mixed"),
  trained_prior = NULL,
  inside_weights = NULL,
  outside_weights = NULL,
  ...
)
```

## Arguments

- z:

  Tested-variant z scores.

- n:

  GWAS sample size.

- R_analysis:

  Tested-variant analysis LD matrix. A `Matrix` sparse matrix invokes
  the scalable `R + dI` backend.

- wgs_ld:

  A cross-LD matrix or list containing `cross_ld`, optional `R_SS`, and
  variant IDs/alleles.

- examined_ids:

  IDs for tested variants.

- effect_prior:

  A prior returned by
  [`susieout_effect_prior()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_effect_prior.md)
  or its name.

- L:

  Maximum number of single-effect components.

- count_prior:

  An L-invariant prior returned by
  [`susieout_count_prior()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_count_prior.md).

- infinitesimal_variance:

  Total local infinitesimal effect variance (local heritability). The
  corresponding z-score covariance multiplier is
  `n * infinitesimal_variance / length(z)`.

- ld_mismatch_variance:

  Fixed independent LD-mismatch variance.

- transform_wgs:

  Apply whitening/recoloring when WGS `R_SS` is supplied.

- ridge:

  Numerical ridge relative to the covariance scale.

- max_iter:

  Maximum IBSS sweeps.

- tol:

  Maximum fitted-mean change required for convergence.

- initialization:

  Variational initialization. `inside` greedily seeds tested variants;
  `mixed` seeds tested variants and OUT states together.

- trained_prior:

  Optional frozen object from
  [`susieout_train_prior()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_train_prior.md).

- inside_weights:

  Optional within-set prior weights.

- outside_weights:

  Optional within-set prior weights.

- ...:

  Reserved; extra arguments are rejected.

## Value

An object of class `susieout_fit`, a list that includes:

- `pip_in`, `pip_out`: PIPs for tested variants and for omitted variants
  (OUT PIPs).

- `alpha_in`, `alpha_out`: per-component assignment probabilities.

- `any_causal_omitted`: posterior probability that at least one sparse
  effect is assigned to an omitted variant (locus OUT support), and
  `candidate_set_complete`, its complement.

- `component_null`, `component_outside`: per-component null and OUT
  probabilities.

- `fitted`, `residual`: fitted and residual z scores.

- `infinitesimal`, `ld_mismatch`: nuisance variance settings.

- `diagnostics`: convergence and numerical diagnostics.

Use
[`susieout_candidate_support()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_candidate_support.md),
[`susieout_credible_sets()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_credible_sets.md)
and
[summary()](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit-methods.md)
to summarize it.

## References

Wang G, Sarkar A, Carbonetto P, Stephens M (2020). A simple new approach
to variable selection in regression, with application to genetic fine
mapping. *Journal of the Royal Statistical Society: Series B* 82,
1273–1300.

## Examples

``` r
locus <- susieout_simulate_locus(1, architecture = "entirely_omitted")
fit <- susieout_fit(
  locus$z, locus$n, locus$R_analysis, locus$wgs_ld, locus$examined_ids,
  L = 10, count_prior = susieout_count_prior(1, out_fraction = 0.10)
)
fit$any_causal_omitted
#> [1] 0.1669364
head(susieout_candidate_support(fit))
#>     variant_id   state        pip maximum_component_probability
#> v40        v40  inside 0.04196349                   0.004277781
#> v28        v28  inside 0.03861921                   0.003930730
#> v29        v29  inside 0.02933064                   0.002972515
#> v41        v41 outside 0.02656404                   0.002688715
#> v44        v44 outside 0.02528099                   0.002557339
#> v43        v43 outside 0.02146891                   0.002167927
```
