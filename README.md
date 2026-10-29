# SuSiEOut

`SuSiEOut` implements SuSiE-OUT, a sparse summary-statistic fine-mapping method
that allows causal variants to fall outside the tested variant set. Each single
effect can be null, assigned to a tested variant, or assigned to an untested
variant represented by its linkage disequilibrium (LD) pattern in an expanded
whole-genome sequencing (WGS) reference.

## Installation

```r
# install.packages("remotes")
remotes::install_github("weinstocklab/SuSiEOut")
```

## Example

```r
library(SuSiEOut)
locus <- susieout_simulate_locus(1, architecture = "entirely_omitted")
fit <- susieout_fit(
    locus$z, locus$n, locus$R_analysis, locus$wgs_ld,
    locus$examined_ids, effect_prior = "robust_mixture", L = 10,
    count_prior = susieout_count_prior(lambda_total = 1,
                                       out_fraction = 0.10))
summary(fit)
susieout_credible_sets(fit)
```

The dense nuisance parameter is total local infinitesimal variance, not the
coefficient placed directly on the z-score covariance. For `p` examined
variants and GWAS size `n`, `infinitesimal_variance = h2_dense` fits

```text
R + (n * h2_dense / p) R^2 + omega I.
```

This scaling lets the same local architecture be compared coherently across
GWAS sample sizes. Fixed fits report both `local_h2` and the resulting
`covariance_multiplier`; fixed variance fractions are not labeled posterior
probabilities.

The count prior is specified at the locus level. Changing `L` changes the
numerical approximation but not the prior probability that the locus contains
any sparse effect. `lambda_total` controls sparse-effect count and
`out_fraction` controls how that count is divided between tested (IN) and
excluded (OUT) candidates.

Heterogeneous banks should not impose one sparse count or one positive dense
variance on every locus. Spike-and-grid averaging combines sparse-count and
nuisance states using a comparable variational evidence:

```r
sparse_prior <- susieout_train_prior(training_loci, em_iterations = 0)
nuisance_prior <- susieout_train_nuisance_grid(
    training_loci, sparse_prior,
    infinitesimal_grid = c(0, 0.002, 0.008, 0.032),
    ld_mismatch_grid = c(0, 0.04, 0.12))
prior_grid <- susieout_prior_grid()
fit <- susieout_fit_grid(
    locus$z, locus$n, locus$R_analysis, locus$wgs_ld,
    locus$examined_ids, trained_prior = sparse_prior,
    prior_grid = prior_grid, nuisance_prior = nuisance_prior)
```

`susieout_train_prior(..., em_update = "nuisance_only")` provides a global EM
diagnostic without feeding fitted OUT responsibilities back into the count
prior. A zero `out_fraction` fits a candidate-complete model without WGS
cross-LD.

For large loci, `R_analysis` may be a `Matrix` sparse matrix. The sparse backend
avoids forming the full examined-by-examined fingerprint matrix and supports the
`R + dI` nuisance covariance (`infinitesimal_variance = 0`). The accompanying
scale audit covers 120 loci with 15,000 examined variants.
