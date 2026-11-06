# Package index

## Fitting

- [`susieout_fit()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit.md)
  : Fit the SuSiE-OUT fine-mapping model
- [`susieout_fit_grid()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit_grid.md)
  : Fit SuSiE-OUT by averaging sparse-count and nuisance states
- [`susieout_fit_nuisance_grid()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit_nuisance_grid.md)
  : Fit a SuSiE-OUT model averaged over local nuisance states
- [`susieout_fit_nuisance_grid_factorized()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit_nuisance_grid_factorized.md)
  : Fit the joint SuSiE-OUT nuisance model from a low-rank genotype
  reference

## Priors

- [`susieout_count_prior()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_count_prior.md)
  : Construct an L-invariant sparse-effect count prior
- [`susieout_effect_prior()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_effect_prior.md)
  : Construct an effect-size mixture prior
- [`susieout_prior_grid()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_prior_grid.md)
  : Construct a spike-and-grid sparse-effect prior

## Training priors

- [`susieout_train_prior()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_train_prior.md)
  : Train and freeze L-invariant SuSiE-OUT hyperparameters
- [`susieout_train_prior_grid()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_train_prior_grid.md)
  : Learn global sparse-grid weights from retained-data evidence
- [`susieout_train_nuisance_grid()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_train_nuisance_grid.md)
  : Learn a spike-and-grid nuisance prior from training loci
- [`susieout_calibrate_prior_grid()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_calibrate_prior_grid.md)
  : Calibrate sparse-grid weights with masked-WGS targets

## Summarizing and plotting fits

- [`summary(`*`<susieout_fit>`*`)`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit-methods.md)
  [`print(`*`<summary.susieout_fit>`*`)`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit-methods.md)
  [`print(`*`<susieout_fit>`*`)`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit-methods.md)
  [`plot(`*`<susieout_fit>`*`)`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit-methods.md)
  [`predict(`*`<susieout_fit>`*`)`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_fit-methods.md)
  : Summarize, print, plot, and predict from a SuSiE-OUT fit
- [`susieout_locuszoom()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_locuszoom.md)
  : LocusZoom-style plot of a SuSiE-OUT fit
- [`susieout_candidate_support()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_candidate_support.md)
  : Variant support under SuSiE-OUT
- [`susieout_credible_sets()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_credible_sets.md)
  : SuSiE-OUT credible sets

## LD fingerprints

- [`susieout_candidate_mask()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_candidate_mask.md)
  : Apply realistic candidate-set masks
- [`susieout_harmonize_cross_ld()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_harmonize_cross_ld.md)
  : Harmonize WGS cross-LD rows to analysis alleles
- [`susieout_transform_fingerprints()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_transform_fingerprints.md)
  : Ridge-stabilized whitening and recoloring of WGS fingerprints
- [`susieout_bootstrap_fingerprints()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_bootstrap_fingerprints.md)
  : Bootstrap WGS cross-LD fingerprints

## Simulation and command line

- [`susieout_simulate_locus()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_simulate_locus.md)
  : Simulate a deterministic SuSiE-OUT locus
- [`susieout_cli()`](https://weinstocklab.github.io/SuSiEOut/reference/susieout_cli.md)
  : Command-line entry point for batch locus analysis
