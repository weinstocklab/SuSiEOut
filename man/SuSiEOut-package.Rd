\name{SuSiEOut}
\alias{SuSiEOut-package}
\alias{susieout_cli}
\alias{susieout_bootstrap_fingerprints}
\alias{susieout_candidate_mask}
\alias{susieout_candidate_support}
\alias{susieout_credible_sets}
\alias{susieout_effect_prior}
\alias{susieout_fit}
\alias{susieout_fit_nuisance_grid}
\alias{susieout_harmonize_cross_ld}
\alias{susieout_simulate_locus}
\alias{susieout_train_prior}
\alias{susieout_train_nuisance_grid}
\alias{susieout_transform_fingerprints}
\alias{plot.susieout_fit}
\alias{predict.susieout_fit}
\alias{print.susieout_fit}
\alias{summary.susieout_fit}
\title{SuSiE-OUT fine-mapping beyond the tested variant set}
\description{
Fits sparse summary-statistic fine-mapping models with explicit null, tested-variant,
and OUT states. OUT states represent variants omitted from the tested set by
their LD fingerprints in a sequencing reference panel. Locus-level count-rate priors preserve their
interpretation when the component limit changes. The package includes evolutionary
effect mixtures, ridge-stabilized reference-panel LD, prior training, diagnostics,
credible sets, simulation utilities, and a batch command-line interface.
}
\details{
See the package README and the accompanying manuscript for mathematical details
and reproducible examples. Large SuSiE-OUT loci may provide a
\code{Matrix} sparse analysis-LD matrix; this invokes the sparse
\code{R + dI} covariance backend without materializing all tested-variant
fingerprints as a dense matrix. Arguments named \code{examined_*} refer to
tested variants and \code{excluded_*} or \code{outside_*} to omitted variants.
}
