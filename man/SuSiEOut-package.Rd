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
Fits sparse summary-statistic fine-mapping models with explicit null, tested,
and WGS-enumerated OUT states. Locus-level count-rate priors preserve their
interpretation when the component limit changes. The package includes evolutionary
effect mixtures, ridge-stabilized WGS LD transport, prior training, diagnostics,
credible sets, simulation utilities, and a batch command-line interface.
}
\details{
See the package README and the accompanying manuscript for mathematical details
and reproducible examples. Large SuSiE-OUT loci may provide a
\code{Matrix} sparse analysis-LD matrix; this invokes the sparse
\code{R + dI} covariance backend without materializing all inside
fingerprints as a dense matrix.
}
