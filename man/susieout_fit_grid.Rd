\name{susieout_fit_grid}
\alias{susieout_fit_grid}
\alias{susieout_train_prior_grid}
\alias{susieout_calibrate_prior_grid}
\title{Fit and train spike-and-grid SuSiE-OUT models}
\usage{
susieout_fit_grid(
  z, n, R_analysis, wgs_ld = NULL, examined_ids = NULL,
  effect_prior = "robust_mixture", L = 10L,
  prior_grid = susieout_prior_grid(), nuisance_prior = NULL,
  trained_prior = NULL, infinitesimal_variance = 0,
  ld_mismatch_variance = 0,
  starts = c("zero", "inside", "mixed"), keep_fits = FALSE, ...
)

susieout_train_prior_grid(
  simulation_bank, prior_grid = susieout_prior_grid(),
  nuisance_prior = NULL, trained_prior = NULL, L = 10L,
  iterations = 50L, anchor_strength = length(simulation_bank),
  tol = 1e-06, ...
)

susieout_calibrate_prior_grid(
  calibration_bank, prior_grid = susieout_prior_grid(),
  dirichlet_strength = NULL,
  iterations = 500L
)
}
\arguments{
\item{z, n, R_analysis, wgs_ld, examined_ids, effect_prior, L}{Locus inputs passed to \code{susieout_fit()}.}
\item{prior_grid}{A grid returned by \code{susieout_prior_grid()}.}
\item{nuisance_prior}{Optional frozen nuisance grid.}
\item{trained_prior}{Optional frozen effect and nuisance prior.}
\item{infinitesimal_variance, ld_mismatch_variance}{Singleton nuisance values.}
\item{starts}{Deterministic initialization strategies.}
\item{keep_fits}{Whether to retain conditional fits.}
\item{simulation_bank}{Full training loci.}
\item{calibration_bank}{Masked conditional fits and reference-PIP targets from unfiltered fits.}
\item{iterations, anchor_strength, dirichlet_strength, tol}{Grid-weight training and calibration controls.}
\item{...}{Additional arguments passed to \code{susieout_fit()}.}
}
\value{A model-averaged fit or learned frozen prior grid.}
\description{
Fits and learns spike-and-grid sparse-count priors using complete variational
evidence and globally anchored empirical-Bayes updates.
}
