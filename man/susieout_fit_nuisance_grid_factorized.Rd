\name{susieout_fit_nuisance_grid_factorized}
\alias{susieout_fit_nuisance_grid_factorized}
\title{Fit the Joint SuSiE-OUT Nuisance Model from a Low-Rank Genotype Reference}
\usage{
susieout_fit_nuisance_grid_factorized(
  z,
  n,
  examined_genotypes,
  outside_genotypes,
  examined_ids = NULL,
  excluded_ids = NULL,
  effect_prior = "robust_mixture",
  L = 10L,
  count_prior = susieout_count_prior(),
  nuisance_prior,
  trained_prior = NULL,
  ridge = 1e-06,
  max_iter = 200L,
  tol = 1e-05
)
}
\arguments{
\item{z, n, examined_ids, effect_prior, L, count_prior, trained_prior}{Arguments analogous to \code{susieout_fit_nuisance_grid()}.}
\item{examined_genotypes}{Reference dosage matrix for examined variants, with samples in rows.}
\item{outside_genotypes}{Reference dosage matrix for excluded variants, with the same samples and ordering.}
\item{excluded_ids}{IDs for columns of \code{outside_genotypes}.}
\item{nuisance_prior}{Frozen grid prior from \code{susieout_train_nuisance_grid()}.}
\item{ridge, max_iter, tol}{Numerical controls.}
}
\value{An \code{susieout_fit} averaged over the learned joint nuisance grid.}
\description{
An exact computational backend for LD matrices constructed from a finite
reference panel. It uses one thin singular-value decomposition of the reference
genotypes while fitting the same null, inside, outside, infinitesimal, and
LD-mismatch states as \code{susieout_fit_nuisance_grid()}.
}
