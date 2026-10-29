\name{susieout_count_prior}
\alias{susieout_count_prior}
\alias{susieout_prior_grid}
\title{Construct L-invariant SuSiE-OUT count priors}
\usage{
susieout_count_prior(lambda_total = 1, out_fraction = 0.1)

susieout_prior_grid(
  lambda_total = c(0, 0.5, 1, 2, 4, 6),
  out_fraction = c(0, 0.1, 0.3, 0.6, 1),
  weights = NULL,
  dirichlet_strength = 5
)
}
\arguments{
\item{lambda_total}{Nonnegative locus-level sparse-effect rate or grid.}
\item{out_fraction}{Fraction of sparse effects assigned to OUT candidates.}
\item{weights}{Optional prior-grid weights.}
\item{dirichlet_strength}{Regularization strength for learned grid weights.}
}
\value{A validated count prior or frozen prior-grid object.}
\description{
Constructs sparse-effect priors whose locus-level meaning is independent of the
computational component limit \code{L}.
}
