#' Construct an L-invariant sparse-effect count prior
#'
#' `lambda_total` is the locus-level sparse-effect rate. `out_fraction` is the
#' expected fraction of sparse effects assigned to variants outside the tested
#' candidate set. The prior is converted to component probabilities only when
#' a model is fitted, so its locus-level interpretation does not change with
#' the computational component limit `L`.
#'
#' @param lambda_total Nonnegative sparse-effect rate for the locus.
#' @param out_fraction Fraction of sparse effects assigned to OUT candidates.
#' @return A `susieout_count_prior` object.
#' @export
susieout_count_prior <- function(lambda_total = 1, out_fraction = 0.10) {
    if (!is.numeric(lambda_total) || length(lambda_total) != 1L ||
        !is.finite(lambda_total) || lambda_total < 0) {
        stop("lambda_total must be a nonnegative finite scalar", call. = FALSE)
    }
    if (!is.numeric(out_fraction) || length(out_fraction) != 1L ||
        !is.finite(out_fraction) || out_fraction < 0 || out_fraction > 1) {
        stop("out_fraction must lie in [0, 1]", call. = FALSE)
    }
    if (lambda_total == 0) out_fraction <- 0
    structure(list(
        lambda_total = as.numeric(lambda_total),
        lambda_in = as.numeric(lambda_total * (1 - out_fraction)),
        lambda_out = as.numeric(lambda_total * out_fraction),
        out_fraction = as.numeric(out_fraction)),
        class = "susieout_count_prior")
}

.susieout_validate_count_prior <- function(count_prior) {
    if (!inherits(count_prior, "susieout_count_prior")) {
        stop("count_prior must be returned by susieout_count_prior()",
             call. = FALSE)
    }
    susieout_count_prior(count_prior$lambda_total, count_prior$out_fraction)
}

.susieout_count_prior_masses <- function(count_prior, L) {
    count_prior <- .susieout_validate_count_prior(count_prior)
    if (!is.numeric(L) || length(L) != 1L || !is.finite(L) ||
        L < 1 || L != as.integer(L)) {
        stop("L must be a positive integer", call. = FALSE)
    }
    null <- exp(-count_prior$lambda_total / L)
    nonnull <- 1 - null
    outside <- nonnull * count_prior$out_fraction
    inside <- nonnull - outside
    list(
        null = null, inside = inside, outside = outside,
        expected_total = L * nonnull,
        expected_inside = L * inside,
        expected_outside = L * outside,
        probability_any_sparse = 1 - null^L,
        probability_any_outside = 1 - (1 - outside)^L,
        L = as.integer(L))
}

#' Construct a spike-and-grid sparse-effect prior
#'
#' The default grid contains a null-locus spike and positive count-rate states
#' spanning candidate-complete through entirely omitted architectures.
#'
#' @param lambda_total Unique nonnegative total-effect rates.
#' @param out_fraction Unique OUT fractions in `[0, 1]`.
#' @param weights Optional weights for the resulting grid rows.
#' @param dirichlet_strength Regularization strength used when learning weights.
#' @return A frozen `susieout_prior_grid` object.
#' @export
susieout_prior_grid <- function(
        lambda_total = c(0, 0.5, 1, 2, 4, 6),
        out_fraction = c(0, 0.1, 0.3, 0.6, 1),
        weights = NULL, dirichlet_strength = 5) {
    validate_grid <- function(x, name, upper = Inf) {
        if (!is.numeric(x) || !length(x) || any(!is.finite(x)) ||
            any(x < 0) || any(x > upper) || anyDuplicated(x)) {
            stop(name, " must contain unique finite values in its allowed range",
                 call. = FALSE)
        }
        sort(as.numeric(x))
    }
    lambda_total <- validate_grid(lambda_total, "lambda_total")
    out_fraction <- validate_grid(out_fraction, "out_fraction", 1)
    positive <- expand.grid(
        lambda_total = lambda_total[lambda_total > 0],
        out_fraction = out_fraction,
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    grid <- if (any(lambda_total == 0)) {
        rbind(data.frame(lambda_total = 0, out_fraction = 0), positive)
    } else positive
    if (!nrow(grid)) stop("prior grid contains no states", call. = FALSE)
    if (is.null(weights)) {
        weights <- numeric(nrow(grid))
        null_state <- grid$lambda_total == 0
        positive_state <- !null_state
        weights[null_state] <- if (any(null_state)) 0.10 else 0
        remaining <- 1 - sum(weights)
        positive_lambdas <- unique(grid$lambda_total[positive_state])
        for (lambda in positive_lambdas) {
            index <- which(grid$lambda_total == lambda)
            zero_out <- index[grid$out_fraction[index] == 0]
            positive_out <- setdiff(index, zero_out)
            lambda_mass <- remaining / length(positive_lambdas)
            if (length(zero_out) && length(positive_out)) {
                weights[zero_out] <- 0.5 * lambda_mass / length(zero_out)
                weights[positive_out] <- 0.5 * lambda_mass /
                    length(positive_out)
            } else {
                weights[index] <- lambda_mass / length(index)
            }
        }
    }
    if (!is.numeric(weights) || length(weights) != nrow(grid) ||
        any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) {
        stop("weights must be nonnegative and match the prior grid", call. = FALSE)
    }
    if (!is.numeric(dirichlet_strength) || length(dirichlet_strength) != 1L ||
        !is.finite(dirichlet_strength) || dirichlet_strength < 0) {
        stop("dirichlet_strength must be nonnegative", call. = FALSE)
    }
    grid$prior_weight <- weights / sum(weights)
    structure(list(
        grid = grid,
        base_weight = grid$prior_weight,
        dirichlet_strength = as.numeric(dirichlet_strength),
        frozen = TRUE), class = "susieout_prior_grid")
}

.susieout_validate_prior_grid <- function(prior_grid) {
    if (!inherits(prior_grid, "susieout_prior_grid") ||
        !isTRUE(prior_grid$frozen) || !is.data.frame(prior_grid$grid)) {
        stop("prior_grid must be returned by susieout_prior_grid()",
             call. = FALSE)
    }
    required <- c("lambda_total", "out_fraction", "prior_weight")
    if (!all(required %in% names(prior_grid$grid))) {
        stop("prior_grid is missing required columns", call. = FALSE)
    }
    grid <- prior_grid$grid
    if (any(!is.finite(grid$lambda_total)) || any(grid$lambda_total < 0) ||
        any(!is.finite(grid$out_fraction)) ||
        any(grid$out_fraction < 0 | grid$out_fraction > 1) ||
        any(!is.finite(grid$prior_weight)) || any(grid$prior_weight < 0) ||
        sum(grid$prior_weight) <= 0 ||
        any(grid$lambda_total == 0 & grid$out_fraction != 0)) {
        stop("prior_grid contains invalid states or weights", call. = FALSE)
    }
    grid$prior_weight <- grid$prior_weight / sum(grid$prior_weight)
    prior_grid$grid <- grid
    prior_grid$base_weight <- prior_grid$base_weight %||% grid$prior_weight
    prior_grid$base_weight <- prior_grid$base_weight /
        sum(prior_grid$base_weight)
    prior_grid$dirichlet_strength <- prior_grid$dirichlet_strength %||% 5
    prior_grid
}
