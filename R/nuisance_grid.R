#' Fit a SuSiE-OUT model averaged over local nuisance states
#'
#' This post-lock extension gives each locus an explicit spike at zero dense
#' variance and positive infinitesimal/mismatch alternatives. Conditional fits
#' are averaged using a covariance-comparable posterior expected-likelihood
#' score, so a
#' global nuisance estimate need not be imposed on candidate-complete sparse
#' loci.
#'
#' @param z,n,R_analysis,wgs_ld,examined_ids,effect_prior,L,count_prior Arguments passed to
#'   [susieout_fit()].
#' @param infinitesimal_grid Grid of total local infinitesimal variances.
#' @param ld_mismatch_grid Grid of independent mismatch variances.
#' @param infinitesimal_prior Optional prior weights over the dense grid.
#' @param ld_mismatch_prior Optional prior weights over the mismatch grid.
#' @param nuisance_prior Optional learned joint grid prior returned by
#'   [susieout_train_nuisance_grid()].
#' @param trained_prior Optional frozen prior; its nuisance values are replaced
#'   by each grid point without modifying the original object.
#' @param ... Further arguments passed to [susieout_fit()].
#' @return An `susieout_fit` whose PIPs and locus quantities are averaged over the
#'   nuisance grid. `diagnostics$nuisance_grid` contains conditional scores and
#'   posterior weights.
#' @export
susieout_fit_nuisance_grid <- function(
    z,
    n,
    R_analysis,
    wgs_ld = NULL,
    examined_ids = NULL,
    effect_prior = "robust_mixture",
    L = 10L,
    count_prior = susieout_count_prior(),
    infinitesimal_grid = c(0, 0.002, 0.008, 0.032),
    ld_mismatch_grid = c(0, 0.04, 0.12),
    infinitesimal_prior = NULL,
    ld_mismatch_prior = NULL,
    nuisance_prior = NULL,
    trained_prior = NULL,
    ...
) {
    validate_grid <- function(values, name) {
        if (
            !is.numeric(values) ||
                !length(values) ||
                any(!is.finite(values)) ||
                any(values < 0) ||
                anyDuplicated(values)
        ) {
            stop(
                name,
                " must contain unique nonnegative finite values",
                call. = FALSE
            )
        }
        as.numeric(values)
    }
    default_spike_prior <- function(values) {
        if (length(values) == 1L) {
            return(1)
        }
        if (any(values == 0)) {
            weights <- rep(0.5 / (length(values) - 1L), length(values))
            weights[values == 0] <- 0.5
            return(weights)
        }
        rep(1 / length(values), length(values))
    }
    validate_prior <- function(weights, values, name) {
        weights <- weights %||% default_spike_prior(values)
        if (
            !is.numeric(weights) ||
                length(weights) != length(values) ||
                any(!is.finite(weights)) ||
                any(weights < 0) ||
                sum(weights) <= 0
        ) {
            stop(name, " must be nonnegative and match its grid", call. = FALSE)
        }
        weights / sum(weights)
    }
    if (!is.null(nuisance_prior)) {
        if (
            !inherits(nuisance_prior, "susieout_nuisance_prior") ||
                !isTRUE(nuisance_prior$frozen)
        ) {
            stop(
                "nuisance_prior must be a frozen susieout_nuisance_prior",
                call. = FALSE
            )
        }
        grid <- nuisance_prior$grid[, c(
            "infinitesimal_h2",
            "ld_mismatch_variance",
            "prior_weight"
        )]
        infinitesimal_grid <- unique(grid$infinitesimal_h2)
        ld_mismatch_grid <- unique(grid$ld_mismatch_variance)
        infinitesimal_prior <- nuisance_prior$infinitesimal_marginal
        ld_mismatch_prior <- nuisance_prior$ld_mismatch_marginal
    } else {
        infinitesimal_grid <- validate_grid(
            infinitesimal_grid,
            "infinitesimal_grid"
        )
        ld_mismatch_grid <- validate_grid(
            ld_mismatch_grid,
            "ld_mismatch_grid"
        )
        infinitesimal_prior <- validate_prior(
            infinitesimal_prior,
            infinitesimal_grid,
            "infinitesimal_prior"
        )
        ld_mismatch_prior <- validate_prior(
            ld_mismatch_prior,
            ld_mismatch_grid,
            "ld_mismatch_prior"
        )
        grid <- expand.grid(
            infinitesimal_h2 = infinitesimal_grid,
            ld_mismatch_variance = ld_mismatch_grid,
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
        )
        grid$prior_weight <- infinitesimal_prior[
            match(grid$infinitesimal_h2, infinitesimal_grid)
        ] *
            ld_mismatch_prior[
                match(grid$ld_mismatch_variance, ld_mismatch_grid)
            ]
    }
    fits <- vector("list", nrow(grid))
    for (index in seq_len(nrow(grid))) {
        conditional_prior <- trained_prior
        if (!is.null(conditional_prior)) {
            conditional_prior$infinitesimal_variance <-
                grid$infinitesimal_h2[index]
            conditional_prior$ld_mismatch_variance <-
                grid$ld_mismatch_variance[index]
        }
        arguments <- list(
            z = z,
            n = n,
            R_analysis = R_analysis,
            wgs_ld = wgs_ld,
            examined_ids = examined_ids,
            effect_prior = effect_prior,
            L = L,
            count_prior = count_prior,
            infinitesimal_variance = grid$infinitesimal_h2[index],
            ld_mismatch_variance = grid$ld_mismatch_variance[index],
            trained_prior = conditional_prior
        )
        fits[[index]] <- do.call(susieout_fit, c(arguments, list(...)))
    }
    grid$expected_nll <- vapply(
        seq_len(nrow(grid)),
        function(index) {
            .susieout_nuisance_expected_nll(
                fits[[index]],
                grid$infinitesimal_h2[index],
                grid$ld_mismatch_variance[index]
            )
        },
        numeric(1)
    )
    grid$expected_loglik_score <- -0.5 * grid$expected_nll
    log_weight <- log(grid$prior_weight) + grid$expected_loglik_score
    log_normalizer <- .susieout_logsumexp(log_weight)
    grid$posterior_weight <- exp(log_weight - log_normalizer)
    grid$converged <- vapply(
        fits,
        function(fit) {
            isTRUE(fit$diagnostics$converged)
        },
        logical(1)
    )
    map <- which.max(grid$posterior_weight)
    result <- fits[[map]]
    weighted <- function(extract) {
        Reduce(
            `+`,
            Map(
                function(fit, weight) weight * extract(fit),
                fits,
                grid$posterior_weight
            )
        )
    }
    result$pip_in <- weighted(function(fit) fit$pip_in)
    result$pip_out <- weighted(function(fit) fit$pip_out)
    result$alpha_in <- weighted(function(fit) fit$alpha_in)
    result$alpha_out <- weighted(function(fit) fit$alpha_out)
    result$component_null <- weighted(function(fit) fit$component_null)
    result$component_outside <- weighted(function(fit) fit$component_outside)
    result$component_mixture <- weighted(function(fit) fit$component_mixture)
    result$component_contribution <- weighted(
        function(fit) fit$component_contribution
    )
    result$fitted <- rowSums(result$component_contribution)
    result$residual <- result$z - result$fitted
    alpha <- cbind(result$alpha_in, result$alpha_out)
    joint_mean <- weighted(function(fit) {
        cbind(fit$alpha_in, fit$alpha_out) *
            fit$component_beta_mean
    })
    joint_second <- weighted(function(fit) {
        cbind(fit$alpha_in, fit$alpha_out) *
            fit$component_beta_second
    })
    result$component_beta_mean <- ifelse(alpha > 0, joint_mean / alpha, 0)
    result$component_beta_second <- ifelse(alpha > 0, joint_second / alpha, 0)
    result$any_causal_omitted_raw <- weighted(
        function(fit) fit$any_causal_omitted
    )
    result$any_causal_omitted <- result$any_causal_omitted_raw
    calibration <- trained_prior$locus_calibration %||% NULL
    if (!is.null(calibration)) {
        raw <- pmin(pmax(result$any_causal_omitted_raw, 1e-8), 1 - 1e-8)
        result$any_causal_omitted <- stats::plogis(
            calibration$intercept + calibration$slope * stats::qlogis(raw)
        )
        result$prior$locus_calibration <- calibration
    }
    result$candidate_set_complete <- 1 - result$any_causal_omitted
    result$null <- list(
        per_effect = result$component_null,
        all_effects = weighted(function(fit) fit$null$all_effects)
    )
    result$infinitesimal <- list(
        local_h2 = sum(grid$posterior_weight * grid$infinitesimal_h2),
        variance = sum(grid$posterior_weight * grid$infinitesimal_h2),
        covariance_multiplier = sum(
            grid$posterior_weight *
                vapply(
                    fits,
                    function(fit) fit$infinitesimal$covariance_multiplier,
                    numeric(1)
                )
        ),
        variance_fraction = sum(
            grid$posterior_weight *
                vapply(
                    fits,
                    function(fit) fit$infinitesimal$variance_fraction,
                    numeric(1)
                )
        ),
        posterior_probability_nonzero = sum(
            grid$posterior_weight[grid$infinitesimal_h2 > 0]
        )
    )
    result$ld_mismatch <- list(
        variance = sum(grid$posterior_weight * grid$ld_mismatch_variance),
        variance_fraction = sum(
            grid$posterior_weight *
                vapply(
                    fits,
                    function(fit) fit$ld_mismatch$variance_fraction,
                    numeric(1)
                )
        ),
        posterior_probability_nonzero = sum(
            grid$posterior_weight[grid$ld_mismatch_variance > 0]
        )
    )
    result$backend <- "openworld_nuisance_grid"
    result$prior$nuisance_grid <- list(
        infinitesimal_prior = infinitesimal_prior,
        ld_mismatch_prior = ld_mismatch_prior
    )
    result$diagnostics$nuisance_grid <- grid
    result$diagnostics$nuisance_log_normalizer <- log_normalizer
    result$diagnostics$nuisance_map_index <- map
    result$diagnostics$converged <- all(grid$converged[
        grid$posterior_weight > 1e-6
    ])
    class(result) <- "susieout_fit"
    result
}

#' Learn a spike-and-grid nuisance prior from training loci
#'
#' Conditional dense/mismatch likelihood scores are computed once per training
#' locus. A finite-mixture EM then estimates joint grid weights without using
#' architecture labels or nuisance truth.
#'
#' @param simulation_bank List of full loci labelled `training`.
#' @param trained_prior Frozen sparse/outside prior used in conditional fits.
#' @param infinitesimal_grid,ld_mismatch_grid Nuisance grids.
#' @param L Number of sparse single-effect components.
#' @param iterations Maximum mixture-weight EM iterations.
#' @param max_loci Maximum number of training loci.
#' @param dirichlet_strength Total pseudocount assigned according to the
#'   default spike prior.
#' @param tol Maximum grid-weight change for convergence.
#' @param seed Training-locus subsampling seed.
#' @param ... Additional [susieout_fit()] arguments.
#' @return A frozen `susieout_nuisance_prior` accepted by
#'   [susieout_fit_nuisance_grid()].
#' @export
susieout_train_nuisance_grid <- function(
    simulation_bank,
    trained_prior,
    infinitesimal_grid = c(0, 0.002, 0.008, 0.032),
    ld_mismatch_grid = c(0, 0.04, 0.12),
    L = 10L,
    iterations = 50L,
    max_loci = 500L,
    dirichlet_strength = 1,
    tol = 1e-6,
    seed = 1L,
    ...
) {
    if (
        !is.list(simulation_bank) ||
            !length(simulation_bank) ||
            !all(vapply(
                simulation_bank,
                function(locus) {
                    is.list(locus) &&
                        !is.null(locus$z) &&
                        !is.null(locus$R_analysis) &&
                        !is.null(locus$wgs_ld)
                },
                logical(1)
            ))
    ) {
        stop(
            "simulation_bank must be a nonempty list of full loci",
            call. = FALSE
        )
    }
    splits <- tolower(vapply(
        simulation_bank,
        function(locus) {
            as.character(locus$split %||% locus$truth$split %||% "training")
        },
        character(1)
    ))
    if (any(splits != "training")) {
        stop("nuisance-grid weights may use training loci only", call. = FALSE)
    }
    if (
        !inherits(trained_prior, "susieout_prior") ||
            !isTRUE(trained_prior$frozen)
    ) {
        stop("trained_prior must be a frozen susieout_prior", call. = FALSE)
    }
    validate <- function(values, name) {
        if (
            !is.numeric(values) ||
                !length(values) ||
                any(!is.finite(values)) ||
                any(values < 0) ||
                anyDuplicated(values)
        ) {
            stop(
                name,
                " must contain unique nonnegative finite values",
                call. = FALSE
            )
        }
        as.numeric(values)
    }
    default_prior <- function(values) {
        if (length(values) == 1L) {
            return(1)
        }
        if (any(values == 0)) {
            result <- rep(0.5 / (length(values) - 1L), length(values))
            result[values == 0] <- 0.5
            return(result)
        }
        rep(1 / length(values), length(values))
    }
    infinitesimal_grid <- validate(
        infinitesimal_grid,
        "infinitesimal_grid"
    )
    ld_mismatch_grid <- validate(
        ld_mismatch_grid,
        "ld_mismatch_grid"
    )
    grid <- expand.grid(
        infinitesimal_h2 = infinitesimal_grid,
        ld_mismatch_variance = ld_mismatch_grid,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    initial_weight <- default_prior(infinitesimal_grid)[match(
        grid$infinitesimal_h2,
        infinitesimal_grid
    )] *
        default_prior(ld_mismatch_grid)[match(
            grid$ld_mismatch_variance,
            ld_mismatch_grid
        )]
    if (!is.numeric(max_loci) || max_loci < 1L) {
        stop("max_loci must be positive", call. = FALSE)
    }
    set.seed(seed)
    loci <- simulation_bank
    if (length(loci) > max_loci) {
        loci <- loci[sample.int(length(loci), max_loci)]
    }
    score <- matrix(NA_real_, length(loci), nrow(grid))
    for (locus_index in seq_along(loci)) {
        locus <- loci[[locus_index]]
        for (grid_index in seq_len(nrow(grid))) {
            conditional_prior <- trained_prior
            conditional_prior$infinitesimal_variance <-
                grid$infinitesimal_h2[grid_index]
            conditional_prior$ld_mismatch_variance <-
                grid$ld_mismatch_variance[grid_index]
            fit <- susieout_fit(
                locus$z,
                locus$n,
                locus$R_analysis,
                locus$wgs_ld,
                locus$examined_ids,
                L = L,
                trained_prior = conditional_prior,
                ...
            )
            score[locus_index, grid_index] <- -0.5 *
                .susieout_nuisance_expected_nll(
                    fit,
                    grid$infinitesimal_h2[grid_index],
                    grid$ld_mismatch_variance[grid_index]
                )
        }
    }
    weights <- initial_weight / sum(initial_weight)
    trace <- vector("list", as.integer(iterations))
    converged <- FALSE
    completed <- 0L
    for (iteration in seq_len(as.integer(iterations))) {
        log_joint <- sweep(score, 2L, log(weights), `+`)
        log_normalizer <- apply(log_joint, 1L, .susieout_logsumexp)
        responsibility <- exp(log_joint - log_normalizer)
        updated <- (colSums(responsibility) +
            dirichlet_strength * initial_weight) /
            (nrow(score) + dirichlet_strength)
        maximum_change <- max(abs(updated - weights))
        weights <- updated / sum(updated)
        completed <- iteration
        trace[[iteration]] <- data.frame(
            iteration = iteration,
            training_log_score = sum(log_normalizer),
            maximum_weight_change = maximum_change,
            stringsAsFactors = FALSE
        )
        if (maximum_change <= tol) {
            converged <- TRUE
            break
        }
    }
    grid$prior_weight <- weights
    result <- list(
        grid = grid,
        infinitesimal_marginal = as.numeric(rowsum(
            grid$prior_weight,
            match(grid$infinitesimal_h2, infinitesimal_grid)
        )[, 1L]),
        ld_mismatch_marginal = as.numeric(rowsum(
            grid$prior_weight,
            match(grid$ld_mismatch_variance, ld_mismatch_grid)
        )[, 1L]),
        training_loci = length(loci),
        iterations = completed,
        converged = converged,
        trace = do.call(rbind, trace[seq_len(completed)]),
        trained_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
        frozen = TRUE
    )
    class(result) <- "susieout_nuisance_prior"
    result
}
