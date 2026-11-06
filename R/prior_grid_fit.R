.susieout_average_fits <- function(fits, weights, grid, backend) {
    map <- which.max(weights)
    result <- fits[[map]]
    weighted <- function(extract) {
        Reduce(
            `+`,
            Map(
                function(fit, weight) weight * extract(fit),
                fits,
                weights
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
        cbind(fit$alpha_in, fit$alpha_out) * fit$component_beta_mean
    })
    joint_second <- weighted(function(fit) {
        cbind(fit$alpha_in, fit$alpha_out) * fit$component_beta_second
    })
    result$component_beta_mean <- ifelse(alpha > 0, joint_mean / alpha, 0)
    result$component_beta_second <- ifelse(alpha > 0, joint_second / alpha, 0)
    result$any_causal_omitted <- weighted(
        function(fit) fit$any_causal_omitted
    )
    result$candidate_set_complete <- 1 - result$any_causal_omitted
    result$null <- list(
        per_effect = result$component_null,
        all_effects = weighted(function(fit) fit$null$all_effects)
    )
    result$infinitesimal <- list(
        local_h2 = sum(weights * grid$infinitesimal_h2),
        variance = sum(weights * grid$infinitesimal_h2),
        covariance_multiplier = weighted(
            function(fit) fit$infinitesimal$covariance_multiplier
        ),
        variance_fraction = weighted(
            function(fit) fit$infinitesimal$variance_fraction
        ),
        posterior_probability_nonzero = sum(
            weights[grid$infinitesimal_h2 > 0]
        )
    )
    result$ld_mismatch <- list(
        variance = sum(weights * grid$ld_mismatch_variance),
        variance_fraction = weighted(
            function(fit) fit$ld_mismatch$variance_fraction
        ),
        posterior_probability_nonzero = sum(
            weights[grid$ld_mismatch_variance > 0]
        )
    )
    result$prior <- list(
        lambda_total_posterior_mean = sum(weights * grid$lambda_total),
        out_fraction_posterior_mean = sum(weights * grid$out_fraction),
        probability_positive_out_rate = sum(
            weights[grid$lambda_total > 0 & grid$out_fraction > 0]
        ),
        grid = grid
    )
    result$backend <- backend
    result$diagnostics$model_grid <- grid
    result$diagnostics$model_map_index <- map
    result$diagnostics$converged <- all(
        grid$converged[weights > 1e-6]
    )
    class(result) <- "susieout_fit"
    result
}

#' Fit SuSiE-OUT by averaging sparse-count and nuisance states
#'
#' Every sparse-count state is combined with every nuisance state. Posterior
#' state weights use the complete variational evidence stored by the
#' conditional fit, including likelihood normalization and variational
#' uncertainty.
#'
#' @param z,n,R_analysis,wgs_ld,examined_ids,effect_prior,L Arguments passed to
#'   [susieout_fit()].
#' @param prior_grid A grid returned by [susieout_prior_grid()].
#' @param nuisance_prior Optional frozen nuisance prior returned by
#'   [susieout_train_nuisance_grid()].
#' @param trained_prior Optional frozen effect/nuisance prior. Its count prior
#'   is replaced by each sparse-grid state.
#' @param infinitesimal_variance,ld_mismatch_variance Singleton nuisance values
#'   used when `nuisance_prior` is `NULL`.
#' @param keep_fits Retain conditional fits in diagnostics.
#' @param starts Deterministic initialization strategies. The highest-evidence
#'   converged fit is retained for each grid state.
#' @param ... Additional arguments passed to [susieout_fit()].
#' @return A model-averaged `susieout_fit`.
#' @export
susieout_fit_grid <- function(
    z,
    n,
    R_analysis,
    wgs_ld = NULL,
    examined_ids = NULL,
    effect_prior = "robust_mixture",
    L = 10L,
    prior_grid = susieout_prior_grid(),
    nuisance_prior = NULL,
    trained_prior = NULL,
    infinitesimal_variance = 0,
    ld_mismatch_variance = 0,
    starts = c("zero", "inside", "mixed"),
    keep_fits = FALSE,
    ...
) {
    prior_grid <- .susieout_validate_prior_grid(prior_grid)
    sparse_grid <- prior_grid$grid[, c(
        "lambda_total",
        "out_fraction",
        "prior_weight"
    )]
    names(sparse_grid)[3L] <- "sparse_prior_weight"
    if (is.null(nuisance_prior)) {
        nuisance_grid <- data.frame(
            infinitesimal_h2 = infinitesimal_variance,
            ld_mismatch_variance = ld_mismatch_variance,
            nuisance_prior_weight = 1
        )
    } else {
        if (
            !inherits(nuisance_prior, "susieout_nuisance_prior") ||
                !isTRUE(nuisance_prior$frozen)
        ) {
            stop(
                "nuisance_prior must be a frozen susieout_nuisance_prior",
                call. = FALSE
            )
        }
        nuisance_grid <- nuisance_prior$grid[, c(
            "infinitesimal_h2",
            "ld_mismatch_variance",
            "prior_weight"
        )]
        names(nuisance_grid)[3L] <- "nuisance_prior_weight"
    }
    grid <- merge(sparse_grid, nuisance_grid, by = NULL)
    grid$prior_weight <- grid$sparse_prior_weight *
        grid$nuisance_prior_weight
    fits <- vector("list", nrow(grid))
    for (index in seq_len(nrow(grid))) {
        conditional_prior <- trained_prior
        if (!is.null(conditional_prior)) {
            if (
                !inherits(conditional_prior, "susieout_prior") ||
                    !isTRUE(conditional_prior$frozen)
            ) {
                stop(
                    "trained_prior must be a frozen susieout_prior",
                    call. = FALSE
                )
            }
            conditional_prior$count_prior <- susieout_count_prior(
                grid$lambda_total[index],
                grid$out_fraction[index]
            )
            conditional_prior$infinitesimal_variance <-
                grid$infinitesimal_h2[index]
            conditional_prior$ld_mismatch_variance <-
                grid$ld_mismatch_variance[index]
        }
        available_starts <- if (inherits(R_analysis, "sparseMatrix")) {
            "zero"
        } else {
            unique(match.arg(
                starts,
                c("zero", "inside", "mixed"),
                several.ok = TRUE
            ))
        }
        start_fits <- lapply(available_starts, function(initialization) {
            susieout_fit(
                z = z,
                n = n,
                R_analysis = R_analysis,
                wgs_ld = wgs_ld,
                examined_ids = examined_ids,
                effect_prior = effect_prior,
                L = L,
                count_prior = susieout_count_prior(
                    grid$lambda_total[index],
                    grid$out_fraction[index]
                ),
                infinitesimal_variance = grid$infinitesimal_h2[index],
                ld_mismatch_variance = grid$ld_mismatch_variance[index],
                initialization = initialization,
                trained_prior = conditional_prior,
                ...
            )
        })
        evidence <- vapply(
            start_fits,
            function(fit) {
                fit$diagnostics$variational_evidence
            },
            numeric(1)
        )
        converged <- vapply(
            start_fits,
            function(fit) {
                isTRUE(fit$diagnostics$converged)
            },
            logical(1)
        )
        eligible <- if (any(converged)) {
            which(converged)
        } else {
            seq_along(start_fits)
        }
        selected <- eligible[which.max(evidence[eligible])]
        fits[[index]] <- start_fits[[selected]]
        selected_pip <- c(
            start_fits[[selected]]$pip_in,
            start_fits[[selected]]$pip_out
        )
        maximum_pip_difference <- vapply(
            start_fits,
            function(fit) {
                max(abs(c(fit$pip_in, fit$pip_out) - selected_pip))
            },
            numeric(1)
        )
        fits[[index]]$diagnostics$initialization_search <- data.frame(
            initialization = available_starts,
            variational_evidence = evidence,
            converged = converged,
            maximum_pip_difference_from_selected = maximum_pip_difference,
            selected = seq_along(start_fits) == selected
        )
    }
    grid$variational_evidence <- vapply(
        fits,
        function(fit) {
            fit$diagnostics$variational_evidence
        },
        numeric(1)
    )
    grid$converged <- vapply(
        fits,
        function(fit) {
            isTRUE(fit$diagnostics$converged)
        },
        logical(1)
    )
    log_weight <- log(grid$prior_weight) + grid$variational_evidence
    log_normalizer <- .susieout_logsumexp(log_weight)
    grid$posterior_weight <- exp(log_weight - log_normalizer)
    result <- .susieout_average_fits(
        fits,
        grid$posterior_weight,
        grid,
        "openworld_prior_grid"
    )
    result$diagnostics$model_log_normalizer <- log_normalizer
    if (isTRUE(keep_fits)) {
        result$diagnostics$conditional_fits <- fits
    }
    result
}

#' Learn global sparse-grid weights from retained-data evidence
#'
#' This sensitivity analysis updates only global grid weights. A Dirichlet
#' anchor prevents weakly identified retained-data likelihoods from replacing
#' the externally specified prior.
#'
#' @param simulation_bank Full training loci.
#' @param prior_grid Starting grid and anchor weights.
#' @param nuisance_prior Optional frozen nuisance grid.
#' @param trained_prior Optional frozen effect/nuisance prior.
#' @param L Component limit.
#' @param iterations Maximum mixture-weight updates.
#' @param anchor_strength Dirichlet effective sample size. By default it equals
#'   the number of training loci.
#' @param tol Weight convergence tolerance.
#' @param ... Additional [susieout_fit()] arguments.
#' @return A learned, frozen `susieout_prior_grid`.
#' @export
susieout_train_prior_grid <- function(
    simulation_bank,
    prior_grid = susieout_prior_grid(),
    nuisance_prior = NULL,
    trained_prior = NULL,
    L = 10L,
    iterations = 50L,
    anchor_strength = length(simulation_bank),
    tol = 1e-6,
    ...
) {
    if (!is.list(simulation_bank) || !length(simulation_bank)) {
        stop(
            "simulation_bank must be a nonempty full-locus list",
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
        stop("prior-grid weights may use training loci only", call. = FALSE)
    }
    prior_grid <- .susieout_validate_prior_grid(prior_grid)
    state_count <- nrow(prior_grid$grid)
    score <- matrix(NA_real_, length(simulation_bank), state_count)
    singleton_nuisance <- if (is.null(nuisance_prior)) NULL else nuisance_prior
    for (locus_index in seq_along(simulation_bank)) {
        locus <- simulation_bank[[locus_index]]
        for (state in seq_len(state_count)) {
            cell <- prior_grid
            cell$grid <- prior_grid$grid[state, , drop = FALSE]
            cell$grid$prior_weight <- 1
            fit <- susieout_fit_grid(
                locus$z,
                locus$n,
                locus$R_analysis,
                locus$wgs_ld,
                locus$examined_ids,
                L = L,
                prior_grid = cell,
                nuisance_prior = singleton_nuisance,
                trained_prior = trained_prior,
                ...
            )
            score[locus_index, state] <-
                fit$diagnostics$model_log_normalizer
        }
    }
    base <- prior_grid$base_weight
    weights <- prior_grid$grid$prior_weight
    trace <- data.frame()
    for (iteration in seq_len(as.integer(iterations))) {
        log_joint <- sweep(score, 2L, log(weights), `+`)
        normalizer <- apply(log_joint, 1L, .susieout_logsumexp)
        responsibility <- exp(log_joint - normalizer)
        updated <- (colSums(responsibility) + anchor_strength * base) /
            (nrow(score) + anchor_strength)
        maximum_change <- max(abs(updated - weights))
        weights <- updated / sum(updated)
        trace <- rbind(
            trace,
            data.frame(
                iteration = iteration,
                training_log_score = sum(normalizer),
                maximum_weight_change = maximum_change
            )
        )
        if (maximum_change <= tol) break
    }
    prior_grid$grid$prior_weight <- weights
    prior_grid$training_loci <- length(simulation_bank)
    prior_grid$anchor_strength <- anchor_strength
    prior_grid$training_trace <- trace
    prior_grid$trained_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
    prior_grid$frozen <- TRUE
    prior_grid
}

#' Calibrate sparse-grid weights with masked-WGS targets
#'
#' Conditional fits use masked statistics only. Reference PIPs from unfiltered fits enter only
#' this across-locus hyperparameter calibration, and callers must keep
#' calibration and evaluation loci disjoint.
#'
#' @param calibration_bank List of loci containing `evidence`, a
#'   variants-by-states `pip_in` matrix, and `target_pip`.
#' @param prior_grid Starting grid and Dirichlet center.
#' @param dirichlet_strength Effective strength of the grid-weight penalty.
#' @param iterations Optimizer iteration limit.
#' @return A learned `susieout_prior_grid` with calibration diagnostics.
#' @export
susieout_calibrate_prior_grid <- function(
    calibration_bank,
    prior_grid = susieout_prior_grid(),
    dirichlet_strength = NULL,
    iterations = 500L
) {
    prior_grid <- .susieout_validate_prior_grid(prior_grid)
    dirichlet_strength <- dirichlet_strength %||%
        prior_grid$dirichlet_strength %||%
        5
    state_count <- nrow(prior_grid$grid)
    if (!is.list(calibration_bank) || !length(calibration_bank)) {
        stop("calibration_bank must be a nonempty locus list", call. = FALSE)
    }
    validate_locus <- function(locus) {
        if (
            !is.list(locus) ||
                !is.numeric(locus$evidence) ||
                length(locus$evidence) != state_count ||
                !is.matrix(locus$pip_in) ||
                ncol(locus$pip_in) != state_count ||
                !is.numeric(locus$target_pip) ||
                length(locus$target_pip) != nrow(locus$pip_in) ||
                any(
                    !is.finite(c(
                        locus$evidence,
                        locus$pip_in,
                        locus$target_pip
                    ))
                )
        ) {
            stop(
                "each calibration locus must contain compatible finite evidence, pip_in, and target_pip",
                call. = FALSE
            )
        }
        if (
            any(locus$pip_in < 0 | locus$pip_in > 1) ||
                any(locus$target_pip < 0 | locus$target_pip > 1)
        ) {
            stop("calibration PIPs must lie in [0, 1]", call. = FALSE)
        }
        TRUE
    }
    invisible(lapply(calibration_bank, validate_locus))
    base <- prior_grid$base_weight
    if (length(base) != state_count) {
        base <- prior_grid$grid$prior_weight
    }
    base <- pmax(base / sum(base), 1e-12)
    softmax <- function(theta) {
        shifted <- theta - max(theta)
        exp(shifted) / sum(exp(shifted))
    }
    objective <- function(theta) {
        weights <- softmax(theta)
        locus_loss <- vapply(
            calibration_bank,
            function(locus) {
                log_state <- log(weights) + locus$evidence
                posterior <- exp(log_state - .susieout_logsumexp(log_state))
                predicted <- as.numeric(locus$pip_in %*% posterior)
                mean((predicted - locus$target_pip)^2)
            },
            numeric(1)
        )
        mean(locus_loss) -
            dirichlet_strength /
                length(calibration_bank) *
                sum(base * log(weights))
    }
    starts <- list(log(base), rep(0, state_count))
    solutions <- lapply(starts, function(start) {
        stats::optim(
            start,
            objective,
            method = "BFGS",
            control = list(maxit = as.integer(iterations), reltol = 1e-10)
        )
    })
    selected <- which.min(vapply(solutions, `[[`, numeric(1), "value"))
    solution <- solutions[[selected]]
    weights <- softmax(solution$par)
    prior_grid$grid$prior_weight <- weights
    prior_grid$calibration <- list(
        method = "locus_balanced_masked_wgs_pip_brier",
        loci = length(calibration_bank),
        dirichlet_strength = dirichlet_strength,
        objective = solution$value,
        optim_converged = identical(solution$convergence, 0L),
        starts = length(starts),
        trained_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
    prior_grid$frozen <- TRUE
    prior_grid
}
