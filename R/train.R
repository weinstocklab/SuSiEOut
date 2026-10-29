.susieout_training_frame <- function(simulation_bank) {
    if (is.data.frame(simulation_bank)) return(simulation_bank)
    if (!is.list(simulation_bank) || !length(simulation_bank)) {
        stop("simulation_bank must be a nonempty data frame or locus list",
             call. = FALSE)
    }
    rows <- lapply(simulation_bank, function(locus) {
        truth <- locus$truth %||% locus
        omitted_count <- truth$omitted_count %||% NULL
        if (is.null(omitted_count)) {
            causal_index <- truth$causal_index %||% integer()
            omitted_count <- if (length(causal_index) &&
                length(locus$examined_ids %||% character())) {
                sum(causal_index > length(locus$examined_ids))
            } else {
                as.numeric(truth$any_omitted %||% truth$omitted %||% FALSE)
            }
        }
        data.frame(
            omitted = as.logical(truth$any_omitted %||% truth$omitted),
            is_null = as.logical(truth$is_null %||% FALSE),
            split = as.character(locus$split %||% truth$split %||% "training"),
            infinitesimal_variance = as.numeric(
                truth$infinitesimal_variance %||% NA_real_),
            ld_mismatch_variance = as.numeric(
                truth$ld_mismatch_variance %||% NA_real_),
            causal_beta = as.numeric((truth$causal_beta %||% NA_real_)[1L]),
            causal_count = length(truth$causal %||% character()),
            omitted_count = as.numeric(omitted_count),
            architecture = as.character(truth$architecture %||% "unknown"),
            stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

.susieout_em_spectrum <- function(fit) {
    factorized <- fit$diagnostics$factorized_spectrum %||% NULL
    if (!is.null(factorized)) {
        projected_contributions <- crossprod(
            factorized$V, fit$component_contribution)
        alpha <- cbind(fit$alpha_in, fit$alpha_out)
        posterior_variance <- numeric(length(factorized$lambda))
        for (component in seq_len(nrow(alpha))) {
            second_weight <- alpha[component, ] *
                fit$component_beta_second[component, ]
            second_moment <- fit$n * as.numeric(
                (factorized$projected_fingerprints^2) %*% second_weight)
            posterior_variance <- posterior_variance + pmax(
                second_moment - projected_contributions[, component]^2, 0)
        }
        projected_residual <- as.numeric(crossprod(
            factorized$V, fit$residual))
        return(list(
            lambda = factorized$lambda,
            residual = projected_residual,
            posterior_variance = posterior_variance,
            residual_null_squared = max(
                sum(fit$residual^2) - sum(projected_residual^2), 0),
            null_dimension = length(fit$residual) -
                length(factorized$lambda),
            n = fit$n, p = length(fit$residual),
            ridge = fit$diagnostics$ridge %||% 1e-6))
    }
    decomposition <- eigen(as.matrix(fit$R_analysis), symmetric = TRUE)
    vectors <- decomposition$vectors
    fingerprints <- if (is.list(fit$fingerprints)) {
        cbind(as.matrix(fit$fingerprints$inside),
              as.matrix(fit$fingerprints$outside))
    } else {
        as.matrix(fit$fingerprints)
    }
    alpha <- cbind(fit$alpha_in, fit$alpha_out)
    projected_fingerprints <- crossprod(vectors, fingerprints)
    projected_contributions <- crossprod(
        vectors, fit$component_contribution)
    posterior_variance <- numeric(nrow(fingerprints))
    for (component in seq_len(nrow(alpha))) {
        second_weight <- alpha[component, ] *
            fit$component_beta_second[component, ]
        second_moment <- fit$n * as.numeric(
            (projected_fingerprints^2) %*% second_weight)
        posterior_variance <- posterior_variance + pmax(
            second_moment - projected_contributions[, component]^2, 0)
    }
    list(
        lambda = pmax(decomposition$values, 0),
        residual = as.numeric(crossprod(vectors, fit$residual)),
        posterior_variance = posterior_variance,
        n = fit$n,
        p = nrow(fit$R_analysis),
        ridge = fit$diagnostics$ridge %||% 1e-6)
}

.susieout_nuisance_expected_nll <- function(fit, local_h2, mismatch) {
    entry <- .susieout_em_spectrum(fit)
    multiplier <- .susieout_infinitesimal_multiplier(
        local_h2, entry$n, entry$p)
    variance <- pmax(entry$lambda + multiplier * entry$lambda^2 +
        mismatch + entry$ridge, 1e-8)
    expected_squared_residual <- entry$residual^2 +
        entry$posterior_variance
    result <- sum(log(variance) + expected_squared_residual / variance)
    if (!is.null(entry$null_dimension) && entry$null_dimension > 0L) {
        null_variance <- pmax(mismatch + entry$ridge, 1e-8)
        result <- result + entry$null_dimension * log(null_variance) +
            entry$residual_null_squared / null_variance
    }
    result
}

.susieout_em_nuisance <- function(fits, starting = c(0, 0)) {
    spectra <- lapply(fits, .susieout_em_spectrum)
    objective_parameters <- function(parameters) {
        local_h2 <- parameters[1L]
        mismatch <- parameters[2L]
        total <- 0
        dimension <- 0L
        for (entry in spectra) {
            multiplier <- .susieout_infinitesimal_multiplier(
                local_h2, entry$n, entry$p)
            variance <- pmax(entry$lambda + multiplier * entry$lambda^2 +
                mismatch + entry$ridge, 1e-8)
            expected_squared_residual <- entry$residual^2 +
                entry$posterior_variance
            total <- total + sum(log(variance) +
                expected_squared_residual / variance)
            dimension <- dimension + length(variance)
        }
        total / dimension
    }
    objective_theta <- function(theta) {
        objective_parameters(pmax(exp(theta) - 1, 0))
    }
    starting <- pmax(as.numeric(starting), 0)
    solution <- stats::optim(
        log1p(starting), objective_theta, method = "L-BFGS-B",
        lower = c(0, 0), upper = log1p(c(1, 1e4)))
    proposed <- pmax(exp(solution$par) - 1, 0)
    objective_before <- objective_parameters(starting)
    objective_after <- objective_parameters(proposed)
    damping <- 1
    while (objective_after > objective_before + 1e-10 && damping > 1 / 128) {
        damping <- damping / 2
        proposed <- starting + damping * (proposed - starting)
        objective_after <- objective_parameters(proposed)
    }
    list(
        parameters = proposed,
        objective_before = objective_before,
        objective_after = objective_after,
        damping = damping,
        optim_converged = identical(solution$convergence, 0L))
}

.susieout_calibration_delta <- function(fits, truth) {
    objective <- function(delta) {
        probability <- vapply(fits, function(fit) {
            component <- stats::plogis(stats::qlogis(pmin(pmax(
                fit$component_outside, 1e-12), 1 - 1e-12)) + delta)
            1 - prod(1 - component)
        }, numeric(1))
        probability <- pmin(pmax(probability, 1e-10), 1 - 1e-10)
        -mean(truth * log(probability) + (1 - truth) * log1p(-probability))
    }
    stats::optimize(objective, c(-8, 8), tol = 1e-6)$minimum
}

#' Train and freeze L-invariant SuSiE-OUT hyperparameters
#'
#' Sparse count rates are estimated once from known training counts. Optional
#' variational EM updates only effect-mixture weights and nuisance covariance;
#' it never feeds inferred OUT responsibilities back into the OUT count prior.
#'
#' @param simulation_bank Training/validation data frame or full-locus list.
#' @param effect_prior Starting effect-size mixture.
#' @param L Component limit used by the fitted training models.
#' @param em_iterations Number of effect/nuisance EM iterations.
#' @param em_max_loci Maximum number of training loci used by EM.
#' @param em_update Update effect weights and nuisance parameters or nuisance
#'   parameters only.
#' @param initial_infinitesimal_variance Initial total local dense variance.
#' @param initial_ld_mismatch_variance Initial independent mismatch variance.
#' @param em_tol Relative parameter convergence tolerance.
#' @param em_pip_tol Maximum PIP-change convergence tolerance.
#' @param seed Training-locus subsampling seed.
#' @param verbose Report EM progress.
#' @return A frozen `susieout_prior` accepted by [susieout_fit()].
#' @export
susieout_train_prior <- function(
        simulation_bank, effect_prior = "robust_mixture", L = 10L,
        em_iterations = 0L, em_max_loci = 500L,
        em_update = c("effect_and_nuisance", "nuisance_only"),
        initial_infinitesimal_variance = 0,
        initial_ld_mismatch_variance = 0,
        em_tol = 1e-4, em_pip_tol = 1e-3,
        seed = 1L, verbose = FALSE) {
    em_update <- match.arg(em_update)
    frame <- .susieout_training_frame(simulation_bank)
    split <- tolower(as.character(frame$split %||% "training"))
    if (any(split %in% c("test", "locked", "locked_test", "heldout",
                         "held_out"))) {
        stop("hyperparameters cannot be fit using held-out/locked test loci",
             call. = FALSE)
    }
    parameter_frame <- frame[split == "training", , drop = FALSE]
    if (!nrow(parameter_frame)) stop("no eligible training loci", call. = FALSE)
    if (!"causal_count" %in% names(parameter_frame) ||
        any(!is.finite(parameter_frame$causal_count)) ||
        !"omitted_count" %in% names(parameter_frame) ||
        any(!is.finite(parameter_frame$omitted_count))) {
        stop("training data must supply causal_count and omitted_count truth",
             call. = FALSE)
    }
    if (any(parameter_frame$causal_count < 0) ||
        any(parameter_frame$omitted_count < 0) ||
        any(parameter_frame$omitted_count > parameter_frame$causal_count)) {
        stop("training effect counts are inconsistent", call. = FALSE)
    }
    lambda_total <- mean(parameter_frame$causal_count)
    out_fraction <- if (sum(parameter_frame$causal_count) > 0) {
        sum(parameter_frame$omitted_count) /
            sum(parameter_frame$causal_count)
    } else 0
    count_prior <- susieout_count_prior(lambda_total, out_fraction)
    starting_prior <- .susieout_validate_effect_prior(effect_prior)
    observed_beta <- parameter_frame$causal_beta[
        is.finite(parameter_frame$causal_beta)]
    if (length(observed_beta) >= 10L) {
        empirical_variance <- mean(observed_beta^2)
        prior_variance <- sum(starting_prior$weights *
            (starting_prior$variances + starting_prior$means^2))
        scale <- min(max(empirical_variance / prior_variance, 0.25), 4)
        starting_prior$variances <- starting_prior$variances * scale
    }
    infinitesimal <- initial_infinitesimal_variance
    mismatch <- initial_ld_mismatch_variance
    em_trace <- data.frame()
    em_converged <- FALSE
    em_completed <- 0L
    training_loci <- list()
    if (em_iterations > 0L) {
        locus_bank <- is.list(simulation_bank) && length(simulation_bank) &&
            all(vapply(simulation_bank, function(locus)
                is.list(locus) && !is.null(locus$z) &&
                    !is.null(locus$R_analysis), logical(1)))
        if (!locus_bank) {
            stop("positive em_iterations requires full simulated loci",
                 call. = FALSE)
        }
        locus_split <- tolower(vapply(simulation_bank, function(locus)
            as.character(locus$split %||% locus$truth$split %||% "training"),
            character(1)))
        training_loci <- simulation_bank[locus_split == "training"]
        set.seed(seed)
        if (length(training_loci) > em_max_loci) {
            training_loci <- training_loci[sample.int(
                length(training_loci), em_max_loci)]
        }
        previous_pips <- NULL
        for (iteration in seq_len(as.integer(em_iterations))) {
            fits <- lapply(training_loci, function(locus) susieout_fit(
                locus$z, locus$n, locus$R_analysis, locus$wgs_ld,
                locus$examined_ids, effect_prior = starting_prior, L = L,
                count_prior = count_prior,
                infinitesimal_variance = infinitesimal,
                ld_mismatch_variance = mismatch))
            old_parameters <- c(starting_prior$weights, infinitesimal, mismatch)
            if (em_update == "effect_and_nuisance") {
                mixture_count <- Reduce(`+`, lapply(fits, function(fit)
                    colSums(fit$component_mixture)))
                starting_prior$weights <- (mixture_count + 1) /
                    sum(mixture_count + 1)
            }
            nuisance <- .susieout_em_nuisance(fits, c(infinitesimal, mismatch))
            infinitesimal <- nuisance$parameters[1L]
            mismatch <- nuisance$parameters[2L]
            current_pips <- unlist(lapply(fits, function(fit)
                c(fit$pip_in, fit$pip_out)), use.names = FALSE)
            maximum_pip_change <- if (is.null(previous_pips)) Inf else
                max(abs(current_pips - previous_pips))
            previous_pips <- current_pips
            new_parameters <- c(starting_prior$weights, infinitesimal, mismatch)
            maximum_relative_parameter_change <- max(
                abs(new_parameters - old_parameters) /
                    pmax(abs(old_parameters), 1e-4))
            em_trace <- rbind(em_trace, data.frame(
                iteration = iteration,
                infinitesimal_variance = infinitesimal,
                ld_mismatch_variance = mismatch,
                nuisance_objective_before = nuisance$objective_before,
                nuisance_objective_after = nuisance$objective_after,
                maximum_pip_change = maximum_pip_change,
                maximum_relative_parameter_change =
                    maximum_relative_parameter_change))
            em_completed <- iteration
            if (isTRUE(verbose)) message("SuSiE-OUT EM iteration ", iteration)
            if (iteration >= 2L && maximum_pip_change <= em_pip_tol &&
                maximum_relative_parameter_change <= em_tol) {
                em_converged <- TRUE
                break
            }
        }
    }
    result <- list(
        count_prior = count_prior,
        count_prior_training = list(
            mean_total_effects = lambda_total,
            mean_outside_effects = mean(parameter_frame$omitted_count),
            total_effects = sum(parameter_frame$causal_count),
            outside_effects = sum(parameter_frame$omitted_count)),
        prior_at_training_L = .susieout_count_prior_masses(count_prior, L),
        training_L = as.integer(L),
        infinitesimal_variance = infinitesimal,
        ld_mismatch_variance = mismatch,
        effect_prior = starting_prior,
        training_loci = nrow(parameter_frame),
        training_splits = "training",
        em_iterations = as.integer(em_iterations),
        em_iterations_completed = em_completed,
        em_converged = if (em_iterations > 0L) em_converged else NA,
        em_update = em_update,
        em_training_loci = length(training_loci),
        em_trace = em_trace,
        trained_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
        frozen = TRUE)
    class(result) <- "susieout_prior"
    result
}
