.susieout_solve_chol <- function(chol_factor, x) {
    backsolve(chol_factor, forwardsolve(t(chol_factor), x))
}

.susieout_infinitesimal_multiplier <- function(local_h2, n, p) {
    n * local_h2 / p
}

.susieout_covariance_fractions <- function(R, infinitesimal_multiplier,
                                       ld_mismatch_variance) {
    sampling_trace <- sum(.susieout_diag(R))
    infinitesimal_trace <- infinitesimal_multiplier * sum(R^2)
    mismatch_trace <- nrow(R) * ld_mismatch_variance
    total_trace <- sampling_trace + infinitesimal_trace + mismatch_trace
    c(
        infinitesimal = infinitesimal_trace / total_trace,
        ld_mismatch = mismatch_trace / total_trace)
}

.susieout_rank_one_update <- function(residual, fingerprints, solved_fingerprints,
                                  quadratic, sqrt_n, mixture, state_prior,
                                  chol_factor = NULL, score = NULL,
                                  multiply = NULL) {
    if (is.null(score)) {
        score <- as.numeric(crossprod(solved_fingerprints, residual))
    }
    candidate_count <- length(quadratic)
    mixture_count <- length(mixture$weights)
    log_joint <- matrix(-Inf, candidate_count, mixture_count)
    log_bayes_factor <- matrix(0, candidate_count, mixture_count)
    posterior_mean <- matrix(0, candidate_count, mixture_count)
    posterior_variance <- matrix(0, candidate_count, mixture_count)
    for (m in seq_len(mixture_count)) {
        variance <- mixture$variances[m]
        prior_mean <- mixture$means[m]
        denominator <- 1 + sqrt_n^2 * variance * quadratic
        centered_score <- score - sqrt_n * prior_mean * quadratic
        log_bf <- -0.5 * log(denominator) +
            sqrt_n * prior_mean * score -
            0.5 * sqrt_n^2 * prior_mean^2 * quadratic +
            0.5 * sqrt_n^2 * variance * centered_score^2 / denominator
        log_joint[, m] <- log(state_prior[seq_len(candidate_count)]) +
            log(mixture$weights[m]) + log_bf
        log_bayes_factor[, m] <- log_bf
        posterior_mean[, m] <-
            (prior_mean + sqrt_n * variance * score) / denominator
        posterior_variance[, m] <- variance / denominator
    }
    null_log_weight <- log(state_prior[1L + candidate_count])
    log_normalizer <- .susieout_logsumexp(c(null_log_weight, as.vector(log_joint)))
    responsibility <- exp(log_joint - log_normalizer)
    null_probability <- exp(null_log_weight - log_normalizer)
    alpha <- rowSums(responsibility)
    joint_beta_mean <- rowSums(responsibility * posterior_mean)
    joint_beta_second <- rowSums(
        responsibility * (posterior_variance + posterior_mean^2))
    contribution <- sqrt_n * if (is.null(multiply)) {
        as.numeric(fingerprints %*% joint_beta_mean)
    } else {
        as.numeric(multiply(joint_beta_mean))
    }
    positive <- responsibility > 0
    categorical_kl <- sum(responsibility[positive] * (
        log(responsibility[positive]) -
        log_joint[positive] +
        log_bayes_factor[positive]))
    if (null_probability > 0) {
        categorical_kl <- categorical_kl + null_probability *
            (log(null_probability) - null_log_weight)
    }
    normal_kl <- 0
    for (m in seq_len(mixture_count)) {
        variance <- mixture$variances[m]
        prior_mean <- mixture$means[m]
        normal_term <- 0.5 * (
            (posterior_variance[, m] +
             (posterior_mean[, m] - prior_mean)^2) / variance - 1 +
            log(variance / posterior_variance[, m]))
        normal_kl <- normal_kl + sum(responsibility[, m] * normal_term)
    }
    expected_second_quadratic <- sqrt_n^2 * sum(
        responsibility * (posterior_variance + posterior_mean^2) *
            matrix(quadratic, candidate_count, mixture_count))
    list(
        alpha = alpha,
        null = null_probability,
        beta_mean = ifelse(alpha > 0, joint_beta_mean / alpha, 0),
        beta_second = ifelse(alpha > 0, joint_beta_second / alpha, 0),
        joint_beta_mean = joint_beta_mean,
        mixture_responsibility = colSums(responsibility),
        contribution = contribution,
        kl = categorical_kl + normal_kl,
        expected_second_quadratic = expected_second_quadratic,
        log_normalizer = log_normalizer,
        score = score
    )
}

.susieout_sparse_world_fit <- function(
        z, n, R_analysis, parsed_wgs, examined_ids, mixture, state_prior,
        L, ridge, max_iter, tol, infinitesimal_variance,
        ld_mismatch_variance, count_prior, prior_masses) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("sparse LD support requires the Matrix package", call. = FALSE)
    }
    if (infinitesimal_variance != 0) {
        stop(paste(
            "sparse LD currently requires infinitesimal_variance = 0;",
            "use a dense matrix for the R + tau R^2 covariance"), call. = FALSE)
    }
    if (ld_mismatch_variance > 1e-4) {
        stop(paste(
            "sparse LD currently requires ld_mismatch_variance <= 1e-4",
            "because the inside quadratic drops an order-d^2 inverse term"),
            call. = FALSE)
    }
    p <- length(examined_ids)
    q <- length(parsed_wgs$excluded_ids)
    cross <- parsed_wgs$transformed
    diagonal_shift <- ld_mismatch_variance + ridge
    covariance <- Matrix::forceSymmetric(
        R_analysis + Matrix::Diagonal(p, diagonal_shift), uplo = "U")
    factor <- tryCatch(
        Matrix::Cholesky(covariance, LDL = FALSE, perm = TRUE),
        error = function(error) NULL)
    if (is.null(factor)) stop("sparse nuisance covariance factorization failed",
                              call. = FALSE)
    solve_covariance <- function(rhs) {
        result <- Matrix::solve(factor, rhs, system = "A")
        if (is.null(dim(rhs))) as.numeric(result) else as.matrix(result)
    }
    # For A = R + dI, R A^-1 R = R - dI + d^2 A^-1. The omitted
    # d^2 diagonal term is below 1e-10 for the default numerical ridge in
    # well-conditioned LD blocks and avoids materializing a dense inverse.
    quadratic_inside <- pmax(as.numeric(.susieout_diag(R_analysis)) - diagonal_shift,
                             .Machine$double.eps)
    solved_outside <- solve_covariance(cross)
    quadratic_outside <- colSums(cross * solved_outside)
    quadratic <- c(quadratic_inside, quadratic_outside)
    multiply <- function(beta) {
        as.numeric(R_analysis %*% beta[seq_len(p)] +
            cross %*% beta[p + seq_len(q)])
    }

    component_count <- as.integer(L)
    candidate_count <- p + q
    alpha <- matrix(0, component_count, candidate_count)
    component_null <- rep(1, component_count)
    beta_mean <- matrix(0, component_count, candidate_count)
    beta_second <- matrix(0, component_count, candidate_count)
    component_mixture <- matrix(0, component_count, length(mixture$weights))
    contributions <- matrix(0, p, component_count)
    elbo <- numeric(max_iter)
    log_normalizers <- rep(NA_real_, component_count)
    component_kl <- rep(NA_real_, component_count)
    component_second_quadratic <- rep(NA_real_, component_count)
    sqrt_n <- sqrt(n)
    converged <- FALSE
    for (iteration in seq_len(max_iter)) {
        maximum_change <- 0
        total_contribution <- rowSums(contributions)
        for (component in seq_len(component_count)) {
            residual <- z - (total_contribution - contributions[, component])
            solved_residual <- solve_covariance(residual)
            score <- c(
                as.numeric(crossprod(R_analysis, solved_residual)),
                as.numeric(crossprod(cross, solved_residual)))
            update <- .susieout_rank_one_update(
                residual, NULL, NULL, quadratic, sqrt_n, mixture, state_prior,
                score = score, multiply = multiply)
            old <- contributions[, component]
            contributions[, component] <- update$contribution
            total_contribution <- total_contribution - old + update$contribution
            maximum_change <- max(maximum_change,
                                  max(abs(update$contribution - old)))
            alpha[component, ] <- update$alpha
            component_null[component] <- update$null
            beta_mean[component, ] <- update$beta_mean
            beta_second[component, ] <- update$beta_second
            component_mixture[component, ] <- update$mixture_responsibility
            log_normalizers[component] <- update$log_normalizer
            component_kl[component] <- update$kl
            component_second_quadratic[component] <-
                update$expected_second_quadratic
        }
        residual <- z - rowSums(contributions)
        solved_residual <- solve_covariance(residual)
        elbo[iteration] <- sum(log_normalizers) -
            0.5 * sum(residual * solved_residual)
        if (maximum_change < tol) {
            converged <- TRUE
            elbo <- elbo[seq_len(iteration)]
            break
        }
    }
    if (!converged) elbo <- elbo[seq_len(max_iter)]
    covariance_logdet <- as.numeric(
        Matrix::determinant(covariance, logarithm = TRUE)$modulus)
    residual <- z - rowSums(contributions)
    solved_residual <- solve_covariance(residual)
    solved_contributions <- solve_covariance(contributions)
    contribution_quadratic <- colSums(contributions * solved_contributions)
    posterior_variance_quadratic <- sum(pmax(
        component_second_quadratic - contribution_quadratic, 0))
    variational_evidence <- -0.5 * (
        p * log(2 * pi) + covariance_logdet +
        sum(residual * solved_residual) + posterior_variance_quadratic) -
        sum(component_kl)
    alpha_in <- alpha[, seq_len(p), drop = FALSE]
    alpha_out <- alpha[, p + seq_len(q), drop = FALSE]
    colnames(alpha_in) <- examined_ids
    colnames(alpha_out) <- parsed_wgs$excluded_ids
    component_outside <- rowSums(alpha_out)
    any_omitted <- 1 - prod(1 - component_outside)
    pip_in <- .susieout_pip(alpha_in)
    pip_out <- .susieout_pip(alpha_out)
    names(pip_in) <- examined_ids
    names(pip_out) <- parsed_wgs$excluded_ids
    fitted <- rowSums(contributions)
    covariance_fractions <- .susieout_covariance_fractions(
        R_analysis, 0, ld_mismatch_variance)
    fit <- list(
        call = match.call(), backend = "openworld_ibss_sparse",
        z = z, n = n, R_analysis = R_analysis,
        examined_ids = examined_ids, excluded_ids = parsed_wgs$excluded_ids,
        pip_in = pip_in, pip_out = pip_out,
        alpha_in = alpha_in, alpha_out = alpha_out,
        component_null = component_null,
        component_outside = component_outside,
        component_beta_mean = beta_mean,
        component_beta_second = beta_second,
        component_mixture = component_mixture,
        component_contribution = contributions,
        fitted = fitted, residual = z - fitted,
        any_causal_omitted = any_omitted,
        candidate_set_complete = 1 - any_omitted,
        null = list(per_effect = component_null,
                    all_effects = prod(component_null)),
        infinitesimal = list(local_h2 = 0, variance = 0,
                             covariance_multiplier = 0,
                             variance_fraction = 0),
        ld_mismatch = list(variance = ld_mismatch_variance,
                           variance_fraction = unname(
                               covariance_fractions["ld_mismatch"])),
        prior = c(list(count = count_prior, state = state_prior), prior_masses),
        effect_prior = mixture,
        fingerprints = list(inside = R_analysis, outside = cross),
        raw_outside_fingerprints = parsed_wgs$raw,
        transformed_outside_fingerprints = parsed_wgs$transformed,
        diagnostics = list(
            converged = converged, iterations = length(elbo), elbo = elbo,
            ridge = ridge,
            covariance_logdet = covariance_logdet,
            variational_evidence = variational_evidence,
            component_kl = component_kl,
            posterior_variance_quadratic = posterior_variance_quadratic,
            effective_rank = .susieout_participation_rank(R_analysis),
            covariance_effective_rank = .susieout_participation_rank(covariance),
            wgs_transform_applied = parsed_wgs$transform_applied,
            wgs_bootstrap = list(available = FALSE),
            sparse_quadratic_omitted_order = diagonal_shift^2,
            maximum_responsibility_error = max(abs(
                rowSums(alpha) + component_null - 1)))
    )
    class(fit) <- "susieout_fit"
    fit
}

.susieout_closed_world_fit <- function(z, n, R_analysis, examined_ids, L,
                                   dots = list()) {
    if (!requireNamespace("susieR", quietly = TRUE)) {
        stop("outside_prior = 0 requires the suggested susieR package",
             call. = FALSE)
    }
    arguments <- c(list(z = z, R = R_analysis, n = n, L = L,
                        estimate_residual_variance = FALSE), dots)
    susie_fit <- do.call(susieR::susie_rss, arguments)
    alpha <- as.matrix(susie_fit$alpha)
    colnames(alpha) <- examined_ids
    pip <- as.numeric(susie_fit$pip)
    names(pip) <- examined_ids
    fit <- list(
        call = match.call(), backend = "susieR_closed_world",
        z = z, n = n, R_analysis = R_analysis,
        examined_ids = examined_ids, excluded_ids = character(),
        pip_in = pip, pip_out = numeric(), alpha_in = alpha,
        alpha_out = matrix(numeric(), nrow(alpha), 0L),
        component_null = pmax(0, 1 - rowSums(alpha)),
        component_outside = rep(0, nrow(alpha)),
        component_beta_mean = susie_fit$mu,
        any_causal_omitted = 0, candidate_set_complete = 1,
        null = list(per_effect = pmax(0, 1 - rowSums(alpha)),
                    all_effects = prod(pmax(0, 1 - rowSums(alpha)))),
        infinitesimal = list(local_h2 = 0, variance = 0,
                             covariance_multiplier = 0,
                             variance_fraction = 0),
        ld_mismatch = list(variance = 0, variance_fraction = 0),
        diagnostics = list(
            converged = isTRUE(susie_fit$converged),
            iterations = length(susie_fit$elbo %||% numeric()),
            elbo = susie_fit$elbo %||% numeric(),
            effective_rank = .susieout_participation_rank(R_analysis),
            wgs_bootstrap = list(available = FALSE)),
        credible_sets_closed = susie_fit$sets,
        susie_fit = susie_fit,
        effect_prior = NULL,
        fingerprints = R_analysis
    )
    class(fit) <- "susieout_fit"
    fit
}

#' Fit the SuSiE-OUT fine-mapping model
#'
#' @param z Tested-variant z scores.
#' @param n GWAS sample size.
#' @param R_analysis Tested-variant analysis LD matrix. A `Matrix` sparse
#'   matrix invokes the scalable `R + dI` backend.
#' @param wgs_ld A cross-LD matrix or list containing `cross_ld`, optional
#'   `R_SS`, and variant IDs/alleles.
#' @param examined_ids IDs for tested variants.
#' @param effect_prior A prior returned by [susieout_effect_prior()] or its name.
#' @param L Maximum number of single-effect components.
#' @param count_prior An L-invariant prior returned by
#'   [susieout_count_prior()].
#' @param infinitesimal_variance Total local infinitesimal effect variance
#'   (local heritability). The corresponding z-score covariance multiplier is
#'   `n * infinitesimal_variance / length(z)`.
#' @param ld_mismatch_variance Fixed independent LD-mismatch variance.
#' @param transform_wgs Apply whitening/recoloring when WGS `R_SS` is supplied.
#' @param ridge Numerical ridge relative to the covariance scale.
#' @param max_iter Maximum IBSS sweeps.
#' @param tol Maximum fitted-mean change required for convergence.
#' @param initialization Variational initialization. `inside` greedily seeds
#'   tested variants; `mixed` seeds tested variants and OUT states together.
#' @param trained_prior Optional frozen object from [susieout_train_prior()].
#' @param inside_weights Optional within-set prior weights.
#' @param outside_weights Optional within-set prior weights.
#' @param ... Reserved; extra arguments are rejected.
#' @return An object of class `susieout_fit`.
#' @export
susieout_fit <- function(z, n, R_analysis, wgs_ld = NULL, examined_ids = NULL,
                     effect_prior = "robust_mixture", L = 10L,
                     count_prior = susieout_count_prior(),
                     infinitesimal_variance = 0,
                     ld_mismatch_variance = 0,
                     transform_wgs = TRUE, ridge = 1e-6,
                     max_iter = 200L, tol = 1e-5,
                     initialization = c("zero", "inside", "mixed"),
                     trained_prior = NULL,
                     inside_weights = NULL, outside_weights = NULL, ...) {
    z <- .susieout_validate_vector(z, "z")
    initialization <- match.arg(initialization)
    p <- length(z)
    examined_ids <- .susieout_validate_ids(examined_ids %||% names(z), p, "variant_")
    R_analysis <- .susieout_validate_matrix(
        R_analysis, "R_analysis", p, p, symmetric = TRUE)
    .susieout_check_psd(R_analysis, ridge)
    if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n <= 1) {
        stop("n must be a finite scalar greater than one", call. = FALSE)
    }
    if (!is.numeric(L) || length(L) != 1L || L < 1 || L != as.integer(L)) {
        stop("L must be a positive integer", call. = FALSE)
    }
    count_prior <- .susieout_validate_count_prior(count_prior)
    dots <- list(...)
    if (!is.null(trained_prior)) {
        if (!inherits(trained_prior, "susieout_prior") ||
            !isTRUE(trained_prior$frozen)) {
            stop("trained_prior must be a frozen susieout_prior", call. = FALSE)
        }
        if (is.null(trained_prior$count_prior)) {
            stop(paste(
                "trained_prior uses the removed component-prior schema;",
                "retrain it with the L-invariant count-prior API"), call. = FALSE)
        }
        count_prior <- .susieout_validate_count_prior(trained_prior$count_prior)
        infinitesimal_variance <- trained_prior$infinitesimal_variance
        ld_mismatch_variance <- trained_prior$ld_mismatch_variance
        effect_prior <- trained_prior$effect_prior
    }
    if (length(dots)) {
        stop("legacy or extra fit arguments are not supported", call. = FALSE)
    }
    for (entry in c("infinitesimal_variance", "ld_mismatch_variance")) {
        value <- get(entry)
        if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
            value < 0) stop(entry, " must be nonnegative", call. = FALSE)
    }
    mixture <- .susieout_validate_effect_prior(effect_prior)
    prior_masses <- .susieout_count_prior_masses(count_prior, L)
    parsed_wgs <- if (is.null(wgs_ld)) {
        if (prior_masses$outside > 0) {
            stop("a positive OUT rate requires WGS cross-LD", call. = FALSE)
        }
        empty <- matrix(numeric(), p, 0L,
                        dimnames = list(examined_ids, character()))
        list(raw = empty, transformed = empty, excluded_ids = character(),
             R_wgs_SS = NULL, transform_applied = FALSE, bootstrap = NULL)
    } else .susieout_parse_wgs_ld(
        wgs_ld, examined_ids, R_analysis, transform_wgs, ridge)
    q <- length(parsed_wgs$excluded_ids)
    if (q == 0L && prior_masses$outside > 0) {
        stop("a positive OUT rate requires omitted variants in the reference LD", call. = FALSE)
    }
    inside_weights <- inside_weights %||% rep(1, p)
    outside_weights <- outside_weights %||% rep(1, q)
    inside_weights <- .susieout_validate_vector(inside_weights, "inside_weights", p)
    if (q > 0L) {
        outside_weights <- .susieout_validate_vector(
            outside_weights, "outside_weights", q)
    }
    if (any(inside_weights < 0) || sum(inside_weights) <= 0 ||
        (q > 0L && (any(outside_weights < 0) || sum(outside_weights) <= 0))) {
        stop("within-set prior weights must be nonnegative with positive sums",
             call. = FALSE)
    }
    state_prior <- c(
        prior_masses$inside * inside_weights / sum(inside_weights),
        if (q > 0L) prior_masses$outside * outside_weights /
            sum(outside_weights) else numeric(),
        prior_masses$null)
    if (inherits(R_analysis, "sparseMatrix")) {
        if (initialization != "zero") {
            stop("sparse LD currently supports initialization = 'zero' only",
                 call. = FALSE)
        }
        return(.susieout_sparse_world_fit(
            z, n, R_analysis, parsed_wgs, examined_ids, mixture, state_prior,
            L, ridge, max_iter, tol, infinitesimal_variance,
            ld_mismatch_variance, count_prior, prior_masses))
    }
    fingerprints <- cbind(R_analysis, parsed_wgs$transformed)
    colnames(fingerprints) <- c(examined_ids, parsed_wgs$excluded_ids)
    infinitesimal_multiplier <- .susieout_infinitesimal_multiplier(
        infinitesimal_variance, n, p)
    covariance <- R_analysis +
        infinitesimal_multiplier * tcrossprod(R_analysis) +
        diag(ld_mismatch_variance + ridge, p)
    covariance <- (covariance + t(covariance)) / 2
    chol_factor <- tryCatch(chol(covariance), error = function(error) NULL)
    if (is.null(chol_factor)) {
        covariance <- covariance + diag(10 * ridge, p)
        chol_factor <- chol(covariance)
    }
    covariance_logdet <- 2 * sum(log(diag(chol_factor)))
    solved_fingerprints <- .susieout_solve_chol(chol_factor, fingerprints)
    quadratic <- colSums(fingerprints * solved_fingerprints)
    if (any(quadratic <= 0)) {
        stop("one or more LD fingerprints have zero covariance-adjusted norm",
             call. = FALSE)
    }
    component_count <- as.integer(L)
    candidate_count <- p + q
    alpha <- matrix(0, component_count, candidate_count)
    component_null <- rep(1, component_count)
    beta_mean <- matrix(0, component_count, candidate_count)
    beta_second <- matrix(0, component_count, candidate_count)
    component_mixture <- matrix(0, component_count, length(mixture$weights))
    contributions <- matrix(0, p, component_count)
    if (initialization != "zero") {
        eligible <- if (initialization == "inside") seq_len(p) else
            seq_len(candidate_count)
        greedy_residual <- z
        used <- integer()
        for (component in seq_len(min(component_count, length(eligible)))) {
            greedy_score <- as.numeric(crossprod(
                solved_fingerprints[, eligible, drop = FALSE],
                greedy_residual))
            greedy_score[eligible %in% used] <- 0
            selected_local <- which.max(abs(greedy_score) /
                sqrt(quadratic[eligible]))
            selected <- eligible[selected_local]
            if (!length(selected) || greedy_score[selected_local] == 0) break
            beta_seed <- greedy_score[selected_local] /
                (sqrt(n) * quadratic[selected])
            contributions[, component] <-
                sqrt(n) * fingerprints[, selected] * beta_seed
            greedy_residual <- greedy_residual - contributions[, component]
            used <- c(used, selected)
        }
    }
    elbo <- numeric(max_iter)
    log_normalizers <- rep(NA_real_, component_count)
    component_kl <- rep(NA_real_, component_count)
    component_second_quadratic <- rep(NA_real_, component_count)
    sqrt_n <- sqrt(n)
    converged <- FALSE
    for (iteration in seq_len(max_iter)) {
        maximum_change <- 0
        total_contribution <- rowSums(contributions)
        for (component in seq_len(component_count)) {
            residual <- z - (total_contribution - contributions[, component])
            update <- .susieout_rank_one_update(
                residual, fingerprints, solved_fingerprints, quadratic,
                sqrt_n, mixture, state_prior, chol_factor)
            old <- contributions[, component]
            contributions[, component] <- update$contribution
            total_contribution <- total_contribution - old + update$contribution
            maximum_change <- max(maximum_change,
                                  max(abs(update$contribution - old)))
            alpha[component, ] <- update$alpha
            component_null[component] <- update$null
            beta_mean[component, ] <- update$beta_mean
            beta_second[component, ] <- update$beta_second
            component_mixture[component, ] <- update$mixture_responsibility
            log_normalizers[component] <- update$log_normalizer
            component_kl[component] <- update$kl
            component_second_quadratic[component] <-
                update$expected_second_quadratic
        }
        residual <- z - rowSums(contributions)
        solved_residual <- .susieout_solve_chol(chol_factor, residual)
        elbo[iteration] <- sum(log_normalizers) -
            0.5 * sum(residual * solved_residual)
        if (maximum_change < tol) {
            converged <- TRUE
            elbo <- elbo[seq_len(iteration)]
            break
        }
    }
    if (!converged) elbo <- elbo[seq_len(max_iter)]
    residual <- z - rowSums(contributions)
    solved_residual <- .susieout_solve_chol(chol_factor, residual)
    solved_contributions <- .susieout_solve_chol(chol_factor, contributions)
    contribution_quadratic <- colSums(contributions * solved_contributions)
    posterior_variance_quadratic <- sum(pmax(
        component_second_quadratic - contribution_quadratic, 0))
    variational_evidence <- -0.5 * (
        p * log(2 * pi) + covariance_logdet +
        sum(residual * solved_residual) + posterior_variance_quadratic) -
        sum(component_kl)
    alpha_in <- alpha[, seq_len(p), drop = FALSE]
    alpha_out <- alpha[, p + seq_len(q), drop = FALSE]
    colnames(alpha_in) <- examined_ids
    colnames(alpha_out) <- parsed_wgs$excluded_ids
    component_outside <- rowSums(alpha_out)
    any_omitted <- 1 - prod(1 - component_outside)
    pip_in <- .susieout_pip(alpha_in)
    pip_out <- .susieout_pip(alpha_out)
    names(pip_in) <- examined_ids
    names(pip_out) <- parsed_wgs$excluded_ids
    bootstrap <- parsed_wgs$bootstrap
    bootstrap_diagnostics <- if (is.null(bootstrap)) {
        list(available = FALSE)
    } else {
        list(available = TRUE, B = bootstrap$B %||% length(bootstrap$cross_ld),
             seed = bootstrap$seed %||% NA_integer_,
             mean_fingerprint_sd = mean(bootstrap$sd %||% NA_real_, na.rm = TRUE))
    }
    covariance_fractions <- .susieout_covariance_fractions(
        R_analysis, infinitesimal_multiplier, ld_mismatch_variance)
    fit <- list(
        call = match.call(), backend = "openworld_ibss",
        z = z, n = n, R_analysis = R_analysis,
        examined_ids = examined_ids, excluded_ids = parsed_wgs$excluded_ids,
        pip_in = pip_in, pip_out = pip_out,
        alpha_in = alpha_in, alpha_out = alpha_out,
        component_null = component_null,
        component_outside = component_outside,
        component_beta_mean = beta_mean,
        component_beta_second = beta_second,
        component_mixture = component_mixture,
        component_contribution = contributions,
        fitted = rowSums(contributions), residual = z - rowSums(contributions),
        any_causal_omitted = any_omitted,
        candidate_set_complete = 1 - any_omitted,
        null = list(per_effect = component_null,
                    all_effects = prod(component_null)),
        infinitesimal = list(local_h2 = infinitesimal_variance,
                             variance = infinitesimal_variance,
                             covariance_multiplier = infinitesimal_multiplier,
                             variance_fraction = unname(
                                 covariance_fractions["infinitesimal"])),
        ld_mismatch = list(variance = ld_mismatch_variance,
                           variance_fraction = unname(
                               covariance_fractions["ld_mismatch"])),
        prior = c(list(count = count_prior, state = state_prior), prior_masses),
        effect_prior = mixture,
        fingerprints = fingerprints,
        raw_outside_fingerprints = parsed_wgs$raw,
        transformed_outside_fingerprints = parsed_wgs$transformed,
        diagnostics = list(
            converged = converged, iterations = length(elbo), elbo = elbo,
            ridge = ridge,
            covariance_logdet = covariance_logdet,
            variational_evidence = variational_evidence,
            component_kl = component_kl,
            posterior_variance_quadratic = posterior_variance_quadratic,
            nuisance_profile_score = utils::tail(elbo, 1L) -
                0.5 * covariance_logdet,
            effective_rank = .susieout_participation_rank(R_analysis),
            covariance_effective_rank = .susieout_participation_rank(covariance),
            wgs_transform_applied = parsed_wgs$transform_applied,
            wgs_bootstrap = bootstrap_diagnostics,
            maximum_responsibility_error = max(abs(
                rowSums(alpha) + component_null - 1)))
    )
    class(fit) <- "susieout_fit"
    fit
}
