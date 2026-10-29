.susieout_standardize_genotypes <- function(genotypes, name) {
    matrix <- as.matrix(genotypes)
    storage.mode(matrix) <- "double"
    if (!length(matrix) || nrow(matrix) < 3L || ncol(matrix) < 1L) {
        stop(name, " must have at least three samples and one variant",
             call. = FALSE)
    }
    for (column in seq_len(ncol(matrix))) {
        missing <- !is.finite(matrix[, column])
        if (all(missing)) stop(name, " contains an entirely missing variant",
                               call. = FALSE)
        if (any(missing)) {
            matrix[missing, column] <- mean(matrix[!missing, column])
        }
    }
    means <- colMeans(matrix)
    matrix <- sweep(matrix, 2L, means, `-`)
    scales <- sqrt(colSums(matrix^2) / (nrow(matrix) - 1L))
    if (any(!is.finite(scales) | scales <= 0)) {
        stop(name, " contains a monomorphic variant", call. = FALSE)
    }
    sweep(matrix, 2L, scales, `/`)
}

.susieout_factorized_genotype_cache <- function(examined_genotypes,
                                             outside_genotypes,
                                             tolerance = 1e-10) {
    examined <- .susieout_standardize_genotypes(
        examined_genotypes, "examined_genotypes")
    outside <- .susieout_standardize_genotypes(
        outside_genotypes, "outside_genotypes")
    if (nrow(examined) != nrow(outside)) {
        stop("tested and omitted genotypes must have the same samples",
             call. = FALSE)
    }
    denominator <- sqrt(nrow(examined) - 1L)
    G <- examined / denominator
    H <- outside / denominator
    decomposition <- svd(G, nu = min(dim(G)), nv = min(dim(G)))
    retained <- decomposition$d >
        max(decomposition$d) * max(tolerance, .Machine$double.eps^0.75)
    singular <- decomposition$d[retained]
    U <- decomposition$u[, retained, drop = FALSE]
    V <- decomposition$v[, retained, drop = FALSE]
    lambda <- singular^2
    projected_inside <- sweep(t(V), 1L, lambda, `*`)
    projected_outside <- sweep(
        crossprod(U, H), 1L, singular, `*`)
    list(
        sample_count = nrow(examined), p = ncol(examined),
        q = ncol(outside), V = V, lambda = lambda,
        projected_fingerprints = cbind(
            projected_inside, projected_outside),
        effective_rank = sum(lambda)^2 / sum(lambda^2))
}

.susieout_factorized_world_fit <- function(
        z, n, cache, examined_ids, excluded_ids, effect_prior, L,
        count_prior, infinitesimal_variance,
        ld_mismatch_variance, ridge, max_iter, tol, trained_prior) {
    z <- .susieout_validate_vector(z, "z", cache$p)
    examined_ids <- .susieout_validate_ids(examined_ids, cache$p, "variant_")
    excluded_ids <- .susieout_validate_ids(
        excluded_ids, cache$q, "outside_variant_")
    if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n <= 1) {
        stop("n must be a finite scalar greater than one", call. = FALSE)
    }
    if (!is.numeric(L) || length(L) != 1L || L < 1 || L != as.integer(L)) {
        stop("L must be a positive integer", call. = FALSE)
    }
    if (!is.null(trained_prior)) {
        if (!inherits(trained_prior, "susieout_prior") ||
            !isTRUE(trained_prior$frozen)) {
            stop("trained_prior must be a frozen susieout_prior", call. = FALSE)
        }
        if (is.null(trained_prior$count_prior)) {
            stop("trained_prior uses the removed component-prior schema",
                 call. = FALSE)
        }
        count_prior <- trained_prior$count_prior
        infinitesimal_variance <- trained_prior$infinitesimal_variance
        ld_mismatch_variance <- trained_prior$ld_mismatch_variance
        effect_prior <- trained_prior$effect_prior
    }
    count_prior <- .susieout_validate_count_prior(count_prior)
    mixture <- .susieout_validate_effect_prior(effect_prior)
    prior_masses <- .susieout_count_prior_masses(count_prior, L)
    state_prior <- c(
        rep(prior_masses$inside / cache$p, cache$p),
        rep(prior_masses$outside / cache$q, cache$q),
        prior_masses$null)

    multiplier <- .susieout_infinitesimal_multiplier(
        infinitesimal_variance, n, cache$p)
    subspace_variance <- pmax(
        cache$lambda + multiplier * cache$lambda^2 +
            ld_mismatch_variance + ridge,
        .Machine$double.eps)
    nullspace_variance <- ld_mismatch_variance + ridge
    projected <- cache$projected_fingerprints
    solved_projected <- sweep(projected, 1L, subspace_variance, `/`)
    quadratic <- colSums(projected * solved_projected)
    if (any(!is.finite(quadratic) | quadratic <= 0)) {
        stop("one or more factorized fingerprints have zero norm",
             call. = FALSE)
    }
    multiply <- function(beta) {
        as.numeric(cache$V %*% (projected %*% beta))
    }

    component_count <- as.integer(L)
    candidate_count <- cache$p + cache$q
    alpha <- matrix(0, component_count, candidate_count)
    component_null <- rep(1, component_count)
    beta_mean <- matrix(0, component_count, candidate_count)
    beta_second <- matrix(0, component_count, candidate_count)
    component_mixture <- matrix(
        0, component_count, length(mixture$weights))
    contributions <- matrix(0, cache$p, component_count)
    log_normalizers <- rep(NA_real_, component_count)
    component_kl <- rep(NA_real_, component_count)
    component_second_quadratic <- rep(NA_real_, component_count)
    elbo <- numeric(max_iter)
    sqrt_n <- sqrt(n)
    converged <- FALSE
    for (iteration in seq_len(max_iter)) {
        maximum_change <- 0
        total_contribution <- rowSums(contributions)
        for (component in seq_len(component_count)) {
            residual <- z - (total_contribution - contributions[, component])
            projected_residual <- as.numeric(crossprod(cache$V, residual))
            score <- as.numeric(crossprod(
                solved_projected, projected_residual))
            update <- .susieout_rank_one_update(
                residual, NULL, NULL, quadratic, sqrt_n, mixture,
                state_prior, score = score, multiply = multiply)
            old <- contributions[, component]
            contributions[, component] <- update$contribution
            total_contribution <- total_contribution - old + update$contribution
            maximum_change <- max(
                maximum_change, max(abs(update$contribution - old)))
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
        projected_residual <- as.numeric(crossprod(cache$V, residual))
        null_residual_squared <- max(
            sum(residual^2) - sum(projected_residual^2), 0)
        quadratic_residual <- sum(
            projected_residual^2 / subspace_variance) +
            null_residual_squared / nullspace_variance
        elbo[iteration] <- sum(log_normalizers) - 0.5 * quadratic_residual
        if (maximum_change < tol) {
            converged <- TRUE
            elbo <- elbo[seq_len(iteration)]
            break
        }
    }
    if (!converged) elbo <- elbo[seq_len(max_iter)]

    alpha_in <- alpha[, seq_len(cache$p), drop = FALSE]
    alpha_out <- alpha[, cache$p + seq_len(cache$q), drop = FALSE]
    colnames(alpha_in) <- examined_ids
    colnames(alpha_out) <- excluded_ids
    component_outside <- rowSums(alpha_out)
    any_omitted <- 1 - prod(1 - component_outside)
    pip_in <- .susieout_pip(alpha_in)
    pip_out <- .susieout_pip(alpha_out)
    names(pip_in) <- examined_ids
    names(pip_out) <- excluded_ids
    covariance_logdet <- sum(log(subspace_variance)) +
        (cache$p - length(cache$lambda)) * log(nullspace_variance)
    residual <- z - rowSums(contributions)
    projected_residual <- as.numeric(crossprod(cache$V, residual))
    null_residual_squared <- max(
        sum(residual^2) - sum(projected_residual^2), 0)
    residual_quadratic <- sum(projected_residual^2 / subspace_variance) +
        null_residual_squared / nullspace_variance
    projected_contribution <- crossprod(cache$V, contributions)
    contribution_quadratic <- colSums(
        sweep(projected_contribution^2, 1L, subspace_variance, `/`))
    posterior_variance_quadratic <- sum(pmax(
        component_second_quadratic - contribution_quadratic, 0))
    variational_evidence <- -0.5 * (
        cache$p * log(2 * pi) + covariance_logdet + residual_quadratic +
        posterior_variance_quadratic) - sum(component_kl)
    sampling_trace <- sum(cache$lambda)
    dense_trace <- multiplier * sum(cache$lambda^2)
    mismatch_trace <- cache$p * ld_mismatch_variance
    total_trace <- sampling_trace + dense_trace + mismatch_trace
    covariance_eigenvalues <- c(
        subspace_variance,
        rep(nullspace_variance, cache$p - length(cache$lambda)))
    fit <- list(
        call = match.call(), backend = "openworld_ibss_factorized",
        z = z, n = n, R_analysis = NULL,
        examined_ids = examined_ids, excluded_ids = excluded_ids,
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
        infinitesimal = list(
            local_h2 = infinitesimal_variance,
            variance = infinitesimal_variance,
            covariance_multiplier = multiplier,
            variance_fraction = dense_trace / total_trace),
        ld_mismatch = list(
            variance = ld_mismatch_variance,
            variance_fraction = mismatch_trace / total_trace),
        prior = c(list(count = count_prior, state = state_prior), prior_masses),
        effect_prior = mixture,
        fingerprints = list(
            factorized = TRUE,
            projected = projected,
            eigenvectors = cache$V,
            eigenvalues = cache$lambda),
        raw_outside_fingerprints = NULL,
        transformed_outside_fingerprints = NULL,
        diagnostics = list(
            converged = converged, iterations = length(elbo), elbo = elbo,
            ridge = ridge, covariance_logdet = covariance_logdet,
            variational_evidence = variational_evidence,
            component_kl = component_kl,
            posterior_variance_quadratic = posterior_variance_quadratic,
            nuisance_profile_score = utils::tail(elbo, 1L) -
                0.5 * covariance_logdet,
            effective_rank = cache$effective_rank,
            covariance_effective_rank =
                sum(covariance_eigenvalues)^2 /
                sum(covariance_eigenvalues^2),
            factorized_spectrum = list(
                lambda = cache$lambda, V = cache$V,
                projected_fingerprints = projected),
            wgs_transform_applied = FALSE,
            wgs_bootstrap = list(available = FALSE),
            maximum_responsibility_error = max(abs(
                rowSums(alpha) + component_null - 1)))
    )
    class(fit) <- "susieout_fit"
    fit
}

#' Fit the joint SuSiE-OUT nuisance model from a low-rank genotype reference
#'
#' This is an exact computational backend for LD matrices constructed from a
#' finite reference panel. It uses a single thin singular-value decomposition
#' of the reference genotypes while fitting the same null, tested-variant, OUT,
#' infinitesimal, and LD-mismatch states as [susieout_fit_nuisance_grid()].
#'
#' @param z,n,examined_ids,effect_prior,L,count_prior,trained_prior Arguments
#'   analogous to [susieout_fit_nuisance_grid()].
#' @param examined_genotypes Reference dosage matrix for tested variants,
#'   with samples in rows.
#' @param outside_genotypes Reference dosage matrix for omitted variants,
#'   with the same samples and ordering.
#' @param excluded_ids IDs for columns of `outside_genotypes`.
#' @param nuisance_prior Frozen grid prior from [susieout_train_nuisance_grid()].
#' @param ridge,max_iter,tol Numerical controls.
#' @return An `susieout_fit` averaged over the learned joint nuisance grid.
#' @export
susieout_fit_nuisance_grid_factorized <- function(
        z, n, examined_genotypes, outside_genotypes,
        examined_ids = NULL, excluded_ids = NULL,
        effect_prior = "robust_mixture", L = 10L,
        count_prior = susieout_count_prior(), nuisance_prior,
        trained_prior = NULL, ridge = 1e-6,
        max_iter = 200L, tol = 1e-5) {
    if (!inherits(nuisance_prior, "susieout_nuisance_prior") ||
        !isTRUE(nuisance_prior$frozen)) {
        stop("nuisance_prior must be a frozen susieout_nuisance_prior",
             call. = FALSE)
    }
    cache <- .susieout_factorized_genotype_cache(
        examined_genotypes, outside_genotypes)
    examined_ids <- .susieout_validate_ids(
        examined_ids %||% colnames(examined_genotypes), cache$p, "variant_")
    excluded_ids <- .susieout_validate_ids(
        excluded_ids %||% colnames(outside_genotypes),
        cache$q, "outside_variant_")
    grid <- nuisance_prior$grid[, c(
        "infinitesimal_h2", "ld_mismatch_variance", "prior_weight")]
    fits <- vector("list", nrow(grid))
    for (index in seq_len(nrow(grid))) {
        conditional_prior <- trained_prior
        if (!is.null(conditional_prior)) {
            conditional_prior$infinitesimal_variance <-
                grid$infinitesimal_h2[index]
            conditional_prior$ld_mismatch_variance <-
                grid$ld_mismatch_variance[index]
        }
        fits[[index]] <- .susieout_factorized_world_fit(
            z = z, n = n, cache = cache, examined_ids = examined_ids,
            excluded_ids = excluded_ids, effect_prior = effect_prior, L = L,
            count_prior = count_prior,
            infinitesimal_variance = grid$infinitesimal_h2[index],
            ld_mismatch_variance = grid$ld_mismatch_variance[index],
            ridge = ridge, max_iter = max_iter, tol = tol,
            trained_prior = conditional_prior)
    }
    grid$expected_nll <- vapply(seq_len(nrow(grid)), function(index)
        .susieout_nuisance_expected_nll(
            fits[[index]], grid$infinitesimal_h2[index],
            grid$ld_mismatch_variance[index]), numeric(1))
    grid$expected_loglik_score <- -0.5 * grid$expected_nll
    log_weight <- log(grid$prior_weight) + grid$expected_loglik_score
    log_normalizer <- .susieout_logsumexp(log_weight)
    grid$posterior_weight <- exp(log_weight - log_normalizer)
    grid$converged <- vapply(fits, function(fit)
        isTRUE(fit$diagnostics$converged), logical(1))
    map <- which.max(grid$posterior_weight)
    result <- fits[[map]]
    weighted <- function(extract) Reduce(`+`, Map(
        function(fit, weight) weight * extract(fit),
        fits, grid$posterior_weight))
    result$pip_in <- weighted(function(fit) fit$pip_in)
    result$pip_out <- weighted(function(fit) fit$pip_out)
    result$alpha_in <- weighted(function(fit) fit$alpha_in)
    result$alpha_out <- weighted(function(fit) fit$alpha_out)
    result$component_null <- weighted(function(fit) fit$component_null)
    result$component_outside <- weighted(function(fit) fit$component_outside)
    result$component_mixture <- weighted(function(fit) fit$component_mixture)
    result$component_contribution <- weighted(
        function(fit) fit$component_contribution)
    result$fitted <- rowSums(result$component_contribution)
    result$residual <- result$z - result$fitted
    alpha <- cbind(result$alpha_in, result$alpha_out)
    joint_mean <- weighted(function(fit)
        cbind(fit$alpha_in, fit$alpha_out) * fit$component_beta_mean)
    joint_second <- weighted(function(fit)
        cbind(fit$alpha_in, fit$alpha_out) * fit$component_beta_second)
    result$component_beta_mean <- ifelse(alpha > 0, joint_mean / alpha, 0)
    result$component_beta_second <- ifelse(
        alpha > 0, joint_second / alpha, 0)
    result$any_causal_omitted_raw <- weighted(
        function(fit) fit$any_causal_omitted)
    result$any_causal_omitted <- result$any_causal_omitted_raw
    calibration <- trained_prior$locus_calibration %||% NULL
    if (!is.null(calibration)) {
        raw <- pmin(pmax(result$any_causal_omitted_raw, 1e-8), 1 - 1e-8)
        result$any_causal_omitted <- stats::plogis(
            calibration$intercept + calibration$slope * stats::qlogis(raw))
        result$prior$locus_calibration <- calibration
    }
    result$candidate_set_complete <- 1 - result$any_causal_omitted
    result$null <- list(
        per_effect = result$component_null,
        all_effects = weighted(function(fit) fit$null$all_effects))
    result$infinitesimal <- list(
        local_h2 = sum(grid$posterior_weight * grid$infinitesimal_h2),
        variance = sum(grid$posterior_weight * grid$infinitesimal_h2),
        covariance_multiplier = weighted(
            function(fit) fit$infinitesimal$covariance_multiplier),
        variance_fraction = weighted(
            function(fit) fit$infinitesimal$variance_fraction),
        posterior_probability_nonzero = sum(
            grid$posterior_weight[grid$infinitesimal_h2 > 0]))
    result$ld_mismatch <- list(
        variance = sum(grid$posterior_weight * grid$ld_mismatch_variance),
        variance_fraction = weighted(
            function(fit) fit$ld_mismatch$variance_fraction),
        posterior_probability_nonzero = sum(
            grid$posterior_weight[grid$ld_mismatch_variance > 0]))
    result$backend <- "openworld_nuisance_grid_factorized"
    result$prior$nuisance_grid <- list(
        infinitesimal_prior = nuisance_prior$infinitesimal_marginal,
        ld_mismatch_prior = nuisance_prior$ld_mismatch_marginal)
    result$diagnostics$nuisance_grid <- grid
    result$diagnostics$nuisance_log_normalizer <- log_normalizer
    result$diagnostics$nuisance_map_index <- map
    result$diagnostics$converged <- all(
        grid$converged[grid$posterior_weight > 1e-6])
    class(result) <- "susieout_fit"
    result
}
