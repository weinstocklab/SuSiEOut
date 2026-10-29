test_that("a known outside fingerprint is recovered", {
    R <- diag(4)
    C <- matrix(c(0.50, -0.45, 0.40, -0.35), 4, 1,
                dimnames = list(paste0("i", 1:4), "outside_truth"))
    n <- 1000000
    beta <- 0.04
    z <- sqrt(n) * C[, 1] * beta
    fit <- susieout_fit(z, n, R, C, paste0("i", 1:4),
                    effect_prior = list(weights = 1, means = 0,
                                        variances = 0.0025),
                    L = 1, count_prior = susieout_count_prior(-log(0.5), 0.5))
    expect_gt(unname(fit$pip_out["outside_truth"]), 0.9)
    expect_gt(fit$any_causal_omitted, 0.9)
    expect_equal(susieout_credible_sets(fit)[[1]]$state[1], "outside")
})

test_that("mixed inside and outside effects can occupy separate components", {
    R <- diag(6)
    C <- matrix(c(0.45, -0.35, 0.30, -0.25, 0.20, -0.15), 6, 1,
                dimnames = list(paste0("i", 1:6), "o1"))
    z <- sqrt(1000000) * (R[, 2] * 0.035 + C[, 1] * 0.045)
    fit <- susieout_fit(z, 1000000, R, C, paste0("i", 1:6),
                    effect_prior = list(weights = 1, means = 0,
                                        variances = 0.0025),
                    L = 2, count_prior = susieout_count_prior(-2 * log(0.4), 0.5))
    expect_gt(unname(fit$pip_in["i2"]), 0.8)
    expect_gt(unname(fit$pip_out["o1"]), 0.8)
    expect_true(fit$diagnostics$converged)
})

test_that("null and nuisance controls do not force outside states", {
    R <- diag(5)
    C <- matrix(c(0.4, -0.3, 0.2, -0.1, 0.05), 5, 1)
    null_fit <- susieout_fit(rep(0, 5), 100000, R, C,
                         paste0("i", 1:5), L = 2,
                         count_prior = susieout_count_prior(-2 * log(0.98), 0.2))
    expect_lt(null_fit$any_causal_omitted, 0.05)
    expect_gt(prod(null_fit$component_null), 0.9)

    dense_fit <- susieout_fit(c(0.8, -0.7, 0.6, -0.5, 0.4), 100000,
                          R, C, paste0("i", 1:5), L = 2,
                          count_prior = susieout_count_prior(-2 * log(0.98), 0.2),
                          infinitesimal_variance = 1e-5)
    mismatch_fit <- susieout_fit(c(0.8, -0.7, 0.6, -0.5, 0.4), 100000,
                             R, C, paste0("i", 1:5), L = 2,
                             count_prior = susieout_count_prior(-2 * log(0.98), 0.2),
                             ld_mismatch_variance = 0.2)
    expect_lt(dense_fit$any_causal_omitted, 0.1)
    expect_lt(mismatch_fit$any_causal_omitted, 0.1)
    expect_equal(dense_fit$infinitesimal$covariance_multiplier, 0.2)
    expect_gt(dense_fit$infinitesimal$variance_fraction, 0)
    expect_gt(mismatch_fit$ld_mismatch$variance_fraction, 0)
})

test_that("infinitesimal covariance uses total local variance on the z scale", {
    p <- 5L
    n <- 200000
    R <- diag(p)
    C <- matrix(seq(0.1, 0.5, length.out = p), p, 1)
    fit <- susieout_fit(rep(0, p), n, R, C, paste0("i", seq_len(p)),
                    L = 1, count_prior = susieout_count_prior(-log(0.99), 0.2),
                    infinitesimal_variance = 0.01)
    expect_equal(fit$infinitesimal$local_h2, 0.01)
    expect_equal(fit$infinitesimal$covariance_multiplier, n * 0.01 / p)
    expect_false("posterior_probability" %in% names(fit$infinitesimal))
})

test_that("nuisance-grid averaging is normalized and singleton-exact", {
    locus <- susieout_simulate_locus(
        411, p = 8, q = 2, n = 100000,
        architecture = "null", infinitesimal_h2 = 0.008,
        oracle_wgs = TRUE)
    arguments <- list(
        z = locus$z, n = locus$n, R_analysis = locus$R_analysis,
        wgs_ld = locus$wgs_ld, examined_ids = locus$examined_ids,
        L = 1, count_prior = susieout_count_prior(-log(0.9), 0.2))
    direct <- do.call(susieout_fit, c(arguments, list(
        infinitesimal_variance = 0.008, ld_mismatch_variance = 0)))
    singleton <- do.call(susieout_fit_nuisance_grid, c(arguments, list(
        infinitesimal_grid = 0.008, ld_mismatch_grid = 0)))
    expect_equal(singleton$pip_in, direct$pip_in, tolerance = 1e-12)
    expect_equal(singleton$any_causal_omitted,
                 direct$any_causal_omitted, tolerance = 1e-12)
    averaged <- do.call(susieout_fit_nuisance_grid, c(arguments, list(
        infinitesimal_grid = c(0, 0.008),
        ld_mismatch_grid = c(0, 0.08))))
    expect_equal(sum(averaged$diagnostics$nuisance_grid$posterior_weight),
                 1, tolerance = 1e-12)
    expect_true(averaged$infinitesimal$posterior_probability_nonzero >= 0 &&
                averaged$infinitesimal$posterior_probability_nonzero <= 1)
    expect_true(all(averaged$pip_in >= 0 & averaged$pip_in <= 1))
})

test_that("factorized finite-panel backend matches the dense joint model", {
    set.seed(902)
    samples <- 14L
    p <- 24L
    q <- 6L
    examined <- matrix(rbinom(samples * p, 2, 0.35), samples, p)
    outside <- matrix(rbinom(samples * q, 2, 0.20), samples, q)
    examined_ids <- paste0("i", seq_len(p))
    outside_ids <- paste0("o", seq_len(q))
    R <- cor(examined)
    C <- cor(examined, outside)
    z <- rnorm(p)
    grid <- expand.grid(
        infinitesimal_h2 = c(0, 0.002),
        ld_mismatch_variance = c(0, 0.04),
        KEEP.OUT.ATTRS = FALSE)
    grid$prior_weight <- c(0.35, 0.15, 0.35, 0.15)
    nuisance <- structure(list(
        grid = grid,
        infinitesimal_marginal = c(0.70, 0.30),
        ld_mismatch_marginal = c(0.50, 0.50),
        frozen = TRUE), class = "susieout_nuisance_prior")
    dense <- susieout_fit_nuisance_grid(
        z, 100000, R,
        list(R_SS = R, cross_ld = C,
             examined_ids = examined_ids, excluded_ids = outside_ids),
        examined_ids, L = 2,
        count_prior = susieout_count_prior(-2 * log(0.9), 0.2),
        nuisance_prior = nuisance, transform_wgs = FALSE)
    factorized <- susieout_fit_nuisance_grid_factorized(
        z, 100000, examined, outside,
        examined_ids, outside_ids, L = 2,
        count_prior = susieout_count_prior(-2 * log(0.9), 0.2),
        nuisance_prior = nuisance)
    expect_identical(factorized$backend,
                     "openworld_nuisance_grid_factorized")
    expect_equal(factorized$pip_in, dense$pip_in, tolerance = 1e-9)
    expect_equal(factorized$pip_out, dense$pip_out, tolerance = 1e-9)
    expect_equal(factorized$any_causal_omitted,
                 dense$any_causal_omitted, tolerance = 1e-9)
    expect_equal(
        factorized$diagnostics$nuisance_grid$posterior_weight,
        dense$diagnostics$nuisance_grid$posterior_weight,
        tolerance = 1e-9)
})

test_that("sparse LD backend matches dense inference and recovers outside signal", {
    p <- 120L
    R_sparse <- Matrix::bandSparse(
        p, k = c(-2L, -1L, 0L, 1L, 2L),
        diagonals = list(rep(0.08, p - 2L), rep(0.20, p - 1L),
                         rep(1, p), rep(0.20, p - 1L),
                         rep(0.08, p - 2L)))
    ids <- paste0("i", seq_len(p))
    C <- matrix(0, p, 2L,
                dimnames = list(ids, c("outside_1", "outside_2")))
    C[55:65, 1L] <- c(0.03, 0.05, 0.08, 0.12, 0.20, 0.30,
                       0.20, 0.12, 0.08, 0.05, 0.03)
    C[85:91, 2L] <- c(0.02, 0.06, 0.12, 0.24, 0.12, 0.06, 0.02)
    z <- sqrt(1e6) * C[, 1L] * 0.05
    arguments <- list(
        z = z, n = 1e6, wgs_ld = C, examined_ids = ids,
        effect_prior = list(weights = 1, means = 0, variances = 0.0025),
        L = 1, count_prior = susieout_count_prior(-log(0.5), 0.5))
    sparse_fit <- do.call(susieout_fit, c(arguments, list(R_analysis = R_sparse)))
    dense_fit <- do.call(susieout_fit, c(arguments, list(R_analysis = as.matrix(R_sparse))))
    expect_identical(sparse_fit$backend, "openworld_ibss_sparse")
    expect_true(sparse_fit$diagnostics$converged)
    expect_gt(unname(sparse_fit$pip_out["outside_1"]), 0.9)
    expect_equal(sparse_fit$pip_in, dense_fit$pip_in, tolerance = 1e-7)
    expect_equal(sparse_fit$pip_out, dense_fit$pip_out, tolerance = 1e-7)
    expect_lt(sparse_fit$diagnostics$maximum_responsibility_error, 1e-10)
    expect_error(do.call(susieout_fit, c(arguments, list(
        R_analysis = R_sparse, ld_mismatch_variance = 0.2))),
        "ld_mismatch_variance <= 1e-4", fixed = TRUE)
})
