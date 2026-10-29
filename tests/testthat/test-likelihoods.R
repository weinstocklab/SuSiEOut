test_that("one-fingerprint integrated likelihood matches a direct normal density", {
    residual <- 1.7
    fingerprint <- matrix(0.65, 1, 1)
    covariance <- matrix(1.2, 1, 1)
    chol_factor <- chol(covariance)
    solved <- SuSiEOut:::.susieout_solve_chol(chol_factor, fingerprint)
    quadratic <- colSums(fingerprint * solved)
    mixture <- list(weights = 1, means = 0, variances = 0.04)
    prior <- c(0.35, 0.65)
    update <- SuSiEOut:::.susieout_rank_one_update(
        residual, fingerprint, solved, quadratic, sqrt(200), mixture,
        prior, chol_factor)
    alternative_sd <- sqrt(covariance[1, 1] + 200 * 0.04 * fingerprint[1, 1]^2)
    bf <- dnorm(residual, 0, alternative_sd) /
        dnorm(residual, 0, sqrt(covariance[1, 1]))
    direct <- prior[1] * bf / (prior[1] * bf + prior[2])
    expect_equal(update$alpha, direct, tolerance = 1e-10)
    expect_equal(update$alpha + update$null, 1, tolerance = 1e-12)
})

test_that("one-inside one-outside responsibilities sum to one", {
    fingerprints <- matrix(c(1, 0.35), 1, 2)
    covariance <- matrix(1, 1, 1)
    chol_factor <- chol(covariance)
    solved <- SuSiEOut:::.susieout_solve_chol(chol_factor, fingerprints)
    mixture <- list(weights = c(0.7, 0.3), means = c(0, 0),
                    variances = c(0.01, 0.1))
    update <- SuSiEOut:::.susieout_rank_one_update(
        2.2, fingerprints, solved, colSums(fingerprints * solved),
        sqrt(1000), mixture, c(0.08, 0.02, 0.90), chol_factor)
    expect_equal(sum(update$alpha) + update$null, 1, tolerance = 1e-12)
    expect_true(all(update$alpha >= 0 & update$alpha <= 1))
})

test_that("fit responsibilities and reported probabilities are coherent", {
    locus <- susieout_simulate_locus(41, p = 12, q = 3, n = 500000,
                                 architecture = "entirely_omitted",
                                 oracle_wgs = TRUE)
    fit <- susieout_fit(locus$z, locus$n, locus$R_analysis, locus$wgs_ld,
                    locus$examined_ids, L = 2,
                    count_prior = susieout_count_prior(1, 0.3))
    expect_equal(rowSums(cbind(fit$alpha_in, fit$alpha_out)) +
                     fit$component_null,
                 rep(1, 2), tolerance = 1e-10)
    expect_equal(fit$any_causal_omitted,
                 1 - prod(1 - fit$component_outside))
    expect_equal(fit$candidate_set_complete, 1 - fit$any_causal_omitted)
    expect_true(all(fit$pip_in >= 0 & fit$pip_in <= 1))
    expect_true(all(fit$pip_out >= 0 & fit$pip_out <= 1))
})
