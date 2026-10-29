test_that("count-prior locus probabilities are invariant to L", {
    prior <- susieout_count_prior(1.7, 0.3)
    summaries <- lapply(c(3L, 10L, 20L), function(L)
        SuSiEOut:::.susieout_count_prior_masses(prior, L))
    expect_equal(vapply(summaries, `[[`, numeric(1), "probability_any_sparse"),
                 rep(1 - exp(-1.7), 3), tolerance = 1e-14)
    expect_true(all(vapply(summaries, function(x)
        abs(x$null + x$inside + x$outside - 1), numeric(1)) < 1e-14))
})

test_that("OUT candidate count does not change total prior mass", {
    R <- diag(4)
    z <- rep(0, 4)
    prior <- susieout_count_prior(1, 0.4)
    one <- susieout_fit(z, 10000, R, matrix(0.1, 4, 1),
                        paste0("i", 1:4), L = 2, count_prior = prior)
    many <- susieout_fit(z, 10000, R, matrix(0.1, 4, 8),
                         paste0("i", 1:4), L = 2, count_prior = prior)
    expect_equal(sum(one$prior$state[5]),
                 sum(many$prior$state[5:12]), tolerance = 1e-14)
    expect_equal(one$prior$probability_any_sparse,
                 many$prior$probability_any_sparse, tolerance = 1e-14)
})

test_that("L=1 variational evidence equals exact marginal likelihood", {
    R <- diag(3)
    C <- matrix(c(0.4, -0.3, 0.2), 3, 1)
    z <- c(2, -1, 0.5)
    n <- 1000
    variance <- 0.002
    fit <- susieout_fit(
        z, n, R, C, paste0("i", 1:3),
        effect_prior = list(weights = 1, means = 0, variances = variance),
        L = 1, count_prior = susieout_count_prior(0.8, 0.3))
    fingerprints <- cbind(R, C)
    covariance <- R + diag(1e-6, 3)
    log_density <- function(V) {
        factor <- chol(V)
        -0.5 * (3 * log(2 * pi) + 2 * sum(log(diag(factor))) +
            sum(backsolve(factor, z, transpose = TRUE)^2))
    }
    values <- c(
        log(fit$prior$state[5]) + log_density(covariance),
        vapply(seq_len(4), function(index)
            log(fit$prior$state[index]) + log_density(
                covariance + n * variance *
                    tcrossprod(fingerprints[, index])), numeric(1)))
    exact <- max(values) + log(sum(exp(values - max(values))))
    expect_equal(fit$diagnostics$variational_evidence, exact,
                 tolerance = 1e-10)
})

test_that("singleton prior grid exactly reproduces its conditional fit", {
    locus <- susieout_simulate_locus(
        717, p = 8, q = 2, architecture = "partially_omitted",
        causal_count = 2, omitted_count = 1, oracle_wgs = TRUE)
    count <- susieout_count_prior(1, 0.3)
    direct <- susieout_fit(
        locus$z, locus$n, locus$R_analysis, locus$wgs_ld,
        locus$examined_ids, L = 2, count_prior = count)
    grid <- susieout_prior_grid(1, 0.3, weights = 1)
    averaged <- susieout_fit_grid(
        locus$z, locus$n, locus$R_analysis, locus$wgs_ld,
        locus$examined_ids, L = 2, prior_grid = grid, starts = "zero")
    expect_equal(averaged$pip_in, direct$pip_in, tolerance = 1e-12)
    expect_equal(averaged$pip_out, direct$pip_out, tolerance = 1e-12)
    expect_equal(averaged$diagnostics$model_grid$posterior_weight, 1)
})

test_that("six sparse signals and explicit omitted counts are supported", {
    locus <- susieout_simulate_locus(
        812, p = 12, q = 8, architecture = "partially_omitted",
        causal_count = 6, omitted_count = 3, oracle_wgs = TRUE)
    expect_equal(length(locus$truth$causal), 6)
    expect_equal(locus$truth$omitted_count, 3)
})

test_that("legacy component-prior arguments are rejected", {
    expect_error(susieout_fit(
        rep(0, 3), 1000, diag(3), matrix(0.1, 3, 1),
        paste0("i", 1:3), outside_prior = 0.1),
        "legacy or extra fit arguments are not supported")
})

test_that("masked-WGS calibration learns normalized global weights", {
    grid <- susieout_prior_grid(c(0.5, 1), c(0, 0.5),
                                dirichlet_strength = 0.1)
    states <- nrow(grid$grid)
    bank <- lapply(seq_len(4), function(index) list(
        evidence = rep(0, states),
        pip_in = cbind(
            c(0.8, 0.1), c(0.7, 0.2),
            c(0.2, 0.6), c(0.1, 0.7)),
        target_pip = c(0.15, 0.65)))
    learned <- susieout_calibrate_prior_grid(bank, grid,
                                             dirichlet_strength = 0.1)
    expect_s3_class(learned, "susieout_prior_grid")
    expect_equal(sum(learned$grid$prior_weight), 1, tolerance = 1e-12)
    expect_true(learned$calibration$optim_converged)
    expect_gt(sum(learned$grid$prior_weight[3:4]),
              sum(grid$grid$prior_weight[3:4]))
})
