test_that("inputs reject inconsistent and missing data", {
    R <- diag(3)
    C <- matrix(c(0.2, 0.1, 0.3), 3, 1)
    expect_error(susieout_fit(c(1, 2), 1000, R, C), "length|rows")
    expect_error(susieout_fit(c(1, NA, 2), 1000, R, C), "missing")
    expect_error(susieout_fit(c(1, 2, 3), 1000, R,
                          matrix(c(0.2, NA, 0.1), 3, 1)), "missing")
    bad_R <- R
    bad_R[1, 2] <- 0.4
    expect_error(susieout_fit(c(1, 2, 3), 1000, bad_R, C), "symmetric")
    indefinite <- matrix(c(1, 2, 2, 1), 2)
    expect_error(susieout_fit(c(1, 2), 1000, indefinite,
                          matrix(c(0.2, 0.1), 2, 1)),
                 "positive semidefinite")
})

test_that("rank-deficient PSD analysis LD is stabilized", {
    R <- matrix(1, 3, 3)
    C <- matrix(c(0.8, 0.8, 0.8), 3, 1)
    fit <- susieout_fit(c(2, 2, 2), 5000, R, C,
                    examined_ids = paste0("i", 1:3), L = 1)
    expect_s3_class(fit, "susieout_fit")
    expect_true(is.finite(fit$diagnostics$effective_rank))
    expect_lt(fit$diagnostics$maximum_responsibility_error, 1e-10)
})

test_that("allele harmonization reorders and flips LD", {
    cross <- matrix(c(0.2, -0.4), 2, 1,
                    dimnames = list(c("b", "a"), "x"))
    wgs <- data.frame(
        id = c("b", "a"), effect_allele = c("A", "C"),
        other_allele = c("G", "T"))
    analysis <- data.frame(
        id = c("a", "b"), effect_allele = c("T", "T"),
        other_allele = c("C", "C"))
    aligned <- susieout_harmonize_cross_ld(cross, wgs, analysis)
    expect_equal(as.numeric(aligned), c(0.4, 0.2))
    broken <- analysis
    broken$effect_allele[1] <- "A"
    expect_error(susieout_harmonize_cross_ld(cross, wgs, broken),
                 "unresolved allele")
})

test_that("whitening and recoloring has the expected identity case", {
    R <- matrix(c(1, 0.3, 0.3, 1), 2)
    C <- matrix(c(0.5, 0.2), 2, 1)
    transformed <- susieout_transform_fingerprints(R, R, C)
    expect_equal(transformed, C, tolerance = 1e-6)

    rank_one <- matrix(1, 2, 2)
    stable <- susieout_transform_fingerprints(R, rank_one, C)
    expect_true(all(is.finite(stable)))
})

test_that("WGS bootstrap is seeded and reproducible", {
    set.seed(9)
    G <- matrix(rnorm(600), 100, 6)
    first <- susieout_bootstrap_fingerprints(G, 1:4, 5:6, B = 8, seed = 17)
    second <- susieout_bootstrap_fingerprints(G, 1:4, 5:6, B = 8, seed = 17)
    expect_equal(first$mean, second$mean)
    expect_equal(first$sd, second$sd)
    expect_equal(first$B, 8L)
})

