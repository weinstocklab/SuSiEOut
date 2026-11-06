test_that("zero OUT rate supports candidate-complete native fitting", {
    R <- matrix(c(1, 0.2, 0.2, 1), 2)
    z <- c(4, 1)
    wrapped <- susieout_fit(
        z,
        10000,
        R,
        wgs_ld = NULL,
        examined_ids = c("a", "b"),
        L = 2,
        count_prior = susieout_count_prior(1, 0)
    )
    expect_true(all(wrapped$pip_in >= 0 & wrapped$pip_in <= 1))
    expect_equal(ncol(wrapped$alpha_out), 0)
    expect_equal(
        rowSums(wrapped$alpha_in) + wrapped$component_null,
        rep(1, 2),
        tolerance = 1e-10
    )
    expect_equal(wrapped$any_causal_omitted, 0)
    expect_equal(wrapped$candidate_set_complete, 1)
    expect_equal(nrow(susieout_candidate_support(wrapped)), 2)
})
