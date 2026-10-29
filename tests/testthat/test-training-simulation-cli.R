test_that("simulation is deterministic and exposes SuSiE-OUT truth", {
    first <- susieout_simulate_locus(2026, p = 10, q = 3,
                                 architecture = "partially_omitted",
                                 causal_count = 2)
    second <- susieout_simulate_locus(2026, p = 10, q = 3,
                                  architecture = "partially_omitted",
                                  causal_count = 2)
    expect_equal(first$z, second$z)
    expect_equal(first$wgs_ld$cross_ld, second$wgs_ld$cross_ld)
    expect_true(first$truth$any_omitted)
    expect_equal(length(first$truth$causal), 2)
})

test_that("dense and sparse effects can coexist with coherent scaling", {
    locus <- susieout_simulate_locus(
        2027, p = 20, q = 4, n = 500000,
        architecture = "entirely_omitted", causal_count = 2,
        infinitesimal_h2 = 0.008, ld_mismatch_variance = 0.04,
        oracle_wgs = TRUE)
    expect_true(locus$truth$any_omitted)
    expect_equal(locus$truth$infinitesimal_h2, 0.008)
    expect_equal(locus$truth$infinitesimal_per_variant_variance, 0.008 / 20)
    expect_equal(locus$truth$infinitesimal_covariance_multiplier,
                 500000 * 0.008 / 20)
    expect_equal(locus$truth$ld_mismatch_variance, 0.04)
    expect_false(locus$truth$is_null)
})

test_that("candidate masks combine all exclusion mechanisms", {
    variants <- data.frame(
        maf = c(0.02, 0.005, 0.02, 0.02, 0.02),
        variant_type = c("SNP", "SNP", "INDEL", "SNP", "SNP"),
        info = c(0.9, 0.9, 0.9, 0.7, 0.9),
        reference_present = c(TRUE, TRUE, TRUE, TRUE, FALSE),
        panel_polymorphic = c(TRUE, TRUE, TRUE, TRUE, TRUE))
    expect_equal(susieout_candidate_mask(variants),
                 c(TRUE, FALSE, FALSE, FALSE, FALSE))
})

test_that("training freezes priors and rejects locked test rows", {
    bank <- data.frame(
        omitted = c(TRUE, FALSE, TRUE, FALSE),
        omitted_count = c(1, 0, 1, 0),
        causal_count = c(1, 1, 1, 0),
        is_null = c(FALSE, FALSE, FALSE, TRUE),
        split = c("training", "training", "validation", "validation"),
        infinitesimal_variance = c(0, 0, 0.1, 0),
        ld_mismatch_variance = c(0, 0.05, 0, 0),
        causal_beta = c(0.02, 0.03, 0.04, NA))
    prior <- susieout_train_prior(bank)
    expect_s3_class(prior, "susieout_prior")
    expect_true(prior$frozen)
    expect_equal(prior$count_prior$lambda_total, 1)
    expect_equal(prior$count_prior$out_fraction, 0.5)
    expect_equal(prior$training_loci, 2L)
    expect_identical(prior$training_splits, "training")
    expect_equal(prior$infinitesimal_variance, 0)
    bank$split[1] <- "locked_test"
    expect_error(susieout_train_prior(bank), "held-out/locked")
})

test_that("full-locus training performs responsibility and nuisance M-steps", {
    bank <- lapply(seq_len(6), function(index) {
        locus <- susieout_simulate_locus(
            900 + index, p = 8, q = 2, n = 100000,
            architecture = if (index %% 2) "contained" else "entirely_omitted",
            oracle_wgs = TRUE)
        locus$split <- if (index <= 4) "training" else "validation"
        locus$truth$split <- locus$split
        locus
    })
    prior <- susieout_train_prior(bank, L = 2, em_iterations = 1,
                              em_max_loci = 4, em_update = "nuisance_only",
                              seed = 11)
    expect_equal(prior$em_iterations, 1L)
    expect_equal(nrow(prior$em_trace), 1L)
    expect_identical(prior$em_update, "nuisance_only")
    expect_true(prior$em_trace$nuisance_objective_after <=
                prior$em_trace$nuisance_objective_before + 1e-8)
    expect_true(all(prior$effect_prior$weights >= 0))
    expect_equal(sum(prior$effect_prior$weights), 1, tolerance = 1e-10)
})

test_that("nuisance-grid weights are learned from training loci only", {
    bank <- lapply(seq_len(4), function(index) {
        locus <- susieout_simulate_locus(
            1200 + index, p = 6, q = 2, n = 100000,
            architecture = "null",
            infinitesimal_h2 = if (index <= 2) 0 else 0.008,
            oracle_wgs = TRUE)
        locus$split <- "training"
        locus$truth$split <- "training"
        locus
    })
    sparse_prior <- susieout_train_prior(bank, L = 1)
    nuisance_prior <- susieout_train_nuisance_grid(
        bank, sparse_prior, infinitesimal_grid = c(0, 0.008),
        ld_mismatch_grid = 0, L = 1, iterations = 3)
    expect_s3_class(nuisance_prior, "susieout_nuisance_prior")
    expect_equal(sum(nuisance_prior$grid$prior_weight), 1,
                 tolerance = 1e-12)
    expect_equal(nuisance_prior$training_loci, 4L)
    bank[[1L]]$split <- "validation"
    bank[[1L]]$truth$split <- "validation"
    expect_error(susieout_train_nuisance_grid(
        bank, sparse_prior, infinitesimal_grid = c(0, 0.008),
        ld_mismatch_grid = 0, L = 1, iterations = 1),
        "training loci only")
})

test_that("CLI writes a fit and candidate-support table", {
    locus <- susieout_simulate_locus(81, p = 8, q = 2, n = 200000,
                                 architecture = "entirely_omitted",
                                 oracle_wgs = TRUE)
    input <- tempfile(fileext = ".rds")
    output <- tempfile(fileext = ".rds")
    support <- tempfile(fileext = ".tsv")
    saveRDS(locus, input)
    expect_equal(susieout_cli(c(
        "--input", input, "--output", output, "--summary", support,
        "--L", "2", "--lambda-total", "1", "--out-fraction", "0.25")), 0L)
    expect_s3_class(readRDS(output), "susieout_fit")
    expect_true(file.exists(support))
    table <- read.delim(support)
    expect_equal(nrow(table), 10)
})
