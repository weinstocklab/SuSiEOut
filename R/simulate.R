#' Apply realistic candidate-set masks
#'
#' @param variants Data frame containing `maf`, `variant_type`, `info`,
#'   `reference_present`, and `panel_polymorphic`.
#' @export
susieout_candidate_mask <- function(variants, maf_min = 0.01,
                                snp_only = TRUE, info_min = 0.8) {
    required <- c("maf", "variant_type", "info", "reference_present",
                  "panel_polymorphic")
    if (!is.data.frame(variants) || !all(required %in% names(variants))) {
        stop("variants must contain maf, variant_type, info, reference_present, and panel_polymorphic",
             call. = FALSE)
    }
    if (anyNA(variants[, required])) stop("candidate-mask fields cannot be missing",
                                         call. = FALSE)
    variants$maf >= maf_min & variants$info >= info_min &
        (!snp_only | toupper(variants$variant_type) == "SNP") &
        as.logical(variants$reference_present) &
        as.logical(variants$panel_polymorphic)
}

.susieout_ar1 <- function(size, rho) {
    index <- seq_len(size)
    rho^abs(outer(index, index, `-`))
}

#' Simulate a deterministic SuSiE-OUT locus
#'
#' This compact truth-known generator supports smoke tests and workflow banks.
#' It represents recent growth through a rare-allele-enriched frequency draw,
#' prospective causal assignment, derived-allele orientation, finite WGS
#' panels, ancestry mismatch, and sparse/null/infinitesimal/LD-error regimes.
#'
#' @export
susieout_simulate_locus <- function(seed, p = 40L, q = 10L, n = 100000,
                                wgs_n = 1000L,
                                architecture = c("contained", "partially_omitted",
                                                 "entirely_omitted", "null",
                                                 "infinitesimal", "ld_mismatch"),
                                evolution = c("neutral", "stabilizing", "directional"),
                                ancestry = c("matched", "mismatched"),
                                causal_count = NULL,
                                omitted_count = NULL,
                                oracle_wgs = FALSE,
                                infinitesimal_h2 = NULL,
                                ld_mismatch_variance = NULL) {
    architecture <- match.arg(architecture)
    evolution <- match.arg(evolution)
    ancestry <- match.arg(ancestry)
    if (p < 2L || q < 1L || n <= 1 || wgs_n < 10L) {
        stop("p, q, n, and wgs_n are outside supported ranges", call. = FALSE)
    }
    set.seed(seed)
    total <- p + q
    population_R <- .susieout_ar1(total, rho = 0.88)
    analysis_R <- population_R[seq_len(p), seq_len(p), drop = FALSE]
    if (ancestry == "mismatched") {
        wgs_population_R <- .susieout_ar1(total, rho = 0.72)
    } else {
        wgs_population_R <- population_R
    }
    if (isTRUE(oracle_wgs)) {
        wgs_R <- wgs_population_R
    } else {
        eig <- .susieout_psd_eigen(wgs_population_R, 1e-10)
        draws <- matrix(stats::rnorm(wgs_n * total), wgs_n, total) %*%
            (sweep(t(eig$vectors), 1L, sqrt(eig$values), `*`))
        wgs_R <- stats::cor(draws)
    }
    maf <- pmin(0.5, pmax(1 / (2 * wgs_n),
                         exp(stats::runif(total, log(0.0005), log(0.5)))))
    variant_type <- ifelse(stats::runif(total) < 0.12, "INDEL", "SNP")
    variants <- data.frame(
        id = paste0("v", seq_len(total)), maf = maf,
        variant_type = variant_type,
        info = pmin(1, pmax(0, 0.75 + 0.35 * stats::runif(total))),
        reference_present = stats::runif(total) > 0.04,
        panel_polymorphic = 2 * wgs_n * maf > 1,
        derived_allele = sample(c("A", "C", "G", "T"), total, replace = TRUE),
        stringsAsFactors = FALSE)
    examined <- seq_len(p)
    excluded <- p + seq_len(q)
    variants$maf[examined] <- pmax(variants$maf[examined], 0.01)
    variants$variant_type[examined] <- "SNP"
    variants$info[examined] <- pmax(variants$info[examined], 0.8)
    variants$reference_present[examined] <- TRUE
    variants$panel_polymorphic[examined] <- TRUE
    exclusion_reason <- rep("examined", total)
    mechanisms <- c("maf_below_1pct", "indel", "info_dropout",
                    "reference_absence", "finite_panel_monomorphism")
    for (index in seq_along(excluded)) {
        variant <- excluded[index]
        mechanism <- mechanisms[(index - 1L) %% length(mechanisms) + 1L]
        exclusion_reason[variant] <- mechanism
        if (mechanism == "maf_below_1pct") variants$maf[variant] <- min(0.009, variants$maf[variant])
        if (mechanism == "indel") variants$variant_type[variant] <- "INDEL"
        if (mechanism == "info_dropout") variants$info[variant] <- min(0.79, variants$info[variant])
        if (mechanism == "reference_absence") variants$reference_present[variant] <- FALSE
        if (mechanism == "finite_panel_monomorphism") variants$panel_polymorphic[variant] <- FALSE
    }
    variants$exclusion_reason <- exclusion_reason
    if (is.null(causal_count)) {
        causal_count <- if (architecture %in%
            c("null", "infinitesimal", "ld_mismatch")) 0L else 1L
    }
    causal_count <- as.integer(causal_count)
    if (causal_count < 0L || causal_count > 6L) {
        stop("causal_count must be between zero and six", call. = FALSE)
    }
    if (is.null(omitted_count)) {
        omitted_count <- switch(
            architecture,
            contained = 0L,
            entirely_omitted = causal_count,
            partially_omitted = if (causal_count > 0L)
                max(causal_count - 1L, 1L) else 0L,
            0L)
    }
    omitted_count <- as.integer(omitted_count)
    if (omitted_count < 0L || omitted_count > causal_count ||
        omitted_count > q || causal_count - omitted_count > p) {
        stop("omitted_count is incompatible with causal_count, p, or q",
             call. = FALSE)
    }
    causal <- integer()
    if (causal_count > 0L) {
        inside_count <- causal_count - omitted_count
        causal <- c(
            if (inside_count > 0L) sample(examined, inside_count) else integer(),
            if (omitted_count > 0L) sample(excluded, omitted_count) else integer())
    }
    beta <- numeric(length(causal))
    if (length(causal)) {
        beta <- switch(
            evolution,
            neutral = stats::rnorm(length(causal), 0, 0.018),
            stabilizing = stats::rnorm(length(causal), 0,
                0.012 * (pmax(variants$maf[causal], 0.001) / 0.01)^(-0.18)),
            directional = stats::rnorm(length(causal), 0.018, 0.010))
        if (evolution == "directional") {
            derived_sign <- ifelse(variants$derived_allele[causal] %in% c("A", "G"), 1, -1)
            beta <- abs(beta) * derived_sign
        }
    }
    mean_z_full <- numeric(total)
    if (length(causal)) {
        mean_z_full <- sqrt(n) * as.numeric(
            population_R[, causal, drop = FALSE] %*% beta)
    }
    if (is.null(infinitesimal_h2)) {
        # The historical dense simulation used per-variant variance 2e-4 at
        # p = 40, corresponding to total local variance h2 = 0.008.
        infinitesimal_h2 <- if (architecture == "infinitesimal") 0.008 else 0
    }
    if (is.null(ld_mismatch_variance)) {
        ld_mismatch_variance <- if (architecture == "ld_mismatch") 0.08 else 0
    }
    if (!is.numeric(infinitesimal_h2) || length(infinitesimal_h2) != 1L ||
        !is.finite(infinitesimal_h2) || infinitesimal_h2 < 0) {
        stop("infinitesimal_h2 must be a nonnegative finite scalar", call. = FALSE)
    }
    if (!is.numeric(ld_mismatch_variance) ||
        length(ld_mismatch_variance) != 1L ||
        !is.finite(ld_mismatch_variance) || ld_mismatch_variance < 0) {
        stop("ld_mismatch_variance must be a nonnegative finite scalar",
             call. = FALSE)
    }
    infinitesimal_per_variant_variance <- infinitesimal_h2 / p
    dense_beta <- if (infinitesimal_h2 > 0)
        stats::rnorm(p, 0, sqrt(infinitesimal_per_variant_variance)) else numeric(p)
    if (length(dense_beta)) {
        mean_z_full[examined] <- mean_z_full[examined] +
            sqrt(n) * as.numeric(analysis_R %*% dense_beta)
    }
    noise_full <- .susieout_mvn_draw(
        population_R + diag(ld_mismatch_variance, total))
    z_full <- mean_z_full + noise_full
    z <- z_full[examined]
    list(
        z = z, n = n, R_analysis = analysis_R,
        wgs_ld = list(
            R_SS = wgs_R[examined, examined, drop = FALSE],
            cross_ld = wgs_R[examined, excluded, drop = FALSE],
            examined_ids = variants$id[examined],
            excluded_ids = variants$id[excluded]),
        examined_ids = variants$id[examined], variants = variants,
        oracle = list(z = z_full, R = population_R, ids = variants$id),
        split = "training",
        truth = list(
            architecture = architecture, evolution = evolution,
            ancestry = ancestry, causal = variants$id[causal],
            causal_index = causal, causal_beta = beta,
            any_omitted = any(causal %in% excluded),
            omitted_count = sum(causal %in% excluded),
            is_null = !length(causal) && infinitesimal_h2 == 0 &&
                ld_mismatch_variance == 0,
            sparse_architecture = architecture,
            infinitesimal_variance = infinitesimal_h2,
            infinitesimal_h2 = infinitesimal_h2,
            infinitesimal_per_variant_variance =
                infinitesimal_per_variant_variance,
            infinitesimal_covariance_multiplier =
                .susieout_infinitesimal_multiplier(infinitesimal_h2, n, p),
            ld_mismatch_variance = ld_mismatch_variance,
            max_proxy_r2 = if (any(causal %in% excluded))
                max(population_R[examined, causal[causal %in% excluded]]^2) else NA_real_,
            effective_rank = .susieout_participation_rank(analysis_R)),
        seed = seed)
}
