#' Harmonize WGS cross-LD rows to analysis alleles
#'
#' @param cross_ld Tested-by-omitted cross-LD matrix.
#' @param wgs_examined Data frame with `id`, `effect_allele`, and `other_allele`.
#' @param analysis_examined Data frame with the same columns in analysis order.
#' @param allow_strand Whether unambiguous strand complements may be aligned.
#' @return The row-aligned, sign-corrected cross-LD matrix.
#' @export
susieout_harmonize_cross_ld <- function(
    cross_ld,
    wgs_examined,
    analysis_examined,
    allow_strand = TRUE
) {
    required <- c("id", "effect_allele", "other_allele")
    if (
        !is.data.frame(wgs_examined) ||
            !all(required %in% names(wgs_examined)) ||
            !is.data.frame(analysis_examined) ||
            !all(required %in% names(analysis_examined))
    ) {
        stop(
            "allele tables must contain id, effect_allele, and other_allele",
            call. = FALSE
        )
    }
    if (anyDuplicated(wgs_examined$id) || anyDuplicated(analysis_examined$id)) {
        stop("allele-table IDs must be unique", call. = FALSE)
    }
    index <- match(analysis_examined$id, wgs_examined$id)
    if (anyNA(index)) {
        stop("reference LD is missing tested variant IDs", call. = FALSE)
    }
    cross_ld <- .susieout_validate_matrix(
        cross_ld,
        "cross_ld",
        nrow = nrow(wgs_examined)
    )
    cross_ld <- cross_ld[index, , drop = FALSE]
    wgs <- wgs_examined[index, , drop = FALSE]
    analysis_effect <- toupper(as.character(analysis_examined$effect_allele))
    analysis_other <- toupper(as.character(analysis_examined$other_allele))
    wgs_effect <- toupper(as.character(wgs$effect_allele))
    wgs_other <- toupper(as.character(wgs$other_allele))
    same <- analysis_effect == wgs_effect & analysis_other == wgs_other
    swapped <- analysis_effect == wgs_other & analysis_other == wgs_effect
    if (allow_strand) {
        strand_same <- analysis_effect == .susieout_complement(wgs_effect) &
            analysis_other == .susieout_complement(wgs_other)
        strand_swapped <- analysis_effect == .susieout_complement(wgs_other) &
            analysis_other == .susieout_complement(wgs_effect)
        same <- same | strand_same
        swapped <- swapped | strand_swapped
    }
    unresolved <- !(same | swapped)
    if (any(unresolved)) {
        stop(
            "unresolved allele mismatch for: ",
            paste(analysis_examined$id[unresolved], collapse = ", "),
            call. = FALSE
        )
    }
    cross_ld[swapped, ] <- -cross_ld[swapped, , drop = FALSE]
    rownames(cross_ld) <- as.character(analysis_examined$id)
    cross_ld
}

#' Ridge-stabilized whitening and recoloring of WGS fingerprints
#'
#' Computes `R_analysis^(1/2) R_wgs_SS^(-1/2) cross_ld` with eigenvalue
#' regularization, including rank-deficient input matrices. Use it when LD
#' among tested variants in the sequencing panel differs from the analysis
#' LD, so that the fingerprints are expressed on the analysis LD scale.
#'
#' @param R_analysis LD matrix among tested variants used in the analysis.
#' @param R_wgs_SS LD matrix among the same tested variants in the
#'   sequencing reference panel.
#' @param cross_ld Reference-panel LD between tested variants (rows) and
#'   omitted variants (columns).
#' @param ridge Positive eigenvalue floor used for regularization.
#' @return A matrix of transported fingerprints with the same dimensions as
#'   `cross_ld`.
#' @export
susieout_transform_fingerprints <- function(
    R_analysis,
    R_wgs_SS,
    cross_ld,
    ridge = 1e-6
) {
    if (
        !is.numeric(ridge) ||
            length(ridge) != 1L ||
            !is.finite(ridge) ||
            ridge <= 0
    ) {
        stop("ridge must be positive", call. = FALSE)
    }
    if (
        inherits(R_analysis, "sparseMatrix") ||
            inherits(R_wgs_SS, "sparseMatrix")
    ) {
        stop(
            paste(
                "whitening/recoloring currently requires dense LD matrices;",
                "supply pre-transported sparse-locus fingerprints without R_SS"
            ),
            call. = FALSE
        )
    }
    p <- nrow(R_analysis)
    R_analysis <- .susieout_validate_matrix(
        R_analysis,
        "R_analysis",
        p,
        p,
        symmetric = TRUE
    )
    R_wgs_SS <- .susieout_validate_matrix(
        R_wgs_SS,
        "R_wgs_SS",
        p,
        p,
        symmetric = TRUE
    )
    cross_ld <- .susieout_validate_matrix(cross_ld, "cross_ld", nrow = p)
    analysis_half <- .susieout_matrix_power(R_analysis, 0.5, ridge)
    wgs_inverse_half <- .susieout_matrix_power(R_wgs_SS, -0.5, ridge)
    transformed <- analysis_half %*% wgs_inverse_half %*% cross_ld
    dimnames(transformed) <- dimnames(cross_ld)
    transformed
}

.susieout_parse_wgs_ld <- function(
    wgs_ld,
    examined_ids,
    R_analysis,
    transform_wgs,
    ridge
) {
    p <- length(examined_ids)
    if (is.matrix(wgs_ld)) {
        cross <- .susieout_validate_matrix(wgs_ld, "wgs_ld", nrow = p)
        excluded_ids <- .susieout_validate_ids(
            colnames(cross),
            ncol(cross),
            "excluded_"
        )
        return(list(
            raw = cross,
            transformed = cross,
            excluded_ids = excluded_ids,
            R_wgs_SS = NULL,
            transform_applied = FALSE,
            bootstrap = NULL
        ))
    }
    if (!is.list(wgs_ld)) {
        stop("wgs_ld must be a cross-LD matrix or a WGS LD list", call. = FALSE)
    }
    cross <- wgs_ld$cross_ld %||% wgs_ld$C
    if (is.null(cross)) {
        stop("wgs_ld must contain cross_ld", call. = FALSE)
    }
    if (!is.null(wgs_ld$wgs_examined) || !is.null(wgs_ld$analysis_examined)) {
        if (is.null(wgs_ld$wgs_examined) || is.null(wgs_ld$analysis_examined)) {
            stop(
                "both WGS and analysis allele tables are required",
                call. = FALSE
            )
        }
        cross <- susieout_harmonize_cross_ld(
            cross,
            wgs_ld$wgs_examined,
            wgs_ld$analysis_examined
        )
    } else if (!is.null(wgs_ld$examined_ids)) {
        index <- match(examined_ids, as.character(wgs_ld$examined_ids))
        if (anyNA(index)) {
            stop(
                "reference LD does not contain every tested variant ID",
                call. = FALSE
            )
        }
        cross <- cross[index, , drop = FALSE]
    }
    cross <- .susieout_validate_matrix(cross, "wgs_ld$cross_ld", nrow = p)
    excluded_ids <- .susieout_validate_ids(
        wgs_ld$excluded_ids %||% colnames(cross),
        ncol(cross),
        "excluded_"
    )
    if (length(intersect(examined_ids, excluded_ids))) {
        stop("tested and omitted variant IDs must be disjoint", call. = FALSE)
    }
    R_wgs_SS <- wgs_ld$R_SS
    transformed <- cross
    applied <- FALSE
    if (isTRUE(transform_wgs) && !is.null(R_wgs_SS)) {
        if (!is.null(wgs_ld$examined_ids)) {
            R_wgs_SS <- R_wgs_SS[index, index, drop = FALSE]
        }
        if (max(abs(R_analysis - R_wgs_SS)) <= 1e-10) {
            transformed <- cross
        } else {
            transformed <- susieout_transform_fingerprints(
                R_analysis,
                R_wgs_SS,
                cross,
                ridge
            )
            applied <- TRUE
        }
    }
    colnames(cross) <- colnames(transformed) <- excluded_ids
    rownames(cross) <- rownames(transformed) <- examined_ids
    list(
        raw = cross,
        transformed = transformed,
        excluded_ids = excluded_ids,
        R_wgs_SS = R_wgs_SS,
        transform_applied = applied,
        bootstrap = wgs_ld$bootstrap %||% NULL
    )
}

#' Bootstrap WGS cross-LD fingerprints
#'
#' @param genotypes Individual-by-variant dosage matrix.
#' @param examined Integer or character indices of tested variants.
#' @param excluded Integer or character indices of omitted variants.
#' @param B Number of bootstrap panels.
#' @param seed Reproducibility seed.
#' @return A list containing bootstrapped cross-LD matrices and uncertainty.
#' @export
susieout_bootstrap_fingerprints <- function(
    genotypes,
    examined,
    excluded,
    B = 100L,
    seed = 1L
) {
    genotypes <- .susieout_validate_matrix(genotypes, "genotypes")
    if (B < 1L) {
        stop("B must be positive", call. = FALSE)
    }
    examined <- if (is.character(examined)) {
        match(examined, colnames(genotypes))
    } else {
        examined
    }
    excluded <- if (is.character(excluded)) {
        match(excluded, colnames(genotypes))
    } else {
        excluded
    }
    if (
        anyNA(c(examined, excluded)) ||
            any(
                c(examined, excluded) < 1L |
                    c(examined, excluded) > ncol(genotypes)
            )
    ) {
        stop(
            "tested and omitted variants must index genotype columns",
            call. = FALSE
        )
    }
    set.seed(seed)
    matrices <- vector("list", B)
    for (b in seq_len(B)) {
        rows <- sample.int(nrow(genotypes), nrow(genotypes), replace = TRUE)
        correlation <- stats::cor(genotypes[
            rows,
            c(examined, excluded),
            drop = FALSE
        ])
        if (anyNA(correlation)) {
            stop("bootstrap produced a monomorphic variant", call. = FALSE)
        }
        matrices[[b]] <- correlation[
            seq_along(examined),
            length(examined) + seq_along(excluded),
            drop = FALSE
        ]
    }
    array <- simplify2array(matrices)
    list(
        cross_ld = matrices,
        mean = apply(array, c(1L, 2L), mean),
        sd = apply(array, c(1L, 2L), stats::sd),
        B = as.integer(B),
        seed = as.integer(seed)
    )
}
