`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

.susieout_logsumexp <- function(x) {
    if (length(x) == 0L) {
        return(-Inf)
    }
    maximum <- max(x)
    if (!is.finite(maximum)) {
        return(maximum)
    }
    maximum + log(sum(exp(x - maximum)))
}

.susieout_validate_vector <- function(x, name, length = NULL) {
    if (!is.numeric(x) || is.matrix(x) || anyNA(x) || any(!is.finite(x))) {
        stop(
            name,
            " must be a finite numeric vector without missing values",
            call. = FALSE
        )
    }
    if (!is.null(length) && base::length(x) != length) {
        stop(name, " must have length ", length, call. = FALSE)
    }
    as.numeric(x)
}

.susieout_validate_matrix <- function(
    x,
    name,
    nrow = NULL,
    ncol = NULL,
    symmetric = FALSE
) {
    sparse <- inherits(x, "sparseMatrix")
    finite <- if (sparse) {
        is.numeric(x@x) && !anyNA(x@x) && all(is.finite(x@x))
    } else {
        is.matrix(x) && is.numeric(x) && !anyNA(x) && all(is.finite(x))
    }
    if (!finite) {
        stop(
            name,
            " must be a finite numeric matrix without missing values",
            call. = FALSE
        )
    }
    if (!is.null(nrow) && base::nrow(x) != nrow) {
        stop(name, " must have ", nrow, " rows", call. = FALSE)
    }
    if (!is.null(ncol) && base::ncol(x) != ncol) {
        stop(name, " must have ", ncol, " columns", call. = FALSE)
    }
    if (symmetric) {
        if (base::nrow(x) != base::ncol(x)) {
            stop(name, " must be square", call. = FALSE)
        }
        transpose <- if (sparse) Matrix::t(x) else t(x)
        asymmetry <- max(abs(x - transpose))
        if (asymmetry > 1e-6) {
            stop(
                name,
                " must be symmetric (maximum asymmetry ",
                signif(asymmetry, 3),
                ")",
                call. = FALSE
            )
        }
        x <- (x + transpose) / 2
        if (sparse) x <- Matrix::drop0(x)
    }
    if (!sparse) {
        storage.mode(x) <- "double"
    }
    x
}

.susieout_validate_ids <- function(ids, expected, name) {
    if (is.null(ids)) {
        ids <- paste0(name, seq_len(expected))
    }
    if (
        length(ids) != expected || anyNA(ids) || any(!nzchar(as.character(ids)))
    ) {
        stop(
            name,
            " must contain one non-missing ID per variant",
            call. = FALSE
        )
    }
    ids <- as.character(ids)
    if (anyDuplicated(ids)) {
        stop(name, " must be unique", call. = FALSE)
    }
    ids
}

.susieout_diag <- function(x) {
    if (inherits(x, "sparseMatrix")) Matrix::diag(x) else diag(x)
}

.susieout_psd_eigen <- function(x, ridge = 1e-8) {
    x <- (x + t(x)) / 2
    decomposition <- eigen(x, symmetric = TRUE)
    scale <- max(1, max(abs(decomposition$values)))
    if (min(decomposition$values) < -1e-6 * scale) {
        stop(
            "LD matrix is not positive semidefinite; minimum eigenvalue is ",
            signif(min(decomposition$values), 4),
            call. = FALSE
        )
    }
    values <- pmax(decomposition$values, ridge * scale)
    list(
        vectors = decomposition$vectors,
        values = values,
        raw_values = decomposition$values
    )
}

.susieout_check_psd <- function(x, ridge = 1e-8, exact_limit = 1000L) {
    if (inherits(x, "sparseMatrix")) {
        scale <- max(1, max(abs(.susieout_diag(x))))
        stabilized <- Matrix::forceSymmetric(
            x + Matrix::Diagonal(nrow(x), ridge * scale),
            uplo = "U"
        )
        factor <- tryCatch(
            Matrix::Cholesky(stabilized, LDL = FALSE, perm = TRUE),
            error = function(error) NULL
        )
        if (is.null(factor)) {
            stop(
                "sparse LD matrix is not positive semidefinite after ridge stabilization",
                call. = FALSE
            )
        }
        return(invisible(TRUE))
    }
    if (nrow(x) <= exact_limit) {
        .susieout_psd_eigen(x, ridge)
        return(invisible(TRUE))
    }
    scale <- max(1, max(abs(diag(x))))
    factor <- tryCatch(
        chol((x + t(x)) / 2 + diag(ridge * scale, nrow(x))),
        error = function(error) NULL
    )
    if (is.null(factor)) {
        stop(
            "LD matrix is not positive semidefinite after ridge stabilization",
            call. = FALSE
        )
    }
    invisible(TRUE)
}

.susieout_matrix_power <- function(x, power, ridge = 1e-8) {
    eig <- .susieout_psd_eigen(x, ridge)
    tcrossprod(sweep(eig$vectors, 2L, eig$values^power, `*`), eig$vectors)
}

.susieout_participation_rank <- function(x) {
    denominator <- sum(x * x)
    if (denominator <= 0) {
        return(0)
    }
    as.numeric(sum(.susieout_diag(x))^2 / denominator)
}

.susieout_mvn_draw <- function(covariance) {
    eig <- .susieout_psd_eigen(covariance, ridge = 1e-10)
    as.numeric(
        eig$vectors %*% (sqrt(eig$values) * stats::rnorm(nrow(covariance)))
    )
}

.susieout_pip <- function(alpha) {
    if (ncol(alpha) == 0L) {
        return(numeric())
    }
    bounded <- pmin(pmax(alpha, 0), 1 - .Machine$double.eps)
    1 - exp(colSums(log1p(-bounded)))
}

.susieout_complement <- function(allele) {
    translate <- c(A = "T", T = "A", C = "G", G = "C")
    result <- unname(translate[toupper(as.character(allele))])
    result[is.na(result)] <- ""
    result
}
