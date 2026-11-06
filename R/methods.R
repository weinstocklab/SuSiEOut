#' Variant support under SuSiE-OUT
#'
#' Lists every tested variant and every omitted variant (OUT state) with its
#' posterior inclusion probability (PIP), sorted from highest to lowest.
#'
#' @param fit A `susieout_fit` object from [susieout_fit()] or one of the
#'   grid-averaged fitting functions.
#' @return A data frame with columns `variant_id`, `state` (`"inside"` for
#'   tested variants, `"outside"` for omitted variants), `pip`, and
#'   `maximum_component_probability`, the largest probability that any single
#'   effect assigns to the variant.
#' @examples
#' locus <- susieout_simulate_locus(1, architecture = "entirely_omitted")
#' fit <- susieout_fit(
#'   locus$z, locus$n, locus$R_analysis, locus$wgs_ld, locus$examined_ids
#' )
#' susieout_candidate_support(fit)[1:5, ]
#' @export
susieout_candidate_support <- function(fit) {
    if (!inherits(fit, "susieout_fit")) {
        stop("fit must be a susieout_fit", call. = FALSE)
    }
    inside <- data.frame(
        variant_id = fit$examined_ids,
        state = rep("inside", length(fit$examined_ids)),
        pip = unname(fit$pip_in),
        maximum_component_probability = if (length(fit$pip_in)) {
            apply(fit$alpha_in, 2L, max)
        } else {
            numeric()
        },
        stringsAsFactors = FALSE
    )
    outside <- data.frame(
        variant_id = fit$excluded_ids,
        state = rep("outside", length(fit$excluded_ids)),
        pip = unname(fit$pip_out),
        maximum_component_probability = if (length(fit$pip_out)) {
            apply(fit$alpha_out, 2L, max)
        } else {
            numeric()
        },
        stringsAsFactors = FALSE
    )
    result <- rbind(inside, outside)
    result[order(result$pip, decreasing = TRUE), , drop = FALSE]
}

#' SuSiE-OUT credible sets
#'
#' Credible sets are formed separately for each single-effect component and may
#' contain tested or omitted variants. Probabilities are normalized
#' conditional on a non-null sparse effect; the returned `nonnull_probability`
#' preserves the unconditional component support.
#'
#' @param fit A `susieout_fit` object.
#' @param coverage Target coverage of each credible set, in (0, 1).
#' @return A list of class `susieout_credible_sets` with one data frame per
#'   single-effect component. Each data frame lists the variants in the set
#'   with their `state` (`"inside"` or `"outside"`), unconditional
#'   `probability`, `conditional_probability` given a non-null effect, and
#'   `cumulative` conditional probability. Components with no non-null support
#'   give empty data frames. The coverage is stored as an attribute.
#' @examples
#' locus <- susieout_simulate_locus(1, architecture = "entirely_omitted")
#' fit <- susieout_fit(
#'   locus$z, locus$n, locus$R_analysis, locus$wgs_ld, locus$examined_ids
#' )
#' sets <- susieout_credible_sets(fit)
#' Filter(nrow, sets)
#' @export
susieout_credible_sets <- function(fit, coverage = 0.95) {
    if (!inherits(fit, "susieout_fit")) {
        stop("fit must be a susieout_fit", call. = FALSE)
    }
    if (
        !is.numeric(coverage) ||
            length(coverage) != 1L ||
            !is.finite(coverage) ||
            coverage <= 0 ||
            coverage >= 1
    ) {
        stop("coverage must lie in (0, 1)", call. = FALSE)
    }
    alpha <- cbind(fit$alpha_in, fit$alpha_out)
    ids <- c(fit$examined_ids, fit$excluded_ids)
    state <- c(
        rep("inside", length(fit$examined_ids)),
        rep("outside", length(fit$excluded_ids))
    )
    sets <- vector("list", nrow(alpha))
    for (component in seq_len(nrow(alpha))) {
        nonnull <- sum(alpha[component, ])
        if (nonnull <= .Machine$double.eps) {
            sets[[component]] <- data.frame(
                component = integer(),
                variant_id = character(),
                state = character(),
                probability = numeric(),
                conditional_probability = numeric(),
                cumulative = numeric(),
                nonnull_probability = numeric(),
                stringsAsFactors = FALSE
            )
            next
        }
        conditional <- alpha[component, ] / nonnull
        order <- order(conditional, decreasing = TRUE)
        cumulative <- cumsum(conditional[order])
        count <- which(cumulative >= coverage)[1L]
        selected <- order[seq_len(count)]
        sets[[component]] <- data.frame(
            component = component,
            variant_id = ids[selected],
            state = state[selected],
            probability = alpha[component, selected],
            conditional_probability = conditional[selected],
            cumulative = cumulative[seq_len(count)],
            nonnull_probability = nonnull,
            stringsAsFactors = FALSE
        )
    }
    structure(
        sets,
        coverage = coverage,
        class = c("susieout_credible_sets", "list")
    )
}

#' Summarize, print, plot, and predict from a SuSiE-OUT fit
#'
#' `summary()` collects locus-level quantities and the ten best-supported
#' variants; `print()` displays that summary. `plot()` draws a locus plot
#' ([susieout_locuszoom()], the default when \pkg{ggplot2} and \pkg{cowplot}
#' are installed), a bar chart of the highest PIPs coloured by tested (blue)
#' or omitted (orange) state, or the null, tested and OUT probabilities of each
#' single-effect component.
#' `predict()` returns fitted z scores, residuals, variant support, or the
#' probability that at least one causal variant was omitted.
#'
#' @param object,x A `susieout_fit` object.
#' @param type For `plot()`, `"locuszoom"`, `"pip"` or `"components"`. For
#'   `predict()`, one of `"z"`, `"residual"`, `"candidate_support"`, or
#'   `"outside_probability"`.
#' @param top Number of variants to show in the `"pip"` bar chart.
#' @param newdata Optional fingerprint matrix with one column per tested and
#'   omitted variant; if supplied with `type = "z"`, the fitted effects are
#'   projected onto it.
#' @param ... For `plot()`, further arguments passed to
#'   [susieout_locuszoom()] (`type = "locuszoom"`, for example `positions` and
#'   `genes`) or to [graphics::barplot()] (the other types). Ignored by the
#'   other methods.
#' @return `summary()` returns an object of class `summary.susieout_fit`.
#'   `print()` and `plot()` return their input invisibly. `predict()` returns
#'   a numeric vector, a data frame (for `"candidate_support"`), or a scalar
#'   (for `"outside_probability"`).
#' @name susieout_fit-methods
#' @examples
#' locus <- susieout_simulate_locus(1, architecture = "entirely_omitted")
#' fit <- susieout_fit(
#'   locus$z, locus$n, locus$R_analysis, locus$wgs_ld, locus$examined_ids
#' )
#' summary(fit)
#' plot(fit, type = "pip", top = 15)
#' predict(fit, type = "outside_probability")
NULL

#' @rdname susieout_fit-methods
#' @export
summary.susieout_fit <- function(object, ...) {
    support <- susieout_candidate_support(object)
    result <- list(
        backend = object$backend,
        examined_variants = length(object$examined_ids),
        excluded_variants = length(object$excluded_ids),
        any_causal_omitted = object$any_causal_omitted,
        candidate_set_complete = object$candidate_set_complete,
        top_variants = utils::head(support, 10L),
        component_outside = object$component_outside,
        component_null = object$component_null,
        infinitesimal = object$infinitesimal,
        ld_mismatch = object$ld_mismatch,
        diagnostics = object$diagnostics
    )
    class(result) <- "summary.susieout_fit"
    result
}

#' @rdname susieout_fit-methods
#' @export
print.summary.susieout_fit <- function(x, ...) {
    cat("SuSiE-OUT fine-mapping summary\n")
    cat("  backend:", x$backend, "\n")
    cat(
        "  variants:",
        x$examined_variants,
        "tested +",
        x$excluded_variants,
        "omitted\n"
    )
    cat(
        "  P(any causal omitted):",
        format(round(x$any_causal_omitted, 4), nsmall = 4),
        "\n"
    )
    cat(
        "  P(candidate set complete):",
        format(round(x$candidate_set_complete, 4), nsmall = 4),
        "\n"
    )
    cat(
        "  converged:",
        x$diagnostics$converged,
        "in",
        x$diagnostics$iterations,
        "iterations\n"
    )
    if (nrow(x$top_variants)) {
        cat("\nTop variant support:\n")
        print(x$top_variants, row.names = FALSE)
    }
    invisible(x)
}

#' @rdname susieout_fit-methods
#' @export
print.susieout_fit <- function(x, ...) {
    print(summary(x), ...)
    invisible(x)
}

#' @rdname susieout_fit-methods
#' @export
plot.susieout_fit <- function(
    x,
    type = c("locuszoom", "pip", "components"),
    top = 30L,
    ...
) {
    if (missing(type)) {
        can_plot <- requireNamespace("ggplot2", quietly = TRUE) &&
            requireNamespace("cowplot", quietly = TRUE)
        type <- if (can_plot) "locuszoom" else "pip"
    } else {
        type <- match.arg(type)
    }
    if (type == "locuszoom") {
        print(susieout_locuszoom(x, ...))
    } else if (type == "pip") {
        support <- utils::head(susieout_candidate_support(x), top)
        colors <- ifelse(support$state == "inside", "#2B6CB0", "#C05621")
        graphics::barplot(
            support$pip,
            names.arg = support$variant_id,
            las = 2,
            col = colors,
            border = NA,
            ylab = "Posterior inclusion probability",
            main = "SuSiE-OUT variant support",
            ...
        )
        graphics::abline(h = 0.9, lty = 3, col = "grey40")
    } else {
        matrix <- rbind(
            outside = x$component_outside,
            null = x$component_null,
            inside = pmax(0, 1 - x$component_outside - x$component_null)
        )
        graphics::barplot(
            matrix,
            beside = FALSE,
            ylim = c(0, 1),
            col = c("#C05621", "grey80", "#2B6CB0"),
            xlab = "Single-effect component",
            ylab = "Posterior probability",
            main = "Component state responsibilities",
            ...
        )
    }
    invisible(x)
}

#' @rdname susieout_fit-methods
#' @export
predict.susieout_fit <- function(
    object,
    newdata = NULL,
    type = c("z", "residual", "candidate_support", "outside_probability"),
    ...
) {
    type <- match.arg(type)
    if (type == "candidate_support") {
        return(susieout_candidate_support(object))
    }
    if (type == "outside_probability") {
        return(object$any_causal_omitted)
    }
    if (type == "residual") {
        return(object$residual %||% object$z - object$fitted)
    }
    if (is.null(newdata)) {
        return(object$fitted %||% object$z)
    }
    if (is.list(object$fingerprints)) {
        expected <- length(object$examined_ids) + length(object$excluded_ids)
    } else {
        expected <- ncol(object$fingerprints)
    }
    newdata <- .susieout_validate_matrix(
        newdata,
        "newdata",
        ncol = expected
    )
    alpha <- cbind(object$alpha_in, object$alpha_out)
    joint_beta <- colSums(alpha * object$component_beta_mean)
    sqrt(object$n) * as.numeric(newdata %*% joint_beta)
}
