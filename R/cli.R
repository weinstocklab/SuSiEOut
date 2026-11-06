.susieout_parse_cli <- function(args) {
    value <- list(
        L = 10L,
        lambda_total = 1,
        out_fraction = 0.10,
        summary = NULL
    )
    index <- 1L
    while (index <= length(args)) {
        key <- args[index]
        if (key %in% c("--help", "-h")) {
            value$help <- TRUE
            index <- index + 1L
            next
        }
        if (index == length(args)) {
            stop("missing value for ", key, call. = FALSE)
        }
        argument <- args[index + 1L]
        if (key == "--input") {
            value$input <- argument
        } else if (key == "--output") {
            value$output <- argument
        } else if (key == "--summary") {
            value$summary <- argument
        } else if (key == "--L") {
            value$L <- as.integer(argument)
        } else if (key == "--lambda-total") {
            value$lambda_total <- as.numeric(argument)
        } else if (key == "--out-fraction") {
            value$out_fraction <- as.numeric(argument)
        } else {
            stop("unknown argument: ", key, call. = FALSE)
        }
        index <- index + 2L
    }
    value
}

#' Command-line entry point for batch locus analysis
#'
#' The input is an RDS list containing the named arguments of [susieout_fit()].
#' The output is a complete fit RDS, with an optional tab-separated support file.
#' The installed script `exec/susieout` calls this function.
#'
#' @param args Character vector of command-line arguments: `--input`,
#'   `--output`, and optionally `--summary`, `--L`, `--lambda-total`, and
#'   `--out-fraction`. Use `--help` for usage.
#' @return Invisibly, `0` on success. The fit and optional support table are
#'   written to the paths given.
#' @examples
#' susieout_cli("--help")
#' @export
susieout_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
    options <- .susieout_parse_cli(args)
    if (isTRUE(options$help)) {
        cat(paste(
            "Usage: susieout --input LOCUS.rds --output FIT.rds",
            paste(
                "[--summary SUPPORT.tsv] [--L 10] [--lambda-total 1]",
                "[--out-fraction 0.10]\n"
            )
        ))
        return(invisible(0L))
    }
    if (is.null(options$input) || is.null(options$output)) {
        stop("--input and --output are required; use --help", call. = FALSE)
    }
    locus <- readRDS(options$input)
    if (!is.list(locus)) {
        stop("input RDS must contain a locus list", call. = FALSE)
    }
    locus$L <- locus$L %||% options$L
    locus$count_prior <- locus$count_prior %||%
        susieout_count_prior(
            options$lambda_total,
            options$out_fraction
        )
    allowed <- names(formals(susieout_fit))
    arguments <- locus[intersect(names(locus), allowed)]
    fit <- do.call(susieout_fit, arguments)
    saveRDS(fit, options$output)
    if (!is.null(options$summary)) {
        support <- susieout_candidate_support(fit)
        support$locus_any_causal_omitted <- fit$any_causal_omitted
        support$candidate_set_complete <- fit$candidate_set_complete
        utils::write.table(
            support,
            options$summary,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )
    }
    invisible(0L)
}
