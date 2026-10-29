#' Construct an effect-size mixture prior
#'
#' The robust mixture is the primary prior and gives equal architecture mass to
#' neutral, stabilizing, and recent directional-optimum-shift components.
#'
#' @param architecture One of `robust_mixture`, `generic`, `neutral`,
#'   `stabilizing`, or `directional`.
#' @return A normalized list of normal-mixture weights, means, and variances.
#' @export
susieout_effect_prior <- function(architecture = c(
        "robust_mixture", "generic", "neutral", "stabilizing", "directional")) {
    architecture <- match.arg(architecture)
    definitions <- switch(
        architecture,
        generic = list(
            weights = c(0.45, 0.35, 0.20),
            means = c(0, 0, 0),
            variances = c(2.5e-5, 2.5e-4, 2.5e-3),
            labels = c("small", "medium", "large")),
        neutral = list(
            weights = c(0.65, 0.30, 0.05),
            means = c(0, 0, 0),
            variances = c(1e-5, 1e-4, 1e-3),
            labels = c("neutral_small", "neutral_medium", "neutral_large")),
        stabilizing = list(
            weights = c(0.35, 0.45, 0.20),
            means = c(0, 0, 0),
            variances = c(2.5e-5, 4e-4, 4e-3),
            labels = c("stabilizing_small", "stabilizing_medium",
                       "stabilizing_large")),
        directional = list(
            weights = c(0.30, 0.35, 0.35),
            means = c(0, 0.012, -0.012),
            variances = c(5e-5, 5e-4, 5e-4),
            labels = c("directional_symmetric", "directional_derived_positive",
                       "directional_derived_negative")),
        robust_mixture = list(
            weights = c(0.65, 0.30, 0.05,
                        0.35, 0.45, 0.20,
                        0.30, 0.35, 0.35) / 3,
            means = c(0, 0, 0, 0, 0, 0, 0, 0.012, -0.012),
            variances = c(1e-5, 1e-4, 1e-3,
                          2.5e-5, 4e-4, 4e-3,
                          5e-5, 5e-4, 5e-4),
            labels = c("neutral_small", "neutral_medium", "neutral_large",
                       "stabilizing_small", "stabilizing_medium",
                       "stabilizing_large", "directional_symmetric",
                       "directional_positive", "directional_negative"))
    )
    definitions$architecture <- architecture
    class(definitions) <- "susieout_effect_prior"
    definitions
}

.susieout_validate_effect_prior <- function(effect_prior) {
    if (is.character(effect_prior) && length(effect_prior) == 1L) {
        effect_prior <- susieout_effect_prior(effect_prior)
    }
    if (!is.list(effect_prior)) {
        stop("effect_prior must be an architecture name or a mixture list",
             call. = FALSE)
    }
    weights <- .susieout_validate_vector(effect_prior$weights, "effect_prior$weights")
    variances <- .susieout_validate_vector(
        effect_prior$variances, "effect_prior$variances", length(weights))
    means <- effect_prior$means %||% rep(0, length(weights))
    means <- .susieout_validate_vector(means, "effect_prior$means", length(weights))
    if (any(weights < 0) || sum(weights) <= 0) {
        stop("effect-prior weights must be nonnegative with positive sum",
             call. = FALSE)
    }
    if (any(variances <= 0)) {
        stop("effect-prior variances must be positive", call. = FALSE)
    }
    labels <- effect_prior$labels %||% paste0("component_", seq_along(weights))
    list(weights = weights / sum(weights), means = means,
         variances = variances, labels = labels,
         architecture = effect_prior$architecture %||% "custom")
}

