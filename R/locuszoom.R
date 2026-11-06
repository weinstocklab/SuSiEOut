#' LocusZoom-style plot of a SuSiE-OUT fit
#'
#' Draws stacked, position-aligned panels for one locus:
#'
#' 1. **Association.** Tested variants at their position, with
#'    \eqn{-\log_{10} p}{-log10(p)} computed from the z score and colour showing
#'    the squared correlation (\eqn{r^2}{r^2}) with the lead variant, as in
#'    LocusZoom. Omitted candidate variants have no z score; ticks below the
#'    axis mark where they lie.
#' 2. **Posterior inclusion probability.** The PIP of every tested (circle)
#'    and omitted (triangle) variant, filled by the 95% credible set that
#'    contains it.
#' 3. **Genes** (optional), when `genes` is supplied.
#'
#' The subtitle reports the posterior probability that at least one causal
#' variant was omitted from the tested set.
#'
#' Requires the suggested packages \pkg{ggplot2} and \pkg{cowplot};
#' \pkg{ggrastr} is used to rasterize large point layers when available.
#'
#' @param fit A `susieout_fit` from [susieout_fit()].
#' @param positions Base-pair positions for every tested and omitted variant,
#'   either a numeric vector named by variant ID or a data frame with columns
#'   `variant_id` and `position`. If `NULL`, positions are parsed from IDs of
#'   the form `chr1:12345:A:G` or `chr1-12345-A-G`; if that is not possible
#'   variants are placed by their order and a message is shown.
#' @param genes Optional data frame of genes with columns `gene_name`, `start`
#'   and `end`, and optionally `chrom`. Genes outside the plotted range are
#'   dropped. Requires genomic positions.
#' @param chrom Chromosome of the locus (for example `"chr1"`), used to select
#'   genes and in the title. Needed when `genes` has several chromosomes and
#'   `positions` was supplied explicitly; otherwise inferred from variant IDs.
#' @param lead ID of the tested variant that anchors the \eqn{r^2}{r^2}
#'   colouring. Defaults to the tested variant with the largest `|z|`.
#' @param ld If `FALSE`, skip LD colouring and draw all tested variants in one
#'   colour.
#' @param coverage Credible-set coverage passed to [susieout_credible_sets()].
#' @param min_probability Minimum probability that a single-effect component
#'   is non-null for its credible set to be shown. Components without
#'   non-null support yield large, uninformative sets that would otherwise
#'   colour most of the locus.
#' @param title Optional plot title. Defaults to the locus range when
#'   positions are genomic.
#' @return A \pkg{ggplot2} object (a \pkg{cowplot} grid of panels).
#' @seealso [plot.susieout_fit()] calls this function by default.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("cowplot", quietly = TRUE)) {
#'     locus <- susieout_simulate_locus(1, architecture = "entirely_omitted")
#'     fit <- susieout_fit(
#'         locus$z, locus$n, locus$R_analysis, locus$wgs_ld,
#'         locus$examined_ids
#'     )
#'     # Simulated loci have no coordinates; supply them for a genomic axis.
#'     ids <- c(fit$examined_ids, fit$excluded_ids)
#'     positions <- stats::setNames(1e6 + 2500 * seq_along(ids), ids)
#'     susieout_locuszoom(fit, positions = positions)
#' }
#' @export
susieout_locuszoom <- function(
    fit,
    positions = NULL,
    genes = NULL,
    chrom = NULL,
    lead = NULL,
    ld = TRUE,
    coverage = 0.95,
    min_probability = 0.5,
    title = NULL
) {
    if (!inherits(fit, "susieout_fit")) {
        stop("fit must be a susieout_fit", call. = FALSE)
    }
    require_plot_packages()

    tested_ids <- fit$examined_ids
    omitted_ids <- fit$excluded_ids
    coords <- resolve_positions(c(tested_ids, omitted_ids), positions)
    if (!is.null(chrom)) {
        coords$chrom <- as.character(chrom)
    }
    pos <- coords$position
    tested_pos <- pos[seq_along(tested_ids)]
    omitted_pos <- pos[length(tested_ids) + seq_along(omitted_ids)]

    # Association panel data ------------------------------------------------
    z <- as.numeric(fit$z)
    neglog10p <- -(stats::pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10)
    lead_index <- if (is.null(lead)) {
        which.max(abs(z))
    } else {
        match(lead, tested_ids)
    }
    if (length(lead_index) != 1L || is.na(lead_index)) {
        stop("lead must be the ID of a tested variant", call. = FALSE)
    }
    r2 <- rep(NA_real_, length(z))
    if (isTRUE(ld) && !is.null(fit$R_analysis)) {
        r2 <- as.numeric(fit$R_analysis[, lead_index])^2
    }
    tested <- data.frame(
        variant_id = tested_ids,
        position = tested_pos,
        neglog10p = neglog10p,
        r2_bin = ld_bin(r2),
        pip = as.numeric(fit$pip_in),
        stringsAsFactors = FALSE
    )
    omitted <- data.frame(
        variant_id = omitted_ids,
        position = omitted_pos,
        pip = as.numeric(fit$pip_out),
        stringsAsFactors = FALSE
    )

    # Credible sets ---------------------------------------------------------
    set_label <- credible_set_labels(
        fit,
        c(tested_ids, omitted_ids),
        coverage,
        min_probability
    )
    tested$credible_set <- set_label[seq_along(tested_ids)]
    omitted$credible_set <- set_label[
        length(tested_ids) + seq_along(omitted_ids)
    ]

    xlim <- range(pos)
    pad <- max(diff(xlim) * 0.02, 1)
    xlim <- xlim + c(-pad, pad)

    p_assoc <- assoc_panel(
        tested,
        omitted,
        tested_pos[lead_index],
        xlim,
        coords
    )
    p_pip <- pip_panel(
        tested,
        omitted,
        xlim,
        coords,
        fit,
        has_genes = !is.null(genes)
    )

    panels <- list(p_assoc, p_pip)
    heights <- c(3, 2)
    if (!is.null(genes)) {
        if (!coords$genomic) {
            stop("genes require genomic positions", call. = FALSE)
        }
        gene_plot <- gene_panel(genes, xlim, coords$chrom)
        panels <- c(panels, list(gene_plot$plot))
        heights <- c(heights, 0.5 + 0.45 * gene_plot$tracks)
    }

    main <- title
    if (is.null(main) && coords$genomic) {
        range_text <- paste(
            format(round(range(pos)), big.mark = ",", scientific = FALSE),
            collapse = "-"
        )
        main <- if (is.na(coords$chrom)) {
            paste(range_text, "bp")
        } else {
            paste0(coords$chrom, ":", range_text)
        }
    }
    header <- cowplot::ggdraw() +
        cowplot::draw_label(
            main %||% "",
            x = 0.01,
            hjust = 0,
            fontface = "bold",
            size = 13
        )
    subtitle <- cowplot::ggdraw() +
        cowplot::draw_label(
            sprintf(
                "P(at least one causal variant omitted) = %.2f",
                fit$any_causal_omitted
            ),
            x = 0.01,
            hjust = 0,
            size = 10,
            colour = "grey30"
        )
    body <- cowplot::plot_grid(
        plotlist = panels,
        ncol = 1,
        align = "v",
        axis = "lr",
        rel_heights = heights
    )
    cowplot::plot_grid(
        header,
        subtitle,
        body,
        ncol = 1,
        rel_heights = c(0.06, 0.05, 1)
    )
}

# ggplot2 supplies the `.data` pronoun inside aes(); it is a suggested
# dependency, so declare the name instead of importing it.
utils::globalVariables(".data")

# Palettes ------------------------------------------------------------------

ld_levels <- c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0")
ld_colors <- c(
    "0.0-0.2" = "#3B4CC0",
    "0.2-0.4" = "#8FB8E6",
    "0.4-0.6" = "#58B368",
    "0.6-0.8" = "#F4A53A",
    "0.8-1.0" = "#D7301F"
)
# Okabe-Ito colours without orange, which marks omitted variants elsewhere.
set_colors <- c(
    "#0072B2",
    "#CC79A7",
    "#009E73",
    "#56B4E9",
    "#D55E00",
    "#F0E442",
    "#999999"
)

ld_bin <- function(r2) {
    bin <- cut(
        r2,
        breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
        labels = ld_levels,
        right = FALSE
    )
    factor(as.character(bin), levels = ld_levels)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

require_plot_packages <- function() {
    for (pkg in c("ggplot2", "cowplot")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop(
                "Package '",
                pkg,
                "' is required for locus plots; install it with ",
                "install.packages(\"",
                pkg,
                "\")",
                call. = FALSE
            )
        }
    }
    invisible(TRUE)
}

# Positions -----------------------------------------------------------------

resolve_positions <- function(ids, positions) {
    if (!is.null(positions)) {
        if (is.data.frame(positions)) {
            if (!all(c("variant_id", "position") %in% names(positions))) {
                stop(
                    "positions data frame needs columns variant_id and position",
                    call. = FALSE
                )
            }
            lookup <- stats::setNames(
                as.numeric(positions$position),
                as.character(positions$variant_id)
            )
        } else {
            if (is.null(names(positions))) {
                stop("positions must be named by variant ID", call. = FALSE)
            }
            lookup <- stats::setNames(as.numeric(positions), names(positions))
        }
        missing_ids <- setdiff(ids, names(lookup))
        if (length(missing_ids)) {
            stop(
                "positions is missing ",
                length(missing_ids),
                " variant(s), including: ",
                paste(utils::head(missing_ids, 3L), collapse = ", "),
                call. = FALSE
            )
        }
        return(list(
            position = unname(lookup[ids]),
            genomic = TRUE,
            chrom = NA_character_
        ))
    }
    pattern <- "^(chr)?([0-9XYM]+)[:_-]([0-9]+)([:_-].*)?$"
    if (all(grepl(pattern, ids))) {
        return(list(
            position = as.numeric(sub(pattern, "\\3", ids)),
            genomic = TRUE,
            chrom = {
                chroms <- unique(sub(pattern, "chr\\2", ids))
                if (length(chroms) == 1L) chroms else NA_character_
            }
        ))
    }
    message(
        "Could not parse positions from variant IDs; placing variants by ",
        "order. Supply `positions` for a genomic axis."
    )
    list(
        position = as.numeric(seq_along(ids)),
        genomic = FALSE,
        chrom = NA_character_
    )
}

# Credible sets -------------------------------------------------------------

credible_set_labels <- function(fit, ids, coverage, min_probability) {
    label <- rep(NA_character_, length(ids))
    sets <- susieout_credible_sets(fit, coverage = coverage)
    shown <- 0L
    for (set in sets) {
        if (!nrow(set) || set$nonnull_probability[1L] < min_probability) {
            next
        }
        shown <- shown + 1L
        members <- match(set$variant_id, ids)
        label[members[is.na(label[members])]] <- paste0("CS", shown)
    }
    label
}

# Panels --------------------------------------------------------------------

panel_theme <- function(show_x) {
    ggplot2::theme_classic(base_size = 11) +
        ggplot2::theme(
            axis.text.x = if (show_x) {
                ggplot2::element_text()
            } else {
                ggplot2::element_blank()
            },
            axis.ticks.x = if (show_x) {
                ggplot2::element_line()
            } else {
                ggplot2::element_blank()
            },
            axis.title.x = if (show_x) {
                ggplot2::element_text()
            } else {
                ggplot2::element_blank()
            },
            legend.key.size = ggplot2::unit(0.4, "cm"),
            legend.title = ggplot2::element_text(size = 9),
            legend.text = ggplot2::element_text(size = 8),
            plot.margin = ggplot2::margin(2, 6, 2, 6)
        )
}

x_scale <- function(xlim, coords) {
    if (coords$genomic) {
        ggplot2::scale_x_continuous(
            limits = xlim,
            expand = c(0, 0),
            labels = function(x) {
                format(x / 1e6, trim = TRUE, scientific = FALSE)
            }
        )
    } else {
        ggplot2::scale_x_continuous(limits = xlim, expand = c(0, 0))
    }
}

x_label <- function(coords) {
    if (coords$genomic) "Position (Mb)" else "Variant order"
}

point_layer <- function(data, mapping, ...) {
    layer <- ggplot2::geom_point(data = data, mapping = mapping, ...)
    if (nrow(data) > 2000L && requireNamespace("ggrastr", quietly = TRUE)) {
        layer <- ggrastr::rasterise(layer, dpi = 300)
    }
    layer
}

assoc_panel <- function(tested, omitted, lead_position, xlim, coords) {
    ymax <- max(tested$neglog10p, 1, na.rm = TRUE)
    tick_y <- -0.05 * ymax
    has_ld <- any(!is.na(tested$r2_bin))
    plain <- tested[is.na(tested$r2_bin), , drop = FALSE]
    colored <- tested[!is.na(tested$r2_bin), , drop = FALSE]
    lead <- tested[tested$position == lead_position, , drop = FALSE][1L, ]

    p <- ggplot2::ggplot()
    if (ymax >= 5) {
        p <- p +
            ggplot2::geom_hline(
                yintercept = -log10(5e-8),
                linetype = "dashed",
                colour = "grey60",
                linewidth = 0.3
            )
    }
    if (nrow(plain)) {
        p <- p +
            point_layer(
                plain,
                ggplot2::aes(x = .data$position, y = .data$neglog10p),
                colour = "grey55",
                size = 1.6,
                alpha = 0.85
            )
    }
    if (nrow(colored)) {
        p <- p +
            point_layer(
                colored,
                ggplot2::aes(
                    x = .data$position,
                    y = .data$neglog10p,
                    colour = .data$r2_bin
                ),
                size = 1.8,
                alpha = 0.95
            )
    }
    p <- p +
        ggplot2::geom_point(
            data = lead,
            ggplot2::aes(x = .data$position, y = .data$neglog10p),
            shape = 23,
            size = 3.4,
            fill = "#7B3294",
            colour = "black",
            stroke = 0.5
        )
    if (nrow(omitted)) {
        p <- p +
            ggplot2::geom_point(
                data = omitted,
                ggplot2::aes(x = .data$position, y = tick_y),
                shape = 17,
                size = 1.8,
                colour = "#C05621"
            )
    }
    p +
        ggplot2::scale_colour_manual(
            values = ld_colors,
            limits = ld_levels,
            drop = FALSE,
            name = expression(r^2 ~ "with lead"),
            guide = if (has_ld) ggplot2::guide_legend(order = 1) else "none"
        ) +
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(mult = c(0.02, 0.08)),
            breaks = function(l) pretty(c(0, l[2L]))[pretty(c(0, l[2L])) >= 0]
        ) +
        x_scale(xlim, coords) +
        ggplot2::labs(y = expression(-log[10](italic(p)))) +
        panel_theme(show_x = FALSE)
}

pip_panel <- function(tested, omitted, xlim, coords, fit, has_genes) {
    tested$state <- "Tested"
    omitted$state <- "Omitted"
    cols <- c("variant_id", "position", "pip", "credible_set", "state")
    both <- rbind(tested[, cols], omitted[, cols])
    sets <- sort(unique(stats::na.omit(both$credible_set)))
    fills <- stats::setNames(set_colors[seq_along(sets)], sets)
    fills <- c(fills, "none" = "white")
    both$fill <- ifelse(is.na(both$credible_set), "none", both$credible_set)
    both$state <- factor(both$state, levels = c("Tested", "Omitted"))
    both <- both[order(both$pip), , drop = FALSE]

    ggplot2::ggplot(
        both,
        ggplot2::aes(
            x = .data$position,
            y = .data$pip,
            fill = .data$fill,
            shape = .data$state,
            colour = .data$state
        )
    ) +
        ggplot2::geom_point(size = 2.2, stroke = 0.5) +
        ggplot2::scale_colour_manual(
            values = c(Tested = "grey20", Omitted = "#C05621"),
            name = "Variant",
            drop = FALSE,
            guide = ggplot2::guide_legend(
                order = 2,
                override.aes = list(fill = "grey85")
            )
        ) +
        ggplot2::scale_shape_manual(
            values = c(Tested = 21, Omitted = 24),
            name = "Variant",
            drop = FALSE,
            guide = ggplot2::guide_legend(order = 2)
        ) +
        ggplot2::scale_fill_manual(
            values = fills,
            breaks = sets,
            name = "Credible set",
            guide = if (length(sets)) {
                ggplot2::guide_legend(
                    order = 3,
                    override.aes = list(shape = 21)
                )
            } else {
                "none"
            }
        ) +
        ggplot2::scale_y_continuous(
            limits = c(0, 1),
            breaks = c(0, 0.5, 1),
            expand = ggplot2::expansion(mult = c(0.03, 0.05))
        ) +
        x_scale(xlim, coords) +
        ggplot2::labs(x = x_label(coords), y = "PIP") +
        panel_theme(show_x = !has_genes)
}

gene_panel <- function(genes, xlim, chrom) {
    needed <- c("gene_name", "start", "end")
    if (!is.data.frame(genes) || !all(needed %in% names(genes))) {
        stop(
            "genes must be a data frame with columns gene_name, start and end",
            call. = FALSE
        )
    }
    if ("chrom" %in% names(genes)) {
        if (is.na(chrom)) {
            if (length(unique(genes$chrom)) > 1L) {
                stop(
                    "genes spans several chromosomes; set `chrom`",
                    call. = FALSE
                )
            }
        } else {
            genes <- genes[genes$chrom == chrom, , drop = FALSE]
        }
    }
    genes <- genes[
        genes$end >= xlim[1L] & genes$start <= xlim[2L],
        ,
        drop = FALSE
    ]
    base <- ggplot2::ggplot() +
        x_scale(xlim, list(genomic = TRUE)) +
        ggplot2::labs(x = "Position (Mb)", y = NULL) +
        ggplot2::theme_classic(base_size = 11) +
        ggplot2::theme(
            axis.line.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(2, 6, 2, 6)
        )
    if (!nrow(genes)) {
        return(list(
            plot = base + ggplot2::scale_y_continuous(limits = c(0, 1)),
            tracks = 1
        ))
    }
    genes$track <- assign_tracks(
        genes$start,
        genes$end,
        gap = diff(xlim) * 0.08
    )
    genes$y <- -genes$track
    genes$mid <- pmin(pmax((genes$start + genes$end) / 2, xlim[1L]), xlim[2L])
    plot <- base +
        ggplot2::geom_segment(
            data = genes,
            ggplot2::aes(
                x = pmax(.data$start, xlim[1L]),
                xend = pmin(.data$end, xlim[2L]),
                y = .data$y,
                yend = .data$y
            ),
            linewidth = 2,
            colour = "#2c3e50",
            lineend = "butt"
        ) +
        ggplot2::geom_text(
            data = genes,
            ggplot2::aes(
                x = .data$mid,
                y = .data$y + 0.4,
                label = .data$gene_name
            ),
            size = 2.8,
            vjust = 0,
            fontface = "italic",
            colour = "#2c3e50"
        ) +
        ggplot2::scale_y_continuous(
            limits = c(-max(genes$track) - 0.2, 0.9),
            expand = c(0, 0)
        )
    list(plot = plot, tracks = max(genes$track))
}

# Greedy interval packing: place each gene on the first track whose last gene
# ends at least `gap` before it starts.
assign_tracks <- function(start, end, gap = 0) {
    order <- order(start)
    track_end <- numeric()
    track <- integer(length(start))
    for (i in order) {
        free <- which(track_end + gap < start[i])
        if (length(free)) {
            chosen <- free[1L]
        } else {
            chosen <- length(track_end) + 1L
            track_end[chosen] <- -Inf
        }
        track[i] <- chosen
        track_end[chosen] <- end[i]
    }
    track
}
