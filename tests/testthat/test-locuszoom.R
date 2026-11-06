simulated_fit <- function(seed = 10, architecture = "entirely_omitted") {
    locus <- susieout_simulate_locus(seed, architecture = architecture)
    fit <- susieout_fit(
        locus$z,
        locus$n,
        locus$R_analysis,
        locus$wgs_ld,
        locus$examined_ids,
        L = 10
    )
    ids <- c(fit$examined_ids, fit$excluded_ids)
    list(
        fit = fit,
        positions = stats::setNames(1e6 + 2500 * seq_along(ids), ids)
    )
}

test_that("locus plot builds for explicit positions and with genes", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("cowplot")
    x <- simulated_fit()
    plot <- susieout_locuszoom(x$fit, positions = x$positions)
    expect_s3_class(plot, "ggplot")

    genes <- data.frame(
        chrom = "chr1",
        gene_name = c("A", "B", "C"),
        start = c(1.01e6, 1.02e6, 1.09e6),
        end = c(1.05e6, 1.04e6, 1.12e6)
    )
    with_genes <- susieout_locuszoom(
        x$fit,
        positions = x$positions,
        genes = genes
    )
    expect_s3_class(with_genes, "ggplot")
    expect_error(
        suppressMessages(susieout_locuszoom(x$fit, genes = genes)),
        "genomic positions"
    )

    # Genes on other chromosomes are dropped when the chromosome is known and
    # explicitly ambiguous otherwise.
    other <- rbind(genes, transform(genes, chrom = "chr2", gene_name = "X"))
    expect_error(
        susieout_locuszoom(x$fit, positions = x$positions, genes = other),
        "set `chrom`"
    )
    expect_s3_class(
        susieout_locuszoom(
            x$fit,
            positions = x$positions,
            genes = other,
            chrom = "chr1"
        ),
        "ggplot"
    )
})

test_that("positions are parsed from IDs, or variants are placed by order", {
    ids <- c("chr1:100:A:G", "chr1:250:C:T", "chr1:400:G:A")
    parsed <- SuSiEOut:::resolve_positions(ids, NULL)
    expect_equal(parsed$position, c(100, 250, 400))
    expect_equal(parsed$chrom, "chr1")
    expect_true(parsed$genomic)

    expect_message(
        ordered <- SuSiEOut:::resolve_positions(c("a", "b"), NULL),
        "placing variants by order"
    )
    expect_false(ordered$genomic)
    expect_equal(ordered$position, c(1, 2))

    frame <- data.frame(variant_id = ids, position = c(5, 6, 7))
    expect_equal(SuSiEOut:::resolve_positions(ids, frame)$position, c(5, 6, 7))
    expect_error(
        SuSiEOut:::resolve_positions(ids, c(chr1 = 1)),
        "missing"
    )
})

test_that("locus plot validates its inputs", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("cowplot")
    x <- simulated_fit()
    expect_error(susieout_locuszoom(list()), "susieout_fit")
    expect_error(
        susieout_locuszoom(x$fit, positions = x$positions, lead = "nope"),
        "tested variant"
    )
})

test_that("only credible sets with non-null support are labelled", {
    label_ids <- function(fit) c(fit$examined_ids, fit$excluded_ids)
    weak <- simulated_fit(seed = 2)$fit
    expect_true(all(is.na(
        SuSiEOut:::credible_set_labels(weak, label_ids(weak), 0.95, 0.5)
    )))
    strong <- simulated_fit(seed = 10, architecture = "contained")$fit
    labels <- SuSiEOut:::credible_set_labels(
        strong,
        label_ids(strong),
        0.95,
        0.5
    )
    expect_true(any(labels == "CS1", na.rm = TRUE))
})

test_that("gene tracks do not overlap", {
    tracks <- SuSiEOut:::assign_tracks(
        start = c(0, 10, 60, 5),
        end = c(50, 40, 90, 30),
        gap = 1
    )
    for (track in unique(tracks)) {
        keep <- tracks == track
        ordered <- order(c(0, 10, 60, 5)[keep])
        starts <- c(0, 10, 60, 5)[keep][ordered]
        ends <- c(50, 40, 90, 30)[keep][ordered]
        expect_true(all(starts[-1L] > ends[-length(ends)]))
    }
})

test_that("plot() dispatches to the locus plot and keeps the legacy types", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("cowplot")
    x <- simulated_fit()
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    expect_identical(
        plot(x$fit, positions = x$positions),
        x$fit
    )
    expect_identical(plot(x$fit, type = "pip", top = 10), x$fit)
    expect_identical(plot(x$fit, type = "components"), x$fit)
})
