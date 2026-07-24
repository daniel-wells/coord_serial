library(ggplot2)

test_that("sort_domains works correctly", {
    # Numeric sorting
    expect_equal(sort_domains(c("10", "1", "2")), c("1", "2", "10"))

    # "chr" prefix handling (case-insensitive)
    expect_equal(sort_domains(c("chr2", "chr1", "chr10")), c("chr1", "chr2", "chr10"))
    expect_equal(sort_domains(c("Chr2", "chr1")), c("chr1", "Chr2"))

    # Genomic domains (X, Y, MT)
    expect_equal(sort_domains(c("2", "X", "1", "Y", "MT")), c("1", "2", "X", "Y", "MT"))
    expect_equal(sort_domains(c("chrMT", "chrY", "chrX")), c("chrX", "chrY", "chrMT"))

    # Combinations and unknown strings
    expect_equal(sort_domains(c("Z", "1", "B", "A")), c("1", "A", "B", "Z"))
})

test_that("compute_serial_layout works correctly", {
    positions <- 1:100
    domains <- rep(c("A", "B"), each = 50)

    layout <- compute_serial_layout(positions, domains, domain_gap = 0.1)

    expect_type(layout, "list")
    expect_named(layout, c("domain_levels", "offsets", "midpoints", "total_range"))

    # "A": positions 1-50, length 49. min=1, max=50
    # "B": positions 51-100, length 49. min=51, max=100
    # Total length = 49 + 49 = 98
    # Gap = 0.1 * 98 = 9.8

    # Offsets
    # A: offset = 0 - min(A) = -1
    # B: offset = (0 + 49 + 9.8) - min(B) = 58.8 - 51 = 7.8
    expect_equal(layout$offsets[["A"]], -1)
    expect_equal(layout$offsets[["B"]], 7.8)

    # Midpoints
    # A: midpoint = 0 + 49/2 = 24.5
    # B: midpoint = 58.8 + 49/2 = 58.8 + 24.5 = 83.3
    expect_equal(layout$midpoints[["A"]], 24.5)
    expect_equal(layout$midpoints[["B"]], 83.3)
})

test_that("compute_serial_layout handles edge cases", {
    # Single point domains
    expect_no_error(
        layout <- compute_serial_layout(c(1, 10), c("A", "B"))
    )
    expect_equal(layout$offsets[["A"]], -1)
    expect_equal(layout$offsets[["B"]], -10) # length A = 0, so cumul for B = 0

    # Manual domain order
    layout_ord <- compute_serial_layout(1:10, rep(c("A", "B"), each = 5), domain_order = c("B", "A"))
    expect_equal(layout_ord$domain_levels, c("B", "A"))

    # Empty input (robustly returns empty layout)
    expect_no_error(layout_empty <- compute_serial_layout(numeric(0), character(0)))
    expect_equal(length(layout_empty$domain_levels), 0)
    expect_equal(layout_empty$total_range, c(0, 0))
})

test_that("coord_serial creates a correct ggproto object", {
    coord <- coord_serial(domain_gap = 0.05, domain_order = c("X", "Y"))
    expect_s3_class(coord, "CoordSerial")
    expect_equal(coord$domain_gap, 0.05)
    expect_equal(coord$domain_order, c("X", "Y"))
})

test_that("ggplot_add.CoordSerial transforms data correctly", {
    df <- data.frame(
        pos = c(1, 10, 1, 5),
        chr = c("1", "1", "2", "2"),
        val = 1:4
    )

    p <- ggplot(df, aes(x = pos, y = val, domain = chr)) +
        geom_point() +
        coord_serial(domain_gap = 0)

    # Building the plot triggers ggplot_add
    pb <- ggplot_build(p)
    data <- pb$data[[1]]

    # chr 1: min=1, max=10, length=9. offset = 0 - 1 = -1
    # chr 2: min=1, max=5, length=4. gap=0. offset = (0+9+0) - 1 = 8

    # chr 1 positions: 1 + (-1) = 0, 10 + (-1) = 9
    # chr 2 positions: 1 + 8 = 9, 5 + 8 = 13

    expect_equal(data$x, c(0, 9, 9, 13))
})

test_that("coord_serial transforms curve endpoints across domains", {
    df <- data.frame(
        pos = c(1, 10, 1, 5),
        chr = c("1", "1", "2", "2"),
        val = 1:4
    )
    links <- data.frame(
        x = 10,
        domain_x = "1",
        xend = 5,
        domain_end = "2",
        y = 1,
        yend = 2
    )

    p <- ggplot(df, aes(x = pos, y = val, domain = chr)) +
        geom_point() +
        geom_curve(
            data = links,
            aes(x = x, xend = xend, y = y, yend = yend),
            inherit.aes = FALSE
        ) +
        coord_serial(domain_gap = 0)

    expect_equal(p$layers[[2]]$data$x, 9)
    expect_equal(p$layers[[2]]$data$xend, 13)
})

test_that("coord_serial transforms rect xmin/xmax across domains", {
    df <- data.frame(
        pos = c(1, 10, 1, 5),
        chr = c("1", "1", "2", "2"),
        val = 1:4
    )
    regions <- data.frame(
        xmin = c(2, 2),
        xmax = c(8, 4),
        domain_xmin = c("1", "2"),
        domain_xmax = c("1", "2"),
        y = c(0.5, 0.5)
    )

    p <- ggplot(df, aes(x = pos, y = val, domain = chr)) +
        geom_point() +
        geom_rect(
            data = regions,
            aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = y),
            inherit.aes = FALSE,
            alpha = 0.2
        ) +
        coord_serial(domain_gap = 0)

    expect_equal(p$layers[[2]]$data$xmin, c(1, 10))
    expect_equal(p$layers[[2]]$data$xmax, c(7, 12))
})

test_that("coord_serial handles missing domain aesthetic gracefully", {
    df <- data.frame(pos = 1:10, val = 1:10)
    p <- ggplot(df, aes(x = pos, y = val)) +
        geom_point() +
        coord_serial()

    # Should not error, just not transform
    pb <- ggplot_build(p)
    expect_equal(pb$data[[1]]$x, 1:10)
})

test_that("coord_serial handles empty data frame", {
    df <- data.frame(pos = numeric(0), val = numeric(0), chr = character(0))
    p <- ggplot(df, aes(x = pos, y = val, domain = chr)) +
        geom_point() +
        coord_serial()

    expect_no_error(ggplot_build(p))
})

test_that("scaffold works correctly", {
    # 1. compute_serial_layout with scaffold
    scaffold <- c("A" = 100, "B" = 200)
    layout <- compute_serial_layout(numeric(0), character(0), scaffold = scaffold, domain_gap = 0)
    
    expect_equal(layout$domain_levels, c("A", "B"))
    expect_equal(layout$offsets[["A"]], 0)
    expect_equal(layout$offsets[["B"]], 100)
    expect_equal(layout$total_range, c(0, 300))
    
    # 2. Reordering scaffold
    layout_ord <- compute_serial_layout(numeric(0), character(0), 
                                        scaffold = scaffold, 
                                        domain_order = c("B", "A"), 
                                        domain_gap = 0)
    expect_equal(layout_ord$domain_levels, c("B", "A"))
    expect_equal(layout_ord$offsets[["B"]], 0)
    expect_equal(layout_ord$offsets[["A"]], 200)

    # 3. ggplot integration
    df <- data.frame(pos = 10, chr = "A", val = 1)
    p <- ggplot(df, aes(x = pos, y = val, domain = chr)) +
        geom_point() +
        coord_serial(scaffold = scaffold, domain_gap = 0)
    
    pb <- ggplot_build(p)
    expect_equal(pb$data[[1]]$x, 10)
    
    # Total range should reflect scaffold
    expect_equal(pb$plot$coordinates$.serial_layout$total_range, c(0, 300))
})
