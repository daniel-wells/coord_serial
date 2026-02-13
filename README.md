# coord_serial - plot grouped serial data along a continuous axis

<!-- badges: start -->
[![R-CMD-check](https://github.com/daniel-wells/coord_serial/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/daniel-wells/coord_serial/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

This is an R package - an extension of the `ggplot2` package - for plotting data along a continuous axis composed of multiple discrete domains (like chromosomes in a manhattan plot), in a grammar of graphics approach:
- no bespoke data pre-processing
- full customisation using the `ggplot2` framework

This is achieved by providing a new coordinate system, `coord_serial()`, that can be used with `ggplot2` to plot data along a single axis.



## Example

While `coord_serial` is general-purpose, it is often used for genetics data like manhattan plots. Here is a concrete example using simulated GWAS data:

```r
library(coord.serial)
library(ggplot2)

# genetics-specific example data
simulated_gwas <- simulate_gwas()

# Create a label for a specific hit
target_hit <- simulated_gwas[simulated_gwas$chrom == "17" & simulated_gwas$causal, ][1, ]
target_hit$label <- "rs9823673"

# plot using the generalized coord_serial()
p <- ggplot(simulated_gwas, aes(x = position,
                                y = log10p,
                                domain = chrom,
                                color = causal)) +
    geom_point() +
    geom_text(data = target_hit, aes(label = label), vjust = -1, size = 5, color = "black") +
    geom_hline(yintercept = 7.3, linetype = "dashed", color = "red") +
    coord_serial()

p + scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey40")) +
    theme_minimal() +
    labs(x = "Chromosome", color = "Causal Variant")
```

![Manhattan Plot](man/figures/manhattan.png)

## genomic domain visualization (Exons)
Often you want to visualize a score or depth across multiple non-contiguous regions like exons. `coord_serial` makes this easy by treating exons as domains and automatically handling the spacing.

```r
# Simulate conservation scores for exons of varying lengths
exon_data <- simulate_serial(
    n = 1000, 
    domains = c("Exon 1", "Exon 2", "Exon 3", "Exon 4"), 
    domain_lengths = c(120, 300, 80, 200)
)

ggplot(exon_data, aes(x = position, y = log10p, domain = domain, fill = domain)) +
    geom_area(show.legend = FALSE) +
    coord_serial() +
    theme_minimal() +
    labs(title = "Conservation Scores across Exons", x = "Position", y = "Score")
```

![Exon Plot](man/figures/exons.png)

## non-genomic time series (Daily Sessions)
`coord_serial` isn't limited to biology. It can be used for any data where you want to stitch together disparate "blocks" of time or space into a single continuous view. Here is an example of monitoring system activity across different daily sessions with widely different durations.

```r
# Simulate activity metrics for sessions with different durations (seconds)
sessions <- simulate_serial(
    n = 2000,
    domains = c("A", "B", "C", "D"),
    domain_lengths = c(3600, 7200, 2400, 300)
)
# Add a random walk per session with randomized starting points
sessions$metric <- ave(rnorm(nrow(sessions)), sessions$domain, FUN = function(x) {
    cumsum(x) + runif(1, -10, 10)
})

ggplot(sessions, aes(x = position, y = metric, domain = domain, color = domain)) +
    geom_line(show.legend = FALSE) +
    coord_serial() +
    theme_light() +
    labs(title = "System Load across Sessions", x = "Time (seconds)", y = "Load Metric")
```

![Session Plot](man/figures/sessions.png)

## Installation

```r
devtools::install_github("daniel-wells/coord_serial")
```
