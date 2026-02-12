# coord_serial - plot grouped serial data along a continuous axis

<!-- badges: start -->
[![R-CMD-check](https://github.com/daniel-wells/coord_serial/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/daniel-wells/coord_serial/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

This is an R package - an extension of the `ggplot2` package - for plotting data along a continuous axis composed of multiple discrete domains (like chromosomes in a manhattan plot), in a grammar of graphics approach without the need for bespoke data pre-processing.

Unlike other packages like fastman, `coord_serial` does not create a new function to plot a manhattan plot with many customisation arguments. Instead, it provides a new coordinate system, `coord_serial()`, that can be used with `ggplot2` to plot data along the genome and so allows full customisation of the plot using the `ggplot2` framework.

## example

While `coord_serial` is general-purpose, it is often used for genetics data like manhattan plots. Here is a concrete example using simulated GWAS data:

```r
library(coord.serial)
library(ggplot2)

# genetics-specific example data
simulated_gwas <- simulate_gwas()

# plot using the generalized coord_serial()
ggplot(simulated_gwas, aes(x = position, y = log10p, domain = chrom, color = causal)) + \
    geom_point(alpha = 0.5) + \
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey40")) + \
    coord_serial() + \
    geom_hline(yintercept = 7.3, linetype = "dashed", color = "red") + \
    theme_minimal() + \
    theme(plot.background = element_rect(fill = "white", colour = NA)) + \
    labs(x = "Chromosome", color = "Causal Variant")
```

![Manhattan Plot](man/figures/manhattan.png)

## installation

```r
devtools::install_github("daniel-wells/coord_serial")
```
