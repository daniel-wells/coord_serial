#' Simulate GWAS data
#'
#' Generates a simulated GWAS dataset with genomic positions and p-values.
#' This is a concrete genetics-specific example of serial data.
#'
#' @param n Total number of variants to simulate.
#' @param seed Random seed for reproducibility.
#' @param p_max Maximum p-value for background distribution (default 1.0).
#'
#' @return A data frame with columns `chrom`, `position`, `log10p`, and `causal`.
#' @importFrom stats runif
#' @export
simulate_gwas <- function(n = 1000000, seed = 123, p_max = 1.0) {
    # Approximate human chromosome lengths (Mb)
    chrom_lengths_mb <- c(
        249, 243, 198, 191, 181, 171, 159, 145, 138, 134,
        135, 133, 114, 107, 102, 90, 83, 80, 59, 64, 47, 51,
        156
    )
    chrom_names <- c(as.character(1:22), "X")

    # Use the generic engine
    df <- simulate_serial(
        n = n,
        domains = chrom_names,
        domain_lengths = chrom_lengths_mb * 1e6,
        hit_domains = c("3", "6", "11", "17"),
        seed = seed,
        p_max = p_max
    )

    # Standardize column names for GWAS
    names(df)[names(df) == "domain"] <- "chrom"

    # Ensure chrom is a factor with standard order
    df$chrom <- factor(df$chrom, levels = chrom_names)

    df
}
