#' Simulate serial data
#'
#' Generates a simulated serial dataset with positions and values across
#' multiple domains.
#'
#' @param n Total number of points to simulate.
#' @param domains Vector of domain names or an integer for number of domains.
#' @param domain_lengths Optional vector of domain lengths. If `NULL`, all
#'   domains get length 1e8.
#' @param hit_domains Optional vector of indices or names where "hits" (higher
#'   values) should be generated.
#' @param seed Random seed for reproducibility.
#'
#' @return A data frame with columns `domain`, `position`, `log10p`, and `causal`.
#' @importFrom stats runif
#' @export
simulate_serial <- function(n = 1000000,
                            domains = 10,
                            domain_lengths = NULL,
                            hit_domains = c(2, 5, 8),
                            seed = 123) {
    set.seed(seed)

    if (is.numeric(domains) && length(domains) == 1) {
        domain_names <- as.character(seq_len(domains))
        if (is.null(domain_lengths)) {
            domain_lengths <- rep(1e8, domains)
        }
    } else {
        domain_names <- as.character(domains)
        if (is.null(domain_lengths)) {
            domain_lengths <- rep(1e8, length(domains))
        }
    }

    # Proportional point counts
    n_per_dom <- round(n * (domain_lengths / sum(domain_lengths)))

    rows <- list()
    for (i in seq_along(domain_lengths)) {
        dom_name <- domain_names[i]
        nn <- n_per_dom[i]

        pos <- sort(sample.int(domain_lengths[i], nn))
        log10p <- -log10(runif(nn))
        causal <- rep(FALSE, nn)

        # Add hits on select domains
        if (i %in% hit_domains || dom_name %in% hit_domains) {
            n_hits <- sample(3:8, 1)
            hit_idx <- sample(nn, n_hits)
            log10p[hit_idx] <- -log10(10^(-runif(n_hits, 5, 10)))
            causal[hit_idx] <- TRUE
        }

        rows[[i]] <- data.frame(
            domain = factor(dom_name, levels = domain_names),
            position = as.integer(pos),
            log10p = round(log10p, 2),
            causal = causal,
            stringsAsFactors = FALSE
        )
    }

    do.call(rbind, rows)
}
