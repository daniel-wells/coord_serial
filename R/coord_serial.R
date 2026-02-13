#' Serial Coordinate System
#'
#' `coord_serial()` provides a coordinate system for plotting data along a
#' continuous axis composed of multiple discrete domains. It takes
#' per-domain positions and maps them to a continuous serial x-axis with
#' domain labels — similar to a manhattan plot — without requiring manual
#' data preprocessing.
#'
#' The domain grouping is specified via the `domain` aesthetic in `aes()`.
#' Internally, `coord_serial()` transforms each layer's x-positions into
#' cumulative coordinates and relabels the x-axis with domain names.
#'
#' @param domain_gap Fractional gap between domains as a proportion of total
#'   length. Default `0.01` (1%).
#' @param domain_order Vector specifying the order of domains. If `NULL`
#'   (the default), domains are sorted naturally.
#' @param ylim Numeric vector of length 2 giving y-axis limits, or `NULL`
#'   (default) for automatic limits.
#' @param expand If `TRUE` (the default), adds a small expansion factor to the
#'   axes.
#' @param clip Should drawing be clipped to the extent of the plot panel?
#'   `"on"` (the default) or `"off"`.
#'
#' @return A ggplot2 `Coord` object that can be added to a ggplot.
#'
#' @examples
#' library(ggplot2)
#'
#' # Generate genetics-specific data on the fly
#' simulated_gwas <- simulate_gwas()
#'
#' ggplot(simulated_gwas, aes(x = position, y = log10p, domain = chrom)) +
#'   geom_point(size = 0.5) +
#'   coord_serial()
#'
#' @export
coord_serial <- function(domain_gap = 0.01,
                         domain_order = NULL,
                         ylim = NULL,
                         expand = TRUE,
                         clip = "on") {
  ggplot2::ggproto(NULL, CoordSerial,
    domain_gap   = domain_gap,
    domain_order = domain_order,
    limits       = list(x = NULL, y = ylim),
    expand       = expand,
    clip         = clip
  )
}


# ========================================================================
# Internal helpers
# ========================================================================

#' Sort domains in natural order
#'
#' Numeric domains first, then known genomic domains (X, Y, MT), then
#' alphabetical.
#'
#' @param domains Character vector of domain names.
#' @return Character vector of unique domains in sorted order.
#' @keywords internal
sort_domains <- function(domains) {
  uc <- unique(as.character(domains))
  stripped <- sub("^chr", "", uc, ignore.case = TRUE)

  # Numeric domains
  num_val <- suppressWarnings(as.numeric(stripped))
  is_num <- !is.na(num_val)
  num_dom <- uc[is_num][order(num_val[is_num])]

  # Known genomic non-numeric
  known <- c("X", "Y", "MT", "M")
  non_num <- uc[!is_num]
  nn_strip <- stripped[!is_num]
  known_idx <- match(toupper(nn_strip), known)
  is_known <- !is.na(known_idx)
  known_dom <- non_num[is_known][order(known_idx[is_known])]

  # Unknown non-numeric
  unknown_dom <- sort(non_num[!is_known])

  c(num_dom, known_dom, unknown_dom)
}


#' Compute cumulative serial layout
#'
#' Given a combined data frame of positions and domains, compute the
#' cumulative offset for each domain.
#'
#' @param positions Numeric vector of positions.
#' @param domains Vector of domain identifiers (same length as `positions`).
#' @param domain_order Optional vector of domain order.
#' @param domain_gap Gap between domains as a fraction of total length.
#' @return A list with elements `domain_levels`, `offsets`, `midpoints`,
#'   `total_range`.
#' @keywords internal
compute_serial_layout <- function(positions, domains, domain_order = NULL,
                                  domain_gap = 0.01) {
  domains <- as.character(domains)

  # Determine domain order
  if (!is.null(domain_order)) {
    domain_levels <- as.character(domain_order)
  } else {
    domain_levels <- sort_domains(domains)
  }

  domain_f <- factor(domains, levels = domain_levels)

  # Per-domain min and max positions
  domain_min <- tapply(positions, domain_f, min, na.rm = TRUE)
  domain_max <- tapply(positions, domain_f, max, na.rm = TRUE)

  domain_min[is.infinite(domain_min)] <- 0
  domain_max[is.infinite(domain_max)] <- 0

  domain_lengths <- domain_max - domain_min
  total_length <- sum(domain_lengths, na.rm = TRUE)
  gap <- domain_gap * total_length

  # Cumulative offsets
  offsets <- numeric(length(domain_levels))
  midpoints <- numeric(length(domain_levels))
  names(offsets) <- names(midpoints) <- domain_levels

  cumul <- 0
  for (i in seq_along(domain_levels)) {
    dom <- domain_levels[i]
    offsets[dom] <- cumul - domain_min[dom]
    midpoints[dom] <- cumul + domain_lengths[dom] / 2
    cumul <- cumul + domain_lengths[dom] + gap
  }

  total_range <- c(0, cumul - gap)

  list(
    domain_levels = domain_levels,
    offsets       = offsets,
    midpoints     = midpoints,
    total_range   = total_range
  )
}


# ========================================================================
# CoordSerial ggproto class
# ========================================================================

#' @title CoordSerial ggproto object
#' @format NULL
#' @usage NULL
#' @export
CoordSerial <- ggplot2::ggproto("CoordSerial", ggplot2::CoordCartesian,

  # --- fields ---------------------------------------------------------------
  domain_gap = 0.01,
  domain_order = NULL,

  # --- setup_panel_params ---------------------------------------------------
  setup_panel_params = function(self, scale_x, scale_y, params = list()) {
    layout <- self$.serial_layout

    if (!is.null(layout)) {
      scale_x$range$range <- layout$total_range
    }

    panel_params <- ggplot2::ggproto_parent(
      ggplot2::CoordCartesian, self
    )$setup_panel_params(scale_x, scale_y, params)

    if (!is.null(layout) && !is.null(panel_params$x)) {
      panel_params$x$break_positions <- function(...) {
        layout$midpoints
      }
      panel_params$x$get_labels <- function(...) {
        layout$domain_levels
      }

      panel_params$x$breaks <- layout$midpoints
      panel_params$x$minor_breaks <- NULL

      # Title can be overridden by labs(x=...)
      panel_params$x$name <- "Domain"
    }

    panel_params
  },

  # --- stashed state --------------------------------------------------------
  .serial_layout = NULL
)


# ========================================================================
# ggplot_add method
# ========================================================================

#' @importFrom ggplot2 ggplot_add
#' @importFrom rlang as_name
#' @method ggplot_add CoordSerial
#' @export
ggplot_add.CoordSerial <- function(object, plot, ...) {
  plot$coordinates <- object

  # 1. Determine the mapping columns from the main plot mapping
  domain_mapping <- plot$mapping$domain
  x_mapping <- plot$mapping$x

  if (is.null(domain_mapping) || is.null(x_mapping)) {
    return(plot)
  }

  domain_col <- rlang::as_name(domain_mapping)
  x_col <- rlang::as_name(x_mapping)

  # 2. Compute layout using the union of ALL data (plot data + layer data)
  # This ensures all domains from all layers are accounted for in the layout.
  all_x <- plot$data[[x_col]]
  all_dom <- as.character(plot$data[[domain_col]])

  for (layer in plot$layers) {
    l_data <- layer$data %||% plot$data
    if (!is.null(l_data) && domain_col %in% names(l_data) && x_col %in% names(l_data)) {
      all_x <- c(all_x, l_data[[x_col]])
      all_dom <- c(all_dom, as.character(l_data[[domain_col]]))
    }
  }

  if (length(all_x) == 0) {
    return(plot)
  }

  layout <- compute_serial_layout(
    positions    = all_x,
    domains      = all_dom,
    domain_order = object$domain_order,
    domain_gap   = object$domain_gap
  )
  plot$coordinates$.serial_layout <- layout

  # 3. Transform the main plot data
  if (!is.null(plot$data) && domain_col %in% names(plot$data) && x_col %in% names(plot$data)) {
    domain_char <- as.character(plot$data[[domain_col]])
    plot$data[[x_col]] <- plot$data[[x_col]] + layout$offsets[domain_char]
  }

  # 4. Transform ALL layer data
  for (i in seq_along(plot$layers)) {
    # If the layer has its own data, we MUST transform it
    if (!is.null(plot$layers[[i]]$data) && !is.waive(plot$layers[[i]]$data)) {
      l_data <- plot$layers[[i]]$data
      if (domain_col %in% names(l_data) && x_col %in% names(l_data)) {
        domain_char <- as.character(l_data[[domain_col]])
        l_data[[x_col]] <- l_data[[x_col]] + layout$offsets[domain_char]
        plot$layers[[i]]$data <- l_data
      }
    }
  }

  plot
}

# Helper to check for waive
is.waive <- function(x) inherits(x, "waiver")


# Small utility: NULL-default operator
`%||%` <- function(a, b) if (!is.null(a)) a else b
