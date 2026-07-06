# Shared helpers for the laser-generic R/reticulate tutorials.
# Source from `.Rmd` chunks: source("helpers.R")
#
# Convention: every helper takes plain R values and handles the
# R <-> Python dtype boundary inside, so the tutorial chunks read like
# ordinary R code.

# ---- module references ---------------------------------------------------

# These are filled by `setup_python()` and reused by everything else. We
# stash them in a per-session environment to avoid repeated `import()` calls.
.lg_env <- new.env(parent = emptyenv())

#' Lazy-import the laser modules we need.
#'
#' Called once at the top of each `.Rmd`. Safe to call repeatedly.
setup_python <- function() {
    if (!is.null(.lg_env$np)) return(invisible(NULL))
    .lg_env$np            <- reticulate::import("numpy",                            convert = FALSE)
    .lg_env$laser_generic <- reticulate::import("laser.generic",                    convert = FALSE)
    .lg_env$laser_core    <- reticulate::import("laser.core",                       convert = FALSE)
    .lg_env$laser_utils   <- reticulate::import("laser.core.utils",                 convert = FALSE)
    .lg_env$laser_random  <- reticulate::import("laser.core.random",                convert = FALSE)
    .lg_env$laser_dist    <- reticulate::import("laser.core.distributions",         convert = FALSE)
    .lg_env$laser_demog   <- reticulate::import("laser.core.demographics",          convert = FALSE)
    .lg_env$laser_propset <- reticulate::import("laser.core.propertyset",           convert = FALSE)
    .lg_env$laser_vd      <- reticulate::import("laser.generic.vitaldynamics",      convert = FALSE)
    .lg_env$laser_gutils  <- reticulate::import("laser.generic.utils",              convert = FALSE)
    invisible(NULL)
}

# Convenience accessors so chunks can write `lg$Model` etc.
lg            <- function() .lg_env$laser_generic
np            <- function() .lg_env$np
PropertySet   <- function(...) .lg_env$laser_propset$PropertySet(...)
ValuesMap     <- function() .lg_env$laser_gutils$ValuesMap

# ---- seeding -------------------------------------------------------------

#' Seed both R and laser's PRNG for reproducible knits.
seed_everything <- function(seed = 271828L) {
    set.seed(seed)
    .lg_env$laser_random$seed(as.integer(seed)) # this also sets NumPy's random seed
    invisible(seed)
}

# ---- scenario factory ----------------------------------------------------

#' Build a single-node (or grid) scenario.
#'
#' Wraps `laser.core.utils.grid()`. For a single-node scenario, all you need
#' is the population. For a grid, pass `M > 1` and/or `N > 1` and either a
#' scalar `pop` (replicated across nodes) or a vector of length M * N.
#'
#' Returns a GeoDataFrame (proxied as a Python object). Set initial
#' compartment sizes via `scenario[["S"]] <- ...`, etc.
make_scenario <- function(pop,
                          M = 1L, N = 1L,
                          longitude = -122.33,
                          latitude  =   47.60) {
    # Allow either scalar or vector populations.
    if (length(pop) == 1L) {
        # `grid()` expects a Python callable; an R closure auto-converts.
        pop_fn <- function(row, col) as.integer(pop)
    } else {
        stopifnot(length(pop) == M * N)
        pops <- as.integer(pop)
        idx_env <- new.env(parent = emptyenv())
        idx_env$i <- 0L
        pop_fn <- function(row, col) {
            v <- pops[idx_env$i + 1L]
            idx_env$i <- idx_env$i + 1L
            as.integer(v)
        }
    }

    .lg_env$laser_utils$grid(
        M = as.integer(M),
        N = as.integer(N),
        population_fn = pop_fn,
        origin_x = longitude,
        origin_y = latitude
    )
}

# ---- compartment extraction ---------------------------------------------

#' Pull per-tick compartment counts out of a model into a tidy data.frame.
#'
#' `model$nodes$<COMPARTMENT>` is a 2D numpy array of shape (nticks + 1, n_nodes).
#' This function flattens it across nodes (summing per tick) so the resulting
#' data.frame has one row per (tick, compartment) — handy for ggplot.
compartments_df <- function(model, states) {
    rows <- lapply(states, function(state) {
        arr <- reticulate::py_to_r(reticulate::py_get_attr(model$nodes, state))
        # arr is an R matrix: rows = ticks, cols = nodes. Total population
        # of this compartment per tick = rowSums (1 column for single-node).
        if (is.matrix(arr)) totals <- rowSums(arr) else totals <- as.numeric(arr)
        data.frame(
            tick        = seq_len(length(totals)) - 1L,
            count       = totals,
            compartment = state,
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, rows)
}

# ---- plotting ------------------------------------------------------------

#' Plot compartment trajectories with sensible defaults.
plot_compartments <- function(df, title = NULL, subtitle = NULL) {
    p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = tick, y = count, color = compartment)
    ) +
        ggplot2::geom_line(linewidth = 0.7) +
        ggplot2::labs(
            x = "Day", y = "Agents",
            color = NULL,
            title = title, subtitle = subtitle
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(legend.position = "top")
    p
}

# ---- spatial helpers -----------------------------------------------------

#' Per-tick, per-node compartment counts.
#'
#' Unlike `compartments_df()` (which sums across nodes), this returns one
#' row per (tick, node, compartment) so faceted / heat-map plots work
#' naturally.
compartments_per_node_df <- function(model, states) {
    rows <- lapply(states, function(state) {
        arr <- reticulate::py_to_r(reticulate::py_get_attr(model$nodes, state))
        # arr: rows = ticks, cols = nodes (single-node arrays come through
        # as a vector; promote to a 1-column matrix so the loop below is uniform).
        if (!is.matrix(arr)) arr <- matrix(arr, ncol = 1L)
        n_ticks <- nrow(arr)
        n_nodes <- ncol(arr)
        do.call(rbind, lapply(seq_len(n_nodes), function(j) {
            data.frame(
                tick        = seq_len(n_ticks) - 1L,
                node        = j - 1L,                # zero-indexed to match Python
                count       = arr[, j],
                compartment = state,
                stringsAsFactors = FALSE
            )
        }))
    })
    do.call(rbind, rows)
}

#' Facet-by-node line plot of per-node compartment trajectories.
plot_compartments_by_node <- function(df, title = NULL, subtitle = NULL,
                                      ncol = NULL) {
    p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = tick, y = count, color = compartment)
    ) +
        ggplot2::geom_line(linewidth = 0.6) +
        ggplot2::facet_wrap(~ node, ncol = ncol, labeller = ggplot2::label_both) +
        ggplot2::labs(
            x = "Day", y = "Agents", color = NULL,
            title = title, subtitle = subtitle
        ) +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(legend.position = "top")
    p
}

#' Heat-map of a 2D matrix (nodes x ticks, or any other layout).
plot_matrix_heatmap <- function(mat, title = NULL, xlab = "column", ylab = "row",
                                fill_label = "value") {
    df <- expand.grid(row = seq_len(nrow(mat)), col = seq_len(ncol(mat)))
    df$value <- as.vector(mat)
    ggplot2::ggplot(df, ggplot2::aes(col, row, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_y_reverse() +
        ggplot2::scale_fill_viridis_c(option = "viridis") +
        ggplot2::labs(title = title, x = xlab, y = ylab, fill = fill_label) +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(panel.grid = ggplot2::element_blank())
}
