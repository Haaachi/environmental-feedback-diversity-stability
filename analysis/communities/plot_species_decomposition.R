community_workspace <- function() {
  override <- Sys.getenv("COMMUNITY_WORKSPACE", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script)) stop("Run with Rscript or set COMMUNITY_WORKSPACE.")
  dirname(normalizePath(sub("^--file=", "", script[[1]]), mustWork = TRUE))
}

# ===========================================================
# plot_species_decomposition_core.R
# Original workspace: /home/hachi/Tem_mortality_workspace
#
# Core species decomposition and turnover figures:
#
# 1. synchrony_phi
# 2. community_cv vs total_biomass_cv
# 3. weighted_mean_pairwise_correlation
# 4. relative_covariance_contribution
# 5. cumulative_excess_richness
# 6. cumulative_excess_shannon
# 7. turnover_events
#
# Plot requirements:
# Draw last3 and long windows separately.
# Draw mortality and temperature separately.
# Use square figures.
# Omit legends from the figures.
# Mortality X-axis exponents are positive: 10^1, 10^2, etc.
#   - raw points + mean ± SEM
#
# Inputs:
#   processed/species_decomposition/
#     last3/community_decomposition.csv
#     long/community_decomposition.csv
#
# Outputs:
#   figures/species_decomposition_core/
# ===========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

BASE_DIR <- community_workspace()
DATA_DIR <- file.path(BASE_DIR, "processed", "species_decomposition")
FIG_DIR  <- file.path(BASE_DIR, "figures", "species_decomposition_core")
STAT_DIR <- file.path(BASE_DIR, "processed", "species_decomposition_core_stats")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(STAT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Plot parameters
# ============================================================

FONT_FAMILY <- "Arial"
FONT_AX     <- 10

BORDER_SIZE <- 0.25
TICK_SIZE   <- 0.20
TICK_LEN    <- -1.0

# Square canvas
W_PLOT <- 75
H_PLOT <- 75

POINT_SIZE  <- 1.35
POINT_ALPHA <- 0.65
MEAN_SIZE   <- 2.0

LINE_WIDTH  <- 0.60
ERR_WIDTH   <- 0.12
ERR_LW      <- 0.40

POINT_COLOR <- "#7A7A7A"
LINE_COLOR  <- "#2E6DA4"
REF_COLOR   <- "black"
SIG_COLOR   <- "black"
SIG_TEXT_SIZE <- 3.0
SIG_LW <- 0.25

# Set TRUE to add a light boxplot behind the points.
# FALSE is recommended for main figures showing raw points and mean +/- SEM.
DRAW_BOXPLOT <- FALSE

BOX_FILL  <- "grey90"
BOX_COLOR <- "grey45"

# ============================================================
# X-axis labels
# Use positive dilution-factor exponents for mortality.
# ============================================================

X_LABELS <- list(
  mortality = c(
    "W1" = expression(10^{1}),
    "W2" = expression(10^{2}),
    "W3" = expression(10^{3}),
    "W4" = expression(10^{4}),
    "W5" = expression(10^{5})
  ),
  temperature = c(
    "W1" = "10\u00b0C",
    "W2" = "20\u00b0C",
    "W3" = "30\u00b0C",
    "W4" = "40\u00b0C",
    "W5" = "50\u00b0C"
  )
)

scale_x_condition <- function(experiment) {
  scale_x_continuous(
    breaks   = 1:5,
    labels   = X_LABELS[[experiment]],
    expand   = expansion(mult = c(0.08, 0.08)),
    sec.axis = dup_axis(labels = NULL, name = NULL)
  )
}

# ============================================================
# Plot theme
# ============================================================

theme_pub <- function() {
  theme_classic(base_size = FONT_AX, base_family = FONT_FAMILY) +
    theme(
      axis.line              = element_blank(),
      panel.border           = element_rect(
        linewidth = BORDER_SIZE,
        color = "black",
        fill = NA
      ),
      axis.ticks             = element_line(
        linewidth = TICK_SIZE,
        color = "black"
      ),
      axis.ticks.x.top       = element_line(
        linewidth = TICK_SIZE,
        color = "black"
      ),
      axis.ticks.y.right     = element_line(
        linewidth = TICK_SIZE,
        color = "black"
      ),
      axis.ticks.length               = unit(TICK_LEN, "mm"),
      axis.ticks.length.x.top         = unit(TICK_LEN, "mm"),
      axis.ticks.length.y.right       = unit(TICK_LEN, "mm"),
      axis.text.x.top        = element_blank(),
      axis.text.y.right      = element_blank(),
      axis.text.x  = element_text(
        size   = FONT_AX,
        color  = "black",
        margin = margin(t = 3)
      ),
      axis.text.y  = element_text(
        size   = FONT_AX,
        color  = "black",
        margin = margin(r = 3)
      ),
      axis.title        = element_blank(),
      legend.position   = "none",
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA),
      plot.margin       = margin(3, 3, 3, 3, "mm")
    )
}

sem <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(NA_real_)
  sd(x) / sqrt(length(x))
}

p_to_sig_label <- function(p) {
  case_when(
    is.na(p)  ~ NA_character_,
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

safe_test_p <- function(expr) {
  tryCatch(expr, error = function(e) NA_real_)
}

add_x_pos <- function(df) {
  df %>%
    mutate(
      condition = as.character(condition),
      x_pos = as.integer(sub("W", "", condition))
    )
}

# ============================================================
# Y-axis limit helpers
# ============================================================

get_line_ylim <- function(values,
                          ymin_values = NULL,
                          ymax_values = NULL,
                          force_zero = TRUE,
                          force_nonnegative = FALSE,
                          symmetric = FALSE,
                          fixed_upper_one = FALSE) {
  
  vals <- c(values, ymin_values, ymax_values)
  vals <- vals[is.finite(vals)]
  
  if (length(vals) == 0) return(c(0, 1))
  
  if (symmetric) {
    m <- max(abs(vals), na.rm = TRUE)
    if (!is.finite(m) || m == 0) m <- 1
    upper <- ceiling(m * 10) / 10
    return(c(-upper, upper))
  }
  
  ymin <- min(vals, na.rm = TRUE)
  ymax <- max(vals, na.rm = TRUE)
  
  if (force_zero) {
    ymin <- min(0, ymin)
    ymax <- max(0, ymax)
  }
  
  if (force_nonnegative) {
    ymin <- 0
  }
  
  if (fixed_upper_one) {
    ymax <- max(1, ymax)
  }
  
  if (!is.finite(ymin) || !is.finite(ymax) || ymin == ymax) {
    ymin <- 0
    ymax <- 1
  }
  
  y_range <- ymax - ymin
  pad <- 0.06 * y_range
  
  if (force_nonnegative) {
    ymin <- 0
    ymax <- ymax + pad
  } else {
    ymin <- ymin - pad
    ymax <- ymax + pad
  }
  
  if (fixed_upper_one) {
    ymax <- max(1, ymax)
  }
  
  c(ymin, ymax)
}

coord_square_condition <- function(ylim) {
  x_range <- 5.4 - 0.6
  y_range <- ylim[2] - ylim[1]
  if (!is.finite(y_range) || y_range <= 0) {
    return(coord_cartesian(clip = "off"))
  }
  coord_fixed(ratio = x_range / y_range, clip = "off")
}

# ============================================================
# Read data.
# ============================================================

read_window <- function(window, filename = "community_decomposition.csv") {
  path <- file.path(DATA_DIR, window, filename)
  if (!file.exists(path)) {
    warning(sprintf("File not found: %s", path))
    return(NULL)
  }
  
  read_csv(path, show_col_types = FALSE) %>%
    mutate(window = window) %>%
    add_x_pos()
}

df_all <- bind_rows(
  read_window("last3"),
  read_window("long")
)

df_fluctuating <- bind_rows(
  read_window("last3", "community_decomposition_fluctuating_cv.csv"),
  read_window("long", "community_decomposition_fluctuating_cv.csv")
)

if (
  is.null(df_fluctuating) ||
  nrow(df_fluctuating) == 0 ||
  !"synchrony_phi" %in% names(df_fluctuating)
) {
  df_fluctuating <- tibble()
} else {
  if (!"synchrony_phi_robust" %in% names(df_fluctuating)) {
    df_fluctuating$synchrony_phi_robust <- NA_real_
  }
  df_fluctuating <- df_fluctuating %>%
    mutate(
      synchrony_phi_fluctuating_cv = synchrony_phi,
      synchrony_phi_robust_fluctuating_cv = synchrony_phi_robust,
      weighted_mean_pairwise_correlation_fluctuating_cv =
        weighted_mean_pairwise_correlation,
      relative_covariance_contribution_fluctuating_cv =
        relative_covariance_contribution,
      cumulative_excess_richness_fluctuating_cv =
        cumulative_excess_richness,
      cumulative_excess_shannon_fluctuating_cv =
        cumulative_excess_shannon,
      turnover_events_fluctuating_cv =
        turnover_events
    )
}

if (is.null(df_all) || nrow(df_all) == 0) {
  stop("No community_decomposition.csv found.")
}

cat(sprintf("community rows: %d\n", nrow(df_all)))
cat(sprintf("fluctuating rows: %d\n", nrow(df_fluctuating)))

STAT_ROWS <- list()

record_temperature_w3_w4_test <- function(d, metric, window, experiment) {
  if (experiment != "temperature") return(tibble())
  
  test_df <- d %>%
    filter(condition %in% c("W3", "W4")) %>%
    transmute(condition, value = .data[[metric]]) %>%
    filter(is.finite(value))
  
  vals_30 <- test_df %>%
    filter(condition == "W3") %>%
    pull(value)
  vals_40 <- test_df %>%
    filter(condition == "W4") %>%
    pull(value)
  
  if (length(vals_30) < 1 || length(vals_40) < 1) {
    return(tibble())
  }
  
  wilcox_p <- safe_test_p(
    wilcox.test(
      vals_30,
      vals_40,
      alternative = "less",
      exact = FALSE
    )$p.value
  )
  welch_p <- if (length(vals_30) >= 2 && length(vals_40) >= 2) {
    safe_test_p(
      t.test(
        vals_30,
        vals_40,
        alternative = "less"
      )$p.value
    )
  } else {
    NA_real_
  }
  
  out <- tibble(
    window = window,
    experiment = experiment,
    metric = metric,
    comparison = "30C_lt_40C",
    alternative = "30C < 40C",
    condition_30C = "W3",
    condition_40C = "W4",
    n_30C = length(vals_30),
    n_40C = length(vals_40),
    mean_30C = mean(vals_30, na.rm = TRUE),
    mean_40C = mean(vals_40, na.rm = TRUE),
    median_30C = median(vals_30, na.rm = TRUE),
    median_40C = median(vals_40, na.rm = TRUE),
    delta_mean_30C_minus_40C = mean_30C - mean_40C,
    wilcox_p_30C_less_40C = wilcox_p,
    welch_t_p_30C_less_40C = welch_p,
    sig_label = p_to_sig_label(wilcox_p)
  )
  
  STAT_ROWS[[length(STAT_ROWS) + 1]] <<- out
  out
}

# ============================================================
# Shared scatter plot with a mean +/- SEM trend
# ============================================================

plot_metric_mean_points <- function(df,
                                    metric,
                                    window,
                                    experiment,
                                    force_nonnegative = FALSE,
                                    symmetric = FALSE,
                                    fixed_upper_one = FALSE,
                                    add_zero_line = FALSE) {
  
  if (!metric %in% names(df)) {
    message(sprintf("[skip] metric not found: %s", metric))
    return(NULL)
  }
  
  d <- df %>%
    filter(
      window == !!window,
      experiment == !!experiment,
      is.finite(.data[[metric]])
    ) %>%
    arrange(condition, community, replica)
  
  if (nrow(d) == 0) {
    message(sprintf("[skip] %s / %s / %s: no data",
                    metric, window, experiment))
    return(NULL)
  }
  
  summary_df <- d %>%
    group_by(condition, x_pos) %>%
    summarise(
      n      = n(),
      mean_v = mean(.data[[metric]], na.rm = TRUE),
      sem_v  = sem(.data[[metric]]),
      .groups = "drop"
    ) %>%
    mutate(
      ymin = mean_v - sem_v,
      ymax = mean_v + sem_v
    ) %>%
    arrange(x_pos)
  
  if (force_nonnegative) {
    summary_df <- summary_df %>%
      mutate(ymin = pmax(0, ymin))
  }
  
  ylim <- get_line_ylim(
    values = d[[metric]],
    ymin_values = summary_df$ymin,
    ymax_values = summary_df$ymax,
    force_zero = TRUE,
    force_nonnegative = force_nonnegative,
    symmetric = symmetric,
    fixed_upper_one = fixed_upper_one
  )
  
  sig_df <- record_temperature_w3_w4_test(d, metric, window, experiment)
  if (nrow(sig_df) > 0 && is.finite(ylim[1]) && is.finite(ylim[2])) {
    y_range <- ylim[2] - ylim[1]
    if (!is.finite(y_range) || y_range <= 0) y_range <- 1
    sig_df <- sig_df %>%
      mutate(
        x1 = 3,
        x2 = 4,
        x_mid = 3.5,
        y = ylim[2] - 0.08 * y_range,
        y_tick = y - 0.04 * y_range,
        y_label = y + 0.015 * y_range
      )
  }
  
  y_breaks <- pretty(ylim, n = 5)
  y_breaks <- y_breaks[y_breaks >= ylim[1] & y_breaks <= ylim[2]]
  
  p <- ggplot()
  
  if (add_zero_line) {
    p <- p +
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        linewidth = 0.35,
        color = REF_COLOR
      )
  }
  
  if (DRAW_BOXPLOT) {
    p <- p +
      geom_boxplot(
        data = d,
        aes(x = x_pos, y = .data[[metric]], group = x_pos),
        width = 0.36,
        fill = BOX_FILL,
        color = BOX_COLOR,
        linewidth = 0.25,
        outlier.shape = NA
      )
  }
  
  p <- p +
    geom_jitter(
      data = d,
      aes(x = x_pos, y = .data[[metric]]),
      width  = 0.08,
      height = 0,
      size   = POINT_SIZE,
      alpha  = POINT_ALPHA,
      color  = POINT_COLOR
    ) +
    geom_line(
      data = summary_df,
      aes(x = x_pos, y = mean_v),
      color     = LINE_COLOR,
      linewidth = LINE_WIDTH
    ) +
    geom_errorbar(
      data = summary_df,
      aes(x = x_pos, ymin = ymin, ymax = ymax),
      width     = ERR_WIDTH,
      linewidth = ERR_LW,
      color     = LINE_COLOR,
      na.rm     = TRUE
    ) +
    geom_point(
      data = summary_df,
      aes(x = x_pos, y = mean_v),
      shape = 16,
      size  = MEAN_SIZE,
      color = LINE_COLOR
    ) +
    scale_x_condition(experiment) +
    scale_y_continuous(
      limits   = ylim,
      breaks   = y_breaks,
      expand   = expansion(mult = c(0, 0)),
      sec.axis = dup_axis(labels = NULL, name = NULL)
    ) +
    labs(x = NULL, y = NULL) +
    coord_square_condition(ylim) +
    theme_pub()
  
  if (nrow(sig_df) > 0 && !is.na(sig_df$sig_label[1])) {
    p <- p +
      geom_segment(
        data = sig_df,
        aes(x = x1, xend = x2, y = y, yend = y),
        inherit.aes = FALSE,
        linewidth = SIG_LW,
        color = SIG_COLOR
      ) +
      geom_segment(
        data = sig_df,
        aes(x = x1, xend = x1, y = y_tick, yend = y),
        inherit.aes = FALSE,
        linewidth = SIG_LW,
        color = SIG_COLOR
      ) +
      geom_segment(
        data = sig_df,
        aes(x = x2, xend = x2, y = y_tick, yend = y),
        inherit.aes = FALSE,
        linewidth = SIG_LW,
        color = SIG_COLOR
      ) +
      geom_text(
        data = sig_df,
        aes(x = x_mid, y = y_label, label = sig_label),
        inherit.aes = FALSE,
        size = SIG_TEXT_SIZE,
        family = FONT_FAMILY,
        color = SIG_COLOR,
        vjust = 0
      )
  }
  
  out_dir <- file.path(FIG_DIR, metric)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  out_path <- file.path(
    out_dir,
    sprintf("%s_%s_%s.pdf", metric, window, experiment)
  )
  
  ggsave(
    filename = out_path,
    plot     = p,
    width    = W_PLOT,
    height   = H_PLOT,
    units    = "mm",
    device   = cairo_pdf
  )
  
  cat(sprintf("  -> %s\n", out_path))
}

# ============================================================
# Figure 2: community_cv versus total_biomass_cv
# ============================================================

plot_cv_relationship <- function(df,
                                 window,
                                 experiment,
                                 output_name = "cv_relationship") {
  
  required <- c("community_cv", "total_biomass_cv")
  if (!all(required %in% names(df))) {
    message("[skip] community_cv or total_biomass_cv not found")
    return(NULL)
  }
  
  d <- df %>%
    filter(
      window == !!window,
      experiment == !!experiment,
      is.finite(community_cv),
      is.finite(total_biomass_cv)
    ) %>%
    arrange(condition, community, replica)
  
  if (nrow(d) == 0) {
    message(sprintf("[skip] cv_relationship / %s / %s: no data",
                    window, experiment))
    return(NULL)
  }
  
  summary_df <- d %>%
    group_by(condition, x_pos) %>%
    summarise(
      n = n(),
      mean_community_cv = mean(community_cv, na.rm = TRUE),
      sem_community_cv  = sem(community_cv),
      mean_total_biomass_cv = mean(total_biomass_cv, na.rm = TRUE),
      sem_total_biomass_cv  = sem(total_biomass_cv),
      .groups = "drop"
    ) %>%
    arrange(x_pos) %>%
    mutate(
      xmin = pmax(0, mean_community_cv - sem_community_cv),
      xmax = mean_community_cv + sem_community_cv,
      ymin = pmax(0, mean_total_biomass_cv - sem_total_biomass_cv),
      ymax = mean_total_biomass_cv + sem_total_biomass_cv
    )
  
  max_xy <- max(
    c(
      d$community_cv,
      d$total_biomass_cv,
      summary_df$xmax,
      summary_df$ymax
    ),
    na.rm = TRUE
  )
  
  max_xy <- ceiling(max_xy * 10) / 10
  max_xy <- max(max_xy, 0.1)
  
  p <- ggplot() +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      linewidth = 0.35,
      color = REF_COLOR
    ) +
    geom_point(
      data = d,
      aes(x = community_cv, y = total_biomass_cv),
      size  = POINT_SIZE,
      alpha = POINT_ALPHA,
      color = POINT_COLOR
    ) +
    geom_errorbar(
      data = summary_df,
      aes(
        x = mean_community_cv,
        ymin = ymin,
        ymax = ymax
      ),
      width = 0,
      linewidth = ERR_LW,
      color = LINE_COLOR,
      na.rm = TRUE
    ) +
    geom_segment(
      data = summary_df,
      aes(
        x = xmin,
        xend = xmax,
        y = mean_total_biomass_cv,
        yend = mean_total_biomass_cv
      ),
      linewidth = ERR_LW,
      color = LINE_COLOR,
      na.rm = TRUE
    ) +
    geom_line(
      data = summary_df,
      aes(x = mean_community_cv, y = mean_total_biomass_cv),
      color = LINE_COLOR,
      linewidth = LINE_WIDTH
    ) +
    geom_point(
      data = summary_df,
      aes(x = mean_community_cv, y = mean_total_biomass_cv),
      shape = 16,
      size  = MEAN_SIZE,
      color = LINE_COLOR
    ) +
    scale_x_continuous(
      limits = c(0, max_xy),
      breaks = pretty(c(0, max_xy), n = 5),
      expand = expansion(mult = c(0, 0)),
      sec.axis = dup_axis(labels = NULL, name = NULL)
    ) +
    scale_y_continuous(
      limits = c(0, max_xy),
      breaks = pretty(c(0, max_xy), n = 5),
      expand = expansion(mult = c(0, 0)),
      sec.axis = dup_axis(labels = NULL, name = NULL)
    ) +
    labs(x = NULL, y = NULL) +
    coord_equal(clip = "off") +
    theme_pub()
  
  out_dir <- file.path(FIG_DIR, output_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  out_path <- file.path(
    out_dir,
    sprintf("%s_%s_%s.pdf", output_name, window, experiment)
  )
  
  ggsave(
    filename = out_path,
    plot     = p,
    width    = W_PLOT,
    height   = H_PLOT,
    units    = "mm",
    device   = cairo_pdf
  )
  
  cat(sprintf("  -> %s\n", out_path))
}

# ============================================================
# Main loop for core figures
# ============================================================

cat("\nPlotting core species decomposition figures...\n")

for (window in c("last3", "long")) {
  for (experiment in c("mortality", "temperature")) {
    
    cat(sprintf("\n--- %s / %s ---\n", window, experiment))
    
    # 1. synchrony_phi
    plot_metric_mean_points(
      df = df_all,
      metric = "synchrony_phi",
      window = window,
      experiment = experiment,
      force_nonnegative = TRUE,
      fixed_upper_one = TRUE,
      add_zero_line = FALSE
    )

    plot_metric_mean_points(
      df = df_fluctuating,
      metric = "synchrony_phi_fluctuating_cv",
      window = window,
      experiment = experiment,
      force_nonnegative = TRUE,
      fixed_upper_one = TRUE,
      add_zero_line = FALSE
    )

    plot_metric_mean_points(
      df = df_all,
      metric = "synchrony_phi_robust",
      window = window,
      experiment = experiment,
      force_nonnegative = TRUE,
      fixed_upper_one = TRUE,
      add_zero_line = FALSE
    )

    plot_metric_mean_points(
      df = df_fluctuating,
      metric = "synchrony_phi_robust_fluctuating_cv",
      window = window,
      experiment = experiment,
      force_nonnegative = TRUE,
      fixed_upper_one = TRUE,
      add_zero_line = FALSE
    )
    
    # 2. community_cv vs total_biomass_cv
    plot_cv_relationship(
      df = df_all,
      window = window,
      experiment = experiment
    )

    plot_cv_relationship(
      df = df_fluctuating,
      window = window,
      experiment = experiment,
      output_name = "cv_relationship_fluctuating_cv"
    )
    
    # 3. weighted_mean_pairwise_correlation
    plot_metric_mean_points(
      df = df_all,
      metric = "weighted_mean_pairwise_correlation",
      window = window,
      experiment = experiment,
      force_nonnegative = FALSE,
      symmetric = TRUE,
      add_zero_line = TRUE
    )

    plot_metric_mean_points(
      df = df_fluctuating,
      metric = "weighted_mean_pairwise_correlation_fluctuating_cv",
      window = window,
      experiment = experiment,
      force_nonnegative = FALSE,
      symmetric = TRUE,
      add_zero_line = TRUE
    )
    
    # 4. relative_covariance_contribution
    plot_metric_mean_points(
      df = df_all,
      metric = "relative_covariance_contribution",
      window = window,
      experiment = experiment,
      force_nonnegative = FALSE,
      symmetric = TRUE,
      add_zero_line = TRUE
    )

    plot_metric_mean_points(
      df = df_fluctuating,
      metric = "relative_covariance_contribution_fluctuating_cv",
      window = window,
      experiment = experiment,
      force_nonnegative = FALSE,
      symmetric = TRUE,
      add_zero_line = TRUE
    )
    
    # 5. cumulative_excess_richness
    plot_metric_mean_points(
      df = df_all,
      metric = "cumulative_excess_richness",
      window = window,
      experiment = experiment,
      force_nonnegative = TRUE,
      fixed_upper_one = FALSE,
      add_zero_line = FALSE
    )

    plot_metric_mean_points(
      df = df_fluctuating,
      metric = "cumulative_excess_richness_fluctuating_cv",
      window = window,
      experiment = experiment,
      force_nonnegative = TRUE,
      fixed_upper_one = FALSE,
      add_zero_line = FALSE
    )
    
    # 6. cumulative_excess_shannon
    plot_metric_mean_points(
      df = df_all,
      metric = "cumulative_excess_shannon",
      window = window,
      experiment = experiment,
      force_nonnegative = TRUE,
      fixed_upper_one = FALSE,
      add_zero_line = FALSE
    )

    plot_metric_mean_points(
      df = df_fluctuating,
      metric = "cumulative_excess_shannon_fluctuating_cv",
      window = window,
      experiment = experiment,
      force_nonnegative = TRUE,
      fixed_upper_one = FALSE,
      add_zero_line = FALSE
    )
    
    # 7. turnover_events
    plot_metric_mean_points(
      df = df_all,
      metric = "turnover_events",
      window = window,
      experiment = experiment,
      force_nonnegative = TRUE,
      fixed_upper_one = FALSE,
      add_zero_line = FALSE
    )

    plot_metric_mean_points(
      df = df_fluctuating,
      metric = "turnover_events_fluctuating_cv",
      window = window,
      experiment = experiment,
      force_nonnegative = TRUE,
      fixed_upper_one = FALSE,
      add_zero_line = FALSE
    )
  }
}

# ============================================================
# Output index
# ============================================================

index_df <- tibble(
  figure_type = c(
    "synchrony_phi",
    "synchrony_phi_fluctuating_cv",
    "synchrony_phi_robust",
    "synchrony_phi_robust_fluctuating_cv",
    "cv_relationship",
    "cv_relationship_fluctuating_cv",
    "weighted_mean_pairwise_correlation",
    "weighted_mean_pairwise_correlation_fluctuating_cv",
    "relative_covariance_contribution",
    "relative_covariance_contribution_fluctuating_cv",
    "cumulative_excess_richness",
    "cumulative_excess_richness_fluctuating_cv",
    "cumulative_excess_shannon",
    "cumulative_excess_shannon_fluctuating_cv",
    "turnover_events",
    "turnover_events_fluctuating_cv"
  ),
  interpretation = c(
    "Final species synchrony index",
    "Species synchrony index only for communities with community_cv >= 0.25 or temporal BC >= 0.15",
    "Robust species synchrony index after excluding near-zero phi denominators",
    "Robust species synchrony index only for fluctuating communities",
    "Geometric relation between species-level fluctuation and total biomass fluctuation",
    "CV relationship only for communities with community_cv >= 0.25 or temporal BC >= 0.15",
    "Weighted mean species-pair correlation",
    "Weighted mean species-pair correlation only for communities with community_cv >= 0.25 or temporal BC >= 0.15",
    "Relative contribution of covariance to biomass variance",
    "Relative covariance contribution only for communities with community_cv >= 0.25 or temporal BC >= 0.15",
    "Cumulative richness minus mean daily richness, based on rel_abund >= 1%",
    "Cumulative richness excess only for communities with community_cv >= 0.25 or temporal BC >= 0.15",
    "Cumulative Shannon minus mean daily Shannon, based on summed absolute biomass across the window",
    "Cumulative Shannon excess only for communities with community_cv >= 0.25 or temporal BC >= 0.15",
    "Total number of within-window presence/absence turnover events, based on rel_abund >= 1%",
    "Turnover events only for communities with community_cv >= 0.25 or temporal BC >= 0.15"
  )
)

write_csv(index_df, file.path(FIG_DIR, "figure_index.csv"))

if (length(STAT_ROWS) > 0) {
  stat_out <- bind_rows(STAT_ROWS) %>%
    arrange(window, metric)
  stat_path <- file.path(STAT_DIR, "temperature_30C_lt_40C_tests.csv")
  write_csv(stat_out, stat_path)
  cat(sprintf("  -> stats: %s\n", stat_path))
}

cat("\nDone!\n")
cat(sprintf("Output directory: %s\n", FIG_DIR))
cat("\nOutput subdirectories:\n")
cat("  synchrony_phi/\n")
cat("  cv_relationship/\n")
cat("  cv_relationship_fluctuating_cv/\n")
cat("  weighted_mean_pairwise_correlation/\n")
cat("  weighted_mean_pairwise_correlation_fluctuating_cv/\n")
cat("  relative_covariance_contribution/\n")
cat("  relative_covariance_contribution_fluctuating_cv/\n")
cat("  cumulative_excess_richness/\n")
cat("  cumulative_excess_richness_fluctuating_cv/\n")
cat("  cumulative_excess_shannon/\n")
cat("  cumulative_excess_shannon_fluctuating_cv/\n")
cat("  turnover_events/\n")
cat("  turnover_events_fluctuating_cv/\n")
cat("\nTo inspect distribution shapes, set DRAW_BOXPLOT <- TRUE at the top of the script.\n")
