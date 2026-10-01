community_workspace <- function() {
  override <- Sys.getenv("COMMUNITY_WORKSPACE", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script)) stop("Run with Rscript or set COMMUNITY_WORKSPACE.")
  dirname(normalizePath(sub("^--file=", "", script[[1]]), mustWork = TRUE))
}

# ===========================================================
# figure4_sum_abs_std.R
#
# Figure 4 sum_abs_std variant:
#   shannon / gamma_shannon vs sum_abs_std
#   x-axis zone threshold:
#     mortality   = 0.33
#     temperature = 0.20
#   diversity zone threshold: Shannon = 0.8
#
# Original community_cv Figure 4 is left unchanged in figure4.R.
# ===========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(cowplot)
})

W_PLOT <- 65
H_PLOT <- 60

FONT_FAMILY <- "Arial"
FONT_AX     <- 10
FONT_LEGEND <- 12

BORDER_SIZE <- 0.25
TICK_SIZE   <- 0.20
TICK_LEN    <- -0.4

DOT_SIZE  <- 1.5
DOT_ALPHA <- 1
N_TICKS   <- 3

SUM_ABS_STD_THRESHOLDS <- c(
  mortality   = 0.25,
  temperature = 0.20
)
DIVERSITY_THRESHOLD <- 0.8

ZONE_COLOR_STABLE    <- "#E8E4F0"
ZONE_COLOR_FLUC_LOW  <- "#FDE8D0"
ZONE_COLOR_FLUC_HIGH <- "#FAD4B5"

COL_STABLE      <- "#9B8EC4"
COL_FLUCTUATION <- "#F4A460"

COND_SHAPES <- c(
  "W1" = 16,
  "W2" = 17,
  "W3" = 15,
  "W4" = 18,
  "W5" = 8
)

COND_LABELS <- list(
  mortality = c(
    "W1" = "10\u00B9",
    "W2" = "10\u00B2",
    "W3" = "10\u00B3",
    "W4" = "10\u2074",
    "W5" = "10\u2075"
  ),
  temperature = c(
    "W1" = "10\u00B0C",
    "W2" = "20\u00B0C",
    "W3" = "30\u00B0C",
    "W4" = "40\u00B0C",
    "W5" = "50\u00B0C"
  )
)

WINDOWS_BY_EXPERIMENT <- list(
  mortality   = c("early", "full"),
  temperature = c("full")
)

X_LIMITS <- list(
  mortality   = c(0, 1.6),
  temperature = c(0, 0.8)
)

Y_LO <- 0
Y_HI <- 1.6

BASE_DIR <- community_workspace()
PROC_DIR <- file.path(BASE_DIR, "processed")
FIG_DIR  <- file.path(BASE_DIR, "figures", "diversity_sum_abs_std_zones")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

make_ticks <- function(lo, hi, n = N_TICKS) {
  seq(lo, hi, length.out = n)
}

fmt_axis <- function(x) {
  ifelse(x == floor(x),
         as.character(as.integer(x)),
         gsub("0+$", "", formatC(x, format = "f", digits = 3)))
}

theme_pub <- function() {
  theme_classic(base_size = FONT_AX, base_family = FONT_FAMILY) +
    theme(
      axis.line         = element_blank(),
      panel.border      = element_rect(linewidth = BORDER_SIZE,
                                       color = "black", fill = NA),
      axis.ticks        = element_line(linewidth = TICK_SIZE,
                                       color = "black"),
      axis.ticks.length = unit(TICK_LEN, "mm"),
      axis.text         = element_text(size = FONT_AX, color = "black"),
      axis.title        = element_blank(),
      legend.position   = "none",
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA),
      plot.margin       = margin(3, 3, 3, 3, "mm")
    )
}

make_zone_scatter <- function(df_plot, x_lim, sum_threshold) {
  x_lo <- x_lim[1]
  x_hi <- x_lim[2]
  x_breaks <- make_ticks(x_lo, x_hi)
  y_breaks <- make_ticks(Y_LO, Y_HI)

  zone_rects <- tibble(
    xmin = c(x_lo,          sum_threshold, sum_threshold),
    xmax = c(sum_threshold, x_hi,          x_hi),
    ymin = c(Y_LO,                  Y_LO,                  DIVERSITY_THRESHOLD),
    ymax = c(Y_HI,                  DIVERSITY_THRESHOLD,   Y_HI),
    zone = c("Stable",
             "Fluctuating + Low diversity",
             "Fluctuating + High diversity")
  )

  zone_colors <- c(
    "Stable"                       = ZONE_COLOR_STABLE,
    "Fluctuating + Low diversity"  = ZONE_COLOR_FLUC_LOW,
    "Fluctuating + High diversity" = ZONE_COLOR_FLUC_HIGH
  )

  ggplot(df_plot, aes(x = x_val, y = value)) +
    geom_rect(
      data = zone_rects,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = zone),
      alpha = 0.5,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = zone_colors) +
    geom_point(aes(color = fluctuating, shape = condition),
               size = DOT_SIZE, alpha = DOT_ALPHA) +
    scale_color_manual(values = c("TRUE" = COL_FLUCTUATION,
                                  "FALSE" = COL_STABLE)) +
    scale_shape_manual(values = COND_SHAPES) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = fmt_axis(x_breaks),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = fmt_axis(y_breaks),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(xlim = x_lim,
                    ylim = c(Y_LO, Y_HI),
                    clip = "on") +
    labs(x = NULL, y = NULL) +
    theme_pub()
}

make_stability_legend_plot <- function() {
  legend_df <- tibble(
    stability = factor(c("Stable", "Fluctuating"),
                       levels = c("Stable", "Fluctuating")),
    x = 1:2,
    y = 1
  )

  ggplot(legend_df, aes(x = x, y = y, color = stability)) +
    geom_point(size = 4) +
    scale_color_manual(
      values = c("Stable" = COL_STABLE, "Fluctuating" = COL_FLUCTUATION),
      name = "Community stability"
    ) +
    theme_void(base_family = FONT_FAMILY) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = FONT_LEGEND, family = FONT_FAMILY),
      legend.text  = element_text(size = FONT_LEGEND, family = FONT_FAMILY),
      legend.key.size = unit(5, "mm"),
      legend.margin = margin(2, 2, 2, 2),
      plot.margin = margin(0, 0, 0, 0)
    )
}

make_condition_shape_legend_plot <- function(experiment) {
  legend_df <- tibble(
    condition = factor(paste0("W", 1:5), levels = paste0("W", 1:5)),
    x = 1:5,
    y = 1
  )

  ggplot(legend_df, aes(x = x, y = y, shape = condition)) +
    geom_point(size = 4, color = "black") +
    scale_shape_manual(values = COND_SHAPES,
                       labels = COND_LABELS[[experiment]]) +
    guides(shape = guide_legend(nrow = 2, byrow = TRUE, title = NULL)) +
    theme_void(base_family = FONT_FAMILY) +
    theme(
      legend.position = "bottom",
      legend.text  = element_text(size = FONT_LEGEND, family = FONT_FAMILY),
      legend.key.size = unit(5, "mm"),
      legend.margin = margin(2, 2, 2, 2),
      plot.margin = margin(0, 0, 0, 0)
    )
}

make_zone_legend_plot <- function(sum_threshold) {
  threshold_label <- sprintf("%.2f", sum_threshold)
  zone_levels <- c(
    sprintf("Stable\n(sum SD < %s)", threshold_label),
    sprintf("Fluctuating + Low diversity\n(sum SD >= %s, D < 0.8)",
            threshold_label),
    sprintf("Fluctuating + High diversity\n(sum SD >= %s, D >= 0.8)",
            threshold_label)
  )
  zone_df <- tibble(
    zone = factor(zone_levels, levels = zone_levels),
    x = 1:3,
    y = 1
  )

  zone_colors <- setNames(
    c(ZONE_COLOR_STABLE, ZONE_COLOR_FLUC_LOW, ZONE_COLOR_FLUC_HIGH),
    zone_levels
  )

  ggplot(zone_df, aes(x = x, y = y, fill = zone)) +
    geom_tile(width = 0.9, height = 0.9,
              color = "grey50", linewidth = 0.3) +
    scale_fill_manual(values = zone_colors, name = "Zone") +
    theme_void(base_family = FONT_FAMILY) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = FONT_LEGEND, family = FONT_FAMILY),
      legend.text  = element_text(size = FONT_LEGEND, family = FONT_FAMILY),
      legend.key.size = unit(6, "mm"),
      legend.margin = margin(2, 2, 2, 2),
      plot.margin = margin(0, 0, 0, 0)
    )
}

save_legend <- function(plot_with_legend, filepath, width_mm, height_mm) {
  legend_grob <- cowplot::get_legend(plot_with_legend)
  legend_plot <- cowplot::ggdraw(legend_grob)
  ggsave(filename = filepath, plot = legend_plot,
         width = width_mm, height = height_mm,
         units = "mm", device = cairo_pdf)
}

for (experiment in names(WINDOWS_BY_EXPERIMENT)) {
  sum_threshold <- SUM_ABS_STD_THRESHOLDS[[experiment]]
  x_lim <- X_LIMITS[[experiment]]

  for (window in WINDOWS_BY_EXPERIMENT[[experiment]]) {
    cat(sprintf("\n========== %s - %s ==========\n", experiment, window))

    fluc_path <- file.path(PROC_DIR, "fluctuations", window, "community_level.csv")
    if (!file.exists(fluc_path)) {
      cat(sprintf("  [SKIP] Fluctuation file not found: %s\n", fluc_path))
      next
    }

    alpha_path <- file.path(PROC_DIR, "diversity", window, "alpha_diversity.csv")
    gamma_path <- file.path(PROC_DIR, "diversity", window, "gamma_diversity.csv")
    if (!file.exists(alpha_path) || !file.exists(gamma_path)) {
      cat("  [SKIP] Diversity file not found\n")
      next
    }

    fluc_df <- read_csv(fluc_path, show_col_types = FALSE) %>%
      filter(experiment == !!experiment, collapsed == FALSE) %>%
      select(condition, community, replica, sum_abs_std)

    alpha_div <- read_csv(alpha_path, show_col_types = FALSE) %>%
      filter(experiment == !!experiment, collapsed == FALSE)

    gamma_div <- read_csv(gamma_path, show_col_types = FALSE) %>%
      filter(experiment == !!experiment, collapsed == FALSE)

    build_plot_data <- function(div_df, div_metric) {
      div_df %>%
        select(condition, community, replica, value = all_of(div_metric)) %>%
        inner_join(
          fluc_df %>%
            select(condition, community, replica, x_val = sum_abs_std),
          by = c("condition", "community", "replica")
        ) %>%
        filter(!is.na(value), !is.na(x_val)) %>%
        mutate(
          condition = factor(condition, levels = paste0("W", 1:5)),
          fluctuating = x_val >= sum_threshold
        )
    }

    df_shannon       <- build_plot_data(alpha_div, "shannon")
    df_gamma_shannon <- build_plot_data(gamma_div, "gamma_shannon")

    if (nrow(df_shannon) >= 1) {
      p1 <- make_zone_scatter(df_shannon, x_lim, sum_threshold)
      ggsave(
        filename = file.path(FIG_DIR, sprintf("shannon_vs_sum_abs_std_%s_%s.pdf",
                                              experiment, window)),
        plot = p1, width = W_PLOT, height = H_PLOT,
        units = "mm", device = cairo_pdf
      )
      cat(sprintf("  -> shannon_vs_sum_abs_std_%s_%s.pdf\n", experiment, window))
    }

    if (nrow(df_gamma_shannon) >= 1) {
      p2 <- make_zone_scatter(df_gamma_shannon, x_lim, sum_threshold)
      ggsave(
        filename = file.path(FIG_DIR,
                             sprintf("gamma_shannon_vs_sum_abs_std_%s_%s.pdf",
                                     experiment, window)),
        plot = p2, width = W_PLOT, height = H_PLOT,
        units = "mm", device = cairo_pdf
      )
      cat(sprintf("  -> gamma_shannon_vs_sum_abs_std_%s_%s.pdf\n",
                  experiment, window))
    }
  }

  save_legend(
    make_condition_shape_legend_plot(experiment),
    file.path(FIG_DIR, sprintf("legend_condition_shape_%s.pdf", experiment)),
    width_mm = 80, height_mm = 20
  )

  save_legend(
    make_zone_legend_plot(sum_threshold),
    file.path(FIG_DIR, sprintf("legend_zones_%s.pdf", experiment)),
    width_mm = 110, height_mm = 25
  )
}

save_legend(
  make_stability_legend_plot(),
  file.path(FIG_DIR, "legend_stability.pdf"),
  width_mm = 60, height_mm = 15
)

cat("\nDone.\n")
cat(sprintf("Output directory: %s\n", FIG_DIR))
cat(sprintf("Thresholds: mortality sum_abs_std=%.2f; temperature sum_abs_std=%.2f; diversity=%.1f\n",
            SUM_ABS_STD_THRESHOLDS[["mortality"]],
            SUM_ABS_STD_THRESHOLDS[["temperature"]],
            DIVERSITY_THRESHOLD))
