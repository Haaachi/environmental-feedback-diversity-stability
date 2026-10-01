community_workspace <- function() {
  override <- Sys.getenv("COMMUNITY_WORKSPACE", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script)) stop("Run with Rscript or set COMMUNITY_WORKSPACE.")
  dirname(normalizePath(sub("^--file=", "", script[[1]]), mustWork = TRUE))
}

# ===========================================================
# plot_community_cv_sum_std_dual_axis.R
#
# From compute_fluctuations.py outputs:
#   processed/fluctuations/<window>/community_level.csv
#
# Draw two mean curves in one panel, following plot_fluctuation_proportion.R:
#   left y-axis  = community_CV
#   right y-axis = sum of stdNi (sum_abs_std)
#
# No raw jitter points. Mean curves keep the same line, point, errorbar,
# tick, border, and canvas logic as plot_fluctuation_proportion.R.
# ===========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(cowplot)
})

BASE_DIR <- community_workspace()
PROC_DIR <- file.path(BASE_DIR, "processed")
FLUC_DIR <- file.path(PROC_DIR, "fluctuations")
FIG_DIR  <- file.path(BASE_DIR, "figures", "community_cv_sum_std_dual_axis")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

SUMMARY_OUT <- file.path(
  FLUC_DIR,
  "community_cv_sum_std_dual_axis_summary.csv"
)

# ============================================================
# Plot parameters, matching the existing R figure style
# ============================================================

FONT_FAMILY <- "Arial"
FONT_AX     <- 10
FONT_LEGEND <- 12
BORDER_SIZE <- 0.25
TICK_SIZE   <- 0.20
TICK_LEN    <- -1.0

W_PLOT <- 80
H_PLOT <- 70
W_LEGEND <- 80
H_LEGEND <- 18

COL_CV      <- "#2E6DA4"  # same as LINE_COLOR in plot_fluctuation_proportion.R
COL_SUM_STD <- "#D95F02"

LINE_WIDTH <- 0.6
ERROR_WIDTH <- 0.12
ERROR_LW <- 0.40
POINT_SIZE <- 1.8

PLOT_CONFIG <- tribble(
  ~experiment,   ~window, ~filename,            ~cv_y_min, ~cv_y_max, ~cv_y_step,
  "temperature", "full",  "temperature_full",   0,         0.3,       0.1,
  "mortality",   "early", "mortality_early",    0,         0.6,       0.2,
  "temperature", "last4", "temperature_last4",  0,         0.3,       0.1,
  "mortality",   "last4", "mortality_last4",   0,         0.6,       0.2
)

X_LABELS <- list(
  temperature = c(
    "W1" = "10\u00B0C",
    "W2" = "20\u00B0C",
    "W3" = "30\u00B0C",
    "W4" = "40\u00B0C",
    "W5" = "50\u00B0C"
  ),
  mortality = c(
    "W1" = expression(10^1),
    "W2" = expression(10^2),
    "W3" = expression(10^3),
    "W4" = expression(10^4),
    "W5" = expression(10^5)
  )
)

theme_pub <- function() {
  theme_classic(base_size = FONT_AX, base_family = FONT_FAMILY) +
    theme(
      axis.line              = element_blank(),
      panel.border           = element_rect(linewidth = BORDER_SIZE,
                                            color = "black", fill = NA),
      axis.ticks             = element_line(linewidth = TICK_SIZE,
                                            color = "black"),
      axis.ticks.x.top       = element_line(linewidth = TICK_SIZE,
                                            color = "black"),
      axis.ticks.y.right     = element_line(linewidth = TICK_SIZE,
                                            color = "black"),
      axis.ticks.length               = unit(TICK_LEN, "mm"),
      axis.ticks.length.x.top         = unit(TICK_LEN, "mm"),
      axis.ticks.length.y.right       = unit(TICK_LEN, "mm"),
      axis.text.x.top        = element_blank(),
      axis.text.y.right      = element_text(size = FONT_AX,
                                            color = "black",
                                            margin = margin(l = 3)),
      axis.text.x            = element_text(size = FONT_AX,
                                            color = "black",
                                            margin = margin(t = 3)),
      axis.text.y            = element_text(size = FONT_AX,
                                            color = "black",
                                            margin = margin(r = 3)),
      axis.title             = element_blank(),
      legend.position        = "none",
      plot.background        = element_rect(fill = "white", color = NA),
      panel.background       = element_rect(fill = "white", color = NA),
      plot.margin            = margin(3, 3, 3, 3, "mm")
    )
}

nice_upper <- function(x) {
  if (!is.finite(x) || x <= 0) return(1)
  step <- case_when(
    x <= 0.10 ~ 0.02,
    x <= 0.50 ~ 0.10,
    x <= 1.00 ~ 0.20,
    TRUE      ~ 0.50
  )
  ceiling(x * 1.12 / step) * step
}

scale_x_condition <- function(experiment) {
  scale_x_continuous(
    breaks = 1:5,
    labels = X_LABELS[[experiment]],
    expand = expansion(mult = c(0.08, 0.08)),
    sec.axis = dup_axis(labels = NULL, name = NULL)
  )
}

read_fluctuation_window <- function(window) {
  csv_path <- file.path(FLUC_DIR, window, "community_level.csv")
  if (!file.exists(csv_path)) {
    stop("Cannot find compute_fluctuations output: ", csv_path)
  }

  read_csv(csv_path, show_col_types = FALSE)
}

summarise_metric_means <- function(df, experiment, window) {
  df %>%
    filter(experiment == !!experiment) %>%
    mutate(x_pos = as.integer(str_remove(condition, "^W"))) %>%
    group_by(experiment, condition, x_pos) %>%
    summarise(
      n = n(),
      community_cv_mean = mean(community_cv, na.rm = TRUE),
      community_cv_sem = sd(community_cv, na.rm = TRUE) / sqrt(n()),
      sum_abs_std_mean = mean(sum_abs_std, na.rm = TRUE),
      sum_abs_std_sem = sd(sum_abs_std, na.rm = TRUE) / sqrt(n()),
      n_collapsed = sum(collapsed, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(window = window, .before = condition) %>%
    arrange(x_pos)
}

plot_dual_axis <- function(summary_df, cfg) {
  experiment <- cfg$experiment
  cv_y_min <- cfg$cv_y_min
  cv_y_max <- cfg$cv_y_max
  cv_y_step <- cfg$cv_y_step

  sum_y_max <- nice_upper(
    max(summary_df$sum_abs_std_mean + summary_df$sum_abs_std_sem,
        na.rm = TRUE)
  )
  scale_factor <- cv_y_max / sum_y_max

  plot_df <- summary_df %>%
    mutate(
      community_cv = community_cv_mean,
      community_cv_ymin = community_cv_mean - community_cv_sem,
      community_cv_ymax = community_cv_mean + community_cv_sem,
      sum_abs_std_scaled = sum_abs_std_mean * scale_factor,
      sum_abs_std_ymin_scaled = pmax(0, sum_abs_std_mean - sum_abs_std_sem) *
        scale_factor,
      sum_abs_std_ymax_scaled = (sum_abs_std_mean + sum_abs_std_sem) *
        scale_factor
    )

  y_breaks_left <- seq(cv_y_min, cv_y_max, by = cv_y_step)

  ggplot(plot_df, aes(x = x_pos)) +
    geom_line(
      aes(y = community_cv, color = "community_CV"),
      linewidth = LINE_WIDTH
    ) +
    geom_errorbar(
      aes(ymin = community_cv_ymin, ymax = community_cv_ymax,
          color = "community_CV"),
      width = ERROR_WIDTH,
      linewidth = ERROR_LW
    ) +
    geom_point(
      aes(y = community_cv, color = "community_CV"),
      shape = 16,
      size = POINT_SIZE
    ) +
    geom_line(
      aes(y = sum_abs_std_scaled, color = "sum of stdNi"),
      linewidth = LINE_WIDTH
    ) +
    geom_errorbar(
      aes(ymin = sum_abs_std_ymin_scaled, ymax = sum_abs_std_ymax_scaled,
          color = "sum of stdNi"),
      width = ERROR_WIDTH,
      linewidth = ERROR_LW
    ) +
    geom_point(
      aes(y = sum_abs_std_scaled, color = "sum of stdNi"),
      shape = 16,
      size = POINT_SIZE
    ) +
    scale_color_manual(
      values = c(
        "community_CV" = COL_CV,
        "sum of stdNi" = COL_SUM_STD
      ),
      breaks = c("community_CV", "sum of stdNi")
    ) +
    scale_x_condition(experiment) +
    scale_y_continuous(
      limits = c(cv_y_min, cv_y_max),
      breaks = y_breaks_left,
      labels = label_number(accuracy = 0.01),
      expand = expansion(mult = c(0, 0)),
      name = NULL,
      sec.axis = sec_axis(
        trans = ~ . / scale_factor,
        name = NULL,
        breaks = pretty(c(0, sum_y_max), n = 4),
        labels = label_number(accuracy = 0.01)
      )
    ) +
    labs(x = NULL, y = NULL) +
    theme_pub() +
    theme(
      legend.position = "none"
    )
}

make_metric_legend <- function() {
  legend_df <- tibble(
    x = rep(1:2, 2),
    y = c(1, 1, 2, 2),
    metric = factor(
      rep(c("community_CV", "sum of stdNi"), each = 2),
      levels = c("community_CV", "sum of stdNi")
    )
  )

  ggplot(legend_df, aes(x = x, y = y, color = metric)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(
      values = c(
        "community_CV" = COL_CV,
        "sum of stdNi" = COL_SUM_STD
      )
    ) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, title = NULL)) +
    theme_void(base_family = FONT_FAMILY) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = FONT_LEGEND,
                                 family = FONT_FAMILY,
                                 color = "black"),
      legend.key.width = unit(10, "mm"),
      legend.key.height = unit(4, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(0, 0, 0, 0)
    )
}

all_summaries <- list()

for (i in seq_len(nrow(PLOT_CONFIG))) {
  cfg <- PLOT_CONFIG[i, ]
  df_window <- read_fluctuation_window(cfg$window)
  summary_df <- summarise_metric_means(
    df = df_window,
    experiment = cfg$experiment,
    window = cfg$window
  )

  if (nrow(summary_df) == 0) next

  all_summaries[[length(all_summaries) + 1]] <- summary_df

  p <- plot_dual_axis(summary_df, cfg)

  out_pdf <- file.path(
    FIG_DIR,
    paste0(cfg$filename, "_community_cv_sum_std_dual_axis.pdf")
  )
  out_png <- file.path(
    FIG_DIR,
    paste0(cfg$filename, "_community_cv_sum_std_dual_axis.png")
  )

  ggsave(
    filename = out_pdf,
    plot = p,
    width = W_PLOT,
    height = H_PLOT,
    units = "mm",
    device = cairo_pdf
  )
  ggsave(
    filename = out_png,
    plot = p,
    width = W_PLOT,
    height = H_PLOT,
    units = "mm",
    dpi = 300
  )

  cat("Saved:", out_pdf, "\n")
}

if (length(all_summaries) > 0) {
  write_csv(bind_rows(all_summaries), SUMMARY_OUT)
  cat("Summary saved:", SUMMARY_OUT, "\n")
}

legend_plot <- cowplot::get_legend(make_metric_legend())
ggsave(
  filename = file.path(FIG_DIR, "community_cv_sum_std_dual_axis_legend.pdf"),
  plot = legend_plot,
  width = W_LEGEND,
  height = H_LEGEND,
  units = "mm",
  device = cairo_pdf
)
ggsave(
  filename = file.path(FIG_DIR, "community_cv_sum_std_dual_axis_legend.png"),
  plot = legend_plot,
  width = W_LEGEND,
  height = H_LEGEND,
  units = "mm",
  dpi = 300
)

cat("Done.\n")
