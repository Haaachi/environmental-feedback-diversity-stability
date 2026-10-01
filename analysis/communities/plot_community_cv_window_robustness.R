community_workspace <- function() {
  override <- Sys.getenv("COMMUNITY_WORKSPACE", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script)) stop("Run with Rscript or set COMMUNITY_WORKSPACE.")
  dirname(normalizePath(sub("^--file=", "", script[[1]]), mustWork = TRUE))
}

# ===========================================================
# plot_community_cv_window_robustness.R
#
# Robustness plots for the choice of community_CV time-window start.
# Input is produced by compute_community_cv_window_robustness.py.
# ===========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

BASE_DIR <- community_workspace()
PROC_DIR <- file.path(BASE_DIR, "processed", "community_cv_window_robustness")
FIG_DIR  <- file.path(BASE_DIR, "figures", "community_cv_window_robustness")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

CV_THRESHOLD <- 0.25

CLASS_COLORS <- c(
  "Fluctuation" = "#F4A261",
  "Stable"      = "#8E6BBE"
)

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

FONT_FAMILY <- "Arial"
FONT_AX <- 10
BORDER_SIZE <- 0.25
TICK_SIZE <- 0.20
TICK_LEN <- -1.0
W_PLOT <- 82
H_PLOT <- 72

theme_pub_axis <- function(show_legend = TRUE) {
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
      axis.text.y.right      = element_blank(),
      axis.text             = element_text(size = FONT_AX, color = "black"),
      axis.title            = element_text(size = FONT_AX, color = "black"),
      legend.position       = if (show_legend) c(0.04, 0.96) else "none",
      legend.justification  = c(0, 1),
      legend.title          = element_blank(),
      legend.text           = element_text(size = FONT_AX, color = "black"),
      legend.background     = element_blank(),
      plot.background       = element_rect(fill = "white", color = NA),
      panel.background      = element_rect(fill = "white", color = NA)
    )
}

plot_window_curve <- function(summary_df, experiment, out_file) {
  d <- summary_df %>%
    filter(experiment == !!experiment,
           !is.na(mean_community_cv)) %>%
    mutate(
      ymin = pmax(0, mean_community_cv - sem_community_cv),
      ymax = mean_community_cv + sem_community_cv,
      fluctuation_class = factor(fluctuation_class,
                                 levels = c("Fluctuation", "Stable"))
    )

  if (nrow(d) == 0) return(invisible(NULL))

  y_upper <- max(d$ymax, na.rm = TRUE)
  y_upper <- max(CV_THRESHOLD * 1.4, ceiling(y_upper * 10) / 10)
  y_breaks <- pretty(c(0, y_upper), n = 5)
  y_breaks <- y_breaks[y_breaks <= y_upper]

  p <- ggplot(d, aes(x = first_day, y = mean_community_cv,
                     color = fluctuation_class,
                     group = fluctuation_class)) +
    geom_hline(yintercept = CV_THRESHOLD, linewidth = 0.25,
               linetype = "dashed", color = "grey55") +
    geom_line(linewidth = 0.62) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax),
                  width = 0.12, linewidth = 0.36, alpha = 0.65) +
    geom_point(size = 1.8) +
    scale_color_manual(values = CLASS_COLORS, drop = FALSE) +
    scale_x_continuous(
      breaks = sort(unique(d$first_day)),
      expand = expansion(mult = c(0.04, 0.06)),
      sec.axis = dup_axis(labels = NULL, name = NULL)
    ) +
    scale_y_continuous(
      limits = c(0, y_upper),
      breaks = y_breaks,
      expand = expansion(mult = c(0, 0.03)),
      sec.axis = dup_axis(labels = NULL, name = NULL)
    ) +
    labs(x = "First day of time window", y = "Community CV") +
    theme_pub_axis(show_legend = FALSE)

  ggsave(out_file, p, width = W_PLOT, height = H_PLOT, units = "mm",
         device = cairo_pdf)
}

plot_window_curve_by_condition <- function(summary_df, experiment, condition) {
  d <- summary_df %>%
    filter(experiment == !!experiment,
           condition == !!condition,
           !is.na(mean_community_cv)) %>%
    mutate(
      ymin = pmax(0, mean_community_cv - sem_community_cv),
      ymax = mean_community_cv + sem_community_cv,
      fluctuation_class = factor(fluctuation_class,
                                 levels = c("Fluctuation", "Stable"))
    )

  if (nrow(d) == 0) return(invisible(NULL))

  y_upper <- max(d$ymax, na.rm = TRUE)
  y_upper <- max(CV_THRESHOLD * 1.4, ceiling(y_upper * 10) / 10)
  y_breaks <- pretty(c(0, y_upper), n = 5)
  y_breaks <- y_breaks[y_breaks <= y_upper]

  p <- ggplot(d, aes(x = first_day, y = mean_community_cv,
                     color = fluctuation_class,
                     group = fluctuation_class)) +
    geom_hline(yintercept = CV_THRESHOLD, linewidth = 0.25,
               linetype = "dashed", color = "grey55") +
    geom_line(linewidth = 0.62) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax),
                  width = 0.12, linewidth = 0.36, alpha = 0.65) +
    geom_point(size = 1.8) +
    scale_color_manual(values = CLASS_COLORS, drop = FALSE) +
    scale_x_continuous(
      breaks = sort(unique(d$first_day)),
      expand = expansion(mult = c(0.04, 0.06)),
      sec.axis = dup_axis(labels = NULL, name = NULL)
    ) +
    scale_y_continuous(
      limits = c(0, y_upper),
      breaks = y_breaks,
      expand = expansion(mult = c(0, 0.03)),
      sec.axis = dup_axis(labels = NULL, name = NULL)
    ) +
    labs(x = "First day of time window", y = "Community CV") +
    theme_pub_axis(show_legend = FALSE)

  cond_dir <- file.path(FIG_DIR, "window_start_by_condition")
  dir.create(cond_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(cond_dir,
                   sprintf("%s_%s_window_start_robustness.pdf",
                           experiment, condition)),
         p, width = W_PLOT, height = H_PLOT, units = "mm",
         device = cairo_pdf)
}

plot_temperature_replicate_distance <- function(dist_df) {
  d <- dist_df %>%
    filter(is.finite(community_cv_max_last3),
           is.finite(mean_replicate_bc_last3)) %>%
    mutate(
      fluctuation_class = factor(fluctuation_class,
                                 levels = c("Fluctuation", "Stable")),
      condition = factor(condition, levels = paste0("W", 1:5))
    )

  if (nrow(d) == 0) return(invisible(NULL))

  p <- ggplot(d, aes(x = community_cv_max_last3,
                     y = mean_replicate_bc_last3,
                     color = fluctuation_class)) +
    geom_vline(xintercept = CV_THRESHOLD, linewidth = 0.25,
               linetype = "dashed", color = "grey55") +
    geom_point(size = 1.9, alpha = 0.82) +
    scale_color_manual(values = CLASS_COLORS, drop = FALSE) +
    scale_x_continuous(
      limits = c(0, max(1, max(d$community_cv_max_last3, na.rm = TRUE) * 1.05)),
      breaks = pretty_breaks(5),
      expand = expansion(mult = c(0, 0.02)),
      sec.axis = dup_axis(labels = NULL, name = NULL)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = expansion(mult = c(0, 0)),
      sec.axis = dup_axis(labels = NULL, name = NULL)
    ) +
    labs(x = "Community CV",
         y = "Replicate distance") +
    theme_pub_axis(show_legend = FALSE)

  ggsave(file.path(FIG_DIR,
                   "temperature_replicate_distance_vs_max_community_CV.pdf"),
         p, width = W_PLOT, height = H_PLOT, units = "mm",
         device = cairo_pdf)

  cond_dir <- file.path(FIG_DIR, "temperature_replicate_distance_by_condition")
  dir.create(cond_dir, showWarnings = FALSE, recursive = TRUE)
  for (cond in paste0("W", 1:5)) {
    d_cond <- d %>% filter(condition == cond)
    if (nrow(d_cond) == 0) next

    p_cond <- ggplot(d_cond, aes(x = community_cv_max_last3,
                                 y = mean_replicate_bc_last3,
                                 color = fluctuation_class)) +
      geom_vline(xintercept = CV_THRESHOLD, linewidth = 0.25,
                 linetype = "dashed", color = "grey55") +
      geom_point(size = 1.9, alpha = 0.82) +
      scale_color_manual(values = CLASS_COLORS, drop = FALSE) +
      scale_x_continuous(
        limits = c(0, max(1, max(d$community_cv_max_last3, na.rm = TRUE) * 1.05)),
        breaks = pretty_breaks(5),
        expand = expansion(mult = c(0, 0.02)),
        sec.axis = dup_axis(labels = NULL, name = NULL)
      ) +
      scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        expand = expansion(mult = c(0, 0)),
        sec.axis = dup_axis(labels = NULL, name = NULL)
      ) +
      labs(x = "Community CV",
           y = "Replicate distance") +
      theme_pub_axis(show_legend = FALSE)

    ggsave(file.path(cond_dir,
                     sprintf("temperature_%s_replicate_distance_vs_max_community_CV.pdf",
                             cond)),
           p_cond, width = W_PLOT, height = H_PLOT, units = "mm",
           device = cairo_pdf)
  }
}

summary_path <- file.path(PROC_DIR, "window_cv_summary.csv")
summary_condition_path <- file.path(PROC_DIR, "window_cv_summary_by_condition.csv")
dist_path <- file.path(PROC_DIR, "temperature_replicate_bc_vs_max_cv.csv")

if (!file.exists(summary_path) || !file.exists(summary_condition_path) ||
    !file.exists(dist_path)) {
  stop("Missing robustness CSVs. Run compute_community_cv_window_robustness.py first.",
       call. = FALSE)
}

summary_df <- read_csv(summary_path, show_col_types = FALSE)
summary_condition_df <- read_csv(summary_condition_path, show_col_types = FALSE)
dist_df <- read_csv(dist_path, show_col_types = FALSE)

cat("Plotting pooled window-start robustness curves...\n")
for (experiment in c("mortality", "temperature")) {
  plot_window_curve(
    summary_df,
    experiment,
    file.path(FIG_DIR,
              sprintf("%s_window_start_robustness_pooled_conditions.pdf",
                      experiment))
  )
}

cat("Plotting condition-specific window-start robustness curves...\n")
for (experiment in c("mortality", "temperature")) {
  for (condition in paste0("W", 1:5)) {
    plot_window_curve_by_condition(summary_condition_df, experiment, condition)
  }
}

cat("Plotting temperature replicate-distance scatter...\n")
plot_temperature_replicate_distance(dist_df)

cat("Done. Output directory: figures/community_cv_window_robustness\n")
