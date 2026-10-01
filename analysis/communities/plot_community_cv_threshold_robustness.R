community_workspace <- function() {
  override <- Sys.getenv("COMMUNITY_WORKSPACE", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script)) stop("Run with Rscript or set COMMUNITY_WORKSPACE.")
  dirname(normalizePath(sub("^--file=", "", script[[1]]), mustWork = TRUE))
}

# ===========================================================
# plot_community_cv_threshold_robustness.R
#
# Legacy robustness script. Binary classes use the main fluctuation criterion:
# community_cv >= 0.25.
# ===========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

BASE_DIR <- community_workspace()
PROC_DIR <- file.path(BASE_DIR, "processed")
FIG_DIR  <- file.path(BASE_DIR, "figures", "community_cv_threshold_robustness")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

COMMUNITY_CV_THRESHOLD <- 0.25

is_fluctuating_composite <- function(experiment, community_cv) {
  !is.na(community_cv) & community_cv >= COMMUNITY_CV_THRESHOLD
}

EXP_WINDOW_MAP <- list(
  mortality   = "early",
  temperature = "full"
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

COLOR_STABLE <- "#9B8EC4"
COLOR_FLUCT  <- "#F4A460"
LINE_COLOR   <- "#2E6DA4"

FONT_FAMILY <- "Arial"
FONT_AX     <- 10
BORDER_SIZE <- 0.25
TICK_SIZE   <- 0.20
TICK_LEN    <- -1.0

W_PLOT <- 80
H_PLOT <- 70

DIVERSITY_INPUTS <- list(
  alpha = list(
    file = "alpha_diversity.csv",
    metrics = c("richness", "shannon", "effective_shannon")
  ),
  mean_daily = list(
    file = "mean_daily_diversity.csv",
    metrics = c("mean_daily_richness",
                "mean_daily_shannon",
                "mean_daily_effective_shannon")
  ),
  gamma = list(
    file = "gamma_diversity.csv",
    metrics = c("gamma_richness",
                "gamma_shannon",
                "gamma_effective_shannon")
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
      axis.text.y.right      = element_blank(),
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

scale_x_condition <- function(experiment) {
  scale_x_continuous(
    breaks   = 1:5,
    labels   = X_LABELS[[experiment]],
    expand   = expansion(mult = c(0.08, 0.08)),
    sec.axis = dup_axis(labels = NULL, name = NULL)
  )
}

read_cv_table <- function(experiment) {
  window <- EXP_WINDOW_MAP[[experiment]]
  path <- file.path(PROC_DIR, "fluctuations", window, "community_level.csv")
  if (!file.exists(path)) {
    warning(sprintf("Missing community_CV table: %s", path))
    return(tibble())
  }

  read_csv(path, show_col_types = FALSE) %>%
    filter(experiment == !!experiment) %>%
    mutate(
      x_pos = as.integer(sub("W", "", condition)),
      threshold = COMMUNITY_CV_THRESHOLD,
      fluctuation_class = if_else(
        is_fluctuating_composite(experiment, community_cv),
        "Fluctuation",
        "Stable"
      ),
      fluctuation_class = factor(fluctuation_class,
                                 levels = c("Fluctuation", "Stable"))
    )
}

plot_threshold_proportion <- function(experiment) {
  out_dir <- file.path(FIG_DIR, "proportion")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  window <- EXP_WINDOW_MAP[[experiment]]
  cv_df <- read_cv_table(experiment) %>%
    filter(is.finite(community_cv))
  if (nrow(cv_df) == 0) return(invisible(NULL))

  summary_prop <- cv_df %>%
    group_by(condition, x_pos) %>%
    summarise(
      n_osc = sum(fluctuation_class == "Fluctuation", na.rm = TRUE),
      n_total = n(),
      .groups = "drop"
    ) %>%
    mutate(
      proportion = n_osc / n_total,
      sem = sqrt(proportion * (1 - proportion) / n_total),
      ymin = pmax(0, proportion - sem),
      ymax = pmin(1, proportion + sem)
    ) %>%
    arrange(x_pos)

  write_csv(
    summary_prop,
    file.path(out_dir,
              sprintf("%s_%s_community_cv_0.25_proportion_summary.csv",
                      experiment, window))
  )

  p <- ggplot(summary_prop, aes(x = x_pos, y = proportion)) +
    geom_line(color = LINE_COLOR, linewidth = 0.6) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax),
                  width = 0.12, linewidth = 0.40, color = LINE_COLOR) +
    geom_point(shape = 16, size = 1.8, color = LINE_COLOR) +
    scale_x_condition(experiment) +
    scale_y_continuous(
      limits = c(0, 1.0),
      breaks = seq(0, 1.0, by = 0.2),
      labels = sprintf("%.1f", seq(0, 1.0, by = 0.2)),
      expand = expansion(mult = c(0, 0)),
      sec.axis = dup_axis(labels = NULL, name = NULL)
    ) +
    labs(x = NULL, y = NULL) +
    theme_pub()

  ggsave(
    file.path(out_dir,
              sprintf("%s_%s_community_cv_0.25_proportion.pdf",
                      experiment, window)),
    p, width = W_PLOT, height = H_PLOT, units = "mm",
    device = cairo_pdf
  )
}

plot_flucstable_diversity_metric <- function(div_df, cv_df, metric,
                                             experiment) {
  joined <- div_df %>%
    left_join(
      cv_df %>%
        select(experiment, condition, community, replica,
               community_cv, fluctuation_class),
      by = c("experiment", "condition", "community", "replica")
    ) %>%
    filter(!is.na(fluctuation_class),
           is.finite(.data[[metric]])) %>%
    mutate(x_pos = as.integer(sub("W", "", condition)))

  if (nrow(joined) == 0) return(NULL)

  summary_df <- joined %>%
    group_by(condition, x_pos, fluctuation_class) %>%
    summarise(
      n = n(),
      mean_val = mean(.data[[metric]], na.rm = TRUE),
      sem_val = sd(.data[[metric]], na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      ymin = pmax(0, mean_val - sem_val),
      ymax = mean_val + sem_val
    ) %>%
    arrange(x_pos, fluctuation_class)

  y_upper <- max(summary_df$ymax, na.rm = TRUE)
  y_upper <- ceiling(y_upper * 10) / 10
  y_upper <- max(y_upper, 0.2)
  y_breaks <- pretty(c(0, y_upper), n = 5)
  y_breaks <- y_breaks[y_breaks <= y_upper]

  p <- ggplot(summary_df,
              aes(x = x_pos, y = mean_val,
                  color = fluctuation_class,
                  group = fluctuation_class)) +
    geom_line(linewidth = 0.58) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax),
                  width = 0.12, linewidth = 0.35, alpha = 0.72) +
    geom_point(size = 1.8) +
    scale_color_manual(values = c("Fluctuation" = COLOR_FLUCT,
                                  "Stable" = COLOR_STABLE),
                       drop = FALSE) +
    scale_x_condition(experiment) +
    scale_y_continuous(
      limits = c(0, y_upper),
      breaks = y_breaks,
      expand = expansion(mult = c(0, 0.03)),
      sec.axis = dup_axis(labels = NULL, name = NULL)
    ) +
    labs(x = NULL, y = NULL) +
    theme_pub()

  attr(p, "joined_data") <- joined
  attr(p, "summary_data") <- summary_df
  p
}

plot_threshold_diversity <- function(experiment) {
  window <- EXP_WINDOW_MAP[[experiment]]
  cv_df <- read_cv_table(experiment)
  if (nrow(cv_df) == 0) return(invisible(NULL))

  out_dir <- file.path(FIG_DIR, "diversity_flucstable", experiment)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  all_joined <- list()
  all_summary <- list()

  for (input_name in names(DIVERSITY_INPUTS)) {
    input_cfg <- DIVERSITY_INPUTS[[input_name]]
    path <- file.path(PROC_DIR, "diversity", window, input_cfg$file)
    if (!file.exists(path)) next

    div_df <- read_csv(path, show_col_types = FALSE) %>%
      filter(experiment == !!experiment)

    input_dir <- file.path(out_dir, input_name)
    dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)

    for (metric in input_cfg$metrics) {
      if (!metric %in% names(div_df)) next

      p <- plot_flucstable_diversity_metric(div_df, cv_df, metric,
                                            experiment)
      if (is.null(p)) next

      ggsave(
        file.path(input_dir,
                  sprintf("%s_community_cv_0.25_flucstable.pdf", metric)),
        p, width = W_PLOT, height = H_PLOT, units = "mm",
        device = cairo_pdf
      )

      joined <- attr(p, "joined_data") %>%
        mutate(input = input_name, metric = metric)
      summary_df <- attr(p, "summary_data") %>%
        mutate(input = input_name, metric = metric)
      all_joined[[length(all_joined) + 1]] <- joined
      all_summary[[length(all_summary) + 1]] <- summary_df
    }
  }

  if (length(all_joined) > 0) {
    write_csv(bind_rows(all_joined),
              file.path(out_dir,
                        "community_cv_0.25_diversity_flucstable_values.csv"))
  }
  if (length(all_summary) > 0) {
    write_csv(bind_rows(all_summary),
              file.path(out_dir,
                        "community_cv_0.25_diversity_flucstable_summary.csv"))
  }
}

plot_community_cv_threshold_robustness <- function() {
  for (experiment in names(EXP_WINDOW_MAP)) {
    cat(sprintf("Plotting species-instability proportion: %s\n",
                experiment))
    plot_threshold_proportion(experiment)

    cat(sprintf("Plotting species-instability diversity fluc/stable: %s\n",
                experiment))
    plot_threshold_diversity(experiment)
  }
}

plot_community_cv_threshold_robustness()

cat("Done. Output directory: figures/community_cv_threshold_robustness\n")
