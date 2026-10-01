community_workspace <- function() {
  override <- Sys.getenv("COMMUNITY_WORKSPACE", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script)) stop("Run with Rscript or set COMMUNITY_WORKSPACE.")
  dirname(normalizePath(sub("^--file=", "", script[[1]]), mustWork = TRUE))
}

# ===========================================================
# plot_diversity_cv_halves.R
#
# Purpose:
#   Plot diversity metrics for the high-community-CV half and the
#   low-community-CV half across the full experiment x window.
#
# This is not a fixed-threshold stable/fluctuation classification.
# Within each experiment x window, communities/replicates are
# ranked by community_CV:
#   - High CV half: top 50%
#   - Low CV half: bottom 50%
# In the temperature experiment, collapsed 50°C communities are excluded before
# ranking and therefore do not enter the corresponding significance tests.
#
# Colors reuse the previous visual language:
#   - High CV half: fluctuation orange
#   - Low CV half: stable purple
#
# Main analysis windows:
#   mortality   -> early
#   temperature -> full
#
# Output:
#   figures/diversity_cv_halves_global_rank/{experiment}/{window}/
# ===========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

BASE_DIR <- community_workspace()
PROC_DIR <- file.path(BASE_DIR, "processed")
FIG_DIR  <- file.path(BASE_DIR, "figures", "diversity_cv_halves_global_rank")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

EXP_WINDOW_MAP <- list(
  mortality   = c("early"),
  temperature = c("full")
)

X_LABELS <- list(
  mortality = c(
    "W1" = expression(10^{-1}),
    "W2" = expression(10^{-2}),
    "W3" = expression(10^{-3}),
    "W4" = expression(10^{-4}),
    "W5" = expression(10^{-5})
  ),
  temperature = c(
    "W1" = "10\u00b0C",
    "W2" = "20\u00b0C",
    "W3" = "30\u00b0C",
    "W4" = "40\u00b0C",
    "W5" = "50\u00b0C"
  )
)

ALPHA_METRICS <- list(
  richness                  = list(ylim = NULL),
  survival_fraction         = list(ylim = c(0, 1)),
  shannon                   = list(ylim = NULL),
  effective_shannon         = list(ylim = NULL),
  simpson                   = list(ylim = c(0, 1)),
  evenness                  = list(ylim = c(0, 1))
)

GAMMA_METRICS <- list(
  gamma_richness            = list(ylim = NULL),
  gamma_shannon             = list(ylim = NULL),
  gamma_effective_shannon   = list(ylim = NULL),
  gamma_simpson             = list(ylim = c(0, 1)),
  gamma_evenness            = list(ylim = c(0, 1))
)

MEAN_DAILY_METRICS <- list(
  mean_daily_richness          = list(ylim = NULL),
  mean_daily_shannon           = list(ylim = NULL),
  mean_daily_effective_shannon = list(ylim = NULL),
  mean_daily_simpson           = list(ylim = c(0, 1)),
  mean_daily_evenness          = list(ylim = c(0, 1)),
  mean_daily_survival_fraction = list(ylim = c(0, 1))
)

COLOR_LOW_CV  <- "#9B8EC4"
COLOR_HIGH_CV <- "#F4A460"
LINE_COLOR    <- "#333333"

FONT_FAMILY <- "Arial"
FONT_AX     <- 10
BORDER_SIZE <- 0.25
TICK_SIZE   <- 0.20
TICK_LEN    <- -1.0

W_PLOT <- 80
H_PLOT <- 70

LINE_WIDTH <- 0.55
ERR_WIDTH  <- 0.12
ERR_LW     <- 0.35
MEAN_SIZE  <- 1.9
SIG_TEXT_SIZE <- 3.0

p_to_sig_label <- function(p) {
  case_when(
    is.na(p)  ~ NA_character_,
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ NA_character_
  )
}

theme_pub <- function() {
  theme_classic(base_size = FONT_AX, base_family = FONT_FAMILY) +
    theme(
      axis.line = element_blank(),
      panel.border = element_rect(
        linewidth = BORDER_SIZE, color = "black", fill = NA
      ),
      axis.ticks = element_line(linewidth = TICK_SIZE, color = "black"),
      axis.ticks.x.top = element_line(linewidth = TICK_SIZE, color = "black"),
      axis.ticks.y.right = element_line(linewidth = TICK_SIZE, color = "black"),
      axis.ticks.length = unit(TICK_LEN, "mm"),
      axis.ticks.length.x.top = unit(TICK_LEN, "mm"),
      axis.ticks.length.y.right = unit(TICK_LEN, "mm"),
      axis.text.x.top = element_blank(),
      axis.text.y.right = element_blank(),
      axis.text = element_text(size = FONT_AX, color = "black"),
      axis.title = element_blank(),
      legend.position = "none",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(3, 3, 3, 3, "mm")
    )
}

ensure_collapsed_col <- function(df) {
  if (!"collapsed" %in% names(df)) {
    df <- df %>% mutate(collapsed = FALSE)
  } else {
    df <- df %>% mutate(collapsed = if_else(is.na(collapsed), FALSE, collapsed))
  }
  df
}

is_valid_combo <- function(experiment, window) {
  allowed <- EXP_WINDOW_MAP[[experiment]]
  !is.null(allowed) && window %in% allowed
}

read_cv_halves <- function(experiment, window) {
  fluc_path <- file.path(PROC_DIR, "fluctuations", window, "community_level.csv")
  if (!file.exists(fluc_path)) return(tibble())

  read_csv(fluc_path, show_col_types = FALSE) %>%
    filter(experiment == !!experiment) %>%
    ensure_collapsed_col() %>%
    mutate(
      community_cv = as.numeric(community_cv),
      x_pos = as.integer(sub("W", "", condition))
    ) %>%
    filter(
      !is.na(community_cv),
      !(experiment == "temperature" & condition == "W5" & collapsed == TRUE)
    ) %>%
    arrange(desc(community_cv)) %>%
    mutate(
      cv_rank = row_number(),
      n_rank_set = n(),
      cv_half = if_else(
        cv_rank <= ceiling(n_rank_set / 2),
        "High CV half",
        "Low CV half"
      ),
      cv_half = factor(cv_half, levels = c("Low CV half", "High CV half"))
    ) %>%
    select(
      experiment, condition, community, replica, community_cv,
      collapsed, x_pos, cv_rank, n_rank_set, cv_half
    )
}

plot_cv_half_metric <- function(df_metric, cv_df, metric_col, experiment,
                                ylim_fixed = NULL) {
  if (nrow(cv_df) == 0) return(NULL)

  join_keys <- intersect(
    c("experiment", "condition", "community", "replica"),
    names(df_metric)
  )

  df_joined <- df_metric %>%
    left_join(cv_df, by = join_keys) %>%
    filter(
      !is.na(cv_half),
      !is.na(.data[[metric_col]])
    )

  if (nrow(df_joined) == 0) return(NULL)

  summary_df <- df_joined %>%
    group_by(condition, x_pos, cv_half) %>%
    summarise(
      n = sum(!is.na(.data[[metric_col]])),
      mean_val = mean(.data[[metric_col]], na.rm = TRUE),
      sem_val = sd(.data[[metric_col]], na.rm = TRUE) / sqrt(n),
      .groups = "drop"
    ) %>%
    mutate(
      ymin = mean_val - sem_val,
      ymax = mean_val + sem_val
    )

  if (nrow(summary_df) == 0) return(NULL)

  if (!is.null(ylim_fixed)) {
    y_lower <- ylim_fixed[1]
    y_upper <- ylim_fixed[2]
    y_breaks <- pretty(c(y_lower, y_upper), n = 4)
    y_breaks <- y_breaks[y_breaks >= y_lower & y_breaks <= y_upper]
  } else {
    y_lower <- 0
    y_upper <- max(summary_df$ymax, df_joined[[metric_col]], na.rm = TRUE) * 1.05
    y_breaks <- pretty(c(y_lower, y_upper), n = 4)
    y_upper <- max(y_breaks)
  }

  sig_df <- test_high_half_greater(df_joined, metric_col, "plot") %>%
    filter(!is.na(wilcox_p_greater), wilcox_p_greater < 0.05) %>%
    mutate(
      x_pos = as.integer(sub("W", "", condition)),
      sig_label = p_to_sig_label(wilcox_p_greater),
      higher_half = if_else(mean_high > mean_low, "High CV half", "Low CV half"),
      color = if_else(higher_half == "High CV half", COLOR_HIGH_CV, COLOR_LOW_CV)
    ) %>%
    filter(!is.na(sig_label))

  y_range <- y_upper - y_lower
  if (!is.finite(y_range) || y_range <= 0) y_range <- 1
  if (nrow(sig_df) > 0) {
    max_per_cond <- summary_df %>%
      group_by(condition, x_pos) %>%
      summarise(max_upper = max(ymax, na.rm = TRUE), .groups = "drop")

    sig_df <- sig_df %>%
      select(condition, x_pos, sig_label, color) %>%
      inner_join(max_per_cond, by = c("condition", "x_pos")) %>%
      mutate(
        label_y = max_upper + y_range * 0.08,
        label_y = pmin(label_y, y_upper - y_range * 0.02)
      )
  }

  p <- ggplot(summary_df, aes(x = x_pos, y = mean_val, color = cv_half)) +
    geom_line(linewidth = LINE_WIDTH) +
    geom_errorbar(
      aes(ymin = ymin, ymax = ymax),
      width = ERR_WIDTH,
      linewidth = ERR_LW
    ) +
    geom_point(shape = 16, size = MEAN_SIZE) +
    scale_color_manual(
      values = c(
        "Low CV half" = COLOR_LOW_CV,
        "High CV half" = COLOR_HIGH_CV
      )
    ) +
    scale_x_continuous(
      breaks = 1:5,
      labels = X_LABELS[[experiment]],
      expand = expansion(mult = c(0.10, 0.15)),
      sec.axis = dup_axis(labels = NULL, name = NULL)
    ) +
    scale_y_continuous(
      limits = c(y_lower, y_upper),
      breaks = y_breaks,
      expand = expansion(mult = c(0, 0)),
      labels = label_number(accuracy = 0.1),
      sec.axis = dup_axis(labels = NULL, name = NULL)
    ) +
    theme_pub()

  if (nrow(sig_df) > 0) {
    p <- p +
      geom_text(
        data = sig_df,
        aes(
          x = x_pos,
          y = label_y,
          label = sig_label,
          color = I(color)
        ),
        inherit.aes = FALSE,
        size = SIG_TEXT_SIZE,
        family = FONT_FAMILY,
        vjust = 0
      )
  }

  attr(p, "joined_data") <- df_joined
  attr(p, "summary_data") <- summary_df
  p
}

test_high_half_greater <- function(df_joined, metric_col, metric_group) {
  if (is.null(df_joined) || nrow(df_joined) == 0) return(tibble())

  needed <- c("condition", "cv_half", metric_col)
  if (!all(needed %in% names(df_joined))) return(tibble())

  test_df <- df_joined %>%
    filter(
      cv_half %in% c("Low CV half", "High CV half"),
      !is.na(.data[[metric_col]])
    )

  if (nrow(test_df) == 0) return(tibble())

  out <- bind_rows(lapply(sort(unique(test_df$condition)), function(cond) {
    cond_df <- test_df %>% filter(condition == cond)
    high_vals <- cond_df %>%
      filter(cv_half == "High CV half") %>%
      select(all_of(metric_col)) %>%
      pull(1)
    low_vals <- cond_df %>%
      filter(cv_half == "Low CV half") %>%
      select(all_of(metric_col)) %>%
      pull(1)

    if (length(high_vals) < 2 || length(low_vals) < 2) return(tibble())

    wilcox_res <- suppressWarnings(
      wilcox.test(
        high_vals,
        low_vals,
        alternative = "greater",
        exact = FALSE
      )
    )
    t_res <- suppressWarnings(
      t.test(high_vals, low_vals, alternative = "greater")
    )

    tibble(
      metric_group = metric_group,
      metric = metric_col,
      condition = cond,
      n_high = length(high_vals),
      n_low = length(low_vals),
      mean_high = mean(high_vals, na.rm = TRUE),
      mean_low = mean(low_vals, na.rm = TRUE),
      diff_high_minus_low = mean_high - mean_low,
      wilcox_p_greater = wilcox_res$p.value,
      t_test_p_greater = t_res$p.value
    )
  }))

  if (nrow(out) == 0) return(out)

  out %>%
    mutate(
      wilcox_p_adj_bh = p.adjust(wilcox_p_greater, method = "BH"),
      t_test_p_adj_bh = p.adjust(t_test_p_greater, method = "BH")
    )
}

save_metric_set <- function(data_path, metric_list, cv_df, experiment,
                            fig_dir, stat_dir, suffix) {
  if (!file.exists(data_path)) return(invisible(NULL))

  df_metric <- read_csv(data_path, show_col_types = FALSE) %>%
    filter(experiment == !!experiment)

  test_results <- list()

  for (metric_name in names(metric_list)) {
    if (!metric_name %in% names(df_metric)) next

    p <- plot_cv_half_metric(
      df_metric = df_metric,
      cv_df = cv_df,
      metric_col = metric_name,
      experiment = experiment,
      ylim_fixed = metric_list[[metric_name]]$ylim
    )
    if (is.null(p)) next

    out_file <- file.path(fig_dir, sprintf("%s_cv_halves.pdf", metric_name))
    ggsave(
      out_file,
      p,
      width = W_PLOT,
      height = H_PLOT,
      units = "mm",
      device = cairo_pdf
    )

    joined <- attr(p, "joined_data")
    summary_df <- attr(p, "summary_data")
    if (!is.null(joined)) {
      write_csv(
        joined,
        file.path(stat_dir, sprintf("%s_%s_joined.csv", metric_name, suffix))
      )
    }
    if (!is.null(summary_df)) {
      write_csv(
        summary_df,
        file.path(stat_dir, sprintf("%s_%s_summary.csv", metric_name, suffix))
      )
    }

    high_tests <- test_high_half_greater(
      joined,
      metric_name,
      suffix
    )
    if (nrow(high_tests) > 0) {
      test_results[[length(test_results) + 1]] <- high_tests
      write_csv(
        high_tests,
        file.path(
          stat_dir,
          sprintf("%s_%s_high_gt_low_tests.csv", metric_name, suffix)
        )
      )
    }

    cat(sprintf("  -> %s_cv_halves.pdf\n", metric_name))
  }

  if (length(test_results) > 0) {
    write_csv(
      bind_rows(test_results),
      file.path(stat_dir, sprintf("%s_high_gt_low_tests_all.csv", suffix))
    )
  }
}

for (experiment in names(EXP_WINDOW_MAP)) {
  for (window in EXP_WINDOW_MAP[[experiment]]) {
    cat(sprintf("\n=== %s / %s ===\n", experiment, window))

    cv_df <- read_cv_halves(experiment, window)
    if (nrow(cv_df) == 0) {
      cat("  [skip] missing or empty fluctuation data\n")
      next
    }

    fig_dir <- file.path(FIG_DIR, experiment, window)
    stat_dir <- file.path(FIG_DIR, "stats", experiment, window)
    dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(stat_dir, showWarnings = FALSE, recursive = TRUE)

    write_csv(
      cv_df,
      file.path(stat_dir, "community_cv_halves_assignment.csv")
    )

    save_metric_set(
      file.path(PROC_DIR, "diversity", window, "alpha_diversity.csv"),
      ALPHA_METRICS, cv_df, experiment, fig_dir, stat_dir, "alpha"
    )

    save_metric_set(
      file.path(PROC_DIR, "diversity", window, "gamma_diversity.csv"),
      GAMMA_METRICS, cv_df, experiment, fig_dir, stat_dir, "gamma"
    )

    save_metric_set(
      file.path(PROC_DIR, "diversity", window, "mean_daily_diversity.csv"),
      MEAN_DAILY_METRICS, cv_df, experiment, fig_dir, stat_dir, "mean_daily"
    )

    save_metric_set(
      file.path(PROC_DIR, "diversity", window, "beta_within.csv"),
      list(mean_bc = list(ylim = c(0, 1))),
      cv_df, experiment, fig_dir, stat_dir, "beta_within"
    )

    save_metric_set(
      file.path(PROC_DIR, "diversity", window, "beta_across.csv"),
      list(bc_vs_reference = list(ylim = c(0, 1))),
      cv_df, experiment, fig_dir, stat_dir, "beta_across"
    )
  }
}

cat("\nDone.\n")
cat(sprintf("Output directory: %s\n", FIG_DIR))
cat("High CV half uses fluctuation orange; low CV half uses stable purple.\n")
