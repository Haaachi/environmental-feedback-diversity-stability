community_workspace <- function() {
  override <- Sys.getenv("COMMUNITY_WORKSPACE", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script)) stop("Run with Rscript or set COMMUNITY_WORKSPACE.")
  dirname(normalizePath(sub("^--file=", "", script[[1]]), mustWork = TRUE))
}

# ===========================================================
# plot_diversity_flucstable_clean.R
# Purpose:
# 1. Plot stable/fluctuating alpha, gamma and beta diversity curves without raw points.
# 2. Show colored stars for within-condition Wilcoxon comparisons.
# 3. Compare fluctuating communities across conditions using brackets and stars:
# Temperature: 30 versus 10 and 30 versus 40.
#    - mortality：W3 vs W1, W3 vs W5
# 4. Omit stable-group Kruskal-Wallis results from the figure.
# 5. Also compute Student's t tests alongside Mann-Whitney tests,
# and stable-group Kruskal-Wallis tests; save results to CSV.
#
# Implementation notes:
# Add ensure_collapsed_col().
# Ensure the collapsed column exists before using it.
# Exclude collapsed temperature W5 communities with richness == 0 first.
# Apply collapse exclusions before the species-instability criterion,
# to avoid incorrectly labeling collapsed communities as fluctuating.
# ===========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

BASE_DIR <- community_workspace()
PROC_DIR <- file.path(BASE_DIR, "processed")
FIG_DIR  <- file.path(BASE_DIR, "figures", "diversity_flucstable")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

COMMUNITY_CV_THRESHOLD <- 0.265

is_fluctuating_composite <- function(experiment, community_cv) {
  !is.na(community_cv) & community_cv >= COMMUNITY_CV_THRESHOLD
}
COLOR_STABLE   <- "#9B8EC4"
COLOR_FLUCT    <- "#F4A460"

FONT_FAMILY <- "Arial"
FONT_AX     <- 10
BORDER_SIZE <- 0.25
TICK_SIZE   <- 0.20
TICK_LEN    <- -1.0

W_PLOT <- 80
H_PLOT <- 70

LINE_WIDTH     <- 0.5
ERR_WIDTH      <- 0.12
ERR_LW         <- 0.35
MEAN_SIZE      <- 1.8
SIG_TEXT_SIZE  <- 3.0

# ============================================================
# Allowed experiment-window combinations
# ============================================================

EXP_WINDOW_MAP <- list(
  mortality   = c("early"),
  temperature = c("full")
)

is_valid_combo <- function(experiment, window) {
  allowed <- EXP_WINDOW_MAP[[experiment]]
  !is.null(allowed) && window %in% allowed
}

# ============================================================
# Helper for a missing collapsed column
# ============================================================

ensure_collapsed_col <- function(df) {
  if (!"collapsed" %in% colnames(df)) {
    df <- df %>%
      mutate(collapsed = FALSE)
  } else {
    df <- df %>%
      mutate(collapsed = if_else(is.na(collapsed), FALSE, collapsed))
  }
  df
}

# ============================================================
# X-axis labels
# ============================================================

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
  richness          = list(ylim = NULL),
  shannon           = list(ylim = c(0, 1.6)),
  effective_shannon = list(ylim = NULL),
  simpson           = list(ylim = c(0, 1)),
  evenness          = list(ylim = c(0, 1)),
  dominance         = list(ylim = c(0, 1)),
  survival_fraction = list(ylim = c(0, 1))
)

GAMMA_METRICS <- list(
  gamma_richness          = list(ylim = NULL),
  gamma_shannon           = list(ylim = NULL),
  gamma_effective_shannon = list(ylim = NULL),
  gamma_simpson           = list(ylim = c(0, 1)),
  gamma_evenness          = list(ylim = c(0, 1))
)

MEAN_DAILY_METRICS <- list(
  mean_daily_richness          = list(ylim = NULL),
  mean_daily_shannon           = list(ylim = NULL),
  mean_daily_effective_shannon = list(ylim = NULL),
  mean_daily_simpson           = list(ylim = c(0, 1)),
  mean_daily_evenness          = list(ylim = c(0, 1)),
  mean_daily_survival_fraction = list(ylim = c(0, 1))
)

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
        size = FONT_AX,
        color = "black",
        margin = margin(t = 3)
      ),
      axis.text.y  = element_text(
        size = FONT_AX,
        color = "black",
        margin = margin(r = 3)
      ),
      axis.title         = element_blank(),
      legend.position    = "none",
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white", color = NA),
      plot.margin        = margin(3, 3, 3, 3, "mm")
    )
}

# ============================================================
# Plotting function
# ============================================================

plot_flucstable_metric <- function(df_metric, cv_df,
                                   metric_col, experiment,
                                   ylim_fixed = NULL,
                                   add_fluct_bracket = TRUE) {
  
  x_labs <- X_LABELS[[experiment]]
  
  # Ensure the collapsed column exists.
  cv_df <- ensure_collapsed_col(cv_df)
  
  join_keys <- intersect(
    c("experiment", "condition", "community", "replica"),
    colnames(df_metric)
  )
  
  cv_join <- cv_df %>%
    select(
      experiment,
      condition,
      community,
      replica,
      community_cv,
      collapsed
    ) %>%
    ensure_collapsed_col()
  
  df_joined <- df_metric %>%
    left_join(cv_join, by = join_keys) %>%
    ensure_collapsed_col() %>%
    mutate(
      x_pos = as.integer(sub("W", "", condition)),
      state = case_when(
        is.na(community_cv)                          ~ "no_data",
        collapsed == TRUE                             ~ "stable",
        
        # Mandatory exclusion:
        # For temperature W5, retain communities 3, 7 and 12 only.
        experiment == "temperature" &
          condition == "W5" &
          !community %in% c(3, 7, 12) ~ "no_data",
        
        experiment == "temperature" &
          condition == "W5" &
          collapsed == TRUE ~ "no_data",
        
        is_fluctuating_composite(experiment,
                                 community_cv)        ~ "fluctuating",
        TRUE                                          ~ "stable"
      )
    ) %>%
    filter(state != "no_data")
  
  if (nrow(df_joined) == 0) return(NULL)
  if (sum(!is.na(df_joined[[metric_col]])) == 0) return(NULL)
  
  summary_df <- df_joined %>%
    filter(!is.na(.data[[metric_col]])) %>%
    group_by(condition, x_pos, state) %>%
    summarise(
      mean_val = mean(.data[[metric_col]], na.rm = TRUE),
      sem_val  = sd(.data[[metric_col]], na.rm = TRUE) /
        sqrt(sum(!is.na(.data[[metric_col]]))),
      .groups  = "drop"
    ) %>%
    mutate(
      ymin = mean_val - sem_val,
      ymax = mean_val + sem_val
    )
  
  data_max <- max(df_joined[[metric_col]], na.rm = TRUE)
  
  if (!is.null(ylim_fixed)) {
    y_lower <- ylim_fixed[1]
    y_upper <- ylim_fixed[2]
    y_breaks <- pretty(c(y_lower, y_upper), n = 4)
    y_breaks <- y_breaks[y_breaks >= y_lower & y_breaks <= y_upper]
  } else {
    y_lower <- 0
    y_upper <- data_max * 1.05
    y_breaks <- pretty(c(y_lower, y_upper), n = 4)
    y_upper <- max(y_breaks)
  }
  
  # ------------------------------------------------------------
  # Within-condition stable/fluctuating Wilcoxon comparisons
  # ------------------------------------------------------------
  
  star_data <- data.frame()
  conditions <- unique(df_joined$condition)
  
  for (cond in conditions) {
    vals_stable <- df_joined %>%
      filter(
        condition == cond,
        state == "stable",
        !is.na(.data[[metric_col]])
      ) %>%
      pull(metric_col)
    
    vals_fluct <- df_joined %>%
      filter(
        condition == cond,
        state == "fluctuating",
        !is.na(.data[[metric_col]])
      ) %>%
      pull(metric_col)
    
    if (length(vals_stable) >= 2 && length(vals_fluct) >= 2) {
      test_res <- wilcox.test(vals_stable, vals_fluct, exact = FALSE)
      p_val <- test_res$p.value
      
      mean_stable <- mean(vals_stable, na.rm = TRUE)
      mean_fluct  <- mean(vals_fluct, na.rm = TRUE)
      higher_group <- ifelse(mean_stable > mean_fluct, "stable", "fluctuating")
      
      sig_star <- case_when(
        p_val < 0.001 ~ "***",
        p_val < 0.01  ~ "**",
        p_val < 0.05  ~ "*",
        TRUE          ~ NA_character_
      )
      
      if (!is.na(sig_star)) {
        star_data <- bind_rows(
          star_data,
          tibble(
            condition = cond,
            x_pos = as.integer(sub("W", "", cond)),
            star = sig_star,
            color = ifelse(
              higher_group == "stable",
              COLOR_STABLE,
              COLOR_FLUCT
            )
          )
        )
      }
    }
  }
  
  if (nrow(star_data) > 0) {
    max_per_cond <- summary_df %>%
      group_by(condition, x_pos) %>%
      summarise(
        max_upper = max(ymax, na.rm = TRUE),
        .groups = "drop"
      )
    
    y_range <- y_upper - y_lower
    label_y_offset <- y_range * 0.08
    
    label_data <- max_per_cond %>%
      inner_join(star_data, by = c("condition", "x_pos")) %>%
      mutate(label_y = max_upper + label_y_offset) %>%
      mutate(label_y = pmin(label_y, y_upper + y_range * 0.02))
  } else {
    label_data <- data.frame()
  }
  
  # ------------------------------------------------------------
  # Base plot
  # ------------------------------------------------------------
  
  p <- ggplot() +
    geom_line(
      data = summary_df %>% filter(state == "stable"),
      aes(x = x_pos, y = mean_val, color = state),
      linewidth = LINE_WIDTH
    ) +
    geom_errorbar(
      data = summary_df %>% filter(state == "stable"),
      aes(x = x_pos, ymin = ymin, ymax = ymax, color = state),
      width = ERR_WIDTH,
      linewidth = ERR_LW
    ) +
    geom_point(
      data = summary_df %>% filter(state == "stable"),
      aes(x = x_pos, y = mean_val, color = state),
      shape = 16,
      size = MEAN_SIZE
    ) +
    geom_line(
      data = summary_df %>% filter(state == "fluctuating"),
      aes(x = x_pos, y = mean_val, color = state),
      linewidth = LINE_WIDTH
    ) +
    geom_errorbar(
      data = summary_df %>% filter(state == "fluctuating"),
      aes(x = x_pos, ymin = ymin, ymax = ymax, color = state),
      width = ERR_WIDTH,
      linewidth = ERR_LW
    ) +
    geom_point(
      data = summary_df %>% filter(state == "fluctuating"),
      aes(x = x_pos, y = mean_val, color = state),
      shape = 16,
      size = MEAN_SIZE
    ) +
    scale_x_continuous(
      breaks = 1:5,
      labels = x_labs,
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
    scale_color_manual(
      values = c(
        "stable" = COLOR_STABLE,
        "fluctuating" = COLOR_FLUCT
      )
    ) +
    theme_pub()
  
  if (nrow(label_data) > 0) {
    p <- p +
      geom_text(
        data = label_data,
        aes(
          x = x_pos,
          y = label_y,
          label = star,
          color = I(color)
        ),
        size = SIG_TEXT_SIZE,
        family = FONT_FAMILY,
        vjust = 0
      )
  }
  
  # ------------------------------------------------------------
  # Between-condition brackets for fluctuating communities
  # ------------------------------------------------------------
  
  if (add_fluct_bracket) {
    
    if (experiment == "temperature") {
      fluct_data <- df_joined %>%
        filter(
          state == "fluctuating",
          condition %in% c("W1", "W3", "W4")
        )
      
      if (n_distinct(fluct_data$condition) == 3) {
        x1 <- 1
        x_mid <- 3
        x2 <- 4
        
        y_vals_10 <- fluct_data[[metric_col]][fluct_data$condition == "W1"]
        y_vals_30 <- fluct_data[[metric_col]][fluct_data$condition == "W3"]
        y_vals_40 <- fluct_data[[metric_col]][fluct_data$condition == "W4"]
        
        if (
          length(y_vals_10) >= 2 &&
          length(y_vals_30) >= 2 &&
          length(y_vals_40) >= 2
        ) {
          test_mid_low <- wilcox.test(y_vals_30, y_vals_10, exact = FALSE)
          test_mid_high <- wilcox.test(y_vals_30, y_vals_40, exact = FALSE)
          
          star1 <- case_when(
            test_mid_low$p.value < 0.001 ~ "***",
            test_mid_low$p.value < 0.01  ~ "**",
            test_mid_low$p.value < 0.05  ~ "*",
            TRUE ~ "ns"
          )
          
          star2 <- case_when(
            test_mid_high$p.value < 0.001 ~ "***",
            test_mid_high$p.value < 0.01  ~ "**",
            test_mid_high$p.value < 0.05  ~ "*",
            TRUE ~ "ns"
          )
          
          y_range <- y_upper - y_lower
          bracket_y <- y_upper + y_range * 0.08
          bracket_y2 <- bracket_y + y_range * 0.03
          
          p <- p +
            geom_segment(
              aes(x = x1, xend = x_mid, y = bracket_y, yend = bracket_y),
              color = COLOR_FLUCT,
              linewidth = 0.3
            ) +
            geom_segment(
              aes(
                x = x1,
                xend = x1,
                y = bracket_y,
                yend = bracket_y - y_range * 0.01
              ),
              color = COLOR_FLUCT,
              linewidth = 0.3
            ) +
            geom_segment(
              aes(
                x = x_mid,
                xend = x_mid,
                y = bracket_y,
                yend = bracket_y - y_range * 0.01
              ),
              color = COLOR_FLUCT,
              linewidth = 0.3
            ) +
            annotate(
              "text",
              x = (x1 + x_mid) / 2,
              y = bracket_y2,
              label = star1,
              color = COLOR_FLUCT,
              size = SIG_TEXT_SIZE
            ) +
            geom_segment(
              aes(x = x_mid, xend = x2, y = bracket_y, yend = bracket_y),
              color = COLOR_FLUCT,
              linewidth = 0.3
            ) +
            geom_segment(
              aes(
                x = x_mid,
                xend = x_mid,
                y = bracket_y,
                yend = bracket_y - y_range * 0.01
              ),
              color = COLOR_FLUCT,
              linewidth = 0.3
            ) +
            geom_segment(
              aes(
                x = x2,
                xend = x2,
                y = bracket_y,
                yend = bracket_y - y_range * 0.01
              ),
              color = COLOR_FLUCT,
              linewidth = 0.3
            ) +
            annotate(
              "text",
              x = (x_mid + x2) / 2,
              y = bracket_y2,
              label = star2,
              color = COLOR_FLUCT,
              size = SIG_TEXT_SIZE
            )
        }
      }
      
    } else if (experiment == "mortality") {
      fluct_data <- df_joined %>%
        filter(
          state == "fluctuating",
          condition %in% c("W1", "W3", "W5")
        )
      
      if (n_distinct(fluct_data$condition) == 3) {
        x_low <- 1
        x_mid <- 3
        x_high <- 5
        
        y_vals_low <- fluct_data[[metric_col]][fluct_data$condition == "W1"]
        y_vals_mid <- fluct_data[[metric_col]][fluct_data$condition == "W3"]
        y_vals_high <- fluct_data[[metric_col]][fluct_data$condition == "W5"]
        
        if (
          length(y_vals_low) >= 2 &&
          length(y_vals_mid) >= 2 &&
          length(y_vals_high) >= 2
        ) {
          test_mid_low <- wilcox.test(y_vals_mid, y_vals_low, exact = FALSE)
          test_mid_high <- wilcox.test(y_vals_mid, y_vals_high, exact = FALSE)
          
          star1 <- case_when(
            test_mid_low$p.value < 0.001 ~ "***",
            test_mid_low$p.value < 0.01  ~ "**",
            test_mid_low$p.value < 0.05  ~ "*",
            TRUE ~ "ns"
          )
          
          star2 <- case_when(
            test_mid_high$p.value < 0.001 ~ "***",
            test_mid_high$p.value < 0.01  ~ "**",
            test_mid_high$p.value < 0.05  ~ "*",
            TRUE ~ "ns"
          )
          
          y_range <- y_upper - y_lower
          bracket_y <- y_upper + y_range * 0.08
          bracket_y2 <- bracket_y + y_range * 0.03
          
          p <- p +
            geom_segment(
              aes(x = x_low, xend = x_mid, y = bracket_y, yend = bracket_y),
              color = COLOR_FLUCT,
              linewidth = 0.3
            ) +
            geom_segment(
              aes(
                x = x_low,
                xend = x_low,
                y = bracket_y,
                yend = bracket_y - y_range * 0.01
              ),
              color = COLOR_FLUCT,
              linewidth = 0.3
            ) +
            geom_segment(
              aes(
                x = x_mid,
                xend = x_mid,
                y = bracket_y,
                yend = bracket_y - y_range * 0.01
              ),
              color = COLOR_FLUCT,
              linewidth = 0.3
            ) +
            annotate(
              "text",
              x = (x_low + x_mid) / 2,
              y = bracket_y2,
              label = star1,
              color = COLOR_FLUCT,
              size = SIG_TEXT_SIZE
            ) +
            geom_segment(
              aes(x = x_mid, xend = x_high, y = bracket_y, yend = bracket_y),
              color = COLOR_FLUCT,
              linewidth = 0.3
            ) +
            geom_segment(
              aes(
                x = x_mid,
                xend = x_mid,
                y = bracket_y,
                yend = bracket_y - y_range * 0.01
              ),
              color = COLOR_FLUCT,
              linewidth = 0.3
            ) +
            geom_segment(
              aes(
                x = x_high,
                xend = x_high,
                y = bracket_y,
                yend = bracket_y - y_range * 0.01
              ),
              color = COLOR_FLUCT,
              linewidth = 0.3
            ) +
            annotate(
              "text",
              x = (x_mid + x_high) / 2,
              y = bracket_y2,
              label = star2,
              color = COLOR_FLUCT,
              size = SIG_TEXT_SIZE
            )
        }
      }
    }
  }
  
  attr(p, "joined_data") <- df_joined
  
  return(p)
}

# ============================================================
# Additional statistical calculations
# ============================================================

compute_and_save_stats <- function(df_joined, metric_col, experiment, stat_dir, metric_name) {
  
  # 1. Within each condition: stable versus fluctuating.
  conditions <- unique(df_joined$condition)
  all_res <- list()
  
  for (cond in conditions) {
    vals_stable <- df_joined %>%
      filter(
        condition == cond,
        state == "stable",
        !is.na(.data[[metric_col]])
      ) %>%
      pull(metric_col)
    
    vals_fluct <- df_joined %>%
      filter(
        condition == cond,
        state == "fluctuating",
        !is.na(.data[[metric_col]])
      ) %>%
      pull(metric_col)
    
    if (length(vals_stable) >= 2 && length(vals_fluct) >= 2) {
      w_res <- wilcox.test(vals_stable, vals_fluct, exact = FALSE)
      t_res <- t.test(vals_stable, vals_fluct)
      
      all_res[[cond]] <- tibble(
        condition = cond,
        n_stable = length(vals_stable),
        n_fluct = length(vals_fluct),
        mean_stable = mean(vals_stable),
        mean_fluct = mean(vals_fluct),
        wilcoxon_p = w_res$p.value,
        t_p = t_res$p.value,
        t_statistic = as.numeric(t_res$statistic),
        t_df = as.numeric(t_res$parameter)
      )
    }
  }
  
  if (length(all_res) > 0) {
    res_df <- bind_rows(all_res)
    write_csv(
      res_df,
      file.path(stat_dir, paste0(metric_name, "_stable_vs_fluct.csv"))
    )
  }
  
  # 2. Selected between-condition comparisons for fluctuating communities.
  if (experiment == "temperature") {
    fluct_data <- df_joined %>%
      filter(
        state == "fluctuating",
        condition %in% c("W1", "W3", "W4")
      )
    
    if (n_distinct(fluct_data$condition) == 3) {
      y_10 <- fluct_data[[metric_col]][fluct_data$condition == "W1"]
      y_30 <- fluct_data[[metric_col]][fluct_data$condition == "W3"]
      y_40 <- fluct_data[[metric_col]][fluct_data$condition == "W4"]
      
      if (length(y_10) >= 2 && length(y_30) >= 2 && length(y_40) >= 2) {
        w_30_10 <- wilcox.test(y_30, y_10, exact = FALSE)
        t_30_10 <- t.test(y_30, y_10)
        w_30_40 <- wilcox.test(y_30, y_40, exact = FALSE)
        t_30_40 <- t.test(y_30, y_40)
        
        fluct_comp <- tibble(
          comparison = c("30°C_vs_10°C", "30°C_vs_40°C"),
          wilcoxon_p = c(w_30_10$p.value, w_30_40$p.value),
          t_p = c(t_30_10$p.value, t_30_40$p.value),
          t_statistic = c(
            as.numeric(t_30_10$statistic),
            as.numeric(t_30_40$statistic)
          ),
          t_df = c(
            as.numeric(t_30_10$parameter),
            as.numeric(t_30_40$parameter)
          )
        )
        
        write_csv(
          fluct_comp,
          file.path(stat_dir, paste0(metric_name, "_fluct_selected.csv"))
        )
      }
    }
    
  } else if (experiment == "mortality") {
    fluct_data <- df_joined %>%
      filter(
        state == "fluctuating",
        condition %in% c("W1", "W3", "W5")
      )
    
    if (n_distinct(fluct_data$condition) == 3) {
      y_low <- fluct_data[[metric_col]][fluct_data$condition == "W1"]
      y_mid <- fluct_data[[metric_col]][fluct_data$condition == "W3"]
      y_high <- fluct_data[[metric_col]][fluct_data$condition == "W5"]
      
      if (length(y_low) >= 2 && length(y_mid) >= 2 && length(y_high) >= 2) {
        w_mid_low <- wilcox.test(y_mid, y_low, exact = FALSE)
        t_mid_low <- t.test(y_mid, y_low)
        w_mid_high <- wilcox.test(y_mid, y_high, exact = FALSE)
        t_mid_high <- t.test(y_mid, y_high)
        
        fluct_comp <- tibble(
          comparison = c("W3_vs_W1", "W3_vs_W5"),
          wilcoxon_p = c(w_mid_low$p.value, w_mid_high$p.value),
          t_p = c(t_mid_low$p.value, t_mid_high$p.value),
          t_statistic = c(
            as.numeric(t_mid_low$statistic),
            as.numeric(t_mid_high$statistic)
          ),
          t_df = c(
            as.numeric(t_mid_low$parameter),
            as.numeric(t_mid_high$parameter)
          )
        )
        
        write_csv(
          fluct_comp,
          file.path(stat_dir, paste0(metric_name, "_fluct_selected.csv"))
        )
      }
    }
  }
  
  # 3. Kruskal-Wallis tests for stable communities.
  if (experiment == "temperature") {
    stable_data <- df_joined %>%
      filter(
        state == "stable",
        condition %in% c("W1", "W2", "W3", "W4")
      )
    
    if (n_distinct(stable_data$condition) >= 3) {
      kw_res <- kruskal.test(stable_data[[metric_col]] ~ stable_data$condition)
      
      kw_df <- tibble(
        metric = metric_name,
        p_value = kw_res$p.value,
        method = "Kruskal-Wallis"
      )
      
      write_csv(
        kw_df,
        file.path(stat_dir, paste0(metric_name, "_stable_kruskal.csv"))
      )
    }
    
  } else if (experiment == "mortality") {
    stable_data <- df_joined %>%
      filter(
        state == "stable",
        condition %in% c("W1", "W2", "W3", "W4", "W5")
      )
    
    if (n_distinct(stable_data$condition) >= 3) {
      kw_res <- kruskal.test(stable_data[[metric_col]] ~ stable_data$condition)
      
      kw_df <- tibble(
        metric = metric_name,
        p_value = kw_res$p.value,
        method = "Kruskal-Wallis"
      )
      
      write_csv(
        kw_df,
        file.path(stat_dir, paste0(metric_name, "_stable_kruskal.csv"))
      )
    }
  }
}

# ============================================================
# Main loop
# ============================================================

for (experiment in c("mortality", "temperature")) {
  for (window in c("full", "early")) {
    
    if (!is_valid_combo(experiment, window)) {
      cat(sprintf(
        "\n[skip] %s × %s (unsupported combination; see EXP_WINDOW_MAP)\n",
        experiment,
        window
      ))
      next
    }
    
    cat(sprintf("\n=== %s / %s ===\n", experiment, window))
    
    fluc_path <- file.path(
      PROC_DIR,
      "fluctuations",
      window,
      "community_level.csv"
    )
    
    if (!file.exists(fluc_path)) {
      cat("  [skip] fluctuations data missing\n")
      next
    }
    
    cv_df <- read_csv(fluc_path, show_col_types = FALSE) %>%
      filter(experiment == !!experiment) %>%
      select(experiment, condition, community, replica,
             community_cv) %>%
      ensure_collapsed_col()
    
    # ------------------------------------------------------------
    # Infer collapse flags.
    # Treat richness == 0 as collapsed.
    # If the alpha file is missing, set collapsed = FALSE for all rows.
    # ------------------------------------------------------------
    
    alpha_path_for_collapse <- file.path(
      PROC_DIR,
      "diversity",
      window,
      "alpha_diversity.csv"
    )
    
    if (file.exists(alpha_path_for_collapse)) {
      collapse_df <- read_csv(alpha_path_for_collapse, show_col_types = FALSE) %>%
        filter(experiment == !!experiment) %>%
        mutate(collapsed = !is.na(richness) & richness == 0) %>%
        select(experiment, condition, community, replica, collapsed) %>%
        ensure_collapsed_col()
      
      cv_df <- cv_df %>%
        select(-collapsed) %>%
        left_join(
          collapse_df,
          by = c("experiment", "condition", "community", "replica")
        ) %>%
        ensure_collapsed_col()
    } else {
      cv_df <- cv_df %>%
        ensure_collapsed_col()
    }
    
    fig_dir <- file.path(FIG_DIR, experiment, window)
    stat_dir <- file.path(FIG_DIR, "stats", experiment, window)
    
    dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(stat_dir, showWarnings = FALSE, recursive = TRUE)
    
    # ------------------------------------------------------------
    # Alpha diversity
    # ------------------------------------------------------------
    
    alpha_path <- file.path(
      PROC_DIR,
      "diversity",
      window,
      "alpha_diversity.csv"
    )
    
    if (file.exists(alpha_path)) {
      alpha_df <- read_csv(alpha_path, show_col_types = FALSE) %>%
        filter(experiment == !!experiment)
      
      for (metric in names(ALPHA_METRICS)) {
        if (!metric %in% colnames(alpha_df)) next
        
        current_ylim <- ALPHA_METRICS[[metric]]$ylim
        
        if (experiment == "temperature" && metric == "survival_fraction") {
          current_ylim <- c(0, 0.8)
        }
        
        p <- plot_flucstable_metric(
          alpha_df,
          cv_df,
          metric,
          experiment,
          current_ylim,
          add_fluct_bracket = TRUE
        )
        
        if (is.null(p)) next
        
        out_file <- file.path(fig_dir, sprintf("%s_flucstable.pdf", metric))
        
        ggsave(
          out_file,
          p,
          width = W_PLOT,
          height = H_PLOT,
          units = "mm",
          device = cairo_pdf
        )
        
        cat(sprintf("  -> Figure: %s_flucstable.pdf\n", metric))
        
        df_joined <- attr(p, "joined_data")
        
        if (!is.null(df_joined)) {
          compute_and_save_stats(
            df_joined,
            metric,
            experiment,
            stat_dir,
            metric
          )
        }
      }
    }
    
    # ------------------------------------------------------------
    # Gamma diversity
    # ------------------------------------------------------------
    
    gamma_path <- file.path(
      PROC_DIR,
      "diversity",
      window,
      "gamma_diversity.csv"
    )
    
    if (file.exists(gamma_path)) {
      gamma_df <- read_csv(gamma_path, show_col_types = FALSE) %>%
        filter(experiment == !!experiment)
      
      for (metric in names(GAMMA_METRICS)) {
        if (!metric %in% colnames(gamma_df)) next
        
        p <- plot_flucstable_metric(
          gamma_df,
          cv_df,
          metric,
          experiment,
          GAMMA_METRICS[[metric]]$ylim,
          add_fluct_bracket = TRUE
        )
        
        if (is.null(p)) next
        
        out_file <- file.path(fig_dir, sprintf("%s_flucstable.pdf", metric))
        
        ggsave(
          out_file,
          p,
          width = W_PLOT,
          height = H_PLOT,
          units = "mm",
          device = cairo_pdf
        )
        
        cat(sprintf("  -> Figure: %s_flucstable.pdf\n", metric))
        
        df_joined <- attr(p, "joined_data")
        
        if (!is.null(df_joined)) {
          compute_and_save_stats(
            df_joined,
            metric,
            experiment,
            stat_dir,
            metric
          )
        }
      }
    }

    mean_daily_path <- file.path(
      PROC_DIR,
      "diversity",
      window,
      "mean_daily_diversity.csv"
    )

    if (file.exists(mean_daily_path)) {
      mean_daily_df <- read_csv(mean_daily_path, show_col_types = FALSE) %>%
        filter(experiment == !!experiment)

      for (metric in names(MEAN_DAILY_METRICS)) {
        if (!metric %in% colnames(mean_daily_df)) next

        p <- plot_flucstable_metric(
          mean_daily_df,
          cv_df,
          metric,
          experiment,
          MEAN_DAILY_METRICS[[metric]]$ylim,
          add_fluct_bracket = TRUE
        )

        if (is.null(p)) next

        out_file <- file.path(fig_dir, sprintf("%s_flucstable.pdf", metric))

        ggsave(
          out_file,
          p,
          width = W_PLOT,
          height = H_PLOT,
          units = "mm",
          device = cairo_pdf
        )

        cat(sprintf("  -> Figure: %s_flucstable.pdf\n", metric))

        df_joined <- attr(p, "joined_data")

        if (!is.null(df_joined)) {
          compute_and_save_stats(
            df_joined,
            metric,
            experiment,
            stat_dir,
            metric
          )
        }
      }
    }
    
    # ------------------------------------------------------------
    # Beta within full only
    # ------------------------------------------------------------
    
    if (window == "full") {
      bw_path <- file.path(
        PROC_DIR,
        "diversity",
        window,
        "beta_within.csv"
      )
      
      if (file.exists(bw_path)) {
        bw_df <- read_csv(bw_path, show_col_types = FALSE) %>%
          filter(experiment == !!experiment)
        
        p <- plot_flucstable_metric(
          bw_df,
          cv_df,
          "mean_bc",
          experiment,
          c(0, 1),
          add_fluct_bracket = TRUE
        )
        
        if (!is.null(p)) {
          out_file <- file.path(fig_dir, "beta_within_flucstable.pdf")
          
          ggsave(
            out_file,
            p,
            width = W_PLOT,
            height = H_PLOT,
            units = "mm",
            device = cairo_pdf
          )
          
          cat("  -> Figure: beta_within_flucstable.pdf\n")
          
          df_joined <- attr(p, "joined_data")
          
          if (!is.null(df_joined)) {
            compute_and_save_stats(
              df_joined,
              "mean_bc",
              experiment,
              stat_dir,
              "beta_within"
            )
          }
        }
      }
    }
    
    # ------------------------------------------------------------
    # Beta across
    # ------------------------------------------------------------
    
    ba_path <- file.path(
      PROC_DIR,
      "diversity",
      window,
      "beta_across.csv"
    )
    
    if (file.exists(ba_path)) {
      ba_df <- read_csv(ba_path, show_col_types = FALSE) %>%
        filter(experiment == !!experiment)
      
      p <- plot_flucstable_metric(
        ba_df,
        cv_df,
        "bc_vs_reference",
        experiment,
        c(0, 1),
        add_fluct_bracket = TRUE
      )
      
      if (!is.null(p)) {
        out_file <- file.path(fig_dir, "beta_across_flucstable.pdf")
        
        ggsave(
          out_file,
          p,
          width = W_PLOT,
          height = H_PLOT,
          units = "mm",
          device = cairo_pdf
        )
        
        cat("  -> Figure: beta_across_flucstable.pdf\n")
        
        df_joined <- attr(p, "joined_data")
        
        if (!is.null(df_joined)) {
          compute_and_save_stats(
            df_joined,
            "bc_vs_reference",
            experiment,
            stat_dir,
            "beta_across"
          )
        }
      }
    }
  }
}

cat("\n===== Done!=====\n")
cat("Figures saved; annotations show within-condition stable/fluctuating stars and between-condition fluctuating-group brackets.\n")
cat("All statistical results (Wilcoxon, t-test, Kruskal-Wallis) saved in the stats subdirectory.\n")
cat("\nExperiment - Window constraints:\n")

for (exp_name in names(EXP_WINDOW_MAP)) {
  cat(sprintf(
    "  %-12s → %s\n",
    exp_name,
    paste(EXP_WINDOW_MAP[[exp_name]], collapse = ", ")
  ))
}
