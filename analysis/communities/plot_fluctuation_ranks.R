community_workspace <- function() {
  override <- Sys.getenv("COMMUNITY_WORKSPACE", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script)) stop("Run with Rscript or set COMMUNITY_WORKSPACE.")
  dirname(normalizePath(sub("^--file=", "", script[[1]]), mustWork = TRUE))
}

# ===========================================================
# plot_fluctuation_rank.R
# Original workspace: /home/hachi/Tem_mortality_workspace
#
# Rank plots of instability metrics:
#   1. community_cv     = Σσᵢ_abs / mean_total_abs
#   2. total_biomass_cv = std(OD_t) / mean(OD_t)
# 3. temporal_bc_mean: mean Bray-Curtis distance over all day pairs in the window.
#   4. weighted_log_sd  = Σ wᵢ sd(log(absᵢ + eps))
#
# Default version: R1 only.
# Additional temperature version: plot R1/R2/R3 together (allrep).
#
# Export rank-only plots; interpret operational guides from the plotted values.
# Use a shared Y upper limit per metric across full/early, conditions and replicates.
#
# Output: figures/fluctuations/rank/
#   <metric>/
#     full/
#       mortality_W1.pdf ~ mortality_W5.pdf
#       mortality_combined.pdf
#       temperature_W1.pdf ~ temperature_W5.pdf
#       temperature_combined.pdf
#       temperature_W1_allrep.pdf ~ temperature_W5_allrep.pdf
#       temperature_combined_allrep.pdf
#     early/
# Same structure as above.
# ===========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

BASE_DIR <- community_workspace()
PROC_DIR <- file.path(BASE_DIR, "processed")
FIG_DIR  <- file.path(BASE_DIR, "figures", "fluctuations", "rank")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Operational thresholds; adjust here.
# ============================================================
RANK_METRICS <- list(
  sum_abs_std       = NA_real_,
  sum_rel_std       = NA_real_,
  temporal_bc_mean  = NA_real_,
  temporal_bc_max   = 0.20,
  weighted_log_sd   = NA_real_,
  community_cv      = 0.25,
  community_cv_rel  = 0.25,
  total_biomass_std = NA_real_,
  total_biomass_cv  = NA_real_
)

threshold_for <- function(metric, experiment) {
  threshold <- RANK_METRICS[[metric]]
  if (!is.null(threshold) && is.finite(threshold)) {
    return(threshold)
  }
  NA_real_
}

# Colors
COL_FLUCTUATION <- "#F4A460"
COL_STABLE      <- "#9B8EC4"
COL_RANK_ONLY   <- "#6F6F6F"

# Plot parameters
FONT_FAMILY  <- "Arial"
FONT_AX      <- 14
FONT_LEGEND  <- 18
BORDER_SIZE  <- 0.25
TICK_SIZE    <- 0.20
TICK_LEN     <- -0.5

W_PLOT <- 90
H_PLOT <- 80

DOT_SIZE      <- 1.6
DOT_ALPHA     <- 1.0
THRESHOLD_LW  <- 0.4
THRESHOLD_LT  <- "dashed"
THRESHOLD_COL <- "black"

# Publication plot theme
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
      axis.text.x  = element_text(size   = FONT_AX,
                                  color  = "black",
                                  margin = margin(t = 3)),
      axis.text.y  = element_text(size   = FONT_AX,
                                  color  = "black",
                                  margin = margin(r = 3)),
      axis.title   = element_blank(),
      legend.position  = "none",
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin      = margin(3, 3, 3, 3, "mm")
    )
}

# Core plotting function
plot_rank <- function(df_input, metric, threshold,
                      y_upper, y_breaks, show_legend = FALSE) {
  has_threshold <- is.finite(threshold)
  
  df_plot <- df_input %>%
    filter(!is.na(.data[[metric]])) %>%
    arrange(desc(.data[[metric]])) %>%
    mutate(
      rank        = row_number(),
      value       = .data[[metric]],
      fluctuating = if (has_threshold) value >= threshold else FALSE
    )
  
  if (nrow(df_plot) == 0) return(NULL)
  
  n        <- nrow(df_plot)
  x_breaks <- unique(c(1, pretty(c(1, n), n = 4)))
  x_breaks <- x_breaks[x_breaks >= 1 & x_breaks <= n]
  
  p <- ggplot(df_plot, aes(x = rank, y = value))

  if (has_threshold) {
    p <- p +
      geom_hline(yintercept = threshold,
                 linetype   = THRESHOLD_LT,
                 color      = THRESHOLD_COL,
                 linewidth  = THRESHOLD_LW)
  }

  if (has_threshold) {
    p <- p +
      geom_point(aes(color = fluctuating),
                 size  = DOT_SIZE,
                 alpha = DOT_ALPHA) +
      scale_color_manual(
        values = c("TRUE"  = COL_FLUCTUATION,
                   "FALSE" = COL_STABLE)
      )
  } else {
    p <- p +
      geom_point(size  = DOT_SIZE,
                 alpha = DOT_ALPHA,
                 color = COL_RANK_ONLY)
  }

  p <- p +
    scale_x_continuous(
      limits   = c(0, n + 1),
      breaks   = x_breaks,
      expand   = expansion(mult = c(0.02, 0.02)),
      sec.axis = dup_axis(labels = NULL, name = NULL)
    ) +
    scale_y_continuous(
      limits   = c(0, y_upper),
      breaks   = y_breaks,
      expand   = expansion(mult = c(0, 0)),
      sec.axis = dup_axis(labels = NULL, name = NULL)
    ) +
    labs(x = NULL, y = NULL) +
    theme_pub()
  
  if (show_legend && has_threshold) {
    p <- p +
      annotate("text",
               x      = n * 0.05,
               y      = y_upper * 0.97,
               label  = "Fluctuation",
               color  = COL_FLUCTUATION,
               size   = FONT_LEGEND / .pt,
               hjust  = 0, vjust = 1,
               family = FONT_FAMILY) +
      annotate("text",
               x      = n * 0.05,
               y      = y_upper * 0.97 - y_upper * 0.09,
               label  = "Stable",
               color  = COL_STABLE,
               size   = FONT_LEGEND / .pt,
               hjust  = 0, vjust = 1,
               family = FONT_FAMILY)
  }
  
  p
}

# Main loop: metrics x windows x experiments x conditions
for (metric in names(RANK_METRICS)) {
  metric_threshold_values <- c(
    threshold_for(metric, "mortality"),
    threshold_for(metric, "temperature")
  )
  has_any_threshold <- any(is.finite(metric_threshold_values))
  
  if (has_any_threshold) {
    cat(sprintf("\n========== %s (experiment-specific threshold) ==========\n",
                metric))
  } else {
    cat(sprintf("\n========== %s (rank only; no threshold) ==========\n",
                metric))
  }
  
  # Read both windows.
  # Use all_data_allrep for shared Y limits and temperature allrep plots.
  # Use all_data_r1 for default R1 plots.
  all_data_allrep <- list()
  all_data_r1     <- list()
  
  for (window in c("full", "early", "last4")) {
    csv_path <- file.path(PROC_DIR, "fluctuations", window,
                          "community_level.csv")
    if (!file.exists(csv_path)) {
      cat(sprintf("  [skip] %s does not exist\n", csv_path))
      next
    }
    
    df_tmp_allrep <- read_csv(csv_path, show_col_types = FALSE)

    if (!metric %in% names(df_tmp_allrep)) {
      cat(sprintf("  [skip] %s does not contain column %s\n",
                  csv_path, metric))
      next
    }

    df_tmp_allrep <- df_tmp_allrep %>%
      filter(
        !is.na(.data[[metric]]),
        !(experiment == "temperature" & window == "early")
      )
    
    df_tmp_r1 <- df_tmp_allrep %>%
      filter(replica == 1)
    
    if (nrow(df_tmp_allrep) == 0) next
    
    all_data_allrep[[window]] <- df_tmp_allrep
    all_data_r1[[window]]     <- df_tmp_r1
    
    cat(sprintf("  %s: %d rows allrep, %d rows R1\n",
                window, nrow(df_tmp_allrep), nrow(df_tmp_r1)))
  }
  
  if (length(all_data_allrep) == 0) {
    cat(sprintf("  [skip] %s has no valid data\n", metric))
    next
  }
  
  # A separate shared Y upper limit for each metric
  # Use allrep data to avoid clipping those figures.
  y_max_global <- max(sapply(all_data_allrep, function(d)
    max(d[[metric]], na.rm = TRUE)))
  y_upper  <- ceiling(y_max_global * 10) / 10
  if (has_any_threshold) {
    y_upper  <- max(y_upper, max(metric_threshold_values, na.rm = TRUE) * 1.5)
  }
  y_breaks <- pretty(c(0, y_upper), n = 5)
  y_breaks <- y_breaks[y_breaks <= y_upper]
  
  cat(sprintf("  Shared Y-axis upper limit: %.2f\n", y_upper))
  
  for (window in c("full", "early", "last4")) {
    if (!window %in% names(all_data_r1)) next
    
    df_window_r1     <- all_data_r1[[window]]
    df_window_allrep <- all_data_allrep[[window]]
    
    fig_dir <- file.path(FIG_DIR, metric, window)
    dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
    
    cat(sprintf("\n  --- %s ---\n", window))
    
    for (experiment in c("mortality", "temperature")) {
      if (experiment == "temperature" && window == "early") next
      threshold <- threshold_for(metric, experiment)
      has_threshold <- is.finite(threshold)
      
      # ======================================================
      # Default version: R1 only.
      # ======================================================
      df_exp <- df_window_r1 %>%
        filter(experiment == !!experiment)
      
      if (nrow(df_exp) > 0) {
        cat(sprintf("    %s / R1 only\n", experiment))
        
        for (cond in paste0("W", 1:5)) {
          df_cond <- df_exp %>%
            filter(condition == cond)
          
          if (nrow(df_cond) == 0) next
          
          n_fluct  <- if (has_threshold) {
            sum(df_cond[[metric]] >= threshold, na.rm = TRUE)
          } else {
            NA_integer_
          }
          n_stable <- if (has_threshold) nrow(df_cond) - n_fluct else NA_integer_
          
          p <- plot_rank(df_cond,
                         metric      = metric,
                         threshold   = threshold,
                         y_upper     = y_upper,
                         y_breaks    = y_breaks,
                         show_legend = (cond == "W1"))
          
          ggsave(
            filename = file.path(fig_dir,
                                 sprintf("%s_%s.pdf", experiment, cond)),
            plot     = p,
            width    = W_PLOT,
            height   = H_PLOT,
            units    = "mm",
            device   = cairo_pdf
          )
          
          if (has_threshold) {
            cat(sprintf("      %s: n=%d  F=%d  S=%d\n",
                        cond, nrow(df_cond), n_fluct, n_stable))
          } else {
            cat(sprintf("      %s: n=%d\n", cond, nrow(df_cond)))
          }
        }
        
        p_combined <- plot_rank(df_exp,
                                metric      = metric,
                                threshold   = threshold,
                                y_upper     = y_upper,
                                y_breaks    = y_breaks,
                                show_legend = TRUE)
        
        ggsave(
          filename = file.path(fig_dir,
                               sprintf("%s_combined.pdf", experiment)),
          plot     = p_combined,
          width    = W_PLOT,
          height   = H_PLOT,
          units    = "mm",
          device   = cairo_pdf
        )
        
        cat(sprintf("      combined: n=%d\n", nrow(df_exp)))
      }
      
      # ======================================================
      # Additional temperature allrep version: R1/R2/R3 together.
      # Export this version for temperature only.
      # ======================================================
      if (experiment == "temperature") {
        
        df_exp_allrep <- df_window_allrep %>%
          filter(experiment == "temperature")
        
        if (nrow(df_exp_allrep) == 0) next
        
        cat(sprintf("    temperature / allrep\n"))
        
        for (cond in paste0("W", 1:5)) {
          df_cond_allrep <- df_exp_allrep %>%
            filter(condition == cond)
          
          if (nrow(df_cond_allrep) == 0) next
          
          n_fluct_allrep  <- if (has_threshold) {
            sum(df_cond_allrep[[metric]] >= threshold, na.rm = TRUE)
          } else {
            NA_integer_
          }
          n_stable_allrep <- if (has_threshold) {
            nrow(df_cond_allrep) - n_fluct_allrep
          } else {
            NA_integer_
          }
          
          p_allrep <- plot_rank(df_cond_allrep,
                                metric      = metric,
                                threshold   = threshold,
                                y_upper     = y_upper,
                                y_breaks    = y_breaks,
                                show_legend = (cond == "W1"))
          
          ggsave(
            filename = file.path(
              fig_dir,
              sprintf("temperature_%s_allrep.pdf", cond)
            ),
            plot     = p_allrep,
            width    = W_PLOT,
            height   = H_PLOT,
            units    = "mm",
            device   = cairo_pdf
          )
          
          if (has_threshold) {
            cat(sprintf("      %s_allrep: n=%d  F=%d  S=%d\n",
                        cond, nrow(df_cond_allrep),
                        n_fluct_allrep, n_stable_allrep))
          } else {
            cat(sprintf("      %s_allrep: n=%d\n",
                        cond, nrow(df_cond_allrep)))
          }
        }
        
        p_combined_allrep <- plot_rank(df_exp_allrep,
                                       metric      = metric,
                                       threshold   = threshold,
                                       y_upper     = y_upper,
                                       y_breaks    = y_breaks,
                                       show_legend = TRUE)
        
        ggsave(
          filename = file.path(fig_dir,
                               "temperature_combined_allrep.pdf"),
          plot     = p_combined_allrep,
          width    = W_PLOT,
          height   = H_PLOT,
          units    = "mm",
          device   = cairo_pdf
        )
        
        cat(sprintf("      combined_allrep: n=%d\n",
                    nrow(df_exp_allrep)))
      }
    }
  }
}

cat("\nDone!\n")
cat(sprintf("Output directory: figures/fluctuations/rank/\n"))
cat("  sum_abs_std/       (mortality threshold=0.20; temperature threshold=0.15)\n")
cat("  sum_rel_std/       (rank only; no threshold)\n")
cat("  temporal_bc_mean/  (rank only; no threshold)\n")
cat("  temporal_bc_max/   (threshold=0.20)\n")
cat("  weighted_log_sd/   (rank only; no threshold)\n")
cat("  community_cv/      (threshold=0.265)\n")
cat("  community_cv_rel/  (rank only; no threshold)\n")
cat("  total_biomass_std/ (rank only; no threshold)\n")
cat("  total_biomass_cv/  (rank only; no threshold)\n")
cat("\nAdditional temperature allrep figures:\n")
cat("  temperature_W1_allrep.pdf ~ temperature_W5_allrep.pdf\n")
cat("  temperature_combined_allrep.pdf\n")
