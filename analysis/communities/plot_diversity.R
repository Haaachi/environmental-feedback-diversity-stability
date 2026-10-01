community_workspace <- function() {
  override <- Sys.getenv("COMMUNITY_WORKSPACE", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script)) stop("Run with Rscript or set COMMUNITY_WORKSPACE.")
  dirname(normalizePath(sub("^--file=", "", script[[1]]), mustWork = TRUE))
}

# ===========================================================
# plot_diversity.R
# 主目录: /home/hachi/Tem_mortality_workspace
#
# temperature：散点橙/紫区分震荡状态，Y轴固定范围
# mortality：  散点统一单色，Y轴上限贴近数据最大值
# 折线+errorbar：所有点的均值±SEM
# Y轴：0附近无padding
#
# 输出: figures/diversity/{experiment}/{window}/
# ===========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

BASE_DIR <- community_workspace()
PROC_DIR <- file.path(BASE_DIR, "processed")
FIG_DIR  <- file.path(BASE_DIR, "figures", "diversity_global_cv_rank")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 判据 ─────────────────────────────────────────────────────
COMMUNITY_CV_THRESHOLD <- 0.265

is_fluctuating_composite <- function(experiment, community_cv) {
  !is.na(community_cv) & community_cv >= COMMUNITY_CV_THRESHOLD
}

ensure_collapsed_col <- function(df) {
  if (!"collapsed" %in% names(df)) {
    df <- df %>% mutate(collapsed = FALSE)
  } else {
    df <- df %>% mutate(collapsed = if_else(is.na(collapsed), FALSE, collapsed))
  }
  df
}

# ── 配色 ─────────────────────────────────────────────────────
COLOR_FLUCT    <- "#F4A460"
COLOR_STABLE   <- "#9B8EC4"
COLOR_NODATA   <- "#BBBBBB"
CV_COLOR_LOW    <- "#0057FF"
CV_COLOR_HIGH   <- "#FF0000"
CV_COLOR_VALUES <- c(COLOR_STABLE, "#F2F2F2", COLOR_FLUCT)
CV_COLOR_LIMITS <- c(0, 1)
COLOR_MORTALITY <- "#7A9FC2"   # mortality单色（中性蓝灰）

# ── 图形参数 ─────────────────────────────────────────────────
FONT_FAMILY <- "Arial"
FONT_AX     <- 10
FONT_LAB    <- 12
BORDER_SIZE <- 0.25
TICK_SIZE   <- 0.20
TICK_LEN    <- -1.0

W_PLOT <- 80
H_PLOT <- 70

JITTER_WIDTH <- 0.15
DOT_SIZE     <- 1.0
DOT_ALPHA    <- 0.8
LINE_COLOR   <- "#2E6DA4"
LINE_WIDTH   <- 0.5
ERR_WIDTH    <- 0.12
ERR_LW       <- 0.35
MEAN_SIZE    <- 2.0

# ── X轴标签 ──────────────────────────────────────────────────
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


# ── 指标定义 ─────────────────────────────────────────────────
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

# ── 学术主题 ─────────────────────────────────────────────────
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
      axis.title         = element_blank(),
      legend.position    = "none",
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white", color = NA),
      plot.margin        = margin(3, 3, 3, 3, "mm")
    )
}

save_community_cv_legend <- function(out_dir) {
  legend_df <- tibble(
    x = 1,
    community_cv_rank = seq(CV_COLOR_LIMITS[1], CV_COLOR_LIMITS[2], length.out = 240)
  )

  p_legend <- ggplot(
    legend_df,
    aes(x = x, y = community_cv_rank, fill = community_cv_rank)
  ) +
    geom_tile(width = 0.28, height = diff(CV_COLOR_LIMITS) / 240) +
    scale_fill_gradientn(
      colors = CV_COLOR_VALUES,
      limits = CV_COLOR_LIMITS,
      oob = squish,
      guide = "none"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.35, 0.35))) +
    scale_y_continuous(
      name = "CV rank\nwithin experiment",
      limits = CV_COLOR_LIMITS,
      breaks = c(0, 0.5, 1),
      labels = c("Low", "Median", "High"),
      position = "right",
      expand = expansion(mult = c(0, 0))
    ) +
    theme_void(base_family = FONT_FAMILY) +
    theme(
      axis.title.y.right = element_text(
        size = 7,
        color = "black",
        angle = 90,
        margin = margin(l = 3)
      ),
      axis.text.y.right = element_text(size = 6, color = "black"),
      axis.ticks.y.right = element_line(linewidth = 0.2, color = "black"),
      axis.ticks.length.y.right = unit(1.0, "mm"),
      plot.margin = margin(2, 4, 2, 1, "mm")
    )

  ggsave(
    file.path(out_dir, "community_CV_rank_color_legend.pdf"),
    p_legend,
    width = 24,
    height = 42,
    units = "mm",
    device = cairo_pdf
  )
}

plot_metric <- function(df_metric, cv_df,
                        metric_col, experiment,
                        ylim_fixed = NULL) {
  
  x_labs <- X_LABELS[[experiment]]
  df_metric <- df_metric %>%
    ensure_collapsed_col() %>%
    filter(!(experiment == "temperature" & condition == "W5" & collapsed == TRUE))
  
  join_keys <- intersect(
    c("experiment", "condition", "community", "replica"),
    colnames(df_metric)
  )
  
  df_plot <- df_metric %>%
    left_join(
      cv_df %>%
        select(experiment, condition, community, replica,
               community_cv, community_cv_rank),
      by = join_keys
    ) %>%
    mutate(
      x_pos = as.integer(sub("W", "", condition)),
      state = case_when(
        is.na(community_cv)                           ~ "no_data",
        is_fluctuating_composite(experiment,
                                 community_cv)        ~ "fluctuating",
        TRUE                                          ~ "stable"
      )
    )
  
  if (nrow(df_plot) == 0) return(NULL)
  if (sum(!is.na(df_plot[[metric_col]])) == 0) return(NULL)
  
  # 均值和SEM
  summary_df <- df_plot %>%
    filter(!is.na(.data[[metric_col]])) %>%
    group_by(condition, x_pos) %>%
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
  
  # ── Y轴范围（统一逻辑）────────────────────────────────────
  data_max <- max(df_plot[[metric_col]], na.rm = TRUE)
  
  # ── Y轴范围（统一逻辑）────────────────────────────────────
  if (!is.null(ylim_fixed)) {
    y_lower <- ylim_fixed[1]
    y_upper <- ylim_fixed[2]
    # 生成刻度后严格截断，防止 pretty() 自动外扩改变 limits
    y_breaks <- pretty(c(y_lower, y_upper), n = 4)
    y_breaks <- y_breaks[y_breaks >= y_lower & y_breaks <= y_upper]
  } else {
    data_max <- max(df_plot[[metric_col]], na.rm = TRUE)
    y_lower <- 0
    y_upper <- data_max * 1.05
    y_breaks <- pretty(c(y_lower, y_upper), n = 4)
    y_upper <- max(y_breaks) # 动态模式下以最大刻度为上限
  }
  set.seed(42)
  
  # ── 散点层（统一颜色映射）────────────────────────────────
  scatter_layer <- geom_jitter(
    data   = df_plot %>%
      filter(!is.na(.data[[metric_col]])),
    aes(x     = x_pos,
        y     = .data[[metric_col]],
        color = community_cv_rank),
    width  = JITTER_WIDTH,
    height = 0,
    size   = DOT_SIZE,
    alpha  = DOT_ALPHA
  )
  
  p <- ggplot() +
    scatter_layer +
    geom_line(data      = summary_df,
              aes(x     = x_pos, y = mean_val),
              color     = LINE_COLOR,
              linewidth = LINE_WIDTH) +
    geom_errorbar(data  = summary_df,
                  aes(x    = x_pos,
                      ymin = ymin,
                      ymax = ymax),
                  width     = ERR_WIDTH,
                  linewidth = ERR_LW,
                  color     = LINE_COLOR) +
    geom_point(data  = summary_df,
               aes(x = x_pos, y = mean_val),
               shape = 16,
               size  = MEAN_SIZE,
               color = LINE_COLOR) +
    scale_x_continuous(
      breaks   = 1:5,
      labels   = x_labs,
      expand   = expansion(mult = c(0.10, 0.10)),
      sec.axis = dup_axis(labels = NULL, name = NULL)
    ) +
    scale_y_continuous(
      limits   = c(y_lower, y_upper),
      breaks   = y_breaks,
      expand   = expansion(mult = c(0, 0)),
      labels   = label_number(accuracy = 0.1),
      sec.axis = dup_axis(labels = NULL, name = NULL)
    ) +
    scale_color_gradientn(
      name     = "CV rank within experiment",
      colors   = CV_COLOR_VALUES,
      limits   = CV_COLOR_LIMITS,
      oob      = squish,
      na.value = COLOR_NODATA,
      guide = "none"
    ) +
    theme_pub()
  
  p
}

# ── 主循环 ───────────────────────────────────────────────────
save_community_cv_legend(FIG_DIR)

for (experiment in c("mortality", "temperature")) {
  for (window in c("full", "early")) {
    # 跳过 mortality / full
    if (experiment == "mortality" && window == "full") next
    if (experiment == "temperature" && window == "early") next
    
    cat(sprintf("\n=== %s / %s ===\n", experiment, window))
    
    fluc_path <- file.path(PROC_DIR, "fluctuations", window,
                           "community_level.csv")
    if (!file.exists(fluc_path)) {
      cat("  [跳过] fluctuations数据不存在\n")
      next
    }
    cv_df <- read_csv(fluc_path, show_col_types = FALSE) %>%
      filter(experiment == !!experiment) %>%
      ensure_collapsed_col() %>%
      mutate(community_cv = as.numeric(community_cv)) %>%
      filter(!(experiment == "temperature" & condition == "W5" & collapsed == TRUE)) %>%
      select(experiment, condition, community, replica,
             community_cv) %>%
      mutate(
        community_cv_rank = if (sum(!is.na(community_cv)) > 1) {
          percent_rank(community_cv)
        } else {
          if_else(is.na(community_cv), NA_real_, 0.5)
        }
      )
    
    fig_dir <- file.path(FIG_DIR, experiment, window)
    dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
    
    # ── Alpha diversity ───────────────────────────────────────
    alpha_path <- file.path(PROC_DIR, "diversity", window,
                            "alpha_diversity.csv")
    if (file.exists(alpha_path)) {
      alpha_df <- read_csv(alpha_path, show_col_types = FALSE) %>%
        filter(experiment == !!experiment)
      
      for (metric in names(ALPHA_METRICS)) {
        if (!metric %in% colnames(alpha_df)) next
        
        # --- 新增：动态判断并覆盖 ylim ---
        current_ylim <- ALPHA_METRICS[[metric]]$ylim
        if (experiment == "temperature" && metric == "survival_fraction") {
          current_ylim <- c(0, 0.8)
        }
        # ---------------------------------
        
        # 将传入的 ALPHA_METRICS[[metric]]$ylim 替换为 current_ylim
        p <- plot_metric(alpha_df, cv_df, metric, experiment,
                         current_ylim)
        
        if (is.null(p)) next
        ggsave(file.path(fig_dir, sprintf("%s.pdf", metric)),
               p, width = W_PLOT, height = H_PLOT,
               units = "mm", device = cairo_pdf)
        cat(sprintf("  -> %s.pdf\n", metric))
      }
    }
    # ── Gamma diversity ───────────────────────────────────────
    gamma_path <- file.path(PROC_DIR, "diversity", window,
                            "gamma_diversity.csv")
    if (file.exists(gamma_path)) {
      gamma_df <- read_csv(gamma_path, show_col_types = FALSE) %>%
        filter(experiment == !!experiment)
      
      for (metric in names(GAMMA_METRICS)) {
        if (!metric %in% colnames(gamma_df)) next
        p <- plot_metric(gamma_df, cv_df, metric, experiment,
                         GAMMA_METRICS[[metric]]$ylim)
        if (is.null(p)) next
        ggsave(file.path(fig_dir, sprintf("%s.pdf", metric)),
               p, width = W_PLOT, height = H_PLOT,
               units = "mm", device = cairo_pdf)
        cat(sprintf("  -> %s.pdf\n", metric))
      }
    }

    mean_daily_path <- file.path(PROC_DIR, "diversity", window,
                                 "mean_daily_diversity.csv")
    if (file.exists(mean_daily_path)) {
      mean_daily_df <- read_csv(mean_daily_path, show_col_types = FALSE) %>%
        filter(experiment == !!experiment)

      for (metric in names(MEAN_DAILY_METRICS)) {
        if (!metric %in% colnames(mean_daily_df)) next
        p <- plot_metric(mean_daily_df, cv_df, metric, experiment,
                         MEAN_DAILY_METRICS[[metric]]$ylim)
        if (is.null(p)) next
        ggsave(file.path(fig_dir, sprintf("%s.pdf", metric)),
               p, width = W_PLOT, height = H_PLOT,
               units = "mm", device = cairo_pdf)
        cat(sprintf("  -> %s.pdf\n", metric))
      }
    }
    
    # ── Beta within（仅full窗口）─────────────────────────────
    if (window == "full") {
      bw_path <- file.path(PROC_DIR, "diversity", window,
                           "beta_within.csv")
      if (file.exists(bw_path)) {
        bw_df <- read_csv(bw_path, show_col_types = FALSE) %>%
          filter(experiment == !!experiment)
        p <- plot_metric(bw_df, cv_df, "mean_bc", experiment,
                         c(0, 1))
        if (!is.null(p)) {
          ggsave(file.path(fig_dir, "beta_within.pdf"),
                 p, width = W_PLOT, height = H_PLOT,
                 units = "mm", device = cairo_pdf)
          cat("  -> beta_within.pdf\n")
        }
      }
    }
    
    # ── Beta across ───────────────────────────────────────────
    ba_path <- file.path(PROC_DIR, "diversity", window,
                         "beta_across.csv")
    if (file.exists(ba_path)) {
      ba_df <- read_csv(ba_path, show_col_types = FALSE) %>%
        filter(experiment == !!experiment)
      p <- plot_metric(ba_df, cv_df, "bc_vs_reference", experiment,
                       c(0, 1))
      if (!is.null(p)) {
        ggsave(file.path(fig_dir, "beta_across.pdf"),
               p, width = W_PLOT, height = H_PLOT,
               units = "mm", device = cairo_pdf)
        cat("  -> beta_across.pdf\n")
      }
    }
  }
}

cat("\n完成！\n")
cat(sprintf("  判据: community_cv >= %.2f\n", COMMUNITY_CV_THRESHOLD))
cat("  temperature: 橙=Fluctuation 紫=Stable 灰=无数据\n")
cat("  mortality:   散点统一蓝灰色，Y轴贴近数据最大值\n")
