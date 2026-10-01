community_workspace <- function() {
  override <- Sys.getenv("COMMUNITY_WORKSPACE", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script)) stop("Run with Rscript or set COMMUNITY_WORKSPACE.")
  dirname(normalizePath(sub("^--file=", "", script[[1]]), mustWork = TRUE))
}

# ===========================================================
# plot_mortality_diversity_cv_zones.R
# 主目录: /home/hachi/Tem_mortality_workspace
#
# mortality 工况 early + full 窗口
# shannon / gamma_shannon vs community_cv
# 背景按稳定性和多样性分区着色（3区域，无分隔线）
# 散点颜色: 震荡/稳定 | 形状: 稀释率
# 主图无图例，图例单独输出
# ===========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(cowplot)
})

# ============================================================
# ══ 绘图参数 ══════════════════════════════════════════════════
# ============================================================

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

N_TICKS <- 3

# ── 分区阈值 ─────────────────────────────────────────────────
CV_THRESHOLD        <- 0.25
DIVERSITY_THRESHOLD <- 0.8

# ── 背景分区颜色 ─────────────────────────────────────────────
ZONE_COLOR_STABLE        <- "#E8E4F0"   # 稳定（浅柔紫）
ZONE_COLOR_FLUC_LOW      <- "#FDE8D0"   # 震荡+低多样性（浅暖沙橙）
ZONE_COLOR_FLUC_HIGH     <- "#FAD4B5"   # 震荡+高多样性（中暖沙橙）

# ── 散点颜色（震荡 vs 稳定） ─────────────────────────────────
COL_STABLE      <- "#9B8EC4"   # 稳定柔紫
COL_FLUCTUATION <- "#F4A460"   # 震荡暖沙橙

# ── 稀释率形状 ───────────────────────────────────────────────
COND_SHAPES <- c(
  "W1" = 16,
  "W2" = 17,
  "W3" = 15,
  "W4" = 18,
  "W5" = 8
)

COND_LABELS <- c(
  "W1" = "10\u00B9",
  "W2" = "10\u00B2",
  "W3" = "10\u00B3",
  "W4" = "10\u2074",
  "W5" = "10\u2075"
)

# ── X/Y 轴范围 ───────────────────────────────────────────────
X_LO <- 0
X_HI <- 1.2
Y_LO <- 0
Y_HI <- 1.6

# ============================================================
# ── 路径配置 ──────────────────────────────────────────────────
# ============================================================
BASE_DIR <- community_workspace()
PROC_DIR <- file.path(BASE_DIR, "processed")
FIG_DIR  <- file.path(BASE_DIR, "figures", "mortality_diversity_cv_zones")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# ── 工具函数 ──────────────────────────────────────────────────
# ============================================================

make_ticks <- function(lo, hi, n = N_TICKS) {
  seq(lo, hi, length.out = n)
}

theme_pub <- function() {
  theme_classic(base_size = FONT_AX, base_family = FONT_FAMILY) +
    theme(
      axis.line          = element_blank(),
      panel.border       = element_rect(linewidth = BORDER_SIZE,
                                        color = "black", fill = NA),
      axis.ticks         = element_line(linewidth = TICK_SIZE,
                                        color = "black"),
      axis.ticks.length  = unit(TICK_LEN, "mm"),
      axis.text          = element_text(size = FONT_AX, color = "black"),
      axis.title         = element_blank(),
      legend.position    = "none",
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white", color = NA),
      plot.margin        = margin(3, 3, 3, 3, "mm")
    )
}

fmt_axis <- function(x) {
  formatC(x, format = "f", digits = 1)
}

axis_hi_1_decimal <- function(x, default_hi) {
  x_max <- suppressWarnings(max(x, na.rm = TRUE))
  if (!is.finite(x_max)) return(default_hi)
  max(default_hi, ceiling(x_max * 10) / 10)
}

parse_days <- function(days_string) {
  as.integer(stringr::str_extract_all(days_string[1], "\\d+")[[1]])
}

read_abs_abundance <- function(experiment_name) {
  purrr::map_dfr(paste0("W", 1:5), function(w) {
    path <- file.path(PROC_DIR, experiment_name, w, "abs_abundance.csv")
    if (!file.exists(path)) return(tibble())
    read_csv(path, show_col_types = FALSE)
  })
}

add_threshold_alpha_richness <- function(alpha_div, abs_df) {
  if (nrow(alpha_div) == 0 || nrow(abs_df) == 0) return(alpha_div)
  
  threshold_richness <- abs_df %>%
    group_by(experiment, condition, community, replica, day) %>%
    summarise(
      richness_1pct    = sum(rel_abund >= 0.01,  na.rm = TRUE),
      richness_0p1pct  = sum(rel_abund >= 0.001, na.rm = TRUE),
      .groups = "drop"
    )
  
  alpha_div %>%
    left_join(
      threshold_richness,
      by = c("experiment", "condition", "community", "replica", "day")
    )
}

add_threshold_gamma_richness <- function(gamma_div, abs_df) {
  if (nrow(gamma_div) == 0 || nrow(abs_df) == 0) return(gamma_div)
  
  window_days <- parse_days(gamma_div$days)
  
  threshold_gamma <- abs_df %>%
    filter(day %in% window_days) %>%
    group_by(experiment, condition, community, replica, taxon) %>%
    summarise(max_rel_abund = max(rel_abund, na.rm = TRUE), .groups = "drop") %>%
    group_by(experiment, condition, community, replica) %>%
    summarise(
      accumulative_richness_1pct   = sum(max_rel_abund >= 0.01,  na.rm = TRUE),
      accumulative_richness_0p1pct = sum(max_rel_abund >= 0.001, na.rm = TRUE),
      .groups = "drop"
    )
  
  gamma_div %>%
    left_join(
      threshold_gamma,
      by = c("experiment", "condition", "community", "replica")
    )
}

# ============================================================
# ── 核心画图函数：主散点图（无图例，无分隔线） ──────────────
# ============================================================
make_zone_scatter <- function(df_plot) {
  
  x_breaks <- make_ticks(X_LO, X_HI)
  y_hi <- axis_hi_1_decimal(df_plot$value, Y_HI)
  y_breaks <- make_ticks(Y_LO, y_hi)
  
  # 分区背景数据
  zone_rects <- data.frame(
    xmin = c(0,             CV_THRESHOLD,  CV_THRESHOLD),
    xmax = c(CV_THRESHOLD,  X_HI,          X_HI),
    ymin = c(0,             0,             DIVERSITY_THRESHOLD),
    ymax = c(y_hi,          DIVERSITY_THRESHOLD, y_hi),
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
    # 背景分区（无分隔线，纯颜色过渡）
    geom_rect(data = zone_rects,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = zone),
              alpha = 0.5, inherit.aes = FALSE) +
    scale_fill_manual(values = zone_colors) +
    
    # 散点：颜色 = 震荡/稳定，形状 = 稀释率
    geom_point(aes(color = fluctuating, shape = condition),
               size = DOT_SIZE, alpha = DOT_ALPHA) +
    scale_color_manual(values = c("TRUE" = COL_FLUCTUATION, 
                                  "FALSE" = COL_STABLE)) +
    scale_shape_manual(values = COND_SHAPES) +
    
    # 坐标轴
    scale_x_continuous(
      breaks   = x_breaks,
      labels   = fmt_axis(x_breaks),
      expand   = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      breaks   = y_breaks,
      labels   = fmt_axis(y_breaks),
      expand   = expansion(mult = c(0, 0))
    ) +
    
    coord_cartesian(xlim = c(X_LO, X_HI),
                    ylim = c(Y_LO, y_hi),
                    clip = "on") +
    
    labs(x = NULL, y = NULL) +
    theme_pub()
}

# ============================================================
# ── 图例构建函数 ──────────────────────────────────────────────
# ============================================================

# 图例1：震荡/稳定颜色
make_stability_legend_plot <- function() {
  legend_df <- data.frame(
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

# 图例2：稀释率形状（无文字标签）
# 图例2：稀释率形状（只显示标签，无标题）
make_dilution_shape_legend_plot <- function() {
  legend_df <- data.frame(
    condition = factor(paste0("W", 1:5), levels = paste0("W", 1:5)),
    x = 1:5,
    y = 1
  )
  
  ggplot(legend_df, aes(x = x, y = y, shape = condition)) +
    geom_point(size = 4, color = "black") +
    scale_shape_manual(
      values = COND_SHAPES,
      labels = COND_LABELS
    ) +
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

# 图例3：背景分区色块
make_zone_legend_plot <- function() {
  zone_df <- data.frame(
    zone = factor(
      c("Stable\n(CV < 0.25)",
        "Fluctuating + Low diversity\n(CV ≥ 0.25, D < 0.8)",
        "Fluctuating + High diversity\n(CV ≥ 0.25, D ≥ 0.8)"),
      levels = c("Stable\n(CV < 0.25)",
                 "Fluctuating + Low diversity\n(CV ≥ 0.25, D < 0.8)",
                 "Fluctuating + High diversity\n(CV ≥ 0.25, D ≥ 0.8)")
    ),
    x = 1:3,
    y = 1
  )
  
  zone_colors <- c(
    "Stable\n(CV < 0.25)"                                  = ZONE_COLOR_STABLE,
    "Fluctuating + Low diversity\n(CV ≥ 0.25, D < 0.8)"   = ZONE_COLOR_FLUC_LOW,
    "Fluctuating + High diversity\n(CV ≥ 0.25, D ≥ 0.8)"  = ZONE_COLOR_FLUC_HIGH
  )
  
  ggplot(zone_df, aes(x = x, y = y, fill = zone)) +
    geom_tile(width = 0.9, height = 0.9, color = "grey50", linewidth = 0.3) +
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
  
  legend_plot <- ggplot() +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0)) +
    annotation_custom(
      grob = legend_grob,
      xmin = -Inf, xmax = Inf,
      ymin = -Inf, ymax = Inf
    )
  
  ggsave(
    filename = filepath,
    plot = legend_plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    device = cairo_pdf
  )
}

# ============================================================
# ── 主流程 ────────────────────────────────────────────────────
# ============================================================

for (experiment in c("mortality", "temperature")) {
for (window in c("early", "full")) {
  if (experiment == "temperature" && window == "early") next
  
  cat(sprintf("\n========== %s - %s ==========\n", experiment, window))
  
  # ── 读取数据 ────────────────────────────────────────────────
  fluc_path <- file.path(PROC_DIR, "fluctuations", window, "community_level.csv")
  if (!file.exists(fluc_path)) {
    cat(sprintf("  [SKIP] Fluctuation file not found: %s\n", fluc_path))
    next
  }
  
  fluc_df <- read_csv(fluc_path, show_col_types = FALSE) %>%
    filter(experiment == !!experiment, collapsed == FALSE) %>%
    select(condition, community, replica, community_cv)
  
  alpha_path <- file.path(PROC_DIR, "diversity", window, "alpha_diversity.csv")
  gamma_path <- file.path(PROC_DIR, "diversity", window, "gamma_diversity.csv")
  
  if (!file.exists(alpha_path) || !file.exists(gamma_path)) {
    cat("  [SKIP] Diversity file not found\n")
    next
  }
  
  alpha_div <- read_csv(alpha_path, show_col_types = FALSE) %>%
    filter(experiment == !!experiment, collapsed == FALSE)
  
  gamma_div <- read_csv(gamma_path, show_col_types = FALSE) %>%
    filter(experiment == !!experiment, collapsed == FALSE)
  
  abs_df <- read_abs_abundance(experiment)
  alpha_div <- add_threshold_alpha_richness(alpha_div, abs_df)
  gamma_div <- add_threshold_gamma_richness(gamma_div, abs_df)
  
  cat(sprintf("  Fluctuation: %d rows\n", nrow(fluc_df)))
  cat(sprintf("  Alpha div: %d rows\n", nrow(alpha_div)))
  cat(sprintf("  Gamma div: %d rows\n", nrow(gamma_div)))
  
  build_plot_data <- function(div_df, div_metric) {
    if (!div_metric %in% names(div_df)) return(tibble())
    div_df %>%
      select(condition, community, replica, value = all_of(div_metric)) %>%
      inner_join(
        fluc_df %>% select(condition, community, replica, x_val = community_cv),
        by = c("condition", "community", "replica")
      ) %>%
      filter(!is.na(value), !is.na(x_val)) %>%
      mutate(
        condition = factor(condition, levels = paste0("W", 1:5)),
        fluctuating = x_val >= CV_THRESHOLD
      )
  }
  
  df_shannon       <- build_plot_data(alpha_div, "shannon")
  df_gamma_shannon <- build_plot_data(gamma_div, "gamma_shannon")
  
  # ── 绘制并保存主图 ──────────────────────────────────────────
  if (nrow(df_shannon) >= 1) {
    p1 <- make_zone_scatter(df_shannon)
    ggsave(
      filename = file.path(FIG_DIR, sprintf("shannon_vs_community_cv_%s_%s.pdf",
                                            experiment, window)),
      plot = p1, width = W_PLOT, height = H_PLOT, units = "mm", device = cairo_pdf
    )
    cat(sprintf("  -> shannon_vs_community_cv_%s_%s.pdf\n", experiment, window))
  }
  
  if (nrow(df_gamma_shannon) >= 1) {
    p2 <- make_zone_scatter(df_gamma_shannon)
    ggsave(
      filename = file.path(FIG_DIR, sprintf("gamma_shannon_vs_community_cv_%s_%s.pdf",
                                            experiment, window)),
      plot = p2, width = W_PLOT, height = H_PLOT, units = "mm", device = cairo_pdf
    )
    cat(sprintf("  -> gamma_shannon_vs_community_cv_%s_%s.pdf\n", experiment, window))
  }
  
  extra_alpha_metrics <- c(
    "richness",
    "richness_1pct",
    "richness_0p1pct"
  )
  
  extra_gamma_metrics <- c(
    "gamma_richness",
    "accumulative_richness_1pct",
    "accumulative_richness_0p1pct"
  )
  
  for (metric in extra_alpha_metrics) {
    df_metric <- build_plot_data(alpha_div, metric)
    if (nrow(df_metric) < 1) next
    p <- make_zone_scatter(df_metric)
    out_name <- sprintf("%s_vs_community_cv_%s_%s.pdf", metric, experiment, window)
    ggsave(
      filename = file.path(FIG_DIR, out_name),
      plot = p, width = W_PLOT, height = H_PLOT, units = "mm", device = cairo_pdf
    )
    cat(sprintf("  -> %s\n", out_name))
  }
  
  for (metric in extra_gamma_metrics) {
    df_metric <- build_plot_data(gamma_div, metric)
    if (nrow(df_metric) < 1) next
    p <- make_zone_scatter(df_metric)
    out_name <- sprintf("%s_vs_community_cv_%s_%s.pdf", metric, experiment, window)
    ggsave(
      filename = file.path(FIG_DIR, out_name),
      plot = p, width = W_PLOT, height = H_PLOT, units = "mm", device = cairo_pdf
    )
    cat(sprintf("  -> %s\n", out_name))
  }
}
}

# ============================================================
# ── 生成并保存图例 ────────────────────────────────────────────
# ============================================================

# 震荡/稳定颜色图例
p_legend_stability <- make_stability_legend_plot()
save_legend(
  p_legend_stability,
  file.path(FIG_DIR, "legend_stability.pdf"),
  width_mm = 60, height_mm = 15
)
cat("  -> legend_stability.pdf\n")

# 稀释率形状图例（两行排列）
p_legend_shape <- make_dilution_shape_legend_plot()
save_legend(
  p_legend_shape,
  file.path(FIG_DIR, "legend_dilution_shape.pdf"),
  width_mm = 80, height_mm = 20
)
cat("  -> legend_dilution_shape.pdf\n")

# 背景分区色块图例
p_legend_zone <- make_zone_legend_plot()
save_legend(
  p_legend_zone,
  file.path(FIG_DIR, "legend_zones.pdf"),
  width_mm = 100, height_mm = 25
)
cat("  -> legend_zones.pdf\n")

# ============================================================
# 完成提示
# ============================================================
cat("\n完成！\n\n")
cat("输出目录：figures/mortality_diversity_cv_zones/\n")
cat("\n生成文件：\n")
cat("  主图：\n")
cat("    shannon_vs_community_cv_mortality_early.pdf\n")
cat("    shannon_vs_community_cv_mortality_full.pdf\n")
cat("    gamma_shannon_vs_community_cv_mortality_early.pdf\n")
cat("    gamma_shannon_vs_community_cv_mortality_full.pdf\n")
cat("  图例：\n")
cat("    legend_stability.pdf        (震荡/稳定颜色)\n")
cat("    legend_dilution_shape.pdf   (稀释率形状)\n")
cat("    legend_zones.pdf            (背景分区色块)\n")
cat("\n设计说明：\n")
cat("  - 背景三色块：无分隔线，纯颜色区分\n")
cat("  - 散点颜色：震荡(暖沙橙) vs 稳定(柔紫)\n")
cat("  - 散点形状：W1-W5 区分稀释率\n")
cat(sprintf("  - 阈值: CV=%.2f, Diversity=%.1f\n", CV_THRESHOLD, DIVERSITY_THRESHOLD))
