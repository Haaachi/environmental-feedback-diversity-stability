community_workspace <- function() {
  override <- Sys.getenv("COMMUNITY_WORKSPACE", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script)) stop("Run with Rscript or set COMMUNITY_WORKSPACE.")
  dirname(normalizePath(sub("^--file=", "", script[[1]]), mustWork = TRUE))
}

# ===========================================================
# plot_fluctuation_proportion.R
# 主目录: /home/hachi/Tem_mortality_workspace
#
# 对震荡指标绘制两类图：
#   图A：震荡比例图（indicator >= threshold → 震荡）
#   图B：指标系综均值图（mean ± SEM，跨群落）
#
# 判据阈值在顶部 THRESHOLDS 统一调整
# 分母使用每个 experiment x window x condition 的实际观测数
# 误差棒：SEM（二项分布用于比例图，普通SEM用于均值图）
#
# 输出: figures/fluctuation_proportion/
#   <metric>/
#     {experiment}_{window}_proportion.pdf   ← 震荡比例
#     {experiment}_{window}_mean.pdf         ← 指标均值
#   oscillating_communities/
#     {metric}_{experiment}_{window}_oscillating.csv
# ===========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

BASE_DIR <- community_workspace()
PROC_DIR <- file.path(BASE_DIR, "processed")
FIG_DIR  <- file.path(BASE_DIR, "figures", "fluctuation_proportion")
OSC_DIR  <- file.path(FIG_DIR, "oscillating_communities")
STAT_DIR <- file.path(FIG_DIR, "proportion_tests")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OSC_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(STAT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# ── 阈值设置（在此处手动调整）───────────────────────────────
# ============================================================
THRESHOLDS <- list(
  community_cv = 0.265
)

BINARY_METRICS <- character(0)
MEAN_ONLY_METRICS <- c("sum_abs_std", "weighted_log_sd", "total_biomass_std")

KEY_CONTRASTS <- list(
  mortality = list(
    c("W1", "W2"),
    c("W2", "W3"),
    c("W3", "W4"),
    c("W4", "W5"),
    c("W1", "W3"),
    c("W3", "W5"),
    c("W1", "W5")
  ),
  temperature = list(
    c("W1", "W2"),
    c("W2", "W3"),
    c("W3", "W4"),
    c("W4", "W5"),
    c("W1", "W3"),
    c("W1", "W4")
  )
)

# ── 分母：每个 experiment × window 的总群落数 ────────────────
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

# ── 图形参数 ─────────────────────────────────────────────────
FONT_FAMILY <- "Arial"
FONT_AX     <- 10
BORDER_SIZE <- 0.25
TICK_SIZE   <- 0.20
TICK_LEN    <- -1.0

W_PLOT <- 80
H_PLOT <- 70

LINE_COLOR <- "#2E6DA4"
POINT_COLOR <- "#777777"
POINT_ALPHA <- 0.55
POINT_SIZE  <- 1.15
JITTER_WIDTH <- 0.055
SIG_LW <- 0.25
SIG_TEXT_SIZE <- 3.0

KEY_MEAN_CONTRASTS <- list(
  mortality = list(c("W1", "W5")),
  temperature = list(c("W3", "W2"), c("W3", "W1"))
)

# ── 学术主题（与原版完全一致）────────────────────────────────
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
      axis.title        = element_blank(),
      legend.position   = "none",
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA),
      plot.margin       = margin(3, 3, 3, 3, "mm")
    )
}

# ── 通用X轴 scale（两类图共用）───────────────────────────────
scale_x_condition <- function(experiment) {
  scale_x_continuous(
    breaks   = 1:5,
    labels   = X_LABELS[[experiment]],
    expand   = expansion(mult = c(0.08, 0.08)),
    sec.axis = dup_axis(labels = NULL, name = NULL)
  )
}

add_composite_metric <- function(df) {
  df %>%
    mutate(
      community_cv_threshold = THRESHOLDS$community_cv,
      community_cv_binary = as.numeric(
        !is.na(community_cv) & community_cv >= community_cv_threshold
      ),
      composite_fluctuation = community_cv_binary
    )
}

proportion_axis_upper <- function(summary_prop, experiment) {
  data_upper <- max(summary_prop$ymax, summary_prop$proportion, na.rm = TRUE)
  target <- if (experiment == "temperature") 0.60 else 0.82
  if (data_upper <= target) return(target)
  ceiling((data_upper + 0.03) * 10) / 10
}

proportion_axis_breaks <- function(y_upper) {
  if (abs(y_upper - 0.82) < 1e-9) return(seq(0, 0.8, by = 0.2))
  if (abs(y_upper - 0.60) < 1e-9) return(seq(0, 0.6, by = 0.2))

  step <- if (y_upper <= 0.6) 0.2 else 0.25
  br <- seq(0, y_upper, by = step)
  br[br <= y_upper]
}

run_pairwise_proportion_tests <- function(df, metric, threshold,
                                          experiment, window) {
  summary_counts <- df %>%
    group_by(condition) %>%
    summarise(
      n_osc = sum(.data[[metric]] >= threshold, na.rm = TRUE),
      n_total = n(),
      proportion = n_osc / n_total,
      .groups = "drop"
    )

  conds <- sort(unique(summary_counts$condition))
  if (length(conds) < 2) return(tibble())

  all_pairs <- combn(conds, 2, simplify = FALSE)
  key_pairs <- KEY_CONTRASTS[[experiment]]

  bind_rows(lapply(all_pairs, function(pair) {
    a <- summary_counts %>% filter(condition == pair[1])
    b <- summary_counts %>% filter(condition == pair[2])
    mat <- matrix(
      c(a$n_osc, a$n_total - a$n_osc,
        b$n_osc, b$n_total - b$n_osc),
      nrow = 2,
      byrow = TRUE
    )

    is_key <- any(vapply(
      key_pairs,
      function(k) identical(pair, k) || identical(rev(pair), k),
      logical(1)
    ))

    tibble(
      metric = metric,
      experiment = experiment,
      window = window,
      condition_a = pair[1],
      condition_b = pair[2],
      n_osc_a = a$n_osc,
      n_total_a = a$n_total,
      proportion_a = a$proportion,
      n_osc_b = b$n_osc,
      n_total_b = b$n_total,
      proportion_b = b$proportion,
      diff_b_minus_a = b$proportion - a$proportion,
      fisher_p = fisher.test(mat)$p.value,
      prop_test_p = suppressWarnings(
        prop.test(c(a$n_osc, b$n_osc), c(a$n_total, b$n_total),
                  correct = FALSE)$p.value
      ),
      key_contrast = is_key
    )
  })) %>%
    mutate(
      fisher_p_adj_BH = p.adjust(fisher_p, method = "BH"),
      prop_test_p_adj_BH = p.adjust(prop_test_p, method = "BH")
    )
}

p_to_sig_label <- function(p) {
  case_when(
    is.na(p)  ~ "NA",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

key_mean_tests <- function(df, metric, experiment, window) {
  key_pairs <- KEY_MEAN_CONTRASTS[[experiment]]
  if (is.null(key_pairs)) return(tibble())

  bind_rows(lapply(key_pairs, function(pair) {
    focal <- pair[1]
    ref <- pair[2]
    x_vals <- df %>%
      filter(condition == focal, !is.na(.data[[metric]])) %>%
      select(all_of(metric)) %>%
      pull(1)
    y_vals <- df %>%
      filter(condition == ref, !is.na(.data[[metric]])) %>%
      select(all_of(metric)) %>%
      pull(1)

    if (length(x_vals) < 2 || length(y_vals) < 2) return(tibble())

    test_res <- suppressWarnings(
      wilcox.test(x_vals, y_vals, alternative = "greater", exact = FALSE)
    )

    tibble(
      metric = metric,
      experiment = experiment,
      window = window,
      condition_focal = focal,
      condition_ref = ref,
      x1 = as.integer(sub("W", "", ref)),
      x2 = as.integer(sub("W", "", focal)),
      n_focal = length(x_vals),
      n_ref = length(y_vals),
      mean_focal = mean(x_vals, na.rm = TRUE),
      mean_ref = mean(y_vals, na.rm = TRUE),
      p_value = test_res$p.value,
      sig_label = p_to_sig_label(test_res$p.value)
    )
  }))
}

add_sig_brackets <- function(p, sig_df, y_lower, y_upper) {
  if (is.null(sig_df) || nrow(sig_df) == 0) return(p)

  sig_df <- sig_df %>%
    filter(!is.na(p_value), p_value < 0.05)

  if (nrow(sig_df) == 0) return(p)

  y_range <- y_upper - y_lower
  if (!is.finite(y_range) || y_range <= 0) y_range <- y_upper
  if (!is.finite(y_range) || y_range <= 0) y_range <- 1

  sig_df <- sig_df %>%
    mutate(
      idx = row_number(),
      y = y_upper - y_range * (0.08 + 0.075 * (idx - 1)),
      tick = y_range * 0.025
    )

  p +
    geom_segment(
      data = sig_df,
      aes(x = x1, xend = x2, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = SIG_LW,
      color = "black"
    ) +
    geom_segment(
      data = sig_df,
      aes(x = x1, xend = x1, y = y, yend = y - tick),
      inherit.aes = FALSE,
      linewidth = SIG_LW,
      color = "black"
    ) +
    geom_segment(
      data = sig_df,
      aes(x = x2, xend = x2, y = y, yend = y - tick),
      inherit.aes = FALSE,
      linewidth = SIG_LW,
      color = "black"
    ) +
    geom_text(
      data = sig_df,
      aes(x = (x1 + x2) / 2, y = y + tick * 0.3, label = sig_label),
      inherit.aes = FALSE,
      size = SIG_TEXT_SIZE,
      family = FONT_FAMILY,
      vjust = 0,
      color = "black"
    )
}

plot_metric_mean <- function(df, summary_mean, metric, experiment, window,
                             y_lower, y_upper, y_breaks) {
  sig_df <- key_mean_tests(df, metric, experiment, window)

  p <- ggplot(summary_mean, aes(x = x_pos, y = mean_v)) +
    geom_jitter(
      data = df,
      aes(x = x_pos, y = .data[[metric]]),
      inherit.aes = FALSE,
      width = JITTER_WIDTH,
      height = 0,
      size = POINT_SIZE,
      alpha = POINT_ALPHA,
      color = POINT_COLOR
    ) +
    geom_line(color = LINE_COLOR, linewidth = 0.6) +
    geom_errorbar(
      aes(ymin = ymin, ymax = ymax),
      width = 0.12,
      linewidth = 0.40,
      color = LINE_COLOR
    ) +
    geom_point(shape = 16, size = 1.8, color = LINE_COLOR) +
    scale_x_condition(experiment) +
    scale_y_continuous(
      limits = c(y_lower, y_upper),
      breaks = y_breaks,
      expand = expansion(mult = c(0, 0)),
      sec.axis = dup_axis(labels = NULL, name = NULL)
    ) +
    labs(x = NULL, y = NULL) +
    theme_pub()

  add_sig_brackets(p, sig_df, y_lower, y_upper)
}

# ============================================================
# 主循环：指标 × 实验 × 窗口
# ============================================================
all_prop_tests <- list()

for (metric in names(THRESHOLDS)) {
  threshold <- THRESHOLDS[[metric]]
  
  metric_dir <- file.path(FIG_DIR, metric)
  dir.create(metric_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat(sprintf("\n========== %s (threshold=%.2f) ==========\n",
              metric, threshold))
  
  for (experiment in c("mortality", "temperature")) {
    for (window in c("full", "early", "last4")) {
      if (experiment == "temperature" && window == "early") next
      
      csv_path <- file.path(PROC_DIR, "fluctuations", window,
                            "community_level.csv")
      if (!file.exists(csv_path)) next
      
      df <- read_csv(csv_path, show_col_types = FALSE) %>%
        filter(experiment == !!experiment)

      df <- add_composite_metric(df)

      df <- df %>%
        filter(!is.na(.data[[metric]]))
      
      if (nrow(df) == 0) next
      
      cat(sprintf("\n  --- %s / %s ---\n", experiment, window))
      
      # ── 震荡群落编号输出 ──────────────────────────────────
      osc_df <- df %>%
        filter(.data[[metric]] >= threshold) %>%
        arrange(condition, community, replica) %>%
        select(condition, community, replica,
               all_of(metric), collapsed) %>%
        mutate(across(all_of(metric), ~ round(.x, 4)))
      
      osc_path <- file.path(
        OSC_DIR,
        sprintf("%s_%s_%s_oscillating.csv", metric, experiment, window)
      )
      write_csv(osc_df, osc_path)
      cat(sprintf("  震荡群落: %d / %d\n", nrow(osc_df), nrow(df)))
      
      # x_pos：W1→1, W2→2, ...
      df <- df %>%
        mutate(x_pos = as.integer(sub("W", "", condition)))
      
      # ======================================================
      # 图A：震荡比例图
      # 分母 n_obs 含 collapsed（collapsed 在 df 中 cv=0，不超阈值）
      # ======================================================
      summary_prop <- df %>%
        group_by(condition, x_pos) %>%
        summarise(
          n_osc = sum(.data[[metric]] >= threshold, na.rm = TRUE),
          n_total = n(),
          .groups = "drop"
        ) %>%
        mutate(
          proportion = n_osc / n_total,
          sem        = sqrt(proportion * (1 - proportion) / n_total),
          ymin       = pmax(0, proportion - sem),
          ymax       = pmin(1, proportion + sem)
        ) %>%
        arrange(x_pos)

      all_prop_tests[[length(all_prop_tests) + 1]] <-
        run_pairwise_proportion_tests(
          df = df,
          metric = metric,
          threshold = threshold,
          experiment = experiment,
          window = window
        )
      
      y_upper_prop <- proportion_axis_upper(summary_prop, experiment)
      y_breaks_prop <- proportion_axis_breaks(y_upper_prop)
      
      p_prop <- ggplot(summary_prop,
                       aes(x = x_pos, y = proportion)) +
        geom_line(color     = LINE_COLOR,
                  linewidth = 0.6) +
        geom_errorbar(aes(ymin = ymin, ymax = ymax),
                      width     = 0.12,
                      linewidth = 0.40,
                      color     = LINE_COLOR) +
        geom_point(shape = 16,
                   size  = 1.8,
                   color = LINE_COLOR) +
        scale_x_condition(experiment) +
        scale_y_continuous(
          limits   = c(0, y_upper_prop),
          breaks   = y_breaks_prop,
          labels   = label_number(accuracy = 0.01),
          expand   = expansion(mult = c(0, 0)),
          sec.axis = dup_axis(labels = NULL, name = NULL)
        ) +
        labs(x = NULL, y = NULL) +
        theme_pub()
      
      ggsave(
        filename = file.path(
          metric_dir,
          sprintf("%s_%s_proportion.pdf", experiment, window)
        ),
        plot   = p_prop,
        width  = W_PLOT, height = H_PLOT, units = "mm",
        device = cairo_pdf
      )
      cat(sprintf("  -> %s_%s_proportion.pdf\n", experiment, window))

      if (metric %in% BINARY_METRICS) {
        cat(sprintf("  -> skip %s_%s_mean.pdf (binary criterion)\n",
                    experiment, window))
        next
      }
      
      # ======================================================
      # 图B：指标系综均值图（mean ± SEM）
      # collapsed 群落指标已在Python中置0，纳入平均
      # 避免50°C等高collapse条件下均值异常偏大
      # ======================================================
      summary_mean <- df %>%
        group_by(condition, x_pos) %>%
        summarise(
          n      = n(),
          mean_v = mean(.data[[metric]], na.rm = TRUE),
          sem_v  = sd(.data[[metric]],   na.rm = TRUE) / sqrt(n()),
          .groups = "drop"
        ) %>%
        mutate(
          ymin = mean_v - sem_v,
          ymax = mean_v + sem_v
        ) %>%
        arrange(x_pos)
      
      # Y轴：基于本实验×窗口数据自适应，上限留10%空白
      y_max_m  <- max(summary_mean$ymax, na.rm = TRUE)
      y_upper_m <- ceiling(y_max_m * 1.18 * 10) / 10
      y_upper_m <- max(y_upper_m, threshold * 1.5)
      y_lower_m <- 0
      if (experiment == "mortality" && metric == "community_cv") {
        y_lower_m <- 0.1
        y_upper_m <- 0.6
      }
      y_breaks_m <- pretty(c(y_lower_m, y_upper_m), n = 5)
      y_breaks_m <- y_breaks_m[
        y_breaks_m >= y_lower_m & y_breaks_m <= y_upper_m
      ]

      p_mean <- plot_metric_mean(
        df = df,
        summary_mean = summary_mean,
        metric = metric,
        experiment = experiment,
        window = window,
        y_lower = y_lower_m,
        y_upper = y_upper_m,
        y_breaks = y_breaks_m
      )
      
      ggsave(
        filename = file.path(
          metric_dir,
          sprintf("%s_%s_mean.pdf", experiment, window)
        ),
        plot   = p_mean,
        width  = W_PLOT, height = H_PLOT, units = "mm",
        device = cairo_pdf
      )
      cat(sprintf("  -> %s_%s_mean.pdf\n", experiment, window))
    }
  }
}

if (length(all_prop_tests) > 0) {
  prop_tests <- bind_rows(all_prop_tests) %>%
    arrange(metric, experiment, window, desc(key_contrast),
            condition_a, condition_b)

  write_csv(
    prop_tests,
    file.path(STAT_DIR, "pairwise_proportion_tests_all.csv")
  )
  write_csv(
    prop_tests %>% filter(key_contrast),
    file.path(STAT_DIR, "pairwise_proportion_tests_key_contrasts.csv")
  )

  cat(sprintf("\n  -> proportion tests: %s\n", STAT_DIR))
}

for (metric in MEAN_ONLY_METRICS) {
  metric_dir <- file.path(FIG_DIR, metric)
  dir.create(metric_dir, showWarnings = FALSE, recursive = TRUE)

  cat(sprintf("\n========== %s (mean only) ==========\n", metric))

  for (experiment in c("mortality", "temperature")) {
    for (window in c("full", "early", "last4")) {
      if (experiment == "temperature" && window == "early") next

      csv_path <- file.path(PROC_DIR, "fluctuations", window,
                            "community_level.csv")
      if (!file.exists(csv_path)) next

      df <- read_csv(csv_path, show_col_types = FALSE)
      if (!metric %in% names(df)) {
        cat(sprintf("  skip %s / %s: missing column %s\n",
                    experiment, window, metric))
        next
      }

      df <- df %>%
        filter(experiment == !!experiment,
               !is.na(.data[[metric]])) %>%
        mutate(x_pos = as.integer(sub("W", "", condition)))

      if (nrow(df) == 0) next

      summary_mean <- df %>%
        group_by(condition, x_pos) %>%
        summarise(
          n      = n(),
          mean_v = mean(.data[[metric]], na.rm = TRUE),
          sem_v  = sd(.data[[metric]],   na.rm = TRUE) / sqrt(n()),
          .groups = "drop"
        ) %>%
        mutate(
          ymin = pmax(0, mean_v - sem_v),
          ymax = mean_v + sem_v
        ) %>%
        arrange(x_pos)

      y_max_m  <- max(summary_mean$ymax, na.rm = TRUE)
      y_upper_m <- ceiling(y_max_m * 1.18 * 100) / 100
      y_upper_m <- max(y_upper_m, 0.01)
      y_lower_m <- 0
      if (experiment == "mortality" && metric == "community_cv") {
        y_lower_m <- 0.1
        y_upper_m <- 0.6
      }
      y_breaks_m <- pretty(c(y_lower_m, y_upper_m), n = 5)
      y_breaks_m <- y_breaks_m[
        y_breaks_m >= y_lower_m & y_breaks_m <= y_upper_m
      ]

      p_mean <- plot_metric_mean(
        df = df,
        summary_mean = summary_mean,
        metric = metric,
        experiment = experiment,
        window = window,
        y_lower = y_lower_m,
        y_upper = y_upper_m,
        y_breaks = y_breaks_m
      )

      ggsave(
        filename = file.path(
          metric_dir,
          sprintf("%s_%s_mean.pdf", experiment, window)
        ),
        plot   = p_mean,
        width  = W_PLOT, height = H_PLOT, units = "mm",
        device = cairo_pdf
      )
      cat(sprintf("  -> %s_%s_mean.pdf\n", experiment, window))
    }
  }
}

cat("\n完成！\n")
cat("figures/fluctuation_proportion/\n")
for (metric in names(THRESHOLDS)) {
  cat(sprintf("  %s/  (threshold=%.2f)\n", metric, THRESHOLDS[[metric]]))
  cat("    {experiment}_{window}_proportion.pdf\n")
  if (!(metric %in% BINARY_METRICS)) {
    cat("    {experiment}_{window}_mean.pdf\n")
  }
}
for (metric in MEAN_ONLY_METRICS) {
  cat(sprintf("  %s/  (mean only)\n", metric))
  cat("    {experiment}_{window}_mean.pdf\n")
}
cat("  oscillating_communities/\n")
cat("    {metric}_{experiment}_{window}_oscillating.csv\n")
cat("  proportion_tests/\n")
cat("    pairwise_proportion_tests_all.csv\n")
cat("    pairwise_proportion_tests_key_contrasts.csv\n")
