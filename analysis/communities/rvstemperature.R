community_workspace <- function() {
  override <- Sys.getenv("COMMUNITY_WORKSPACE", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script)) stop("Run with Rscript or set COMMUNITY_WORKSPACE.")
  dirname(normalizePath(sub("^--file=", "", script[[1]]), mustWork = TRUE))
}

# ===========================================================
# plot_growth_rate.R
# 主目录: /home/hachi/Tem_mortality_workspace
# ===========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
})

BASE_DIR <- community_workspace()
DATA_PATH <- file.path(BASE_DIR, "data", "rK_expfit_skip0.xlsx")
FIG_DIR   <- file.path(BASE_DIR, "figures", "growth_rate")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 参数 ─────────────────────────────────────────────────────
FONT_FAMILY <- "Arial"
FONT_AX     <- 8
BORDER_SIZE <- 0.25
TICK_SIZE   <- 0.20
TICK_LEN    <- -1.0

DOT_SIZE  <- 1.8
JIT_SIZE  <- 1.0
JIT_ALPHA <- 0.45
JIT_WIDTH <- 0.08
ERR_LW    <- 0.45
ERR_W     <- 0.10
W_3T <- 35; H_3T <- 45
W_2T <- 35;   H_2T <- 45

COL_TEMP <- c("20°C" = "#85C1E9",
              "30°C" = "#2E6DA4",
              "40°C" = "#C0392B")

# ── 主题（移除Y轴标题占位）───────────────────────────────────
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
      axis.text.x  = element_text(size = FONT_AX, color = "black",
                                  margin = margin(t = 3)),
      axis.text.y  = element_text(size = FONT_AX, color = "black",
                                  margin = margin(r = 3)),
      axis.title.x  = element_blank(),
      axis.title.y  = element_blank(),  # 👈 彻底清空 Y 轴标题区域
      legend.position  = "none",
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin      = margin(3, 3, 3, 3, "mm")
    )
}

# ── 显著性检验 ─────────────────────────────────────────────
auto_sig_label <- function(x1, x2) {
  x1 <- x1[!is.na(x1)]; x2 <- x2[!is.na(x2)]
  normal1 <- if (length(x1) >= 3) shapiro.test(x1)$p.value > 0.05 else FALSE
  normal2 <- if (length(x2) >= 3) shapiro.test(x2)$p.value > 0.05 else FALSE
  if (normal1 && normal2) {
    res <- t.test(x1, x2, var.equal = FALSE)
    cat("  检验方法：Welch t检验\n")
  } else {
    res <- wilcox.test(x1, x2, exact = FALSE)
    cat("  检验方法：Mann-Whitney U\n")
  }
  p <- res$p.value
  cat(sprintf("  p = %.4f\n", p))
  if      (p > 0.05)  "ns"
  else if (p > 0.01)  "*"
  else if (p > 0.001) "**"
  else                "***"
}

add_sig_annotation <- function(plt, label, y_bar, y_top,
                               x1 = 1, x2 = 2,
                               tick_h = NULL,
                               font_size = FONT_AX) {
  tick_h <- tick_h %||% (y_top * 0.02)
  seg_df <- data.frame(
    x    = c(x1, x1, x2, x2),
    xend = c(x1, x2, x2, x2),
    y    = c(y_bar - tick_h, y_bar, y_bar, y_bar - tick_h),
    yend = c(y_bar, y_bar, y_bar, y_bar - tick_h)
  )
  plt +
    geom_segment(data = seg_df,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE,
                 linewidth = 0.3, color = "black") +
    annotate("text", x = (x1 + x2) / 2, y = y_bar + tick_h * 0.5,
             label = label,
             size  = font_size / ggplot2::.pt,
             vjust = 0, family = FONT_FAMILY)
}

# ── 读取数据 ─────────────────────────────────────────────────
df_all <- read_excel(DATA_PATH) %>%
  mutate(Temperature = factor(Temperature,
                              levels = c(20, 30, 40),
                              labels = c("20°C", "30°C", "40°C")))

# ── 通用作图函数（固定Y轴上限 1.5，支持多组显著性）───────────
plot_r <- function(df, x_breaks, x_labels, x_limits, w, h, filename,
                   sig_list = NULL) {
  smry <- df %>%
    group_by(Temperature) %>%
    summarise(
      x   = cur_group_id(),
      mu  = mean(r, na.rm = TRUE),
      sem = sd(r, na.rm = TRUE) / sqrt(sum(!is.na(r))),
      .groups = "drop"
    )
  
  jit <- df %>% mutate(x = as.integer(Temperature))
  
  # ▼▼▼ 固定 Y 轴上限为 1.5 ▼▼▼
  y_max <- 1.5
  
  p <- ggplot() +
    geom_jitter(data = jit,
                aes(x = x, y = r, color = Temperature),
                width = JIT_WIDTH, size = JIT_SIZE,
                alpha = JIT_ALPHA, shape = 16) +
    geom_errorbar(data = smry,
                  aes(x = x, ymin = mu - sem, ymax = mu + sem,
                      color = Temperature),
                  width = ERR_W * 2, linewidth = ERR_LW) +
    geom_point(data = smry,
               aes(x = x, y = mu, color = Temperature),
               shape = 21, size = DOT_SIZE,
               fill = "white", stroke = 0.7) +
    scale_color_manual(values = COL_TEMP) +
    scale_x_continuous(breaks = x_breaks,
                       labels = x_labels,
                       limits = x_limits,
                       expand = expansion(0),
                       sec.axis = dup_axis(labels = NULL, name = NULL)) +
    scale_y_continuous(limits = c(0, y_max),
                       breaks = seq(0, y_max, by = 0.5),  # 显式设置刻度
                       expand = expansion(mult = c(0, 0)),
                       sec.axis = dup_axis(labels = NULL, name = NULL)) +
    labs(x = NULL, y = NULL) +
    theme_pub()
  
  # 循环添加显著性标注（自动垂直错开）
  if (!is.null(sig_list)) {
    for (i in seq_along(sig_list)) {
      s <- sig_list[[i]]
      # 高度基于固定 y_max=1.5 计算，避免重叠
      y_pos <- y_max * (0.82 + (i - 1) * 0.06)
      p <- add_sig_annotation(p,
                              label  = s$label,
                              y_bar  = y_pos,
                              y_top  = y_max,
                              x1     = s$x1,
                              x2     = s$x2,
                              tick_h = y_max * 0.025)
    }
  }
  
  ggsave(file.path(FIG_DIR, filename), p,
         width = w, height = h,
         units = "mm", device = cairo_pdf)
  cat(sprintf("-> %s\n", filename))
}

# ── 预计算所有组间显著性 ────────────────────────────────────
r20 <- df_all %>% filter(Temperature == "20°C") %>% pull(r)
r30 <- df_all %>% filter(Temperature == "30°C") %>% pull(r)
r40 <- df_all %>% filter(Temperature == "40°C") %>% pull(r)

cat("growth rate 20 vs 30 显著性检验：\n")
sig_20_30 <- auto_sig_label(r20, r30)

cat("growth rate 30 vs 40 显著性检验：\n")
sig_30_40 <- auto_sig_label(r30, r40)

# ── 图1：三温度（20 / 30 / 40）+ 双显著性标注 ───────────────
plot_r(df       = df_all,
       x_breaks = c(1, 2, 3),
       x_labels = c("20°C", "30°C", "40°C"),
       x_limits = c(0.6, 3.4),
       w        = W_3T, h = H_3T,
       filename = "growth_rate_r_3temp.pdf",
       sig_list = list(
         list(label = sig_20_30, x1 = 1, x2 = 2),
         list(label = sig_30_40, x1 = 2, x2 = 3)
       ))

# ── 图2：两温度（30 / 40）+ 单显著性标注 ─────────────────────
df_2t <- df_all %>%
  filter(Temperature %in% c("30°C", "40°C")) %>%
  mutate(Temperature = droplevels(Temperature))

plot_r(df       = df_2t,
       x_breaks = c(1, 2),
       x_labels = c("30°C", "40°C"),
       x_limits = c(0.6, 2.4),
       w        = W_2T, h = H_2T,
       filename = "growth_rate_r_30v40.pdf",
       sig_list = list(list(label = sig_30_40, x1 = 1, x2 = 2)))

cat("\n完成！figures/growth_rate/\n")
cat("  growth_rate_r_3temp.pdf   三温度 + 双显著性 (Y轴固定 0-1.5)\n")
cat("  growth_rate_r_30v40.pdf   30 vs 40 + 显著性 (Y轴固定 0-1.5)\n")