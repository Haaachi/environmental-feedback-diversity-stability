community_workspace <- function() {
  override <- Sys.getenv("COMMUNITY_WORKSPACE", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script)) stop("Run with Rscript or set COMMUNITY_WORKSPACE.")
  dirname(normalizePath(sub("^--file=", "", script[[1]]), mustWork = TRUE))
}

# ===========================================================
# plot_taxon_heatmap.R
# 主目录: /home/hachi/Tem_mortality_workspace
#
# 每个taxon在各工况下的平均abs_abund热图
# X轴：工况（W1-W5）
# Y轴：taxon（按总biomass排序）
# 颜色：该工况下所有群落×replicate的mean abs_abund
#       只统计rel_abund > 1%的存活记录（低于检测限视为0）
#
# 时间点：full=D10，early=D6
#
# 输出: figures/taxon_heatmap/
#   {experiment}_{window}_taxon_heatmap.pdf
# ===========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

BASE_DIR <- community_workspace()
PROC_DIR <- file.path(BASE_DIR, "processed")
FIG_DIR  <- file.path(BASE_DIR, "figures", "taxon_heatmap")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ── 参数 ─────────────────────────────────────────────────────
REL_THRESHOLD <- 1.0   # rel_abund > 1% 才视为存活

WINDOWS <- list(
  full  = list(day = 10, reps = c(1, 2, 3)),
  early = list(day = 6,  reps = c(1))
)

# ── X轴标签 ──────────────────────────────────────────────────
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

# ── 图形参数 ─────────────────────────────────────────────────
FONT_FAMILY <- "Arial"
FONT_AX     <- 7
FONT_TAXON  <- 6      # taxon名字号
BORDER_SIZE <- 0.25
TICK_SIZE   <- 0.20
TICK_LEN    <- -1.0

W_PLOT <- 60    # mm，X轴5个工况
# H_PLOT按taxon数量动态计算

# 配色：白→深蓝，零值为白色
HEATMAP_LOW  <- "white"
HEATMAP_HIGH <- "#1A4E8A"
HEATMAP_NA   <- "grey95"   # 完全缺失的格子

# ── 学术主题 ─────────────────────────────────────────────────
theme_pub <- function() {
  theme_classic(base_size = FONT_AX, base_family = FONT_FAMILY) +
    theme(
      axis.line              = element_blank(),
      panel.border           = element_rect(linewidth = BORDER_SIZE,
                                            color = "black", fill = NA),
      axis.ticks.x           = element_line(linewidth = TICK_SIZE,
                                            color = "black"),
      axis.ticks.y           = element_blank(),
      axis.ticks.x.top       = element_line(linewidth = TICK_SIZE,
                                            color = "black"),
      axis.ticks.y.right     = element_blank(),
      axis.ticks.length               = unit(TICK_LEN, "mm"),
      axis.ticks.length.x.top         = unit(TICK_LEN, "mm"),
      axis.text.x.top        = element_blank(),
      axis.text.y.right      = element_blank(),
      axis.text.x  = element_text(size   = FONT_AX,
                                  color  = "black",
                                  margin = margin(t = 3)),
      axis.text.y  = element_text(size   = FONT_TAXON,
                                  color  = "black",
                                  hjust  = 1,
                                  margin = margin(r = 2)),
      axis.title        = element_blank(),
      legend.position   = "right",
      legend.key.width  = unit(2, "mm"),
      legend.key.height = unit(8, "mm"),
      legend.text       = element_text(size = FONT_AX,
                                       color = "black"),
      legend.title      = element_text(size = FONT_AX,
                                       color = "black"),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA),
      plot.margin       = margin(3, 3, 3, 3, "mm")
    )
}

# ── 主循环 ───────────────────────────────────────────────────
for (experiment in c("mortality", "temperature")) {
  for (window_name in c("full", "early")) {
    if (experiment == "temperature" && window_name == "early") next
    
    cfg     <- WINDOWS[[window_name]]
    last_day <- cfg$day
    use_reps <- cfg$reps
    
    cat(sprintf("\n=== %s / %s (D%d, R%s) ===\n",
                experiment, window_name, last_day,
                paste(use_reps, collapse = "+")))
    
    # 读取所有工况数据
    records <- list()
    
    for (cond in paste0("W", 1:5)) {
      csv_path <- file.path(PROC_DIR, experiment, cond,
                            "abs_abundance.csv")
      if (!file.exists(csv_path)) next
      
      df <- read_csv(csv_path, show_col_types = FALSE) %>%
        filter(day     == last_day,
               replica %in% use_reps)
      
      if (nrow(df) == 0) next
      
      # 只保留存活记录（rel_abund > 1%）
      # 不存活的abs_abund视为0，不纳入均值计算
      df_alive <- df %>%
        filter(rel_abund > REL_THRESHOLD) %>%
        mutate(condition = cond)
      
      records[[cond]] <- df_alive
    }
    
    if (length(records) == 0) next
    
    df_all <- bind_rows(records)
    
    # ── 计算每个taxon × 工况的mean abs_abund ────────────────
    # 分母：该工况下所有群落×replicate的数量（存活与否都计入分母）
    # 即：mean = sum(abs_abund_alive) / n_total_samples
    # 这样不存活的样本贡献0，反映真实平均
    
    # 每个工况的总样本数（群落×replicate）
    n_samples <- map_int(paste0("W", 1:5), function(cond) {
      csv_path <- file.path(PROC_DIR, experiment, cond,
                            "abs_abundance.csv")
      if (!file.exists(csv_path)) return(0L)
      read_csv(csv_path, show_col_types = FALSE) %>%
        filter(day == last_day, replica %in% use_reps) %>%
        distinct(community, replica) %>%
        nrow()
    }) %>% setNames(paste0("W", 1:5))
    
    heatmap_df <- df_all %>%
      group_by(condition, taxon) %>%
      summarise(sum_abs = sum(abs_abund, na.rm = TRUE),
                .groups = "drop") %>%
      mutate(
        n_total  = n_samples[condition],
        mean_abs = sum_abs / n_total
      )
    
    # ── Taxon排序：按总biomass降序 ───────────────────────────
    taxon_order <- heatmap_df %>%
      group_by(taxon) %>%
      summarise(total = sum(mean_abs), .groups = "drop") %>%
      arrange(desc(total)) %>%
      pull(taxon)
    
    # ── 补全缺失格子（某工况下某taxon完全未存活）────────────
    all_conds <- paste0("W", 1:5)
    heatmap_full <- expand_grid(
      condition = all_conds,
      taxon     = taxon_order
    ) %>%
      left_join(heatmap_df %>% select(condition, taxon, mean_abs),
                by = c("condition", "taxon")) %>%
      replace_na(list(mean_abs = 0)) %>%
      mutate(
        condition = factor(condition, levels = all_conds),
        taxon     = factor(taxon, levels = rev(taxon_order))
      )
    
    n_taxa <- length(taxon_order)
    cat(sprintf("  存活taxon数: %d\n", n_taxa))
    
    # H动态计算：每个taxon 3mm + 上下边距
    H_PLOT <- n_taxa * 3 + 10
    
    # ── 绘图 ─────────────────────────────────────────────────
    p <- ggplot(heatmap_full,
                aes(x    = condition,
                    y    = taxon,
                    fill = mean_abs)) +
      geom_tile(color = "white", linewidth = 0.3) +
      scale_fill_gradient(
        low      = HEATMAP_LOW,
        high     = HEATMAP_HIGH,
        na.value = HEATMAP_NA,
        name     = "Mean\nabs abund",
        labels   = label_number(accuracy = 0.001)
      ) +
      scale_x_discrete(
        labels   = X_LABELS[[experiment]],
        expand   = expansion(mult = c(0, 0))
      ) +
      scale_y_discrete(
        expand = expansion(mult = c(0, 0))
      ) +
      theme_pub()
    
    # 移除上方重复的sec.axis（热图用不到）
    p <- p + guides(fill = guide_colorbar(
      barwidth  = unit(2,  "mm"),
      barheight = unit(15, "mm"),
      ticks     = TRUE
    ))
    
    out_path <- file.path(FIG_DIR,
                          sprintf("%s_%s_taxon_heatmap.pdf",
                                  experiment, window_name))
    ggsave(
      filename = out_path,
      plot     = p,
      width    = W_PLOT,
      height   = H_PLOT,
      units    = "mm",
      device   = cairo_pdf
    )
    cat(sprintf("  -> %s_%s_taxon_heatmap.pdf\n",
                experiment, window_name))
  }
}

cat("\n完成！\n")
cat("figures/taxon_heatmap/\n")
cat("  mortality_full_taxon_heatmap.pdf\n")
cat("  mortality_early_taxon_heatmap.pdf\n")
cat("  temperature_full_taxon_heatmap.pdf\n")
cat(sprintf("\n  存活判定: rel_abund > %.1f%%\n", REL_THRESHOLD))
cat("  颜色: mean abs_abund（不存活样本贡献0）\n")
cat("  Y轴: taxon按总biomass降序排列\n")
