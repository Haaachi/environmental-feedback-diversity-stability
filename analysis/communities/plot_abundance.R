community_workspace <- function() {
  override <- Sys.getenv("COMMUNITY_WORKSPACE", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script)) stop("Run with Rscript or set COMMUNITY_WORKSPACE.")
  dirname(normalizePath(sub("^--file=", "", script[[1]]), mustWork = TRUE))
}

# ===========================================================
# plot_abundance.R
# 主目录: /home/hachi/Tem_mortality_workspace
#
# 统一输出所有丰度堆叠柱状图：
#   abs_full              abs_abund, D1-10,  R1
#   abs_early             abs_abund, D1-6,   R1
#   rel_full              rel_abund, D1-10,  R1
#   rel_early             rel_abund, D1-6,   R1
#   abs_last3_replicates  abs_abund, D8-10,  R1+R2+R3 横排
#   rel_last3_replicates  rel_abund, D8-10,  R1+R2+R3 横排
#
# Scope (按 experiment 限制任务，避免无谓出图):
#   mortality   : 仅 early 版本           (abs_early, rel_early)
#   temperature : 仅 full 相关             (abs_full, rel_full,
#                                          abs_last3_replicates,
#                                          rel_last3_replicates)
#
# 输出: figures/{SERIES_NAME}/{task_suffix}/{experiment}/{cond}/
#         {experiment}_{cond}_C{nn}.pdf          (单replicate任务)
#         {experiment}_{cond}_C{nn}_R1R2R3.pdf   (多replicate任务)
#
# SERIES_NAME 可通过环境变量 ABUNDANCE_SERIES 覆盖
# (默认 "series_01")，方便不同批次出图互不覆盖。
#
# ── 修复版变更说明 ────────────────────────────────────────────
# 1. COLORS 由全局命名向量改为按 experiment 拆分的列表
#    COLORS_BY_EXP，避免 mortality / temperature 同名 taxon
#    互相覆盖颜色。
# 2. 色板参考图改用复合 key (experiment::taxon) 染色，
#    facet 内的颜色与实际画图严格一致。
# 3. abs_abund 路径同样调用 complete()，保证不同 community 间
#    的堆叠顺序一致。
# ===========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

# ── 路径 ─────────────────────────────────────────────────────
BASE_DIR <- community_workspace()
PROC_DIR    <- file.path(BASE_DIR, "processed")
SERIES_NAME <- Sys.getenv("ABUNDANCE_SERIES", unset = "series_01")
FIG_DIR     <- file.path(BASE_DIR, "figures", SERIES_NAME)

# ── 全局图形参数 ─────────────────────────────────────────────
FONT_FAMILY <- "Arial"
FONT_AX     <- 7      # pt 刻度文字
FONT_LAB    <- 8      # pt 轴标签（保留备用）
BORDER_SIZE <- 0.25   # mm 边框线宽
TICK_SIZE   <- 0.20   # mm 刻度线宽
TICK_LEN    <- -1.0   # mm 刻度线长（负值=向内）

W_SINGLE <- 45   # mm 单子图宽度
H_SINGLE <- 42   # mm 单子图高度

# ── 任务定义 ─────────────────────────────────────────────────
# abund_col : 使用的丰度列名
# days      : 纳入的天数
# reps      : 纳入的replicate编号（长度>1时横向拼图）
# suffix    : 输出子目录名
TASKS <- list(
  abs_full             = list(abund_col = "abs_abund", days = 1:10,      reps = 1,   suffix = "abs_full"),
  abs_early            = list(abund_col = "abs_abund", days = 1:6,       reps = 1,   suffix = "abs_early"),
  rel_full             = list(abund_col = "rel_abund", days = 1:10,      reps = 1,   suffix = "rel_full"),
  rel_early            = list(abund_col = "rel_abund", days = 1:6,       reps = 1,   suffix = "rel_early"),
  abs_last3_replicates = list(abund_col = "abs_abund", days = c(8,9,10), reps = 1:3, suffix = "abs_last3_replicates"),
  rel_last3_replicates = list(abund_col = "rel_abund", days = c(8,9,10), reps = 1:3, suffix = "rel_last3_replicates")
)

# ── 每个 experiment 只跑哪些 task（关键提速开关）─────────────
# mortality   : 只有 D1-6 早期数据   → 仅 early
# temperature : 只有 full / last3   → 跳过 early
TASK_SCOPE <- list(
  mortality   = c("abs_early", "rel_early"),
  temperature = c("abs_full", "rel_full",
                  "abs_last3_replicates", "rel_last3_replicates")
)

# 可选：通过环境变量 ABUNDANCE_TASK_PATTERN 进一步筛选 task
TASK_PATTERN <- Sys.getenv("ABUNDANCE_TASK_PATTERN", unset = "")
if (nzchar(TASK_PATTERN)) {
  TASK_SCOPE <- map(TASK_SCOPE, ~ keep(.x, ~ str_detect(.x, TASK_PATTERN)))
}

# ── Y轴配置 ──────────────────────────────────────────────────
# abs：按实验分别设定（OD量级不同）
# rel：全局统一0-1
Y_CONFIG_ABS <- list(
  mortality   = list(limits = c(0, 2.2),
                     breaks = c(0, 1.1, 2.2),
                     labels = c("0", "1.1", "2.2")),
  temperature = list(limits = c(0, 1.2),
                     breaks = c(0, 0.6, 1.2),
                     labels = c("0", "0.6", "1.2"))
)
Y_CONFIG_REL <- list(
  limits = c(0, 1.001),
  breaks = c(0, 0.5, 1),
  labels = c("0", "0.5", "1")
)

# ── X轴刻度辅助函数 ──────────────────────────────────────────
x_breaks_labels <- function(days) {
  d <- sort(days)
  if (length(d) == 3) {
    # last3_replicates：显示全部3个天数
    list(breaks = as.character(d),
         labels = as.character(d))
  } else if (max(d) == 10) {
    list(breaks = c("1", "5", "10"),
         labels = c("1", "5", "10"))
  } else {
    list(breaks = c("1", "6"),
         labels = c("1", "6"))
  }
}

# ── 学术主题（单一定义）──────────────────────────────────────
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
      plot.title        = element_blank(),
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA),
      plot.margin       = margin(3, 3, 3, 3, "mm")
    )
}

# ── 核心画图函数（单个replicate的单张图）────────────────────
plot_one <- function(df, abund_col, experiment, days, rep_id) {
  
  if (abund_col == "abs_abund") {
    y_cfg <- Y_CONFIG_ABS[[experiment]]
  } else {
    y_cfg <- Y_CONFIG_REL
  }
  
  df_plot <- df %>%
    filter(replica == rep_id, day %in% days)
  
  if (nrow(df_plot) == 0) {
    return(
      ggplot() + theme_void() +
        theme(plot.background = element_rect(fill = "white", color = NA))
    )
  }
  
  # 统一从全局 taxon levels 取顺序（保证不同 community 堆叠一致）
  all_taxa <- COLOR_LEVELS[[experiment]]
  if (is.null(all_taxa)) all_taxa <- unique(df$taxon)
  
  # rel / abs 都补全缺失行 → 同样的因子层级、同样的堆叠顺序
  if (abund_col == "rel_abund") {
    df_plot <- df_plot %>%
      complete(taxon = all_taxa,
               day   = sort(days),
               fill  = list(rel_abund = 0)) %>%
      mutate(plot_val = replace_na(rel_abund, 0) / 100)
  } else {
    df_plot <- df_plot %>%
      complete(taxon = all_taxa,
               day   = sort(days),
               fill  = list(abs_abund = 0)) %>%
      mutate(plot_val = replace_na(abs_abund, 0))
  }
  
  if (nrow(df_plot) == 0 ||
      !any(is.finite(df_plot$plot_val)) ||
      sum(df_plot$plot_val, na.rm = TRUE) == 0) {
    return(
      ggplot() + theme_void() +
        theme(plot.background = element_rect(fill = "white", color = NA))
    )
  }
  
  df_plot <- df_plot %>% mutate(taxon = factor(taxon, levels = all_taxa))
  
  y_max_actual <- 0
  y_totals <- df_plot %>%
    group_by(day) %>%
    summarise(total = sum(plot_val, na.rm = TRUE), .groups = "drop") %>%
    pull(total)
  
  if (length(y_totals) > 0 && any(is.finite(y_totals))) {
    y_max_actual <- max(y_totals[is.finite(y_totals)], na.rm = TRUE)
  }
  
  if (y_max_actual > y_cfg$limits[2]) {
    warning(sprintf(
      "[%s R%d] Y轴超出上限: 实际最大值 %.4f > 设定上限 %.4f",
      experiment, rep_id, y_max_actual, y_cfg$limits[2]
    ))
  }
  
  x_cfg    <- x_breaks_labels(days)
  x_levels <- as.character(sort(days))
  
  # ★ 按 experiment 取颜色，避免 mortality/temperature 同名 taxon 互相覆盖
  exp_colors  <- COLORS_BY_EXP[[experiment]]
  used_colors <- exp_colors[names(exp_colors) %in% levels(droplevels(df_plot$taxon))]
  
  ggplot(df_plot,
         aes(x    = factor(day, levels = x_levels),
             y    = plot_val,
             fill = taxon)) +
    geom_col(position = "stack", width = 0.9, color = NA) +
    scale_fill_manual(values = used_colors, drop = FALSE) +
    scale_y_continuous(
      breaks   = y_cfg$breaks,
      labels   = y_cfg$labels,
      expand   = expansion(mult = c(0, 0)),
      sec.axis = dup_axis(labels = NULL, name = NULL)
    ) +
    scale_x_discrete(
      breaks = x_cfg$breaks,
      labels = x_cfg$labels,
      drop   = FALSE
    ) +
    coord_cartesian(
      ylim = y_cfg$limits,
      clip = "off"
    ) +
    theme_pub()
}

# ── 读取全局色板 ─────────────────────────────────────────────
color_df <- read_csv(file.path(PROC_DIR, "color_map.csv"),
                     show_col_types = FALSE)
if (!"display_label" %in% names(color_df)) {
  color_df <- color_df %>% mutate(display_label = taxon)
}

# ★ 按 experiment 拆分 taxon → color 映射
#   每个 experiment 独立的命名向量，键名是 taxon
COLORS_BY_EXP <- color_df %>%
  split(.$experiment) %>%
  map(~ setNames(.x$color, .x$taxon))

# 按 experiment 拆分 taxon 顺序（用于 factor levels 和堆叠顺序）
COLOR_LEVELS <- split(color_df$taxon, color_df$experiment)

# ── 主循环 ───────────────────────────────────────────────────
cat(sprintf("输出根目录: %s\n", FIG_DIR))
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

for (experiment in c("mortality", "temperature")) {
  cat(sprintf("\n=== %s ===\n", experiment))
  
  # 选取该 experiment 需要跑的 task
  scope_keys <- TASK_SCOPE[[experiment]]
  if (is.null(scope_keys) || length(scope_keys) == 0) {
    cat("  (跳过：scope 中无任务)\n")
    next
  }
  scoped_tasks <- TASKS[scope_keys]
  cat(sprintf("  待跑任务: %s\n", paste(scope_keys, collapse = ", ")))
  
  for (cond in paste0("W", 1:5)) {
    csv_path <- file.path(PROC_DIR, experiment, cond, "abs_abundance.csv")
    if (!file.exists(csv_path)) next
    
    df_cond <- read_csv(csv_path, show_col_types = FALSE)
    cat(sprintf("  %s\n", cond))
    
    for (task in scoped_tasks) {
      
      fig_dir <- file.path(FIG_DIR, task$suffix, experiment, cond)
      dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
      
      multi_rep <- length(task$reps) > 1
      
      for (comm_num in sort(unique(df_cond$community))) {
        df_comm <- df_cond %>% filter(community == comm_num)
        
        if (multi_rep) {
          # 多replicate：横向拼图
          plots <- map(task$reps, function(r) {
            plot_one(df_comm, task$abund_col, experiment,
                     task$days, rep_id = r)
          })
          combined <- wrap_plots(plots, nrow = 1) &
            theme(plot.margin = margin(3, 3, 3, 3, "mm"))
          
          fname <- sprintf("%s_%s_C%02d_R1R2R3.pdf",
                           experiment, cond, comm_num)
          ggsave(
            filename = file.path(fig_dir, fname),
            plot     = combined,
            width    = W_SINGLE * length(task$reps),
            height   = H_SINGLE,
            units    = "mm",
            device   = cairo_pdf
          )
        } else {
          # 单replicate：单张输出
          p <- plot_one(df_comm, task$abund_col, experiment,
                        task$days, rep_id = task$reps)
          
          fname <- sprintf("%s_%s_C%02d.pdf",
                           experiment, cond, comm_num)
          ggsave(
            filename = file.path(fig_dir, fname),
            plot     = p,
            width    = W_SINGLE,
            height   = H_SINGLE,
            units    = "mm",
            device   = cairo_pdf
          )
        }
      }
      cat(sprintf("    -> %s/%s/%s/ 完成\n",
                  task$suffix, experiment, cond))
    }
  }
}

# ── 色板参考图 ───────────────────────────────────────────────
cat("\n生成色板参考图...\n")

# ★ 用复合 key (experiment::taxon) 作为 fill aesthetic
#   这样两个 facet 即使 taxon 同名，也会分别拿到自己的颜色
color_df_for_palette <- color_df %>%
  mutate(
    fill_key = paste(experiment, taxon, sep = "::"),
    order    = row_number()
  )

palette_fill_lookup <- setNames(
  color_df_for_palette$color,
  color_df_for_palette$fill_key
)

p_pal <- color_df_for_palette %>%
  ggplot(aes(x = 1, y = reorder(taxon, -order), fill = fill_key)) +
  geom_tile(width = 0.55, height = 0.88,
            color = "white", linewidth = 0.15) +
  geom_text(aes(x = 1.38, label = display_label),
            hjust = 0, size = 5 / .pt,
            family = FONT_FAMILY) +
  facet_wrap(~experiment, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = palette_fill_lookup) +
  scale_x_continuous(limits = c(0.65, 5.9)) +
  theme_void(base_family = FONT_FAMILY) +
  theme(
    legend.position = "none",
    strip.text      = element_text(size = 7, face = "bold",
                                   margin = margin(b = 3)),
    panel.spacing   = unit(6, "mm"),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(5, 5, 5, 5, "mm")
  )

ggsave(
  filename = file.path(FIG_DIR, "color_palette_reference.pdf"),
  plot     = p_pal,
  width    = 120, height = 170,
  units    = "mm",
  device   = cairo_pdf
)

cat("完成！\n\n")
cat("输出结构:\n")
cat(sprintf("figures/%s/\n", SERIES_NAME))
cat("  color_palette_reference.pdf\n")
for (experiment in names(TASK_SCOPE)) {
  scope_keys <- TASK_SCOPE[[experiment]]
  for (key in scope_keys) {
    task    <- TASKS[[key]]
    rep_tag <- if (length(task$reps) > 1) "R1R2R3横排" else "R1"
    cat(sprintf("  %s/%s/  [%s, Day%s, %s]\n",
                task$suffix, experiment,
                task$abund_col,
                paste(range(task$days), collapse = "-"),
                rep_tag))
  }
}