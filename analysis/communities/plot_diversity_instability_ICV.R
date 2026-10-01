community_workspace <- function() {
  override <- Sys.getenv("COMMUNITY_WORKSPACE", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script)) stop("Run with Rscript or set COMMUNITY_WORKSPACE.")
  dirname(normalizePath(sub("^--file=", "", script[[1]]), mustWork = TRUE))
}

# ===========================================================
# plot_diversity_instability_ICV.R
# Robust diversity-instability plots.
#
# X metrics:
#   log_invariability_comm  = -log10(community_cv)
#   log_invariability_total = -log10(total_biomass_cv)
#
# Y metrics:
#   effective_species       = exp(shannon)
#   gamma_effective_species = exp(gamma_shannon)
#
# Points are plotted in a single color; no threshold-based classes are used.
#
# Main analysis windows:
#   mortality   -> early  (last 3 days, R1 only)
#   temperature -> full   (last 3 days, all replicates)
#
# Output:
#   figures/diversity_instability_robust/
#     {pearson|spearman}/{experiment}/{window}/{fluc_metric}_vs_{div_metric}/
# ===========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

# -----------------------------------------------------------
# Paths
# -----------------------------------------------------------
BASE_DIR <- community_workspace()
PROC_DIR <- file.path(BASE_DIR, "processed")
FIG_DIR  <- file.path(BASE_DIR, "figures", "diversity_instability_robust")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------
# Plot parameters
# -----------------------------------------------------------
W_PLOT <- 35
H_PLOT <- 30

FONT_FAMILY <- "Arial"
FONT_AX     <- 7
FONT_ANNOT  <- 6.5

BORDER_SIZE <- 0.25
TICK_SIZE   <- 0.20
TICK_LEN    <- -0.4

DOT_SIZE  <- 1.2
DOT_ALPHA <- 0.75
REG_LW    <- 0.35

P_THRESHOLD <- 0.05
COMMUNITY_CV_THRESHOLD <- 0.265
LOG_INV_COMM_THRESHOLD <- -log10(COMMUNITY_CV_THRESHOLD)
CV_FLOOR <- 1e-3
N_TICKS <- 3

COL_SINGLE   <- "#777777"
COLOR_STABLE <- "#9B8EC4"
COLOR_FLUCT  <- "#F4A460"
COL_POS_SIG  <- "#C75B4E"
COL_NEG_SIG  <- "#4A7FB5"
COL_NS       <- "#B8B8B8"
ANNOT_COL_SIG <- "#5A5A5A"
ANNOT_COL_NS  <- "#AAAAAA"

# -----------------------------------------------------------
# Analysis configuration
# -----------------------------------------------------------
EXP_WINDOW_MAP <- list(
  mortality   = c("early"),
  temperature = c("full")
)

FLUC_XLIM <- list(
  log_invariability_comm  = NULL,
  log_invariability_total = NULL
)

FLUC_METRICS <- list(
  log_invariability_comm = list(
    source_col = "community_cv",
    transform = function(x) -log10(pmax(x, CV_FLOOR)),
    use_color = FALSE
  ),
  log_invariability_total = list(
    source_col = "total_biomass_cv",
    transform = function(x) -log10(pmax(x, CV_FLOOR)),
    use_color = FALSE
  )
)

DIVERSITY_METRICS <- list(
  effective_species = list(
    file = "alpha_diversity.csv",
    source_col = "shannon",
    transform = function(x) exp(x),
    ylim = NULL
  ),
  gamma_effective_species = list(
    file = "gamma_diversity.csv",
    source_col = "gamma_shannon",
    transform = function(x) exp(x),
    ylim = NULL
  )
)

# -----------------------------------------------------------
# Helper functions
# -----------------------------------------------------------
is_valid_combo <- function(experiment, window) {
  allowed <- EXP_WINDOW_MAP[[experiment]]
  !is.null(allowed) && window %in% allowed
}

make_ticks <- function(lo, hi, n = N_TICKS) {
  seq(lo, hi, length.out = n)
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
      axis.text.x = element_text(
        size = FONT_AX, color = "black", margin = margin(t = 3)
      ),
      axis.text.y = element_text(
        size = FONT_AX, color = "black", margin = margin(r = 3)
      ),
      axis.title = element_blank(),
      legend.position = "none",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(3, 3, 3, 3, "mm")
    )
}

fmt_axis <- function(x) {
  ifelse(
    x == floor(x),
    as.character(as.integer(x)),
    gsub("0+$", "", formatC(x, format = "f", digits = 3))
  )
}

nice_range_signed <- function(vals) {
  vals <- vals[!is.na(vals) & is.finite(vals)]
  if (length(vals) == 0) return(list(lo = 0, hi = 1))

  lo <- min(vals)
  hi <- max(vals)
  if (lo == hi) {
    pad <- max(abs(lo) * 0.1, 0.1)
    return(list(lo = lo - pad, hi = hi + pad))
  }

  pad <- (hi - lo) * 0.05
  br <- pretty(c(lo - pad, hi + pad), n = 4)
  list(lo = min(br), hi = max(br))
}

nice_range_positive <- function(vals) {
  vals <- vals[!is.na(vals) & is.finite(vals)]
  if (length(vals) == 0) return(list(lo = 0, hi = 1))

  hi <- max(vals)
  if (!is.finite(hi) || hi <= 0) hi <- 1
  br <- pretty(c(0, hi * 1.05), n = 4)
  list(lo = 0, hi = max(br))
}

get_xlim <- function(fluc_metric, experiment) {
  xlim_cfg <- FLUC_XLIM[[fluc_metric]]
  if (is.null(xlim_cfg)) return(NULL)

  if (is.list(xlim_cfg) && !is.null(names(xlim_cfg))) {
    if (experiment %in% names(xlim_cfg)) return(xlim_cfg[[experiment]])
    warning(sprintf(
      "No xlim for experiment '%s' in fluc_metric '%s'; using auto range",
      experiment, fluc_metric
    ))
    return(NULL)
  }

  if (is.numeric(xlim_cfg) && length(xlim_cfg) == 2) return(xlim_cfg)
  NULL
}

calc_corr_single <- function(x, y, method) {
  df <- data.frame(x = x, y = y) %>%
    filter(!is.na(x), !is.na(y), is.finite(x), is.finite(y))
  n <- nrow(df)
  if (n < 3) return(list(r = NA_real_, p = NA_real_, n = n))

  res <- tryCatch(
    cor.test(
      df$x, df$y,
      method = method,
      exact = if (method == "spearman") FALSE else NULL
    ),
    error = function(e) NULL
  )

  list(
    r = if (!is.null(res)) round(as.numeric(res$estimate), 2) else NA_real_,
    p = if (!is.null(res)) res$p.value else NA_real_,
    n = n
  )
}

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("p<0.001")
  if (p < 0.01) return(sprintf("p=%.3f", p))
  sprintf("p=%.2f", p)
}

fmt_annot_single <- function(corr, method) {
  prefix <- if (method == "pearson") "r" else "rho"
  if (is.na(corr$r)) return(sprintf("%s=NA", prefix))
  sprintf("%s=%.2f, %s", prefix, corr$r, fmt_p(corr$p))
}

is_sig <- function(corr) {
  !is.na(corr$p) && corr$p < P_THRESHOLD
}

reg_style <- function(corr) {
  if (!is_sig(corr)) return(list(color = COL_NS, linetype = "dashed"))
  list(
    color = if (!is.na(corr$r) && corr$r > 0) COL_POS_SIG else COL_NEG_SIG,
    linetype = "solid"
  )
}

classify_by_community_cv <- function(community_cv) {
  if_else(
    !is.na(community_cv) & community_cv >= COMMUNITY_CV_THRESHOLD,
    "Fluctuation",
    "Stable"
  )
}

make_scatter <- function(df_plot, fluc_cfg, corr_method, x_lo, x_hi, y_lo, y_hi) {
  x_breaks <- make_ticks(x_lo, x_hi)
  y_breaks <- make_ticks(y_lo, y_hi)

  corr <- calc_corr_single(df_plot$x_val, df_plot$value, corr_method)
  annot <- fmt_annot_single(corr, corr_method)
  style <- reg_style(corr)
  annot_col <- if (is_sig(corr)) ANNOT_COL_SIG else ANNOT_COL_NS

  point_layer <- if (isTRUE(fluc_cfg$use_color)) {
    geom_point(aes(color = fluctuation_class), size = DOT_SIZE, alpha = DOT_ALPHA)
  } else {
    geom_point(color = COL_SINGLE, size = DOT_SIZE, alpha = DOT_ALPHA)
  }

  color_scale <- if (isTRUE(fluc_cfg$use_color)) {
    scale_color_manual(
      values = c("Stable" = COLOR_STABLE, "Fluctuation" = COLOR_FLUCT),
      drop = FALSE
    )
  } else {
    NULL
  }

  ggplot(df_plot, aes(x = x_val, y = value)) +
    point_layer +
    geom_smooth(
      method = "lm",
      se = FALSE,
      linewidth = REG_LW,
      linetype = style$linetype,
      color = style$color,
      fullrange = FALSE
    ) +
    annotate(
      "text",
      x = x_hi,
      y = y_hi * 0.98,
      label = annot,
      hjust = 1,
      vjust = 1,
      size = FONT_ANNOT / .pt,
      family = FONT_FAMILY,
      lineheight = 0.9,
      color = annot_col
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = fmt_axis(x_breaks),
      expand = expansion(mult = c(0, 0)),
      sec.axis = dup_axis(labels = NULL, name = NULL)
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = fmt_axis(y_breaks),
      expand = expansion(mult = c(0, 0)),
      sec.axis = dup_axis(labels = NULL, name = NULL)
    ) +
    coord_cartesian(xlim = c(x_lo, x_hi), ylim = c(y_lo, y_hi), clip = "on") +
    labs(x = NULL, y = NULL) +
    color_scale +
    theme_pub()
}

# -----------------------------------------------------------
# Main loop
# -----------------------------------------------------------
for (corr_method in c("pearson", "spearman")) {
  cat(sprintf("\n========== %s ==========\n", toupper(corr_method)))

  for (experiment in c("mortality", "temperature")) {
    for (window in c("full", "early")) {
      if (!is_valid_combo(experiment, window)) {
        cat(sprintf("  [skip] %s x %s: not in EXP_WINDOW_MAP\n", experiment, window))
        next
      }

      fluc_path <- file.path(PROC_DIR, "fluctuations", window, "community_level.csv")
      if (!file.exists(fluc_path)) {
        cat(sprintf("  [skip] missing fluctuation file: %s\n", fluc_path))
        next
      }

      fluc_raw <- read_csv(fluc_path, show_col_types = FALSE) %>%
        filter(experiment == !!experiment, collapsed == FALSE)

      if (nrow(fluc_raw) == 0) {
        cat(sprintf("  [skip] %s x %s: no fluctuation rows\n", experiment, window))
        next
      }

      needed_source_cols <- unique(c(
        sapply(FLUC_METRICS, `[[`, "source_col"),
        "community_cv"
      ))

      fluc_df <- fluc_raw %>%
        select(condition, community, replica, all_of(needed_source_cols))

      for (fluc_metric in names(FLUC_METRICS)) {
        cfg <- FLUC_METRICS[[fluc_metric]]
        fluc_df[[fluc_metric]] <- cfg$transform(fluc_df[[cfg$source_col]])
      }

      div_cache <- list()
      for (fname in unique(sapply(DIVERSITY_METRICS, `[[`, "file"))) {
        div_path <- file.path(PROC_DIR, "diversity", window, fname)
        if (file.exists(div_path)) {
          div_cache[[fname]] <- read_csv(div_path, show_col_types = FALSE) %>%
            filter(experiment == !!experiment, collapsed == FALSE)
        }
      }

      for (fluc_metric in names(FLUC_METRICS)) {
        fluc_cfg <- FLUC_METRICS[[fluc_metric]]

        for (div_metric in names(DIVERSITY_METRICS)) {
          div_cfg <- DIVERSITY_METRICS[[div_metric]]
          div_data <- div_cache[[div_cfg$file]]
          if (is.null(div_data)) next
          if (!div_cfg$source_col %in% colnames(div_data)) next

          div_data_local <- div_data %>%
            mutate(value = div_cfg$transform(.data[[div_cfg$source_col]]))

          df_all <- div_data_local %>%
            select(condition, community, replica, value) %>%
            inner_join(
              fluc_df %>%
                mutate(x_val = .data[[fluc_metric]]) %>%
                select(condition, community, replica, x_val, community_cv),
              by = c("condition", "community", "replica")
            ) %>%
            filter(!is.na(value), !is.na(x_val), is.finite(value), is.finite(x_val)) %>%
            mutate(
              condition = factor(condition, levels = paste0("W", 1:5)),
              fluctuation_class = classify_by_community_cv(community_cv),
              fluctuation_class = factor(
                fluctuation_class,
                levels = c("Stable", "Fluctuation")
              )
            )

          if (nrow(df_all) == 0) next

          xlim_fixed <- get_xlim(fluc_metric, experiment)
          if (!is.null(xlim_fixed)) {
            x_lo <- xlim_fixed[1]
            x_hi <- xlim_fixed[2]
          } else {
            rng <- nice_range_signed(df_all$x_val)
            x_lo <- rng$lo
            x_hi <- rng$hi
          }

          if (!is.null(div_cfg$ylim)) {
            y_lo <- div_cfg$ylim[1]
            y_hi <- div_cfg$ylim[2]
          } else {
            rng <- nice_range_positive(df_all$value)
            y_lo <- rng$lo
            y_hi <- rng$hi
          }

          combo_name <- sprintf("%s_vs_%s", fluc_metric, div_metric)
          fig_dir <- file.path(FIG_DIR, corr_method, experiment, window, combo_name)
          dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

          for (cond in levels(df_all$condition)) {
            df_cond <- df_all %>% filter(condition == cond)
            if (nrow(df_cond) < 3) next

            p <- make_scatter(df_cond, fluc_cfg, corr_method, x_lo, x_hi, y_lo, y_hi)
            ggsave(
              filename = file.path(fig_dir, sprintf("%s_%s.pdf", experiment, cond)),
              plot = p,
              width = W_PLOT,
              height = H_PLOT,
              units = "mm",
              device = cairo_pdf
            )
          }

          if (nrow(df_all) >= 3) {
            p_combined <- make_scatter(
              df_all, fluc_cfg, corr_method, x_lo, x_hi, y_lo, y_hi
            )
            ggsave(
              filename = file.path(fig_dir, sprintf("%s_combined.pdf", experiment)),
              plot = p_combined,
              width = W_PLOT,
              height = H_PLOT,
              units = "mm",
              device = cairo_pdf
            )
          }

          cat(sprintf(
            "  -> %s/%s/%s/%s/ (xlim: [%.2f, %.2f])\n",
            corr_method, experiment, window, combo_name, x_lo, x_hi
          ))
        }
      }
    }
  }
}

cat("\nDone.\n")
cat(sprintf("Output directory: %s\n", FIG_DIR))
cat("Points are not colored by stable/fluctuation class.\n")
cat("Valid experiment-window combinations:\n")
for (exp_name in names(EXP_WINDOW_MAP)) {
  cat(sprintf("  %-12s -> %s\n", exp_name, paste(EXP_WINDOW_MAP[[exp_name]], collapse = ", ")))
}
