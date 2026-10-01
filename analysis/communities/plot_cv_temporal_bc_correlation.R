community_workspace <- function() {
  override <- Sys.getenv("COMMUNITY_WORKSPACE", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script)) stop("Run with Rscript or set COMMUNITY_WORKSPACE.")
  dirname(normalizePath(sub("^--file=", "", script[[1]]), mustWork = TRUE))
}

# ===========================================================
# plot_cv_temporal_bc_correlation.R
#
# Diagnostic correlation plots between the two current species-instability
# axes:
#   x = sum_abs_std      = sum_i SD(N_i)
#   y = temporal_bc_mean = mean temporal Bray-Curtis distance
#
# Thresholds are configurable below. Current working thresholds:
#   sum_abs_std: mortality = 0.20, temperature = 0.15
#   temporal_bc_mean = 0.17
# ===========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

BASE_DIR <- community_workspace()
PROC_DIR <- file.path(BASE_DIR, "processed")
FIG_DIR  <- file.path(BASE_DIR, "figures", "sum_abs_std_temporal_bc_correlation")
OUT_DIR  <- file.path(PROC_DIR, "sum_abs_std_temporal_bc_correlation")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

MAIN_WINDOWS <- tribble(
  ~experiment,   ~window,
  "mortality",   "early",
  "temperature", "full"
)

WINDOWS <- c("early", "full", "last4")

SUM_ABS_STD_THRESHOLDS <- c(
  mortality   = 0.33,
  temperature = 0.2
)
TEMPORAL_BC_THRESHOLD <- 0.16

SUM_ABS_STD_X_LIMITS <- list(
  mortality   = c(0, 1.6),
  temperature = c(0, 0.8)
)

sum_abs_std_threshold_for <- function(experiment) {
  threshold <- unname(SUM_ABS_STD_THRESHOLDS[as.character(experiment)])
  threshold[is.na(threshold)] <- NA_real_
  threshold
}

CLASS_COLORS <- c(
  "Stable core" = "#9B8EC4",
  "Sum SD high" = "#E6A157",
  "BC high"     = "#4C78A8",
  "Both high"   = "#55A868"
)
LINE_COL <- "#3F3F3F"

FONT_FAMILY <- "Arial"
FONT_AX     <- 10
BORDER_SIZE <- 0.25
TICK_SIZE   <- 0.20
TICK_LEN    <- -1.0
DOT_SIZE    <- 1.6
DOT_ALPHA   <- 0.78
LINE_SIZE   <- 0.28

W_SINGLE <- 110
H_SINGLE <- 92
W_FACET  <- 210
H_FACET  <- 96

theme_pub_axis <- function() {
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
      axis.text              = element_text(size = FONT_AX, color = "black"),
      axis.title             = element_text(size = FONT_AX, color = "black"),
      aspect.ratio           = 1,
      strip.background       = element_blank(),
      strip.text             = element_text(size = FONT_AX, color = "black"),
      legend.position        = "right",
      legend.title           = element_blank(),
      legend.text            = element_text(size = FONT_AX - 1,
                                            color = "black"),
      plot.background        = element_rect(fill = "white", color = NA),
      panel.background       = element_rect(fill = "white", color = NA)
    )
}

read_fluctuation_window <- function(window) {
  path <- file.path(PROC_DIR, "fluctuations", window, "community_level.csv")
  if (!file.exists(path)) return(tibble())

  read_csv(path, show_col_types = FALSE) %>%
    mutate(window = window)
}

clean_fluctuation_df <- function(df) {
  df %>%
    mutate(
      sum_abs_std = as.numeric(sum_abs_std),
      temporal_bc_mean = as.numeric(temporal_bc_mean),
      collapsed = if_else(is.na(collapsed), FALSE, collapsed),
      sum_abs_std_threshold = map_dbl(
        experiment, sum_abs_std_threshold_for
      ),
      temporal_bc_threshold = TEMPORAL_BC_THRESHOLD,
      sum_abs_std_high = sum_abs_std >= sum_abs_std_threshold,
      temporal_bc_high = temporal_bc_mean >= TEMPORAL_BC_THRESHOLD,
      species_instability_class = case_when(
        sum_abs_std_high & temporal_bc_high ~ "Both high",
        sum_abs_std_high                    ~ "Sum SD high",
        temporal_bc_high                    ~ "BC high",
        TRUE                                ~ "Stable core"
      ),
      species_instability_class = factor(
        species_instability_class,
        levels = c("Stable core", "Sum SD high", "BC high", "Both high")
      )
    ) %>%
    filter(
      is.finite(sum_abs_std),
      is.finite(temporal_bc_mean),
      !(experiment == "temperature" & window == "early")
    )
}

axis_limits <- function(x, floor_zero = TRUE) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(0, 1))
  hi <- max(x, na.rm = TRUE)
  if (!is.finite(hi) || hi <= 0) hi <- 1
  pretty_hi <- max(pretty(c(0, hi * 1.05), n = 5), na.rm = TRUE)
  if (floor_zero) c(0, pretty_hi) else range(x, na.rm = TRUE)
}

fmt_p <- function(p) {
  if (is.na(p)) return("p=NA")
  if (p < 0.001) return("p<0.001")
  sprintf("p=%.3f", p)
}

spearman_result <- function(df) {
  d <- df %>%
    filter(is.finite(sum_abs_std), is.finite(temporal_bc_mean))
  if (nrow(d) < 3) {
    return(tibble(n = nrow(d), rho = NA_real_, p_value = NA_real_))
  }

  res <- suppressWarnings(
    cor.test(d$sum_abs_std, d$temporal_bc_mean,
             method = "spearman", exact = FALSE)
  )
  tibble(
    n = nrow(d),
    rho = as.numeric(res$estimate),
    p_value = res$p.value
  )
}

corr_label <- function(df) {
  res <- spearman_result(df)
  sprintf("\u03c1=%.2f, %s, n=%d",
          res$rho, fmt_p(res$p_value), res$n)
}

x_limits_for_experiment <- function(experiment) {
  limits <- SUM_ABS_STD_X_LIMITS[[as.character(experiment)[1]]]
  if (is.null(limits)) c(0, 0.8) else limits
}

make_annotation <- function(df, x_limits, y_limits, facet = FALSE) {
  if (facet) {
    df %>%
      group_by(experiment) %>%
      group_modify(~ tibble(label = corr_label(.x))) %>%
      ungroup() %>%
      rowwise() %>%
      mutate(x = x_limits_for_experiment(experiment)[2],
             y = y_limits[2]) %>%
      ungroup()
  } else {
    tibble(
      x = x_limits[2],
      y = y_limits[2],
      label = corr_label(df)
    )
  }
}

plot_sum_std_temporal_bc <- function(df, x_limits, y_limits, facet = FALSE) {
  if (nrow(df) == 0) return(NULL)

  ann <- make_annotation(df, x_limits, y_limits, facet = facet)
  facet_limits <- if (facet) {
    bind_rows(lapply(unique(df$experiment), function(exp_i) {
      limits <- x_limits_for_experiment(exp_i)
      tibble(experiment = exp_i, x = limits, y = y_limits[1])
    }))
  } else {
    tibble(experiment = character(), x = numeric(), y = numeric())
  }
  threshold_lines <- df %>%
    distinct(experiment, sum_abs_std_threshold, temporal_bc_threshold) %>%
    filter(is.finite(sum_abs_std_threshold),
           is.finite(temporal_bc_threshold))

  p <- ggplot(df, aes(x = sum_abs_std, y = temporal_bc_mean)) +
    geom_blank(data = facet_limits,
               aes(x = x, y = y),
               inherit.aes = FALSE) +
    geom_vline(data = threshold_lines,
               aes(xintercept = sum_abs_std_threshold),
               inherit.aes = FALSE,
               linetype = "dashed", linewidth = LINE_SIZE,
               color = "grey45") +
    geom_hline(data = threshold_lines,
               aes(yintercept = temporal_bc_threshold),
               inherit.aes = FALSE,
               linetype = "dashed", linewidth = LINE_SIZE,
               color = "grey45") +
    geom_point(aes(color = species_instability_class),
               size = DOT_SIZE, alpha = DOT_ALPHA) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                linewidth = LINE_SIZE, color = LINE_COL) +
    geom_text(
      data = ann,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1.02, vjust = 1.25,
      size = (FONT_AX - 1) / .pt,
      family = FONT_FAMILY,
      color = "grey20"
    ) +
    scale_color_manual(values = CLASS_COLORS, drop = FALSE) +
    scale_x_continuous(
      limits = if (facet) NULL else x_limits,
      breaks = pretty_breaks(5),
      labels = label_number(accuracy = 0.1),
      expand = expansion(mult = c(0, 0)),
      sec.axis = dup_axis(labels = NULL, name = NULL)
    ) +
    scale_y_continuous(
      limits = y_limits,
      breaks = pretty_breaks(5),
      labels = label_number(accuracy = 0.1),
      expand = expansion(mult = c(0, 0)),
      sec.axis = dup_axis(labels = NULL, name = NULL)
    ) +
    labs(x = "Sum species SD", y = "Temporal BC") +
    coord_cartesian(xlim = if (facet) NULL else x_limits,
                    ylim = y_limits, clip = "on") +
    theme_pub_axis()

  if (facet) {
    p <- p + facet_wrap(~ experiment, nrow = 1, scales = "free_x")
  }

  p
}

all_df <- bind_rows(lapply(WINDOWS, read_fluctuation_window)) %>%
  clean_fluctuation_df()

if (nrow(all_df) == 0) {
  stop("No valid fluctuation rows with sum_abs_std and temporal_bc_mean found.")
}

main_df <- all_df %>%
  semi_join(MAIN_WINDOWS, by = c("experiment", "window"))

x_limits <- c(0, 0.8)
y_limits <- axis_limits(all_df$temporal_bc_mean)
y_limits[2] <- max(y_limits[2], TEMPORAL_BC_THRESHOLD * 1.15)

correlation_results <- bind_rows(
  main_df %>%
    group_by(experiment) %>%
    group_modify(~ spearman_result(.x)) %>%
    ungroup() %>%
    mutate(scope = "main", window = "main"),
  all_df %>%
    group_by(window, experiment) %>%
    group_modify(~ spearman_result(.x)) %>%
    ungroup() %>%
    mutate(scope = "window")
) %>%
  select(scope, window, experiment, n, rho, p_value)

write_csv(
  correlation_results,
  file.path(OUT_DIR, "sum_abs_std_temporal_bc_spearman.csv")
)

write_csv(
  main_df %>%
    arrange(experiment, condition, community, replica) %>%
    select(window, experiment, condition, community, replica,
           sum_abs_std, temporal_bc_mean, community_cv,
           total_biomass_std, total_biomass_cv,
           sum_abs_std_high, temporal_bc_high,
           sum_abs_std_threshold, temporal_bc_threshold,
           species_instability_class, mean_total_abs, collapsed),
  file.path(OUT_DIR, "main_sum_abs_std_temporal_bc_values.csv")
)

classification_summary <- all_df %>%
  count(window, experiment, condition, species_instability_class,
        name = "n") %>%
  arrange(window, experiment, condition, species_instability_class)

write_csv(
  classification_summary,
  file.path(OUT_DIR, "species_instability_class_summary.csv")
)

p_main <- plot_sum_std_temporal_bc(main_df, x_limits, y_limits, facet = TRUE)
if (!is.null(p_main)) {
  ggsave(
    file.path(FIG_DIR, "main_sum_abs_std_temporal_bc_by_experiment.pdf"),
    p_main, width = W_FACET, height = H_FACET, units = "mm",
    device = cairo_pdf
  )
}

p_active <- plot_sum_std_temporal_bc(
  main_df %>% filter(collapsed == FALSE),
  x_limits, y_limits, facet = TRUE
)
if (!is.null(p_active)) {
  ggsave(
    file.path(FIG_DIR, "main_sum_abs_std_temporal_bc_active_only_by_experiment.pdf"),
    p_active, width = W_FACET, height = H_FACET, units = "mm",
    device = cairo_pdf
  )
}

for (i in seq_len(nrow(MAIN_WINDOWS))) {
  exp_i <- MAIN_WINDOWS$experiment[i]
  win_i <- MAIN_WINDOWS$window[i]
  d <- main_df %>% filter(experiment == exp_i, window == win_i)
  p <- plot_sum_std_temporal_bc(
    d, x_limits_for_experiment(exp_i), y_limits, facet = FALSE
  )
  if (is.null(p)) next

  ggsave(
    file.path(FIG_DIR, sprintf("%s_%s_sum_abs_std_temporal_bc.pdf",
                              exp_i, win_i)),
    p, width = W_SINGLE, height = H_SINGLE, units = "mm",
    device = cairo_pdf
  )
}

last4_df <- all_df %>% filter(window == "last4")
p_last4 <- plot_sum_std_temporal_bc(last4_df, x_limits, y_limits, facet = TRUE)
if (!is.null(p_last4)) {
  ggsave(
    file.path(FIG_DIR, "last4_sum_abs_std_temporal_bc_by_experiment.pdf"),
    p_last4, width = W_FACET, height = H_FACET, units = "mm",
    device = cairo_pdf
  )
}

cat("Done.\n")
cat(sprintf("Figures: %s\n", FIG_DIR))
cat(sprintf("Tables:  %s\n", OUT_DIR))
