community_workspace <- function() {
  override <- Sys.getenv("COMMUNITY_WORKSPACE", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script)) stop("Run with Rscript or set COMMUNITY_WORKSPACE.")
  dirname(normalizePath(sub("^--file=", "", script[[1]]), mustWork = TRUE))
}

# ===========================================================
# plot_summary_24h_temperature.R
# Main directory: /home/hachi/Tem_mortality_workspace
#
# Summary.xlsx contains single-strain measurements across one 24 h cycle.
# This script plots the 24 h OD and CFU across 20, 30, and 40 C.
#
# CFU values in the workbook are exponents n for counts estimated from
# 5 uL. CFU/mL is therefore 10^n * 200, and the plotted value is
# log10(CFU/mL) = n + log10(200).
# ===========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(readxl)
})

PROJECT_DIR <- community_workspace()
DATA_PATH   <- file.path(PROJECT_DIR, "data", "Summary.xlsx")
if (!file.exists(DATA_PATH)) {
  DATA_PATH <- file.path(PROJECT_DIR, "Summary.xlsx")
}
if (!file.exists(DATA_PATH)) {
  stop(sprintf(
    "Summary.xlsx was not found. Checked: %s and %s",
    file.path(PROJECT_DIR, "data", "Summary.xlsx"),
    file.path(PROJECT_DIR, "Summary.xlsx")
  ))
}
FIG_DIR <- file.path(PROJECT_DIR, "figures", "summary_24h_temperature")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Plot parameters -------------------------------------------------------

BG <- 0.035
CFU_VOLUME_FACTOR <- 200

TEMP_LEVELS <- c("20C", "30C", "40C")
TEMP_LABELS <- c("20C" = "20C", "30C" = "30C", "40C" = "40C")
TEMP_COLORS <- c(
  "20C" = "#4E79A7",
  "30C" = "#2E6DA4",
  "40C" = "#C0392B"
)

FONT_FAMILY <- "Arial"
FONT_AX     <- 8
BORDER_SIZE <- 0.25
TICK_SIZE   <- 0.15
TICK_LEN    <- -0.8

W_ST <- 42
H_ST <- 45

LINE_ALPHA <- 0.22
LINE_LW    <- 0.28
JIT_SIZE   <- 1.15
JIT_ALPHA  <- 0.55
JIT_WIDTH  <- 0.08
DOT_SIZE   <- 2.0
ERR_LW     <- 0.45
ERR_W      <- 0.12

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
      axis.ticks.length      = unit(TICK_LEN, "mm"),
      axis.ticks.length.x.top = unit(TICK_LEN, "mm"),
      axis.ticks.length.y.right = unit(TICK_LEN, "mm"),
      axis.text.x.top        = element_blank(),
      axis.text.y.right      = element_blank(),
      axis.text.x            = element_text(size = FONT_AX, color = "black",
                                            margin = margin(t = 3)),
      axis.text.y            = element_text(size = FONT_AX, color = "black",
                                            margin = margin(r = 3)),
      axis.title.x           = element_blank(),
      axis.title.y           = element_blank(),
      legend.position        = "none",
      plot.background        = element_rect(fill = "white", color = NA),
      panel.background       = element_rect(fill = "white", color = NA),
      plot.margin            = margin(3, 3, 3, 3, "mm")
    )
}

fmt_one_decimal <- function(x) {
  formatC(x, format = "f", digits = 1)
}

find_time_column <- function(df, target_time) {
  target_time <- as.character(target_time)
  cleaned_names <- trimws(names(df))
  idx <- which(cleaned_names == target_time)
  if (length(idx) == 0) {
    numeric_names <- suppressWarnings(as.numeric(cleaned_names))
    idx <- which(numeric_names == as.numeric(target_time))
  }
  if (length(idx) == 0) {
    stop(sprintf(
      "Column '%s' not found. Available columns: %s",
      target_time,
      paste(names(df), collapse = ", ")
    ))
  }
  names(df)[idx[1]]
}

finite_max <- function(x, default_value) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(default_value)
  max(x, na.rm = TRUE)
}

finite_min <- function(x, default_value) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(default_value)
  min(x, na.rm = TRUE)
}

read_24h_sheet <- function(temp, metric) {
  sheet_name <- paste0(temp, metric)
  df <- read_excel(DATA_PATH, sheet = sheet_name, .name_repair = "unique")
  names(df)[1:2] <- c("species_id", "species")
  time_col <- find_time_column(df, 24)
  df %>%
    transmute(
      species_id = as.character(species_id),
      species = as.character(species),
      temp = factor(paste0(temp, "C"), levels = TEMP_LEVELS),
      x = match(temp, c(20, 30, 40)),
      raw_24h = suppressWarnings(as.numeric(.data[[time_col]]))
    )
}

read_24h_metric <- function(metric) {
  bind_rows(
    read_24h_sheet(20, metric),
    read_24h_sheet(30, metric),
    read_24h_sheet(40, metric)
  )
}

od_24h <- read_24h_metric("OD") %>%
  mutate(value = pmax(raw_24h - BG, 0)) %>%
  filter(is.finite(value))

cfu_24h <- read_24h_metric("CFU") %>%
  mutate(value = raw_24h + log10(CFU_VOLUME_FACTOR)) %>%
  filter(is.finite(value))

if (nrow(od_24h) == 0) {
  stop("No finite OD values were found at 24 h.")
}
if (nrow(cfu_24h) == 0) {
  stop("No finite CFU values were found at 24 h.")
}

cat(sprintf(
  "Reading data from: %s\n",
  DATA_PATH
))
cat(sprintf(
  "Writing figures to: %s\n",
  FIG_DIR
))
cat(sprintf(
  "OD 24 h: n=%d, range=%.3f-%.3f\n",
  nrow(od_24h),
  finite_min(od_24h$value, NA_real_),
  finite_max(od_24h$value, NA_real_)
))
cat(sprintf(
  "log10(CFU/mL) 24 h: n=%d, range=%.3f-%.3f\n",
  nrow(cfu_24h),
  finite_min(cfu_24h$value, NA_real_),
  finite_max(cfu_24h$value, NA_real_)
))

summarise_temperature <- function(df) {
  df %>%
    group_by(temp, x) %>%
    summarise(
      mu = mean(value, na.rm = TRUE),
      sem = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
      .groups = "drop"
    )
}

save_panel <- function(plot, name) {
  pdf_device <- if (capabilities("cairo")) cairo_pdf else "pdf"
  ggsave(
    file.path(FIG_DIR, paste0(name, ".pdf")),
    plot,
    width = W_ST,
    height = H_ST,
    units = "mm",
    device = pdf_device
  )
  ggsave(
    file.path(FIG_DIR, paste0(name, ".png")),
    plot,
    width = W_ST,
    height = H_ST,
    units = "mm",
    dpi = 600
  )
}

plot_temperature_24h <- function(df, y_limits, y_breaks, y_labels,
                                 mean_shape = 21, connect_species = TRUE) {
  smry <- summarise_temperature(df)
  
  base <- ggplot(df, aes(x = x, y = value))
  if (connect_species) {
    base <- base +
      geom_line(
        aes(group = species_id),
        color = "grey45",
        linewidth = LINE_LW,
        alpha = LINE_ALPHA
      )
  }
  
  base +
    geom_jitter(
      aes(color = temp),
      width = JIT_WIDTH,
      size = JIT_SIZE,
      alpha = JIT_ALPHA,
      shape = 16
    ) +
    geom_errorbar(
      data = smry,
      aes(x = x, ymin = mu - sem, ymax = mu + sem, color = temp),
      width = ERR_W * 2,
      linewidth = ERR_LW,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = smry,
      aes(x = x, y = mu, color = temp),
      shape = mean_shape,
      size = DOT_SIZE,
      fill = "white",
      stroke = 0.7,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = TEMP_COLORS) +
    scale_x_continuous(
      breaks = 1:3,
      labels = TEMP_LABELS[TEMP_LEVELS],
      limits = c(0.6, 3.4),
      expand = expansion(0)
    ) +
    scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      labels = y_labels,
      expand = expansion(mult = c(0, 0))
    ) +
    labs(x = NULL, y = NULL) +
    theme_pub()
}

od_y_hi <- ceiling(finite_max(od_24h$value, 1.6) * 10) / 10
od_y_hi <- max(1.6, od_y_hi)
od_breaks <- seq(0, od_y_hi, length.out = 3)

p_od <- plot_temperature_24h(
  od_24h,
  y_limits = c(0, od_y_hi),
  y_breaks = od_breaks,
  y_labels = fmt_one_decimal(od_breaks),
  mean_shape = 21
)

cfu_y_lo <- floor(finite_min(cfu_24h$value, 3))
cfu_y_hi <- ceiling(finite_max(cfu_24h$value, 10))
if (cfu_y_hi <= cfu_y_lo) cfu_y_hi <- cfu_y_lo + 1
cfu_breaks <- seq(cfu_y_lo, cfu_y_hi, by = 2)
if (tail(cfu_breaks, 1) < cfu_y_hi) cfu_breaks <- c(cfu_breaks, cfu_y_hi)

p_cfu <- plot_temperature_24h(
  cfu_24h,
  y_limits = c(cfu_y_lo, cfu_y_hi),
  y_breaks = cfu_breaks,
  y_labels = function(x) parse(text = paste0("10^", x)),
  mean_shape = 23
)

save_panel(p_od, "OD_24h_by_temperature")
save_panel(p_cfu, "CFU_24h_by_temperature")

summary_out <- bind_rows(
  od_24h %>% mutate(metric = "OD_24h_bg_corrected"),
  cfu_24h %>% mutate(metric = "log10_CFU_per_mL_24h")
) %>%
  select(metric, temp, species_id, species, raw_24h, value)

write_csv(summary_out, file.path(FIG_DIR, "summary_24h_values.csv"))

cat("\nDone: figures/summary_24h_temperature/\n")
cat("  OD_24h_by_temperature.pdf/png\n")
cat("  CFU_24h_by_temperature.pdf/png\n")
cat("  summary_24h_values.csv\n")
