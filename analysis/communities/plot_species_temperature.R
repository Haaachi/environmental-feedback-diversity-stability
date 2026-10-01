community_workspace <- function() {
  override <- Sys.getenv("COMMUNITY_WORKSPACE", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script)) stop("Run with Rscript or set COMMUNITY_WORKSPACE.")
  dirname(normalizePath(sub("^--file=", "", script[[1]]), mustWork = TRUE))
}

# ===========================================================
# plot_species_temperature.R
# Original workspace: /home/hachi/Tem_mortality_workspace
#
# Figure 1a: OD trajectories at 30 degrees C, translucent isolates and mean +/- SEM.
# Figure 1b: OD trajectories at 40 degrees C.
# Figure 2a: OD_max comparisons with individual points and automatic significance.
# Figure 2b: log10(CFU_max) comparisons with individual points and significance.
# Figure 3: ensemble CV.
#
# Statistical test selection:
# Shapiro-Wilk p > 0.05: Welch t test.
# Shapiro-Wilk p <= 0.05: Mann-Whitney U test.
# Subtract OD background 0.035 and clip at zero.
# CFU 0/NA → 1e6 CFU/mL
# ===========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(readxl)
})

BASE_DIR <- community_workspace()
DATA_PATH <- file.path(BASE_DIR, "data", "Speciestemperature.xlsx")
FIG_DIR   <- file.path(BASE_DIR, "figures", "species_temperature")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Parameters
BG       <- 0.035
CFU_ZERO <- 1e6

COL_30  <- "#2E6DA4"
COL_40  <- "#C0392B"
COL_CFU <- "#E07B39"

FONT_FAMILY <- "Arial"
FONT_AX     <- 8
BORDER_SIZE <- 0.25
TICK_SIZE   <- 0.15
TICK_LEN    <- -0.8

W_TS <- 55; H_TS <- 45
W_ST <- 30; H_ST <- 45
W_CV <- 40; H_CV <- 40

LINE_ALPHA  <- 0.25
LINE_LW     <- 0.30
MEAN_LW     <- 0.80
DOT_SIZE    <- 2.0
JIT_SIZE    <- 1.0      
JIT_ALPHA   <- 0.45     
JIT_WIDTH   <- 0.08     
ERR_LW      <- 0.45
ERR_W       <- 0.10

# Publication plot theme
theme_pub <- function() {
  theme_classic(base_size = FONT_AX,
                base_family = FONT_FAMILY) +
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
      axis.text.x  = element_text(size   = FONT_AX, color = "black",
                                  margin = margin(t = 3)),
      axis.text.y  = element_text(size   = FONT_AX, color = "black",
                                  margin = margin(r = 3)),
      axis.title.x  = element_blank(),
      axis.title.y  = element_blank(),  # Remove the Y-axis title area together with labs(NULL).
      legend.position   = "none",
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA),
      plot.margin       = margin(3, 3, 3, 3, "mm")
    )
}

# Automatic statistical test selection
auto_sig_label <- function(x1, x2) {
  x1 <- x1[!is.na(x1)]; x2 <- x2[!is.na(x2)]
  cat(sprintf("  n = %d vs %d\n", length(x1), length(x2)))
  
  normal1 <- if (length(x1) >= 3) shapiro.test(x1)$p.value > 0.05 else FALSE
  normal2 <- if (length(x2) >= 3) shapiro.test(x2)$p.value > 0.05 else FALSE
  
  if (normal1 && normal2) {
    test_result <- t.test(x1, x2, var.equal = FALSE)
    cat("  Test: Welch t test\n")
  } else {
    test_result <- wilcox.test(x1, x2, exact = FALSE)
    cat("  Test: Mann-Whitney U\n")
  }
  
  p <- test_result$p.value
  cat(sprintf("  p = %.4f\n", p))
  
  if      (p > 0.05)  "ns"
  else if (p > 0.01)  "*"
  else if (p > 0.001) "**"
  else                "***"
}

# Significance annotation helper
add_sig_annotation <- function(plt, label, y_bar, y_top,
                               x1 = 1, x2 = 2,
                               tick_h = NULL,
                               font_size = FONT_AX) {
  tick_h  <- tick_h %||% (y_top * 0.03)
  seg_df  <- data.frame(
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

# Read data.
day_cols <- paste0("Day", 1:5)

od30 <- read_excel(DATA_PATH, sheet = "30OD") %>%
  rename(species = 1) %>%
  mutate(across(all_of(day_cols), as.numeric)) %>%
  mutate(across(all_of(day_cols), ~pmax(. - BG, 0)))

od40 <- read_excel(DATA_PATH, sheet = "40OD") %>%
  rename(species = 1) %>%
  mutate(across(all_of(day_cols), as.numeric)) %>%
  mutate(across(all_of(day_cols), ~pmax(. - BG, 0)))

cfu30 <- read_excel(DATA_PATH, sheet = "30CFU") %>%
  rename(species = 1) %>%
  mutate(across(all_of(day_cols),
                ~{v <- suppressWarnings(as.numeric(.)) * 100
                if_else(is.na(v) | v == 0, CFU_ZERO, v)}))

cfu40 <- read_excel(DATA_PATH, sheet = "40CFU") %>%
  rename(species = 1) %>%
  mutate(across(c("Day3","Day4","Day5"),
                ~{v <- suppressWarnings(as.numeric(.)) * 100
                if_else(is.na(v) | v == 0, CFU_ZERO, as.numeric(v))}))

# Figure 1: OD time series
plot_timeseries <- function(od_df, color) {
  df_long <- od_df %>%
    pivot_longer(all_of(day_cols), names_to = "day", values_to = "OD") %>%
    mutate(day = as.integer(str_extract(day, "\\d+")))
  
  smry <- df_long %>%
    group_by(day) %>%
    summarise(mu  = mean(OD, na.rm = TRUE),
              sem = sd(OD, na.rm = TRUE) / sqrt(sum(!is.na(OD))),
              .groups = "drop")
  
  y_br <- pretty(c(0, max(df_long$OD, na.rm = TRUE) * 1.05), n = 3)
  
  ggplot() +
    geom_line(data = df_long, aes(x = day, y = OD, group = species),
              color = color, linewidth = LINE_LW, alpha = LINE_ALPHA) +
    geom_line(data = smry, aes(x = day, y = mu),
              color = color, linewidth = MEAN_LW) +
    geom_errorbar(data = smry, aes(x = day, ymin = mu - sem, ymax = mu + sem),
                  width = ERR_W, linewidth = ERR_LW, color = color) +
    geom_point(data = smry, aes(x = day, y = mu),
               shape = 21, size = DOT_SIZE,
               color = color, fill = "white", stroke = 0.6) +
    scale_x_continuous(breaks = c(1, 3, 5), limits = c(0.8, 5.2),
                       expand = expansion(0),
                       sec.axis = dup_axis(labels = NULL, name = NULL)) +
    scale_y_continuous(limits = c(0, max(y_br)), breaks = y_br,
                       expand = expansion(mult = c(0, 0)),
                       sec.axis = dup_axis(labels = NULL, name = NULL)) +
    labs(x = NULL, y = NULL) +   # Remove axis titles.
    theme_pub()
}

p1a <- plot_timeseries(od30, COL_30)
p1b <- plot_timeseries(od40, COL_40)
ggsave(file.path(FIG_DIR, "OD_timeseries_30.pdf"), p1a,
       width = W_TS, height = H_TS, units = "mm", device = cairo_pdf)
ggsave(file.path(FIG_DIR, "OD_timeseries_40.pdf"), p1b,
       width = W_TS, height = H_TS, units = "mm", device = cairo_pdf)
cat("-> OD_timeseries_30/40.pdf\n")

# Figure 2a: OD_max comparisons
odmax30 <- od30 %>% rowwise() %>%
  mutate(v = max(c_across(all_of(day_cols)), na.rm = TRUE)) %>% pull(v)
odmax40 <- od40 %>% rowwise() %>%
  mutate(v = max(c_across(all_of(day_cols)), na.rm = TRUE)) %>% pull(v)

stat_od <- tibble(
  temp = factor(c("30°C", "40°C"), levels = c("30°C", "40°C")),
  x    = c(1, 2),
  mu   = c(mean(odmax30, na.rm = TRUE), mean(odmax40, na.rm = TRUE)),
  sem  = c(sd(odmax30, na.rm = TRUE) / sqrt(sum(!is.na(odmax30))),
           sd(odmax40, na.rm = TRUE) / sqrt(sum(!is.na(odmax40))))
)

jit_od <- bind_rows(
  tibble(x = 1, val = odmax30, temp = "30°C"),
  tibble(x = 2, val = odmax40, temp = "40°C")
) %>% mutate(temp = factor(temp, levels = c("30°C", "40°C")))

cat("OD_max significance test:\n")
sig_od <- auto_sig_label(odmax30, odmax40)
cat(sprintf("  Annotation:%s\n", sig_od))
cat(sprintf("  Means  30: %.3f±%.3f  40: %.3f±%.3f\n",
            stat_od$mu[1], stat_od$sem[1],
            stat_od$mu[2], stat_od$sem[2]))

od_y_bar <- 0.92
od_y_top <- 1.0

p2a <- ggplot() +
  geom_jitter(data = jit_od,
              aes(x = x, y = val, color = temp),
              width = JIT_WIDTH, size = JIT_SIZE, alpha = JIT_ALPHA,
              shape = 16) +
  geom_errorbar(data = stat_od,
                aes(x = x, ymin = mu - sem, ymax = mu + sem, color = temp),
                width = ERR_W * 2, linewidth = ERR_LW) +
  geom_point(data = stat_od,
             aes(x = x, y = mu, color = temp),
             shape = 21, size = DOT_SIZE, fill = "white", stroke = 0.7) +
  scale_color_manual(values = c("30°C" = COL_30, "40°C" = COL_40)) +
  scale_x_continuous(breaks = c(1, 2), labels = c("30°C", "40°C"),
                     limits = c(0.6, 2.4), expand = expansion(0)) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0,0.5, 1.0),
    expand = expansion(mult = c(0, 0)),
    sec.axis = dup_axis(labels = NULL, name = NULL)
  ) +
  labs(x = NULL, y = NULL) +   # Remove axis titles.
  theme_pub()

p2a <- add_sig_annotation(p2a, sig_od, y_bar = od_y_bar, y_top = od_y_top)

ggsave(file.path(FIG_DIR, "stat_OD.pdf"), p2a,
       width = W_ST, height = H_ST, units = "mm", device = cairo_pdf)
cat("-> stat_OD.pdf\n")

# Figure 2b: log10(CFU_max) comparisons
log_cfumax30 <- cfu30 %>% rowwise() %>%
  mutate(v = log10(max(c_across(c("Day3","Day4","Day5")), na.rm = TRUE))) %>%
  pull(v)
log_cfumax40 <- cfu40 %>% rowwise() %>%
  mutate(v = log10(max(c_across(c("Day3","Day4","Day5")), na.rm = TRUE))) %>%
  pull(v)

stat_cfu <- tibble(
  temp = factor(c("30°C", "40°C"), levels = c("30°C", "40°C")),
  x    = c(1, 2),
  mu   = c(mean(log_cfumax30, na.rm = TRUE), mean(log_cfumax40, na.rm = TRUE)),
  sem  = c(sd(log_cfumax30, na.rm = TRUE) / sqrt(sum(!is.na(log_cfumax30))),
           sd(log_cfumax40, na.rm = TRUE) / sqrt(sum(!is.na(log_cfumax40))))
)

jit_cfu <- bind_rows(
  tibble(x = 1, val = log_cfumax30, temp = "30°C"),
  tibble(x = 2, val = log_cfumax40, temp = "40°C")
) %>% mutate(temp = factor(temp, levels = c("30°C", "40°C")))

cat("CFU_max significance test:\n")
sig_cfu <- auto_sig_label(log_cfumax30, log_cfumax40)
cat(sprintf("  Annotation:%s\n", sig_cfu))
cat(sprintf("  logCFU  30: %.2f±%.2f  40: %.2f±%.2f\n",
            stat_cfu$mu[1], stat_cfu$sem[1],
            stat_cfu$mu[2], stat_cfu$sem[2]))

cfu_lo <- 5.5; cfu_hi <- 9
cfu_y_bar <- cfu_hi - 0.3
cfu_y_top <- cfu_hi

p2b <- ggplot() +
  geom_jitter(data = jit_cfu,
              aes(x = x, y = val, color = temp),
              width = JIT_WIDTH, size = JIT_SIZE, alpha = JIT_ALPHA,
              shape = 16) +
  geom_errorbar(data = stat_cfu,
                aes(x = x, ymin = mu - sem, ymax = mu + sem, color = temp),
                width = ERR_W * 2, linewidth = ERR_LW) +
  geom_point(data = stat_cfu,
             aes(x = x, y = mu, color = temp),
             shape = 23, size = DOT_SIZE,
             fill = "white", stroke = 0.7) +
  scale_color_manual(values = c("30°C" = COL_30, "40°C" = COL_40)) +
  scale_x_continuous(breaks = c(1, 2), labels = c("30°C", "40°C"),
                     limits = c(0.6, 2.4), expand = expansion(0)) +
  scale_y_continuous(
    limits = c(cfu_lo, cfu_hi),
    breaks = seq(6, cfu_hi, by = 1),  # Use positive CFU ticks; fixes an earlier seq(-6, ...) typo.
    labels = function(x) parse(text = paste0("10^", x)),
    expand = expansion(mult = c(0, 0)),
    sec.axis = dup_axis(labels = NULL, name = NULL)
  ) +
  labs(x = NULL, y = NULL) +   # Remove axis titles.
  theme_pub()

p2b <- add_sig_annotation(p2b, sig_cfu,
                          y_bar = cfu_y_bar, y_top = cfu_y_top,
                          tick_h = 0.08)

ggsave(file.path(FIG_DIR, "stat_CFU.pdf"), p2b,
       width = W_ST, height = H_ST, units = "mm", device = cairo_pdf)
cat("-> stat_CFU.pdf\n")

# Figure 3: ensemble CV
calc_cv <- function(od_df) {
  std_i  <- od_df %>% rowwise() %>%
    mutate(s = sd(c_across(all_of(day_cols)), na.rm = TRUE)) %>% pull(s)
  mean_i <- od_df %>% rowwise() %>%
    mutate(m = mean(c_across(all_of(day_cols)), na.rm = TRUE)) %>% pull(m)
  mean(std_i, na.rm = TRUE) / mean(mean_i, na.rm = TRUE)
}

cv30 <- calc_cv(od30)
cv40 <- calc_cv(od40)
cat(sprintf("Ensemble CV 30°C: %.4f\n", cv30))
cat(sprintf("Ensemble CV 40°C: %.4f\n", cv40))

cv_df <- tibble(
  temp = factor(c("30°C", "40°C"), levels = c("30°C", "40°C")),
  x    = c(1, 2),
  cv   = c(cv30, cv40)
)

y_br_cv <- pretty(c(0, max(cv_df$cv) * 1.2), n = 3)

p3 <- ggplot(cv_df, aes(x = x, y = cv)) +
  geom_point(aes(color = temp), shape = 16, size = DOT_SIZE + 0.5) +
  scale_color_manual(values = c("30°C" = COL_30, "40°C" = COL_40)) +
  scale_x_continuous(breaks = c(1, 2), labels = c("30°C", "40°C"),
                     limits = c(0.6, 2.4), expand = expansion(0)) +
  scale_y_continuous(limits = c(0, max(y_br_cv)), breaks = y_br_cv,
                     expand = expansion(mult = c(0, 0)),
                     sec.axis = dup_axis(labels = NULL, name = NULL)) +
  labs(x = NULL, y = NULL) +   # Remove axis titles.
  theme_pub()

ggsave(file.path(FIG_DIR, "ensemble_cv.pdf"), p3,
       width = W_CV, height = H_CV, units = "mm", device = cairo_pdf)
cat("-> ensemble_cv.pdf\n")

cat("\nDone!figures/species_temperature/\n")
cat("  OD_timeseries_30.pdf\n")
cat("  OD_timeseries_40.pdf\n")
cat("  stat_OD.pdf   individual points + mean +/- SEM + automatic significance (no axis titles)\n")
cat("  stat_CFU.pdf  individual points + mean +/- SEM + automatic significance (no axis titles)\n")
cat("  ensemble_cv.pdf\n")
