community_workspace <- function() {
  override <- Sys.getenv("COMMUNITY_WORKSPACE", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script)) stop("Run with Rscript or set COMMUNITY_WORKSPACE.")
  dirname(normalizePath(sub("^--file=", "", script[[1]]), mustWork = TRUE))
}

# ===========================================================
# plot_diversity_distribution_fluctuating.R
# Minimal Shannon / accumulated Shannon distributions for fluctuating communities.
# ===========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

BASE_DIR <- community_workspace()
PROC_DIR <- file.path(BASE_DIR, "processed")
FIG_DIR  <- file.path(BASE_DIR, "figures", "diversity_distribution")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

experiment <- "mortality"
COMMUNITY_CV_THRESHOLD <- 0.25

is_fluctuating_composite <- function(experiment, community_cv) {
  !is.na(community_cv) & community_cv >= COMMUNITY_CV_THRESHOLD
}

for (window in c("early", "full")) {
  
  fluc_path <- file.path(PROC_DIR, "fluctuations", window, "community_level.csv")
  alpha_path <- file.path(PROC_DIR, "diversity", window, "alpha_diversity.csv")
  gamma_path <- file.path(PROC_DIR, "diversity", window, "gamma_diversity.csv")
  
  fluc_df <- read_csv(fluc_path, show_col_types = FALSE) %>%
    filter(experiment == !!experiment, collapsed == FALSE) %>%
    select(condition, community, replica, community_cv) %>%
    mutate(fluctuating = is_fluctuating_composite(experiment, community_cv))
  
  alpha_div <- read_csv(alpha_path, show_col_types = FALSE) %>%
    filter(experiment == !!experiment, collapsed == FALSE) %>%
    select(condition, community, replica, shannon) %>%
    inner_join(fluc_df, by = c("condition", "community", "replica")) %>%
    filter(fluctuating)
  
  gamma_div <- read_csv(gamma_path, show_col_types = FALSE) %>%
    filter(experiment == !!experiment, collapsed == FALSE) %>%
    select(condition, community, replica, gamma_shannon) %>%
    inner_join(fluc_df, by = c("condition", "community", "replica")) %>%
    filter(fluctuating)
  
  # Shannon
  p1 <- ggplot(alpha_div, aes(x = shannon)) +
    geom_histogram(bins = 15, fill = "#4A7FB5", alpha = 0.6, boundary = 0) +
    geom_density(alpha = 0.4, linewidth = 0.5, color = "#2A5F8F") +
    theme_void()
  
  ggsave(file.path(FIG_DIR, sprintf("shannon_fluctuating_%s_%s.pdf", experiment, window)),
         p1, width = 40, height = 25, units = "mm", device = cairo_pdf)
  
  # Gamma Shannon
  p2 <- ggplot(gamma_div, aes(x = gamma_shannon)) +
    geom_histogram(bins = 15, fill = "#C75B4E", alpha = 0.6, boundary = 0) +
    geom_density(alpha = 0.4, linewidth = 0.5, color = "#A04030") +
    theme_void()
  
  ggsave(file.path(FIG_DIR, sprintf("gamma_shannon_fluctuating_%s_%s.pdf", experiment, window)),
         p2, width = 40, height = 25, units = "mm", device = cairo_pdf)
  
  cat(sprintf("%s: Shannon n=%d | Gamma n=%d\n", window, nrow(alpha_div), nrow(gamma_div)))
}

cat("Done\n")
