community_workspace <- function() {
  override <- Sys.getenv("COMMUNITY_WORKSPACE", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script)) stop("Run with Rscript or set COMMUNITY_WORKSPACE.")
  dirname(normalizePath(sub("^--file=", "", script[[1]]), mustWork = TRUE))
}

# ===========================================================
# plot_taxon_heatmap.R
# Original workspace: /home/hachi/Tem_mortality_workspace
#
# Heatmap of mean abs_abund for each taxon across conditions.
# X-axis: conditions W1-W5.
# Y-axis: taxa ranked by total biomass.
# Color: mean abs_abund across all communities and replicates in a condition.
# Only rel_abund > 1% counts as present; below-detection abundance contributes zero.
#
# Endpoint: full = D10; early = D6.
#
# Output: figures/taxon_heatmap/
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

# Parameters
REL_THRESHOLD <- 1.0   # Count taxa as present only when rel_abund > 1%.

WINDOWS <- list(
  full  = list(day = 10, reps = c(1, 2, 3)),
  early = list(day = 6,  reps = c(1))
)

# X-axis labels
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

# Plot parameters
FONT_FAMILY <- "Arial"
FONT_AX     <- 7
FONT_TAXON  <- 6      # Taxon label font size
BORDER_SIZE <- 0.25
TICK_SIZE   <- 0.20
TICK_LEN    <- -1.0

W_PLOT <- 60    # Plot width in mm for five X-axis conditions
# Compute height dynamically from the number of taxa.

# White-to-dark-blue palette; zero abundance is white.
HEATMAP_LOW  <- "white"
HEATMAP_HIGH <- "#1A4E8A"
HEATMAP_NA   <- "grey95"   # Cells with no measurements

# Publication plot theme
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

# Main loop
for (experiment in c("mortality", "temperature")) {
  for (window_name in c("full", "early")) {
    if (experiment == "temperature" && window_name == "early") next
    
    cfg     <- WINDOWS[[window_name]]
    last_day <- cfg$day
    use_reps <- cfg$reps
    
    cat(sprintf("\n=== %s / %s (D%d, R%s) ===\n",
                experiment, window_name, last_day,
                paste(use_reps, collapse = "+")))
    
    # Read all conditions.
    records <- list()
    
    for (cond in paste0("W", 1:5)) {
      csv_path <- file.path(PROC_DIR, experiment, cond,
                            "abs_abundance.csv")
      if (!file.exists(csv_path)) next
      
      df <- read_csv(csv_path, show_col_types = FALSE) %>%
        filter(day     == last_day,
               replica %in% use_reps)
      
      if (nrow(df) == 0) next
      
      # Retain present records only: rel_abund > 1%.
      # Absent taxa contribute zero abundance, with no term in the numerator.
      df_alive <- df %>%
        filter(rel_abund > REL_THRESHOLD) %>%
        mutate(condition = cond)
      
      records[[cond]] <- df_alive
    }
    
    if (length(records) == 0) next
    
    df_all <- bind_rows(records)
    
    # Compute mean abs_abund for each taxon and condition.
    # Denominator: all communities and replicates, including absent observations.
    # mean = sum(abs_abund_alive) / n_total_samples
    # Absent observations contribute zero to the true condition mean.
    
    # Total sample count per condition: communities x replicates.
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
    
    # Order taxa by decreasing total biomass.
    taxon_order <- heatmap_df %>%
      group_by(taxon) %>%
      summarise(total = sum(mean_abs), .groups = "drop") %>%
      arrange(desc(total)) %>%
      pull(taxon)
    
    # Fill cells where a taxon was absent throughout a condition.
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
    cat(sprintf("  Present taxa: %d\n", n_taxa))
    
    # Dynamic height: 3 mm per taxon plus top/bottom margins.
    H_PLOT <- n_taxa * 3 + 10
    
    # Draw the heatmap.
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
    
    # Remove the duplicate upper secondary axis.
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

cat("\nDone!\n")
cat("figures/taxon_heatmap/\n")
cat("  mortality_full_taxon_heatmap.pdf\n")
cat("  mortality_early_taxon_heatmap.pdf\n")
cat("  temperature_full_taxon_heatmap.pdf\n")
cat(sprintf("\n  Presence criterion: rel_abund > %.1f%%\n", REL_THRESHOLD))
cat("  Color: mean abs_abund (absent observations contribute zero)\n")
cat("  Y-axis: taxa ordered by decreasing total biomass\n")
