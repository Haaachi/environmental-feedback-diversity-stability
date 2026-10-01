community_workspace <- function() {
  override <- Sys.getenv("COMMUNITY_WORKSPACE", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script)) stop("Run with Rscript or set COMMUNITY_WORKSPACE.")
  dirname(normalizePath(sub("^--file=", "", script[[1]]), mustWork = TRUE))
}

# ===========================================================
# plot_biomass.R
# Original workspace: /home/hachi/Tem_mortality_workspace
#
# Total community biomass across experimental conditions.
# Biomass: mean daily sum of all taxon abs_abund within the window.
# Lines, SEM error bars and hollow points.
#
# Output: figures/biomass/{experiment}/{window}/biomass.pdf
# ===========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

BASE_DIR <- community_workspace()
PROC_DIR <- file.path(BASE_DIR, "processed")
FIG_DIR  <- file.path(BASE_DIR, "figures", "biomass")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Parameters
LAST_DAYS <- list(full = c(8, 9, 10), early = c(4, 5, 6))
USE_REPS  <- list(full = c(1, 2, 3),  early = c(1))

# Dark green distinguishes biomass from blue fluctuation fractions.
LINE_COLOR <- "#2A7F5E"
LINE_WIDTH <- 0.55
ERR_LW     <- 0.40

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
BORDER_SIZE <- 0.25
TICK_SIZE   <- 0.20
TICK_LEN    <- -0.5

W_PLOT <- 40
H_PLOT <- 35

# Publication plot theme
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

# Main loop
for (experiment in c("mortality", "temperature")) {
  for (window in c("full", "early")) {
    if (experiment == "temperature" && window == "early") next
    
    last_days <- LAST_DAYS[[window]]
    use_reps  <- USE_REPS[[window]]
    
    biomass_records <- list()
    
    for (cond in paste0("W", 1:5)) {
      csv_path <- file.path(PROC_DIR, experiment, cond,
                            "abs_abundance.csv")
      if (!file.exists(csv_path)) next
      
      df <- read_csv(csv_path, show_col_types = FALSE) %>%
        filter(day %in% last_days, replica %in% use_reps)
      
      df_biomass <- df %>%
        group_by(condition, community, replica, day) %>%
        summarise(daily_biomass = sum(abs_abund, na.rm = TRUE),
                  .groups = "drop") %>%
        group_by(condition, community, replica) %>%
        summarise(biomass = mean(daily_biomass, na.rm = TRUE),
                  .groups = "drop")
      
      biomass_records[[cond]] <- df_biomass
    }
    
    if (length(biomass_records) == 0) next
    
    df_all <- bind_rows(biomass_records) %>%
      mutate(x_pos = as.integer(sub("W", "", condition)))
    
    summary_df <- df_all %>%
      group_by(condition, x_pos) %>%
      summarise(
        mean_val = mean(biomass, na.rm = TRUE),
        sem_val  = sd(biomass,  na.rm = TRUE) /
          sqrt(sum(!is.na(biomass))),
        .groups  = "drop"
      ) %>%
      mutate(
        ymin = mean_val - sem_val,
        ymax = mean_val + sem_val
      )
    
    # Y-axis
    if (experiment == "temperature") {
      y_upper  <- 1.0
      y_breaks <- c(0, 0.5, 1)
      y_labels <- c("0", "0.5", "1")
    } else {
      # Use fixed Y-axis limits for mortality.
      y_upper  <- 2.0
      y_breaks <- c(0, 1, 2)
      y_labels <- c("0", "1", "2")
    }
    
    # X-axis
    # Show only W1, W3 and W5 on the X-axis.
    x_breaks     <- c(1, 3, 5)
    x_labels_sub <- X_LABELS[[experiment]][c("W1", "W3", "W5")]
    
    p <- ggplot(summary_df, aes(x = x_pos, y = mean_val)) +
      geom_line(color     = LINE_COLOR,
                linewidth = LINE_WIDTH) +
      geom_errorbar(aes(ymin = ymin, ymax = ymax),
                    width     = 0.25,
                    linewidth = ERR_LW,
                    color     = LINE_COLOR) +
      geom_point(shape  = 21,
                 size   = 1.6,
                 color  = LINE_COLOR,
                 fill   = "white",
                 stroke = 0.55) +
      scale_x_continuous(
        breaks   = x_breaks,
        labels   = x_labels_sub,
        expand   = expansion(mult = c(0.12, 0.12)),
        sec.axis = dup_axis(labels = NULL, name = NULL)
      ) +
      scale_y_continuous(
        limits   = c(0, y_upper),
        breaks   = y_breaks,
        labels   = y_labels,
        expand   = expansion(mult = c(0, 0.02)), # No padding below zero.
        sec.axis = dup_axis(labels = NULL, name = NULL)
      ) +
      theme_pub()
    
    fig_dir  <- file.path(FIG_DIR, experiment, window)
    dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
    out_path <- file.path(fig_dir, "biomass.pdf")
    
    ggsave(out_path, p,
           width  = W_PLOT, height = H_PLOT,
           units  = "mm",   device = cairo_pdf)
    
    cat(sprintf("  -> %s/%s/biomass.pdf\n", experiment, window))
  }
}

cat("\nDone!\n")
