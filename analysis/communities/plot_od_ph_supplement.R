community_workspace <- function() {
  override <- Sys.getenv("COMMUNITY_WORKSPACE", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script)) stop("Run with Rscript or set COMMUNITY_WORKSPACE.")
  dirname(normalizePath(sub("^--file=", "", script[[1]]), mustWork = TRUE))
}

# ===========================================================
# plot_od_ph_supplement.R
#
# Supplementary OD/pH trajectory and replicate-divergence figures.
# This script follows the existing compact publication style, but keeps
# axis labels because these are supplementary diagnostic figures.
# ===========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(readxl)
  library(scales)
})

BASE_DIR <- community_workspace()
DATA_DIR <- file.path(BASE_DIR, "data")
PROC_DIR <- file.path(BASE_DIR, "processed")
FIG_DIR  <- file.path(BASE_DIR, "figures", "supplement_od_ph")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

OD_FILES <- list(
  mortality   = file.path(DATA_DIR, "Mortality_OD_Unified.xlsx"),
  temperature = file.path(DATA_DIR, "Temperature_OD_Unified.xlsx")
)

PH_FILES <- list(
  mortality   = file.path(DATA_DIR, "Mortality_pH_Unified.xlsx"),
  temperature = file.path(DATA_DIR, "Temperature_pH_Unified.xlsx")
)

OD_BACKGROUND <- c(mortality = 0.040, temperature = 0.035)
OD_LIMITS <- list(
  mortality   = c(0, 2.4),
  temperature = c(0, 1.2)
)
PH_LIMITS <- c(2, 10)
MAX_DAY <- c(mortality = 6, temperature = 10)
LAST3_DAYS <- list(mortality = 4:6, temperature = 8:10)
CV_WINDOW <- c(mortality = "early", temperature = "full")
COMMUNITY_CV_THRESHOLD <- 0.25
DIVERSITY_METRIC <- "shannon"

is_fluctuating_composite <- function(experiment, community_cv) {
  !is.na(community_cv) & community_cv >= COMMUNITY_CV_THRESHOLD
}

COL_FLUCTUATION <- "#F4A460"
COL_STABLE      <- "#9B8EC4"
LINE_COLS <- c("1" = "#4C78A8", "2" = "#72B7B2", "3" = "#E07B67")

FONT_FAMILY <- "Arial"
FONT_AX     <- 12
FONT_LEGEND <- 11
BORDER_SIZE <- 0.25
TICK_SIZE   <- 0.20
TICK_LEN    <- -1.0

W_SMALL <- 92
H_SMALL <- 78
W_SCAT  <- 82
H_SCAT  <- 72
W_TABLE <- 250
H_TABLE <- 180

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
      axis.text             = element_text(size = FONT_AX, color = "black"),
      axis.title            = element_text(size = FONT_AX, color = "black"),
      legend.text           = element_text(size = FONT_LEGEND, color = "black"),
      legend.title          = element_text(size = FONT_LEGEND, color = "black"),
      legend.position       = "none",
      plot.background       = element_rect(fill = "white", color = NA),
      panel.background      = element_rect(fill = "white", color = NA)
    )
}

plot_replicate_color_key <- function() {
  key_df <- tibble(
    replica = factor(names(LINE_COLS), levels = names(LINE_COLS)),
    x = seq_along(LINE_COLS),
    y = 1
  )

  p <- ggplot(key_df, aes(x = x, y = y, color = replica)) +
    geom_segment(aes(xend = x + 0.55, yend = y), linewidth = 1.2,
                 lineend = "round") +
    geom_point(size = 2.2) +
    geom_text(aes(x = x + 0.75, label = paste0("R", replica)),
              hjust = 0, size = FONT_LEGEND / 2.845,
              family = FONT_FAMILY, color = "black") +
    scale_color_manual(values = LINE_COLS, drop = FALSE) +
    coord_cartesian(xlim = c(0.8, length(LINE_COLS) + 1.4), ylim = c(0.8, 1.2),
                    clip = "off") +
    theme_void(base_family = FONT_FAMILY) +
    theme(legend.position = "none",
          plot.background = element_rect(fill = "white", color = NA))

  ggsave(file.path(FIG_DIR, "replicate_color_key.pdf"), p,
         width = 70, height = 18, units = "mm", device = cairo_pdf)
}

parse_unified_workbook <- function(path, experiment, value_col, max_day,
                                   background = 0) {
  empty <- tibble(
    experiment = character(),
    condition  = character(),
    community  = integer(),
    replica    = integer(),
    day        = integer()
  )
  empty[[value_col]] <- numeric()

  if (!file.exists(path)) {
    warning(sprintf("Workbook not found: %s", path))
    return(empty)
  }

  rows <- list()
  for (sheet in paste0("W", 1:5)) {
    raw <- read_excel(path, sheet = sheet, col_names = FALSE)
    current_comm <- NA_integer_

    for (i in seq_len(nrow(raw))) {
      c0 <- raw[[1]][i]
      c0_chr <- ifelse(is.na(c0), "", as.character(c0))

      if (startsWith(c0_chr, "Community") && c0_chr != "Community/Replicate") {
        current_comm <- as.integer(str_extract(c0_chr, "\\d+"))
      }

      rep_label <- if (ncol(raw) >= 13) as.character(raw[[13]][i]) else NA_character_
      if (is.na(rep_label) || !startsWith(rep_label, "R") || is.na(current_comm)) next

      replica <- as.integer(str_extract(rep_label, "\\d+"))
      for (day in seq_len(max_day)) {
        col_i <- day + 2
        if (col_i > ncol(raw)) next
        val <- suppressWarnings(as.numeric(raw[[col_i]][i]))
        if (is.na(val)) next
        if (value_col == "OD") {
          val <- max(0, val - background)
        }
        rows[[length(rows) + 1]] <- tibble(
          experiment = experiment,
          condition  = sheet,
          community  = current_comm,
          replica    = replica,
          day        = day,
          value      = val
        )
      }
    }
  }

  if (length(rows) == 0) {
    warning(sprintf("No %s rows parsed from: %s", value_col, path))
    return(empty)
  }

  bind_rows(rows) %>% rename(!!value_col := value)
}

load_all_trajectories <- function(files, value_col, background = NULL) {
  out <- bind_rows(lapply(names(files), function(exp) {
    bg <- if (is.null(background)) 0 else background[[exp]]
    parse_unified_workbook(files[[exp]], exp, value_col, MAX_DAY[[exp]], bg)
  }))

  if (nrow(out) == 0) {
    stop(
      sprintf(
        "No %s data loaded. Checked files:\n  %s",
        value_col,
        paste(unlist(files), collapse = "\n  ")
      ),
      call. = FALSE
    )
  }

  out
}

plot_series_set <- function(df, value_col, y_label, subdir) {
  out_root <- file.path(FIG_DIR, subdir)
  dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

  for (exp in unique(df$experiment)) {
    for (cond in paste0("W", 1:5)) {
      df_cond <- df %>% filter(experiment == exp, condition == cond)
      if (nrow(df_cond) == 0) next

      cond_dir <- file.path(out_root, exp, cond)
      dir.create(cond_dir, showWarnings = FALSE, recursive = TRUE)

      for (comm in sort(unique(df_cond$community))) {
        d <- df_cond %>% filter(community == comm)
        if (nrow(d) == 0) next

        y_scale <- if (value_col == "OD") {
          scale_y_continuous(
            limits = OD_LIMITS[[exp]],
            breaks = pretty_breaks(5),
            expand = expansion(mult = c(0, 0))
          )
        } else if (value_col == "pH") {
          scale_y_continuous(
            limits = PH_LIMITS,
            breaks = pretty_breaks(5),
            expand = expansion(mult = c(0, 0))
          )
        } else {
          scale_y_continuous(
            breaks = pretty_breaks(5),
            expand = expansion(mult = c(0.02, 0.08))
          )
        }

        p <- ggplot(d, aes(x = day, y = .data[[value_col]],
                           color = factor(replica), group = replica)) +
          geom_line(linewidth = 0.68, alpha = 0.96) +
          geom_point(size = 1.55, alpha = 0.96) +
          scale_color_manual(values = LINE_COLS, drop = FALSE) +
          scale_x_continuous(
            limits = c(1, MAX_DAY[[exp]]),
            breaks = pretty_breaks(5),
            expand = expansion(mult = c(0, 0))
          ) +
          y_scale +
          labs(x = "Day", y = y_label) +
          theme_pub_axis()

        ggsave(
          file.path(cond_dir, sprintf("C%02d_%s.pdf", comm, value_col)),
          p, width = W_SMALL, height = H_SMALL, units = "mm",
          device = cairo_pdf
        )
      }
    }
  }
}

plot_od_ph_scatter <- function(od_df, ph_df) {
  joined <- od_df %>%
    inner_join(ph_df, by = c("experiment", "condition", "community",
                             "replica", "day"))

  out_all <- file.path(FIG_DIR, "od_ph_scatter_all_days")
  out_last <- file.path(FIG_DIR, "od_ph_scatter_final_day")
  dir.create(out_all, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_last, showWarnings = FALSE, recursive = TRUE)

  for (exp in unique(joined$experiment)) {
    for (cond in paste0("W", 1:5)) {
      d <- joined %>% filter(experiment == exp, condition == cond)
      if (nrow(d) == 0) next

      p_all <- ggplot(d, aes(x = pH, y = OD, color = day)) +
        geom_point(size = 1.4, alpha = 0.72) +
        scale_color_viridis_c(option = "D", end = 0.9) +
        scale_x_continuous(
          limits = PH_LIMITS,
          breaks = pretty_breaks(5),
          expand = expansion(mult = c(0, 0))
        ) +
        scale_y_continuous(
          limits = OD_LIMITS[[exp]],
          breaks = pretty_breaks(5),
          expand = expansion(mult = c(0, 0))
        ) +
        labs(x = "pH", y = "OD") +
        theme_pub_axis() +
        theme(legend.position = "right")

      ggsave(file.path(out_all, sprintf("%s_%s_OD_pH_all_days.pdf", exp, cond)),
             p_all, width = W_SCAT, height = H_SCAT, units = "mm",
             device = cairo_pdf)

      final_day <- MAX_DAY[[exp]]
      d_last <- d %>% filter(day == final_day)
      if (nrow(d_last) == 0) next

      p_last <- ggplot(d_last, aes(x = pH, y = OD)) +
        geom_point(size = 1.6, alpha = 0.82, color = "#2E6DA4") +
        scale_x_continuous(
          limits = PH_LIMITS,
          breaks = pretty_breaks(5),
          expand = expansion(mult = c(0, 0))
        ) +
        scale_y_continuous(
          limits = OD_LIMITS[[exp]],
          breaks = pretty_breaks(5),
          expand = expansion(mult = c(0, 0))
        ) +
        labs(x = "pH", y = "OD") +
        theme_pub_axis()

      ggsave(file.path(out_last, sprintf("%s_%s_OD_pH_final_day.pdf", exp, cond)),
             p_last, width = W_SCAT, height = H_SCAT, units = "mm",
             device = cairo_pdf)
    }
  }
}

compute_od_divergence <- function(od_df) {
  bind_rows(lapply(unique(od_df$experiment), function(exp) {
    days <- LAST3_DAYS[[exp]]
    od_df %>%
      filter(experiment == exp, day %in% days) %>%
      group_by(experiment, condition, community) %>%
      group_modify(function(d, key) {
        reps <- sort(unique(d$replica))
        pair_vals <- c()
        if (length(reps) >= 2) {
          for (pair in combn(reps, 2, simplify = FALSE)) {
            a <- d %>% filter(replica == pair[1]) %>% select(day, OD)
            b <- d %>% filter(replica == pair[2]) %>% select(day, OD)
            ab <- inner_join(a, b, by = "day", suffix = c("_a", "_b"))
            if (nrow(ab) >= 2) {
              pair_vals <- c(pair_vals, sqrt(sum((ab$OD_a - ab$OD_b)^2)))
            }
          }
        }
        tibble(
          od_divergence = ifelse(length(pair_vals) > 0,
                                 mean(pair_vals, na.rm = TRUE), NA_real_),
          n_replica_pairs = length(pair_vals)
        )
      }) %>%
      ungroup()
  }))
}

load_cv_r1 <- function(exp) {
  path <- file.path(PROC_DIR, "fluctuations", CV_WINDOW[[exp]],
                    "community_level.csv")
  empty <- tibble(
    experiment = character(),
    condition = character(),
    community = integer(),
    community_cv = numeric(),
    sum_abs_std = numeric(),
    temporal_bc_mean = numeric(),
    fluctuating = logical()
  )
  if (!file.exists(path)) {
    warning(sprintf("CV file not found: %s", path))
    return(empty)
  }
  read_csv(path, show_col_types = FALSE) %>%
    filter(experiment == exp, replica == 1) %>%
    select(experiment, condition, community,
           community_cv, sum_abs_std, temporal_bc_mean) %>%
    mutate(
      fluctuating = is_fluctuating_composite(experiment, community_cv)
    )
}

load_diversity_r1 <- function(exp) {
  path <- file.path(PROC_DIR, "diversity", CV_WINDOW[[exp]],
                    "alpha_diversity.csv")
  empty <- tibble(
    experiment = character(),
    condition = character(),
    community = integer(),
    diversity_r1 = numeric()
  )
  if (!file.exists(path)) {
    warning(sprintf("Alpha diversity file not found: %s", path))
    return(empty)
  }
  read_csv(path, show_col_types = FALSE) %>%
    filter(experiment == exp, replica == 1) %>%
    select(experiment, condition, community,
           diversity_r1 = all_of(DIVERSITY_METRIC))
}

plot_divergence_panels <- function(od_df) {
  div_df <- compute_od_divergence(od_df)
  out_dir <- file.path(FIG_DIR, "replicate_od_divergence")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  for (exp in c("mortality", "temperature")) {
    cv_df <- load_cv_r1(exp)
    if (nrow(cv_df) == 0) {
      warning(sprintf("Skip %s divergence plots: no R1 community_cv data.", exp))
      next
    }
    d <- div_df %>%
      filter(experiment == exp) %>%
      inner_join(cv_df, by = c("experiment", "condition", "community")) %>%
      filter(is.finite(od_divergence), is.finite(community_cv))

    if (nrow(d) > 0) {
      p_cv <- ggplot(d, aes(x = community_cv, y = od_divergence,
                            color = fluctuating)) +
        geom_point(size = 1.8, alpha = 0.9) +
        scale_color_manual(values = c("TRUE" = COL_FLUCTUATION,
                                      "FALSE" = COL_STABLE)) +
        scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
        scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
        labs(x = "Community CV (R1)", y = "Replicate OD divergence") +
        theme_pub_axis()

      ggsave(file.path(out_dir, sprintf("%s_divergence_vs_community_cv.pdf", exp)),
             p_cv, width = W_SCAT, height = H_SCAT, units = "mm",
             device = cairo_pdf)
    }

    div_r1 <- load_diversity_r1(exp)
    if (nrow(div_r1) == 0) {
      warning(sprintf("Skip %s divergence-vs-R1-shannon plot: no R1 alpha diversity data.", exp))
      next
    }
    d2 <- d %>%
      inner_join(div_r1, by = c("experiment", "condition", "community")) %>%
      filter(is.finite(diversity_r1))

    if (nrow(d2) > 0) {
      p_div <- ggplot(d2, aes(x = diversity_r1, y = od_divergence,
                              color = fluctuating)) +
        geom_point(size = 1.8, alpha = 0.9) +
        scale_color_manual(values = c("TRUE" = COL_FLUCTUATION,
                                      "FALSE" = COL_STABLE)) +
        scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
        scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
        labs(x = "R1 Shannon diversity", y = "Replicate OD divergence") +
        theme_pub_axis()

      ggsave(file.path(out_dir, sprintf("%s_divergence_vs_r1_shannon.pdf", exp)),
             p_div, width = W_SCAT, height = H_SCAT, units = "mm",
             device = cairo_pdf)
    } else {
      warning(sprintf(
        "Skip %s divergence-vs-R1-shannon plot: joined data has zero finite rows.",
        exp
      ))
    }
  }
}

format_metric_value <- function(x, digits = 2) {
  ifelse(is.na(x), "",
         formatC(x, digits = digits, format = "f"))
}

load_alpha_last_day <- function() {
  alpha_paths <- c(
    mortality   = file.path(PROC_DIR, "diversity", "early", "alpha_diversity.csv"),
    temperature = file.path(PROC_DIR, "diversity", "full", "alpha_diversity.csv")
  )

  bind_rows(lapply(names(alpha_paths), function(exp) {
    if (!file.exists(alpha_paths[[exp]])) return(tibble())
    read_csv(alpha_paths[[exp]], show_col_types = FALSE) %>%
      filter(experiment == exp) %>%
      select(experiment, condition, community, replica,
             final_day = day,
             final_shannon = shannon,
             final_richness = richness)
  }))
}

load_last3_decomposition <- function() {
  path <- file.path(PROC_DIR, "species_decomposition", "last3",
                    "community_decomposition.csv")
  if (!file.exists(path)) return(tibble())

  read_csv(path, show_col_types = FALSE) %>%
    select(experiment, condition, community, replica,
           community_cv,
           synchrony_phi,
           turnover_events,
           turnover_event_rate,
           cumulative_shannon,
           mean_total_biomass)
}

load_main_fluctuations <- function() {
  fluc_paths <- c(
    mortality   = file.path(PROC_DIR, "fluctuations", "early",
                            "community_level.csv"),
    temperature = file.path(PROC_DIR, "fluctuations", "full",
                            "community_level.csv")
  )

  bind_rows(lapply(names(fluc_paths), function(exp) {
    if (!file.exists(fluc_paths[[exp]])) return(tibble())
    read_csv(fluc_paths[[exp]], show_col_types = FALSE) %>%
      filter(experiment == exp) %>%
      mutate(
        temporal_bc_mean = as.numeric(temporal_bc_mean),
        community_cv = as.numeric(community_cv),
        sum_abs_std = as.numeric(sum_abs_std),
        fluctuating_composite = is_fluctuating_composite(
          experiment, community_cv
        )
      ) %>%
      select(experiment, condition, community, replica,
             community_cv_main = community_cv,
             total_biomass_std,
             total_biomass_cv,
             sum_abs_std,
             temporal_bc_mean)
  }))
}

make_metric_summary_table <- function(od_df) {
  out_dir <- file.path(FIG_DIR, "metric_summary_tables")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  alpha_df <- load_alpha_last_day()
  decomp_df <- load_last3_decomposition()
  fluc_df <- load_main_fluctuations()

  final_od <- bind_rows(lapply(unique(od_df$experiment), function(exp) {
    od_df %>%
      filter(experiment == exp, day == MAX_DAY[[exp]]) %>%
      transmute(experiment, condition, community, replica, final_OD = OD)
  }))

  mean_last3_od <- bind_rows(lapply(unique(od_df$experiment), function(exp) {
    od_df %>%
      filter(experiment == exp, day %in% LAST3_DAYS[[exp]]) %>%
      group_by(experiment, condition, community, replica) %>%
      summarise(mean_last3_OD = mean(OD, na.rm = TRUE), .groups = "drop")
  }))

  summary_df <- od_df %>%
    distinct(experiment, condition, community, replica) %>%
    left_join(final_od,
              by = c("experiment", "condition", "community", "replica")) %>%
    left_join(mean_last3_od,
              by = c("experiment", "condition", "community", "replica")) %>%
    left_join(alpha_df,
              by = c("experiment", "condition", "community", "replica")) %>%
    left_join(decomp_df,
              by = c("experiment", "condition", "community", "replica")) %>%
    left_join(fluc_df,
              by = c("experiment", "condition", "community", "replica")) %>%
    mutate(
      community_cv_display = coalesce(community_cv_main, community_cv)
    ) %>%
    arrange(experiment, condition, community, replica)

  write_csv(summary_df, file.path(out_dir, "metric_summary_table_values.csv"))

  metric_order <- c(
    "Shannon", "Richness", "OD", "Sum species SD", "Temporal BC",
    "Community CV", "Total biomass SD", "Total biomass CV", "Phi",
    "Turnover rate", "Cumulative Shannon"
  )

  table_long <- summary_df %>%
    transmute(
      experiment, condition, community,
      replica = factor(replica, levels = 1:3),
      Shannon = final_shannon,
      Richness = final_richness,
      OD = final_OD,
      `Sum species SD` = sum_abs_std,
      `Temporal BC` = temporal_bc_mean,
      `Community CV` = community_cv_display,
      `Total biomass SD` = total_biomass_std,
      `Total biomass CV` = total_biomass_cv,
      Phi = synchrony_phi,
      `Turnover rate` = turnover_event_rate,
      `Cumulative Shannon` = cumulative_shannon
    ) %>%
    pivot_longer(cols = all_of(metric_order),
                 names_to = "metric", values_to = "value",
                 values_transform = list(value = as.character)) %>%
    mutate(
      metric = factor(metric, levels = metric_order),
      row_label = factor(
        paste0(metric, "  R", replica),
        levels = rev(as.vector(outer(metric_order, paste0("R", 1:3),
                                     paste, sep = "  ")))
      ),
      value_numeric = suppressWarnings(as.numeric(value)),
      value_label = case_when(
        metric == "Richness" ~ ifelse(
          is.na(value_numeric), "", as.character(round(value_numeric))
        ),
        TRUE ~ format_metric_value(value_numeric, 2)
      )
    )

  for (exp in unique(table_long$experiment)) {
    for (cond in paste0("W", 1:5)) {
      d <- table_long %>% filter(experiment == exp, condition == cond)
      if (nrow(d) == 0) next

      p <- ggplot(d, aes(x = factor(community), y = row_label)) +
        geom_tile(fill = "white", color = "grey82", linewidth = 0.25) +
        geom_text(aes(label = value_label),
                  size = 2.15, family = FONT_FAMILY, color = "black",
                  na.rm = TRUE) +
        scale_x_discrete(position = "top", drop = FALSE) +
        labs(x = "Community", y = NULL) +
        theme_minimal(base_size = 8, base_family = FONT_FAMILY) +
        theme(
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 8, color = "black"),
          axis.text.y = element_text(size = 7, color = "black"),
          axis.title.x = element_text(size = 9, color = "black"),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(4, 4, 4, 4, "mm")
        )

      ggsave(file.path(out_dir, sprintf("%s_%s_metric_summary_table.pdf",
                                        exp, cond)),
             p, width = W_TABLE, height = H_TABLE, units = "mm",
             device = cairo_pdf)
    }
  }

  summary_df
}

spearman_label <- function(d, x_col, y_col) {
  d <- d %>% filter(is.finite(.data[[x_col]]), is.finite(.data[[y_col]]))
  if (nrow(d) < 3) return("Spearman rho = NA\nP = NA")
  ct <- suppressWarnings(cor.test(d[[x_col]], d[[y_col]], method = "spearman",
                                  exact = FALSE))
  sprintf("Spearman rho = %.2f\nP = %.3g", unname(ct$estimate), ct$p.value)
}

plot_spearman_panel <- function(d, x_col, y_col, x_lab, y_lab, out_file,
                                color_col = NULL) {
  d <- d %>% filter(is.finite(.data[[x_col]]), is.finite(.data[[y_col]]))
  if (nrow(d) < 3) return(invisible(NULL))

  label <- spearman_label(d, x_col, y_col)
  label_x <- min(d[[x_col]], na.rm = TRUE)
  label_y <- max(d[[y_col]], na.rm = TRUE)

  p <- ggplot(d, aes(x = .data[[x_col]], y = .data[[y_col]]))
  if (!is.null(color_col)) {
    p <- p + geom_point(aes(color = .data[[color_col]]),
                        size = 1.8, alpha = 0.82)
  } else {
    p <- p + geom_point(size = 1.8, alpha = 0.82, color = "#4C78A8")
  }

  p <- p +
    annotate("text", x = label_x, y = label_y, label = label,
             hjust = 0, vjust = 1, size = 3.0, family = FONT_FAMILY) +
    labs(x = x_lab, y = y_lab) +
    theme_pub_axis() +
    theme(legend.position = ifelse(is.null(color_col), "none", "right"))

  if (!is.null(color_col)) {
    p <- p + scale_color_brewer(palette = "Set2", name = "Condition")
  }

  ggsave(out_file, p, width = W_SCAT, height = H_SCAT, units = "mm",
         device = cairo_pdf)
}

plot_spearman_analyses <- function(summary_df) {
  out_dir <- file.path(FIG_DIR, "spearman_analysis")
  by_cond_dir <- file.path(out_dir, "by_condition")
  by_exp_dir <- file.path(out_dir, "by_experiment")
  dir.create(by_cond_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(by_exp_dir, showWarnings = FALSE, recursive = TRUE)

  corr_input <- summary_df %>%
    mutate(condition = factor(condition, levels = paste0("W", 1:5)))

  corr_rows <- list()
  for (exp in unique(corr_input$experiment)) {
    for (cond in paste0("W", 1:5)) {
      d <- corr_input %>% filter(experiment == exp, condition == cond)
      if (nrow(d) == 0) next

      plot_spearman_panel(
        d, "final_shannon", "final_OD",
        "Final-day Shannon diversity", "Final-day OD",
        file.path(by_cond_dir,
                  sprintf("%s_%s_final_shannon_vs_final_OD.pdf", exp, cond))
      )
      plot_spearman_panel(
        d, "community_cv", "mean_last3_OD",
        "Community CV (last 3 days)", "Mean OD (last 3 days)",
        file.path(by_cond_dir,
                  sprintf("%s_%s_community_CV_vs_mean_last3_OD.pdf", exp, cond))
      )

      for (pair in list(
        c("final_shannon", "final_OD", "final_shannon_vs_final_OD"),
        c("community_cv", "mean_last3_OD", "community_CV_vs_mean_last3_OD")
      )) {
        dd <- d %>% filter(is.finite(.data[[pair[1]]]),
                           is.finite(.data[[pair[2]]]))
        if (nrow(dd) >= 3) {
          ct <- suppressWarnings(cor.test(dd[[pair[1]]], dd[[pair[2]]],
                                          method = "spearman", exact = FALSE))
          corr_rows[[length(corr_rows) + 1]] <- tibble(
            experiment = exp, condition = cond, analysis = pair[3],
            n = nrow(dd), rho = unname(ct$estimate), p_value = ct$p.value
          )
        }
      }
    }

    d_exp <- corr_input %>% filter(experiment == exp)
    plot_spearman_panel(
      d_exp, "final_shannon", "final_OD",
      "Final-day Shannon diversity", "Final-day OD",
      file.path(by_exp_dir,
                sprintf("%s_final_shannon_vs_final_OD_all_conditions.pdf", exp)),
      color_col = "condition"
    )
    plot_spearman_panel(
      d_exp, "community_cv", "mean_last3_OD",
      "Community CV (last 3 days)", "Mean OD (last 3 days)",
      file.path(by_exp_dir,
                sprintf("%s_community_CV_vs_mean_last3_OD_all_conditions.pdf", exp)),
      color_col = "condition"
    )

    for (pair in list(
      c("final_shannon", "final_OD", "final_shannon_vs_final_OD"),
      c("community_cv", "mean_last3_OD", "community_CV_vs_mean_last3_OD")
    )) {
      dd <- d_exp %>% filter(is.finite(.data[[pair[1]]]),
                             is.finite(.data[[pair[2]]]))
      if (nrow(dd) >= 3) {
        ct <- suppressWarnings(cor.test(dd[[pair[1]]], dd[[pair[2]]],
                                        method = "spearman", exact = FALSE))
        corr_rows[[length(corr_rows) + 1]] <- tibble(
          experiment = exp, condition = "all_conditions", analysis = pair[3],
          n = nrow(dd), rho = unname(ct$estimate), p_value = ct$p.value
        )
      }
    }
  }

  if (length(corr_rows) > 0) {
    write_csv(bind_rows(corr_rows),
              file.path(out_dir, "spearman_results.csv"))
  }
}

cat("Loading OD and pH workbooks...\n")
od_df <- load_all_trajectories(OD_FILES, "OD", OD_BACKGROUND)
ph_df <- load_all_trajectories(PH_FILES, "pH")

cat("Plotting replicate color key...\n")
plot_replicate_color_key()

cat("Plotting OD trajectories...\n")
plot_series_set(od_df, "OD", "OD", "od_series")

cat("Plotting pH trajectories...\n")
plot_series_set(ph_df, "pH", "pH", "pH_series")

cat("Plotting OD-pH scatter plots...\n")
plot_od_ph_scatter(od_df, ph_df)

cat("Plotting replicate OD divergence panels...\n")
plot_divergence_panels(od_df)

cat("Plotting metric summary tables...\n")
summary_df <- make_metric_summary_table(od_df)

cat("Plotting Spearman analyses...\n")
plot_spearman_analyses(summary_df)

cat("Done. Output directory: figures/supplement_od_ph\n")
