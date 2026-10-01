community_workspace <- function() {
  override <- Sys.getenv("COMMUNITY_WORKSPACE", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  script <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(script)) stop("Run with Rscript or set COMMUNITY_WORKSPACE.")
  dirname(normalizePath(sub("^--file=", "", script[[1]]), mustWork = TRUE))
}

# ===========================================================
# calculate_systematic_significance_tests.R
#
# Systematic significance tests for fluctuation proportions,
# fluctuation metrics, and diversity metrics.
#
# Primary biological hypotheses:
#   mortality:
#     W1 is more fluctuating than W5; fluctuation decreases from W1 to W5.
#   temperature:
#     W3 is more fluctuating than W1/W2/W4/W5.
#
# Temperature early-window tests are intentionally skipped.
#
# Outputs:
#   processed/statistical_tests/
#     fluctuation_proportion_pairwise_tests.csv
#     fluctuation_proportion_trend_tests.csv
#     fluctuation_metric_pairwise_tests.csv
#     fluctuation_metric_trend_tests.csv
#     diversity_pairwise_tests.csv
#     species_decomposition_pairwise_tests.csv
# ===========================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

BASE_DIR <- community_workspace()
PROC_DIR <- file.path(BASE_DIR, "processed")
OUT_DIR  <- file.path(PROC_DIR, "statistical_tests")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

COMMUNITY_CV_THRESHOLD <- 0.25
TOTAL_BIOMASS_CV_THRESHOLD <- 0.10

FLUCTUATION_WINDOWS <- tribble(
  ~experiment,   ~window, ~analysis_set,
  "mortality",   "early", "main",
  "mortality",   "last4", "last4",
  "temperature", "full",  "main",
  "temperature", "last4", "last4"
)

DIVERSITY_WINDOWS <- tribble(
  ~experiment,   ~window, ~analysis_set,
  "mortality",   "early", "main",
  "temperature", "full",  "main"
)

SPECIES_DECOMP_WINDOWS <- tribble(
  ~experiment,   ~window, ~analysis_set,
  "mortality",   "last3", "last3",
  "mortality",   "long",  "long",
  "temperature", "last3", "last3",
  "temperature", "long",  "long"
)

FLUCTUATION_METRICS <- c(
  "community_cv",
  "temporal_bc_mean",
  "total_biomass_cv",
  "weighted_log_sd",
  "sum_abs_std",
  "total_biomass_std"
)

BINARY_METRICS <- c(
  "composite_fluctuation",
  "community_cv_binary",
  "total_biomass_cv_binary"
)

DIVERSITY_INPUTS <- list(
  alpha = list(
    file = "alpha_diversity.csv",
    metrics = c("richness", "survival_fraction", "shannon",
                "effective_shannon", "simpson", "evenness")
  ),
  gamma = list(
    file = "gamma_diversity.csv",
    metrics = c("gamma_richness", "gamma_shannon",
                "gamma_effective_shannon", "gamma_simpson",
                "gamma_evenness")
  ),
  mean_daily = list(
    file = "mean_daily_diversity.csv",
    metrics = c("mean_daily_richness", "mean_daily_shannon",
                "mean_daily_effective_shannon", "mean_daily_simpson",
                "mean_daily_evenness", "mean_daily_survival_fraction")
  )
)

SPECIES_DECOMP_METRICS <- c(
  "synchrony_phi",
  "synchrony_phi_robust",
  "weighted_mean_pairwise_correlation",
  "relative_covariance_contribution"
)

# For one-sided fluctuation tests:
# alternative = "greater" means condition_a > condition_b.
# alternative = "two.sided" is retained for neutral pairwise checks.
FLUCTUATION_CONTRASTS <- list(
  mortality = tribble(
    ~condition_a, ~condition_b, ~alternative, ~hypothesis,
    "W1", "W5", "greater", "W1_more_fluctuating_than_W5",
    "W1", "W2", "greater", "decreasing_with_dilution_adjacent",
    "W2", "W3", "greater", "decreasing_with_dilution_adjacent",
    "W3", "W4", "greater", "decreasing_with_dilution_adjacent",
    "W4", "W5", "greater", "decreasing_with_dilution_adjacent",
    "W1", "W3", "greater", "decreasing_with_dilution_broad",
    "W3", "W5", "greater", "decreasing_with_dilution_broad"
  ),
  temperature = tribble(
    ~condition_a, ~condition_b, ~alternative, ~hypothesis,
    "W3", "W1", "greater", "W3_peak",
    "W3", "W2", "greater", "W3_peak",
    "W3", "W4", "greater", "W3_peak",
    "W3", "W5", "greater", "W3_peak",
    "W1", "W2", "two.sided", "W1_vs_W2_pairwise",
    "W4", "W5", "two.sided", "W4_vs_W5_pairwise"
  )
)

DIVERSITY_CONTRASTS <- list(
  mortality = tribble(
    ~condition_a, ~condition_b, ~alternative, ~hypothesis,
    "W1", "W5", "two.sided", "W1_vs_W5",
    "W1", "W2", "two.sided", "W1_vs_W2",
    "W2", "W3", "two.sided", "W2_vs_W3",
    "W3", "W4", "two.sided", "W3_vs_W4",
    "W4", "W5", "two.sided", "W4_vs_W5"
  ),
  temperature = tribble(
    ~condition_a, ~condition_b, ~alternative, ~hypothesis,
    "W3", "W1", "two.sided", "W3_vs_W1",
    "W3", "W2", "two.sided", "W3_vs_W2",
    "W3", "W4", "two.sided", "W3_vs_W4",
    "W3", "W5", "two.sided", "W3_vs_W5",
    "W1", "W2", "two.sided", "W1_vs_W2",
    "W4", "W5", "two.sided", "W4_vs_W5"
  )
)

SPECIES_DECOMP_CONTRASTS <- list(
  mortality = tribble(
    ~condition_a, ~condition_b, ~alternative, ~hypothesis,
    "W5", "W1", "less", "complexity_driven_asynchrony_W5_lt_W1",
    "W5", "W3", "less", "complexity_driven_asynchrony_W5_lt_W3",
    "W3", "W1", "two.sided", "intermediate_vs_W1"
  ),
  temperature = tribble(
    ~condition_a, ~condition_b, ~alternative, ~hypothesis,
    "W3", "W1", "less", "complexity_driven_asynchrony_W3_trough",
    "W3", "W2", "less", "complexity_driven_asynchrony_W3_trough",
    "W3", "W4", "less", "complexity_driven_asynchrony_W3_trough",
    "W3", "W5", "less", "complexity_driven_asynchrony_W3_trough",
    "W1", "W2", "two.sided", "low_temperature_pairwise",
    "W4", "W5", "two.sided", "high_temperature_pairwise"
  )
)

condition_pos <- function(condition) {
  as.integer(sub("W", "", condition))
}

add_binary_fluctuation_metrics <- function(df) {
  df %>%
    mutate(
      community_cv_threshold = COMMUNITY_CV_THRESHOLD,
      community_cv_binary = as.numeric(
        !is.na(community_cv) &
          community_cv >= community_cv_threshold
      ),
      total_biomass_cv_binary = as.numeric(
        !is.na(total_biomass_cv) &
          total_biomass_cv >= TOTAL_BIOMASS_CV_THRESHOLD
      ),
      composite_fluctuation = as.numeric(
        community_cv_binary == 1
      ),
      x_pos = condition_pos(condition)
    )
}

read_fluctuation_table <- function(window) {
  path <- file.path(PROC_DIR, "fluctuations", window, "community_level.csv")
  if (!file.exists(path)) return(tibble())
  read_csv(path, show_col_types = FALSE) %>%
    mutate(window = window) %>%
    add_binary_fluctuation_metrics()
}

read_diversity_table <- function(window, input_name, input_cfg) {
  path <- file.path(PROC_DIR, "diversity", window, input_cfg$file)
  if (!file.exists(path)) return(tibble())

  read_csv(path, show_col_types = FALSE) %>%
    mutate(
      window = window,
      input = input_name,
      x_pos = condition_pos(condition)
    )
}

read_species_decomp_table <- function(window) {
  path <- file.path(PROC_DIR, "species_decomposition",
                    window, "community_decomposition.csv")
  if (!file.exists(path)) return(tibble())

  read_csv(path, show_col_types = FALSE) %>%
    mutate(
      window = window,
      x_pos = condition_pos(condition)
    )
}

one_sided_from_two_sided <- function(estimate, p_two, direction) {
  if (is.na(estimate) || is.na(p_two)) return(NA_real_)
  if (direction == "less") {
    if (estimate < 0) p_two / 2 else 1 - p_two / 2
  } else if (direction == "greater") {
    if (estimate > 0) p_two / 2 else 1 - p_two / 2
  } else {
    p_two
  }
}

safe_fisher <- function(mat, alternative) {
  tryCatch(
    fisher.test(mat, alternative = alternative)$p.value,
    error = function(e) NA_real_
  )
}

safe_prop <- function(success, total, alternative) {
  tryCatch(
    suppressWarnings(
      prop.test(success, total, alternative = alternative,
                correct = FALSE)$p.value
    ),
    error = function(e) NA_real_
  )
}

safe_t_test <- function(x, y, alternative) {
  if (length(x) < 2 || length(y) < 2) return(list(p = NA_real_))
  tryCatch(
    list(p = t.test(x, y, alternative = alternative)$p.value),
    error = function(e) list(p = NA_real_)
  )
}

safe_wilcox <- function(x, y, alternative) {
  if (length(x) < 1 || length(y) < 1) return(list(p = NA_real_))
  tryCatch(
    list(p = suppressWarnings(
      wilcox.test(x, y, alternative = alternative, exact = FALSE)$p.value
    )),
    error = function(e) list(p = NA_real_)
  )
}

pairwise_binary_tests <- function(df, metric, experiment, window,
                                  analysis_set) {
  contrasts <- FLUCTUATION_CONTRASTS[[experiment]]
  if (is.null(contrasts) || nrow(df) == 0) return(tibble())

  counts <- df %>%
    filter(!is.na(.data[[metric]])) %>%
    group_by(condition) %>%
    summarise(
      n_osc = sum(.data[[metric]] >= 0.5, na.rm = TRUE),
      n_total = n(),
      proportion = n_osc / n_total,
      .groups = "drop"
    )

  bind_rows(lapply(seq_len(nrow(contrasts)), function(i) {
    con <- contrasts[i, ]
    a <- counts %>% filter(condition == con$condition_a)
    b <- counts %>% filter(condition == con$condition_b)
    if (nrow(a) == 0 || nrow(b) == 0) return(tibble())

    mat <- matrix(
      c(a$n_osc, a$n_total - a$n_osc,
        b$n_osc, b$n_total - b$n_osc),
      nrow = 2, byrow = TRUE
    )

    tibble(
      analysis_type = "fluctuation_proportion",
      metric = metric,
      experiment = experiment,
      window = window,
      analysis_set = analysis_set,
      condition_a = con$condition_a,
      condition_b = con$condition_b,
      alternative = con$alternative,
      hypothesis = con$hypothesis,
      n_osc_a = a$n_osc,
      n_total_a = a$n_total,
      proportion_a = a$proportion,
      n_osc_b = b$n_osc,
      n_total_b = b$n_total,
      proportion_b = b$proportion,
      diff_a_minus_b = a$proportion - b$proportion,
      fisher_p = safe_fisher(mat, con$alternative),
      prop_test_p = safe_prop(
        c(a$n_osc, b$n_osc),
        c(a$n_total, b$n_total),
        con$alternative
      )
    )
  }))
}

pairwise_continuous_tests <- function(df, metric, experiment, window,
                                      analysis_set,
                                      contrasts,
                                      analysis_type) {
  if (is.null(contrasts) || nrow(df) == 0 || !metric %in% names(df)) {
    return(tibble())
  }

  bind_rows(lapply(seq_len(nrow(contrasts)), function(i) {
    con <- contrasts[i, ]
    x <- df %>%
      filter(condition == con$condition_a) %>%
      pull(all_of(metric)) %>%
      as.numeric()
    y <- df %>%
      filter(condition == con$condition_b) %>%
      pull(all_of(metric)) %>%
      as.numeric()
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
    if (length(x) == 0 || length(y) == 0) return(tibble())

    t_res <- safe_t_test(x, y, con$alternative)
    w_res <- safe_wilcox(x, y, con$alternative)

    tibble(
      analysis_type = analysis_type,
      metric = metric,
      experiment = experiment,
      window = window,
      analysis_set = analysis_set,
      condition_a = con$condition_a,
      condition_b = con$condition_b,
      alternative = con$alternative,
      hypothesis = con$hypothesis,
      n_a = length(x),
      n_b = length(y),
      mean_a = mean(x),
      mean_b = mean(y),
      median_a = median(x),
      median_b = median(y),
      diff_mean_a_minus_b = mean(x) - mean(y),
      t_test_p = t_res$p,
      wilcox_p = w_res$p
    )
  }))
}

trend_binary_tests <- function(df, metric, experiment, window, analysis_set) {
  if (experiment != "mortality") return(tibble())
  d <- df %>%
    filter(!is.na(.data[[metric]]), is.finite(x_pos))
  if (nrow(d) < 3 || length(unique(d$x_pos)) < 3) return(tibble())

  fit <- tryCatch(
    glm(as.formula(paste(metric, "~ x_pos")), data = d,
        family = binomial()),
    error = function(e) NULL
  )
  if (is.null(fit)) return(tibble())
  co <- summary(fit)$coefficients
  if (!"x_pos" %in% rownames(co)) return(tibble())

  estimate <- co["x_pos", "Estimate"]
  p_two <- co["x_pos", "Pr(>|z|)"]

  tibble(
    analysis_type = "fluctuation_proportion_trend",
    metric = metric,
    experiment = experiment,
    window = window,
    analysis_set = analysis_set,
    model = "binomial_glm",
    hypothesis = "decreases_from_W1_to_W5",
    slope = estimate,
    p_two_sided = p_two,
    p_one_sided_decrease = one_sided_from_two_sided(
      estimate, p_two, "less"
    ),
    n = nrow(d)
  )
}

trend_continuous_tests <- function(df, metric, experiment, window,
                                   analysis_set, analysis_type) {
  if (experiment != "mortality" || !metric %in% names(df)) return(tibble())
  d <- df %>%
    filter(is.finite(.data[[metric]]), is.finite(x_pos))
  if (nrow(d) < 3 || length(unique(d$x_pos)) < 3) return(tibble())

  fit <- tryCatch(
    lm(as.formula(paste(metric, "~ x_pos")), data = d),
    error = function(e) NULL
  )
  if (is.null(fit)) return(tibble())
  co <- summary(fit)$coefficients
  if (!"x_pos" %in% rownames(co)) return(tibble())

  estimate <- co["x_pos", "Estimate"]
  p_two <- co["x_pos", "Pr(>|t|)"]

  tibble(
    analysis_type = analysis_type,
    metric = metric,
    experiment = experiment,
    window = window,
    analysis_set = analysis_set,
    model = "linear_model",
    hypothesis = "decreases_from_W1_to_W5",
    slope = estimate,
    p_two_sided = p_two,
    p_one_sided_decrease = one_sided_from_two_sided(
      estimate, p_two, "less"
    ),
    n = nrow(d)
  )
}

# ===========================================================
# Fluctuation tests
# ===========================================================

fluct_rows <- bind_rows(lapply(unique(FLUCTUATION_WINDOWS$window),
                               read_fluctuation_table))

fluct_prop_pairwise <- list()
fluct_prop_trend <- list()
fluct_metric_pairwise <- list()
fluct_metric_trend <- list()

for (i in seq_len(nrow(FLUCTUATION_WINDOWS))) {
  cfg <- FLUCTUATION_WINDOWS[i, ]
  d <- fluct_rows %>%
    filter(experiment == cfg$experiment, window == cfg$window)
  if (nrow(d) == 0) next

  for (metric in BINARY_METRICS) {
    fluct_prop_pairwise[[length(fluct_prop_pairwise) + 1]] <-
      pairwise_binary_tests(
        d, metric, cfg$experiment, cfg$window, cfg$analysis_set
      )
    fluct_prop_trend[[length(fluct_prop_trend) + 1]] <-
      trend_binary_tests(
        d, metric, cfg$experiment, cfg$window, cfg$analysis_set
      )
  }

  for (metric in FLUCTUATION_METRICS) {
    fluct_metric_pairwise[[length(fluct_metric_pairwise) + 1]] <-
      pairwise_continuous_tests(
        d, metric, cfg$experiment, cfg$window, cfg$analysis_set,
        FLUCTUATION_CONTRASTS[[cfg$experiment]],
        "fluctuation_metric"
      )
    fluct_metric_trend[[length(fluct_metric_trend) + 1]] <-
      trend_continuous_tests(
        d, metric, cfg$experiment, cfg$window, cfg$analysis_set,
        "fluctuation_metric_trend"
      )
  }
}

fluct_prop_pairwise_df <- bind_rows(fluct_prop_pairwise) %>%
  group_by(metric, experiment, window) %>%
  mutate(
    fisher_p_adj_BH = p.adjust(fisher_p, method = "BH"),
    prop_test_p_adj_BH = p.adjust(prop_test_p, method = "BH")
  ) %>%
  ungroup()

fluct_prop_trend_df <- bind_rows(fluct_prop_trend) %>%
  group_by(metric, experiment, window) %>%
  mutate(p_one_sided_decrease_adj_BH = p.adjust(
    p_one_sided_decrease, method = "BH"
  )) %>%
  ungroup()

fluct_metric_pairwise_df <- bind_rows(fluct_metric_pairwise) %>%
  group_by(metric, experiment, window) %>%
  mutate(
    t_test_p_adj_BH = p.adjust(t_test_p, method = "BH"),
    wilcox_p_adj_BH = p.adjust(wilcox_p, method = "BH")
  ) %>%
  ungroup()

fluct_metric_trend_df <- bind_rows(fluct_metric_trend) %>%
  group_by(metric, experiment, window) %>%
  mutate(p_one_sided_decrease_adj_BH = p.adjust(
    p_one_sided_decrease, method = "BH"
  )) %>%
  ungroup()

write_csv(
  fluct_prop_pairwise_df,
  file.path(OUT_DIR, "fluctuation_proportion_pairwise_tests.csv")
)
write_csv(
  fluct_prop_trend_df,
  file.path(OUT_DIR, "fluctuation_proportion_trend_tests.csv")
)
write_csv(
  fluct_metric_pairwise_df,
  file.path(OUT_DIR, "fluctuation_metric_pairwise_tests.csv")
)
write_csv(
  fluct_metric_trend_df,
  file.path(OUT_DIR, "fluctuation_metric_trend_tests.csv")
)

# ===========================================================
# Diversity tests
# ===========================================================

diversity_pairwise <- list()

for (i in seq_len(nrow(DIVERSITY_WINDOWS))) {
  cfg <- DIVERSITY_WINDOWS[i, ]
  contrasts <- DIVERSITY_CONTRASTS[[cfg$experiment]]

  for (input_name in names(DIVERSITY_INPUTS)) {
    input_cfg <- DIVERSITY_INPUTS[[input_name]]
    d <- read_diversity_table(cfg$window, input_name, input_cfg) %>%
      filter(experiment == cfg$experiment)
    if (nrow(d) == 0) next

    for (metric in input_cfg$metrics) {
      if (!metric %in% names(d)) next
      diversity_pairwise[[length(diversity_pairwise) + 1]] <-
        pairwise_continuous_tests(
          d, metric, cfg$experiment, cfg$window, cfg$analysis_set,
          contrasts,
          paste0("diversity_", input_name)
        )
    }
  }
}

diversity_pairwise_df <- bind_rows(diversity_pairwise) %>%
  group_by(analysis_type, metric, experiment, window) %>%
  mutate(
    t_test_p_adj_BH = p.adjust(t_test_p, method = "BH"),
    wilcox_p_adj_BH = p.adjust(wilcox_p, method = "BH")
  ) %>%
  ungroup()

write_csv(
  diversity_pairwise_df,
  file.path(OUT_DIR, "diversity_pairwise_tests.csv")
)

# ===========================================================
# Species decomposition / synchrony mechanism tests
# ===========================================================

species_rows <- bind_rows(lapply(unique(SPECIES_DECOMP_WINDOWS$window),
                                 read_species_decomp_table))
species_pairwise <- list()

for (i in seq_len(nrow(SPECIES_DECOMP_WINDOWS))) {
  cfg <- SPECIES_DECOMP_WINDOWS[i, ]
  d <- species_rows %>%
    filter(experiment == cfg$experiment, window == cfg$window)
  if (nrow(d) == 0) next

  for (metric in SPECIES_DECOMP_METRICS) {
    if (!metric %in% names(d)) next
    species_pairwise[[length(species_pairwise) + 1]] <-
      pairwise_continuous_tests(
        d, metric, cfg$experiment, cfg$window, cfg$analysis_set,
        SPECIES_DECOMP_CONTRASTS[[cfg$experiment]],
        "species_decomposition"
      )
  }
}

species_pairwise_df <- bind_rows(species_pairwise) %>%
  group_by(metric, experiment, window) %>%
  mutate(
    t_test_p_adj_BH = p.adjust(t_test_p, method = "BH"),
    wilcox_p_adj_BH = p.adjust(wilcox_p, method = "BH")
  ) %>%
  ungroup()

write_csv(
  species_pairwise_df,
  file.path(OUT_DIR, "species_decomposition_pairwise_tests.csv")
)

cat("Done.\n")
cat(sprintf("Output directory: %s\n", OUT_DIR))
cat("Primary tests:\n")
cat("  mortality: one-sided W1 > W5 and decreasing W1->W5 trend.\n")
cat("  temperature: one-sided W3 > W1/W2/W4/W5 for fluctuation metrics.\n")
cat("  diversity: same condition pairs, two-sided by default.\n")
cat("  species decomposition: one-sided tests for complexity-driven asynchrony/low-phi hypotheses.\n")
