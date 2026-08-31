# ==============================================================================
# sos_score.R
#
# Comprehensive scoring of manual vs PELT SOS, matching the project's existing
# EOS/SOF validation approach.
#
# Reports:
# - n, bias, median |error|, mean |error|, % within 8d/16d, Lin's CCC, Bland-Altman LOA
# - Three nested sets: (a) all 44, (b) excluding structural/snow, (c) also excluding bimodal
# - Breakdown by vegetation type
# - One-way ANOVA of SOS error on vegetation type
# - Comparison to established EOS/SOF agreement metrics
#
# Sign convention: diff = manual - PELT (same as project-wide)
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

nk <- function(x) gsub("[-/]", "_", trimws(x))

BASE_DIR <- "C:/Users/jmusinsky/Documents/R_Scripts/MODIS"
OUT_DIR  <- file.path(BASE_DIR, "pelt_out")

# ------------------------------------------------------------ Lin's CCC function
# Matches the formula from sos_mech.R and the main PELT script

ccc_fun <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 2) return(NA_real_)
  2 * cov(x, y) * (n - 1) / n /
    (var(x) * (n - 1) / n + var(y) * (n - 1) / n + (mean(x) - mean(y))^2)
}

# --------------------------------------------------- Bland-Altman limits of agreement

bland_altman_loa <- function(manual, pelt) {
  diff <- manual - pelt
  ok <- is.finite(diff)
  diff <- diff[ok]
  if (length(diff) < 2) return(c(lower = NA_real_, upper = NA_real_))
  bias <- mean(diff)
  sd_diff <- sd(diff)
  c(lower = bias - 1.96 * sd_diff, upper = bias + 1.96 * sd_diff)
}

# -------------------------------------------------------------------- load data

manual <- read_csv(file.path(BASE_DIR, "cmp/manual_sos_TEMPLATE.csv"),
                   show_col_types = FALSE) %>%
  filter(!is.na(sos_doy), nzchar(trimws(site))) %>%
  mutate(site_key = nk(site)) %>%
  select(site_key, manual_sos = sos_doy, manual_baseline_evi = baseline_evi)

pelt_mean <- read_csv(file.path(OUT_DIR, "pelt_transitions_mean_modis.csv"),
                      show_col_types = FALSE) %>%
  mutate(site_key = nk(site_key)) %>%
  select(site_key, pelt_sos = sos_doy, pelt_baseline_evi = baseline_evi,
         veg_type)

# ------------------------------------------------------------------- join data

joined <- manual %>%
  inner_join(pelt_mean, by = "site_key") %>%
  mutate(
    diff = manual_sos - pelt_sos,
    abs_diff = abs(diff)
  )

if (nrow(joined) != nrow(manual)) {
  stop("Join mismatch: manual has ", nrow(manual), " rows, joined has ", nrow(joined))
}

# -------------------------------------------------------------- exclusion sets

# Structural / snow-affected sites where methods legitimately diverge
structural_snow <- c("SJER", "SRER", "BARR", "TOOL", "RMNP", "BONA", "DEJU")

# Known bimodal sites
bimodal <- c("MOAB", "STER")

# Three nested sets
all_sites <- joined
ex_structural <- joined %>% filter(!site_key %in% structural_snow)
ex_both <- joined %>% filter(!site_key %in% c(structural_snow, bimodal))

# ---------------------------------------------------------- scoring function

score_set <- function(df, label) {
  n <- nrow(df)
  v <- df$diff[is.finite(df$diff)]
  n_valid <- length(v)

  if (n_valid == 0) {
    cat(sprintf("%s: n=0 (no valid pairs)\n\n", label))
    return(NULL)
  }

  bias <- mean(v)
  med_abs <- median(abs(v))
  mean_abs <- mean(abs(v))
  w8 <- sum(abs(v) <= 8)
  w16 <- sum(abs(v) <= 16)
  pct_w8 <- 100 * w8 / n_valid
  pct_w16 <- 100 * w16 / n_valid

  ccc <- ccc_fun(df$pelt_sos, df$manual_sos)
  loa <- bland_altman_loa(df$manual_sos, df$pelt_sos)

  cat(sprintf("%s\n", label))
  cat(sprintf("  n:                 %d\n", n_valid))
  cat(sprintf("  bias (man - PELT): %+.2f d\n", bias))
  cat(sprintf("  median |error|:    %.2f d\n", med_abs))
  cat(sprintf("  mean |error|:      %.2f d\n", mean_abs))
  cat(sprintf("  within 8 d:        %d (%.0f%%)\n", w8, pct_w8))
  cat(sprintf("  within 16 d:       %d (%.0f%%)\n", w16, pct_w16))
  cat(sprintf("  Lin's CCC:         %.3f\n", ccc))
  cat(sprintf("  Bland-Altman LOA:  [%.2f, %.2f] d\n", loa["lower"], loa["upper"]))
  cat("\n")

  data.frame(
    set = label,
    n = n_valid,
    bias = bias,
    median_abs_error = med_abs,
    mean_abs_error = mean_abs,
    within_8d = w8,
    pct_within_8d = pct_w8,
    within_16d = w16,
    pct_within_16d = pct_w16,
    ccc = ccc,
    loa_lower = loa["lower"],
    loa_upper = loa["upper"],
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------ main report

cat("==============================================================================\n")
cat("COMPREHENSIVE SOS SCORING: Manual vs. PELT\n")
cat("==============================================================================\n\n")

cat("Sign convention: diff = manual - PELT\n")
cat("  (positive bias = PELT earlier than manual)\n\n")

cat("------------------------------------------------------------------------------\n")
cat("SET A: ALL SITES (n=", nrow(all_sites), ")\n")
cat("------------------------------------------------------------------------------\n\n")
result_a <- score_set(all_sites, "All sites")

cat("------------------------------------------------------------------------------\n")
cat("SET B: EXCLUDING STRUCTURAL/SNOW SITES\n")
cat("  (excluded: ", paste(structural_snow, collapse = ", "), ")\n")
cat("  n=", nrow(ex_structural), "\n")
cat("------------------------------------------------------------------------------\n\n")
result_b <- score_set(ex_structural, "Ex structural/snow")

cat("------------------------------------------------------------------------------\n")
cat("SET C: EXCLUDING STRUCTURAL/SNOW + BIMODAL\n")
cat("  (excluded: ", paste(c(structural_snow, bimodal), collapse = ", "), ")\n")
cat("  n=", nrow(ex_both), "\n")
cat("------------------------------------------------------------------------------\n\n")
result_c <- score_set(ex_both, "Ex structural/snow + bimodal")

# -------------------------------------------------------- breakdown by veg type

cat("==============================================================================\n")
cat("BREAKDOWN BY VEGETATION TYPE\n")
cat("==============================================================================\n\n")

veg_breakdown <- ex_both %>%
  group_by(veg_type) %>%
  summarise(
    n = n(),
    bias = mean(diff, na.rm = TRUE),
    median_abs_error = median(abs(diff), na.rm = TRUE),
    mean_abs_error = mean(abs(diff), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(veg_type)

cat(sprintf("%-8s %3s %10s %16s %14s\n",
            "veg_type", "n", "bias", "median|error|", "mean|error|"))
cat(paste(rep("-", 60), collapse = ""), "\n")
for (i in seq_len(nrow(veg_breakdown))) {
  r <- veg_breakdown[i, ]
  cat(sprintf("%-8s %3d %+10.2f %16.2f %14.2f\n",
              r$veg_type, r$n, r$bias, r$median_abs_error, r$mean_abs_error))
}

cat("\n")

# ------------------------------------------------------------------ ANOVA on veg type

cat("==============================================================================\n")
cat("ONE-WAY ANOVA: SOS error ~ vegetation type (Set C only)\n")
cat("==============================================================================\n\n")

if (nrow(ex_both) > 0 && length(unique(ex_both$veg_type)) > 1) {
  anova_data <- ex_both %>% filter(is.finite(diff))
  if (nrow(anova_data) >= 3) {
    anova_fit <- aov(diff ~ veg_type, data = anova_data)
    anova_summary <- summary(anova_fit)
    print(anova_summary)
    cat("\n")

    p_value <- anova_summary[[1]][["Pr(>F)"]][1]
    cat(sprintf("Interpretation: %s\n",
                ifelse(is.na(p_value), "NA",
                       ifelse(p_value < 0.05,
                              sprintf("SIGNIFICANT vegetation-type effect (p = %.4f)", p_value),
                              sprintf("NO significant vegetation-type effect (p = %.4f)", p_value)))))
    cat("\n")
  } else {
    cat("Insufficient data for ANOVA (n < 3)\n\n")
  }
} else {
  cat("Insufficient vegetation types for ANOVA\n\n")
}

# -------------------------------------------------- comparison to EOS/SOF benchmarks

cat("==============================================================================\n")
cat("COMPARISON TO ESTABLISHED EOS/SOF AGREEMENT\n")
cat("==============================================================================\n\n")

cat("Established benchmarks (mean-curve arm, ex-outliers):\n")
cat("  EOS: median |error| 0.59 d,  CCC 0.998\n")
cat("  SOF: median |error| 2.79 d,  CCC 0.868 (0.751 ex-outliers)\n\n")

sos_med <- result_c$median_abs_error
sos_ccc <- result_c$ccc

cat("SOS (this analysis, Set C = ex structural/snow + bimodal):\n")
cat(sprintf("  median |error|: %.2f d\n", sos_med))
cat(sprintf("  CCC:            %.3f\n\n", sos_ccc))

cat("Relative performance:\n")
cat(sprintf("  SOS vs EOS: median |e| %.2fx larger,  CCC %s\n",
            sos_med / 0.59,
            ifelse(sos_ccc > 0.998, "better", ifelse(sos_ccc < 0.998, "worse", "equal"))))
cat(sprintf("  SOS vs SOF: median |e| %.2fx %s, CCC %s\n",
            abs(sos_med / 2.79),
            ifelse(sos_med < 2.79, "smaller", "larger"),
            ifelse(sos_ccc > 0.868, "better", ifelse(sos_ccc < 0.868, "worse", "equal"))))

cat("\n")

# ------------------------------------------------------------------ summary table

cat("==============================================================================\n")
cat("SUMMARY TABLE (all three sets)\n")
cat("==============================================================================\n\n")

summary_table <- bind_rows(result_a, result_b, result_c)
cat(sprintf("%-30s %4s %7s %11s %9s %8s %9s %6s\n",
            "set", "n", "bias", "median|e|", "mean|e|", "w/in 8d", "w/in 16d", "CCC"))
cat(paste(rep("-", 100), collapse = ""), "\n")
for (i in seq_len(nrow(summary_table))) {
  r <- summary_table[i, ]
  cat(sprintf("%-30s %4d %+7.2f %11.2f %9.2f %7.0f%% %8.0f%% %6.3f\n",
              r$set, r$n, r$bias, r$median_abs_error, r$mean_abs_error,
              r$pct_within_8d, r$pct_within_16d, r$ccc))
}

cat("\n")
cat("Done.\n")
