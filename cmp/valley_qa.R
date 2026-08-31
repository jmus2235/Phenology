# ==============================================================================
# valley_qa.R
#
# QA check for manual SOS dates that land BEFORE the pre-peak valley minimum.
# Discovered at TEAK and OSBS where the 2nd-derivative method found a SOS DOY
# that was 2-3 composites too early (before the bottom of the dormant valley),
# though the baseline EVI value was nearly correct because the valley is flat.
#
# For each site, report:
# - valley_doy: DOY of minimum EVI before the peak
# - manual_sos, PELT sos_doy
# - sos_minus_valley for both (negative = defect signature)
# - flatness metrics: EVI at manual SOS, at valley bottom, and their difference
#
# CAVEAT: Winter-rainfall/snow sites (SJER, BARR, TOOL, RMNP, BONA, DEJU) have
# dormant valleys in the PREVIOUS calendar year, so pre-peak valley is not
# meaningful in DOY-indexed space. These are flagged but not counted as defects.
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

nk <- function(x) gsub("[-/]", "_", trimws(x))

BASE_DIR <- "C:/Users/jmusinsky/Documents/R_Scripts/MODIS"
OUT_DIR  <- file.path(BASE_DIR, "pelt_out")

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
         peak_doy, veg_type)

curves <- read_csv(file.path(OUT_DIR, "pelt_mean_curves_modis.csv"),
                   show_col_types = FALSE) %>%
  mutate(site_key = nk(site_key))

# ----------------------------------------------------------------- join & check

if (!all(manual$site_key %in% pelt_mean$site_key)) {
  stop("manual sites not in PELT: ",
       paste(setdiff(manual$site_key, pelt_mean$site_key), collapse = ", "))
}
if (!all(manual$site_key %in% curves$site_key)) {
  stop("manual sites not in curves: ",
       paste(setdiff(manual$site_key, curves$site_key), collapse = ", "))
}

sites <- manual %>%
  inner_join(pelt_mean, by = "site_key") %>%
  arrange(site_key)

# -------------------------------------------------------- valley analysis by site

valley_results <- lapply(sites$site_key, function(sk) {
  row <- sites[sites$site_key == sk, ]
  crv <- curves %>% filter(site_key == sk) %>% arrange(doy)

  peak_doy_val <- row$peak_doy
  manual_sos_val <- row$manual_sos
  pelt_sos_val <- row$pelt_sos

  # Find pre-peak valley: minimum EVI before peak_doy
  pre_peak <- crv %>% filter(doy < peak_doy_val)
  if (nrow(pre_peak) == 0) {
    return(data.frame(
      site_key = sk,
      veg_type = row$veg_type,
      valley_doy = NA_real_,
      valley_evi = NA_real_,
      manual_sos = manual_sos_val,
      pelt_sos = pelt_sos_val,
      peak_doy = peak_doy_val,
      manual_sos_minus_valley = NA_real_,
      pelt_sos_minus_valley = NA_real_,
      evi_at_manual_sos = NA_real_,
      evi_at_valley = NA_real_,
      flatness_diff = NA_real_,
      note = "no_pre_peak_obs",
      stringsAsFactors = FALSE
    ))
  }

  valley_idx <- which.min(pre_peak$fitted)
  valley_doy_val <- pre_peak$doy[valley_idx]
  valley_evi_val <- pre_peak$fitted[valley_idx]

  # Find EVI at manual SOS by linear interpolation if needed
  evi_at_manual <- if (manual_sos_val %in% crv$doy) {
    crv$fitted[crv$doy == manual_sos_val]
  } else {
    # Linear interpolation
    lo <- crv %>% filter(doy < manual_sos_val) %>% arrange(desc(doy)) %>% slice(1)
    hi <- crv %>% filter(doy > manual_sos_val) %>% arrange(doy) %>% slice(1)
    if (nrow(lo) == 0 || nrow(hi) == 0) NA_real_
    else lo$fitted + (manual_sos_val - lo$doy) * (hi$fitted - lo$fitted) / (hi$doy - lo$doy)
  }

  data.frame(
    site_key = sk,
    veg_type = row$veg_type,
    valley_doy = valley_doy_val,
    valley_evi = valley_evi_val,
    manual_sos = manual_sos_val,
    pelt_sos = pelt_sos_val,
    peak_doy = peak_doy_val,
    manual_sos_minus_valley = manual_sos_val - valley_doy_val,
    pelt_sos_minus_valley = pelt_sos_val - valley_doy_val,
    evi_at_manual_sos = evi_at_manual,
    evi_at_valley = valley_evi_val,
    flatness_diff = evi_at_manual - valley_evi_val,
    note = "",
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

# -------------------------------------------------- flag structural/snow sites
# Winter-edge/snow sites where the pre-peak valley is not meaningful in
# DOY-indexed space (dormant valley in previous calendar year)
structural_snow <- c("SJER", "SRER", "BARR", "TOOL", "RMNP", "BONA", "DEJU")

valley_results <- valley_results %>%
  mutate(
    is_structural_snow = site_key %in% structural_snow,
    flag = case_when(
      note == "no_pre_peak_obs" ~ "no_pre_peak_obs",
      is_structural_snow ~ "winter_edge_not_interpretable",
      manual_sos_minus_valley < 0 ~ "DEFECT_manual_before_valley",
      TRUE ~ ""
    )
  ) %>%
  arrange(manual_sos_minus_valley)

# ---------------------------------------------------------------- report output

cat("============================================================\n")
cat("VALLEY QA: Manual SOS vs. pre-peak valley minimum\n")
cat("============================================================\n\n")

cat("Full table (ranked by manual_sos_minus_valley ascending):\n\n")
cat(sprintf("%-10s %-4s %6s %9s | %6s %6s %6s | %8s %8s | %11s %11s | %9s | %s\n",
            "site", "veg", "valley", "valley_EVI", "man_SOS", "PELT", "peak",
            "man-vly", "pelt-vly", "EVI@man", "EVI@vly", "flat_diff", "flag"))
cat(paste(rep("-", 150), collapse = ""), "\n")

for (i in seq_len(nrow(valley_results))) {
  r <- valley_results[i, ]
  cat(sprintf("%-10s %-4s %6.0f %9.4f | %6.0f %6.0f %6.0f | %+8.0f %+8.0f | %11.4f %11.4f | %+9.4f | %s\n",
              r$site_key,
              ifelse(is.na(r$veg_type), "NA", r$veg_type),
              r$valley_doy,
              r$valley_evi,
              r$manual_sos,
              r$pelt_sos,
              r$peak_doy,
              r$manual_sos_minus_valley,
              r$pelt_sos_minus_valley,
              r$evi_at_manual_sos,
              r$evi_at_valley,
              r$flatness_diff,
              r$flag))
}

cat("\n")
cat("============================================================\n")
cat("SUMMARY\n")
cat("============================================================\n\n")

n_total <- nrow(valley_results)
n_structural <- sum(valley_results$is_structural_snow, na.rm = TRUE)
n_defect <- sum(valley_results$flag == "DEFECT_manual_before_valley", na.rm = TRUE)
n_interpretable <- n_total - n_structural

cat(sprintf("Total sites:                               %d\n", n_total))
cat(sprintf("  Winter-edge/snow (not interpretable):    %d (%s)\n",
            n_structural, paste(structural_snow, collapse = ", ")))
cat(sprintf("  Interpretable sites:                     %d\n", n_interpretable))
cat(sprintf("  DEFECT flagged (manual before valley):   %d\n", n_defect))

if (n_defect > 0) {
  defect_sites <- valley_results %>%
    filter(flag == "DEFECT_manual_before_valley") %>%
    select(site_key, manual_sos, valley_doy, manual_sos_minus_valley, flatness_diff)
  cat("\nDefect sites:\n")
  print(as.data.frame(defect_sites), row.names = FALSE)
}

cat("\n")
cat("Interpretable-only stats (excluding winter-edge/snow sites):\n")
interp <- valley_results %>% filter(!is_structural_snow)
cat(sprintf("  median manual_sos_minus_valley:  %+.0f d\n",
            median(interp$manual_sos_minus_valley, na.rm = TRUE)))
cat(sprintf("  median pelt_sos_minus_valley:    %+.0f d\n",
            median(interp$pelt_sos_minus_valley, na.rm = TRUE)))
cat(sprintf("  median flatness_diff:            %+.4f EVI\n",
            median(interp$flatness_diff, na.rm = TRUE)))

cat("\nDone.\n")
