# ============================================================================
# pelt_sos_transitions_modis.R
#
# METHOD COMPARISON, not a production pipeline.  This is a straight PORT of
# the ALGORITHM in
#   C:/Users/jmusinsky/Documents/R_Scripts/HLS/pelt_sos_transitions.R
# onto MODIS EVI curves, so the single question it answers is: does PELT
# reproduce the author's manual second-derivative EOS/SOF dates on MODIS, the
# way it already does on HLS?  All peak-finding / changepoint / threshold-
# crossing code below is unchanged from the HLS script (same functions, same
# constants) -- only the INPUT LOADING and VEGETATION-TYPE LOOKUP differ,
# because the MODIS input is a different shape than the HLS "Detailed" files,
# and a VALIDATION section is added because HLS had no per-site-year manual
# reference to check against.
#
# Definitions (identical to HLS; author decision, 2026-08-25):
#   SOS       = DOY of the first PELT greenup break on the rising limb
#   baseline  = MEDIAN of EVI over the pre-break segment (phenocamr
#               convention -- robust to exactly where in the dormant plateau
#               the break lands)
#   amplitude = peak - baseline
#   EOS       = first rising-limb crossing of baseline + 0.90 * amplitude
#               (start of the peak-greenness plateau)
#   SOF       = first post-peak crossing below baseline + f * amplitude
#               (end of the peak-greenness plateau),
#               f = 0.80 for DB/EB/EN, 0.90 for AG/GR/SH/TN
# EOS and SOF are linearly interpolated between 8-day composites.  SOS is NOT
# -- it is quantised to the 8-day grid.
#
# Seasonal-prior peak selection (identical to HLS; OFF BY DEFAULT -- screened,
# implemented, and scored NET WORSE on EOS/SOF; see the toggle block below for
# the numbers and the mechanism.  Opt in with PELT_SEASONAL_PRIOR=1):
#   find_peaks() identifies prominent, well-separated candidates (MIN_PEAK_SEP
#   and MIN_PEAK_PROM unchanged).  Instead of taking the GLOBAL MAXIMUM, take
#   the candidate at smallest CIRCULAR DOY distance to the site's own
#   leave-one-out circular-median peak DOY.  No new tuning constant -- the
#   prior is the site's own central tendency.  LEAVE-ONE-OUT prevents
#   self-confirmation: the year being judged does not contribute to its own
#   prior.  Circular distance is required because DOY 5 and DOY 360 are 10 days
#   apart, not 355.  The prior is built from GLOBAL-MAX peaks (not
#   already-re-picked peaks) so it reflects the data, not the rule's own
#   feedback.  At a tie (two candidates equally distant), keep the higher-EVI
#   one (stable to current behaviour).  Inert when n_peaks == 1 (the rule can
#   only choose among detected candidates; the 5 known cases where the true
#   peak was never DETECTED are a prominence problem, out of scope).  Reverts
#   to global max if fewer than 2 other years' peaks exist, and -- by default --
#   because the rule is disabled; enable per arm with PELT_SEASONAL_PRIOR=1
#   (annual arm) or PELT_SEASONAL_PRIOR_MEAN=1 (mean curve only).
#
# INPUT DIFFERENCE FROM HLS (see task brief for the audit that justifies this):
#   MODIS\data_out\ (528 files / 184 ambiguous site-units) is NOT used.
#   Input is the PhenoFlight DB export
#     .../PhenoFlight_DB_Exports/modis_evi_timeseries.csv
#   which is already LONG format (site, year, doy, evi_value 0-1 scaled,
#   multiyear_mean_evi), one curve per site, 44 sites x 22 yr (2003-2024) x
#   46 DOY (8-day cadence), zero duplicate (site, year, doy) keys.  No file
#   discovery / selection logic is needed, unlike HLS.
#
# DATA-QUALITY DIFFERENCE FROM HLS:
#   The export is a perfect rectangle (every DOY x every year), so there is no
#   NA structure that would show which DOYs were actually observed vs.
#   LOESS-extrapolated.  HLS's TRIM_EXTRAP (which needs EVI_RAW to find the
#   observed span) CANNOT be ported -- there is no EVI_RAW here.  This is a
#   KNOWN LIMITATION, left unaddressed by design (join to MODIS\data\ raw GEE
#   exports to reconstruct observation spans is explicitly out of scope).
#   Instead the per-DOY robust filter (ROBUST_K / ROBUST_MAD_FLOOR), which
#   compares each year's value at a DOY against the across-year median AT
#   THAT DOY and needs no raw data, is the PRIMARY guard here (in the HLS
#   script it is an optional secondary guard, off by default, because
#   TRIM_EXTRAP did the main work there).  A known extrapolation artefact is
#   present in this export: 0.30% of values are < -0.2, all at the winter
#   edges (DOY <=65 or >=297); BARR's mean curve alone carries 30 negative
#   DOYs (min -0.303), TOOL 27, HEAL 11, DEJU 8, YELL 1.
#
# Usage:
#   Rscript pelt_sos_transitions_modis.R                 # all sites
#   Rscript pelt_sos_transitions_modis.R --site=ABBY      # one site
#
# Outputs (pelt_out/):
#   pelt_transitions_mean_modis.csv    one row per site, mean-curve arm
#   pelt_transitions_annual_modis.csv  one row per site-year (2003-2024)
#   pelt_modis_comparison.csv          both validation arms, joined, w/ diffs
# ============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(changepoint)
  library(ggplot2)
})

# ---------------------------------------------------------------- configuration

BASE_DIR   <- "C:/Users/jmusinsky/Documents/R_Scripts/MODIS"
OUT_DIR    <- if (nzchar(Sys.getenv("PELT_OUT_DIR"))) Sys.getenv("PELT_OUT_DIR") else
                file.path(BASE_DIR, "pelt_out")
EXPORT_DIR <- file.path("C:/Users/jmusinsky/Documents/Flight Operations",
                        "Flight Season Planning/Phenology/PhenoFlight_DB_Exports")
TS_FILE     <- file.path(EXPORT_DIR, "modis_evi_timeseries.csv")
MEAN_REF    <- file.path(EXPORT_DIR, "modis_mean_transitions.csv")
ANNUAL_REF  <- file.path(EXPORT_DIR, "modis_annual_transitions.csv")

# ---- seasonal-prior peak selection ------------------------------------------
# Toggles for the seasonal-prior rule.  Both default OFF -- opt in with
# PELT_SEASONAL_PRIOR=1 / PELT_SEASONAL_PRIOR_MEAN=1.  Identical to HLS.
#
# The rule passed its peak-SELECTION screen but FAILED when scored on EOS and
# SOF: mean absolute error against the manual record rose on 3 of 4 sensor x
# limb measures (MODIS SOF 7.74 -> 8.26, HLS EOS 4.59 -> 5.17, HLS SOF
# 7.58 -> 8.04; only MODIS EOS improved, 4.42 -> 4.34).  The cause is
# structural: at the crossing_doy() calls in derive_transitions() the peak index
# is a SEARCH-INTERVAL BOUND, not a scalar offset, so on a curve with two maxima
# the same threshold is crossed four times and the peak alone decides which pair
# is reported.  Date error is therefore quantised to the mode separation, and
# wins and losses use the identical operation -- no tuning keeps one and drops
# the other.  Kept here to reproduce the negative result.
SEASONAL_PRIOR      <- !identical(Sys.getenv("PELT_SEASONAL_PRIOR", "0"), "0")
SEASONAL_PRIOR_MEAN <- !identical(Sys.getenv("PELT_SEASONAL_PRIOR_MEAN", "0"), "0")
YEAR_LEN <- 366   # circular distance / circular median operate on this circle

COMPOSITE_DAYS <- 8           # cadence of the input series (unchanged from HLS)
MIN_SEG_OBS    <- 3           # PELT minseglen, in observations (24 d)
PEN_DEFAULT    <- 0.2         # non-DB veg types
PEN_DB         <- 0.5         # DB, per project convention
MIN_PEAK_SEP   <- 45          # days
MIN_PEAK_PROM  <- 0.15        # a 2nd mode must rise 15% of the curve's range
MIN_YEAR_OBS   <- 12          # per-year arm only; mean curve always pools all years
GREENUP_RISE_FRAC <- 0.10
MIN_LIMB_OBS   <- 6

# ---- fitted-value plausibility mask (ports unchanged from HLS) ------------
# Kept as a mechanism (disabled by default, same as HLS) in case an absolute
# bound is ever wanted on top of the per-DOY robust filter below.
EVI_MIN_PLAUSIBLE <- as.numeric(Sys.getenv("PELT_EVI_MIN", "-Inf"))
EVI_MAX_PLAUSIBLE <- as.numeric(Sys.getenv("PELT_EVI_MAX", "Inf"))

# Reject a year's value at a DOY if it sits more than ROBUST_K scaled-MADs
# from the across-year median at that DOY.  THIS IS THE PRIMARY GUARD for
# MODIS (see header comment) -- default is ENABLED here, unlike HLS where it
# defaulted to Inf (off) because TRIM_EXTRAP did the work instead.
ROBUST_K          <- as.numeric(Sys.getenv("PELT_ROBUST_K", "3"))
ROBUST_MAD_FLOOR  <- 0.05
MIN_YRS_PER_DOY   <- as.integer(Sys.getenv("PELT_MIN_YRS_PER_DOY", "3"))
DROP_YEARS        <- {
  v <- Sys.getenv("PELT_DROP_YEARS", "")
  if (nzchar(v)) as.integer(strsplit(v, ",")[[1]]) else integer(0)
}

PG_STRT_FRAC   <- 0.90        # EOS threshold; 90 at every site
SOF_FRAC_BY_VEG <- c(DB = 0.80, EB = 0.80, EN = 0.80,
                     AG = 0.90, GR = 0.90, SH = 0.90, TN = 0.90)

# Satellite ROI vegetation differs from the PhenoCam-side label at CLBJ only
# (export says GR because that reflects the PhenoCam, not the satellite ROI;
# author confirmed DB is correct for the MODIS ROI).  Unlike HLS, NIWO/TOOL
# need NO override here -- the export's own veg_type is already correct for
# those.
VEG_TYPE_OVERRIDES <- c(CLBJ = "DB")

KNOWN_OUTLIER_SITES <- c("STER", "SJER", "LAJA")  # context only, see header

# ---- Calendar-year framing guard -------------------------------------------
# Author, 2026-09-01: "SJER is an exception where PELT should not be used in
# its current code formulation ... the baseline for SOS calculation occurs in
# the prior year to the peak greenness window."
#
# Everything downstream of the changepoint -- baseline, amplitude, and both
# thresholds -- assumes the dormant floor preceding a peak is present in the
# same DOY 1-365 slice.  At a Mediterranean winter-annual site it is not: the
# floor belonging to the spring peak fell in the previous December, so the
# detector substitutes a mid-season value for dormancy, inflating the baseline
# and halving the amplitude.  This is a coordinate-system failure, not a
# detector failure, and the 2nd-derivative method fails on it identically.
#
# Detection (see cmp/yearframe_screen_modis.R, and the HLS twin for the two
# rejected alternatives -- both of which flagged the Arctic/boreal snow sites,
# where winter fitted EVI is negative, while MISSING SJER):
#
#   baseline_gap = ( min(EVI over DOY 1..peak) - min(EVI over the year) )
#                  / ( max - min )
#
# 0 at a conventional site.  Measured 0.481 at SJER here and 0.522 on HLS --
# the only site over 0.25 on either sensor, next highest 0.202 and 0.238.
# Above the cut the transition dates are withheld (NA + flag), not silently
# emitted: PELT as currently formulated has no defensible answer here.
BASELINE_GAP_MAX <- as.numeric(Sys.getenv("PELT_BASELINE_GAP_MAX", "0.25"))

# ---- one-off experiment: does matching the mean curve's pooling window to
# the manual reference's 2015-2024 window change headline-arm agreement? ----
# modis_mean_transitions.csv carries num_years=10 (author's dates are built
# from 2015-2024 only); the DEFAULT mean curve above pools the full
# 2003-2024 record instead. This parameter lets that window be matched
# without touching the full-record behaviour, which remains the default
# reported in pelt_transitions_mean_modis.csv / pelt_modis_comparison.csv.
MEAN_CURVE_YEARS <- as.integer(Sys.getenv("PELT_MEAN_CURVE_YEARS_START", "2015")):
                     as.integer(Sys.getenv("PELT_MEAN_CURVE_YEARS_END",   "2024"))

# ---------------------------------------------------------------------- helpers

#' Circular distance between two DOYs, in days, on a YEAR_LEN circle.
#' Verbatim from HLS/cmp/seasonal_prior_test.R (YEAR_LEN 366).
circ_dist <- function(a, b) {
  d <- abs(a - b) %% YEAR_LEN
  pmin(d, YEAR_LEN - d)
}

#' Circular median: the candidate value minimising total circular distance to
#' the sample.  Exact (searches the sample itself), robust, and needs no
#' distributional assumption.  The arithmetic median is wrong on a circle
#' because DOY 5 and DOY 360 are 10 days apart, not 355.
#' Verbatim from HLS/cmp/seasonal_prior_test.R.
circ_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  if (length(x) == 1) return(x[1])
  tot <- vapply(x, function(c0) sum(circ_dist(x, c0)), numeric(1))
  x[which.min(tot)]
}

args        <- commandArgs(trailingOnly = TRUE)
opt_plots   <- !("--no-plots" %in% args)
opt_dump    <- "--dump-curves" %in% args
opt_site    <- sub("^--site=", "", grep("^--site=", args, value = TRUE))
if (length(opt_site) == 0) opt_site <- NULL

msg <- function(...) cat(sprintf(...), "\n", sep = "")
normalize_key <- function(x) gsub("[-/]", "_", x)

#' Local maxima separated by >= min_sep days AND topographically prominent.
#' Verbatim port from HLS pelt_sos_transitions.R.
find_peaks <- function(doy, y, min_sep = MIN_PEAK_SEP,
                       min_prom_frac = MIN_PEAK_PROM) {
  n <- length(y)
  if (n < 3) return(list(primary_idx = which.max(y), n_peaks = 1L,
                         second_doy = NA_real_, second_evi = NA_real_))

  loc <- which(vapply(seq_len(n), function(i) {
    lo <- max(1, i - 1); hi <- min(n, i + 1)
    y[i] >= y[lo] && y[i] >= y[hi]
  }, logical(1)))
  if (length(loc) == 0) loc <- which.max(y)

  rng <- diff(range(y, na.rm = TRUE))
  prom <- vapply(loc, function(i) {
    higher <- loc[y[loc] > y[i]]
    if (length(higher) == 0) return(y[i] - min(y, na.rm = TRUE))
    l <- higher[higher < i]; r <- higher[higher > i]
    lcol <- if (length(l)) min(y[max(l):i]) else min(y[1:i])
    rcol <- if (length(r)) min(y[i:min(r)]) else min(y[i:n])
    y[i] - max(lcol, rcol)
  }, numeric(1))

  cand <- loc[is.finite(prom) & prom >= min_prom_frac * rng]
  if (length(cand) == 0) cand <- loc[which.max(y[loc])]

  cand <- cand[order(y[cand], decreasing = TRUE)]
  picked <- integer(0)
  for (i in cand) {
    if (all(abs(doy[i] - doy[picked]) >= min_sep)) picked <- c(picked, i)
  }
  list(
    primary_idx = picked[1],
    n_peaks     = length(picked),
    second_doy  = if (length(picked) >= 2) doy[picked[2]] else NA_real_,
    second_evi  = if (length(picked) >= 2) y[picked[2]]   else NA_real_,
    cand_idx    = picked   # full candidate set, decreasing EVI order
  )
}

#' First PELT greenup break on a rising limb.  Verbatim port from HLS.
pelt_break <- function(y, pen_value, min_seg = MIN_SEG_OBS) {
  n <- length(y)
  if (n < MIN_LIMB_OBS || !is.finite(sd(y)) || sd(y) == 0) {
    return(list(idx = NA_integer_, n_breaks = 0L, ok = FALSE))
  }
  z <- (y - mean(y)) / sd(y)
  fit <- try(changepoint::cpt.mean(
    as.numeric(z),
    method          = "PELT",
    test.stat       = "Normal",
    penalty         = "Manual",
    pen.value       = pen_value,
    minseglen       = min_seg,
    param.estimates = TRUE
  ), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(list(idx = NA_integer_, n_breaks = 0L, ok = FALSE))
  }
  brk <- changepoint::cpts(fit)
  brk <- brk[brk > 0 & brk < n]
  if (length(brk) == 0) return(list(idx = NA_integer_, n_breaks = 0L, ok = FALSE))

  segs  <- c(0, brk, n)
  means <- vapply(seq_len(length(segs) - 1),
                  function(k) mean(y[(segs[k] + 1):segs[k + 1]]), numeric(1))
  dormant <- means[1]
  thr     <- dormant + GREENUP_RISE_FRAC * (max(y) - dormant)
  k       <- which(means[-1] >= thr)[1]
  if (is.na(k)) k <- which.max(diff(means))

  list(idx = as.integer(brk[k]), n_breaks = length(brk),
       n_trimmed = as.integer(k - 1L), ok = TRUE)
}

#' Linearly interpolated DOY crossing.  Verbatim port from HLS.
crossing_doy <- function(doy, y, thr, from, to, dir = c("up", "down")) {
  dir <- match.arg(dir)
  if (is.na(from) || is.na(to) || to <= from) return(NA_real_)
  idx <- from:to
  hit <- if (dir == "up") which(y[idx] >= thr) else which(y[idx] <= thr)
  if (length(hit) == 0) return(NA_real_)
  j <- idx[hit[1]]
  if (j == from) return(doy[from])
  y0 <- y[j - 1]; y1 <- y[j]
  if (!is.finite(y0) || !is.finite(y1) || y1 == y0) return(doy[j])
  doy[j - 1] + (thr - y0) / (y1 - y0) * (doy[j] - doy[j - 1])
}

#' Full SOS -> baseline -> amplitude -> EOS/SOF chain for one curve.
#' Verbatim port from HLS (plus seasonal-prior additions).
derive_transitions <- function(doy, y, sof_frac, pen_value,
                               min_seg = MIN_SEG_OBS, prior_doy = NA_real_) {
  flags <- character(0)
  # `extra` OVERRIDES the NA template.  It used to be c()'d onto it, which
  # appends rather than replaces: the result carried two `peak_doy` entries and
  # as.data.frame() renamed the second to `peak_doy.1`, so every diagnostic an
  # early return tried to preserve was silently dropped in favour of the NA.
  na_row <- function(extra = list()) {
    base <- list(peak_doy = NA_real_, peak_evi = NA_real_, n_peaks = NA_integer_,
           second_peak_doy = NA_real_, sos_doy = NA_real_,
           sos_break_index = NA_integer_, n_breaks = NA_integer_,
           n_breaks_trimmed = NA_integer_,
           baseline_evi = NA_real_, baseline_evi_at_sos = NA_real_,
           baseline_n_obs = NA_integer_, amplitude = NA_real_,
           thr_eos_evi = NA_real_, eos_doy = NA_real_,
           thr_sof_evi = NA_real_, sof_doy = NA_real_,
           window_days = NA_real_, baseline_gap = NA_real_,
           prior_doy = NA_real_, peak_doy_globalmax = NA_real_,
           prior_switched = NA_integer_, n_peak_cands = NA_integer_)
    base[names(extra)] <- extra
    base
  }

  if (length(y) < MIN_LIMB_OBS) {
    return(na_row(list(flags = "too_few_obs")))
  }
  if (any(y > EVI_MAX_PLAUSIBLE, na.rm = TRUE)) {
    flags <- c(flags, sprintf("evi_out_of_range_max%.2f", max(y, na.rm = TRUE)))
  }

  pk <- find_peaks(doy, y)
  pidx_globalmax <- pk$primary_idx
  if (is.na(pidx_globalmax)) return(na_row(list(flags = "no_peak",
                                                 prior_doy = prior_doy,
                                                 n_peak_cands = length(pk$cand_idx))))

  # Seasonal-prior peak selection: choose the candidate nearest to the prior.
  # Only switch from global max if: (1) the rule is enabled, (2) a finite prior
  # exists, and (3) there are at least 2 candidates to choose from.  which.min
  # returns the FIRST minimum, and cand_idx is ordered by decreasing EVI, so an
  # exact tie keeps the higher-EVI candidate (stable to current behaviour).
  pidx <- pidx_globalmax   # default: unchanged
  if (SEASONAL_PRIOR && is.finite(prior_doy) && length(pk$cand_idx) >= 2) {
    pidx <- pk$cand_idx[which.min(circ_dist(doy[pk$cand_idx], prior_doy))]
  }
  prior_switched <- as.integer(!identical(pidx, pidx_globalmax))
  if (pk$n_peaks > 1) flags <- c(flags, sprintf("multimodal_n%d", pk$n_peaks))

  if (pidx < MIN_LIMB_OBS) {
    return(na_row(list(flags = paste(c(flags, "rising_limb_too_short"),
                                     collapse = ";"),
                       prior_doy = prior_doy,
                       peak_doy_globalmax = doy[pidx_globalmax],
                       prior_switched = prior_switched,
                       n_peak_cands = length(pk$cand_idx))))
  }

  # Is the dormant floor belonging to this peak inside the calendar year at all?
  rng      <- diff(range(y, na.rm = TRUE))
  base_gap <- if (is.finite(rng) && rng > 0)
                (min(y[1:pidx], na.rm = TRUE) - min(y, na.rm = TRUE)) / rng
              else NA_real_
  if (is.finite(base_gap) && base_gap >= BASELINE_GAP_MAX) {
    return(na_row(list(peak_doy = doy[pidx], peak_evi = y[pidx],
                       n_peaks = pk$n_peaks, second_peak_doy = pk$second_doy,
                       baseline_gap = base_gap,
                       flags = paste(c(flags, sprintf("baseline_in_prior_year_gap%.2f",
                                                      base_gap)),
                                     collapse = ";"),
                       prior_doy = prior_doy,
                       peak_doy_globalmax = doy[pidx_globalmax],
                       prior_switched = prior_switched,
                       n_peak_cands = length(pk$cand_idx))))
  }

  limb <- 1:pidx
  br <- pelt_break(y[limb], pen_value, min_seg)
  if (!br$ok) {
    return(na_row(list(peak_doy = doy[pidx], peak_evi = y[pidx],
                       n_peaks = pk$n_peaks, second_peak_doy = pk$second_doy,
                       flags = paste(c(flags, "no_pelt_break"), collapse = ";"),
                       prior_doy = prior_doy,
                       peak_doy_globalmax = doy[pidx_globalmax],
                       prior_switched = prior_switched,
                       n_peak_cands = length(pk$cand_idx))))
  }

  sos_idx  <- br$idx
  baseline <- median(y[1:sos_idx])
  base_pt  <- y[sos_idx]
  peak_evi <- y[pidx]
  amp      <- peak_evi - baseline
  if (!is.finite(amp) || amp <= 0) flags <- c(flags, "nonpositive_amplitude")

  thr_eos <- baseline + PG_STRT_FRAC * amp
  thr_sof <- baseline + sof_frac     * amp

  eos <- crossing_doy(doy, y, thr_eos, sos_idx, pidx, "up")
  sof <- crossing_doy(doy, y, thr_sof, pidx, length(y), "down")
  if (is.na(eos)) flags <- c(flags, "no_eos_crossing")
  if (is.na(sof)) flags <- c(flags, "no_sof_crossing")

  list(
    peak_doy = doy[pidx], peak_evi = peak_evi, n_peaks = pk$n_peaks,
    second_peak_doy = pk$second_doy,
    sos_doy = doy[sos_idx], sos_break_index = sos_idx, n_breaks = br$n_breaks,
    n_breaks_trimmed = br$n_trimmed,
    baseline_evi = baseline, baseline_evi_at_sos = base_pt,
    baseline_n_obs = sos_idx, amplitude = amp,
    thr_eos_evi = thr_eos, eos_doy = eos,
    thr_sof_evi = thr_sof, sof_doy = sof,
    window_days = if (is.na(eos) || is.na(sof)) NA_real_ else sof - eos,
    baseline_gap = base_gap,
    prior_doy = prior_doy,
    peak_doy_globalmax = doy[pidx_globalmax],
    prior_switched = prior_switched,
    n_peak_cands = length(pk$cand_idx),
    flags = paste(flags, collapse = ";")
  )
}

#' Load and mask one site's per-year curves.  This REPLACES HLS's
#' load_site_curves(): the input is already tidy long-format, so there is no
#' file parsing, no FITTED/EVI_RAW split, and no TRIM_EXTRAP step (no raw
#' observation span to trim to -- see header comment).  The value-level range
#' mask and per-DOY robust mask are otherwise identical in spirit to HLS.
load_site_curves_modis <- function(site_df) {
  d <- site_df %>%
    select(year, doy, fitted = evi_value) %>%
    filter(!is.na(doy), !is.na(fitted)) %>%
    arrange(year, doy)

  n_in <- nrow(d)
  if (length(DROP_YEARS)) d <- d %>% filter(!year %in% DROP_YEARS)
  n_yr <- nrow(d)

  d <- d %>% filter(fitted >= EVI_MIN_PLAUSIBLE, fitted <= EVI_MAX_PLAUSIBLE)
  n_rng <- nrow(d)

  if (is.finite(ROBUST_K)) {
    d <- d %>% group_by(doy) %>%
      mutate(.med = median(fitted, na.rm = TRUE),
             .mad = mad(fitted, na.rm = TRUE)) %>%
      ungroup() %>%
      filter(abs(fitted - .med) <= ROBUST_K * pmax(.mad, ROBUST_MAD_FLOOR)) %>%
      select(-.med, -.mad)
  }

  attr(d, "n_dropped_year")   <- n_in - n_yr
  attr(d, "n_dropped_range")  <- n_yr - n_rng
  attr(d, "n_dropped_robust") <- n_rng - nrow(d)
  d
}

# ---------------------------------------------------------------------- plotting

plot_site <- function(site_key, curves, mean_curve, res, annual, out_path) {
  ann_lab <- function(v) if (is.na(v)) "NA" else sprintf("%.0f", v)
  sub <- sprintf(
    "%s | veg %s | SOF thr %.0f%% | pen %.2f, minseglen %d obs (%d d)  ||  SOS %s  EOS %s  SOF %s",
    site_key, res$veg_type, res$pg_end_thrshld, res$pen_value,
    res$min_seg_obs, res$min_seg_obs * COMPOSITE_DAYS,
    ann_lab(res$sos_doy), ann_lab(res$eos_doy), ann_lab(res$sof_doy))

  p <- ggplot() +
    geom_line(data = curves,
              aes(x = doy, y = fitted, group = year),
              colour = "grey62", linewidth = 0.3, alpha = 0.8) +
    geom_line(data = mean_curve,
              aes(x = doy, y = fitted),
              colour = "black", linewidth = 1.3)

  if (is.finite(res$baseline_evi)) {
    p <- p + geom_hline(yintercept = res$baseline_evi,
                        colour = "grey40", linetype = "dotted", linewidth = 0.4)
  }
  if (is.finite(res$thr_eos_evi)) {
    p <- p + geom_hline(yintercept = res$thr_eos_evi,
                        colour = "#2166ac", linetype = "dashed", linewidth = 0.4)
  }
  if (is.finite(res$thr_sof_evi) &&
      !isTRUE(all.equal(res$thr_sof_evi, res$thr_eos_evi))) {
    p <- p + geom_hline(yintercept = res$thr_sof_evi,
                        colour = "#b2182b", linetype = "dashed", linewidth = 0.4)
  }

  vl <- data.frame(
    doy   = c(res$sos_doy, res$eos_doy, res$sof_doy),
    label = c("SOS (PELT)", "EOS 90%", sprintf("SOF %.0f%%", res$pg_end_thrshld)),
    col   = c("#1b7837", "#2166ac", "#b2182b"),
    stringsAsFactors = FALSE
  )
  vl <- vl[is.finite(vl$doy), ]
  if (nrow(vl) > 0) {
    p <- p +
      geom_vline(data = vl, aes(xintercept = doy), colour = vl$col,
                 linewidth = 0.7) +
      geom_text(data = vl,
                aes(x = doy, y = -Inf,
                    label = sprintf("%s\n%.0f", label, doy)),
                colour = vl$col, vjust = -0.3, hjust = -0.05,
                size = 2.9, lineheight = 0.9)
  }

  # per-year PELT SOS values, so the interannual spread is visible
  ann_sos <- annual$sos_doy[is.finite(annual$sos_doy)]
  if (length(ann_sos) > 0) {
    p <- p + geom_rug(data = data.frame(doy = ann_sos),
                      aes(x = doy), sides = "b",
                      colour = "#1b7837", alpha = 0.8, linewidth = 0.6)
  }

  p <- p +
    scale_x_continuous(breaks = seq(0, 360, 30), limits = c(0, 366),
                       expand = expansion(mult = c(0.01, 0.01))) +
    labs(title = sprintf("%s - PELT-derived phenometrics (MODIS EVI, 8-day composites, 2003-2024)",
                         site_key),
         subtitle = sub,
         x = "Day of year", y = "EVI (fitted)",
         caption = paste("Thin grey = annual fitted curves; thick black = mean curve.",
                         "Green rug = per-year PELT SOS.",
                         "Dotted = baseline (pre-break segment median);",
                         "dashed = EOS/SOF thresholds.")) +
    theme_bw(base_size = 10) +
    theme(plot.subtitle = element_text(size = 8),
          plot.caption  = element_text(size = 7, colour = "grey35", hjust = 0))

  ggsave(out_path, p, width = 10, height = 6, dpi = 150)
}

# -------------------------------------------------------------------------- main

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

msg("Reading %s ...", TS_FILE)
ts_all <- read_csv(TS_FILE, show_col_types = FALSE, progress = FALSE)
sites  <- sort(unique(ts_all$site))
msg("  %d rows -> %d sites, %d-%d, %d distinct DOY",
    nrow(ts_all), length(sites), min(ts_all$year), max(ts_all$year),
    n_distinct(ts_all$doy))
if (!is.null(opt_site)) sites <- intersect(sites, opt_site)
if (length(sites) == 0) stop("no sites selected")

mean_rows   <- list()
annual_rows <- list()
curve_rows  <- list()
prior_audit_rows <- list()

for (site_key in sites) {
  site_df <- ts_all %>% filter(site == site_key)
  curves  <- load_site_curves_modis(site_df)
  if (nrow(curves) == 0) { msg("  %-10s SKIP (no usable years)", site_key); next }

  # ---- veg type and thresholds ---------------------------------------------
  # Unlike HLS, veg_type comes straight from the export (constant per site,
  # verified in the audit) -- no separate pheno_veg_types.csv join is needed
  # or performed.  Only CLBJ is overridden (satellite ROI is DB; the export
  # value reflects the PhenoCam FOV, GR).
  veg_type   <- site_df$veg_type[1]
  veg_source <- "timeseries"
  if (site_key %in% names(VEG_TYPE_OVERRIDES)) {
    veg_type   <- VEG_TYPE_OVERRIDES[[site_key]]
    veg_source <- "override"
  }
  sof_frac  <- if (!is.na(veg_type) && veg_type %in% names(SOF_FRAC_BY_VEG)) {
    SOF_FRAC_BY_VEG[[veg_type]]
  } else 0.90
  pen_value <- if (!is.na(veg_type) && veg_type == "DB") PEN_DB else PEN_DEFAULT

  # ---- mean curve (pooled across all masked years, per DOY) -----------------
  mc <- curves %>% group_by(doy) %>%
    summarise(fitted = mean(fitted, na.rm = TRUE), n_yr = n(), .groups = "drop") %>%
    arrange(doy)
  n_doy_all <- nrow(mc)
  mc <- mc %>% filter(n_yr >= MIN_YRS_PER_DOY)

  # ---- store the exact masked mean curve for --dump-curves -------------------
  if (opt_dump) {
    curve_rows[[site_key]] <- data.frame(site_key = site_key, doy = mc$doy,
                                         fitted = mc$fitted, stringsAsFactors = FALSE)
  }

  # ---- seasonal-prior two-pass: collect global-max peaks ---------------------
  yrs <- sort(unique(curves$year))
  pass1_peaks <- vapply(yrs, function(yy) {
    cy <- curves %>% filter(year == yy) %>% arrange(doy)
    if (nrow(cy) < MIN_YEAR_OBS) return(NA_real_)
    pk <- find_peaks(cy$doy, cy$fitted)
    if (is.na(pk$primary_idx)) return(NA_real_)
    cy$doy[pk$primary_idx]
  }, numeric(1))
  names(pass1_peaks) <- as.character(yrs)

  # mean-curve prior: circular median of ALL years' global-max peaks
  prior_mean <- if (SEASONAL_PRIOR_MEAN) circ_median(pass1_peaks) else NA_real_

  res <- derive_transitions(mc$doy, mc$fitted, sof_frac, pen_value,
                            prior_doy = prior_mean)

  # ---- per-year arm (pass 2: use leave-one-out priors) ----------------------
  # (2003-2024; joined against the 2015-2024 manual annual reference downstream)
  ann <- lapply(yrs, function(yy) {
    cy <- curves %>% filter(year == yy) %>% arrange(doy)
    # leave-one-out prior: circ_median of OTHER years' peaks
    prior_yr <- circ_median(pass1_peaks[names(pass1_peaks) != as.character(yy)])
    r  <- if (nrow(cy) >= MIN_YEAR_OBS) {
      derive_transitions(cy$doy, cy$fitted, sof_frac, pen_value,
                         prior_doy = prior_yr)
    } else {
      derive_transitions(numeric(0), numeric(0), sof_frac, pen_value,
                         prior_doy = prior_yr)
    }
    c(list(site_key = site_key, year = yy, veg_type = veg_type,
           pg_end_thrshld = sof_frac * 100, n_obs = nrow(cy),
           pen_value = pen_value, min_seg_obs = MIN_SEG_OBS,
           sos_resolution_days = COMPOSITE_DAYS), r)
  })
  ann_df <- bind_rows(lapply(ann, as.data.frame, stringsAsFactors = FALSE))
  annual_rows[[site_key]] <- ann_df

  # ---- seasonal-prior audit (one row per site-year + one for mean curve) ----
  for (j in seq_len(nrow(ann_df))) {
    a <- ann_df[j, ]
    cands_doy <- if (is.finite(a$n_peaks) && a$n_peaks >= 1 && is.finite(a$peak_doy_globalmax)) {
      paste(c(a$peak_doy_globalmax,
              if (is.finite(a$second_peak_doy)) a$second_peak_doy else NULL),
            collapse = ";")
    } else NA_character_
    n_yrs_prior <- sum(is.finite(pass1_peaks[names(pass1_peaks) != as.character(a$year)]))
    prior_audit_rows[[length(prior_audit_rows) + 1]] <- data.frame(
      site_key = site_key,
      year = a$year,
      arm = "annual",
      n_peak_cands = a$n_peak_cands,
      cand_doys = cands_doy,
      prior_doy = a$prior_doy,
      n_years_in_prior = n_yrs_prior,
      peak_doy_globalmax = a$peak_doy_globalmax,
      peak_doy_chosen = a$peak_doy,
      prior_switched = a$prior_switched,
      move_days = if (is.finite(a$peak_doy) && is.finite(a$peak_doy_globalmax))
                    circ_dist(a$peak_doy, a$peak_doy_globalmax) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  # mean curve row
  mean_cands_doy <- if (is.finite(res$n_peaks) && res$n_peaks >= 1 && is.finite(res$peak_doy_globalmax)) {
    paste(c(res$peak_doy_globalmax,
            if (is.finite(res$second_peak_doy)) res$second_peak_doy else NULL),
          collapse = ";")
  } else NA_character_
  prior_audit_rows[[length(prior_audit_rows) + 1]] <- data.frame(
    site_key = site_key,
    year = NA_integer_,
    arm = "mean",
    n_peak_cands = res$n_peak_cands,
    cand_doys = mean_cands_doy,
    prior_doy = res$prior_doy,
    n_years_in_prior = sum(is.finite(pass1_peaks)),
    peak_doy_globalmax = res$peak_doy_globalmax,
    peak_doy_chosen = res$peak_doy,
    prior_switched = res$prior_switched,
    move_days = if (is.finite(res$peak_doy) && is.finite(res$peak_doy_globalmax))
                  circ_dist(res$peak_doy, res$peak_doy_globalmax) else NA_real_,
    stringsAsFactors = FALSE
  )

  res$site_key            <- site_key
  res$veg_type            <- veg_type
  res$veg_source          <- veg_source
  res$pg_strt_thrshld     <- PG_STRT_FRAC * 100
  res$pg_end_thrshld      <- sof_frac * 100
  res$n_years             <- length(yrs)
  res$n_obs_mean_curve    <- nrow(mc)
  res$years               <- paste(range(yrs), collapse = "-")
  res$pen_value           <- pen_value
  res$min_seg_obs         <- MIN_SEG_OBS
  res$sos_resolution_days <- COMPOSITE_DAYS
  res$sos_annual_mean     <- mean(ann_df$sos_doy, na.rm = TRUE)
  res$sos_annual_sd       <- sd(ann_df$sos_doy,   na.rm = TRUE)
  res$n_years_sos         <- sum(is.finite(ann_df$sos_doy))
  res$n_pts_dropped_year  <- attr(curves, "n_dropped_year")
  res$n_pts_dropped_range <- attr(curves, "n_dropped_range")
  res$n_pts_dropped_robust<- attr(curves, "n_dropped_robust")
  res$min_yrs_per_doy     <- if (nrow(mc)) min(mc$n_yr) else NA_integer_
  res$n_doy_dropped_thin  <- n_doy_all - nrow(mc)
  res$yrs_at_baseline_doy <- if (is.finite(res$sos_doy))
                               min(mc$n_yr[mc$doy <= res$sos_doy]) else NA_integer_
  mean_rows[[site_key]] <- as.data.frame(res, stringsAsFactors = FALSE)

  if (opt_plots) {
    tryCatch(
      plot_site(site_key, curves, mc, res, ann_df,
                file.path(OUT_DIR, sprintf("%s_pelt.png", site_key))),
      error = function(e) {
        warning(sprintf("Plot failed for %s: %s", site_key, conditionMessage(e)))
      }
    )
  }

  msg("  %-10s veg %-4s (%s) SOF %2.0f%%  SOS %-5s EOS %-5s SOF %-5s  (annual SOS sd %.1f d) %s",
      site_key, ifelse(is.na(veg_type), "NA", veg_type), veg_source, sof_frac * 100,
      ifelse(is.na(res$sos_doy), "NA", sprintf("%.0f", res$sos_doy)),
      ifelse(is.na(res$eos_doy), "NA", sprintf("%.0f", res$eos_doy)),
      ifelse(is.na(res$sof_doy), "NA", sprintf("%.0f", res$sof_doy)),
      ifelse(is.finite(res$sos_annual_sd), res$sos_annual_sd, NA_real_),
      ifelse(nzchar(res$flags), paste0("[", res$flags, "]"), ""))
}

mean_out <- bind_rows(mean_rows) %>%
  select(site_key, veg_type, veg_source, pg_strt_thrshld, pg_end_thrshld,
         n_years, years, n_obs_mean_curve, sos_doy, sos_break_index,
         n_breaks, n_breaks_trimmed,
         sos_annual_mean, sos_annual_sd, n_years_sos, sos_resolution_days,
         baseline_evi, baseline_evi_at_sos, baseline_n_obs,
         peak_doy, peak_evi, n_peaks, second_peak_doy, amplitude, baseline_gap,
         thr_eos_evi, eos_doy, thr_sof_evi, sof_doy, window_days,
         pen_value, min_seg_obs,
         n_pts_dropped_year, n_pts_dropped_range, n_pts_dropped_robust,
         min_yrs_per_doy, n_doy_dropped_thin, yrs_at_baseline_doy,
         prior_doy, peak_doy_globalmax, prior_switched, n_peak_cands,
         flags) %>%
  arrange(site_key)
annual_out <- bind_rows(annual_rows) %>% arrange(site_key, year)

write_csv(mean_out,   file.path(OUT_DIR, "pelt_transitions_mean_modis.csv"))
write_csv(annual_out, file.path(OUT_DIR, "pelt_transitions_annual_modis.csv"))
prior_audit <- bind_rows(prior_audit_rows) %>%
  arrange(site_key, desc(arm == "mean"), year)
write_csv(prior_audit, file.path(OUT_DIR, "pelt_seasonal_prior_audit.csv"))

if (opt_dump) {
  curves_out <- bind_rows(curve_rows) %>% arrange(site_key, doy)
  write_csv(curves_out, file.path(OUT_DIR, "pelt_mean_curves_modis.csv"))
}

msg("")
msg("=== INPUT / PROCESSING SUMMARY ===")
msg("Sites processed        : %d", nrow(mean_out))
msg("  SOS resolved         : %d", sum(is.finite(mean_out$sos_doy)))
msg("  EOS resolved         : %d", sum(is.finite(mean_out$eos_doy)))
msg("  SOF resolved         : %d", sum(is.finite(mean_out$sof_doy)))
msg("  flagged              : %d", sum(nzchar(mean_out$flags)))
if (any(nzchar(mean_out$flags))) {
  bad <- mean_out %>% filter(nzchar(flags)) %>% select(site_key, flags)
  for (i in seq_len(nrow(bad))) msg("    %-10s %s", bad$site_key[i], bad$flags[i])
}

# ---- QA cross-check: our masked mean curve vs. the export's own precomputed
# multiyear_mean_evi at the sites flagged in the audit as carrying winter-edge
# extrapolation artefacts.  Informational only -- not built into the pipeline.
msg("")
msg("=== QA: robust-masked mean curve vs. precomputed multiyear_mean_evi (winter-edge artefact sites) ===")
artefact_sites <- c("BARR", "TOOL", "HEAL", "DEJU", "YELL")
for (sk in intersect(artefact_sites, mean_out$site_key)) {
  raw_min  <- ts_all %>% filter(site == sk) %>% pull(multiyear_mean_evi) %>% min(na.rm = TRUE)
  our_curve <- ts_all %>% filter(site == sk) %>% load_site_curves_modis() %>%
    group_by(doy) %>% summarise(fitted = mean(fitted, na.rm = TRUE), .groups = "drop")
  msg("  %-6s raw multiyear_mean_evi min %+.3f  |  our masked mean-curve min %+.3f",
      sk, raw_min, min(our_curve$fitted, na.rm = TRUE))
}

# ============================================================================
# VALIDATION
# Mirrors HLS/cmp/accept2.R (headline acceptance stats) and HLS/cmp/ccc.R
# (Lin's CCC), extended with a second arm because MODIS -- unlike HLS -- has
# a genuine per-site-year manual reference.  Sign convention throughout:
#   diff = manual - PELT   (same convention as EOS_SOF_comparison.csv)
# ============================================================================

ccc_fun <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 2) return(NA_real_)
  2 * cov(x, y) * (n - 1) / n /
    (var(x) * (n - 1) / n + var(y) * (n - 1) / n + (mean(x) - mean(y))^2)
}

score_block <- function(dE, dS, label) {
  msg("--- %s ---", label)
  for (nm in c("dE", "dS")) {
    v   <- if (nm == "dE") dE else dS
    lab <- if (nm == "dE") "EOS" else "SOF"
    v   <- v[is.finite(v)]
    n   <- length(v)
    if (n == 0) { msg("%s: no resolved pairs", lab); next }
    msg("%s: n=%-3d bias %+.2f  median|e| %.2f  mean|e| %.2f  w/in8 %d (%.0f%%)  w/in3 %d (%.0f%%)",
        lab, n, mean(v), median(abs(v)), mean(abs(v)),
        sum(abs(v) <= 8), 100 * sum(abs(v) <= 8) / n,
        sum(abs(v) <= 3), 100 * sum(abs(v) <= 3) / n)
  }
}

# ---------------------------------------------------------------- ARM 1: mean
msg("")
msg("############################################################")
msg("# HEADLINE ARM -- site-level mean curve vs modis_mean_transitions.csv")
msg("############################################################")

man_mean <- read_csv(MEAN_REF, show_col_types = FALSE, progress = FALSE) %>%
  mutate(site_key = normalize_key(site_id),
         man_eos  = as.integer(format(as.Date(mean_evi_start), "%j")),
         man_sof  = as.integer(format(as.Date(mean_evi_end),   "%j"))) %>%
  select(site_key, man_num_years = num_years, man_eos, man_sof)
# mean_evi_start / mean_evi_end are the EOS/SOF analogues: they mark the
# start and end of the peak-greenness plateau, exactly like HLS's own
# EOS/SOF columns in cmp/EOS_SOF_comparison.csv (col 2 = EOS = plateau
# start, col 3 = SOF = plateau end) -- same naming convention, just spelled
# "mean_evi_start/end" instead of "EOS/SOF" in this export.  veg_type is
# dropped here (kept only from the PELT side below) since both frames carry
# a column of that name and we already trust the export's own veg_type.

pelt_mean_keyed <- mean_out %>% mutate(site_key_norm = normalize_key(site_key))

matched_mean  <- intersect(man_mean$site_key, pelt_mean_keyed$site_key_norm)
dropped_man   <- setdiff(man_mean$site_key, pelt_mean_keyed$site_key_norm)
dropped_pelt  <- setdiff(pelt_mean_keyed$site_key_norm, man_mean$site_key)
msg("JOIN [mean arm]: manual rows %d, PELT sites %d -> matched %d",
    nrow(man_mean), nrow(pelt_mean_keyed), length(matched_mean))
msg("  in manual, NOT in PELT (no MODIS timeseries in export): %s",
    if (length(dropped_man)) paste(dropped_man, collapse = ", ") else "none")
msg("  in PELT, NOT in manual reference: %s",
    if (length(dropped_pelt)) paste(dropped_pelt, collapse = ", ") else "none")
if ("GUAN" %in% dropped_man || "LAJA" %in% dropped_man) {
  msg("  (GUAN/LAJA are expected here -- they carry num_years=0 in the")
  msg("   reference and have no rows at all in modis_evi_timeseries.csv)")
}

a_mean <- pelt_mean_keyed %>%
  inner_join(man_mean, by = c("site_key_norm" = "site_key")) %>%
  mutate(dE = man_eos - eos_doy, dS = man_sof - sof_doy)

msg("")
msg("--- ALL matched sites (n=%d) ---", nrow(a_mean))
score_block(a_mean$dE, a_mean$dS, "ALL sites")
msg("  CCC EOS: %.3f   CCC SOF: %.3f",
    ccc_fun(a_mean$eos_doy, a_mean$man_eos), ccc_fun(a_mean$sof_doy, a_mean$man_sof))

b_mean <- a_mean %>% filter(!site_key %in% KNOWN_OUTLIER_SITES)
msg("")
msg("--- excluding known outlier sites (%s), n=%d ---",
    paste(intersect(KNOWN_OUTLIER_SITES, a_mean$site_key), collapse = ", "), nrow(b_mean))
score_block(b_mean$dE, b_mean$dS, "ex-outliers")
msg("  CCC EOS: %.3f   CCC SOF: %.3f",
    ccc_fun(b_mean$eos_doy, b_mean$man_eos), ccc_fun(b_mean$sof_doy, b_mean$man_sof))

msg("")
msg("span mismatch note: num_years in modis_mean_transitions.csv reflects a")
msg("2015-2024 window (10 yr), while the mean-curve PELT run above pools the")
msg("full 2003-2024 record (up to 22 yr).  Not corrected here, flagged only.")
msg("num_years distribution in the reference:")
print(table(man_mean$man_num_years))

wr <- a_mean %>% mutate(tot = pmax(abs(dE), abs(dS), na.rm = TRUE)) %>%
  arrange(desc(tot)) %>% select(site_key, veg_type, dE, dS) %>% head(8)
msg("\n--- worst residuals (mean arm) ---")
print(as.data.frame(wr), digits = 3, row.names = FALSE)

# -------------------------------------------------------------- ARM 2: annual
msg("")
msg("############################################################")
msg("# SECOND ARM -- per site-year vs modis_annual_transitions.csv")
msg("############################################################")

man_ann <- read_csv(ANNUAL_REF, show_col_types = FALSE, progress = FALSE) %>%
  mutate(site_key = normalize_key(site_id))
# str_day / end_day are already DOY integers in this reference -- no date
# conversion needed (verified: str_day matches yday(evi_start) in 439/440
# rows; the two are redundant so str_day/end_day are used directly).

# Value-level plausibility guard on the REFERENCE itself (mirrors the
# EVI_MIN/MAX_PLAUSIBLE masking philosophy used on the input curves above):
# a DOY must be in [1,366]. One row (NOGP 2018) carries str_day = 1567,
# traced to a corrupted evi_start date ("2022-04-16") in the source export --
# a reference data-entry defect, not a PELT failure. Masked at the value
# level only (str_day -> NA) so the row's otherwise-valid end_day/SOF is
# not thrown away.
n_bad_str <- sum(man_ann$str_day < 1 | man_ann$str_day > 366, na.rm = TRUE)
n_bad_end <- sum(man_ann$end_day < 1 | man_ann$end_day > 366, na.rm = TRUE)
if (n_bad_str > 0 || n_bad_end > 0) {
  bad_rows <- man_ann %>% filter(str_day < 1 | str_day > 366 | end_day < 1 | end_day > 366) %>%
    select(site_id, survey_year, str_day, end_day)
  msg("NOTE: %d reference str_day / %d reference end_day value(s) outside",
      n_bad_str, n_bad_end)
  msg("plausible DOY range [1,366] -- masked to NA (reference data defect, not scored):")
  print(as.data.frame(bad_rows), row.names = FALSE)
}
man_ann <- man_ann %>% mutate(str_day = ifelse(str_day < 1 | str_day > 366, NA, str_day),
                               end_day = ifelse(end_day < 1 | end_day > 366, NA, end_day))

pelt_ann_keyed <- annual_out %>% mutate(site_key_norm = normalize_key(site_key))

matched_ann  <- pelt_ann_keyed %>%
  inner_join(man_ann, by = c("site_key_norm" = "site_key", "year" = "survey_year"))
msg("JOIN [annual arm]: manual site-years %d, PELT site-years %d -> matched %d",
    nrow(man_ann), nrow(pelt_ann_keyed), nrow(matched_ann))
dropped_ann_man <- setdiff(paste(man_ann$site_key, man_ann$survey_year),
                           paste(pelt_ann_keyed$site_key_norm, pelt_ann_keyed$year))
msg("  manual site-years with no PELT counterpart: %d%s",
    length(dropped_ann_man),
    if (length(dropped_ann_man)) paste0(" (e.g. ", paste(head(dropped_ann_man, 5), collapse = "; "), ")") else "")

a_ann <- matched_ann %>% mutate(dE = str_day - eos_doy, dS = end_day - sof_doy)
n_pelt_unresolved_eos <- sum(is.na(a_ann$eos_doy))
n_pelt_unresolved_sof <- sum(is.na(a_ann$sof_doy))
msg("  of matched site-years, PELT EOS unresolved: %d ; PELT SOF unresolved: %d",
    n_pelt_unresolved_eos, n_pelt_unresolved_sof)

msg("")
msg("--- ALL matched site-years (n=%d) ---", nrow(a_ann))
score_block(a_ann$dE, a_ann$dS, "ALL site-years")
msg("  CCC EOS: %.3f   CCC SOF: %.3f",
    ccc_fun(a_ann$eos_doy, a_ann$str_day), ccc_fun(a_ann$sof_doy, a_ann$end_day))

b_ann <- a_ann %>% filter(!site_key %in% KNOWN_OUTLIER_SITES)
msg("")
msg("--- excluding known outlier sites (%s), n=%d ---",
    paste(intersect(KNOWN_OUTLIER_SITES, a_ann$site_key), collapse = ", "), nrow(b_ann))
score_block(b_ann$dE, b_ann$dS, "ex-outliers")
msg("  CCC EOS: %.3f   CCC SOF: %.3f",
    ccc_fun(b_ann$eos_doy, b_ann$str_day), ccc_fun(b_ann$sof_doy, b_ann$end_day))

# ---------------------------------------------------------- per-site SD check
msg("")
msg("=== per-site interannual SD: manual vs PELT (2015-2024 overlap, one line) ===")
sd_tbl <- a_ann %>% group_by(site_key) %>%
  summarise(sd_man_eos = sd(str_day, na.rm = TRUE), sd_pelt_eos = sd(eos_doy, na.rm = TRUE),
            sd_man_sof = sd(end_day, na.rm = TRUE), sd_pelt_sof = sd(sof_doy, na.rm = TRUE),
            .groups = "drop")
msg("median per-site SD -- EOS: manual %.1f d vs PELT %.1f d | SOF: manual %.1f d vs PELT %.1f d",
    median(sd_tbl$sd_man_eos, na.rm = TRUE), median(sd_tbl$sd_pelt_eos, na.rm = TRUE),
    median(sd_tbl$sd_man_sof, na.rm = TRUE), median(sd_tbl$sd_pelt_sof, na.rm = TRUE))

# ---------------------------------------------------------------- write joins
cmp_mean <- a_mean %>%
  transmute(arm = "mean_curve", site_key, year = NA_integer_, veg_type,
            man_eos, pelt_eos = eos_doy, dE,
            man_sof, pelt_sof = sof_doy, dS,
            note = flags)
cmp_ann <- a_ann %>%
  transmute(arm = "annual", site_key, year, veg_type,
            man_eos = str_day, pelt_eos = eos_doy, dE,
            man_sof = end_day, pelt_sof = sof_doy, dS,
            note = flags)
comparison_out <- bind_rows(cmp_mean, cmp_ann) %>%
  mutate(abs_dE = abs(dE), abs_dS = abs(dS),
         is_known_outlier = site_key %in% KNOWN_OUTLIER_SITES)
write_csv(comparison_out, file.path(OUT_DIR, "pelt_modis_comparison.csv"))

msg("")
msg("Outputs written to %s", OUT_DIR)
msg("  pelt_transitions_mean_modis.csv    (%d rows)", nrow(mean_out))
msg("  pelt_transitions_annual_modis.csv  (%d rows)", nrow(annual_out))
msg("  pelt_modis_comparison.csv          (%d rows: %d mean-arm + %d annual-arm)",
    nrow(comparison_out), nrow(cmp_mean), nrow(cmp_ann))
if (opt_dump) {
  msg("  pelt_mean_curves_modis.csv         (%d rows, %d sites)",
      nrow(curves_out), length(unique(curves_out$site_key)))
}

# ============================================================================
# EXPERIMENT: does matching the mean-curve pooling window to the manual
# reference's window (2015-2024, MEAN_CURVE_YEARS) change headline-arm
# agreement?  Requested follow-up only -- does NOT touch the annual arm or
# the default full-record (2003-2024) outputs above, which are unchanged.
# Rebuilds the mean curve per site from `curves` restricted to
# MEAN_CURVE_YEARS, re-runs derive_transitions(), and re-scores against the
# same modis_mean_transitions.csv reference already loaded as man_mean.
# ============================================================================
msg("")
msg("############################################################")
msg("# EXPERIMENT -- mean curve restricted to %d-%d (matches manual window)",
    min(MEAN_CURVE_YEARS), max(MEAN_CURVE_YEARS))
msg("############################################################")

mean_rows_win <- list()
n_years_thin  <- list()
for (site_key in sites) {
  site_df <- ts_all %>% filter(site == site_key)
  curves  <- load_site_curves_modis(site_df)
  if (nrow(curves) == 0) next

  veg_type   <- site_df$veg_type[1]
  veg_source <- "timeseries"
  if (site_key %in% names(VEG_TYPE_OVERRIDES)) {
    veg_type   <- VEG_TYPE_OVERRIDES[[site_key]]
    veg_source <- "override"
  }
  sof_frac  <- if (!is.na(veg_type) && veg_type %in% names(SOF_FRAC_BY_VEG)) {
    SOF_FRAC_BY_VEG[[veg_type]]
  } else 0.90
  pen_value <- if (!is.na(veg_type) && veg_type == "DB") PEN_DB else PEN_DEFAULT

  curves_win <- curves %>% filter(year %in% MEAN_CURVE_YEARS)
  n_yrs_avail <- length(unique(curves_win$year))
  n_years_thin[[site_key]] <- n_yrs_avail

  # Rebuild the prior from this window's peaks (matching the mean curve's span)
  pass1_peaks_win <- vapply(MEAN_CURVE_YEARS, function(yy) {
    cy <- curves_win %>% filter(year == yy) %>% arrange(doy)
    if (nrow(cy) < MIN_YEAR_OBS) return(NA_real_)
    pk <- find_peaks(cy$doy, cy$fitted)
    if (is.na(pk$primary_idx)) return(NA_real_)
    cy$doy[pk$primary_idx]
  }, numeric(1))
  prior_mean_win <- if (SEASONAL_PRIOR_MEAN) circ_median(pass1_peaks_win) else NA_real_

  mc <- curves_win %>% group_by(doy) %>%
    summarise(fitted = mean(fitted, na.rm = TRUE), n_yr = n(), .groups = "drop") %>%
    arrange(doy) %>% filter(n_yr >= MIN_YRS_PER_DOY)

  res <- derive_transitions(mc$doy, mc$fitted, sof_frac, pen_value,
                            prior_doy = prior_mean_win)
  res$site_key   <- site_key
  res$veg_type   <- veg_type
  res$veg_source <- veg_source
  res$n_years    <- n_yrs_avail
  res$years      <- if (n_yrs_avail) paste(range(curves_win$year), collapse = "-") else NA_character_
  mean_rows_win[[site_key]] <- as.data.frame(res, stringsAsFactors = FALSE)
}

mean_out_win <- bind_rows(mean_rows_win) %>%
  select(site_key, veg_type, veg_source, n_years, years,
         sos_doy, baseline_evi, peak_doy, peak_evi, amplitude,
         eos_doy, sof_doy, window_days, flags) %>%
  arrange(site_key)
write_csv(mean_out_win, file.path(OUT_DIR, "pelt_transitions_mean_modis_2015_2024.csv"))

thin_tbl <- tibble(site_key = names(n_years_thin), n_years_2015_2024 = unlist(n_years_thin))
n_thin <- sum(thin_tbl$n_years_2015_2024 < length(MEAN_CURVE_YEARS))
msg("")
msg("Q4: sites with fewer than %d usable years within %d-%d after masking: %d of %d",
    length(MEAN_CURVE_YEARS), min(MEAN_CURVE_YEARS), max(MEAN_CURVE_YEARS), n_thin, nrow(thin_tbl))
if (n_thin > 0) {
  print(as.data.frame(thin_tbl %>% filter(n_years_2015_2024 < length(MEAN_CURVE_YEARS))), row.names = FALSE)
} else {
  msg("  (all %d sites retain the full %d years -- the restriction does not thin the curve anywhere)",
      nrow(thin_tbl), length(MEAN_CURVE_YEARS))
}

pelt_mean_win_keyed <- mean_out_win %>% mutate(site_key_norm = normalize_key(site_key))
a_mean_win <- pelt_mean_win_keyed %>%
  inner_join(man_mean, by = c("site_key_norm" = "site_key")) %>%
  mutate(dE = man_eos - eos_doy, dS = man_sof - sof_doy)

msg("")
msg("--- matched-window (%d-%d) ALL sites (n=%d) ---",
    min(MEAN_CURVE_YEARS), max(MEAN_CURVE_YEARS), nrow(a_mean_win))
score_block(a_mean_win$dE, a_mean_win$dS, "matched-window ALL sites")
ccc_eos_win_all <- ccc_fun(a_mean_win$eos_doy, a_mean_win$man_eos)
ccc_sof_win_all <- ccc_fun(a_mean_win$sof_doy, a_mean_win$man_sof)
msg("  CCC EOS: %.3f   CCC SOF: %.3f", ccc_eos_win_all, ccc_sof_win_all)

b_mean_win <- a_mean_win %>% filter(!site_key %in% KNOWN_OUTLIER_SITES)
msg("")
msg("--- matched-window ex-outliers (%s), n=%d ---",
    paste(intersect(KNOWN_OUTLIER_SITES, a_mean_win$site_key), collapse = ", "), nrow(b_mean_win))
score_block(b_mean_win$dE, b_mean_win$dS, "matched-window ex-outliers")
ccc_eos_win_ex <- ccc_fun(b_mean_win$eos_doy, b_mean_win$man_eos)
ccc_sof_win_ex <- ccc_fun(b_mean_win$sof_doy, b_mean_win$man_sof)
msg("  CCC EOS: %.3f   CCC SOF: %.3f", ccc_eos_win_ex, ccc_sof_win_ex)

# ---- side-by-side summary table (both variants, both outlier treatments) --
qstat <- function(dE, dS, ce, cs) {
  v_e <- dE[is.finite(dE)]; v_s <- dS[is.finite(dS)]
  list(EOS = c(bias = mean(v_e), med = median(abs(v_e)), mean = mean(abs(v_e)),
               w8 = 100*sum(abs(v_e)<=8)/length(v_e), w3 = 100*sum(abs(v_e)<=3)/length(v_e), ccc = ce),
       SOF = c(bias = mean(v_s), med = median(abs(v_s)), mean = mean(abs(v_s)),
               w8 = 100*sum(abs(v_s)<=8)/length(v_s), w3 = 100*sum(abs(v_s)<=3)/length(v_s), ccc = cs))
}
full_all <- qstat(a_mean$dE, a_mean$dS, ccc_fun(a_mean$eos_doy, a_mean$man_eos), ccc_fun(a_mean$sof_doy, a_mean$man_sof))
full_ex  <- qstat(b_mean$dE, b_mean$dS, ccc_fun(b_mean$eos_doy, b_mean$man_eos), ccc_fun(b_mean$sof_doy, b_mean$man_sof))
win_all  <- qstat(a_mean_win$dE, a_mean_win$dS, ccc_eos_win_all, ccc_sof_win_all)
win_ex   <- qstat(b_mean_win$dE, b_mean_win$dS, ccc_eos_win_ex,  ccc_sof_win_ex)

fmt_row <- function(variant, outl, var, s) {
  sprintf("%-28s %-11s %-4s %+6.2f %6.2f %6.2f %5.0f%% %5.0f%% %6.3f",
          variant, outl, var, s["bias"], s["med"], s["mean"], s["w8"], s["w3"], s["ccc"])
}
msg("")
msg("=== SIDE-BY-SIDE: mean-curve window comparison ===")
msg("%-28s %-11s %-4s %6s %6s %6s %6s %6s %6s",
    "variant", "outliers", "var", "bias", "med|e|", "mean|e|", "<=8d", "<=3d", "CCC")
msg("%s", fmt_row("mean curve 2003-2024", "all sites",  "EOS", full_all$EOS))
msg("%s", fmt_row("mean curve 2003-2024", "all sites",  "SOF", full_all$SOF))
msg("%s", fmt_row("mean curve 2015-2024", "all sites",  "EOS", win_all$EOS))
msg("%s", fmt_row("mean curve 2015-2024", "all sites",  "SOF", win_all$SOF))
msg("%s", fmt_row("mean curve 2003-2024", "ex-outliers","EOS", full_ex$EOS))
msg("%s", fmt_row("mean curve 2003-2024", "ex-outliers","SOF", full_ex$SOF))
msg("%s", fmt_row("mean curve 2015-2024", "ex-outliers","EOS", win_ex$EOS))
msg("%s", fmt_row("mean curve 2015-2024", "ex-outliers","SOF", win_ex$SOF))

# ---- Q3: MOAB's SOF error specifically, both variants ----------------------
msg("")
msg("Q3: MOAB SOF diff (man - PELT) -- full 2003-2024: %s | matched 2015-2024: %s",
    { v <- a_mean$dS[a_mean$site_key == "MOAB"]; if (length(v)) sprintf("%+.2f d", v) else "NA" },
    { v <- a_mean_win$dS[a_mean_win$site_key == "MOAB"]; if (length(v)) sprintf("%+.2f d", v) else "NA" })
moab_full_sof_doy <- mean_out$sof_doy[mean_out$site_key == "MOAB"]
moab_win_sof_doy  <- mean_out_win$sof_doy[mean_out_win$site_key == "MOAB"]
msg("    MOAB PELT sof_doy -- full: %.1f | matched: %.1f  (peak_doy full: %.0f, matched: %.0f)",
    moab_full_sof_doy, moab_win_sof_doy,
    mean_out$peak_doy[mean_out$site_key == "MOAB"], mean_out_win$peak_doy[mean_out_win$site_key == "MOAB"])

# ---- joined diffs for this variant, written alongside the original --------
cmp_mean_win <- a_mean_win %>%
  transmute(arm = "mean_curve_2015_2024", site_key, veg_type,
            man_eos, pelt_eos = eos_doy, dE,
            man_sof, pelt_sof = sof_doy, dS,
            note = flags, abs_dE = abs(dE), abs_dS = abs(dS),
            is_known_outlier = site_key %in% KNOWN_OUTLIER_SITES)
write_csv(cmp_mean_win, file.path(OUT_DIR, "pelt_modis_comparison_2015_2024.csv"))

msg("")
msg("Additional outputs written (this experiment only, alongside originals):")
msg("  pelt_transitions_mean_modis_2015_2024.csv  (%d rows)", nrow(mean_out_win))
msg("  pelt_modis_comparison_2015_2024.csv        (%d rows)", nrow(cmp_mean_win))
