nk <- function(x) gsub("[-/]", "_", trimws(x))
m <- read.csv("C:/Users/jmusinsky/Documents/R_Scripts/MODIS/cmp/manual_sos_TEMPLATE.csv", stringsAsFactors=FALSE)
m <- m[!is.na(m$sos_doy) & nzchar(trimws(m$site)), ]; m$k <- nk(m$site)
p <- read.csv("C:/Users/jmusinsky/Documents/R_Scripts/MODIS/pelt_out/pelt_transitions_mean_modis.csv", stringsAsFactors=FALSE); p$k <- nk(p$site_key)
cp <- read.csv("C:/Users/jmusinsky/Documents/R_Scripts/MODIS/pelt_out/pelt_modis_comparison.csv", stringsAsFactors=FALSE)
cp <- cp[cp$arm=="mean_curve",]; cp$k <- nk(cp$site_key)
j <- merge(merge(m,p,by="k"), cp[,c("k","dE","dS")], by="k")
j$dSOS <- j$sos_doy.x - j$sos_doy.y
j$dWin <- j$dS - j$dE                      # manual window - PELT window

cat("=== Is PELT's SOS later, yet its BASELINE lower? ===\n")
cat(sprintf("PELT SOS later than manual at %d of %d sites\n", sum(j$dSOS<0), nrow(j)))
cat(sprintf("PELT segment-median baseline BELOW manual point-baseline at %d of %d\n", sum(j$baseline_evi.y < j$baseline_evi.x), nrow(j)))
cat(sprintf("PELT curve value AT ITS OWN SOS above manual baseline at %d of %d\n\n", sum(j$baseline_evi_at_sos > j$baseline_evi.x), nrow(j)))

cat("=== Worked example, 4 sites ===\n")
ex <- j[j$k %in% c("SCBI","UKFS","BLAN","NOGP"),]
for (i in seq_len(nrow(ex))) with(ex[i,], {
  f <- pg_end_thrshld/100
  te_m <- baseline_evi.x + 0.90*(peak_evi-baseline_evi.x); ts_m <- baseline_evi.x + f*(peak_evi-baseline_evi.x)
  te_p <- baseline_evi.y + 0.90*(peak_evi-baseline_evi.y); ts_p <- baseline_evi.y + f*(peak_evi-baseline_evi.y)
  cat(sprintf("%-6s peak %.4f | SOS man %3d / PELT %3d (PELT LATER by %d d)\n", k, peak_evi, sos_doy.x, sos_doy.y, -dSOS))
  cat(sprintf("        baseline: manual(point@SOS) %.4f  vs  PELT(median of pre-break segment) %.4f  -> PELT LOWER by %.4f\n",
              baseline_evi.x, baseline_evi.y, baseline_evi.x-baseline_evi.y))
  cat(sprintf("        [PELT's curve value at its OWN SOS would have been %.4f - HIGHER, but it is not used]\n", baseline_evi_at_sos))
  cat(sprintf("        thr_EOS: man %.4f > PELT %.4f -> PELT crosses rising limb EARLIER (dE=%+.1f d)\n", te_m, te_p, dE))
  cat(sprintf("        thr_SOF: man %.4f > PELT %.4f -> PELT crosses falling limb LATER  (dS=%+.1f d)\n", ts_m, ts_p, dS))
  cat(sprintf("        => PELT window WIDER by %.1f d\n\n", -dWin))
})
cat(sprintf("=== Network: PELT window wider at %d of %d sites; mean %+.2f d ===\n",
            sum(j$dWin<0), nrow(j), -mean(j$dWin)))
