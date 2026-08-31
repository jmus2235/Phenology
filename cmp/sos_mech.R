nk <- function(x) gsub("[-/]", "_", trimws(x))
m <- read.csv("C:/Users/jmusinsky/Documents/R_Scripts/MODIS/cmp/manual_sos_TEMPLATE.csv", stringsAsFactors=FALSE)
m <- m[!is.na(m$sos_doy) & nzchar(trimws(m$site)), ]; m$k <- nk(m$site)
p <- read.csv("C:/Users/jmusinsky/Documents/R_Scripts/MODIS/pelt_out/pelt_transitions_mean_modis.csv", stringsAsFactors=FALSE); p$k <- nk(p$site_key)
cp <- read.csv("C:/Users/jmusinsky/Documents/R_Scripts/MODIS/pelt_out/pelt_modis_comparison.csv", stringsAsFactors=FALSE)
cp <- cp[cp$arm=="mean_curve", ]; cp$k <- nk(cp$site_key)

j <- merge(merge(m, p, by="k"), cp[,c("k","dE","dS")], by="k")
j$dSOS <- j$sos_doy.x - j$sos_doy.y
j$base_gap <- j$baseline_evi.x - j$baseline_evi.y   # manual(point) - PELT(segment median)

cat("=== CHECK 1: where dSOS==0, does manual baseline == PELT curve value at that DOY? ===\n")
z <- j[j$dSOS==0, c("k","sos_doy.x","baseline_evi.x","baseline_evi_at_sos")]
z$resid <- z$baseline_evi.x - z$baseline_evi_at_sos
print(z, row.names=FALSE, digits=12)

cat("\n=== CHECK 2: sign/direction of the baseline gap ===\n")
cat(sprintf("manual point-baseline ABOVE PELT segment-median at %d of %d sites; median gap %+.4f EVI\n",
            sum(j$base_gap>0), nrow(j), median(j$base_gap)))

cat("\n=== CHECK 3: predicted vs actual SOF shift from the baseline-rule difference ===\n")
# thr_sof = base + f*(peak-base)  ->  d(thr)/d(base) = (1-f)
j$f <- j$pg_end_thrshld/100
j$dthr_sof <- (1-j$f) * j$base_gap          # manual thr - PELT thr (EVI units)
j$dthr_eos <- (1-0.90)   * j$base_gap
ex <- !(j$k %in% c("SJER","SRER","MOAB","STER","BARR","TOOL","RMNP","BONA","DEJU"))
for (lab in c("all","ex_structural")) {
  s <- if (lab=="all") j else j[ex,]
  ct <- cor.test(s$dthr_sof, s$dS, method="pearson")
  cte <- cor.test(s$dthr_eos, s$dE, method="pearson")
  cat(sprintf("%-14s n=%2d  SOF: r=%+.3f p=%.4f | EOS: r=%+.3f p=%.4f\n",
              lab, nrow(s), ct$estimate, ct$p.value, cte$estimate, cte$p.value))
}
cat("\n=== CHECK 4: bias direction restated ===\n")
cat(sprintf("mean dE = %+.2f  -> PELT EOS is %s than manual\n", mean(j$dE), ifelse(mean(j$dE)>0,"EARLIER","LATER")))
cat(sprintf("mean dS = %+.2f  -> PELT SOF is %s than manual\n", mean(j$dS), ifelse(mean(j$dS)>0,"EARLIER","LATER")))
cat(sprintf("PELT window vs manual: %+.2f d (positive = PELT WIDER)\n", mean(j$dE) - mean(j$dS)))
