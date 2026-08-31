nk <- function(x) gsub("[-/]", "_", trimws(x))
m <- read.csv("C:/Users/jmusinsky/Documents/R_Scripts/MODIS/cmp/manual_sos_TEMPLATE.csv",
              stringsAsFactors = FALSE)
m <- m[!is.na(m$sos_doy) & nzchar(trimws(m$site)), ]
m$k <- nk(m$site)
p <- read.csv("C:/Users/jmusinsky/Documents/R_Scripts/MODIS/pelt_out/pelt_transitions_mean_modis.csv",
              stringsAsFactors = FALSE)
p$k <- nk(p$site_key)

cat("manual n=", nrow(m), "  pelt n=", nrow(p), "\n")
cat("in manual not pelt:", paste(setdiff(m$k, p$k), collapse=", "), "\n")
cat("in pelt not manual:", paste(setdiff(p$k, m$k), collapse=", "), "\n\n")

j <- merge(m, p, by = "k")
j$dSOS  <- j$sos_doy.x - j$sos_doy.y
j$dBase_seg <- j$baseline_evi.x - j$baseline_evi.y          # vs segment median
j$dBase_pt  <- j$baseline_evi.x - j$baseline_evi_at_sos      # vs point-at-PELT-SOS
j <- j[order(-abs(j$dSOS)), ]
cat(sprintf("%-10s %4s %4s %7s | %8s %8s %8s | %s\n",
            "site","man","pelt","dSOS","man_base","pelt_seg","pelt_@sos","veg"))
for (i in seq_len(nrow(j))) with(j[i,], cat(sprintf(
  "%-10s %4d %4d %+7d | %8.4f %8.4f %8.4f | %s\n",
  k, sos_doy.x, sos_doy.y, dSOS, baseline_evi.x, baseline_evi.y, baseline_evi_at_sos, veg_type)))
cat(sprintf("\nn=%d  dSOS: median %+0.1f  mean %+0.1f  median|e| %.1f  <=8d %.0f%%  <=16d %.0f%%\n",
    nrow(j), median(j$dSOS), mean(j$dSOS), median(abs(j$dSOS)),
    100*mean(abs(j$dSOS)<=8), 100*mean(abs(j$dSOS)<=16)))
