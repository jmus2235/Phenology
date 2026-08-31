# Would the ORIGINAL three manual values have been flagged by a rule that NEVER
# looks at PELT?  The valley rule uses only (a) the mean curve and (b) the manual
# SOS.  If yes, the revisions are not comparator-driven.
nk <- function(x) gsub("[-/]","_",trimws(x))
cv <- read.csv("C:/Users/jmusinsky/Documents/R_Scripts/MODIS/pelt_out/pelt_mean_curves_modis.csv", stringsAsFactors=FALSE); cv$k <- nk(cv$site_key)
p  <- read.csv("C:/Users/jmusinsky/Documents/R_Scripts/MODIS/pelt_out/pelt_transitions_mean_modis.csv", stringsAsFactors=FALSE); p$k <- nk(p$site_key)
orig <- data.frame(k=c("TEAK","OSBS","NIWO"), sos_orig=c(57,33,105), sos_new=c(81,49,121))
for (i in seq_len(nrow(orig))) { k <- orig$k[i]
  g <- cv[cv$k==k,]; pk <- p$peak_doy[p$k==k]; pre <- g[g$doy<pk,]
  vd <- pre$doy[which.min(pre$fitted)]
  cat(sprintf("%-5s valley_doy=%3d | ORIGINAL sos=%3d margin=%+3d %s | REVISED sos=%3d margin=%+3d | PELT sos=%3d\n",
      k, vd, orig$sos_orig[i], orig$sos_orig[i]-vd,
      ifelse(orig$sos_orig[i]-vd < 0, "<-FLAGGED (before valley)",
      ifelse(orig$sos_orig[i]-vd == 0, "<-FLAGGED (at valley)   ", "   ok                    ")),
      orig$sos_new[i], orig$sos_new[i]-vd, p$sos_doy[p$k==k])) }
