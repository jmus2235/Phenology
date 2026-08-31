nk <- function(x) gsub("[-/]","_",trimws(x))
cv <- read.csv("C:/Users/jmusinsky/Documents/R_Scripts/MODIS/pelt_out/pelt_mean_curves_modis.csv", stringsAsFactors=FALSE)
cv$k <- nk(cv$site_key)
m <- read.csv("C:/Users/jmusinsky/Documents/R_Scripts/MODIS/cmp/manual_sos_TEMPLATE.csv", stringsAsFactors=FALSE)
m <- m[!is.na(m$sos_doy) & nzchar(trimws(m$site)),]; m$k <- nk(m$site)
p <- read.csv("C:/Users/jmusinsky/Documents/R_Scripts/MODIS/pelt_out/pelt_transitions_mean_modis.csv", stringsAsFactors=FALSE); p$k <- nk(p$site_key)
STRUCT <- c("SJER","SRER","BARR","TOOL","RMNP","BONA","DEJU")
res <- do.call(rbind, lapply(split(cv, cv$k), function(g){
  k <- g$k[1]; pk <- p$peak_doy[p$k==k]; ms <- m$sos_doy[m$k==k]
  pre <- g[g$doy < pk,]; if(!nrow(pre)) return(NULL)
  vd <- pre$doy[which.min(pre$fitted)]
  data.frame(k, valley_doy=vd, manual_sos=ms, sos_minus_valley=ms-vd,
             pelt_sos=p$sos_doy[p$k==k], stringsAsFactors=FALSE)}))
res$struct <- res$k %in% STRUCT
cat("=== sites where manual SOS lands BEFORE the pre-peak minimum (the TEAK/OSBS defect) ===\n")
bad <- res[res$sos_minus_valley < 0,]
if(!nrow(bad)) cat("  none\n") else print(bad, row.names=FALSE)
cat("\n=== TEST SENSITIVITY: it can only fire where the valley bottom is LATE ===\n")
i <- res[!res$struct,]
cat(sprintf("interpretable n=%d; valley_doy median=%d\n", nrow(i), median(i$valley_doy)))
cat(sprintf("valley at DOY<=17 (test blind, manual SOS cannot precede it): %d of %d\n", sum(i$valley_doy<=17), nrow(i)))
cat(sprintf("valley at DOY>=33 (test has power here):                      %d of %d\n", sum(i$valley_doy>=33), nrow(i)))
cat("\nsites where the test HAS power, sorted by margin:\n")
pw <- i[i$valley_doy>=33,c("k","valley_doy","manual_sos","sos_minus_valley","pelt_sos")]
print(pw[order(pw$sos_minus_valley),], row.names=FALSE)
