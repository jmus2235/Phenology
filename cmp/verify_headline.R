d <- read.csv("C:/Users/jmusinsky/Documents/R_Scripts/MODIS/pelt_out/pelt_modis_comparison.csv",
              stringsAsFactors = FALSE)
ccc <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  2 * cov(x, y) * (length(x) - 1) / length(x) /
    (var(x) * (length(x) - 1) / length(x) + var(y) * (length(y) - 1) / length(y) +
       (mean(x) - mean(y))^2)
}
for (a in unique(d$arm)) {
  s <- d[d$arm == a, ]
  cat(sprintf("--- arm=%s  n=%d\n", a, nrow(s)))
  for (v in c("E", "S")) {
    dd <- s[[paste0("d", v)]]
    mm <- s[[paste0("man_", if (v == "E") "eos" else "sof")]]
    pp <- s[[paste0("pelt_", if (v == "E") "eos" else "sof")]]
    ok <- is.finite(dd)
    cat(sprintf("  %s  n=%3d bias=%+6.2f  med|e|=%5.2f  <=8d=%5.1f%%  <=3d=%5.1f%%  CCC=%.4f\n",
                if (v == "E") "EOS" else "SOF", sum(ok), mean(dd[ok]), median(abs(dd[ok])),
                100 * mean(abs(dd[ok]) <= 8), 100 * mean(abs(dd[ok]) <= 3), ccc(mm, pp)))
  }
}
