nk <- function(x) gsub("[-/]", "_", trimws(x))
m <- read.csv("C:/Users/jmusinsky/Documents/R_Scripts/MODIS/cmp/manual_sos_TEMPLATE.csv", stringsAsFactors=FALSE)
m <- m[!is.na(m$sos_doy) & nzchar(trimws(m$site)), ]; m$k <- nk(m$site)
p <- read.csv("C:/Users/jmusinsky/Documents/R_Scripts/MODIS/pelt_out/pelt_transitions_mean_modis.csv", stringsAsFactors=FALSE); p$k <- nk(p$site_key)
j <- merge(m,p,by="k"); j$d <- j$sos_doy.x - j$sos_doy.y
ccc <- function(x,y){n<-length(x); 2*cov(x,y)*(n-1)/n/(var(x)*(n-1)/n+var(y)*(n-1)/n+(mean(x)-mean(y))^2)}
STRUCT <- c("SJER","SRER","BARR","TOOL","RMNP","BONA","DEJU"); BIMOD <- c("MOAB","STER")
sets <- list(A=j, B=j[!j$k %in% STRUCT,], C=j[!j$k %in% c(STRUCT,BIMOD),])
for (nm in names(sets)) { s <- sets[[nm]]
  cat(sprintf("SET %s n=%2d bias=%+7.2f med|e|=%5.2f mean|e|=%5.2f <=8d=%4.0f%% <=16d=%4.0f%% CCC=%+.3f LoA=[%+.1f,%+.1f]\n",
    nm,nrow(s),mean(s$d),median(abs(s$d)),mean(abs(s$d)),100*mean(abs(s$d)<=8),100*mean(abs(s$d)<=16),
    ccc(s$sos_doy.x,s$sos_doy.y), mean(s$d)-1.96*sd(s$d), mean(s$d)+1.96*sd(s$d))) }
cat("\n-- ANOVA SOS error ~ veg (Set C) --\n")
s <- sets$C; print(summary(aov(d ~ factor(veg_type), data=s)))
cat("\n-- veg cell sizes (Set C) --\n"); print(table(s$veg_type))
cat("\n-- DYNAMIC RANGE, why cross-variable CCC comparison is unsafe --\n")
cp <- read.csv("C:/Users/jmusinsky/Documents/R_Scripts/MODIS/pelt_out/pelt_modis_comparison.csv", stringsAsFactors=FALSE)
cp <- cp[cp$arm=="mean_curve",]
cat(sprintf("manual SOS: sd=%5.1f d range=%3d-%3d (SetC n=%d)\n", sd(s$sos_doy.x), min(s$sos_doy.x), max(s$sos_doy.x), nrow(s)))
cat(sprintf("manual EOS: sd=%5.1f d range=%3.0f-%3.0f (n=%d)\n", sd(cp$man_eos), min(cp$man_eos), max(cp$man_eos), nrow(cp)))
cat(sprintf("manual SOF: sd=%5.1f d range=%3.0f-%3.0f (n=%d)\n", sd(cp$man_sof), min(cp$man_sof), max(cp$man_sof), nrow(cp)))
cat("\n-- how many dSOS are exact multiples of 8? --\n")
cat(sprintf("%d of %d (%.0f%%)\n", sum(j$d %% 8 == 0), nrow(j), 100*mean(j$d %% 8 == 0)))
