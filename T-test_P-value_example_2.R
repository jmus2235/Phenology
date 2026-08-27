library(dplyr)

# Example data
year <- c(2003	,2004	,2005	,2006	,2007	,2008	,2009	,2010	,2011 ,2012	,2013	,2014	,2015	,2016	,2017	,2018	,2019 ,2020	,2021	,2022, 2023)
# y <- c(0.333559701, 0.319608102, 0.325745957, 0.340294167, 0.325607349, 0.337607904, 0.323931862, 0.331472699, 0.33307083, 0.348654962, 0.324740774, 0.321847658, 0.324415468, 0.298920836, 0.301237884, 0.326103498, 0.318569884, 0.315714466, 0.294759569, 0.319702449)
DOY <- c(196, 221, 210, 208, 192, 213, 222, 220, 216, 214, 213, 214, 198, 201, 205, 210, 211, 211, 224, 211, 206)

# Create a merged dataframe
meanEVI_DOY <- data.frame(year,DOY)

# Perform linear regression
model <- lm(meanEVI_DOY$DOY ~ meanEVI_DOY$year)

# Extract the coefficient for the slope
slope <- coef(model)[[2]]

# Perform t-test
t_test <- summary(model)$coefficients["meanEVI_DOY$year", ]

# Extract the t-value and p-value
t_value <- t_test["t value"]
p_value <- t_test["Pr(>|t|)"]

# Compute the confidence interval
conf_interval <- confint(model)["meanEVI_DOY$year", ]

# Check if the slope is significantly different from zero
if (p_value < 0.05) {
  conclusion <- "The slope is significantly different from zero."
} else {
  conclusion <- "The slope is NOT significantly different from zero."
}

meanEVI_DOY_filt1 <- filter(meanEVI_DOY, year <2014)
meanEVI_DOY_filt1$GROUP <- 'A'

meanEVI_DOY_filt2 <- filter(meanEVI_DOY, year >2013)
meanEVI_DOY_filt2 <- filter(meanEVI_DOY_filt2, year <2024)
meanEVI_DOY_filt2$GROUP <- 'B'

meanEVI_DOY2 <- rbind(meanEVI_DOY_filt1, meanEVI_DOY_filt2)

# Compute two-sample t-test (assuming unequal variances) between two EVI time series (2003-2012 and 2013-2022)
Ttest <- t.test(DOY ~ GROUP, data = meanEVI_DOY2, paired = FALSE, var.equal= FALSE)
# Convert "p.value" field to dataframe
PTtest <- data.frame("Two_sample.P_value"=Ttest$p.value)

# Check if the P.value from Two-sample T-test is significant
if (PTtest$Two_sample.P_value < 0.05) {
  conclusion2 <- "The mean EOS from the earlier years is significantly different than that of the later years."
} else {
  conclusion2 <- "The mean EOS from the earlier years is NOT significantly different than the later years."
}

# Print the results
cat("Slope Test All Years", "\n")
cat("Slope:", slope, "\n")
cat("t-value:", t_value, "\n")
cat("p-value:", p_value, "\n")
cat("95% Confidence Interval:", conf_interval, "\n")
cat("Conclusion:", conclusion, "\n", "\n")

cat("Two-sample T-test", "\n")
cat("p-value:", PTtest$Two_sample.P_value, "\n")
cat("Conclusion:", conclusion2)

