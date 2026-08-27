# Example data
x <- c(2003	,2004	,2005	,2006	,2007	,2008	,2009	,2010	,2011	,2012	,2013	,2014	,2015	,2016	,2017	,2018	,2019	,2020	,2021	,2022)
# y <- c(0.333559701, 0.319608102, 0.325745957, 0.340294167, 0.325607349, 0.337607904, 0.323931862, 0.331472699, 0.33307083, 0.348654962, 0.324740774, 0.321847658, 0.324415468, 0.298920836, 0.301237884, 0.326103498, 0.318569884, 0.315714466, 0.294759569, 0.319702449)
y <- c(204, 216, 199, 191, 216, 209, 220, 201, 230, 195, 211, 218, 210, 228, 208, 186, 247, 231, 284, 211)
# Perform linear regression
model <- lm(y ~ x)

# Extract the coefficient for the slope
slope <- coef(model)[[2]]

# Perform t-test
t_test <- summary(model)$coefficients["x", ]

# Extract the t-value and p-value
t_value <- t_test["t value"]
p_value <- t_test["Pr(>|t|)"]

# Compute the confidence interval
conf_interval <- confint(model)["x", ]

# Check if the slope is significantly different from zero
if (p_value < 0.05) {
  conclusion <- "The slope is significantly different from zero."
} else {
  conclusion <- "The slope is not significantly different from zero."
}

# Print the results
cat("Slope:", slope, "\n")
cat("t-value:", t_value, "\n")
cat("p-value:", p_value, "\n")
cat("95% Confidence Interval:", conf_interval, "\n")
cat("Conclusion:", conclusion)


df_EVI_DOY_filt1 <- filter(df_EVI_DOY, YEAR <2013)
df_EVI_DOY_filt1$GROUP <- 'A'

df_EVI_DOY_filt2 <- filter(df_EVI_DOY, YEAR >2012)
df_EVI_DOY_filt2 <- filter(df_EVI_DOY_filt2, YEAR <2023)
df_EVI_DOY_filt2$GROUP <- 'B'

# Compute two-sample t-test (assuming unequal variances) between two EVI time series (2003-2012 and 2013-2022)
PTtest <- t.test(FITTED ~ GROUP, data = df_EVI_DOY, paired = FALSE, var.equal= FALSE)
# Convert "p.value" field to dataframe
PTtest2 <- data.frame("Two_sample.P_value"=PTtest$p.value)
#Add DOY_to_process to new column in data frame
PTtest2$DOY <- DOY_to_process