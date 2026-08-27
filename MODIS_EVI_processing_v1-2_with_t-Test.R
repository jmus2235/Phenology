## ----install-packages, message = FALSE-----------------------------------
# uncomment for install
# install.packages("lubridate")
# install.packages("changepoint")
# install.packages("readr")
# install.packages("ggplot2")
# install.packages("tidyverse")
# install.packages("dplyr")
# install.packages("ggpubr")
# install.packages("PairedData")
#install.packages("emuR")

##
## Script for processing smoothed EVI 8-day composite time series from high-quality (QA-bit=0) MODIS Terra and Aqua data 
## produced from Earth Engine script "Mean 2003-2019 MODIS VIIRS EVI per Site - JM v4B" 
## Input CSV table must contain three columns of data in following format: 
##
## system:index	date	      mean
## 2003_01_01   1.04138E+12	885.0715865
## 2003_01_17	  1.04276E+12	867.7747245
## 2003_02_02	  1.04414E+12	869.3154614
## ...

library(readr)
library(ggplot2)
library(changepoint)
library(lubridate)
library(tidyverse)
library(dplyr)
#library(emuR)
# library(ggpubr)
# library(PairedData)

# Set working directory
setwd("~/R_Scripts/MODIS/data")

#### USER-SPECIFIED PARAMETERS
EVIFileName = 'D01_HARV_TOS_DB'
QAbits = '01'
startyear <- 2006
sensor = 'TERRA_AQUA'
#sensor = 'TERRA'
#sensor = 'AQUA'
year_to_process <- '2003'
DOY_EVI = '1'

loessSPAN <- 0.25
interpolate <- 0 # 0 for 'No', 1 for 'Yes'
####

# Read csv file
df <- read.csv(file=sprintf("%s_MODIS_%s_EVI_QAbit_%s.csv",EVIFileName, sensor, QAbits), stringsAsFactors = FALSE)

# Re-scale EVI data
df$MEAN_EVI <- df$mean*0.0001

# Separate YEAR, MONTH and DAY and convert to numeric
df <- df %>% separate(system.index, c("YEAR", "MONTH", "DAY"), sep="_")
df$YEAR <- as.numeric(df$YEAR)
df$MONTH <- as.numeric(df$MONTH)
df$DAY <- as.numeric(df$DAY)

# Sort data frame by YEAR, MONTH, DAY
df <- df[order(df$YEAR, df$MONTH, df$DAY),]

# Remove old date column
df <- df[-4]

# Add a new reformatted DATE column
df$DATE <- with(df, ymd(sprintf('%04d%02d%02d', YEAR, MONTH, DAY)))

# Add a  DOY column
df$DOY <- yday(df$DATE)

# Rename "mean" column to EVI
df$EVI <- df$mean

#Remove column
df <- df[-4]

# Reorder columns
df <- df %>% select(MEAN_EVI, everything())
df <- df %>% select(EVI, everything())
df <- df %>% select(DOY, everything())
df <- df %>% select(DATE, everything())
df <- df %>% select(DAY, everything())
df <- df %>% select(MONTH, everything())
df <- df %>% select(YEAR, everything())

if(interpolate == 1){
  # Create a list of 365 days and merge with dataframe, adding NA's to days without MODIS EVI data  
  day_list <- as.data.frame(list(c(1:365)))
  names(day_list) <- c("Day")
  df_day_list <- merge(df, day_list, by = "Day", all=TRUE)
  df_day_list$ID <- 1:nrow(df_day_list)
  df <- df_day_list
  
  # Rename column headings
  # data.table::setnames(df, old=c('Day', 'X2002', 'X2003', 'X2004', 'X2005', 'X2006', 'X2007', 'X2008', 'X2009', 'X2010', 'X2011', 'X2012', 'X2013', 'X2014', 'X2015', 'X2016', 'X2017', 'X2018'),
  #                     new=c('DOY', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018'))
  
  # Re-format data table so that row-labels are Years
  df_gather <- df %>%
    gather(key=YEAR, value = MEAN_EVI, -DOY, -ID)
  df_gather$YEAR <- as.numeric(df_gather$YEAR)
  df_gather$DOY <- as.numeric(df_gather$DOY)
} else {
  
  # Rename column headings
  # data.table::setnames(df, old=c('Day', 'X2002', 'X2003', 'X2004', 'X2005', 'X2006', 'X2007', 'X2008', 'X2009', 'X2010', 'X2011', 'X2012', 'X2013', 'X2014', 'X2015', 'X2016', 'X2017', 'X2018'),
  #                      new=c('DOY', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018'))
  # Add an ID column
  df$ID <- 1:nrow(df)
  
  # Re-format data table so that row-labels are Years
  # df_gather <- df %>%
  #   gather(key=YEAR, value = MEAN_EVI, -DOY, -ID)
  df_gather <- df[-c(2:4,6)]
  df_gather$YEAR <- as.numeric(df_gather$YEAR)
  df_gather$DOY <- as.numeric(df_gather$DOY)
}

# Perform cubic spline interpolation to convert NA values to "EVI" increments 
df_approx <- zoo::na.spline(df_gather, na.rm = FALSE)
df_gather <- as.data.frame(df_approx)

# Remve ID field
df_gather <- df_gather[-4]
### Re-format table back to original layout in preparation for inserting days or smoothing
df_gather <- df_gather %>% group_by(YEAR) %>% mutate(ID = row_number())

# Create an empty data frame
df_EVI_final <- data.frame(ID=numeric(),
                           YEAR=numeric(),
                           DOY=numeric(),
                           MEAN_EVI=numeric(),
                           FITTED=numeric(),
                           SE=numeric(),
                           stringsAsFactors=FALSE)

# Create an empty data frame
df_PTest_final <- data.frame(
                           #ID=numeric(),
                           DOY=numeric(),
                           PTEST=numeric(),
                           stringsAsFactors=FALSE)

# Set up FOR loop to process each year sequentially for LOESS smoothing
for (YEAR in c(startyear:2019)) {
  # Filter to specified year
  year_to_process = YEAR
  df_year <- df_gather %>% 
    filter(YEAR == year_to_process)
  #Remove old ID column
  df_year <- df_year[-4]
  # Add an ID column
  df_year$ID <- as.numeric(1:nrow(df_year))
  
  # Smooth time series with LOESS
  df_year_loess <- loess(df_year$MEAN_EVI ~ df_year$DOY, span=loessSPAN)
  df_year_fitted_Pred <- predict(df_year_loess, se=T)
  
  # Convert "fitted" and "SE" fields to dataframe
  df_year_fitted <- data.frame("FITTED"=df_year_fitted_Pred$fit, "SE"=df_year_fitted_Pred$se.fit)
  
  # Add an ID column
  df_year_fitted$ID <- as.numeric(1:nrow(df_year_fitted))
  
  # Merge  YEAR/ID dataframe with FITTED dataframe
  df_year_final <- merge(df_year, df_year_fitted, by="ID")

  # Append dataframes from different years together into one output dataframe
  df_EVI_final <- bind_rows(df_EVI_final, df_year_final)
}  

df_EVI_final_MOD <- df_EVI_final

# Set up FOR loop to process each DOY sequentially for T-test
for (DOY in c(1,17,25,33,41,49,57,65,73,81,89,97,105,113,121,129,137,145,153,161,169,177,185,193,201,209,217,225,233,241,249,257,265,273,281,289,297,305,313,321,329,337,345,353,361)) {
 # Filter to specified year
  #DOY <- 17 # temporary bit to test FOR loop code -- remove when done
  DOY_to_process = DOY
  df_EVI_DOY <- df_EVI_final %>% 
    filter(DOY == DOY_to_process)
  
  #Remove old ID column
  df_EVI_DOY <- df_EVI_DOY[-1]
  # Add an ID column
  df_EVI_DOY$ID <- as.numeric(1:nrow(df_EVI_DOY))
  
  df_EVI_DOY <- df_EVI_DOY[-c(1),]
  
  # Stratify into two groups, before 2012 and after 2012
  df_EVI_DOY_filt1 <- filter(df_EVI_DOY, YEAR <2012)
  df_EVI_DOY_filt1$GROUP <- 'A'
  
  df_EVI_DOY_filt2 <- filter(df_EVI_DOY, YEAR >2011)
  df_EVI_DOY_filt2$GROUP <- 'B'
  
  # Merge two groups into single dataframe
  df_EVI_DOY <- bind_rows(df_EVI_DOY_filt1, df_EVI_DOY_filt2)

  ## Run basic stats 
  # group_by(df_EVI_DOY, GROUP) %>%
  #   summarise(
  #     count = n(),
  #     mean = mean(MEAN_EVI, na.rm = TRUE),
  #     sd = sd(MEAN_EVI, na.rm = TRUE)
  #   )
  
  # Subset MEAN_EVI data into before 2012
  # before <- subset(df_EVI_DOY,  GROUP == "A", MEAN_EVI,
  #                drop = TRUE)
  # subset MEAN_EVI data into after 2012
  # after <- subset(df_EVI_DOY,  GROUP == "B", MEAN_EVI,
  #                drop = TRUE)
  
  ## Plot paired data
  # pd <- PairedData::paired(before, after)
  # plot(pd, type = "profile") + theme_bw()
  
  ## Test for normality
  # d <- with(df_EVI_DOY, 
  #           MEAN_EVI[GROUP == "A"] - MEAN_EVI[GROUP == "B"])
  # 
  # ST <- shapiro.test(d)
  
  # Compute t-test for unequal variances between two time series (2003-2011 and 2012-2019)
  PTtest <- t.test(FITTED ~ GROUP, data = df_EVI_DOY, paired = FALSE, var.equal= FALSE)
  # Convert "p.value" field to dataframe
  PTtest2 <- data.frame("T.test"=PTtest$p.value)
  #Add DOY_to_process to new column in data frame
  PTtest2$DOY <- DOY_to_process 

  # Compute linear model for each DOY across all years (2003-2019)
  TTest = lm(FITTED ~ YEAR, df_EVI_DOY)
  ## PTest = lm(FITTED ~ YEAR -1, df_EVI_DOY_filt2) # Omit intercept
  # Summarize LM results
  TTest2 <- summary(TTest, correlation = FALSE, symbolic.cor = FALSE) 
  # Write coefficients to dataframe
  TTest3 <- data.frame("TTEST"=TTest2$coefficients) 
  # Remove row (used when intercept is not omitted)
  TTest3 <- TTest3[-c(1),] 
  #Save P.vale in new dataframe
  TTest4 <- data.frame("P.value"=TTest3$TTEST.Pr...t..) 
  #Add P.value to new column in data frame
  PTtest2$P.value <- TTest4$P.value 
  
  # Append Pdataframes from different years together into one output dataframe
  df_PTest_final <- bind_rows(df_PTest_final, PTtest2)
  
  }  

#Remove ID column
df_EVI_final2 <- subset(df_EVI_final, select = -c(ID, MEAN_EVI, SE))

### Re-format table back to original layout
df_EVI_final2 <- df_EVI_final2 %>% group_by(YEAR) %>% mutate(ID = row_number())

# Create a data frame with only DOY
df_DOY <- as.data.frame(df_year_final[,3, drop=FALSE])

#Add ID field to table
df_DOY$ID <- 1:nrow(df_DOY)

# Use SPREAD to reformat data table without the duplicated DOY 
df_EVI_final_spread <- df_EVI_final2 %>% 
  group_by(YEAR) %>% 
  select(-DOY) %>% 
  spread(key = YEAR, value = FITTED)

# Merge the DOY data table with the SPREAD data table
df_EVI_final_spread <- merge(df_EVI_final_spread, df_DOY, by = "ID", all=TRUE)

# Move DOY to first column
df_EVI_final_spread <- df_EVI_final_spread %>%
  select("DOY", everything())

# Remove ID column
df_EVI_final_spread <- df_EVI_final_spread[,-2]

# Merge with df_PTest_final
df_EVI_final_spread<- merge(df_EVI_final_spread, df_PTest_final, by="DOY")

# write output to a .csv file
write.csv(df_EVI_final_spread, file = (sprintf("~/R_Scripts/MODIS/data_out/%s_MODIS_%s_EVI_QAbit_%s_Final%s_SPAN_%s.csv",EVIFileName, sensor, QAbits, YEAR, loessSPAN)))


