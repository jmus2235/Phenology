##############################################################################################
#' @MODIS_EVI_processing_v1-2.R
#
#' @author
#' John Musinsky 1 \emailjmusinsky@battelleecology.org}
#
#' @description 
## Script for processing smoothed EVI 8-day composite time series from high-quality (QA-bit=0) MODIS Terra and Aqua data 
## produced for ROI's from Earth Engine script "GEE Mean 2003-2019 MODIS VIIRS EVI per Site - JM v4B"  
#
#' @param inp01 
## Input CSV table must contain three columns of data in following format: 
##
## system:index	date	      mean
## 2003_01_01   1.04138E+12	885.0715865
## 2003_01_17	  1.04276E+12	867.7747245
## 2003_02_02	  1.04414E+12	869.3154614
## ...
#
#' @return  
## Returns CSV table that includes the DOY, and smoothed, gap-filled EVI values for each year in the time-series 
#	DOY	2003	      2004	      2005	      2006	      2007	      2008	      2009	      2010	      2011	      2012	      2013	      2014	      2015	      2016	      2017	      2018	      2019
# 1	  0.087721705	0.087105282	0.091742409	0.095520278	0.111711909	0.091755388	0.099825213	0.082403387	0.098211333	0.084425081	0.087809931	0.101072331	0.090517611	0.083372446	0.088578561	0.101498778	0.087828407
#	9	  0.088112414	0.087409603	0.091953517	0.096347006	0.108610876	0.092284827	0.099889071	0.081321816	0.098195721	0.084413471	0.088865979	0.100576407	0.08887948	0.084609433	0.087328277	0.100147609	0.089110182
#	17	0.088655321	0.08760886	0.092985998	0.09695912	0.106588358	0.092669246	0.100009514	0.08076124	0.098156793	0.084527096	0.0898208	0.100191322	0.087696089	0.085865051	0.086641346	0.098806365	0.090082507
#	25	0.089353914	0.087657481	0.094852792	0.097384861	0.105627903	0.092933187	0.100168628	0.080665651	0.09810681	0.08479587	0.09065464	0.099920512	0.086953476	0.087106351	0.086530068	0.097485657	0.090747864
#	33	0.090233421	0.087576017	0.097519494	0.097666464	0.105734205	0.093154606	0.100365457	0.08101555	0.098059797	0.085239836	0.091365778	0.099790358	0.086640967	0.088331317	0.087079731	0.096274064	0.091261946
#
#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007
#
#' @keywords phenology, EVI, MODIS
#
#' @seealso "GEE Mean 2003-2019 MODIS VIIRS EVI per Site - JM v4B" GEE script
#
#' @export
#
# changelog and author contributions / copyrights
#   Author 1 (2017-02-14)
#     original creation
##############################################################################################

## ----install-packages, message = FALSE-----------------------------------
# uncomment for install
# install.packages("lubridate")
# install.packages("changepoint")
# install.packages("readr")
# install.packages("ggplot2")
# install.packages("tidyverse")
##
## Script for processing smoothed EVI 8-day composite time series from high-quality (QA-bit=0) MODIS Terra and Aqua data 
## produced from Earth Engine script "GEE Mean 2003-2019 MODIS VIIRS EVI per Site - JM v4B" 
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

# Set working directory
setwd("~/R_Scripts/MODIS/data")

#### USER-SPECIFIED PARAMETERS
EVIFileName = 'D14_JORN_FB'
QAbits = '0'
sensor = 'TERRA_AQUA'
#sensor = 'TERRA'
#sensor = 'AQUA'
year_to_process <- '2003'
loessSPAN <- 0.20
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

# Set up FOR loop to process each year sequentially
for (YEAR in c(2003:2019)) {
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

#Remove ID column
df_EVI_final2 <- subset(df_EVI_final, select = -c(ID, MEAN_EVI, SE))

### Re-format table back to original layout in preparation for inserting days or smoothing
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

# write output to a .csv file
write.csv(df_EVI_final_spread, file = (sprintf("~/R_Scripts/MODIS/data_out/%s_MODIS_%s_EVI_QAbit_%s_Final%s_SPAN_%s.csv",EVIFileName, sensor, QAbits, YEAR, loessSPAN)))

