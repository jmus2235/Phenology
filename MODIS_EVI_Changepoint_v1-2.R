## ----install-packages, message = FALSE-----------------------------------
# uncomment for install
# install.packages("lubridate")
# install.packages("changepoint")
# install.packages("readr")
# install.packages("ggplot2")
# install.packages("tidyverse")
##
## Script for producing smoothed EVI phenophase change-point transitions from 8-day composite time series MODIS Terra and Aqua data 
## produced from Earth Engine script "GEE_Mean_MODIS_EVI_per_site" 
## Input CSV table must contain three columns of data in following format: 
##
## Day	2002	    2003	    2004	    2005	    2006	    2007	    2008	    2009	    2010	    2011	    2012	    2013	    2014	    2015	    2016	    2017	    2018
## 1	  0.224136	0.169765	0.248993	0.228492	0.236967	0.237593	0.211651	0.186319	0.250911	0.297741	0.241411	0.223491	0.277903	0.241074	0.236391	0.225876	0.216436
## 9	  0.219643	0.180689	0.225491	0.220534	0.236496	0.23885	  0.206401	0.185824	0.234282	0.244129	0.247069	0.22869	  0.257781	0.226136	0.259102	0.215989	0.213553
## 17	  0.22022	  0.186139	0.204644	0.213796	0.239709	0.236334	0.198796	0.18518	  0.220915	0.205435	0.250088	0.230804	0.23639	  0.212631	0.2714	  0.210892	0.212738
## ...

library(readr)
library(ggplot2)
library(changepoint)
library(lubridate)
library(tidyverse)

# Set working directory
setwd("~/R_Scripts/phenocamr-master/data")

#### USER-SPECIFIED PARAMETERS
EVIFileName = "D01_HARV"
year_to_process <- '2002'
penalty <- 0.005
seg_length <- 14
interpolate <- 0 # 0 for 'No', 1 for 'Yes'
####

## Standard Parameters
# penalty <- 0.005
# seg_length <- 4
# interpolate <- 0
# 
# penalty <- 0.005
# seg_length <- 3
# interpolate <- 0
# 
# penalty <- 0.005
# seg_length <- 14
# interpolate <- 1

# Read csv file
df <- read.csv(file=sprintf("%s_DB_EVI_smoothed.csv",EVIFileName), stringsAsFactors = FALSE)

if(interpolate == 1){
  # Create a list of 365 days and merge with dataframe, adding NA's to days without MODIS EVI data  
  day_list <- as.data.frame(list(c(1:365)))
  names(day_list) <- c("Day")
  df_day_list <- merge(df, day_list, by = "Day", all=TRUE)
  df_day_list$ID <- 1:nrow(df_day_list)
  df <- df_day_list
  
  # Rename column headings
  data.table::setnames(df, old=c('Day', 'X2002', 'X2003', 'X2004', 'X2005', 'X2006', 'X2007', 'X2008', 'X2009', 'X2010', 'X2011', 'X2012', 'X2013', 'X2014', 'X2015', 'X2016', 'X2017', 'X2018'),
                       new=c('DOY', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018'))
  
  # Re-format data table so that row-labels are Years
  df_gather <- df %>%
    gather(key=YEAR, value = MEAN_EVI, -DOY, -ID)
  df_gather$YEAR <- as.numeric(df_gather$YEAR)
  df_gather$DOY <- as.numeric(df_gather$DOY)
  } else {
    
    # Rename column headings
    data.table::setnames(df, old=c('Day', 'X2002', 'X2003', 'X2004', 'X2005', 'X2006', 'X2007', 'X2008', 'X2009', 'X2010', 'X2011', 'X2012', 'X2013', 'X2014', 'X2015', 'X2016', 'X2017', 'X2018'),
                         new=c('DOY', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018'))
    # Add an ID column
    df$ID <- 1:nrow(df)
    
    # Re-format data table so that row-labels are Years
    df_gather <- df %>%
      gather(key=YEAR, value = MEAN_EVI, -DOY, -ID)
    df_gather$YEAR <- as.numeric(df_gather$YEAR)
    df_gather$DOY <- as.numeric(df_gather$DOY)
}

if(interpolate == 1){
# Perform cubic spline interpolation to convert NA values to 1 day "EVI" increments 
    df_approx <- zoo::na.spline(df_gather, na.rm = FALSE)
    df_gather <- as.data.frame(df_approx)
  } else {
}

# Calculate dates from Day of Year (DOY) and Year
df <- data.frame(dayyear = c(1, 100, 300, 1), year = c(2013, 2013, 2013, 2014))
df_gather$ORIGIN <- as.Date(paste0(df_gather$YEAR, "-01-01"),tz = "UTC") - days(1)
df_gather$DATE <- as.Date(df_gather$DOY, origin = df_gather$ORIGIN, tz = "UTC")  

if(interpolate == 1){
  # Delete "origin" column
  df_gather <- df_gather[,-5]
  } else {
    df_gather <- df_gather[,-5]
}

# Set up FOR loop to process each year sequentially
for (YEAR in c(2002:2018)) {
year_to_process = YEAR
df_year <- df_gather %>% 
  filter(YEAR == year_to_process)

# Filter to specified year
df_cpt <- df_gather %>%
  filter(YEAR == year_to_process)

# Perform change.point analysis
cpt_man = cpt.mean(
  as.numeric(df_cpt$MEAN_EVI),
  method = 'PELT',
  penalty='Manual',
  pen.value = penalty,
  minseglen = seg_length)

# Save change.point locations to file
cpts_loc <- c(cpts(cpt_man))

# Convert change.point locations list to a data frame - RISING
list_cpts <- as.list(cpts_loc)
df_from_list <- data.frame(matrix(unlist(list_cpts), nrow=length(list_cpts), byrow=T), stringsAsFactors = FALSE)
colnames(df_from_list) <- c('ID')

# Add 1 to change.point location to calculate correct start of mean
df_from_list$NEW_ID <- df_from_list$ID + 0

# Filter the original MEAN_EVI data table to just the change.point dates
df_year_rising <- df_year %>% 
  filter(df_year$ID %in% df_from_list$NEW_ID)

# Prepare header info
if(interpolate == 1){
cpts_header1 = paste(
  "# Interpolated daily EVI data\n",
  sprintf("# Change.points penalty = %s\n",penalty),
  sprintf("# Segment length = %s\n",seg_length),
  "# EOS - Rising\n",
  "#",
  sep='')
} else {
  cpts_header1 = paste(
    "#Un-interpolated 8-day EVI data\n",
    sprintf("# Change.points penalty = %s\n",penalty),
    sprintf("# Segment length = %s\n",seg_length),
    "# EOS - Rising\n",
    "#",
    sep='')
  }

# Specify output filename
writeFile1 = sprintf("~/R_Scripts/phenocamr-master/data_out/%s_DB_EVI_CPTS_%s_rising_penalty%s_segmentlength%s.csv",EVIFileName, year_to_process, penalty, seg_length)

# Create function to write .csv file with header
write.cpts <- function(x, file, header, f = write.csv, ...){
  # create and open the file connection
  datafile <- file(file, open = 'wt')
  # close on exit
  on.exit(close(datafile))
  # if a header is defined, write it to the file 
  if(!missing(header)) writeLines(header,con=datafile)
  # write the file using the defined function and required additional arguments  
  f(x, datafile,...)
}

# Write filtered change.point data frame to a .csv file with header
write.cpts(df_year_rising, writeFile1, cpts_header1)

# Convert change.point locations list to a data frame - FALLING
list_cpts <- as.list(cpts_loc)
df_from_list <- data.frame(matrix(unlist(list_cpts), nrow=length(list_cpts), byrow=T), stringsAsFactors = FALSE)
colnames(df_from_list) <- c('ID')

# Subtract 1 to change.point location to calculate correct start of mean
df_from_list$NEW_ID <- df_from_list$ID - 0

# Filter the original MEAN_EVI data table to just the change.point dates
df_year_falling <- df_year %>% 
  filter(df_year$ID %in% df_from_list$NEW_ID)

# Prepare header info
if(interpolate == 1){
  cpts_header2 = paste(
    "# Interpolated daily EVI data\n",
    sprintf("# Change.points penalty = %s\n",penalty),
    sprintf("# Segment length = %s\n",seg_length),
    "# SOF - Falling\n",
    "#",
    sep='')
} else {
  cpts_header2 = paste(
    "#Un-interpolated 8-day EVI data\n",
    sprintf("# Change.points penalty = %s\n",penalty),
    sprintf("# Segment length = %s\n",seg_length),
    "# SOF - Falling\n",
    "#",
    sep='')
}

# Specify output filename
writeFile2 = sprintf("~/R_Scripts/phenocamr-master/data_out/%s_DB_EVI_CPTS_%s_falling_penalty%s_segmentlength%s.csv",EVIFileName, year_to_process, penalty, seg_length)

# Create function to write .csv file with header
write.cpts <- function(x, file, header, f = write.csv, ...){
  # create and open the file connection
  datafile <- file(file, open = 'wt')
  # close on exit
  on.exit(close(datafile))
  # if a header is defined, write it to the file
  if(!missing(header)) writeLines(header,con=datafile)
  # write the file using the defined function and required additional arguments  
  f(x, datafile,...)
}

# Write filtered change.point data frame to a .csv file with header
write.cpts(df_year_falling, writeFile2, cpts_header2)

# Plot and save plot to .jpg file
plot(cpt_man,cpt.width=.5)
jpeg(file=sprintf("~/R_Scripts/phenocamr-master/data_out/%s_DB_EVI_CPTS_%s_penalty%s_segmentlength%s.jpg",EVIFileName, year_to_process, penalty, seg_length))
plot(cpt_man,cpt.width=.5, main=year_to_process)
dev.off()

}

# Display latest year rising change.point locations
df_year_rising
