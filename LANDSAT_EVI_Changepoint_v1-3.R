##############################################################################################
#' @LANDSAT_EVI_Changepoint_v1-3.R

#' @author
#' John Musinsky 1 \emailjmusinsky@battelleecology.org}

#' @description 
## Script for producing smoothed EVI phenophase change-point transitions from 8-day composite time series Landsat 8 Collection 1  
## Tier 1 8-Day EVI Composite produced from Earth Engine script "GEE_Mean 2018 Landsat 8 EVI per Masked Site - JM" 

#' @param inp01 
## Input CSV table must contain three columns of data in following format: 
##
## Day	2002	    2003	    2004	    2005	    2006	    2007	    2008	    2009	    2010	    2011	    2012	    2013	    2014	    2015	    2016	    2017	    2018
## 1	  0.224136	0.169765	0.248993	0.228492	0.236967	0.237593	0.211651	0.186319	0.250911	0.297741	0.241411	0.223491	0.277903	0.241074	0.236391	0.225876	0.216436
## 9	  0.219643	0.180689	0.225491	0.220534	0.236496	0.23885	  0.206401	0.185824	0.234282	0.244129	0.247069	0.22869	  0.257781	0.226136	0.259102	0.215989	0.213553
## 17	  0.22022	  0.186139	0.204644	0.213796	0.239709	0.236334	0.198796	0.18518	  0.220915	0.205435	0.250088	0.230804	0.23639	  0.212631	0.2714	  0.210892	0.212738
## ... 

#' @return Description of what the function returns [units] 

#' @references
#' License: (example) GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @keywords keyword 1, keyword 2

#' @examples Example code for using your function

#' @seealso References to other related functions

#' @export

# changelog and author contributions / copyrights
#   Awesome Author 1 (2017-02-14)
#     original creation
#   Awesome Author 2 (2017-02-15)
#     Made some awesome changes
##############################################################################################
## ----install-packages, message = FALSE-----------------------------------
# uncomment for install
# install.packages("devtools")
# install.packages("maps")
# install.packages("tidyverse")

##


library(readr)
library(ggplot2)
library(changepoint)
library(lubridate)
library(tidyverse)
library(dplyr)

# Set working directory
setwd("~/R_Scripts/phenocamr-master/data")

#### USER-SPECIFIED PARAMETERS
EVIFileName = "D01_BARTgt500mask"
#year_to_process <- '2015'
loess_span = 0.06
penalty <- .05
seg_length <- 14
interpolate <- 1 # 0 for 'No', 1 for 'Yes

# Read csv file
dfa <- read.csv(file=sprintf("%s_L8_EVI_statistics_2015-2018.csv",EVIFileName), stringsAsFactors = FALSE)

# Rename column headings
data.table::setnames(dfa, old=c('Day', 'X2015', 'X2016', 'X2017', 'X2018'),
                     new=c('DOY', '2015', '2016', '2017', '2018'))

# Re-format data table so that row-labels are Years
dfa_gather <- dfa %>%
  gather(key=YEAR, value = MEAN_EVI, -DOY)
dfa_gather$YEAR <- as.numeric(dfa_gather$YEAR)
dfa_gather$DOY <- as.numeric(dfa_gather$DOY)

# Perform cubic spline interpolation to convert NA values to 8 day "EVI" increments 
dfa_approx <- zoo::na.spline(dfa_gather, na.rm = FALSE)
dfa_gather <- as.data.frame(dfa_approx)

### Re-format table back to original layout in preparation for inserting days or smoothing
dfa_gather <- dfa_gather %>% group_by(YEAR) %>% mutate(ID = row_number())

# Create a data frame with only DOY
dfa_DOY <- select(dfa, -"2015", -"2016", -"2017", -"2018")
#Add ID field to table
dfa_DOY$ID <- 1:nrow(dfa_DOY)

# Use SPREAD to reformat data table without the duplicated DOY 9
dfa_spread <- dfa_gather %>% 
  group_by(YEAR) %>% 
  select(-DOY) %>% 
  spread(key = YEAR, value = MEAN_EVI)

# Merge the DOY data table with the SPREAD data table
dfa_spread_final <- merge(dfa_spread, dfa_DOY, by = "ID", all=TRUE)

#### Write 8-day NA interpolated pre-smoothed data for processing in Matlab fLOESS
write.csv(dfa_spread_final, file=(sprintf("~/R_Scripts/phenocamr-master/data_out/%s_L8_EVI_statistics_2015-2018_8day_unsmoothed_filled.csv",EVIFileName)))

if(interpolate == 1){

# Create a list of 365 days and merge with dataframe, adding NA's to days without MODIS EVI data  
day_list <- as.data.frame(list(c(1:365)))
names(day_list) <- c("DOY")
dfa_day_list <- merge(dfa_spread_final, day_list, by = "DOY", all=TRUE)
dfa_day_list$ID <- 1:nrow(dfa_day_list)
dfa_360 <- dfa_day_list

# Re-format data table so that row-labels are Years
dfa_gather2 <- dfa_360 %>%
  gather(key=YEAR, value = MEAN_EVI, -DOY, -ID)
dfa_gather2$YEAR <- as.numeric(dfa_gather2$YEAR)
dfa_gather2$DOY <- as.numeric(dfa_gather2$DOY)

# Calculate dates from Day of Year (DOY) and Year
df <- data.frame(dayyear = c(1, 100, 300, 1), year = c(2013, 2013, 2013, 2014))
dfa_gather2$ORIGIN <- as.Date(paste0(dfa_gather2$YEAR, "-01-01"),tz = "UTC") - days(1)
dfa_gather2$DATE <- as.Date(dfa_gather2$DOY, origin = dfa_gather2$ORIGIN, tz = "UTC")  

# Delete "origin" column
dfa_gather2 <- dfa_gather2[,-2]
dfa_gather2 <- dfa_gather2[,-4]
dfa_gather2$ID <- 1:nrow(dfa_gather2)

# Perform cubic spline interpolation to convert NA values to 1 day "EVI" increments 
dfa_approx2 <- zoo::na.spline(dfa_gather2$MEAN_EVI, na.rm = FALSE)
dfa_approx2 <- as.data.frame(dfa_approx2)
dfa_approx2$ID <- 1:nrow(dfa_approx2)
dfa_gather3 <- merge(dfa_gather2, dfa_approx2, by = "ID", all=FALSE)
dfa_gather3 <- dfa_gather3[,-4]

# Rename column headings
data.table::setnames(dfa_gather3, old=c('dfa_approx2'),new=c('MEAN_EVI'))
dfa_gather <- dfa_gather3
} else {
#   ### Re-format table back to original layout in preparation for inserting days or smoothing
#   
#   # Use SPREAD to reformat data table without the duplicated DOY 9
#   dfa_gather$ID <- 1:nrow(dfa_gather)
#   
  dfa_gather <- dfa_gather[-4]
  dfa_gather$ID <- 1:nrow(dfa_gather)
  # 
#   # Create a data frame with only DOY
#   dfa_DOY <- select(dfa, -"2015", -"2016", -"2017", -"2018")
#   #Add ID field to table
#   dfa_DOY$ID <- 1:nrow(dfa_DOY)
#   
#   dfa_spread <- dfa_gather %>% 
#     group_by(YEAR) %>% 
#     select(-DOY) %>% 
#     spread(key = YEAR, value = MEAN_EVI)
#   
#     #### Write 8-day NA interpolated pre-smoothed data for processing in Matlab fLOESS
#   write.csv(dfa_gather, file=(sprintf("~/R_Scripts/phenocamr-master/data_out/%s_L8_EVI_statistics_2015-2018_8day_unsmoothed_filled.csv",EVIFileName)))
}

# Plot unsmoothed data
#ggplot(dfa_plot, aes(DOY, MEAN_EVI))+
ggplot(dfa_gather, aes(DOY, MEAN_EVI))+
  geom_point(aes(color = YEAR))

# Run LOESS smoothing algorithm
loess_dfa <- loess(MEAN_EVI ~ ID, data=dfa_gather, normalize = TRUE, span = loess_span)

# Predict LOESS
smoothed20 <- predict(loess_dfa)

# Convert smoothed data to a data frame
smoothed20 <- as.data.frame(smoothed20, stringsAsFactors = default.stringsAsFactors())

# Adda a column header
names(smoothed20) <- c("SMOOTHED_MEAN_EVI")

#Add ID fields to table
smoothed20$ID <- 1:nrow(smoothed20)

# Merge tables
#dfa_gather <- merge(dfa_gather, smoothed20, by = "ID", all=TRUE) %>% filter(YEAR==2018)
dfa_gather <- merge(dfa_gather, smoothed20, by = "ID", all=TRUE)

if(interpolate == 1){

#### Reformat and write 365-day LOESS smoothed data forcomparing to unsmoothed data processed in Matlab fLOESS
dfa_gather_reformat <- dfa_gather %>% 
  group_by(YEAR) %>% 
  select(-ID, -MEAN_EVI, -DATE) %>% 
  spread(key = YEAR, value = SMOOTHED_MEAN_EVI)
#dfa_gather_reformat <- merge(dfa_gather_reformat, dfa_DOY, by = "ID", all=TRUE)
#dfa_gather_reformat <- dfa_gather_reformat[,-1]
write.csv(dfa_gather_reformat, file=(sprintf("~/R_Scripts/phenocamr-master/data_out/%s_L8_EVI_statistics_2015-2018_365day_LOESS_smoothed.csv",EVIFileName)))
} else {
    dfa_gather_reformat <- dfa_gather %>% 
      group_by(YEAR) %>% 
      select(-MEAN_EVI, -ID) %>% 
      spread(key = YEAR, value = SMOOTHED_MEAN_EVI)
  write.csv(dfa_gather_reformat, file=(sprintf("~/R_Scripts/phenocamr-master/data_out/%s_L8_EVI_statistics_2015-2018_8day_LOESS_smoothed.csv",EVIFileName)))
}

# Plot smoothed data
ggplot(dfa_gather, aes(DOY, SMOOTHED_MEAN_EVI))+
  geom_point(aes(color = YEAR))

# Set up FOR loop to process each year sequentially
for (YEAR in c(2015:2018)) {
  year_to_process = YEAR
dfa_year <- dfa_gather %>% 
  filter(YEAR == year_to_process)

# Filter to specified year
dfa_cpt <- dfa_gather %>%
  filter(YEAR == year_to_process)

# Perform change.point analysis
cpt_man = cpt.mean(
  as.numeric(dfa_cpt$SMOOTHED_MEAN_EVI),
  method = 'PELT',
  penalty='Manual',
  pen.value = penalty,
  minseglen = seg_length)

# Save change.point locations to file
cpts_loc <- c(cpts(cpt_man))

# Convert change.point locations list to a data frame - RISING
list_cpts <- as.list(cpts_loc)
dfa_from_list <- data.frame(matrix(unlist(list_cpts), nrow=length(list_cpts), byrow=T), stringsAsFactors = FALSE)
colnames(dfa_from_list) <- c('ID')
#colnames(dfa_from_list) <- c('DOY')

# Add 1 to change.point location to calculate correct start of mean
dfa_from_list$NEW_ID <- dfa_from_list$ID + 0

# Filter the original MEAN_EVI data table to just the change.point dates
dfa_year_rising <- dfa_year %>% 
  filter(dfa_year$DOY %in% dfa_from_list$NEW_ID)

# Prepare header info
if(interpolate == 1){
  cpts_header1 = paste(
    "# Interpolated daily EVI data\n",
    sprintf("# LOESS Span = %s\n", loess_span),
    sprintf("# Change.points penalty = %s\n",penalty),
    sprintf("# Segment length = %s\n",seg_length),
    "# EOS - Rising\n",
    "#",
    sep='')
} else {
  cpts_header1 = paste(
    "#Un-interpolated 8-day EVI data\n",
    sprintf("# LOESS Span = %s\n", loess_span),
    sprintf("# Change.points penalty = %s\n",penalty),
    sprintf("# Segment length = %s\n",seg_length),
    "# EOS - Rising\n",
    "#",
    sep='')
}

# Specify output filename
writeFile1 = sprintf("~/R_Scripts/phenocamr-master/data_out/%s_L8_EVI_statistics_%s_rising_penalty%s_segmentlength%s_span%s.csv",EVIFileName, year_to_process, penalty, seg_length, loess_span)

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
write.cpts(dfa_year_rising, writeFile1, cpts_header1)

# Convert change.point locations list to a data frame - FALLING
list_cpts <- as.list(cpts_loc)
dfa_from_list <- data.frame(matrix(unlist(list_cpts), nrow=length(list_cpts), byrow=T), stringsAsFactors = FALSE)
colnames(dfa_from_list) <- c('ID')

# Subtract 1 to change.point location to calculate correct start of mean
dfa_from_list$NEW_ID <- dfa_from_list$ID - 0

# Filter the original MEAN_EVI data table to just the change.point dates
dfa_year_falling <- dfa_year %>% 
  filter(dfa_year$ID %in% dfa_from_list$NEW_ID)

# Prepare header info
if(interpolate == 1){
  cpts_header2 = paste(
    "# Interpolated daily EVI data\n",
    sprintf("# LOESS Span = %s\n", loess_span),
    sprintf("# Change.points penalty = %s\n",penalty),
    sprintf("# Segment length = %s\n",seg_length),
    "# SOF - Falling\n",
    "#",
    sep='')
} else {
  cpts_header2 = paste(
    "#Un-interpolated 8-day EVI data\n",
    sprintf("# LOESS Span = %s\n", loess_span),
    sprintf("# Change.points penalty = %s\n",penalty),
    sprintf("# Segment length = %s\n",seg_length),
    "# SOF - Falling\n",
    "#",
    sep='')
}

# Specify output filename
writeFile2 = sprintf("~/R_Scripts/phenocamr-master/data_out/%s_L8_EVI_statistics_%s_falling_penalty%s_segmentlength%s_span%s.csv",EVIFileName, year_to_process, penalty, seg_length, loess_span)

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
write.cpts(dfa_year_falling, writeFile2, cpts_header2)

# Plot and save plot to .jpg file
plot(cpt_man,cpt.width=.5)
jpeg(file=sprintf("~/R_Scripts/phenocamr-master/data_out/%s_L8_EVI_statistics_%s_rising_penalty%s_segmentlength%s_span%s.jpg",EVIFileName, year_to_process, penalty, seg_length, loess_span))
plot(cpt_man,cpt.width=.5, main=year_to_process)
dev.off()

if(year_to_process == '2015'){
  dfa_2015 <- dfa_year_rising
} else if(year_to_process == '2016'){
  dfa_2016 <- dfa_year_rising
}  else if(year_to_process == '2017'){
    dfa_2017 <- dfa_year_rising
} else {dfa_2018 <- dfa_year_rising}

} # End FOR loop

# Display rising change.point locations
dfa_2015
dfa_2016
dfa_2017
dfa_2018

