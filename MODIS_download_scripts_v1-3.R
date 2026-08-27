# load the library
library(MODISTools)
library(dplyr)
# library(ggplot2)

# Set working directory
setwd("~/R_Scripts/MODIS/data")

#USER DEFINED PARAMETERS
BoundBox = 0

# McdBand = "MidGreenup.Num_Modes_01"
# PhenoPhase = "MOS"

# McdBand = "Maturity.Num_Modes_01"
# PhenoPhase = "EOS"

McdBand = "Dormancy.Num_Modes_01"
PhenoPhase = "EOF"

# McdBand = "MidGreendown.Num_Modes_01"
# PhenoPhase = "MOF"

# Read csv file
FltboxCP <- read.csv(file="Flightbox_Centerpoints.csv", stringsAsFactors = FALSE)
FltboxCP$ID <- 1:nrow(FltboxCP)

df_final <- data.frame(DomainName=character(),
                       SiteName=character(),
                       meanDOY=integer(),
                       meanDOYdate=character(),
                       FltBoxX=integer(),
                       FltBoxY=integer(),
                       stringsAsFactors=FALSE)

for (i in 1:nrow(FltboxCP)) {

#PRE-DEFINED PARAMETERS
  DomainName <-  data.frame(FltboxCP[i,1])
  SiteName <-  data.frame(FltboxCP[i,2])
  Latitude <- data.frame(FltboxCP[i,3])
  Longitude <-  data.frame(FltboxCP[i,4])
  FltBoxX <-  data.frame(FltboxCP[i,5])
  FltBoxY <-  data.frame(FltboxCP[i,6])
  FltBoxID <-  data.frame(FltboxCP[i,7])

    

# See which MODIS products are available
# products <- mt_products()
# head(products)

# View available bands in MCD12Q2
# bands <- mt_bands(product = "MCD12Q2")
# print(str(bands))

# List available dates
# dates <- mt_dates(product = "MCD12Q2", lat = 44.036956, lon = -71.271877)
# head(dates)


# download data
subset <- mt_subset(product = "MCD12Q2",
                    lat = Latitude,
                    lon = Longitude,
                    band = McdBand,
                    start = "2001-01-01",
                    end = "2017-01-01",
                    km_lr = FltBoxX,
                    km_ab = FltBoxY,
                    site_name = SiteName,
                    internal = TRUE,
                    progress = FALSE)
print(str(subset))

# convert the coordinates
lat_lon <- sin_to_ll(subset$xllcorner, subset$yllcorner)

# bind with the original dataframe
subset <- cbind(subset, lat_lon)

head(subset)

subset$origin <- (as.Date(as.POSIXct(subset$value, origin="1970-01-01")))
subset$maturity_date <- as.Date(subset$value, origin = subset$origin, tz = "UTC")
#subset <- subset[,-24]

if (BoundBox == 1) {
  
# convert to bounding box
bb <- apply(subset, 1, function(x){
  mt_bbox(xllcorner = x['xllcorner'],
          yllcorner = x['yllcorner'],
          cellsize = x['cellsize'],
          nrows = x['nrows'],
          ncols = x['ncols'])
})

# plot one bounding box
plot(bb[[1]])

# add the location of the queried coordinate within the polygon
points(subset$longitude[1],
       subset$latitude[1],
       pch = 20,
       col = "red")
} else {
}

# create data frame with all MODIS year calandar dates and associated maturity dates per pixel in UNIX year zero format
subset2 <- subset["calendar_date"]
subset2$value <- subset$value
subset2$calendar_date <- as.Date(subset2$calendar_date)
# subset2$calendar_date <- as.character(subset2$calendar_date)

subset2$origin <- (as.Date(as.POSIXct(subset2$value, origin="1970-01-01")))
subset2$value_date <- as.Date(subset2$value, origin = subset2$origin, tz = "UTC")

# add a DOY column
subset2$valueDate_DOY <- lubridate::yday(as.Date(subset2$value_date,format='%d/%m/%Y'))

# create data frame with a single mean maturity date for all pixels
df_subset2 <- subset2 %>% 
  group_by(calendar_date) %>% 
  summarize(mean_value = mean(valueDate_DOY, na.rm = TRUE))

# another way to get mean values using the aggregate function
# aggdata <- aggregate(subset2, by=list(subset2$valueDate_DOY), FUN=mean, na.rm=TRUE)

# round the result of the calculated mean
df_subset2$matDate_DOY <- round(df_subset2$mean_value)

# delete column
df_subset <- df_subset2[,-2]

# Add a new column with a readable date format for the maturity date
df_subset$origin <- (as.Date(as.POSIXct(df_subset$matDate_DOY, origin="1970-01-01")))
df_subset$maturity_date <- as.Date(df_subset$matDate_DOY, origin = df_subset$origin, tz = "UTC")

# delete column
df_subset <- df_subset[,-3]

# add a DOY column
#df_subset$matDate_DOY <- lubridate::yday(as.Date(df_subset$maturity_date,format='%d/%m/%Y'))

# calculate mean DOY and Date for printing to screen
meanDOY <- round(mean(df_subset$matDate_DOY))
meanDOYdate <- as.Date(meanDOY, origin = "2020-01-01", tz = "UTC")

df_interum <- data.frame(DomainName, SiteName, meanDOY, meanDOYdate, FltBoxX, FltBoxY)

df_final <- rbind(df_final, df_interum)
  
# Print result
print(SiteName)
print(df_subset)
print(paste0("Mean DOY: ", meanDOY))
print(paste0("Mean DOY Date: ", meanDOYdate))

  }

# write output to a .csv file
write.csv(df_final, file = (sprintf("~/R_Scripts/MODIS/data_out/MCD12Q2_Transition_Dates_%s.csv", PhenoPhase)))




