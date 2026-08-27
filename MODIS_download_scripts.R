# load the library
library(MODISTools)

# Set working directory
setwd("~/R_Scripts/MODIS/data")

# See which MODIS products are available
products <- mt_products()
head(products)

# View available bands in MCD12Q2
bands <- mt_bands(product = "MCD12Q2")
print(str(bands))

# List available dates
dates <- mt_dates(product = "MCD12Q2", lat = 44.036956, lon = -71.271877)
head(dates)

# download data
subset <- mt_subset(product = "MCD12Q2",
                    lat = 44.036960,
                    lon = -71.271692,
                    band = "Maturity.Num_Modes_01",
                    start = "2001-01-01",
                    end = "2017-01-01",
                    km_lr = 5,
                    km_ab = 5,
                    site_name = "BART",
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

# create data frame with all MODIS year calandar dates and associated maturity dates per pixel in UNIX year zero format
subset2 <- subset["calendar_date"]
subset2$value <- subset$value
subset2$calendar_date <- as.Date(subset2$calendar_date)

# create data frame with a single mean maturity date for all pixels
df_subset2 <- subset2 %>% 
  group_by(calendar_date) %>% 
  summarize(mean_value = mean(value, na.rm = TRUE))

# round the result of the calculated mean
df_subset2$new_mean_value <- round(df_subset2$mean_value)

# delete column
df_subset <- df_subset2[,-2]

# Add a new column with a readable date format for the maturity date
df_subset$origin <- (as.Date(as.POSIXct(df_subset$new_mean_value, origin="1970-01-01")))
df_subset$maturity_date <- as.Date(df_subset$new_mean_value, origin = df_subset$origin, tz = "UTC")

# delete column
df_subset <- df_subset[,-3]

# Print result
print(df_subset)

# UnpackDetailedQA <- function(x){
#   bits <- as.integer(intToBits(x))
#   quals <- sapply(seq(1, 16, by=2), function(i) sum(bits[i:(i+1)] * 2^c(0, 1)))[1:7]
#   return(quals)
# }
