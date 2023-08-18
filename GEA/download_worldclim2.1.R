# download bioclim2.1 data

# Clean environment

rm(list = ls())

# Get libraries to extract Data:
library(raster)
library(rgdal)

setwd("/bioclim/wc0.5/")


# Set Link where to find the data:
bio_0.5 <-  "https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_bio.zip"
tmin <- "https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_tmin.zip"
tmax <- "http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_tmax.zip"
tavg <- "http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_tavg.zip"
prec <- "http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_prec.zip"
srad <- "http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_srad.zip"
wind <- "http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_wind.zip"
vapr <- "http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_vapr.zip"

# Downlaod the necessary zip files:
download.file(bio_0.5, destfile = "wc2.1_30s_bio.zip")
download.file(tmin, destfile = "wc2.1_30s_tmin.zip")
download.file(tmax, destfile = "wc2.1_30s_tmax.zip")
download.file(tavg, destfile = "wc2.1_30s_tavg.zip")
download.file(prec, destfile = "wc2.1_30s_prec.zip")
download.file(srad, destfile = "wc2.1_30s_srad.zip")
download.file(wind, destfile = "wc2.1_30s_wind.zip")
download.file(vapr, destfile = "wc2.1_30s_vapr.zip")

## Unzip all your files

bio_0.5  <- unzip("wc2.1_30s_bio.zip")
tmin <- unzip("wc2.1_30s_tmin.zip")
tmax <- unzip("wc2.1_30s_tmax.zip")
tavg <- unzip("wc2.1_30s_tavg.zip")
prec <- unzip("wc2.1_30s_prec.zip")
srad <- unzip("wc2.1_30s_srad.zip")
wind <- unzip("wc2.1_30s_wind.zip")
vapr <- unzip("wc2.1_30s_vapr.zip")

# Read file with coordinates 

t <- read.csv("bdis332_current_bioclim.csv")
head(t)

# data frame with the coordinates:
coords <- data.frame(longitude = t$longitude, latitude = t$latitude)

bioclim2 = c('Annual Mean Temperature',
             'Mean Diurnal Range',
             'Isothermality',
             'Temperature Seasonality',
             'Max Temperature of Warmest Month',
             'Min Temperature of Coldest Month',
             'Temperature Annual Range',
             'Mean Temperature of Wettest Quarter',
             'Mean Temperature of Driest Quarter',
             'Mean Temperature of Warmest Quarter',
             'Mean Temperature of Coldest Quarter',
             'Annual Precipitation',
             'Precipitation of Wettest Month',
             'Precipitation of Driest Month',
             'Precipitation Seasonality',
             'Precipitation of Wettest Quarter',
             'Precipitation of Driest Quarter',
             'Precipitation of Warmest Quarter',
             'Precipitation of Coldest Quarter',
             'minimum temperature',
             'maximum temperature',
             'average temperature',
             'precipitation (mm)',
             'solar radiation (kJ m-2 day-1)',
             'wind speed (m s-1)',
             'water vapor pressure (kPa)'
)


list.files()
files <- list.files(pattern = '\\.tif$', recursive = T)

allraster <- lapply(files,raster)
allraster[1]

# add bioclim variables to data frame 
i <- length(t) + 1
for (f in files) {
  r <- raster(f)
  name = names(r)
  
  # Use Coordinates to create a SpatialPoint object:
  points <- SpatialPoints(coords, proj4string = r@crs)
  
  # Use Extract with cbind and coords to get the desired data.frame:
  values <- extract(r, points, df = FALSE) 
  
  t <- cbind(t, values)  # t new dataframe with bdis and WorldClim2
  colnames(t)[i] <- name
  i <- i + 1
}

summary(t)
write.table(t, file = "/bioclim/wc0.5/my_data_0.5sec.csv",
            quote = F, sep = '\t', row.names = F)

head(t)
# mutate manually
dftest <- t %>%
  mutate(wc2.1_30s_prec_spring = (wc2.1_30s_prec_03 + wc2.1_30s_prec_04 + wc2.1_30s_prec_05 + wc2.1_30s_prec_06)/4) %>%
  mutate(wc2.1_30s_srad_spring = (wc2.1_30s_srad_03 + wc2.1_30s_srad_04 + wc2.1_30s_srad_05 + wc2.1_30s_srad_06)/4) %>%
  mutate(wc2.1_30s_tavg_spring = (wc2.1_30s_tavg_03 + wc2.1_30s_tavg_04 + wc2.1_30s_tavg_05 + wc2.1_30s_tavg_06)/4) %>%
  mutate(wc2.1_30s_tmax_spring = (wc2.1_30s_tmax_03 + wc2.1_30s_tmax_04 + wc2.1_30s_tmax_05 + wc2.1_30s_tmax_06)/4) %>%
  mutate(wc2.1_30s_tmin_spring = (wc2.1_30s_tmin_03 + wc2.1_30s_tmin_04 + wc2.1_30s_tmin_05 + wc2.1_30s_tmin_06)/4) %>%
  mutate(wc2.1_30s_prec_winter = (wc2.1_30s_prec_01 + wc2.1_30s_prec_02 + wc2.1_30s_prec_11 + wc2.1_30s_prec_12)/4) %>%
  mutate(wc2.1_30s_srad_winter = (wc2.1_30s_srad_01 + wc2.1_30s_srad_02 + wc2.1_30s_srad_11 + wc2.1_30s_srad_12)/4) %>%
  mutate(wc2.1_30s_tavg_winter = (wc2.1_30s_tavg_01 + wc2.1_30s_tavg_02 + wc2.1_30s_tavg_11 + wc2.1_30s_tavg_12)/4) %>%
  mutate(wc2.1_30s_tmax_winter = (wc2.1_30s_tmax_01 + wc2.1_30s_tmax_02 + wc2.1_30s_tmax_11 + wc2.1_30s_tmax_12)/4) %>%
  mutate(wc2.1_30s_tmin_winter = (wc2.1_30s_tmin_01 + wc2.1_30s_tmin_02 + wc2.1_30s_tmin_11 + wc2.1_30s_tmin_12)/4)

write.csv(dftest, "/bioclim/wc0.5/bdis332_bioclim.csv")

