###############################################################################
# Project: Causal black carbon on birth weight in MA                          #
# Code: get geoinfo from latitude and longitude                               #
# Input: mabirths_02NOV18.csv (from Anna)                                     #
# Output: "birth_geocode`year`.csv" geocodes for whole MA birth data          #
# Output: "birth_geocode0010.csv" for whole MA birth data                     #
# Output: "birth_geocode
# Author: Shuxin Dong                                                         #
# Date: Oct 21, 2020                                                          #
###############################################################################

############################# 0. Setup ########################################
rm(list = ls())
gc()

library(data.table)
library(doParallel)
library(foreach)
n_cores <- detectCores()-1

setwd("/media/gate/Shuxin/")
dir_input <- "/media/gate/Shuxin/MAbirth/"
dir_output <- "/media/gate/Shuxin/MAbirth/MAbirth_SES/"

## load data
birth <- fread(file = paste0(dir_input, "mabirths_02NOV18.csv"))
birth <- birth[,.(uniqueid_yr, year, long, lat)]
birth <- na.omit(birth)
dim(birth)
# > dim(birth)
# [1] 1174107       4
# birth <- birth[,.SD[sample(.N, min(.N, floor(0.0001*dim(birth)[1])))]] # 0.01% sample

############################# 1. geocode function ##############################
## benchmark as current to get 2010 census tract
## benchmark as 2010 to get 2000 census tract
## the 'get_geoid_latlon' function is to get the corresponding census tract
get_geoid_latlon <- function (lat, lon, benchmark, vintage) {
  if (missing(benchmark)) {
    benchmark <- "Public_AR_Current"
  }
  else {
    benchmark <- benchmark
  }
  if (missing(vintage)) {
    vintage <- "Census2010_Current"
  }
  else {
    vintage <- vintage
  }
  if(benchmark == "Public_AR_Current"){
    benchmark <- 4
  }else if(benchmark == "Public_AR_Census2010"){
    benchmark <- 9
  }
  if(vintage == "Census2010_Current"){
    vintage <- 410
  }else if(vintage == "ACS2013_Current"){
    vintage <- 413
  }else if(vintage == "ACS2014_Current"){
    vintage <- 414
  }else if(vintage == "ACS2015_Current"){
    vintage <- 415
  }else if(vintage == "Census2000_Census2010"){
    vintage <- 900
  }
  call_start <- "https://geocoding.geo.census.gov/geocoder/geographies/coordinates?"
  url0 <- paste0("x=", lon, "&y=", lat)
  benchmark0 <- paste0("&benchmark=", benchmark)
  vintage0 <- paste0("&vintage=", vintage, "&format=json")
  url_full <- paste0(call_start, url0, benchmark0, vintage0)
  r <- httr::GET(url_full)
  httr::stop_for_status(r)
  response <- httr::content(r)
  if (length(response$result$geographies$`Census Tracts`[[1]]$GEOID) == 
      0) {
    message(paste0("Lat/lon (", lat, ", ", lon, ") returned no geocodes. An NA was returned."))
    return(NA_character_)
  }
  else {
    if (length(response$result$geographies$`Census Tracts`[[1]]$GEOID) > 
        1) {
      message(paste0("Lat/lon (", lat, ", ", lon, ") returned more than geocode. The first match was returned."))
    }
    return(response$result$geographies$`Census Tracts`[[1]]$GEOID)
  }
}

############################# 2. get geocode ##################################
############# get geocode for 2000 and 2010 for whole MA birth ################
cl <- makeCluster(n_cores, type="FORK") 
registerDoParallel(cl)
## get geocode for 2010
start_time_2010 <- Sys.time()
geoid2010 <- foreach(birth = iter(birth, by = "row"), 
                     .combine = c, .packages = "data.table") %dopar% 
  get_geoid_latlon(lat = birth[["lat"]], lon = birth[['long']],
                   benchmark = "Public_AR_Current",
                   vintage = "Census2010_Current")
end_time_2010 <- Sys.time()
time_2010 <- end_time_2010 - start_time_2010
print(time_2010)
## get geocode for 2000
start_time_2000 <- Sys.time()
geoid2000 <- foreach(birth = iter(birth, by = "row"), 
                     .combine = c, .packages = "data.table") %dopar% 
  get_geoid_latlon(lat = birth[["lat"]], lon = birth[['long']],
                   benchmark = "Public_AR_Census2010",
                   vintage = "Census2000_Census2010")
end_time_2000 <- Sys.time()
time_2000 <- end_time_2000 - start_time_2000
print(time_2000)
birth <- cbind(birth, geoid2000, geoid2010)
stopCluster(cl)

setDT(birth)
fwrite(birth, paste0(dir_output,"birth_geocode0010.csv"))
rm(birth)
gc()
############# only fetch geocodes for those born in certain year ##############
registerDoParallel(cl)
########## get geocode for 2013 for those born in 2013 (year==2013) ###########
birth_2013 <- birth[year==2013]
start_time_2013 <- Sys.time()
geoid2013 <- foreach(birth = iter(birth_2013, by = "row"),
                     .combine = c, .packages = "data.table") %dopar%
  get_geoid_latlon(lat = birth[["lat"]], lon = birth[['long']],
                   benchmark = "Public_AR_Current",
                   vintage = "ACS2013_Current")
end_time_2013 <- Sys.time()
time_2013 <- end_time_2013 - start_time_2013
print(time_2013)
birth_2013 <- cbind(birth_2013, geoid2013)

setDT(birth_2013)
fwrite(birth_2013, paste0(dir_output,"birth_geocode13.csv"))
rm(birth_2013)
gc()

########## get geocode for 2014 for those born in 2014 (year==2014) ###########
birth_2014 <- birth[year==2014]
start_time_2014 <- Sys.time()
geoid2014 <- foreach(birth = iter(birth_2014, by = "row"),
                     .combine = c, .packages = "data.table") %dopar%
  get_geoid_latlon(lat = birth[["lat"]], lon = birth[['long']],
                   benchmark = "Public_AR_Current",
                   vintage = "ACS2014_Current")
end_time_2014 <- Sys.time()
time_2014 <- end_time_2014 - start_time_2014
print(time_2014)
birth_2014 <- cbind(birth_2014, geoid2014)

setDT(birth_2014)
fwrite(birth_2014, paste0(dir_output,"birth_geocode14.csv"))
rm(birth_2014)
gc()

########## get geocode for 2015 for those born in 2015 (year==2015) ###########
birth_2015 <- birth[year==2015]
start_time_2015 <- Sys.time()
geoid2015 <- foreach(birth = iter(birth_2015, by = "row"),
                     .combine = c, .packages = "data.table") %dopar%
  get_geoid_latlon(lat = birth[["lat"]], lon = birth[['long']],
                   benchmark = "Public_AR_Current",
                   vintage = "ACS2015_Current")
end_time_2015 <- Sys.time()
time_2015 <- end_time_2015 - start_time_2015
print(time_2015)
birth_2015 <- cbind(birth_2015, geoid2015)

setDT(birth_2015)
write.csv(birth_2015, file = paste0(dir_output,"birth_geocode15.csv"))
rm(birth_2015)
gc()

stopCluster(cl)
