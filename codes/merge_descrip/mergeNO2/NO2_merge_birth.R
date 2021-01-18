###############################################################################
# Project: Black Carbon and MA birth                                          #
# Analysis: Link 1km grid NO2                                                 # 
# Machine: 164 RStudio Server                                                 #
# Input: "mabirths_02NOV18.csv" (from Anna)                                   #
# Input: NO2 daily prediction on QNAP3                                        #
# Output: "no2_birth.csv" contains 280 days exposure prior to dob             #
# Output: "no2_birth_average.csv" (30d, 31-90d, 91-280d averages)             #
# Author: Shuxin Dong (modified from Yaguang's codes)                         #
# Date: 2021-01-06                                                            #
###############################################################################

############################# 0. Setup ########################################
rm(list=ls())
gc()

library(data.table)
library(nabor)

setDTthreads(threads = 0)
setwd("/media/gate/Shuxin")

dir_birth <- "/media/gate/Shuxin/"
dir_no2 <- "/media/qnap3/Exposure modeling/3 National scale/USA/5 NO2 v1.2000_16/7 Final predictions 1 km/Daily/"
dir_output <- "/media/gate/Shuxin/"

## read in birth data and convert formats
birth <- fread(paste0(dir_birth, "mabirths_02NOV18.csv"))
names(birth)
dat <- birth[,.(uniqueid_yr, bdob, long, lat)]
rm(birth)
dat <- na.omit(dat)
dat$id <- 1:nrow(dat)
dat$bdob <- as.Date(dat$bdob, tryFormats = "%m/%d/%Y")

## split raw birth data for memory saving
dim(dat)[1]/2
dat_1 <- dat[1:floor(dim(dat)[1]/2),]
dat_2 <- dat[ceiling(dim(dat)[1]/2):dim(dat)[1], ]

gc()
##################### 1. extract NO2 data files ###############################
### extract exposures
no2_files_2000 <- list.files(path = dir_no2, pattern = "^PredictionStep2_NO2_USGrid_200002(.*)rds$")
no2_files_2001 <- list.files(path = dir_no2, pattern = "^PredictionStep2_NO2_USGrid_2001(.*)rds$")
no2_files_2002 <- list.files(path = dir_no2, pattern = "^PredictionStep2_NO2_USGrid_2002(.*)rds$")
no2_files_2003 <- list.files(path = dir_no2, pattern = "^PredictionStep2_NO2_USGrid_2003(.*)rds$")
no2_files_2004 <- list.files(path = dir_no2, pattern = "^PredictionStep2_NO2_USGrid_2004(.*)rds$")
no2_files_2005 <- list.files(path = dir_no2, pattern = "^PredictionStep2_NO2_USGrid_2005(.*)rds$")
no2_files_2006 <- list.files(path = dir_no2, pattern = "^PredictionStep2_NO2_USGrid_2006(.*)rds$")
no2_files_2007 <- list.files(path = dir_no2, pattern = "^PredictionStep2_NO2_USGrid_2007(.*)rds$")
no2_files_2008 <- list.files(path = dir_no2, pattern = "^PredictionStep2_NO2_USGrid_2008(.*)rds$")
no2_files_2009 <- list.files(path = dir_no2, pattern = "^PredictionStep2_NO2_USGrid_2009(.*)rds$")
no2_files_2010 <- list.files(path = dir_no2, pattern = "^PredictionStep2_NO2_USGrid_2010(.*)rds$")
no2_files_2011 <- list.files(path = dir_no2, pattern = "^PredictionStep2_NO2_USGrid_2011(.*)rds$")
no2_files_2012 <- list.files(path = dir_no2, pattern = "^PredictionStep2_NO2_USGrid_2012(.*)rds$")
no2_files_2013 <- list.files(path = dir_no2, pattern = "^PredictionStep2_NO2_USGrid_2013(.*)rds$")
no2_files_2014 <- list.files(path = dir_no2, pattern = "^PredictionStep2_NO2_USGrid_2014(.*)rds$")
no2_files_2015 <- list.files(path = dir_no2, pattern = "^PredictionStep2_NO2_USGrid_2015(.*)rds$")
no2_files <- c(no2_files_2000, no2_files_2001, no2_files_2002, no2_files_2003, no2_files_2004, no2_files_2005,
               no2_files_2006, no2_files_2007, no2_files_2008, no2_files_2009, no2_files_2010,
               no2_files_2011, no2_files_2012, no2_files_2013, no2_files_2014, no2_files_2015)

##################### 2. link NO2 to birth data ###############################
sitecode_no2 <- readRDS(paste0(dir_no2, "USGridSite.rds"))

########################## 1st part of birth data #############################
## identify grid id
link_no2 <- nabor::knn(sitecode_no2[, c("Lon","Lat")],
                       dat_1[, c("long","lat")],
                       k = 1, radius = 2*sqrt(2)/10*0.1)
dat_1 <- dat_1[, grid_id_no2 := link_no2$nn.idx]
## merge
no2 <- data.table(matrix(vector(),0,6))
names(no2) <- c("lon", "lat", "id", "bdob", "no2", "date")
no2$bdob <- as.Date(no2$bdob)
no2$date <- as.Date(no2$date)
for (i in 1:length(no2_files)) { # 
  no2_temp <- readRDS(paste0(dir_no2,no2_files[i]))
  no2_temp <- no2_temp[1, dat_1$grid_id_no2]
  
  no2_temp <- data.table(cbind(dat_1[, .(long, lat, id, bdob)], no2_temp))
  no2_temp$date <- as.Date(i + as.Date('2000-02-01') - 1)
  names(no2_temp) <- c("lon","lat", "id", "bdob", "no2", "date")
  no2_temp <- no2_temp[(date >= bdob-280)&(date < bdob)] # restrict to 280 days prior to dob
  no2 <- rbind(no2, no2_temp)
  
  rm(no2_temp)
  gc()
  
  if (i%%10 == 0) cat(paste0('Processed ', i, ' of ', length(no2_files),'\n'))
}
## save long format
no2_dat_1 <- no2
setDT(no2_dat_1)
no2_dat_1 <- no2_dat_1[order(id,-date)]
fwrite(no2_dat_1, paste0(dir_output, "no2_dat_1_long.csv"))
## transform from long format to wide format
no2_dat_1[, pre_day:= bdob - date][]
no2_dat_1_wide <- dcast(no2_dat_1, id ~ pre_day, value.var = "no2") # may change no2 var name
no2_dat_1_wide <- merge(no2_dat_1_wide, dat[,.(uniqueid_yr,id)], by = "id")
## save wide format part 1
# fwrite(no2_dat_1_wide, paste0(dir_output, "no2_dat_1_wide.csv"))

########################## 2nd part of birth data #############################
## identify grid id
link_no2 <- nabor::knn(sitecode_no2[, c("Lon","Lat")],
                       dat_2[, c("long","lat")],
                       k = 1, radius = 2*sqrt(2)/10*0.1)
dat_2 <- dat_2[, grid_id_no2 := link_no2$nn.idx]
## merge
no2 <- data.table(matrix(vector(),0,6))
names(no2) <- c("lon", "lat", "id", "bdob", "no2", "date")
no2$bdob <- as.Date(no2$bdob)
no2$date <- as.Date(no2$date)
for (i in 1:length(no2_files)) { # 
  no2_temp <- readRDS(paste0(dir_no2,no2_files[i]))
  no2_temp <- no2_temp[1, dat_2$grid_id_no2]
  
  no2_temp <- data.table(cbind(dat_2[, .(long, lat, id, bdob)], no2_temp))
  no2_temp$date <- as.Date(i + as.Date('2000-02-01') - 1)
  names(no2_temp) <- c("lon","lat", "id", "bdob", "no2", "date")
  no2_temp <- no2_temp[(date >= bdob-280)&(date < bdob)] # restrict to 280 days prior to dob
  no2 <- rbind(no2, no2_temp)
  
  rm(no2_temp)
  gc()
  
  if (i%%10 == 0) cat(paste0('Processed ', i, ' of ', length(no2_files),'\n'))
}
## save long format
no2_dat_2 <- no2
setDT(no2_dat_2)
no2_dat_2 <- no2_dat_2[order(id,-date)]
fwrite(no2_dat_2, paste0(dir_output, "no2_dat_2_long.csv"))
## transform from long format to wide format
no2_dat_2[, pre_day:= bdob - date][]
no2_dat_2_wide <- dcast(no2_dat_2, id ~ pre_day, value.var = "no2") # may change no2 var name
no2_dat_2_wide
no2_dat_2_wide <- merge(no2_dat_2_wide, dat[,.(uniqueid_yr,id)], by = "id")
## save wide format part 1
# fwrite(no2_dat_2_wide, paste0(dir_output, "no2_dat_2_wide.csv"))

##################### 3. generate the final NO2 dataset ######################
no2_birth <- rbind(no2_dat_1_wide, no2_dat_2_wide)
setDT(no2_birth)
fwrite(no2_birth, paste0(dir_output, "no2_birth.csv"))

no2_birth[, ":=" (no2_30d = rowMeans(no2_birth[,2:31], na.rm = TRUE),
                  no2_3090d = rowMeans(no2_birth[,32:91], na.rm = TRUE),
                  no2_90280d = rowMeans(no2_birth[,92:281], na.rm = TRUE))][]
no2_birth_average <- na.omit(no2_birth[, 282:285])
fwrite(no2_birth_average, paste0(dir_output, "no2_birth_average.csv"))
