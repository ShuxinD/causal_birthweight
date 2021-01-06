###############################################################################
# Project: Black Carbon and MA birth                                          #
# Analysis: Link 1km grid NO2                                                 # 
# Machine: RStudio Server                                                     #
###############################################################################

######################### 0. Setup #########################
rm(list=ls())
gc()

require(dplyr)
require(magrittr)
require(lubridate)
require(data.table)
library(nabor)
library(fst)

dir_dat <- '/media/qnap3/Yaguang/Abbie/'
dir_pm <- '/media/qnap3/Exposure modeling/3 National scale/USA/4 PM2.5 v2.2000_16/7 Final predictions 1 km/Daily/'
dir_o3 <- '/media/qnap3/Exposure modeling/3 National scale/USA/6 O3 vx.2000_16/7 Final predictions 1 km/Daily/'
dir_no2 <- '/media/qnap3/Exposure modeling/3 National scale/USA/5 NO2 v1.2000_16/7 Final predictions 1 km/Daily/'
dir_temp <- '/media/qnap3/Met_Pred/12km_met_data/TMP/'
dir_humid <- '/media/qnap3/Met_Pred/12km_met_data/SPFH/'
dir_save_GA <- '/media/qnap3/Yaguang/Abbie/GA/'

dat <- read.csv(paste0(dir_dat,'georgia_unique.csv'))
dat$id <- 1:nrow(dat)



######################### 1. pm ####################
### identify grid id
sitecode_pm <- readRDS(paste0(dir_pm,'USGridSite.rds'))
link_pm <- nabor::knn(sitecode_pm[,c("Lon","Lat")],dat[,c("LONG","LAT")],k=1,radius=2*sqrt(2)/10*0.1)
dat <- cbind.data.frame(dat,link_pm$nn.idx)
names(dat)[4] <- "grid_id_pm"

### extract exposures
pm_files_200612 <- list.files(path=dir_pm,pattern = "^PredictionStep2_PM25_USGrid_200612(.*)rds$")  # we need lags so I extract those lagged up to a month
pm_files_2007 <- list.files(path=dir_pm,pattern = "^PredictionStep2_PM25_USGrid_2007(.*)rds$")
pm_files_2008 <- list.files(path=dir_pm,pattern = "^PredictionStep2_PM25_USGrid_2008(.*)rds$")
pm_files_2009 <- list.files(path=dir_pm,pattern = "^PredictionStep2_PM25_USGrid_2009(.*)rds$")
pm_files_2010 <- list.files(path=dir_pm,pattern = "^PredictionStep2_PM25_USGrid_2010(.*)rds$")
pm_files_2011 <- list.files(path=dir_pm,pattern = "^PredictionStep2_PM25_USGrid_2011(.*)rds$")
pm_files_201201 <- list.files(path=dir_pm,pattern = "^PredictionStep2_PM25_USGrid_201201(.*)rds$")  # I also extract those leaded up to a month
pm_files <- c(pm_files_200612,pm_files_2007,pm_files_2008,pm_files_2009,pm_files_2010,pm_files_2011,pm_files_201201)

### merge
pm <- data.frame(matrix(vector(),0,5),stringsAsFactors = FALSE)
names(pm) <- c("lon","lat",'id','pm25','date')
for (i in 1:length(pm_files)) {
  pm_temp <- readRDS(paste0(dir_pm,pm_files[i]))
  pm_temp <- pm_temp[1,dat$grid_id_pm]
  
  pm_temp <- data.frame(cbind(dat[,c("LONG","LAT","id")],pm_temp))
  pm_temp$date <- as.Date(i+as.Date('2006-12-01')-1)
  names(pm_temp) <- c("lon","lat",'id','pm25','date')
  pm <- rbind(pm,pm_temp)
  
  rm(pm_temp)
  gc()
  
  if (i%%10 == 0) cat(paste0('Processed ', i, ' of ', length(pm_files),'\n'))
}

### lagging and leading
pm <- as.data.frame(pm)
pm <- pm[order(pm$id,pm$date),]
pm <- as.data.table(pm)

pm = pm[, pm25_lag1 := shift(pm25,type='lag'),.(id)]
pm = pm[, pm25_lag2 := shift(pm25_lag1,type='lag'),.(id)]
pm = pm[, pm25_lag3 := shift(pm25_lag2,type='lag'),.(id)]
pm = pm[, pm25_lag4 := shift(pm25_lag3,type='lag'),.(id)]
pm = pm[, pm25_lag5 := shift(pm25_lag4,type='lag'),.(id)]
pm = pm[, pm25_lag6 := shift(pm25_lag5,type='lag'),.(id)]
pm = pm[, pm25_lag7 := shift(pm25_lag6,type='lag'),.(id)]
pm = pm[, pm25_lag8 := shift(pm25_lag7,type='lag'),.(id)]
pm = pm[, pm25_lag9 := shift(pm25_lag8,type='lag'),.(id)]
pm = pm[, pm25_lag10 := shift(pm25_lag9,type='lag'),.(id)]
pm = pm[, pm25_lag11 := shift(pm25_lag10,type='lag'),.(id)]
pm = pm[, pm25_lag12 := shift(pm25_lag11,type='lag'),.(id)]
pm = pm[, pm25_lag13 := shift(pm25_lag12,type='lag'),.(id)]
pm = pm[, pm25_lag14 := shift(pm25_lag13,type='lag'),.(id)]
pm = pm[, pm25_lag15 := shift(pm25_lag14,type='lag'),.(id)]
pm = pm[, pm25_lag16 := shift(pm25_lag15,type='lag'),.(id)]
pm = pm[, pm25_lag17 := shift(pm25_lag16,type='lag'),.(id)]
pm = pm[, pm25_lag18 := shift(pm25_lag17,type='lag'),.(id)]
pm = pm[, pm25_lag19 := shift(pm25_lag18,type='lag'),.(id)]
pm = pm[, pm25_lag20 := shift(pm25_lag19,type='lag'),.(id)]

pm$pm25_lag713 <- (pm$pm25_lag7+pm$pm25_lag8+pm$pm25_lag9+pm$pm25_lag10+pm$pm25_lag11+
                     pm$pm25_lag12+pm$pm25_lag13)/7
pm$pm25_lag1420 <- (pm$pm25_lag14+pm$pm25_lag15+pm$pm25_lag16+pm$pm25_lag17+pm$pm25_lag18+
                      pm$pm25_lag19+pm$pm25_lag20)/7

pm = pm[, pm25_lead1 := shift(pm25,type='lead'),.(id)]
pm = pm[, pm25_lead2 := shift(pm25_lead1,type='lead'),.(id)]
pm = pm[, pm25_lead3 := shift(pm25_lead2,type='lead'),.(id)]
pm = pm[, pm25_lead4 := shift(pm25_lead3,type='lead'),.(id)]
pm = pm[, pm25_lead5 := shift(pm25_lead4,type='lead'),.(id)]
pm = pm[, pm25_lead6 := shift(pm25_lead5,type='lead'),.(id)]
pm = pm[, pm25_lead7 := shift(pm25_lead6,type='lead'),.(id)]
pm = pm[, pm25_lead8 := shift(pm25_lead7,type='lead'),.(id)]
pm = pm[, pm25_lead9 := shift(pm25_lead8,type='lead'),.(id)]
pm = pm[, pm25_lead10 := shift(pm25_lead9,type='lead'),.(id)]
pm = pm[, pm25_lead11 := shift(pm25_lead10,type='lead'),.(id)]
pm = pm[, pm25_lead12 := shift(pm25_lead11,type='lead'),.(id)]
pm = pm[, pm25_lead13 := shift(pm25_lead12,type='lead'),.(id)]
pm = pm[, pm25_lead14 := shift(pm25_lead13,type='lead'),.(id)]
pm = pm[, pm25_lead15 := shift(pm25_lead14,type='lead'),.(id)]
pm = pm[, pm25_lead16 := shift(pm25_lead15,type='lead'),.(id)]
pm = pm[, pm25_lead17 := shift(pm25_lead16,type='lead'),.(id)]
pm = pm[, pm25_lead18 := shift(pm25_lead17,type='lead'),.(id)]
pm = pm[, pm25_lead19 := shift(pm25_lead18,type='lead'),.(id)]
pm = pm[, pm25_lead20 := shift(pm25_lead19,type='lead'),.(id)]

pm$pm25_lead713 <- (pm$pm25_lead7+pm$pm25_lead8+pm$pm25_lead9+pm$pm25_lead10+pm$pm25_lead11+
                      pm$pm25_lead12+pm$pm25_lead13)/7
pm$pm25_lead1420 <- (pm$pm25_lead14+pm$pm25_lead15+pm$pm25_lead16+pm$pm25_lead17+pm$pm25_lead18+
                       pm$pm25_lead19+pm$pm25_lead20)/7

pm <- pm[,c("lon","lat","date","pm25","pm25_lag1","pm25_lag2","pm25_lag3","pm25_lag4","pm25_lag5","pm25_lag6","pm25_lag713",
            "pm25_lag1420","pm25_lead1","pm25_lead2","pm25_lead3","pm25_lead4","pm25_lead5","pm25_lead6","pm25_lead713",
            "pm25_lead1420")]

### save
saveRDS(pm,file=paste0(dir_save_GA,'pm_GA.rds'))



######################### 2. o3 ####################
### identify grid id
sitecode_o3 <- readRDS(paste0(dir_o3,'USGridSite_O3.rds'))
link_o3 <- nabor::knn(sitecode_o3[,c("Lon","Lat")],dat[,c("LONG","LAT")],k=1,radius=2*sqrt(2)/10*0.1)
dat <- cbind.data.frame(dat,link_o3$nn.idx)
names(dat)[4] <- "grid_id_o3"

### extract exposures
o3_files_200612 <- list.files(path=dir_o3,pattern = "^Ensemble_predictions_O3_USGrid_200612(.*)rds$")
o3_files_2007 <- list.files(path=dir_o3,pattern = "^Ensemble_predictions_O3_USGrid_2007(.*)rds$")
o3_files_2008 <- list.files(path=dir_o3,pattern = "^Ensemble_predictions_O3_USGrid_2008(.*)rds$")
o3_files_2009 <- list.files(path=dir_o3,pattern = "^Ensemble_predictions_O3_USGrid_2009(.*)rds$")
o3_files_2010 <- list.files(path=dir_o3,pattern = "^Ensemble_predictions_O3_USGrid_2010(.*)rds$")
o3_files_2011 <- list.files(path=dir_o3,pattern = "^Ensemble_predictions_O3_USGrid_2011(.*)rds$")
o3_files_201201 <- list.files(path=dir_o3,pattern = "^Ensemble_predictions_O3_USGrid_201201(.*)rds$")
o3_files <- c(o3_files_200612,o3_files_2007,o3_files_2008,o3_files_2009,o3_files_2010,o3_files_2011,o3_files_201201)

### merge
o3 <- data.frame(matrix(vector(),0,5),stringsAsFactors = FALSE)
names(o3) <- c("lon","lat",'id','o3','date')
for (i in 1:length(o3_files)) {
  o3_temp <- readRDS(paste0(dir_o3,o3_files[i]))
  o3_temp <- o3_temp[dat$grid_id_o3]
  
  o3_temp <- data.frame(cbind(dat[,c("LONG","LAT","id")],o3_temp))
  o3_temp$date <- as.Date(i+as.Date('2006-12-01')-1)
  names(o3_temp) <- c("lon","lat",'id','o3','date')
  o3 <- rbind(o3,o3_temp)
  
  rm(o3_temp)
  gc()
  
  if (i%%10 == 0) cat(paste0('Processed ', i, ' of ', length(o3_files),'\n'))
}

### lagging and leading
o3 <- as.data.frame(o3)
o3 <- o3[order(o3$id,o3$date),]
o3 <- as.data.table(o3)

o3 = o3[, o3_lag1 := shift(o3,type='lag'),.(id)]
o3 = o3[, o3_lag2 := shift(o3_lag1,type='lag'),.(id)]
o3 = o3[, o3_lag3 := shift(o3_lag2,type='lag'),.(id)]
o3 = o3[, o3_lag4 := shift(o3_lag3,type='lag'),.(id)]
o3 = o3[, o3_lag5 := shift(o3_lag4,type='lag'),.(id)]
o3 = o3[, o3_lag6 := shift(o3_lag5,type='lag'),.(id)]
o3 = o3[, o3_lag7 := shift(o3_lag6,type='lag'),.(id)]
o3 = o3[, o3_lag8 := shift(o3_lag7,type='lag'),.(id)]
o3 = o3[, o3_lag9 := shift(o3_lag8,type='lag'),.(id)]
o3 = o3[, o3_lag10 := shift(o3_lag9,type='lag'),.(id)]
o3 = o3[, o3_lag11 := shift(o3_lag10,type='lag'),.(id)]
o3 = o3[, o3_lag12 := shift(o3_lag11,type='lag'),.(id)]
o3 = o3[, o3_lag13 := shift(o3_lag12,type='lag'),.(id)]
o3 = o3[, o3_lag14 := shift(o3_lag13,type='lag'),.(id)]
o3 = o3[, o3_lag15 := shift(o3_lag14,type='lag'),.(id)]
o3 = o3[, o3_lag16 := shift(o3_lag15,type='lag'),.(id)]
o3 = o3[, o3_lag17 := shift(o3_lag16,type='lag'),.(id)]
o3 = o3[, o3_lag18 := shift(o3_lag17,type='lag'),.(id)]
o3 = o3[, o3_lag19 := shift(o3_lag18,type='lag'),.(id)]
o3 = o3[, o3_lag20 := shift(o3_lag19,type='lag'),.(id)]

o3$o3_lag713 <- (o3$o3_lag7+o3$o3_lag8+o3$o3_lag9+o3$o3_lag10+o3$o3_lag11+
                   o3$o3_lag12+o3$o3_lag13)/7
o3$o3_lag1420 <- (o3$o3_lag14+o3$o3_lag15+o3$o3_lag16+o3$o3_lag17+o3$o3_lag18+
                    o3$o3_lag19+o3$o3_lag20)/7

o3 = o3[, o3_lead1 := shift(o3,type='lead'),.(id)]
o3 = o3[, o3_lead2 := shift(o3_lead1,type='lead'),.(id)]
o3 = o3[, o3_lead3 := shift(o3_lead2,type='lead'),.(id)]
o3 = o3[, o3_lead4 := shift(o3_lead3,type='lead'),.(id)]
o3 = o3[, o3_lead5 := shift(o3_lead4,type='lead'),.(id)]
o3 = o3[, o3_lead6 := shift(o3_lead5,type='lead'),.(id)]
o3 = o3[, o3_lead7 := shift(o3_lead6,type='lead'),.(id)]
o3 = o3[, o3_lead8 := shift(o3_lead7,type='lead'),.(id)]
o3 = o3[, o3_lead9 := shift(o3_lead8,type='lead'),.(id)]
o3 = o3[, o3_lead10 := shift(o3_lead9,type='lead'),.(id)]
o3 = o3[, o3_lead11 := shift(o3_lead10,type='lead'),.(id)]
o3 = o3[, o3_lead12 := shift(o3_lead11,type='lead'),.(id)]
o3 = o3[, o3_lead13 := shift(o3_lead12,type='lead'),.(id)]
o3 = o3[, o3_lead14 := shift(o3_lead13,type='lead'),.(id)]
o3 = o3[, o3_lead15 := shift(o3_lead14,type='lead'),.(id)]
o3 = o3[, o3_lead16 := shift(o3_lead15,type='lead'),.(id)]
o3 = o3[, o3_lead17 := shift(o3_lead16,type='lead'),.(id)]
o3 = o3[, o3_lead18 := shift(o3_lead17,type='lead'),.(id)]
o3 = o3[, o3_lead19 := shift(o3_lead18,type='lead'),.(id)]
o3 = o3[, o3_lead20 := shift(o3_lead19,type='lead'),.(id)]

o3$o3_lead713 <- (o3$o3_lead7+o3$o3_lead8+o3$o3_lead9+o3$o3_lead10+o3$o3_lead11+
                    o3$o3_lead12+o3$o3_lead13)/7
o3$o3_lead1420 <- (o3$o3_lead14+o3$o3_lead15+o3$o3_lead16+o3$o3_lead17+o3$o3_lead18+
                     o3$o3_lead19+o3$o3_lead20)/7

o3 <- o3[,c("lon","lat","date","o3","o3_lag1","o3_lag2","o3_lag3","o3_lag4","o3_lag5","o3_lag6","o3_lag713",
            "o3_lag1420","o3_lead1","o3_lead2","o3_lead3","o3_lead4","o3_lead5","o3_lead6","o3_lead713",
            "o3_lead1420")]

### save
saveRDS(o3,file=paste0(dir_save_GA,'o3_GA.rds'))



######################### 3. no2 ####################
### identify grid id
sitecode_no2 <- readRDS(paste0(dir_no2,'USGridSite.rds'))
link_no2 <- nabor::knn(sitecode_no2[,c("Lon","Lat")],dat[,c("LONG","LAT")],k=1,radius=2*sqrt(2)/10*0.1)
dat <- cbind.data.frame(dat,link_no2$nn.idx)
names(dat)[4] <- "grid_id_no2"

### extract exposures
no2_files_200612 <- list.files(path=dir_no2,pattern = "^PredictionStep2_NO2_USGrid_200612(.*)rds$")
no2_files_2007 <- list.files(path=dir_no2,pattern = "^PredictionStep2_NO2_USGrid_2007(.*)rds$")
no2_files_2008 <- list.files(path=dir_no2,pattern = "^PredictionStep2_NO2_USGrid_2008(.*)rds$")
no2_files_2009 <- list.files(path=dir_no2,pattern = "^PredictionStep2_NO2_USGrid_2009(.*)rds$")
no2_files_2010 <- list.files(path=dir_no2,pattern = "^PredictionStep2_NO2_USGrid_2010(.*)rds$")
no2_files_2011 <- list.files(path=dir_no2,pattern = "^PredictionStep2_NO2_USGrid_2011(.*)rds$")
no2_files_201201 <- list.files(path=dir_no2,pattern = "^PredictionStep2_NO2_USGrid_201201(.*)rds$")
no2_files <- c(no2_files_200612,no2_files_2007,no2_files_2008,no2_files_2009,no2_files_2010,no2_files_2011,no2_files_201201)

### merge
no2 <- data.frame(matrix(vector(),0,5),stringsAsFactors = FALSE)
names(no2) <- c("lon","lat",'id','no2','date')
for (i in 1:length(no2_files)) {
  no2_temp <- readRDS(paste0(dir_no2,no2_files[i]))
  no2_temp <- no2_temp[1,dat$grid_id_no2]
  
  no2_temp <- data.frame(cbind(dat[,c("LONG","LAT","id")],no2_temp))
  no2_temp$date <- as.Date(i+as.Date('2006-12-01')-1)
  names(no2_temp) <- c("lon","lat",'id','no2','date')
  no2 <- rbind(no2,no2_temp)
  
  rm(no2_temp)
  gc()
  
  if (i%%10 == 0) cat(paste0('Processed ', i, ' of ', length(no2_files),'\n'))
}

### lagging and leading
no2 <- as.data.frame(no2)
no2 <- no2[order(no2$id,no2$date),]
no2 <- as.data.table(no2)

no2 = no2[, no2_lag1 := shift(no2,type='lag'),.(id)]
no2 = no2[, no2_lag2 := shift(no2_lag1,type='lag'),.(id)]
no2 = no2[, no2_lag3 := shift(no2_lag2,type='lag'),.(id)]
no2 = no2[, no2_lag4 := shift(no2_lag3,type='lag'),.(id)]
no2 = no2[, no2_lag5 := shift(no2_lag4,type='lag'),.(id)]
no2 = no2[, no2_lag6 := shift(no2_lag5,type='lag'),.(id)]
no2 = no2[, no2_lag7 := shift(no2_lag6,type='lag'),.(id)]
no2 = no2[, no2_lag8 := shift(no2_lag7,type='lag'),.(id)]
no2 = no2[, no2_lag9 := shift(no2_lag8,type='lag'),.(id)]
no2 = no2[, no2_lag10 := shift(no2_lag9,type='lag'),.(id)]
no2 = no2[, no2_lag11 := shift(no2_lag10,type='lag'),.(id)]
no2 = no2[, no2_lag12 := shift(no2_lag11,type='lag'),.(id)]
no2 = no2[, no2_lag13 := shift(no2_lag12,type='lag'),.(id)]
no2 = no2[, no2_lag14 := shift(no2_lag13,type='lag'),.(id)]
no2 = no2[, no2_lag15 := shift(no2_lag14,type='lag'),.(id)]
no2 = no2[, no2_lag16 := shift(no2_lag15,type='lag'),.(id)]
no2 = no2[, no2_lag17 := shift(no2_lag16,type='lag'),.(id)]
no2 = no2[, no2_lag18 := shift(no2_lag17,type='lag'),.(id)]
no2 = no2[, no2_lag19 := shift(no2_lag18,type='lag'),.(id)]
no2 = no2[, no2_lag20 := shift(no2_lag19,type='lag'),.(id)]

no2$no2_lag713 <- (no2$no2_lag7+no2$no2_lag8+no2$no2_lag9+no2$no2_lag10+no2$no2_lag11+
                     no2$no2_lag12+no2$no2_lag13)/7
no2$no2_lag1420 <- (no2$no2_lag14+no2$no2_lag15+no2$no2_lag16+no2$no2_lag17+no2$no2_lag18+
                      no2$no2_lag19+no2$no2_lag20)/7

no2 = no2[, no2_lead1 := shift(no2,type='lead'),.(id)]
no2 = no2[, no2_lead2 := shift(no2_lead1,type='lead'),.(id)]
no2 = no2[, no2_lead3 := shift(no2_lead2,type='lead'),.(id)]
no2 = no2[, no2_lead4 := shift(no2_lead3,type='lead'),.(id)]
no2 = no2[, no2_lead5 := shift(no2_lead4,type='lead'),.(id)]
no2 = no2[, no2_lead6 := shift(no2_lead5,type='lead'),.(id)]
no2 = no2[, no2_lead7 := shift(no2_lead6,type='lead'),.(id)]
no2 = no2[, no2_lead8 := shift(no2_lead7,type='lead'),.(id)]
no2 = no2[, no2_lead9 := shift(no2_lead8,type='lead'),.(id)]
no2 = no2[, no2_lead10 := shift(no2_lead9,type='lead'),.(id)]
no2 = no2[, no2_lead11 := shift(no2_lead10,type='lead'),.(id)]
no2 = no2[, no2_lead12 := shift(no2_lead11,type='lead'),.(id)]
no2 = no2[, no2_lead13 := shift(no2_lead12,type='lead'),.(id)]
no2 = no2[, no2_lead14 := shift(no2_lead13,type='lead'),.(id)]
no2 = no2[, no2_lead15 := shift(no2_lead14,type='lead'),.(id)]
no2 = no2[, no2_lead16 := shift(no2_lead15,type='lead'),.(id)]
no2 = no2[, no2_lead17 := shift(no2_lead16,type='lead'),.(id)]
no2 = no2[, no2_lead18 := shift(no2_lead17,type='lead'),.(id)]
no2 = no2[, no2_lead19 := shift(no2_lead18,type='lead'),.(id)]
no2 = no2[, no2_lead20 := shift(no2_lead19,type='lead'),.(id)]

no2$no2_lead713 <- (no2$no2_lead7+no2$no2_lead8+no2$no2_lead9+no2$no2_lead10+no2$no2_lead11+
                      no2$no2_lead12+no2$no2_lead13)/7
no2$no2_lead1420 <- (no2$no2_lead14+no2$no2_lead15+no2$no2_lead16+no2$no2_lead17+no2$no2_lead18+
                       no2$no2_lead19+no2$no2_lead20)/7

no2 <- no2[,c("lon","lat","date","no2","no2_lag1","no2_lag2","no2_lag3","no2_lag4","no2_lag5","no2_lag6","no2_lag713",
              "no2_lag1420","no2_lead1","no2_lead2","no2_lead3","no2_lead4","no2_lead5","no2_lead6","no2_lead713",
              "no2_lead1420")]

### save
saveRDS(no2,file=paste0(dir_save_GA,'no2_GA.rds'))



######################### 4. temperature ####################
# the temp_all data consists of 55403 grids (with 12km resolution) repeated for each day. So I take the first 55403 rows
# to identify the grids linking to the health data, loop to extract those identified grids for each day, and bind them. 

temp_2006 <- read_fst(paste0(dir_temp,"2006_TMP_summary.fst"))
temp_2007 <- read_fst(paste0(dir_temp,"2007_TMP_summary.fst"))
temp_2008 <- read_fst(paste0(dir_temp,"2008_TMP_summary.fst"))
temp_2009 <- read_fst(paste0(dir_temp,"2009_TMP_summary.fst"))
temp_2010 <- read_fst(paste0(dir_temp,"2010_TMP_summary.fst"))
temp_2011 <- read_fst(paste0(dir_temp,"2011_TMP_summary.fst"))
temp_2012 <- read_fst(paste0(dir_temp,"2012_TMP_summary.fst"))
temp_all <- bind_rows(temp_2006,temp_2007,temp_2008,temp_2009,temp_2010,temp_2011,temp_2012)

temp_all$Date_final2 <- as.Date(temp_all$Date_final2)

#CHANGE
temp_all <- temp_all[temp_all$Date_final2>=as.Date('2006-12-01')&temp_all$Date_final2<=as.Date('2012-01-31'),c('x','y','Date_final2','mean')]
names(temp_all) <- c("lon","lat",'date','temp')

### identify grid id
link_temp <- nabor::knn(temp_all[1:55403,c("lon","lat")],dat[,c("LONG","LAT")],k=1)
dat <- cbind.data.frame(dat,link_temp$nn.idx)
names(dat)[4] <- "grid_id_temp"

#CHANGE
dates <- as.Date('2006-12-01')+c(0:(as.Date('2012-01-31')-as.Date('2006-12-01')))

temp <- data.frame(matrix(vector(),0,5),stringsAsFactors = FALSE)
names(temp) <- c("lon","lat","id",'temp','date')
for (i in 1:length(dates)){
  temp_temp <- temp_all[temp_all$date==dates[i],'temp']
  temp_temp <- temp_temp[dat$grid_id_temp]
  
  temp_temp <- data.frame(cbind(dat[,c("LONG","LAT","id")],temp_temp))
  #CHANGE
  temp_temp$date <- as.Date(i+as.Date('2006-12-01')-1)
  names(temp_temp) <- c("lon","lat",'id','temp','date')
  temp <- rbind(temp,temp_temp)
  
  rm(temp_temp)
  gc()
  
  if (i%%10 == 0) cat(paste0('Processed ', i, ' of ', length(dates),'\n'))
}

temp <- as.data.frame(temp)
temp = temp[order(temp$id,temp$date),]
temp <- as.data.table(temp)

temp = temp[, temp_lag1 := shift(temp,type='lag'),.(id)]
temp = temp[, temp_lag2 := shift(temp_lag1,type='lag'),.(id)]
temp = temp[, temp_lag3 := shift(temp_lag2,type='lag'),.(id)]
temp = temp[, temp_lag4 := shift(temp_lag3,type='lag'),.(id)]
temp = temp[, temp_lag5 := shift(temp_lag4,type='lag'),.(id)]
temp = temp[, temp_lag6 := shift(temp_lag5,type='lag'),.(id)]
temp = temp[, temp_lag7 := shift(temp_lag6,type='lag'),.(id)]
temp = temp[, temp_lag8 := shift(temp_lag7,type='lag'),.(id)]
temp = temp[, temp_lag9 := shift(temp_lag8,type='lag'),.(id)]
temp = temp[, temp_lag10 := shift(temp_lag9,type='lag'),.(id)]
temp = temp[, temp_lag11 := shift(temp_lag10,type='lag'),.(id)]
temp = temp[, temp_lag12 := shift(temp_lag11,type='lag'),.(id)]
temp = temp[, temp_lag13 := shift(temp_lag12,type='lag'),.(id)]
temp = temp[, temp_lag14 := shift(temp_lag13,type='lag'),.(id)]
temp = temp[, temp_lag15 := shift(temp_lag14,type='lag'),.(id)]
temp = temp[, temp_lag16 := shift(temp_lag15,type='lag'),.(id)]
temp = temp[, temp_lag17 := shift(temp_lag16,type='lag'),.(id)]
temp = temp[, temp_lag18 := shift(temp_lag17,type='lag'),.(id)]
temp = temp[, temp_lag19 := shift(temp_lag18,type='lag'),.(id)]
temp = temp[, temp_lag20 := shift(temp_lag19,type='lag'),.(id)]

temp$temp_lag713 <- (temp$temp_lag7+temp$temp_lag8+temp$temp_lag9+temp$temp_lag10+temp$temp_lag11+
                       temp$temp_lag12+temp$temp_lag13)/7
temp$temp_lag1420 <- (temp$temp_lag14+temp$temp_lag15+temp$temp_lag16+temp$temp_lag17+temp$temp_lag18+
                        temp$temp_lag19+temp$temp_lag20)/7

temp = temp[, temp_lead1 := shift(temp,type='lead'),.(id)]
temp = temp[, temp_lead2 := shift(temp_lead1,type='lead'),.(id)]
temp = temp[, temp_lead3 := shift(temp_lead2,type='lead'),.(id)]
temp = temp[, temp_lead4 := shift(temp_lead3,type='lead'),.(id)]
temp = temp[, temp_lead5 := shift(temp_lead4,type='lead'),.(id)]
temp = temp[, temp_lead6 := shift(temp_lead5,type='lead'),.(id)]
temp = temp[, temp_lead7 := shift(temp_lead6,type='lead'),.(id)]
temp = temp[, temp_lead8 := shift(temp_lead7,type='lead'),.(id)]
temp = temp[, temp_lead9 := shift(temp_lead8,type='lead'),.(id)]
temp = temp[, temp_lead10 := shift(temp_lead9,type='lead'),.(id)]
temp = temp[, temp_lead11 := shift(temp_lead10,type='lead'),.(id)]
temp = temp[, temp_lead12 := shift(temp_lead11,type='lead'),.(id)]
temp = temp[, temp_lead13 := shift(temp_lead12,type='lead'),.(id)]
temp = temp[, temp_lead14 := shift(temp_lead13,type='lead'),.(id)]
temp = temp[, temp_lead15 := shift(temp_lead14,type='lead'),.(id)]
temp = temp[, temp_lead16 := shift(temp_lead15,type='lead'),.(id)]
temp = temp[, temp_lead17 := shift(temp_lead16,type='lead'),.(id)]
temp = temp[, temp_lead18 := shift(temp_lead17,type='lead'),.(id)]
temp = temp[, temp_lead19 := shift(temp_lead18,type='lead'),.(id)]
temp = temp[, temp_lead20 := shift(temp_lead19,type='lead'),.(id)]

temp$temp_lead713 <- (temp$temp_lead7+temp$temp_lead8+temp$temp_lead9+temp$temp_lead10+temp$temp_lead11+
                        temp$temp_lead12+temp$temp_lead13)/7
temp$temp_lead1420 <- (temp$temp_lead14+temp$temp_lead15+temp$temp_lead16+temp$temp_lead17+temp$temp_lead18+
                         temp$temp_lead19+temp$temp_lead20)/7

temp <- temp[,c("lon","lat","date","temp","temp_lag1","temp_lag2","temp_lag3","temp_lag4","temp_lag5","temp_lag6","temp_lag713",
                "temp_lag1420","temp_lead1","temp_lead2","temp_lead3","temp_lead4","temp_lead5","temp_lead6","temp_lead713",
                "temp_lead1420")]
saveRDS(temp,file=paste0(dir_save_GA,'temp_GA.rds'))



######################### 5. specific humidity ####################
humid_2006 <- read_fst(paste0(dir_humid,"2006_SPFH_summary.fst"))
humid_2007 <- read_fst(paste0(dir_humid,"2007_SPFH_summary.fst"))
humid_2008 <- read_fst(paste0(dir_humid,"2008_SPFH_summary.fst"))
humid_2009 <- read_fst(paste0(dir_humid,"2009_SPFH_summary.fst"))
humid_2010 <- read_fst(paste0(dir_humid,"2010_SPFH_summary.fst"))
humid_2011 <- read_fst(paste0(dir_humid,"2011_SPFH_summary.fst"))
humid_2012 <- read_fst(paste0(dir_humid,"2012_SPFH_summary.fst"))
humid_all <- bind_rows(humid_2006,humid_2007,humid_2008,humid_2009,humid_2010,humid_2011,humid_2012)
humid_all$Date_final2 <- as.Date(humid_all$Date_final2)
#CHANGE
humid_all <- humid_all[humid_all$Date_final2>=as.Date('2006-12-01')&humid_all$Date_final2<=as.Date('2012-01-31'),c('x','y','Date_final2','mean')]
names(humid_all) <- c("lon","lat",'date','humid')

### identify grid id
link_humid <- nabor::knn(humid_all[1:55403,c("lon","lat")],dat[,c("LONG","LAT")],k=1)
dat <- cbind.data.frame(dat,link_humid$nn.idx)
names(dat)[4] <- "grid_id_humid"

#CHANGE
dates <- as.Date('2006-12-01')+c(0:(as.Date('2012-01-31')-as.Date('2006-12-01')))

humid <- data.frame(matrix(vector(),0,5),stringsAsFactors = FALSE)
names(humid) <- c("lon","lat","id",'humid','date')
for (i in 1:length(dates)){
  humid_temp <- humid_all[humid_all$date==dates[i],'humid']
  humid_temp <- humid_temp[dat$grid_id_humid]
  
  humid_temp <- data.frame(cbind(dat[,c("LONG","LAT","id")],humid_temp))
  #CHANGE
  humid_temp$date <- as.Date(i+as.Date('2006-12-01')-1)
  names(humid_temp) <- c("lon","lat",'id','humid','date')
  humid <- rbind(humid,humid_temp)
  
  rm(humid_temp)
  gc()
  
  if (i%%10 == 0) cat(paste0('Processed ', i, ' of ', length(dates),'\n'))
}

humid <- as.data.frame(humid)
humid = humid[order(humid$id,humid$date),]
humid <- as.data.table(humid)

humid = humid[, humid_lag1 := shift(humid,type='lag'),.(id)]
humid = humid[, humid_lag2 := shift(humid_lag1,type='lag'),.(id)]
humid = humid[, humid_lag3 := shift(humid_lag2,type='lag'),.(id)]
humid = humid[, humid_lag4 := shift(humid_lag3,type='lag'),.(id)]
humid = humid[, humid_lag5 := shift(humid_lag4,type='lag'),.(id)]
humid = humid[, humid_lag6 := shift(humid_lag5,type='lag'),.(id)]
humid = humid[, humid_lag7 := shift(humid_lag6,type='lag'),.(id)]
humid = humid[, humid_lag8 := shift(humid_lag7,type='lag'),.(id)]
humid = humid[, humid_lag9 := shift(humid_lag8,type='lag'),.(id)]
humid = humid[, humid_lag10 := shift(humid_lag9,type='lag'),.(id)]
humid = humid[, humid_lag11 := shift(humid_lag10,type='lag'),.(id)]
humid = humid[, humid_lag12 := shift(humid_lag11,type='lag'),.(id)]
humid = humid[, humid_lag13 := shift(humid_lag12,type='lag'),.(id)]
humid = humid[, humid_lag14 := shift(humid_lag13,type='lag'),.(id)]
humid = humid[, humid_lag15 := shift(humid_lag14,type='lag'),.(id)]
humid = humid[, humid_lag16 := shift(humid_lag15,type='lag'),.(id)]
humid = humid[, humid_lag17 := shift(humid_lag16,type='lag'),.(id)]
humid = humid[, humid_lag18 := shift(humid_lag17,type='lag'),.(id)]
humid = humid[, humid_lag19 := shift(humid_lag18,type='lag'),.(id)]
humid = humid[, humid_lag20 := shift(humid_lag19,type='lag'),.(id)]

humid$humid_lag713 <- (humid$humid_lag7+humid$humid_lag8+humid$humid_lag9+humid$humid_lag10+humid$humid_lag11+
                         humid$humid_lag12+humid$humid_lag13)/7
humid$humid_lag1420 <- (humid$humid_lag14+humid$humid_lag15+humid$humid_lag16+humid$humid_lag17+humid$humid_lag18+
                          humid$humid_lag19+humid$humid_lag20)/7

humid = humid[, humid_lead1 := shift(humid,type='lead'),.(id)]
humid = humid[, humid_lead2 := shift(humid_lead1,type='lead'),.(id)]
humid = humid[, humid_lead3 := shift(humid_lead2,type='lead'),.(id)]
humid = humid[, humid_lead4 := shift(humid_lead3,type='lead'),.(id)]
humid = humid[, humid_lead5 := shift(humid_lead4,type='lead'),.(id)]
humid = humid[, humid_lead6 := shift(humid_lead5,type='lead'),.(id)]
humid = humid[, humid_lead7 := shift(humid_lead6,type='lead'),.(id)]
humid = humid[, humid_lead8 := shift(humid_lead7,type='lead'),.(id)]
humid = humid[, humid_lead9 := shift(humid_lead8,type='lead'),.(id)]
humid = humid[, humid_lead10 := shift(humid_lead9,type='lead'),.(id)]
humid = humid[, humid_lead11 := shift(humid_lead10,type='lead'),.(id)]
humid = humid[, humid_lead12 := shift(humid_lead11,type='lead'),.(id)]
humid = humid[, humid_lead13 := shift(humid_lead12,type='lead'),.(id)]
humid = humid[, humid_lead14 := shift(humid_lead13,type='lead'),.(id)]
humid = humid[, humid_lead15 := shift(humid_lead14,type='lead'),.(id)]
humid = humid[, humid_lead16 := shift(humid_lead15,type='lead'),.(id)]
humid = humid[, humid_lead17 := shift(humid_lead16,type='lead'),.(id)]
humid = humid[, humid_lead18 := shift(humid_lead17,type='lead'),.(id)]
humid = humid[, humid_lead19 := shift(humid_lead18,type='lead'),.(id)]
humid = humid[, humid_lead20 := shift(humid_lead19,type='lead'),.(id)]

humid$humid_lead713 <- (humid$humid_lead7+humid$humid_lead8+humid$humid_lead9+humid$humid_lead10+humid$humid_lead11+
                          humid$humid_lead12+humid$humid_lead13)/7
humid$humid_lead1420 <- (humid$humid_lead14+humid$humid_lead15+humid$humid_lead16+humid$humid_lead17+humid$humid_lead18+
                           humid$humid_lead19+humid$humid_lead20)/7

humid <- humid[,c("lon","lat","date","humid","humid_lag1","humid_lag2","humid_lag3","humid_lag4","humid_lag5","humid_lag6","humid_lag713",
                  "humid_lag1420","humid_lead1","humid_lead2","humid_lead3","humid_lead4","humid_lead5","humid_lead6","humid_lead713",
                  "humid_lead1420")]
saveRDS(humid,file=paste0(dir_save_GA,'humid_GA.rds'))
