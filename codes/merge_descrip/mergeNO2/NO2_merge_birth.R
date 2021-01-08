###############################################################################
# Project: Black Carbon and MA birth                                          #
# Analysis: Link 1km grid NO2                                                 # 
# Machine: RStudio Server                                                     #
# Input: mabirths_02NOV18.csv (from Anna)                                     #
# Input: NO2 daily prediction on QNAP3                                        #
# Output: #
# Author: Shuxin Dong                                                         #
# Date: 2021-01-06                                                            #
###############################################################################

############################# 0. Setup ########################################
rm(list=ls())
gc()

require(dplyr)
require(magrittr)
require(lubridate)
require(data.table)
library(nabor)
library(fst)

dir_birth <- "/media/qnap3/Shuxin/"
dir_no2 <- "/media/qnap3/Exposure modeling/3 National scale/USA/5 NO2 v1.2000_16/7 Final predictions 1 km/Daily/"
dir_output <- "/media/qnap3/Shuxin/"

dat <- read.csv(paste0(dir_birth, "mabirths_02NOV18.csv"))
# dat$id <- 1:nrow(dat)


############################# 1. NO2 ########################################
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