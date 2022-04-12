###############################################################################
# Project: Causal black carbon on birth weight in MA                          #
# Code: merge census data into birth registry by tract                        #
# Input: birth_geoid`year`.csv and census (.csv) files                        #
# Output: "birth_SES.csv" SES date which could be merged into MA birth data   #
# Author: Shuxin Dong                                                         #
# Date: Oct 21, 2020                                                          #
###############################################################################

############################# 0. Setup ########################################
rm(list = ls())
gc()

library(data.table)

dir_input_geocode <- "/media/qnap3/Shuxin/bc_no2_MAbirth_causal/data/SES/"
dir_input_census <- "/media/qnap3/Covariates/Census/"
dir_output_birth <- "/media/qnap3/Shuxin/bc_no2_MAbirth_causal/data/SES/"

######################### 1. Load and check birth data #########################
## load data
birth_13 <- fread(paste0(dir_input_geocode, "birth_geocode13.csv"),
                  colClasses = c(rep(NA, 5), rep("character",1)),
                  na.strings = "NA")
birth_14 <- fread(paste0(dir_input_geocode, "birth_geocode14.csv"),
                  colClasses = c(rep(NA, 5), rep("character",1)),
                  na.strings = "NA")
birth_15 <- fread(paste0(dir_input_geocode, "birth_geocode15.csv"),
                  colClasses = c(rep(NA, 5), rep("character",1)),
                  na.strings = "NA")
birth_0010 <- fread(paste0(dir_input_geocode, "birth_geocode0010.csv"),
                    drop = c("long", "lat"),
                    colClasses = c(rep(NA, 5), rep("character",2)),
                    na.strings = "NA")

# birth_0010 %>% filter(is.na(geoid2000)|is.na(geoid2010))
# birth_13 %>% filter(is.na(geoid2013))
# birth_14 %>% filter(is.na(geoid2014))
# birth_15 %>% filter(is.na(geoid2015))

# ## found the birth 2000 has missing value, go back.
# birth_0010missing <- birth_0010 %>% filter(is.na(geoid2000)|is.na(geoid2010))
# write.csv(birth_0010missing, file = paste0(dir_input_geocode, "birth_0010missing.csv"))

################## 2. Merge SES variables for 2013-2015########################
################################# 2013 ########################################
## median household income
mhincome13 <- fread(paste0(dir_input_census, 
                           "MedianHouseholdIncome_tract_US_2000_2018/clean/",
                           "medianhouseholdincome_2013.csv"),
                    select = c("GEO_ID", "mhincome"),
                    colClasses = c(NA,NA,NA,NA,"character"))
setDT(mhincome13)
mhincome13 <- na.omit(mhincome13)
mhincome13[, geoid2013 := substr(GEO_ID,10,20)][]
## median house value
mhvalue13 <- fread(paste0(dir_input_census, 
                         "MedianValueofHouse_tract_US_2000_2018/clean/",
                         "medianvalueofhousevalue_2013.csv"),
                  colClasses = c(NA,NA,NA,NA,"character"),
                  select = c("GEO_ID", "mhvalue"))
setDT(mhvalue13)
mhvalue13[mhvalue=="1,000,000+"]$mhvalue <- NA
mhvalue13 <- na.omit(mhvalue13)
mhvalue13[, mhvalue := as.numeric(mhvalue)]
mhvalue13[, geoid2013 := substr(GEO_ID,10,20)]
## percent below poverty
poverty13 <- fread(paste0(dir_input_census, 
                          "PercentbelowPoverty_tract_US_2000_2018/clean/",
                          "percentPovertypoverty_2013.csv"),
                   colClasses = c(NA,NA,NA,NA,"character"),
                   select = c("GEO_ID", "percentPoverty"))
setDT(poverty13)
poverty13[percentPoverty=="-"]$percentPoverty <- NA
poverty13$percentPoverty <- as.numeric(poverty13$percentPoverty)
poverty13 <- na.omit(poverty13)
poverty13[, geoid2013 := substr(GEO_ID,10,20)]
## merge
birth_13 <- merge(birth_13, mhincome13, by = "geoid2013")
birth_13 <- merge(birth_13, mhvalue13, by="geoid2013")
birth_13 <- merge(birth_13, poverty13, by="geoid2013")
birth13_SES <- birth_13[,.(uniqueid_yr, year, mhincome, mhvalue, percentPoverty)]
# fwrite(birth13_SES, paste0(dir_output_birth, "birth13_SES.csv"))

## for 2014
mhincome14 <- fread(paste0(dir_input_census, 
                           "MedianHouseholdIncome_tract_MA_2000_2015/clean/",
                           "medianhouseholdincome_2014.csv"),
                    select = c("GEO_ID", "mhincome"),
                    colClasses = c(NA,NA,NA,NA,"character"))
setDT(mhincome14)
mhincome14 <- na.omit(mhincome14)
mhincome14[, geoid2014 := substr(GEO_ID,10,20)][]

mhvalue14 <- fread(paste0(dir_input_census, 
                          "MedianValueofHouse_tract_2000_2015/clean/",
                          "medianvalueofhousevalue_2014.csv"),
                   colClasses = c(NA,NA,NA,NA,"character"),
                   select = c("GEO_ID", "mhvalue"),)
setDT(mhvalue14)
mhvalue14[mhvalue=="1,000,000+"]$mhvalue <- NA
mhvalue14[, mhvalue := as.numeric(mhvalue)]
mhvalue14 <- na.omit(mhvalue14)
mhvalue14[, geoid2014 := substr(GEO_ID,10,20)][]

poverty14 <- fread(paste0(dir_input_census, 
                          "PercentbelowPoverty_tract_2000_2015/clean/",
                          "percentPovertypoverty_2014.csv"),
                   colClasses = c(NA,NA,NA,NA,"character"),
                   select = c("GEO_ID", "percentPoverty"))
setDT(poverty14)
poverty14[percentPoverty=="-"]$percentPoverty <- NA
poverty14[, percentPoverty := as.numeric(percentPoverty)][]
poverty14 <- na.omit(poverty14)
poverty14[, geoid2014 := substr(GEO_ID,10,20)][]

birth_14 <- merge(birth_14, mhincome14, by="geoid2014")
birth_14 <- merge(birth_14, mhvalue14, by="geoid2014")
birth_14 <- merge(birth_14, poverty14, by="geoid2014")
birth14_SES <- birth_14[,.(uniqueid_yr, year, mhincome, mhvalue, percentPoverty)]
# fwrite(birth14_SES, file = paste0(dir_output_birth, "birth14_SES.csv"))

## for 2015
mhincome15 <- fread(paste0(dir_input_census, 
                           "MedianHouseholdIncome_tract_MA_2000_2015/clean/",
                           "medianhouseholdincome_2015.csv"),
                    colClasses = c(NA,NA,NA,"numeric","character"),
                    select = c("mhincome", "GEO_ID"))
setDT(mhincome15)
mhincome15[!is.na(mhincome)][is.na(as.numeric(mhincome15[!is.na(mhincome)]$mhincome))]
mhincome15[mhincome=="null"]$mhincome <- NA
mhincome15[, mhincome := as.numeric(mhincome)]
mhincome15 <- na.omit(mhincome15)
mhincome15[, geoid2015 := substr(GEO_ID,10,20)]
mhincome15 <- mhincome15[nchar(mhincome15$geoid2015)==11]
mhincome15
setDT(mhincome15)

mhvalue15 <- fread(paste0(dir_input_census, 
                          "MedianValueofHouse_tract_2000_2015/clean/",
                          "medianvalueofhousevalue_2015.csv"),
                   colClasses = c(NA,NA,NA,NA,"character"),
                   select = c("mhvalue", "GEO_ID"))
mhvalue15$mhvalue
mhvalue15[mhvalue=="null"]$mhvalue <- NA
mhvalue15$mhvalue <- as.numeric(mhvalue15$mhvalue)
mhvalue15[, geoid2015 := substr(GEO_ID,10,20)]
mhvalue15 <- mhvalue15[nchar(mhvalue15$geoid2015)==11]
mhvalue15

poverty15 <- fread(paste0(dir_input_census, 
                          "PercentbelowPoverty_tract_2000_2015/clean/",
                          "percentPovertypoverty_2015.csv"),
                   colClasses = c(NA,NA,NA,NA,"character"),
                   select = c("percentPoverty","GEO_ID"))
poverty15[percentPoverty=="-"]
poverty15[percentPoverty=="-"]$percentPoverty <- NA
poverty15$percentPoverty <- as.numeric(poverty15$percentPoverty)
poverty15[, geoid2015 := substr(GEO_ID,10,20)]
poverty15 <- poverty15[nchar(poverty15$geoid2015)==11]
poverty15

birth_15 <- left_join(birth_15, mhincome15, by="geoid2015")
birth_15 <- left_join(birth_15, mhvalue15, by="geoid2015")
birth_15 <- left_join(birth_15, poverty15, by="geoid2015")
birth15_SES <- birth_15[,.(uniqueid_yr, year, mhincome, mhvalue, percentPoverty)]
# fwrite(birth15_SES, file = paste0(dir_output_birth, "birth15_SES.csv"))

################## 3. Merge SES variables for 2011-2012#########################
## get subset
birth_12 <- birth_0010[year==2012]
birth_12[, geoid2000 := NULL]
birth_12
birth_11 <- birth_0010[year==2011]
birth_11[, geoid2000 := NULL]
birth_11
rm(birth_0010)
gc()

## for 2012
mhincome12 <- fread(paste0(dir_input_census, 
                           "MedianHouseholdIncome_tract_MA_2000_2015/clean/",
                           "medianhouseholdincome_2012.csv"),
                    colClasses = c(NA,NA,NA,"numeric","character"),
                    select = c("mhincome", "GEO_ID"))
summary(mhincome12)
mhincome12 <- na.omit(mhincome12)
mhincome12[, geoid2010 := substr(GEO_ID,10,20)]
mhincome12 <- mhincome12[nchar(mhincome12$geoid2010)==11]
setDT(mhincome12)

mhvalue12 <- fread(paste0(dir_input_census, 
                          "MedianValueofHouse_tract_2000_2015/clean/",
                          "medianvalueofhousevalue_2012.csv"),
                   colClasses = c(NA,NA,NA,NA,"character"),
                   select = c("mhvalue", "GEO_ID"))
mhvalue12$mhvalue
mhvalue12[mhvalue=="1,000,000+"]$mhvalue <- NA
mhvalue12$mhvalue <- as.numeric(mhvalue12$mhvalue)
mhvalue12[, geoid2010 := substr(GEO_ID,10,20)]
mhvalue12 <- mhvalue12[nchar(mhvalue12$geoid2010)==11]
mhvalue12 <- na.omit(mhvalue12)
setDT(mhvalue12)

poverty12 <- fread(paste0(dir_input_census, 
                          "PercentbelowPoverty_tract_2000_2015/clean/",
                          "percentPovertypoverty_2012.csv"),
                   colClasses = c(NA,NA,NA,NA,"character"),
                   select = c("percentPoverty","GEO_ID"))
poverty12[percentPoverty=="-"]
poverty12[percentPoverty=="-"]$percentPoverty <- NA
poverty12$percentPoverty <- as.numeric(poverty12$percentPoverty)
poverty12 <- na.omit(poverty12)
poverty12[, geoid2010 := substr(GEO_ID,10,20)]
poverty12 <- poverty12[nchar(poverty12$geoid2010)==11]
setDT(poverty12)

birth_12 <- merge(birth_12, mhincome12, by="geoid2010")
birth_12 <- merge(birth_12, mhvalue12, by="geoid2010")
birth_12 <- merge(birth_12, poverty12, by="geoid2010")
birth12_SES <- birth_12[,.(uniqueid_yr, year, mhincome, mhvalue, percentPoverty)]
# fwrite(birth12_SES, file = paste0(dir_output_birth, "birth12_SES.csv"))

## for 2011
mhincome11 <- fread(paste0(dir_input_census, 
                           "MedianHouseholdIncome_tract_MA_2000_2015/clean/",
                           "medianhouseholdincome_2011.csv"),
                    colClasses = c(NA,NA,NA,"numeric","character"),
                    select = c("mhincome", "GEO_ID"))
summary(mhincome11)
mhincome11[, geoid2010 := substr(GEO_ID,10,20)]
mhincome11 <- mhincome11[nchar(mhincome11$geoid2010)==11]
mhincome11 <- na.omit(mhincome11)
setDT(mhincome11)

mhvalue11 <- fread(paste0(dir_input_census, 
                          "MedianValueofHouse_tract_2000_2015/clean/",
                          "medianvalueofhousevalue_2011.csv"),
                   colClasses = c(NA,NA,NA,NA,"character"),
                   select = c("mhvalue", "GEO_ID"))
summary(mhvalue11)
mhvalue11[mhvalue=="1,000,000+"]$mhvalue <- NA
mhvalue11$mhvalue <- as.numeric(mhvalue11$mhvalue)
mhvalue11[, geoid2010 := substr(GEO_ID,10,20)]
mhvalue11 <- mhvalue11[nchar(mhvalue11$geoid2010)==11]
mhvalue11 <- na.omit(mhvalue11)
setDT(mhvalue11)

poverty11 <- fread(paste0(dir_input_census, 
                          "PercentbelowPoverty_tract_2000_2015/clean/",
                          "percentPovertypoverty_2011.csv"),
                   colClasses = c(NA,NA,NA,NA,NA,NA, "character"),
                   select = c("percentPoverty","GISJOIN"))
poverty11[percentPoverty=="-"]
poverty11[, geoid2010 := paste0(substr(GISJOIN,2,3),substr(GISJOIN,5,7),substr(GISJOIN,9,14))]
poverty11 <- poverty11[nchar(poverty11$geoid2010)==11]
poverty11 <- na.omit(poverty11)
setDT(poverty11)

birth_11 <- merge(birth_11, mhincome11, by="geoid2010")
birth_11 <- merge(birth_11, mhvalue11, by="geoid2010")
birth_11 <- merge(birth_11, poverty11, by="geoid2010")
birth11_SES <- birth_11[,.(uniqueid_yr, year, mhincome, mhvalue, percentPoverty)]
summary(birth11_SES)
fwrite(birth11_SES, file = paste0(dir_output_birth, "birth11_SES.csv"))

################## 4. Merge SES variables for 2000-2010#########################
birth_0010 <- birth_0010[year!=2012&year!=2011]
birth_0010
summary(birth_0010$year)

## median household income
mhincome00 <- fread(paste0(dir_input_census, 
                           "MedianHouseholdIncome_tract_MA_2000_2015/clean/",
                           "medianhouseholdincome_2000.csv"),
                    colClasses = c(NA,NA,NA,"numeric","character"),
                    select = c("mhincome", "GISJOIN"),
                    col.names = c("mhincome2000", "GISJOIN"),
                    )
mhincome00
mhincome00[, geoid2000 := paste0(substr(GISJOIN,2,3),substr(GISJOIN,5,7),substr(GISJOIN,9,14))]
mhincome00
mhincome00[, GISJOIN := NULL]
mhincome00 <- na.omit(mhincome00)

mhincome10 <- fread(paste0(dir_input_census, 
                           "MedianHouseholdIncome_tract_MA_2000_2015/clean/",
                           "medianhouseholdincome_2010.csv"),
                    colClasses = c(NA,NA,NA,"numeric","character"),
                    select = c("mhincome", "GEO_ID"),
                    col.names = c("mhincome2010", "GEO_ID"),
                    )
summary(mhincome10)
mhincome10[, geoid2010 := substr(GEO_ID,10,20)]
mhincome10
mhincome10 <- mhincome10[nchar(mhincome10$geoid2010)==11]
mhincome10 <- na.omit(mhincome10)
mhincome10[, GEO_ID :=NULL]
mhincome10

## median value of houses
mhvalue00 <- fread(paste0(dir_input_census, 
                          "MedianValueofHouse_tract_2000_2015/clean/",
                          "medianvalueofhousevalue_2000.csv"),
                   colClasses = c(NA,NA,NA,NA,"character"),
                   select = c("mhvalue", "GISJOIN"),
                   col.names = c("mhvalue2000", "GISJOIN"),
                   )
summary(mhvalue00)
mhvalue00[, geoid2000 := paste0(substr(GISJOIN,2,3),substr(GISJOIN,5,7),substr(GISJOIN,9,14))]
mhvalue00
mhvalue00[, GISJOIN := NULL]
mhvalue00 <- na.omit(mhvalue00)
mhvalue10 <- fread(paste0(dir_input_census, 
                          "MedianValueofHouse_tract_2000_2015/clean/",
                          "medianvalueofhousevalue_2010.csv"),
                   colClasses = c(NA,NA,NA,NA,"character"),
                   select = c("mhvalue", "GEO_ID"),
                   col.names = c("mhvalue2010", "GEO_ID"),
                   )
summary(mhvalue10)
mhvalue10$mhvalue2010
mhvalue10[mhvalue2010=="1,000,000+"]$mhvalue2010 <- NA
mhvalue10$mhvalue2010 <- as.numeric(mhvalue10$mhvalue2010)
mhvalue10[, geoid2010 := substr(GEO_ID,10,20)]
mhvalue10
mhvalue10 <- mhvalue10[nchar(mhvalue10$geoid2010)==11]
mhvalue10 <- na.omit(mhvalue10)
mhvalue10[, GEO_ID :=NULL]
mhvalue10

## poverty
poverty00 <- fread(paste0(dir_input_census, 
                          "PercentbelowPoverty_tract_2000_2015/clean/",
                          "percentPovertypoverty_2000.csv"),
                   colClasses = c(NA,"character",NA,NA,NA, "character"),
                   select = c("percentPoverty","areakey"),
                   col.names = c("percentPoverty00", "geoid2000"),
                   )
summary(poverty00)
poverty00 <- na.omit(poverty00)
poverty00
poverty10 <- fread(paste0(dir_input_census, 
                          "PercentbelowPoverty_tract_2000_2015/clean/",
                          "percentPovertypoverty_2010.csv"),
                   colClasses = c(NA,NA,NA,NA,NA,NA, "character"),
                   select = c("percentPoverty","GISJOIN"),
                   col.names = c("percentPoverty10", "GISJOIN"),
                   )
summary(poverty10)
poverty10[, geoid2010 := paste0(substr(GISJOIN,2,3),substr(GISJOIN,5,7),substr(GISJOIN,9,14))]
poverty10
poverty10[, GISJOIN := NULL]
poverty10 <- na.omit(poverty10)
poverty10
birth_0010 <- merge(birth_0010, mhincome00, by = "geoid2000")
birth_0010 <- merge(birth_0010, mhvalue00, by = "geoid2000")
birth_0010 <- merge(birth_0010, poverty00, by = "geoid2000")
birth_0010 <- merge(birth_0010, mhincome10, by = "geoid2010")
birth_0010 <- merge(birth_0010, mhvalue10, by = "geoid2010")
birth_0010 <- merge(birth_0010, poverty10, by = "geoid2010")
birth_0010

birth_00 <- birth_0010[year==2000]
birth_00 <- birth_00[,.(uniqueid_yr,year,mhincome2000,mhvalue2000,percentPoverty00)]
names(birth_00) <- c("uniqueid_yr","year","mhincome","mhvalue","percentPoverty")
fwrite(birth_00, file = paste0(dir_output_birth, "birth00_SES.csv"))

birth_10 <- birth_0010[year==2010]
birth_10 <- birth_10[,.(uniqueid_yr,year,mhincome2010,mhvalue2010,percentPoverty10)]
names(birth_10) <- c("uniqueid_yr","year","mhincome","mhvalue","percentPoverty")
fwrite(birth_10, file = paste0(dir_output_birth, "birth10_SES.csv"))

birth_0109 <- birth_0010[year!=2000&year!=2010]
birth_0109
birth_0109[,c("mhincome","mhvalue","percentPoverty") := 
             list(mhincome = (mhincome2010-mhincome2000)/10*(year-2000)+mhincome2000,
                  mhvalue = (mhvalue2010-mhvalue2000)/10*(year-2000)+mhvalue2000,
                  percentPoverty = (percentPoverty10-percentPoverty00)/10*(year-2000)+percentPoverty00)]
birth0109_SES <- birth_0109[,.(uniqueid_yr,year,mhincome,mhvalue,percentPoverty)]
fwrite(birth0109_SES, file = paste0(dir_output_birth, "birth0109_SES.csv"))

######################### 5. Merge with birth data #############################
# birth <- fread(file = "/Users/shuxind/Desktop/BC_birthweight_data/MAbirth_02Nov18/mabirths_02NOV18.csv")
# birth00 <- fread(file = paste0(dir_output_birth, "birth00_SES.csv"),
#                  drop = "V1")
# birth0109 <- fread(file = paste0(dir_output_birth, "birth0109_SES.csv"),
#                  drop = "V1")
# birth10 <- fread(file = paste0(dir_output_birth, "birth10_SES.csv"),
#                  drop = "V1",
#                  colClasses = c(NA, "character", rep(NA,4)))
# birth11 <- fread(file = paste0(dir_output_birth, "birth11_SES.csv"),
#                  colClasses = c(NA, "character", rep(NA,4)),
#                  drop = "V1")
# birth12 <- fread(file = paste0(dir_output_birth, "birth12_SES.csv"),
#                  colClasses = c(NA, "character", rep(NA,4)),
#                  drop = "V1")
# birth13 <- fread(file = paste0(dir_output_birth, "birth13_SES.csv"),
#                  colClasses = c(NA, "character", rep(NA,4)),
#                  drop = "V1")
# birth14 <- fread(file = paste0(dir_output_birth, "birth14_SES.csv"),
#                  colClasses = c(NA, "character", rep(NA,4)),
#                  drop = "V1")
# birth15 <- fread(file = paste0(dir_output_birth, "birth15_SES.csv"),
#                  colClasses = c(NA, "character", rep(NA,4)),
#                  drop = "V1")
birth_SES <- rbind(birth00, birth0109, birth10, birth11, birth12, birth13, birth14, birth15)
setDT(birth_SES)
fwrite(birth_SES, file = paste0(dir_output_birth, "birth_SES.csv"))
