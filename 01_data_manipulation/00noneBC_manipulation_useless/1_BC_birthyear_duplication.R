###############################################################################
# Project: Causal black carbon on birth weight in MA                          #
# Code: Generate BC data by birth year and detect duplication prob            #
# Author: Shuxin Dong                                                         #
# Date: Sep 11, 2020                                                          #
###############################################################################

############################# 0. Setup ########################################
rm(list = ls())
gc()

library(dplyr)

dir_rawdata <- "/Users/shuxind/Desktop/BC_birthweight_data/" 
dir_newdata <- "/Users/shuxind/Desktop/" 

############################# 1. get the ID list for each birth year ##########
ref <- read.csv(paste0(dir_rawdata, "mabirths_02NOV18.csv"))
b_year_id <- ref[, c("uniqueid_yr", "year")]
rm(ref)
gc()
b_year_id <- b_year_id[b_year_id$year!=2000, ]

############################# 2. generate BC datasets for each birth year #####
bad_id_all <- NULL

# # take b_year 2001 as an example ##############################################
# bad_id <- NULL
# fixed_id <- NULL
# # get the complete exposure history for those who are born in certain year
# raw <- rbind(readRDS(paste0(dir_rawdata, "y2001.rds")),
#              readRDS(paste0(dir_rawdata, "y2000.rds")))
# raw$uniqueid_yr <- as.character(raw$uniqueid_yr)
# 
# # 2.1 select people based on ID by birth year
# raw <- raw %>% 
#   filter(uniqueid_yr %in% b_year_id$uniqueid_yr[b_year_id$year==2001]) %>% 
#   select(-OBJECTI)
# 
# # 2.2 detect problematic IDs (either missing or duplicated)
# prob_id <- raw %>% count(uniqueid_yr) %>% filter(n!=365) %>% .$uniqueid_year
# # *2.3 fix results
# if (!is.null(prob_id)){
#   fix_result <- raw %>% filter(uniqueid_yr %in% prob_id) %>% distinct() %>% 
#     count(uniqueid_yr) 
#   fixed_id <- fix_result %>% filter(n==365) %>% .$uniqueid_yr
#   bad_id <- fix_result %>% filter(n!=365) %>% .$uniqueid_yr
#   if (!is.null(bad_id)){
#     bad_id_all <- rbind(bad_id_all, data.frame(year = 2001, bad_id))
#   }
# }
# 
# # 2.4 generate clean BC datasets by birth year
# if (!is.null(bad_id)){
#   raw <- raw %>% filter(!uniqueid_yr %in% bad_id) # delete obs with bad id
# }
# if (!is.null(fixed_id)){
#   new <- rbind(raw %>% filter(uniqueid_yr %in% fixed_id) %>% distinct(),
#                raw %>% filter(!uniqueid_yr %in% fixed_id))
# } else{
#   new <- raw
# }
# # 2.5 save new dataset
# save(new, "b_y2001.rds")

for (yr in 2001:2015){
  bad_id <- NULL
  fixed_id <- NULL
  # get the complete exposure history 
  raw <- rbind(readRDS(paste0(dir_rawdata, "y", yr, ".rds")),
               readRDS(paste0(dir_rawdata, "y", yr-1, ".rds")))
  raw$uniqueid_yr <- as.character(raw$uniqueid_yr)
  print(paste0("finish importing complete exposure history of ", yr))

  # 2.1 select people based on ID by birth year
  raw <- raw %>%
    filter(uniqueid_yr %in% b_year_id$uniqueid_yr[b_year_id$year==yr]) %>%
    select(-OBJECTI)
  print(paste0("finish selecting ", yr, "births"))
  # 2.2 detect problematic IDs (either missing or duplicated)
  prob_id <- raw %>% count(uniqueid_yr) %>% filter(n!=365) %>% .$uniqueid_year
  # *2.3 fix results
  if (!is.null(prob_id)){
    fix_result <- raw %>% filter(uniqueid_yr %in% prob_id) %>% distinct() %>%
      count(uniqueid_yr)
    fixed_id <- fix_result %>% filter(n==365) %>% .$uniqueid_yr
    bad_id <- fix_result %>% filter(n!=365) %>% .$uniqueid_yr
    if (!is.null(bad_id)){
      bad_id_all <- rbind(bad_id_all, data.frame(year = 2001, bad_id))
    }
  }
  # 2.4 generate clean BC datasets by birth year
  if (!is.null(bad_id)){
    raw <- raw %>% filter(!uniqueid_yr %in% bad_id) # delete obs with bad id
  }
  if (!is.null(fixed_id)){
    new <- rbind(raw %>% filter(uniqueid_yr %in% fixed_id) %>% distinct(),
                 raw %>% filter(!uniqueid_yr %in% fixed_id))
  } else{
    new <- raw
  }
  # 2.5 save new dataset
  save(new, file = paste0(dir_newdata, "b_y", yr, ".rds"))
}
save(bad_id_all, file = paste0(dir_newdata, "bad_id.rds"))