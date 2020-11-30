###############################################################################
# Project: Causal black carbon on birthweight in MA - linear                  #
# Code: Estimate GPS using gradient boosting                                  #
# Machine: RCE                                                                #
# Author: Yaguang Wei / Shuxin Dong                                           #
###############################################################################

############################# 0. Setup ##############################
rm(list=ls())
gc()

library(lubridate)
library(dplyr)
library(h2o)
library(data.table)

h2o.init()

dir_data <- "/Users/dongshuxin/Desktop/MA_births/Reformat_Descrp_Tfm/"
dir_results <- "/Users/dongshuxin/Desktop/MA_births/"

## load covariates
birth <- readRDS(paste0(dir_data, "birth_tfm.rds"))
# names(birth)
# [1] "uniqueid_yr"    "year"           "sex"            "married"        "mage"          
# [6] "mrace"          "m_edu"          "cigdpp"         "cigddp"         "clinega"       
# [11] "kotck"          "pncgov"         "bwg"            "rf_db_gest"     "rf_db_other"   
# [16] "rf_hbp_chronic" "rf_hbp_pregn"   "rf_cervix"      "rf_prev_4kg"    "rf_prev_sga"   
# [21] "med_hs_inc"     "bc_30d"         "bc_90d"         "bc_280d"        "parit_cat"     
# [26] "m_wg_cat"       "log_med_hs_inc"

############################# 1. data manipulation ##############################
birth <- birth[,c("uniqueid_yr","year","sex","married","mage","mrace","m_edu","cigdpp","cigddp",
                  "clinega","kotck","pncgov", "bwg","rf_db_gest","rf_db_other","rf_hbp_chronic",
                  "rf_hbp_pregn","rf_cervix","rf_prev_4kg","rf_prev_sga","bc_30d","bc_90d",
                  "bc_280d","parit_cat","m_wg_cat","log_med_hs_inc")]

#mysample <- birth[sample(1:nrow(birth), 200, replace=FALSE),]

############################# 2. predicting using gradiant boosting ##############################
#birth <- as.h2o(mysample)
birth <- as.h2o(birth)

## bc_30d
xvar_bc_30d <- c("year","sex","married","mage","mrace","m_edu","cigdpp","cigddp",
                 "clinega","kotck","pncgov","rf_db_gest","rf_db_other","rf_hbp_chronic",
                 "rf_hbp_pregn","rf_cervix","rf_prev_4kg","rf_prev_sga","parit_cat",
                 "m_wg_cat","log_med_hs_inc")
yvar_bc_30d <- "bc_30d"
mod_gbm_bc_30d <-  h2o.gbm(x=xvar_bc_30d,y=yvar_bc_30d, 
                           ntrees=200, seed=1, max_depth = 10,
                           training_frame = birth)
pred_gbm_bc_30d <- h2o.predict(mod_gbm_bc_30d, birth)
h2o:::.h2o.garbageCollect()
# birth$resid_gbm_bc_30d <- birth$bc_30d - pred_gbm_bc_30d
# sigma_gbm_bc_30d <- h2o.sd(birth$resid_gbm_bc_30d)

# ## bc_90d
# xvar_bc_90d <- c("year","sex","married","mage","mrace","m_edu","cigdpp","cigddp",
#                  "clinega","kotck","pncgov","rf_db_gest","rf_db_other","rf_hbp_chronic",
#                  "rf_hbp_pregn","rf_cervix","rf_prev_4kg","rf_prev_sga","parit_cat",
#                  "m_wg_cat","log_med_hs_inc")
# yvar_bc_90d <- "bc_90d"
# mod_gbm_bc_90d <-  h2o.gbm(x=xvar_bc_90d,y=yvar_bc_90d, 
#                            ntrees=200, seed=1, max_depth = 10,
#                            training_frame = birth)
# pred_gbm_bc_90d <- h2o.predict(mod_gbm_bc_90d, birth)
# h2o:::.h2o.garbageCollect()
# birth$resid_gbm_bc_90d <- birth$bc_90d - pred_gbm_bc_90d
# sigma_gbm_bc_90d <- h2o.sd(birth$resid_gbm_bc_90d)
# 
# ## bc_280d
# xvar_bc_280d <- c("year","sex","married","mage","mrace","m_edu","cigdpp","cigddp",
#                   "clinega","kotck","pncgov","rf_db_gest","rf_db_other","rf_hbp_chronic",
#                   "rf_hbp_pregn","rf_cervix","rf_prev_4kg","rf_prev_sga","parit_cat",
#                   "m_wg_cat","log_med_hs_inc")
# yvar_bc_280d <- "bc_280d"
# mod_gbm_bc_280d <-  h2o.gbm(x=xvar_bc_280d,y=yvar_bc_280d, 
#                             ntrees=200, seed=1, max_depth = 10,
#                             training_frame = birth)
# pred_gbm_bc_280d <- h2o.predict(mod_gbm_bc_280d, birth)
# h2o:::.h2o.garbageCollect()
# birth$resid_gbm_bc_280d <- birth$bc_280d - pred_gbm_bc_280d
# sigma_gbm_bc_280d <- h2o.sd(birth$resid_gbm_bc_280d)

############################# 3. estimating GPS and save ##############################
birth <- as.data.frame(birth)
pred_gbm_bc_30d <- as.data.frame(pred_gbm_bc_30d)
birth$pred_bc_30d <- pred_gbm_bc_30d$predict

saveRDS(birth, paste0(dir_results,"GBMplotCheck.rds"))
pdf(paste0(dir_results,"GBMpred_obs.pdf"))
plot(birth$bc_30d, birth$pred_bc_30d, main="GBM(depth5_tree200)")
lines(lowess(birth$bc_30d, birth$pred_bc_30d), col="red", lwd=2)
dev.off()

# birth$gps_bc_30d <- dnorm(birth$resid_gbm_bc_30d, mean=0, sd=sigma_gbm_bc_30d)
# birth$gps_bc_90d <- dnorm(birth$resid_gbm_bc_90d, mean=0, sd=sigma_gbm_bc_90d)
# birth$gps_bc_280d <- dnorm(birth$resid_gbm_bc_280d, mean=0, sd=sigma_gbm_bc_280d)
# 
# gps_gbm_birth <- birth[,c("uniqueid_yr", "bwg", "bc_30d", "bc_90d", "bc_280d",
#                           "gps_bc_30d", "gps_bc_90d","gps_bc_280d")]
# 
# saveRDS(gps_gbm_birth, file = paste0(dir_results,'gps_gbm_birth_depth10_tree200.rds')) # change
# beepr::beep(3)
