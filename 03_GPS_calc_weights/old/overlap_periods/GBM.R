###############################################################################
# Project: Causal black carbon on birthweight in MA - linear                  #
# Code: Estimate GPS using gradient boosting                                  #
# Author: Shuxin Dong                                                         #
# Date: May 08, 2020
###############################################################################

############################# 0. Setup ##############################
rm(list=ls())
gc()

library(gbm)
library(dplyr)

dir_data <- "/Users/dongshuxin/Desktop/MA_births/Reformat_Descrp_Tfm/"
dir_results <- "/Users/dongshuxin/Desktop/MA_births/GPS_result/"

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
bin_cols.name <- c("pncgov","rf_db_gest","rf_db_other","rf_hbp_chronic",
                   "rf_hbp_pregn","rf_cervix","rf_prev_4kg","rf_prev_sga")
birth[bin_cols.name] <- sapply(birth[bin_cols.name],as.numeric)
# birth <- birth[sample(1:nrow(birth), 400, replace=FALSE),] ## sample
birth <- birth %>% filter(bc_30d<0.6 & bc_90d<0.6 & bc_280d<0.6) # sensitivity analysis

############################# 2. predicting using gradiant boosting ##############################
## bc_30d
mod_gbm_bc_30d <- gbm(bc_30d ~ year + sex + married + mage + mrace + m_edu + cigdpp + cigddp + 
                        clinega + kotck + pncgov + rf_db_gest + rf_db_other + rf_hbp_chronic + 
                        rf_hbp_pregn + rf_cervix + rf_prev_4kg + rf_prev_sga + parit_cat + 
                        m_wg_cat + log_med_hs_inc, distribution = "gaussian", data = birth,
                      n.trees = 200, interaction.depth = 2)
birth$pred_gbm_bc_30d <- predict.gbm(mod_gbm_bc_30d, birth, n.trees = 200, single.tree = FALSE)
gc()
birth$resid_gbm_bc_30d <- birth$bc_30d - birth$pred_gbm_bc_30d
birth$sigma_gbm_bc_30d <- sd(birth$resid_gbm_bc_30d)

## bc_90d
mod_gbm_bc_90d <- gbm(bc_90d ~ year + sex + married + mage + mrace + m_edu + cigdpp + cigddp + 
                        clinega + kotck + pncgov + rf_db_gest + rf_db_other + rf_hbp_chronic + 
                        rf_hbp_pregn + rf_cervix + rf_prev_4kg + rf_prev_sga + parit_cat + 
                        m_wg_cat + log_med_hs_inc, distribution = "gaussian", data = birth,
                      n.trees = 200, interaction.depth = 2)
birth$pred_gbm_bc_90d <- predict.gbm(mod_gbm_bc_90d, birth, n.trees = 200, single.tree = FALSE)
gc()
birth$resid_gbm_bc_90d <- birth$bc_90d - birth$pred_gbm_bc_90d
birth$sigma_gbm_bc_90d <- sd(birth$resid_gbm_bc_90d)

## bc_280d
mod_gbm_bc_280d <- gbm(bc_280d ~ year + sex + married + mage + mrace + m_edu + cigdpp + cigddp + 
                        clinega + kotck + pncgov + rf_db_gest + rf_db_other + rf_hbp_chronic + 
                        rf_hbp_pregn + rf_cervix + rf_prev_4kg + rf_prev_sga + parit_cat + 
                        m_wg_cat + log_med_hs_inc, distribution = "gaussian", data = birth,
                      n.trees = 200, interaction.depth = 2)
birth$pred_gbm_bc_280d <- predict.gbm(mod_gbm_bc_280d, birth, n.trees = 200, single.tree = FALSE)
gc()
birth$resid_gbm_bc_280d <- birth$bc_280d - birth$pred_gbm_bc_280d
birth$sigma_gbm_bc_280d <- sd(birth$resid_gbm_bc_280d)

############################# 3. estimating GPS and save ##############################
birth$gps_bc_30d <- with(birth, dnorm(resid_gbm_bc_30d, mean=0, sd=sigma_gbm_bc_30d))
birth$gps_bc_90d <- with(birth, dnorm(resid_gbm_bc_90d, mean=0, sd=sigma_gbm_bc_90d))
birth$gps_bc_280d <- with(birth, dnorm(resid_gbm_bc_280d, mean=0, sd=sigma_gbm_bc_280d))

gps_gbm_birth <- birth[,c("uniqueid_yr", "bwg", "bc_30d", "bc_90d", "bc_280d",
                          "gps_bc_30d", "gps_bc_90d","gps_bc_280d")]

saveRDS(gps_gbm_birth, file = paste0(dir_results,'gps_gbm_sens.rds')) # change


