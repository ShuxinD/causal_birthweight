###############################################################################
# Project: Causal black carbon on birth weight in MA                          #
# Code: Detect differences of BC effects in different time periods            #
# Code : (GBM example), fit GPS model for three periods and calculate weights #
# Input: clean birth data and bc data                                         #
# Output: bc averages and corresponding weights                               #
# Author: Shuxin Dong                                                         #
# Date: Sep 17, 2020                                                          #
###############################################################################

############################# 0. Setup ########################################
rm(list=ls())
gc()

library(gbm)
library(dplyr)

dir_input <- "/Users/shuxind/Desktop/BC_birthweight_data/Reformat_Descrp_Tfm/" # change
dir_output <- "/Users/shuxind/Desktop/BC_birthweight_analysis/Fit_psm_code/result" # change
new_data_name <- "gbm_bc_wt.rds"

n_cores <- 8
n_trees <- 200
interaction_depth <- 2

## load covariates
birth <- readRDS(paste0(dir_input, "birth_tfm.rds"))
# names(birth)
# [1] "uniqueid_yr"    "year"           "sex"            "married"        "mage"
# [6] "mrace"          "m_edu"          "cigdpp"         "cigddp"         "clinega"
# [11] "kotck"          "pncgov"         "bwg"            "rf_db_gest"     "rf_db_other"
# [16] "rf_hbp_chronic" "rf_hbp_pregn"   "rf_cervix"      "rf_prev_4kg"    "rf_prev_sga"
# [21] "med_hs_inc"     "bc_30d"         "bc_90d"         "bc_280d"        "parit_cat"
# [26] "m_wg_cat"       "log_med_hs_inc"

############################# 1. data manipulation ############################
birth <- birth[, c("uniqueid_yr","year","sex","married","mage","mrace","m_edu","cigdpp","cigddp",
                   "clinega","kotck","pncgov", "bwg","rf_db_gest","rf_db_other","rf_hbp_chronic",
                   "rf_hbp_pregn","rf_cervix","rf_prev_4kg","rf_prev_sga","bc_30d","bc_90d",
                   "bc_280d","parit_cat","m_wg_cat","log_med_hs_inc")]
bin_cols.name <- c("pncgov","rf_db_gest","rf_db_other","rf_hbp_chronic",
                   "rf_hbp_pregn","rf_cervix","rf_prev_4kg","rf_prev_sga")
birth[bin_cols.name] <- sapply(birth[bin_cols.name],as.numeric)
# birth <- birth[sample(1:nrow(birth), 400, replace=FALSE),] ## sample
# birth <- birth %>% filter(bc_30d<0.6 & bc_90d<0.6 & bc_280d<0.6) # sensitivity analysis

birth$bc_3090d <- (birth$bc_90d*90 - birth$bc_30d*30) / (90-30)
birth$bc_90280d <- (birth$bc_280d*280 - birth$bc_90d*90) / (280-90)

############################# 2. predicting using gradient boosting ###########
## bc_30d: 0-30 days prior to the delivery date
mod_gbm_bc_30d <- gbm(bc_30d ~ year + sex + married + mage + mrace + m_edu + cigdpp + cigddp + 
                        clinega + kotck + pncgov + rf_db_gest + rf_db_other + rf_hbp_chronic + 
                        rf_hbp_pregn + rf_cervix + rf_prev_4kg + rf_prev_sga + parit_cat + 
                        m_wg_cat + log_med_hs_inc + bc_3090d + bc_90280d, 
                      distribution = "gaussian", data = birth,
                      n.trees = n_trees, interaction.depth = interaction_depth,
                      n.cores = n_cores)
birth$pred_gbm_bc_30d <- predict.gbm(mod_gbm_bc_30d, birth, n.trees = n_trees, 
                                     single.tree = FALSE)
gc()
birth$resid_gbm_bc_30d <- birth$bc_30d - birth$pred_gbm_bc_30d
birth$sigma_gbm_bc_30d <- sd(birth$resid_gbm_bc_30d)

## bc_3090d: 30-90 days prior to the delivery date
mod_gbm_bc_3090d <- gbm(bc_3090d ~ year + sex + married + mage + mrace + m_edu + cigdpp + cigddp + 
                        clinega + kotck + pncgov + rf_db_gest + rf_db_other + rf_hbp_chronic + 
                        rf_hbp_pregn + rf_cervix + rf_prev_4kg + rf_prev_sga + parit_cat + 
                        m_wg_cat + log_med_hs_inc + bc_30d + bc_90280d, 
                      distribution = "gaussian", data = birth,
                      n.trees = n_trees, interaction.depth = interaction_depth,
                      n.cores = n_cores)
birth$pred_gbm_bc_3090d <- predict.gbm(mod_gbm_bc_3090d, birth, n.trees = 200, 
                                       single.tree = FALSE)
gc()
birth$resid_gbm_bc_3090d <- birth$bc_3090d - birth$pred_gbm_bc_3090d
birth$sigma_gbm_bc_3090d <- sd(birth$resid_gbm_bc_3090d)

## bc_90280d: 90-280 days prior to the delivery date
mod_gbm_bc_90280d <- gbm(bc_90280d ~ year + sex + married + mage + mrace + m_edu + cigdpp + cigddp + 
                         clinega + kotck + pncgov + rf_db_gest + rf_db_other + rf_hbp_chronic + 
                         rf_hbp_pregn + rf_cervix + rf_prev_4kg + rf_prev_sga + parit_cat + 
                         m_wg_cat + log_med_hs_inc + bc_30d + bc_3090d,
                         distribution = "gaussian", data = birth,
                         n.trees = n_trees, interaction.depth = interaction_depth,
                         n.cores = n_cores)
birth$pred_gbm_bc_90280d <- predict.gbm(mod_gbm_bc_90280d, birth, 
                                        n.trees = 200, single.tree = FALSE)
gc()
birth$resid_gbm_bc_90280d <- birth$bc_90280d - birth$pred_gbm_bc_90280d
birth$sigma_gbm_bc_90280d <- sd(birth$resid_gbm_bc_90280d)

############################# 3. estimate GPS, calculate weights and save #####
# calculate generalized propensity scores from three models
birth$gps_bc_30d <- with(birth, dnorm(resid_gbm_bc_30d, mean=0, sd=sigma_gbm_bc_30d))
birth$gps_bc_3090d <- with(birth, dnorm(resid_gbm_bc_3090d, mean=0, sd=sigma_gbm_bc_3090d))
birth$gps_bc_90280d <- with(birth, dnorm(resid_gbm_bc_90280d, mean=0, sd=sigma_gbm_bc_90280d))

# subset the original data to bc related variables - save space
gps <- birth[,c("uniqueid_yr", "bwg", "bc_30d", "bc_3090d", "bc_90280d",
                "gps_bc_30d", "gps_bc_3090d","gps_bc_90280d")]

# calculate weights
gps$numrt_bc_30d <- with(gps, dnorm(bc_30d, mean = mean(bc_30d), sd=sd(bc_30d))) # numerator of stabilized weight
gps$wt_bc_30d_r <- with(gps, numrt_bc_30d/gps_bc_30d) # raw stabilized weight (without truncation)
gps$numrt_bc_3090d <- with(gps, dnorm(bc_3090d, mean = mean(bc_3090d), sd=sd(bc_3090d)))
gps$wt_bc_3090d_r <- with(gps, numrt_bc_3090d/gps_bc_3090d)
gps$numrt_bc_90280d <- with(gps, dnorm(bc_90280d, mean = mean(bc_90280d), sd=sd(bc_90280d)))
gps$wt_bc_90280d_r <- with(gps, numrt_bc_90280d/gps_bc_90280d)

# truncate
gps <- gps %>% mutate(wt_bc_30d = ifelse(wt_bc_30d_r<quantile(wt_bc_30d_r,0.025),quantile(wt_bc_30d_r,0.025), wt_bc_30d_r),
                      wt_bc_30d = ifelse(wt_bc_30d_r>quantile(wt_bc_30d_r,0.975),quantile(wt_bc_30d_r,0.975), wt_bc_30d_r),
                      wt_bc_3090d = ifelse(wt_bc_3090d_r<quantile(wt_bc_3090d_r,0.025),quantile(wt_bc_3090d_r,0.025), wt_bc_3090d_r),
                      wt_bc_3090d = ifelse(wt_bc_3090d_r>quantile(wt_bc_3090d_r,0.975),quantile(wt_bc_3090d_r,0.975), wt_bc_3090d_r),
                      wt_bc_90280d = ifelse(wt_bc_90280d_r<quantile(wt_bc_90280d_r,0.025),quantile(wt_bc_90280d_r,0.025), wt_bc_90280d_r),
                      wt_bc_90280d = ifelse(wt_bc_90280d_r>quantile(wt_bc_90280d_r,0.975),quantile(wt_bc_90280d_r,0.975), wt_bc_90280d_r))

# new dataset and save
ipw <- gps %>% select(uniqueid_yr, bwg, bc_30d, bc_3090d, bc_90280d, 
                      wt_bc_30d, wt_bc_3090d, wt_bc_90280d)
saveRDS(ipw, paste0(dir_output, new_data_name))
summary(ipw)

gc()