###############################################################################
# Project: Causal black carbon on birth weight in MA                          #
# Code: GBM to estimate GPS with balance stopping rule                        #
# Code : (GBM example) from three non-overlapping periods                     #
# Input: clean birth data, bc averages                                        #
# Output: best model                                                          #
# Author: Shuxin Dong                                                         #
# Date: Sep 28, 2020                                                          #
###############################################################################

############################# 0. Setup ########################################
rm(list = ls())
gc()

library(WeightIt)
library(dplyr)

dir_input <- "/Users/shuxind/Desktop/BC_birthweight_data/"
# dir_input <- "/nfs/home/S/shd968/shared_space/ci3_shd968_proj/"

## load data
birth <- readRDS(paste0(dir_input, "birth_tfm.rds"))
birth <- sample_n(birth, floor(0.1*dim(birth)[1])) # sample as an example

## set hyperparameteres range
n_trees <- 3000
interaction_depth <- 2:5

############################# 1. data manipulation ############################
birth <- birth[, c("uniqueid_yr","year","sex","married","mage","mrace","m_edu",
                   "cigdpp","cigddp", "clinega","kotck","pncgov", "bwg",
                   "rf_db_gest","rf_db_other","rf_hbp_chronic",
                   "rf_hbp_pregn","rf_cervix","rf_prev_4kg","rf_prev_sga",
                   "bc_30d","bc_90d", "bc_280d",
                   "parit_cat","m_wg_cat","log_med_hs_inc")]
bin_cols.name <- c("pncgov","rf_db_gest","rf_db_other","rf_hbp_chronic",
                   "rf_hbp_pregn","rf_cervix","rf_prev_4kg","rf_prev_sga")
birth[bin_cols.name] <- sapply(birth[bin_cols.name],as.numeric)

birth$bc_3090d <- (birth$bc_90d*90 - birth$bc_30d*30) / (90-30)
birth$bc_90280d <- (birth$bc_280d*280 - birth$bc_90d*90) / (280-90)
birth <- birth %>% select (-bc_90d, -bc_280d)

############################# 2. predicting using gradient boosting ###########
## Tuning hyperparameters
gbm_30d <- weightit(bc_30d ~ year + sex + married + mage + mrace + m_edu + cigdpp + cigddp + 
                      clinega + kotck + pncgov + rf_db_gest + rf_db_other + rf_hbp_chronic + 
                      rf_hbp_pregn + rf_cervix + rf_prev_4kg + rf_prev_sga + parit_cat + 
                      m_wg_cat + log_med_hs_inc + bc_3090d + bc_90280d,
                    data = birth,
                    method = "gbm", 
                    stop.method = "p.max",
                    interaction.depth = interaction_depth,
                    n.trees = n_trees,
                    trim.at = .97,
                    verbose = TRUE)
gbm_30d$info$best.tune

gbm_3090d <- weightit(bc_3090d ~ year + sex + married + mage + mrace + m_edu + cigdpp + cigddp + 
                      clinega + kotck + pncgov + rf_db_gest + rf_db_other + rf_hbp_chronic + 
                      rf_hbp_pregn + rf_cervix + rf_prev_4kg + rf_prev_sga + parit_cat + 
                      m_wg_cat + log_med_hs_inc + bc_30d + bc_90280d,
                    data = birth,
                    method = "gbm", 
                    stop.method = "p.max",
                    interaction.depth = interaction_depth,
                    n.trees = n_trees,
                    trim.at = .97,
                    verbose = TRUE)
gbm_3090d$info$best.tune

gbm_90280d <- weightit(bc_3090d ~ year + sex + married + mage + mrace + m_edu + cigdpp + cigddp + 
                        clinega + kotck + pncgov + rf_db_gest + rf_db_other + rf_hbp_chronic + 
                        rf_hbp_pregn + rf_cervix + rf_prev_4kg + rf_prev_sga + parit_cat + 
                        m_wg_cat + log_med_hs_inc + bc_30d + bc_3090d,
                      data = birth,
                      method = "gbm", 
                      stop.method = "p.max",
                      interaction.depth = interaction_depth,
                      n.trees = n_trees,
                      trim.at = .97,
                      verbose = TRUE)
gbm_90280d$info$best.tune

gc()