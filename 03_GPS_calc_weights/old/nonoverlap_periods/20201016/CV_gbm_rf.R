###############################################################################
# Project: Causal black carbon on birth weight in MA                          #
# Code: Detect differences of BC effects in different non overlap time periods#
# Code : fit GPS model with CV GBM and RF                                     #
# Input: clean birth data and bc data                                         #
# Output: bc averages and corresponding weights                               #
# Author: Shuxin Dong                                                         #
# Date: Sep 21, 2020                                                          #
###############################################################################

############################# 0. Setup ########################################
rm(list=ls())
gc()

library(caret)
library(dplyr)
library(doParallel)
cl <- makeCluster(detectCores())
registerDoParallel(cl)

# dir_input <- "/Users/shuxind/Desktop/BC_birthweight_data/" # change
dir_input <- "/nfs/home/S/shd968/shared_space/ci3_shd968_proj/"
# dir_output <- "/Users/shuxind/Desktop/BC_birthweight_analysis/03_GPS_calc_weights/for_time_period/result" # change
dir_output <- "/nfs/home/S/shd968/shared_space/ci3_shd968_proj/"
new_data_name1 <- "gbm_bc_wt.rds"
new_data_name2 <- "rf_bc_wt.rds"

n_cores <- 2

## load covariates
birth <- readRDS(paste0(dir_input, "birth_tfm.rds"))
# birth <- sample_n(birth, 400) # sample as an example

############################# 1. data manipulation ############################
birth <- birth[, c("uniqueid_yr","year","sex","married","mage","mrace","m_edu","cigdpp","cigddp",
                   "clinega","kotck","pncgov", "bwg","rf_db_gest","rf_db_other","rf_hbp_chronic",
                   "rf_hbp_pregn","rf_cervix","rf_prev_4kg","rf_prev_sga","bc_30d","bc_90d",
                   "bc_280d","parit_cat","m_wg_cat","log_med_hs_inc")]
bin_cols.name <- c("pncgov","rf_db_gest","rf_db_other","rf_hbp_chronic",
                   "rf_hbp_pregn","rf_cervix","rf_prev_4kg","rf_prev_sga")
birth[bin_cols.name] <- sapply(birth[bin_cols.name],as.numeric)

birth$bc_3090d <- (birth$bc_90d*90 - birth$bc_30d*30) / (90-30)
birth$bc_90280d <- (birth$bc_280d*280 - birth$bc_90d*90) / (280-90)
birth <- birth %>% select (-bc_90d, -bc_280d)

############################# 2. predicting using gradient boosting ###########
set.seed(998)
inTraining <- createDataPartition(birth$bwg, p = .75, list = FALSE)
training <- birth[ inTraining,]
testing  <- birth[-inTraining,]

fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10)

gbmGrid <-  expand.grid(interaction.depth = c(2, 3, 4, 5), 
                        n.trees = c(150, 200, 250, 300, 350), 
                        shrinkage = 0.05,
                        n.minobsinnode = 5) # change this
nrow(gbmGrid)

## bc_30d: 0-30 days prior to the delivery date
gbmFit1 <- train(bc_30d ~ year + sex + married + mage + mrace + m_edu + cigdpp + cigddp + 
                   clinega + kotck + pncgov + rf_db_gest + rf_db_other + rf_hbp_chronic + 
                   rf_hbp_pregn + rf_cervix + rf_prev_4kg + rf_prev_sga + parit_cat + 
                   m_wg_cat + log_med_hs_inc + bc_3090d + bc_90280d, 
                 data = training, 
                 method = "gbm", 
                 trControl = fitControl, 
                 verbose = FALSE, 
                 tuneGrid = gbmGrid)
gbmFit1

trellis.par.set(caretTheme())
plot(gbmFit1) 

stopCluster(cl)