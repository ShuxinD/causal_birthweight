###############################################################################
# Project: Causal black carbon on birthweight in MA -  linear                 #
# Code: Estimate GPS using m-out-of-n random forest                           #
# Machine: RCE                                                                #
# Author: Yaguang Wei / Shuxin Dong                                           #
###############################################################################

##### Setting up -------------------------------------------------------------------------
rm(list=ls())
gc()

library(readr)
library(parallel)
library(randomForest)
library(lubridate)
library(doParallel)

dir_data <- "/Users/dongshuxin/Desktop/MA_births/Reformat_Descrp_Tfm/" ## Anna change here
dir_results <- "/Users/dongshuxin/Desktop/MA_births/" ## Anna change here

##### Loading the data --------------------------------------------------------
birth <- readRDS(paste0(dir_data, "birth_tfm.rds"))
# names(birth)
# [1] "uniqueid_yr"    "year"           "sex"            "married"        "mage"          
# [6] "mrace"          "m_edu"          "cigdpp"         "cigddp"         "clinega"       
# [11] "kotck"          "pncgov"         "bwg"            "rf_db_gest"     "rf_db_other"   
# [16] "rf_hbp_chronic" "rf_hbp_pregn"   "rf_cervix"      "rf_prev_4kg"    "rf_prev_sga"   
# [21] "med_hs_inc"     "bc_30d"         "bc_90d"         "bc_280d"        "parit_cat"     
# [26] "m_wg_cat"       "log_med_hs_inc"
birth <- birth[,c("uniqueid_yr","year","sex","married","mage","mrace","m_edu","cigdpp","cigddp",
                  "clinega","kotck","pncgov", "bwg","rf_db_gest","rf_db_other","rf_hbp_chronic",
                  "rf_hbp_pregn","rf_cervix","rf_prev_4kg","rf_prev_sga","bc_30d","bc_90d",
                  "bc_280d","parit_cat","m_wg_cat","log_med_hs_inc")]

simulated <- birth[,c("year","sex","married","mage","mrace","m_edu","cigdpp","cigddp",
                      "clinega","kotck","pncgov", "bwg","rf_db_gest","rf_db_other","rf_hbp_chronic",
                      "rf_hbp_pregn","rf_cervix","rf_prev_4kg","rf_prev_sga","bc_30d","bc_90d",
                      "bc_280d","parit_cat","m_wg_cat","log_med_hs_inc")]

##### Drawing subsamples and building the forest (in parallel) -----------------
propSample <- 1e-2 # fraction of the total number of observations used to build subsamples
total_nbtrees <- 200 # total number of trees in the forest
ncores <- 8 # number of cores
max_nodes <- 100 # maximum number of terminal nodes trees can have

moonRF <- mclapply(1:total_nbtrees, function(ind) {
  mySample = sample(1:nrow(simulated),round(propSample*nrow(simulated)), replace = FALSE)
  sampSimulated = simulated[mySample, ]
  res_bc_30d = randomForest(bc_30d ~ ., data=subset(sampSimulated, select=-c(bc_90d, bc_280d)), 
                            ntree=1, maxnodes=max_nodes,replace=FALSE, sampsize=nrow(sampSimulated))
  # res_bc_90d = randomForest(bc_90d ~., data=subset(sampSimulated, select=-c(bc_30d, bc_280d)), 
  #                           ntree=1, maxnodes=max_nodes,replace=FALSE, sampsize=nrow(sampSimulated))
  # res_bc_280d = randomForest(bc_280d ~., data=subset(sampSimulated, select=-c(bc_30d, bc_90d)), 
  #                            ntree=1, maxnodes=max_nodes,replace=FALSE, sampsize=nrow(sampSimulated))
  # return(list(forest_bc_30d=res_bc_30d,
  #             forest_bc_90d=res_bc_90d,
  #             forest_bc_280d=res_bc_280d))
  return(list(forest_bc_30d=res_bc_30d))
}, mc.cores=ncores) 

##### Predicting -----------------------------------------------------
cl = makeCluster(10,outfile='')
registerDoParallel(cl)

pred_bc_30d <- foreach(i = 1:total_nbtrees, .packages=c('randomForest'), .combine=cbind) %dopar% {
  print(i)
  pred_bc_30d_temp <- predict(moonRF[[i]]$forest_bc_30d, simulated, type="response")
  return(pred_bc_30d_temp)
  rm(pred_bc_30d_temp)
  gc()
}
pred_bc_30d <- unlist(pred_bc_30d)
pred_bc_30d <- as.data.frame(pred_bc_30d)

# pred_bc_90d <- foreach(i = 1:total_nbtrees, .packages=c('randomForest'), .combine=cbind) %dopar% {
#   print(i)
#   pred_bc_90d_temp <- predict(moonRF[[i]]$forest_bc_90d, simulated, type="response")
#   return(pred_bc_90d_temp)
#   rm(pred_bc_90d_temp)
#   gc()
# }
# pred_bc_90d <- unlist(pred_bc_90d)
# pred_bc_90d <- as.data.frame(pred_bc_90d)
# 
# pred_bc_280d <- foreach(i = 1:total_nbtrees, .packages=c('randomForest'), .combine=cbind) %dopar% {
#   print(i)
#   pred_bc_280d_temp <- predict(moonRF[[i]]$forest_bc_280d, simulated, type="response")
#   return(pred_bc_280d_temp)
#   rm(pred_bc_280d_temp)
#   gc()
# }
# pred_bc_280d <- unlist(pred_bc_280d)
# pred_bc_280d <- as.data.frame(pred_bc_280d)

stopCluster(cl)

# ##### Estimating GPS -----------------------------------------------------
birth$pred_bc_30d <- rowMeans(pred_bc_30d)
# birth$resid_bc_30d <- birth$bc_30d - birth$pred_bc_30d
# sigma_bc_30d <- sd(birth$resid_bc_30d)
# birth$gps_bc_30d <- dnorm(birth$resid_bc_30d, mean=0, sd=sigma_bc_30d)
# 
# birth$pred_bc_90d <- rowMeans(pred_bc_90d)
# birth$resid_bc_90d <- birth$bc_90d - birth$pred_bc_90d
# sigma_bc_90d <- sd(birth$resid_bc_90d)
# birth$gps_bc_90d <- dnorm(birth$resid_bc_90d, mean=0, sd=sigma_bc_90d)
# 
# birth$pred_bc_280d <- rowMeans(pred_bc_280d)
# birth$resid_bc_280d <- birth$bc_280d - birth$pred_bc_280d
# sigma_bc_280d <- sd(birth$resid_bc_280d)
# birth$gps_bc_280d <- dnorm(birth$resid_bc_280d, mean=0, sd=sigma_bc_280d)
# 
# gps_moonRF_birth <- birth[,c("uniqueid_yr", "bwg", "bc_30d", "bc_90d", "bc_280d",
#                              "gps_bc_30d", "gps_bc_90d","gps_bc_280d")]
# 
# saveRDS(gps_moonRF_birth, file = paste0(dir_results,'gps_moonRF_birth_maxnodes100_tree300.rds')) # change
# beepr::beep(3)

birth <- readRDS(paste0(dir_results,"moonRFplotCheck.rds"))
pdf(paste0(dir_results,"moonRFpred_obs.pdf"))
plot(birth$bc_30d, birth$pred_bc_30d, main="moonRF(maxnodes100_tree200)")
lines(lowess(birth$bc_30d, birth$pred_bc_30d), col="red", lwd=2)
dev.off()



