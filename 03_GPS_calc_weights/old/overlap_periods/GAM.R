##### Setting up -------------------------------------------------------------------------
rm(list=ls())
gc()

library(readr)
library(mgcv)

dir_data <- "/Users/dongshuxin/Desktop/MA_births/Reformat_Descrp_Tfm/" 
dir_results <- "/Users/dongshuxin/Desktop/MA_births/" 

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
# sample <- birth[1:2000, ]
##### GAM model ---------------------------------------------------------------
mod_gam_bc_30d <- gam(bc_30d ~ s(mage,fx=FALSE,bs='cr') + s(cigdpp,fx=FALSE,bs='cr') + 
                        s(cigddp,fx=FALSE,bs='cr') + s(clinega,k=4, fx=FALSE,bs='cr') + 
                        s(log_med_hs_inc,fx=FALSE,bs='cr') +
                        year + sex + married + mrace + m_edu + kotck + pncgov + 
                        rf_db_gest + rf_db_other + rf_hbp_chronic + rf_hbp_pregn + rf_cervix + rf_prev_4kg + rf_prev_sga +
                        parit_cat + m_wg_cat, data = birth, family = gaussian())
pred_bc_30d <- predict.gam(mod_gam_bc_30d, birth)
pred_bc_30d <- as.numeric(pred_bc_30d)
birth$pred_bc_30d <- pred_bc_30d
saveRDS(birth, paste0(dir_results,"GAM.rds"))
