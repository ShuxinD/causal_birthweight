#' Project: causal_birthweight
#' Code: get IPW and balance plot
#' Input: "MAbirth_for_analyses.csv" and grid search results
#' Output: "MAbirth_ipw.csv"
#' Author: Shuxin Dong
#' Create Date: 2021-09-07

## 0. Setup -----
rm(list = ls())
gc()

library(data.table)
library(polycor)
library(h2o)
library(parallel)
n_cores <- detectCores() - 1 

setwd("/media/gate/Shuxin")
dir_in_birth <- "/media/qnap3/Shuxin/airPollution_MAbirth/"
dir_in_gridsearch <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/1GridSearchResults/"
dir_out_ipwplots <- "/media/gate/Shuxin/airPollution_MAbirth/causal_birthweight/results/2ipwPlots/"

#' set default parameters for H2O
min.rows <- 10
learn.rate <- 0.01

#‘ hyperparameter search range
n.trees <- 1000
max.depth <- c(3,4,5,6,8)
col.sample.rate <- c(0.8, 0.9, 1.0)

## 01. data related ----
birth <- fread(paste0(dir_in_birth, "MAbirth_for_analyses.csv"))
# names(birth)
birth$year <- as.factor(birth$year)
birth$m_edu <- as.factor(birth$m_edu)
birth$kotck <- as.factor(birth$kotck)
birth$m_wg_cat <- as.factor(birth$m_wg_cat)

birth <- fastDummies::dummy_cols(birth, select_columns = c("m_edu", "kotck","m_wg_cat"))
names(birth)
# [1] "uniqueid_yr"    "year"           "sex"            "married"        "mage"           "mrace"          "m_edu"          "cigdpp"        
# [9] "cigddp"         "clinega"        "kotck"          "pncgov"         "bdob"           "bwg"            "rf_db_gest"     "rf_db_other"   
# [17] "rf_hbp_chronic" "rf_hbp_pregn"   "rf_cervix"      "rf_prev_4kg"    "rf_prev_sga"    "mhincome"       "mhvalue"        "percentPoverty"
# [25] "bc_30d"         "bc_3090d"       "bc_90280d"      "no2_30d"        "no2_3090d"      "no2_90280d"     "bc_all"         "no2_all"       
# [33] "lbw"            "firstborn"      "m_wg_cat"       "smoker_ddp"     "smoker_dpp"     "mrace_1"        "mrace_2"        "mrace_3"       
# [41] "mrace_4"        "log_mhincome"   "log_mhvalue"    "b_spring"       "b_summer"       "b_autumn"       "b_winter"       "m_edu_1"       
# [49] "m_edu_2"        "m_edu_3"        "m_edu_4"        "m_edu_5"        "kotck_1"        "kotck_2"        "kotck_3"        "kotck_4"       
# [57] "m_wg_cat_1"     "m_wg_cat_2"     "m_wg_cat_3"     "m_wg_cat_4"     "m_wg_cat_5"  

ps_exposures <- c("bc_30d","bc_3090d", "bc_90280d", 
                  "no2_30d", "no2_3090d", "no2_90280d")
ps_vars <- c("year","sex","married","mage", "cigdpp","cigddp",
             "clinega","pncgov", 
             # "rf_db_gest","rf_db_other", "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg", "rf_prev_sga", 
             # "smoker_ddp", "smoker_dpp",
             "mrace_1", "mrace_2", "mrace_3", "mrace_4",
             "m_edu_1", "m_edu_2", "m_edu_3","m_edu_4", "m_edu_5",
             "kotck_1","kotck_2","kotck_3","kotck_4",
             "m_wg_cat_1","m_wg_cat_2","m_wg_cat_3","m_wg_cat_4","m_wg_cat_5",
             # "log_mhvalue", "log_mhincome",
             "percentPoverty", "firstborn")

## 02. specify response ----
response <- "bc_30d" 

## 1.1 setup grid-searched hyperparameters ----
gs <- read.csv(paste0(dir_in_gridsearch,"GridSearchResults_",response,".csv"))
names(gs)
# [1] "X"         "minimum"   "objective"
gs[gs$objective==min(gs$objective),]
# X  minimum  objective
# 12 bc_30ddepth6rate1 5099.688 0.06034888
n.trees <- 5100 # change
max.depth <- 6 # change
rate <- 1 # change

## 1.2 get GPS ----
predictor <- c(ps_vars, ps_exposures[ps_exposures!=response])
h2o.init(nthreads = n_cores, min_mem_size = "200G", port = 54345)
dt <- birth[ , c(ps_exposures,ps_vars), with = F]
dt[, T := get(response)]
dt[, (response) := NULL]

## fit GBM
birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
cat("fitting..grandient boosting model... \n")
gbm <- h2o.gbm(y = "T",
               x = predictor,
               training_frame = birth.hex,
               ntrees = n.trees,
               max_depth = max.depth,
               col_sample_rate = rate,
               min_rows = min.rows,
               learn_rate = learn.rate, 
               distribution = "gaussian")
cat("predicting the response...")
pred_gbm <- h2o.predict(object = gbm, newdata = birth.hex)
GBM.fitted <- as.vector(pred_gbm) # GBM.fitted_30d
h2o.removeAll()
h2o:::.h2o.garbageCollect()
h2o:::.h2o.garbageCollect()
h2o:::.h2o.garbageCollect()
h2o.shutdown(prompt = F)

## 2.1 functions ----
generate_ipw <- function(model_fitted_value, birth_raw_data, exposure_of_interest){
  numerator <- dnorm((birth_raw_data[,get(exposure_of_interest)]-mean(birth_raw_data[,get(exposure_of_interest)]))/sd(birth_raw_data[,get(exposure_of_interest)]),0,1)
  gps <- dnorm((birth_raw_data[,get(exposure_of_interest)]-model_fitted_value)/sd(birth_raw_data[,get(exposure_of_interest)]-model_fitted_value),0,1)
  ipw <- numerator/gps
  return(ipw)
}

truncate_ipw <- function(ipw_raw, upper_bound_percentile, lower_bound_percentile){
  up_bound <- quantile(ipw_raw, upper_bound_percentile)
  low_bound <- quantile(ipw_raw, lower_bound_percentile)
  condition <- rep(0, length(ipw_raw))
  condition[ipw_raw>up_bound] <- 1
  condition[ipw_raw<low_bound] <- -1
  ipw <- ifelse(condition==1, up_bound, ifelse(condition==-1, low_bound, ipw_raw))
  return(ipw)
}

## 2.2 get IPW ----
ipw_raw <- generate_ipw(GBM.fitted, birth, "bc_30d")
summary(ipw_raw)
hist(ipw_raw)

## truncate
ipw.9992 <- truncate_ipw(ipw_raw, 0.9992, 0.0008)
summary(ipw.9992)

## 3. balance ----
ac <- matrix(NA,length(predictor),2)
for (j in 1:length(predictor)){ ## unweighted
  ac[j,1] <-  ifelse(!is.factor(birth[,..predictor][[j]]), 
                     cor(birth[,get(response)], birth[,..predictor][[j]], method = "pearson"),
                     polyserial(birth[,get(response)], birth[,..predictor][[j]]))
  ac[j,2] <-  ifelse(!is.factor(birth[,..predictor][[j]]), 
                     weightedCorr(birth[,get(response)], birth[,..predictor][[j]], method = "pearson", weights = ipw.9992),
                     weightedCorr(birth[,get(response)], birth[,..predictor][[j]], method = "polyserial", weights = ipw.9992))
  cat("finish column", j,"/", length(predictor), "\n")
}

balance.9992 <- data.frame(name.x=predictor,
                           corr=c(ac[,1],ac[,2]),
                           weighted = c(rep("unweighted", length(predictor)),rep("weighted", length(predictor))))

plot.9992 <- ggplot(balance.9992, aes(x = name.x, y = corr, group = weighted)) +
  # geom_line(aes(colour = weighted), size = 0.5, linetype = "dashed")+
  geom_point(aes(colour = weighted), size = 2) +
  geom_hline(yintercept=0, size=0.2) +
  geom_hline(yintercept=-0.1, size=0.2, linetype = "dashed") +
  geom_hline(yintercept=0.1, size=0.2, linetype = "dashed") +
  geom_hline(yintercept=-0.2, size=0.1, linetype = "dashed") +
  geom_hline(yintercept=0.2, size=0.1, linetype = "dashed") +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  coord_flip() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        text = element_text(size=16))
plot.9992
