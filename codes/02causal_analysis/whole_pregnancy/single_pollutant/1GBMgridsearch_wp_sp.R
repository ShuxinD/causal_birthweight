#' Project: causal_birthweight
#' Code: GBM to estimate GPS with balance stopping rule
#' Input: "MAbirth_for_analyses.csv"
#' Output: "GridSearchResults_wholepregnancy.csv"
#' Output: "BestCombinations_gbm.csv"
#' Author: Shuxin Dong
#' First Ceate Date: 2021-01-21

## 0. Setup -----
rm(list = ls())
gc()

library(h2o)
library(dplyr)
library(data.table)
library(polycor)
library(funique)
library(parallel)
library(wCorr)
n_cores <- detectCores() - 1 

setwd("/media/gate/Shuxin/")
dir_in <- "/media/qnap3/Shuxin/bc_no2_MAbirth_causal/data/"
dir_gridsearch <- "/media/qnap3/Shuxin/bc_no2_MAbirth_causal/github_repo/results/1GridSearchResults/whole_pregnancy/two_pollutant/"

# 0.1 hyperparameter range ----
## set parameters for h2o.gbm model
min.rows <- 10 
learn.rate <- 0.01
## hyperparameters range for h2o.gbm model
n.trees <- 1000
max.depth <- c(3,4,5,6,7)
col.sample.rate <- c(0.8,0.9,1.0)

## 0.2  functions ----
truncate_ipw <- function(ipw_raw, upper_bound_percentile, lower_bound_percentile){
  #' ipw_raw: raw stabilized ipw
  #' upper_bound_percentile: truncation upper bound limit usually 0.99
  up_bound <- quantile(ipw_raw, upper_bound_percentile)
  low_bound <- quantile(ipw_raw, lower_bound_percentile)
  condition <- rep(0, length(ipw_raw))
  condition[ipw_raw>up_bound] <- 1
  condition[ipw_raw<low_bound] <- -1
  ipw <- ifelse(condition==1, up_bound, ifelse(condition==-1, low_bound, ipw_raw))
  return(ipw) # output truncated ipw
}

F.aac.iter <- function(col_i, data, ipw_num, ps.model.pred) {
  #' col_i: number of iterations (number of trees) 
  #' data: dataset containing the treatment and the covariates not in h2o structure.
  #' ps.model.pred: the staged prediction results of boosting model to estimate (p(T_iX_i)) 
  #' ipw_num: the estimated p(E_i) 
  GBM.fitted <- as.numeric(as.vector(ps.model.pred[, floor(col_i)]))
  GPS <- dnorm((data$T - GBM.fitted)/sd(data$T - GBM.fitted), 0, 1)
  ipw_raw <- ipw_num/GPS
  # upper <- quantile(wt, 0.9992)
  # lower <- quantile(wt, 0.0008)
  ipw <- truncate_ipw(ipw_raw, 0.975, 0.025)
  
  x <- data %>% select(-T)
  ac <- matrix(NA,dim(x)[2],1)
  x <- as.data.frame(x)
  for (j in 1:dim(x)[2]){
    ac[j] = ifelse (!is.factor(x[,j]), 
                    weightedCorr(data$T, x[,j], method = "pearson", weights = ipw), 
                    weightedCorr(data$T, x[,j], method = "polyserial", weights = ipw))
  }
  aac <- mean(abs(ac), na.rm = TRUE)
  return(aac)
}

## data manipulation ----
birth <- fread(paste0(dir_in, "MAbirth_for_analyses.csv"))
names(birth)
birth$year <- as.factor(birth$year)

ps_exposures <- c("bc_all", "no2_all")
ps_vars <- c("year","sex","married","mage", # "cigdpp","cigddp",
             "clinega","pncgov", 
             "rf_db_gest","rf_db_other", "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg", "rf_prev_sga", 
             "smoker_ddp", # "smoker_dpp",
             "mrace_1", "mrace_2", "mrace_3", # "mrace_4",
             "m_edu_1", "m_edu_2", "m_edu_3","m_edu_4", # "m_edu_5",
             "kotck_1","kotck_2","kotck_3","kotck_4",
             "m_wg_cat_1","m_wg_cat_2","m_wg_cat_3","m_wg_cat_4",# "m_wg_cat_5",
             # "log_mhvalue", "log_mhincome",
             "percentPoverty", "firstborn")
summary(birth[,..ps_vars])

## bc_all ----
exposure <- "bc_all"

h2o.init(nthreads = n_cores, min_mem_size = "250G", port = 54345)
cat("FITTING ", exposure, "\n")
dt <- birth[ , c(ps_vars, exposure) , with = F]
dt[, T := get(exposure)]
dt[, (exposure) := NULL]
independent <- c(ps_vars)
model.num = lm(T~1, data = dt) 
ps.num <- dnorm((dt$T-model.num$fitted)/(summary(model.num))$sigma,0,1) # p(E_i)

opt_aac <- NULL
for (max.depth_i in max.depth){
  for (col.sample.rate_i in col.sample.rate){
    opt_aac_i <- NULL
    birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
    gbm_exposure <- h2o.gbm(y = "T",
                            x = independent,
                            training_frame = birth.hex,
                            ntrees = n.trees, 
                            max_depth = max.depth_i,
                            min_rows = min.rows,
                            learn_rate = learn.rate, 
                            col_sample_rate = col.sample.rate_i,
                            distribution = "gaussian")
    cat(paste0("predicting... with rate of ", col.sample.rate_i, "depth ", max.depth_i, "\n"))
    pred.gbm_exposure <- h2o.staged_predict_proba(object = gbm_exposure, newdata = birth.hex)
    h2o:::.h2o.garbageCollect()
    h2o:::.h2o.garbageCollect()
    h2o:::.h2o.garbageCollect()
    cat(paste0("optimizing... with rate of ", col.sample.rate_i, "depth ", max.depth_i, "\n"))
    opt_aac_i  <- optimize(F.aac.iter, interval = c(50, ncol(pred.gbm_exposure)), data = dt,
                           ipw_num = ps.num, 
                           ps.model.pred = pred.gbm_exposure)
    cat(paste0("finished... with rate of ", col.sample.rate_i, "depth ", max.depth_i, "\n"))
    print(opt_aac_i)
    h2o:::.h2o.garbageCollect()
    opt_aac_i <- as.data.frame(opt_aac_i)
    opt_aac_i$rate <- col.sample.rate_i
    opt_aac_i$depth <- max.depth_i
    opt_aac_i$ntree <- round(as.numeric(opt_aac_i$minimum))
    row.names(opt_aac_i) <- paste0(exposure, "depth", max.depth_i, "rate", col.sample.rate_i)
    opt_aac <- rbind(opt_aac, opt_aac_i)
    # rm(opt_aac_i)
    h2o.removeAll()
    gc()
  }
}
print(opt_aac)
write.csv(opt_aac, paste0(dir_gridsearch, "GridSearch_bc_all_single.csv"))
h2o.shutdown(prompt = FALSE)

## no2_all ----
exposure <- "no2_all"

h2o.init(nthreads = n_cores, min_mem_size = "250G", port = 54345)
cat("FITTING ", exposure, "\n")
dt <- birth[ , c(ps_vars, exposure) , with = F]
dt[, T := get(exposure)]
dt[, (exposure) := NULL]
independent <- c(ps_vars)
model.num = lm(T~1, data = dt) 
ps.num <- dnorm((dt$T-model.num$fitted)/(summary(model.num))$sigma,0,1) # p(E_i)

opt_aac <- NULL
for (max.depth_i in max.depth){
  for (col.sample.rate_i in col.sample.rate){
    opt_aac_i <- NULL
    birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
    gbm_exposure <- h2o.gbm(y = "T",
                            x = independent,
                            training_frame = birth.hex,
                            ntrees = n.trees, 
                            max_depth = max.depth_i,
                            min_rows = min.rows,
                            learn_rate = learn.rate, 
                            col_sample_rate = col.sample.rate_i,
                            distribution = "gaussian")
    cat(paste0("predicting... with rate of ", col.sample.rate_i, "depth ", max.depth_i, "\n"))
    pred.gbm_exposure <- h2o.staged_predict_proba(object = gbm_exposure, newdata = birth.hex)
    h2o:::.h2o.garbageCollect()
    h2o:::.h2o.garbageCollect()
    h2o:::.h2o.garbageCollect()
    cat(paste0("optimizing... with rate of ", col.sample.rate_i, "depth ", max.depth_i, "\n"))
    opt_aac_i  <- optimize(F.aac.iter, interval = c(50, ncol(pred.gbm_exposure)), data = dt,
                           ipw_num = ps.num, 
                           ps.model.pred = pred.gbm_exposure)
    cat(paste0("finished... with rate of ", col.sample.rate_i, "depth ", max.depth_i, "\n"))
    print(opt_aac_i)
    h2o:::.h2o.garbageCollect()
    opt_aac_i <- as.data.frame(opt_aac_i)
    opt_aac_i$rate <- col.sample.rate_i
    opt_aac_i$depth <- max.depth_i
    opt_aac_i$ntree <- round(as.numeric(opt_aac_i$minimum))
    row.names(opt_aac_i) <- paste0(exposure, "depth", max.depth_i, "rate", col.sample.rate_i)
    opt_aac <- rbind(opt_aac, opt_aac_i)
    # rm(opt_aac_i)
    h2o.removeAll()
    gc()
  }
}
print(opt_aac)
write.csv(opt_aac, paste0(dir_gridsearch, "GridSearch_no2_all_single.csv"))
h2o.shutdown(prompt = FALSE)


