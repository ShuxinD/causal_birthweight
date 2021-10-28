#' Project: causal_birthweight
#' Code: GBM to estimate GPS with balance stopping rule
#' Code : from three non-overlapping periods
#' Input: "MAbirth_for_analyses.csv"
#' Output: "GridSearchResults.csv"
#' Output: "BestCombinations_gbm.csv"
#' Author: Shuxin Dong
#' Create Date: 2021-01-21

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
dir_in <- "/media/qnap3/Shuxin/airPollution_MAbirth/"
dir_out_gridsearch_temp <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/1GridSearchResults/"

# 0.1 hyperparameter range ----
## set parameters for h2o.gbm model
min.rows <- 10 
learn.rate <- 0.01
rep.num <- 100
## hyperparameters range for h2o.gbm model
n.trees <- 4000
max.depth <- c(4,5,6,8)
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

F.aac.iter <- function(col_i, data, ipw_num, ps.model.pred, rep) {
  #' col_i: number of iterations (number of trees) 
  #' data: dataset containing the treatment and the covariates not in h2o structure.
  #' ps.model.pred: the staged prediction results of boosting model to estimate (p(T_iX_i)) 
  #' ipw_num: the estimated p(E_i) 
  #' rep: number of replications in bootstrap 
  GBM.fitted <- as.vector(ps.model.pred[, floor(col_i)])
  GPS <- dnorm((data$T - GBM.fitted)/sd(data$T - GBM.fitted), 0, 1)
  ipw_raw <- ipw_num/GPS
  # upper <- quantile(wt, 0.9992)
  # lower <- quantile(wt, 0.0008)
  # ipw <- truncate_ipw(ipw_raw, 0.999, 0.001)
  ipw <- ipw_raw
  aac_iter <- rep(NA,rep) 
  for (i in 1:rep){
    set.seed(i)
    bo <- sample(1:dim(data)[1], size = floor((dim(data)[1])^0.7), replace = TRUE, prob = ipw)
    newsample <- data[bo,]
    # newsample <- Filter(function(x)(length(funique(x))>1), newsample)
    # j.drop <- match(c("T"), names(data))
    # j.drop <- j.drop[!is.na(j.drop)]
    x <- newsample %>% select(-T)
    ac <- matrix(NA,dim(x)[2],1)
    x <- as.data.frame(x)
    for (j in 1:dim(x)[2]){
      ac[j] = ifelse (!is.factor(x[,j]), stats::cor(newsample$T, x[,j], method = "pearson"), polyserial(newsample$T, x[,j]))
    }
    aac_iter[i] <- mean(abs(ac), na.rm = TRUE)
  }
  aac <- mean(aac_iter)
  return(aac)
}

## data manipulation ----
birth <- fread(paste0(dir_in, "MAbirth_for_analyses.csv"))
names(birth)
birth$year <- as.factor(birth$year)

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

## sample run bc_30d----
opt_aac <- NULL

bc_exposures_i <- "bc_90280d" # set sample

cat("FITTING ", bc_exposures_i, "\n")
dt <- birth[ , c(ps_vars, ps_exposures) , with = F]
dt[, T := get(bc_exposures_i)]
dt[, (bc_exposures_i) := NULL]
independent <- c(ps_vars, ps_exposures[ps_exposures!=bc_exposures_i])
model.num = lm(T~1, data = dt) 
ps.num <- dnorm((dt$T-model.num$fitted)/(summary(model.num))$sigma,0,1)

h2o.init(nthreads = n_cores, min_mem_size = "250G", port = 54345)
cat("fitting for ", bc_exposures_i, "\n")
opt_aac_bc_tri_i <- NULL
# for (max.depth_i in max.depth){
max.depth_i <- 8
cat("fitting depth ", max.depth_i, "\n")
# for (col.sample.rate_i in col.sample.rate){
col.sample.rate_i <- 1
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
cat(paste0("predicting... with rate of ", col.sample.rate_i, "\n"))
pred.gbm_exposure <- h2o.staged_predict_proba(object = gbm_exposure, newdata = birth.hex)
h2o:::.h2o.garbageCollect()
h2o:::.h2o.garbageCollect()
h2o:::.h2o.garbageCollect()
cat(paste0("optimizing... with rate of ", col.sample.rate_i, "\n"))
opt_aac_i  <- optimize(F.aac.iter, interval = c(500, ncol(pred.gbm_exposure)), data = dt,
                       ipw_num = ps.num, 
                       ps.model.pred = pred.gbm_exposure,
                       rep = rep.num)
opt_aac_i <- as.data.frame(opt_aac_i)
opt_aac_i$rate <- col.sample.rate_i
opt_aac_i$depth <- max.depth_i
opt_aac_i
row.names(opt_aac_i) <- paste0(bc_exposures_i, "depth", max.depth_i, "rate", col.sample.rate_i)
opt_aac <- rbind(opt_aac, opt_aac_i)
rm(opt_aac_i)
opt_aac
opt_aac$exposure <- c("bc_30d", "bc_3090d", "no2_30d", "no2_3090d", "no2_90280d", "bc_90280d")
opt_aac$ntree <- round(opt_aac$minimum)
fwrite(opt_aac, file=paste0(dir_out_gridsearch_temp, "GridSearchResults1028.csv"))

h2o.shutdown(prompt = FALSE)

## bc three periods----
bc_exposures <- c("bc_30d", "bc_3090d", "bc_90280d")
h2o.init(nthreads = n_cores, min_mem_size = "250G", port = 54345)
opt_aac_bc_tri <- NULL
for (bc_exposures_i in bc_exposures) {
  cat("FITTING ", bc_exposures_i, "\n")
  dt <- birth[ , c(ps_vars, ps_exposures) , with = F]
  dt[, T := get(bc_exposures_i)]
  dt[, (bc_exposures_i) := NULL]
  independent <- c(ps_vars, ps_exposures[ps_exposures!=bc_exposures_i])
  model.num = lm(T~1, data = dt) 
  ps.num <- dnorm((dt$T-model.num$fitted)/(summary(model.num))$sigma,0,1) # p(E_i)
  
  cat("fitting for ", bc_exposures_i, "\n")
  opt_aac_bc_tri_i <- NULL
  for (max.depth_i in max.depth){
    cat("fitting depth ", max.depth_i, "\n")
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
      cat(paste0("predicting... with rate of ", col.sample.rate_i, "\n"))
      pred.gbm_exposure <- h2o.staged_predict_proba(object = gbm_exposure, newdata = birth.hex)
      h2o:::.h2o.garbageCollect()
      h2o:::.h2o.garbageCollect()
      h2o:::.h2o.garbageCollect()
      cat(paste0("optimizing... with rate of ", col.sample.rate_i, "\n"))
      opt_aac_i  <- optimize(F.aac.iter, interval = c(800, ncol(pred.gbm_exposure)), data = dt,
                             ipw_num = ps.num, 
                             ps.model.pred = pred.gbm_exposure,
                             rep = rep.num)
      opt_aac_i
      h2o:::.h2o.garbageCollect()
      opt_aac_i <- as.data.frame(opt_aac_i)
      opt_aac_i$rate <- col.sample.rate_i
      opt_aac_i$depth <- max.depth_i
      row.names(opt_aac_i) <- paste0(bc_exposures_i, "depth", max.depth_i, "rate", col.sample.rate_i)
      opt_aac_bc_tri_i <- rbind(opt_aac_bc_tri_i, opt_aac_i)
      rm(opt_aac_i)
      h2o.removeAll()
      gc()
    }
  }
  opt_aac_bc_tri <- rbind(opt_aac_bc_tri, opt_aac_bc_tri_i)
  write.csv(opt_aac_bc_tri_i, paste0(dir_out_gridsearch_temp, "GridSearch", bc_exposures_i, ".csv"))
  rm(opt_aac_bc_tri_i)
}
opt_aac_bc_tri

h2o.shutdown(prompt = FALSE)

## no2 three periods ----
no2_exposures <- c("no2_30d", "no2_3090d", "no2_90280d")
h2o.init(nthreads = n_cores, min_mem_size = "250G", port = 54345)
opt_aac_bc_tri <- NULL
for (no2_exposures_i in no2_exposures) {
  cat("FITTING ", no2_exposures, "\n")
  dt <- birth[ , c(ps_vars, ps_exposures) , with = F]
  dt[, T := get(no2_exposures)]
  dt[, (no2_exposures) := NULL]
  independent <- c(ps_vars, ps_exposures[ps_exposures!=no2_exposures])
  model.num = lm(T~1, data = dt) 
  ps.num <- dnorm((dt$T-model.num$fitted)/(summary(model.num))$sigma,0,1) # p(E_i)
  
  cat("fitting for ", no2_exposures, "\n")
  opt_aac_bc_tri_i <- NULL
  for (max.depth_i in max.depth){
    cat("fitting depth ", max.depth_i, "\n")
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
      cat(paste0("predicting... with rate of ", col.sample.rate_i, "\n"))
      pred.gbm_exposure <- h2o.staged_predict_proba(object = gbm_exposure, newdata = birth.hex)
      h2o:::.h2o.garbageCollect()
      h2o:::.h2o.garbageCollect()
      h2o:::.h2o.garbageCollect()
      cat(paste0("optimizing... with rate of ", col.sample.rate_i, "\n"))
      opt_aac_i  <- optimize(F.aac.iter, interval = c(800, ncol(pred.gbm_exposure)), data = dt,
                             ipw_num = ps.num, 
                             ps.model.pred = pred.gbm_exposure,
                             rep = rep.num)
      opt_aac_i
      h2o:::.h2o.garbageCollect()
      opt_aac_i <- as.data.frame(opt_aac_i)
      opt_aac_i$rate <- col.sample.rate_i
      opt_aac_i$depth <- max.depth_i
      row.names(opt_aac_i) <- paste0(bc_exposures_i, "depth", max.depth_i, "rate", col.sample.rate_i)
      opt_aac_bc_tri_i <- rbind(opt_aac_bc_tri_i, opt_aac_i)
      rm(opt_aac_i)
      h2o.removeAll()
      gc()
    }
  }
  opt_aac_bc_tri <- rbind(opt_aac_bc_tri, opt_aac_bc_tri_i)
  write.csv(opt_aac_bc_tri_i, paste0(dir_out_gridsearch_temp, "GridSearch", bc_exposures_i, ".csv"))
  rm(opt_aac_bc_tri_i)
}
opt_aac_bc_tri

h2o.shutdown(prompt = FALSE)
