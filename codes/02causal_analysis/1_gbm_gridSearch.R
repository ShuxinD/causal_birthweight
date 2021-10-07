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
n_cores <- detectCores() - 1 

setwd("/media/gate/Shuxin/")
dir_in <- "/media/qnap3/Shuxin/airPollution_MAbirth/"
dir_out_gridsearch <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/1GridSearchResults/"

#' set parameters for h2o.gbm model
min.rows <- 10
learn.rate <- 0.005
rep.bootstrap <- 50

#' set hyperparameters range for h2o.gbm model
n.trees <- 12000
max.depth <- c(3,4,5,6,8)
col.sample.rate <- c(0.8, 0.9, 1.0)

## 1. functions to compute AAC ----
#' change gbm package to h2o.gbm packge
#'  use data.table to manipulate the data instead of default
#'  use subsample of bootstrapping instead of whole sample 
F.aac.iter <- function(i, data, data.hex, ps.num, ps.model.pred, rep) {
  # i: number of iterations (number of trees) 
  # data: dataset containing the treatment and the covariates not in h2o structure.
  # data.hex: dataset containing the treatment and the covariates in h2o env.
  # ps.model.pred: the staged prediction results of boosting model to estimate (p(T_iX_i)) 
  # ps.num: the estimated p(T_i) 
  # rep: number of replications in bootstrap 
  GBM.fitted <- as.vector(ps.model.pred[, floor(i)])
  ps.den <- dnorm((data$T - GBM.fitted)/sd(data$T - GBM.fitted), 0, 1)
  wt <- ps.num/ps.den
  upper <- quantile(wt, 0.9992)
  lower <- quantile(wt, 0.0008)
  wt <- fifelse(wt>upper, upper, wt)
  wt <- fifelse(wt<lower, lower, wt)
  aac_iter <- rep(NA,rep) 
  for (i in 1:rep){
    # set.seed(i)
    bo <- sample(1:dim(data)[1], size = floor((dim(data)[1])^0.7), 
                 replace = TRUE, prob = wt) # change
    newsample <- data[bo,]
    newsample <- Filter(function(x)(length(funique(x))>1), newsample)
    # j.drop <- match(c("T"), names(data))
    # j.drop <- j.drop[!is.na(j.drop)]
    x <- newsample %>% select(-T)
    ac <- matrix(NA,dim(x)[2],1)
    x <- as.data.frame(x)
    for (j in 1:dim(x)[2]){
      ac[j] = ifelse (!is.factor(x[,j]), stats::cor(newsample$T, x[,j],
                                                    method = "pearson"),
                      polyserial(newsample$T, x[,j]))
    }
    aac_iter[i] <- mean(abs(ac), na.rm = TRUE)
  }
  aac <- mean(aac_iter)
  return(aac)
}

F.mac.iter <- function(i, data, data.hex, ps.num, ps.model.pred, rep) {
  # i: number of iterations (number of trees) 
  # data: dataset containing the treatment and the covariates not in h2o structure.
  # data.hex: dataset containing the treatment and the covariates in h2o env.
  # ps.model.pred: the staged prediction results of boosting model to estimate (p(T_iX_i)) 
  # ps.num: the estimated p(T_i) 
  # rep: number of replications in bootstrap 
  GBM.fitted <- as.vector(ps.model.pred[, floor(i)])
  ps.den <- dnorm((data$T - GBM.fitted)/sd(data$T - GBM.fitted), 0, 1)
  wt <- ps.num/ps.den
  upper <- quantile(wt, 0.99)
  lower <- quantile(wt, 0.01)
  wt <- fifelse(wt>upper, upper, wt)
  wt <- fifelse(wt<lower, lower, wt)
  mac_iter <- rep(NA,rep) 
  for (i in 1:rep){
    # set.seed(i)
    bo <- sample(1:dim(data)[1], size = floor((dim(data)[1])^0.7), 
                 replace = TRUE, prob = wt) # change
    newsample <- data[bo,]
    newsample <- Filter(function(x)(length(funique(x))>1), newsample)
    # j.drop <- match(c("T"), names(data))
    # j.drop <- j.drop[!is.na(j.drop)]
    x <- newsample %>% select(-T)
    ac <- matrix(NA,dim(x)[2],1)
    x <- as.data.frame(x)
    for (j in 1:dim(x)[2]){
      ac[j] = ifelse (!is.factor(x[,j]), stats::cor(newsample$T, x[,j],
                                                    method = "pearson"),
                      polyserial(newsample$T, x[,j]))
    }
    mac_iter[i] <- max(abs(ac), na.rm = TRUE)
  }
  mac <- mean(mac_iter)
  return(mac)
}

## 2. data manipulation ----
birth <- fread(paste0(dir_in, "MAbirth_for_analyses.csv"))
names(birth)
# [1] "uniqueid_yr"    "year"           "sex"            "married"        "mage"           "mrace"         
# [7] "m_edu"          "cigdpp"         "cigddp"         "clinega"        "kotck"          "pncgov"        
# [13] "bdob"           "bwg"            "rf_db_gest"     "rf_db_other"    "rf_hbp_chronic" "rf_hbp_pregn"  
# [19] "rf_cervix"      "rf_prev_4kg"    "rf_prev_sga"    "mhincome"       "mhvalue"        "percentPoverty"
# [25] "bc_30d"         "bc_3090d"       "bc_90280d"      "no2_30d"        "no2_3090d"      "no2_90280d"    
# [31] "bc_all"         "no2_all"        "lbw"            "firstborn"      "m_wg_cat"       "smoker_ddp"    
# [37] "smoker_dpp"     "mrace_1"        "mrace_2"        "mrace_3"        "mrace_4"        "log_mhincome"  
# [43] "log_mhvalue"    "b_spring"       "b_summer"       "b_autumn"       "b_winter" 
birth$year <- as.factor(birth$year)
birth$m_edu <- as.factor(birth$m_edu)
birth$kotck <- as.factor(birth$kotck)
birth$m_wg_cat <- as.factor(birth$m_wg_cat)

ps_exposures <- c("bc_30d","bc_3090d", "bc_90280d", 
                  "no2_30d", "no2_3090d", "no2_90280d")
ps_vars <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
         "clinega", "kotck","pncgov", "rf_db_gest","rf_db_other",
         "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
         "rf_prev_sga", "percentPoverty",
         "firstborn","m_wg_cat", "smoker_ddp", "smoker_dpp",
         "mrace_1", "mrace_2", "mrace_3", "mrace_4",
         "log_mhincome", "log_mhvalue",
         "b_spring", "b_summer", "b_autumn", "b_winter")

## 3. grid search for whole pregnancy period exposure ----
#' run each chunk seperately for each exposure of interest
## 3.01 bc_all, with no2_all ----
dt <- birth[ , c(ps_vars, "bc_all", "no2_all"), with = F]
dt[, T := bc_all]
dt[, bc_all := NULL]
independent <- c(ps_vars, "no2_all")
model.num = lm(T~1, data = dt) 
ps.num <- dnorm((dt$T-model.num$fitted)/(summary(model.num))$sigma,0,1)

h2o.init(nthreads = n_cores, min_mem_size = "400G", port = 54345)
cat("fitting for bc_all\n")
opt_aac_bc_all <- NULL
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
    opt_aac_i  <- optimize(F.aac.iter, interval = c(800, n.trees), data = dt,
                           data.hex = birth.hex,
                           ps.num = ps.num, 
                           ps.model.pred = pred.gbm_exposure, # change
                           rep = rep.bootstrap)
    opt_aac_i <- as.data.frame(opt_aac_i)
    row.names(opt_aac_i) <- paste0("depth", max.depth_i, "rate", col.sample.rate_i)
    opt_aac_bc_all <- rbind(opt_aac_bc_all, opt_aac_i)
    rm(opt_aac_i)
    h2o.removeAll()
    gc()
  }
}
h2o.shutdown(prompt = FALSE)
opt_aac_bc_all
write.csv(opt_aac_bc_all, file = paste0(dir_out_gridsearch, "GridSearchResults_bc_all.csv"))

## 3.02 no2_all, with bc_all ----
dt <- birth[ , c(ps_vars, "no2_all", "bc_all"), with = F]
dt[, T := no2_all]
dt[, no2_all := NULL]
independent <- c(ps_vars, "bc_all")
model.num = lm(T~1, data = dt) 
ps.num <- dnorm((dt$T-model.num$fitted)/(summary(model.num))$sigma,0,1)

h2o.init(nthreads = n_cores, min_mem_size = "450G", port = 54345)
cat("fitting for no2_all\n")
opt_aac_no2_all <- NULL
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
    opt_aac_i  <- optimize(F.aac.iter, interval = c(800, n.trees), data = dt,
                           data.hex = birth.hex,
                           ps.num = ps.num, 
                           ps.model.pred = pred.gbm_exposure, # change
                           rep = rep.bootstrap)
    opt_aac_i <- as.data.frame(opt_aac_i)
    row.names(opt_aac_i) <- paste0("depth", max.depth_i, "rate", col.sample.rate_i)
    opt_aac_no2_all <- rbind(opt_aac_no2_all, opt_aac_i)
    rm(opt_aac_i)
    h2o.removeAll()
    gc()
  }
}
h2o.shutdown(prompt = FALSE)
opt_aac_no2_all
write.csv(opt_aac_no2_all, file = paste0(dir_out_gridsearch, "GridSearchResults_no2_all.csv"))

## 3.1 three periods of BC with no2_... ----
bc_exposures <- c("bc_30d", "bc_3090d", "bc_90280d")
opt_aac_bc_tri <- NULL
for (bc_exposures_i in bc_exposures) {
  cat("FITTING ", bc_exposures_i, "\n")
  dt <- birth[ , c(ps_vars, ps_exposures) , with = F]
  dt[, T := get(bc_exposures_i)]
  dt[, (bc_exposures_i) := NULL]
  independent <- c(ps_vars, ps_exposures[ps_exposures!=bc_exposures_i])
  model.num = lm(T~1, data = dt) 
  ps.num <- dnorm((dt$T-model.num$fitted)/(summary(model.num))$sigma,0,1)
  
  h2o.init(nthreads = n_cores, min_mem_size = "400G", port = 54345)
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
      opt_aac_i  <- optimize(F.aac.iter, interval = c(800, n.trees), data = dt,
                             data.hex = birth.hex,
                             ps.num = ps.num, 
                             ps.model.pred = pred.gbm_exposure, # change
                             rep = rep.bootstrap)
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
  rm(opt_aac_bc_tri_i)
}
h2o.shutdown(prompt = FALSE)
opt_aac_bc_tri
write.csv(opt_aac_bc_tri, file = paste0(dir_out_gridsearch, "GridSearchResults_bc_tri.csv"))

## 3.2 three periods of NO2 with bc_... ----
no2_exposures <- c("no2_30d", "no2_3090d", "no2_90280d")
opt_aac_no2_tri <- NULL
for (no2_exposures_i in no2_exposures) {
  cat("FITTING ", no2_exposures_i, "\n")
  dt <- birth[ , c(ps_vars, ps_exposures) , with = F]
  dt[, T := get(no2_exposures_i)]
  dt[, (no2_exposures_i) := NULL]
  independent <- c(ps_vars, ps_exposures[ps_exposures!=no2_exposures_i])
  model.num = lm(T~1, data = dt) 
  ps.num <- dnorm((dt$T-model.num$fitted)/(summary(model.num))$sigma,0,1)
  
  h2o.init(nthreads = n_cores, min_mem_size = "400G", port = 54345)
  cat("fitting for ", no2_exposures_i, "\n")
  opt_aac_no2_tri_i <- NULL
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
      opt_aac_i  <- optimize(F.aac.iter, interval = c(800, n.trees), data = dt,
                             data.hex = birth.hex,
                             ps.num = ps.num, 
                             ps.model.pred = pred.gbm_exposure, # change
                             rep = rep.bootstrap)
      opt_aac_i <- as.data.frame(opt_aac_i)
      row.names(opt_aac_i) <- paste0(no2_exposures_i, "depth", max.depth_i, "rate", col.sample.rate_i)
      opt_aac_no2_tri_i <- rbind(opt_aac_no2_tri_i, opt_aac_i)
      rm(opt_aac_i)
      h2o.removeAll()
      gc()
    }
  }
  opt_aac_no2_tri <- rbind(opt_aac_no2_tri, opt_aac_no2_tri_i)
  rm(opt_aac_no2_tri_i)
}
h2o.shutdown(prompt = FALSE)
opt_aac_no2_tri
write.csv(opt_aac_no2_tri, file = paste0(dir_out_gridsearch, "GridSearchResults_no2_tri.csv"))