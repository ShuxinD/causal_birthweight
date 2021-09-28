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
n.trees <- 25000
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
  upper <- quantile(wt, 0.99)
  lower <- quantile(wt, 0.01)
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
write(opt_aac_bc_tri, file = paste0(dir_out_gridsearch, "GridSearchResults_bc_tri.csv"))


############################## 3.4 no2_30d ######################################
## construct data
dt <- birth[ , var, with = F]
dt[, T := no2_30d]
dt[, no2_30d := NULL]

independent <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
                 "clinega", "kotck","pncgov", "rf_db_gest","rf_db_other",
                 "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
                 "rf_prev_sga", "percentPoverty",
                 "bc_30d", "bc_3090d", "bc_90280d", 
                 # "no2_30d",
                 "no2_3090d",
                 "no2_90280d",
                 "firstborn","m_wg_cat", "smoker_ddp", "smoker_dpp",
                 "mrace_1", "mrace_2", "mrace_3", "mrace_4",
                 "log_mhincome", "log_mhvalue")

## model
model.num = lm(T~1, data = dt) 
ps.num <- dnorm((dt$T-model.num$fitted)/(summary(model.num))$sigma,0,1)

h2o.init(nthreads = n_cores, min_mem_size = "400G", port = 54345)
##### depth 6 ######
opt_aac_depth6_no2_30d <- NULL
for (rate in col.sample.rate){
  birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
  gbm_30d <- h2o.gbm(y = "T",
                     x = independent,
                     training_frame = birth.hex,
                     ntrees = n.trees, 
                     max_depth = max.depth[1], # change
                     min_rows = min.rows,
                     learn_rate = learn.rate, 
                     col_sample_rate = rate,
                     distribution = "gaussian")
  cat(paste0("predicting... with rate of ", rate, "\n"))
  pred.gbm_30d <- h2o.staged_predict_proba(object = gbm_30d, newdata = birth.hex)
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
  cat(paste0("optimizing... with rate of ", rate, "\n"))
  opt_aac_i  <- optimize(F.aac.iter, interval = c(1, n.trees), data = dt,
                         data.hex = birth.hex,
                         ps.num = ps.num, 
                         ps.model.pred = pred.gbm_30d, # change
                         rep = rep.bootstrap)
  opt_aac_i <- as.data.frame(opt_aac_i)
  opt_aac_depth6_no2_30d <- rbind(opt_aac_depth6_no2_30d, opt_aac_i)
  rm(opt_aac_i)
  h2o.removeAll()
  gc()
}
row.names(opt_aac_depth6_no2_30d) <- col.sample.rate

##### depth 8 ######
opt_aac_depth8_no2_30d <- NULL
# h2o.init(nthreads = n_cores, min_mem_size = "460G")
for (rate in col.sample.rate){
  birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
  gbm_30d <- h2o.gbm(y = "T",
                     x = independent,
                     training_frame = birth.hex,
                     ntrees = n.trees, 
                     max_depth = max.depth[2], # change
                     min_rows = min.rows,
                     learn_rate = learn.rate, 
                     col_sample_rate = rate,
                     distribution = "gaussian")
  cat(paste0("predicting... with rate of ", rate, "\n"))
  pred.gbm_30d <- h2o.staged_predict_proba(object = gbm_30d, newdata = birth.hex)
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
  cat(paste0("optimizing... with rate of ", rate, "\n"))
  opt_aac_i  <- optimize(F.aac.iter, interval = c(1, n.trees), data = dt,
                         data.hex = birth.hex,
                         ps.num = ps.num, 
                         ps.model.pred = pred.gbm_30d, # change
                         rep = rep.bootstrap)
  opt_aac_i <- as.data.frame(opt_aac_i)
  opt_aac_depth8_no2_30d <- rbind(opt_aac_depth8_no2_30d, opt_aac_i)
  rm(opt_aac_i)
  h2o.removeAll()
  gc()
}
row.names(opt_aac_depth8_no2_30d) <- col.sample.rate

###### depth 10 #######
opt_aac_depth10_no2_30d <- NULL
for (rate in col.sample.rate){
  birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
  gbm_30d <- h2o.gbm(y = "T",
                     x = independent,
                     training_frame = birth.hex,
                     ntrees = n.trees, 
                     max_depth = max.depth[3], # change
                     min_rows = min.rows,
                     learn_rate = learn.rate, 
                     col_sample_rate = rate,
                     distribution = "gaussian")
  cat(paste0("predicting... with rate of ", rate, "\n"))
  pred.gbm_30d <- h2o.staged_predict_proba(object = gbm_30d, newdata = birth.hex)
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
  cat(paste0("optimizing... with rate of ", rate, "\n"))
  opt_aac_i  <- optimize(F.aac.iter, interval = c(1, n.trees), data = dt,
                         data.hex = birth.hex,
                         ps.num = ps.num, 
                         ps.model.pred = pred.gbm_30d, # change
                         rep = rep.bootstrap)
  opt_aac_i <- as.data.frame(opt_aac_i)
  opt_aac_depth10_no2_30d <- rbind(opt_aac_depth10_no2_30d, opt_aac_i)
  rm(opt_aac_i)
  h2o.removeAll()
  gc()
}
row.names(opt_aac_depth10_no2_30d) <- col.sample.rate
h2o.shutdown(prompt = FALSE)

############################## 3.5 no2_3090d ######################################
## construct data
dt <- birth[ , var, with = F]
dt[, T := no2_3090d]
dt[, no2_3090d := NULL]

independent <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
                 "clinega", "kotck","pncgov", "rf_db_gest","rf_db_other",
                 "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
                 "rf_prev_sga", "percentPoverty",
                 "bc_30d", "bc_3090d", "bc_90280d", 
                 "no2_30d",
                 # "no2_3090d",
                 "no2_90280d",
                 "firstborn","m_wg_cat", "smoker_ddp", "smoker_dpp",
                 "mrace_1", "mrace_2", "mrace_3", "mrace_4",
                 "log_mhincome", "log_mhvalue")
## model
model.num = lm(T~1, data = dt) 
ps.num <- dnorm((dt$T-model.num$fitted)/(summary(model.num))$sigma,0,1)

h2o.init(nthreads = n_cores, min_mem_size = "410G", port = 54345)

##### depth 6 ######
opt_aac_depth6_no2_3090d <- NULL
for (rate in col.sample.rate){
  birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
  gbm_30d <- h2o.gbm(y = "T",
                     x = independent,
                     training_frame = birth.hex,
                     ntrees = n.trees, 
                     max_depth = max.depth[1], # change
                     min_rows = min.rows,
                     learn_rate = learn.rate, 
                     col_sample_rate = rate,
                     distribution = "gaussian")
  cat(paste0("predicting... with rate of ", rate, "\n"))
  pred.gbm_30d <- h2o.staged_predict_proba(object = gbm_30d, newdata = birth.hex)
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
  cat(paste0("optimizing... with rate of ", rate, "\n"))
  opt_aac_i  <- optimize(F.aac.iter, interval = c(1, n.trees), data = dt,
                         data.hex = birth.hex,
                         ps.num = ps.num, 
                         ps.model.pred = pred.gbm_30d, # change
                         rep = rep.bootstrap)
  opt_aac_i <- as.data.frame(opt_aac_i)
  opt_aac_depth6_no2_3090d <- rbind(opt_aac_depth6_no2_3090d, opt_aac_i)
  rm(opt_aac_i)
  h2o.removeAll()
  gc()
}
row.names(opt_aac_depth6_no2_3090d) <- col.sample.rate

##### depth 8 ######
opt_aac_depth8_no2_3090d <- NULL
for (rate in col.sample.rate){
  birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
  gbm_30d <- h2o.gbm(y = "T",
                     x = independent,
                     training_frame = birth.hex,
                     ntrees = n.trees, 
                     max_depth = max.depth[2], # change
                     min_rows = min.rows,
                     learn_rate = learn.rate, 
                     col_sample_rate = rate,
                     distribution = "gaussian")
  cat(paste0("predicting... with rate of ", rate, "\n"))
  pred.gbm_30d <- h2o.staged_predict_proba(object = gbm_30d, newdata = birth.hex)
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
  cat(paste0("optimizing... with rate of ", rate, "\n"))
  opt_aac_i  <- optimize(F.aac.iter, interval = c(1, n.trees), data = dt,
                         data.hex = birth.hex,
                         ps.num = ps.num, 
                         ps.model.pred = pred.gbm_30d, # change
                         rep = rep.bootstrap)
  opt_aac_i <- as.data.frame(opt_aac_i)
  opt_aac_depth8_no2_3090d <- rbind(opt_aac_depth8_no2_3090d, opt_aac_i)
  rm(opt_aac_i)
  h2o.removeAll()
  gc()
}
row.names(opt_aac_depth8_no2_3090d) <- col.sample.rate

###### depth 10 #######
opt_aac_depth10_no2_3090d <- NULL
for (rate in col.sample.rate){
  birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
  gbm_30d <- h2o.gbm(y = "T",
                     x = independent,
                     training_frame = birth.hex,
                     ntrees = n.trees, 
                     max_depth = max.depth[3], # change
                     min_rows = min.rows,
                     learn_rate = learn.rate, 
                     col_sample_rate = rate,
                     distribution = "gaussian")
  cat(paste0("predicting... with rate of ", rate, "\n"))
  pred.gbm_30d <- h2o.staged_predict_proba(object = gbm_30d, newdata = birth.hex)
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
  cat(paste0("optimizing... with rate of ", rate, "\n"))
  opt_aac_i  <- optimize(F.aac.iter, interval = c(1, n.trees), data = dt,
                         data.hex = birth.hex,
                         ps.num = ps.num, 
                         ps.model.pred = pred.gbm_30d, # change
                         rep = rep.bootstrap)
  opt_aac_i <- as.data.frame(opt_aac_i)
  opt_aac_depth10_no2_3090d <- rbind(opt_aac_depth10_no2_3090d, opt_aac_i)
  rm(opt_aac_i)
  h2o.removeAll()
  gc()
}
row.names(opt_aac_depth10_no2_3090d) <- col.sample.rate
h2o.shutdown(prompt = FALSE)

############################## 3.6 no2_90280d ######################################
## construct data
dt <- birth[ , var, with = F]
dt[, T := no2_90280d]
dt[, no2_90280d := NULL]

independent <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
                 "clinega", "kotck","pncgov", "rf_db_gest","rf_db_other",
                 "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
                 "rf_prev_sga", "percentPoverty",
                 "bc_30d", "bc_3090d", "bc_90280d", 
                 "no2_30d",
                 "no2_3090d",
                 # "no2_90280d",
                 "firstborn","m_wg_cat", "smoker_ddp", "smoker_dpp",
                 "mrace_1", "mrace_2", "mrace_3", "mrace_4",
                 "log_mhincome", "log_mhvalue")
## model
model.num = lm(T~1, data = dt) 
ps.num <- dnorm((dt$T-model.num$fitted)/(summary(model.num))$sigma,0,1)

h2o.init(nthreads = n_cores, min_mem_size = "410G", port = 54345)

##### depth 6 ######
opt_aac_depth6_no2_90280d <- NULL
for (rate in col.sample.rate){
  birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
  gbm_30d <- h2o.gbm(y = "T",
                     x = independent,
                     training_frame = birth.hex,
                     ntrees = n.trees, 
                     max_depth = max.depth[1], # change
                     min_rows = min.rows,
                     learn_rate = learn.rate, 
                     col_sample_rate = rate,
                     distribution = "gaussian")
  cat(paste0("predicting... with rate of ", rate, "\n"))
  pred.gbm_30d <- h2o.staged_predict_proba(object = gbm_30d, newdata = birth.hex)
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
  cat(paste0("optimizing... with rate of ", rate, "\n"))
  opt_aac_i  <- optimize(F.aac.iter, interval = c(1, n.trees), data = dt,
                         data.hex = birth.hex,
                         ps.num = ps.num, 
                         ps.model.pred = pred.gbm_30d, # change
                         rep = rep.bootstrap)
  opt_aac_i <- as.data.frame(opt_aac_i)
  opt_aac_depth6_no2_90280d <- rbind(opt_aac_depth6_no2_90280d, opt_aac_i)
  rm(opt_aac_i)
  h2o.removeAll()
  gc()
}
row.names(opt_aac_depth6_no2_90280d) <- col.sample.rate

##### depth 8 ######
opt_aac_depth8_no2_90280d <- NULL
for (rate in col.sample.rate){
  birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
  gbm_30d <- h2o.gbm(y = "T",
                     x = independent,
                     training_frame = birth.hex,
                     ntrees = n.trees, 
                     max_depth = max.depth[2], # change
                     min_rows = min.rows,
                     learn_rate = learn.rate, 
                     col_sample_rate = rate,
                     distribution = "gaussian")
  cat(paste0("predicting... with rate of ", rate, "\n"))
  pred.gbm_30d <- h2o.staged_predict_proba(object = gbm_30d, newdata = birth.hex)
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
  cat(paste0("optimizing... with rate of ", rate, "\n"))
  opt_aac_i  <- optimize(F.aac.iter, interval = c(1, n.trees), data = dt,
                         data.hex = birth.hex,
                         ps.num = ps.num, 
                         ps.model.pred = pred.gbm_30d, # change
                         rep = rep.bootstrap)
  opt_aac_i <- as.data.frame(opt_aac_i)
  opt_aac_depth8_no2_90280d <- rbind(opt_aac_depth8_no2_90280d, opt_aac_i)
  rm(opt_aac_i)
  h2o.removeAll()
  gc()
}
row.names(opt_aac_depth8_no2_90280d) <- col.sample.rate

###### depth 10 #######
opt_aac_depth10_no2_90280d <- NULL
for (rate in col.sample.rate){
  birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
  gbm_30d <- h2o.gbm(y = "T",
                     x = independent,
                     training_frame = birth.hex,
                     ntrees = n.trees, 
                     max_depth = max.depth[3], # change
                     min_rows = min.rows,
                     learn_rate = learn.rate, 
                     col_sample_rate = rate,
                     distribution = "gaussian")
  cat(paste0("predicting... with rate of ", rate, "\n"))
  pred.gbm_30d <- h2o.staged_predict_proba(object = gbm_30d, newdata = birth.hex)
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
  cat(paste0("optimizing... with rate of ", rate, "\n"))
  opt_aac_i  <- optimize(F.aac.iter, interval = c(1, n.trees), data = dt,
                         data.hex = birth.hex,
                         ps.num = ps.num, 
                         ps.model.pred = pred.gbm_30d, # change
                         rep = rep.bootstrap)
  opt_aac_i <- as.data.frame(opt_aac_i)
  opt_aac_depth10_no2_90280d <- rbind(opt_aac_depth10_no2_90280d, opt_aac_i)
  rm(opt_aac_i)
  h2o.removeAll()
  gc()
}
row.names(opt_aac_depth10_no2_90280d) <- col.sample.rate
h2o.shutdown(prompt = FALSE)

############################## 4. summarize results ###########################
result_bc_30d <-
  rbind(
    opt_aac_depth6_30d %>% mutate(depth = 6,col_rate = col.sample.rate),
    opt_aac_depth8_30d %>% mutate(depth = 8,col_rate = col.sample.rate),
    opt_aac_depth10_30d %>% mutate(depth = 10,col_rate = col.sample.rate))
result_bc_30d$label <- "bc_30d"

result_bc_3090d <-
  rbind(
    opt_aac_depth6_3090d %>% mutate(depth = 6,col_rate = col.sample.rate),
    opt_aac_depth8_3090d %>% mutate(depth = 8,col_rate = col.sample.rate),
    opt_aac_depth10_3090d %>% mutate(depth = 10,col_rate = col.sample.rate))
result_bc_3090d$label <- "bc_3090d"

result_bc_90280d <-
  rbind(
    opt_aac_depth6_90280d %>% mutate(depth = 6,col_rate = col.sample.rate),
    opt_aac_depth8_90280d %>% mutate(depth = 8,col_rate = col.sample.rate),
    opt_aac_depth10_90280d %>% mutate(depth = 10,col_rate = col.sample.rate))
result_bc_90280d$label <- "bc_90280d"

result_no2_30d <-
  rbind(
    opt_aac_depth6_no2_30d %>% mutate(depth = 6,col_rate = col.sample.rate),
    opt_aac_depth8_no2_30d %>% mutate(depth = 8,col_rate = col.sample.rate),
    opt_aac_depth10_no2_30d %>% mutate(depth = 10,col_rate = col.sample.rate))
result_no2_30d$label <- "no2_30d"

result_no2_3090d <-
  rbind(
    opt_aac_depth6_no2_3090d %>% mutate(depth = 6,col_rate = col.sample.rate),
    opt_aac_depth8_no2_3090d %>% mutate(depth = 8,col_rate = col.sample.rate),
    opt_aac_depth10_no2_3090d %>% mutate(depth = 10,col_rate = col.sample.rate))
result_no2_3090d$label <- "no2_3090d"

result_no2_90280d <-
  rbind(
    opt_aac_depth6_no2_90280d %>% mutate(depth = 6,col_rate = col.sample.rate),
    opt_aac_depth8_no2_90280d %>% mutate(depth = 8,col_rate = col.sample.rate),
    opt_aac_depth10_no2_90280d %>% mutate(depth = 10,col_rate = col.sample.rate))
result_no2_90280d$label <- "no2_90280d"

results <- rbindlist(list(result_bc_30d, result_bc_3090d, result_bc_90280d,
               result_no2_30d, result_no2_3090d, result_no2_90280d))

colnames(results) <- c("optimal_ntrees", "AAC", "depth", "col_rate", "label")
setDT(results)
results

fwrite(results, paste0(dir_gridsearch, "GridSearchResults.csv"))

results <- fread(paste0(dir_gridsearch, "GridSearchResults.csv"))
best <- results[, .(bestAAC = min(AAC)), by = label]
best <- merge(best, results, by.x = c("label", "bestAAC"), by.y = c("label", "AAC"))
best[, optimal_ntrees := round(optimal_ntrees)][]
fwrite(best, paste0(dir_gridsearch, "BestCombinations_gbm.csv"))

gc()

#######
class(opt_aac_depth6_bc_all)
opt_aac_depth6_bc_all$depth <- 6
opt_aac_depth8_bc_all$depth <- 8
opt_aac_depth10_bc_all$depth <- 10

opt_aac_depth6_bc_all$col_rate <- rownames(opt_aac_depth6_bc_all)
opt_aac_depth8_bc_all$col_rate <- rownames(opt_aac_depth8_bc_all)
opt_aac_depth10_bc_all$col_rate <- rownames(opt_aac_depth10_bc_all)

opt_aac_depth6_bc_all$label <- "bc_all"
opt_aac_depth8_bc_all$label <- "bc_all"
opt_aac_depth10_bc_all$label <- "bc_all"

results_bc_all <- rbind(opt_aac_depth6_bc_all, opt_aac_depth8_bc_all, opt_aac_depth10_bc_all)

dir_gridsearch <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/1GridSearchResults/"
setDT(results_bc_all)
fwrite(results_bc_all, file = paste0(dir_gridsearch, "GridSearchResults_bc_all.csv"))

class(opt_mac_depth6_bc_all)
opt_mac_depth6_bc_all$depth <- 6
opt_mac_depth8_bc_all$depth <- 8
opt_mac_depth10_bc_all$depth <- 10

opt_mac_depth6_bc_all$col_rate <- rownames(opt_mac_depth6_bc_all)
opt_mac_depth8_bc_all$col_rate <- rownames(opt_mac_depth8_bc_all)
opt_mac_depth10_bc_all$col_rate <- rownames(opt_mac_depth10_bc_all)

opt_mac_depth6_bc_all$label <- "bc_all"
opt_mac_depth8_bc_all$label <- "bc_all"
opt_mac_depth10_bc_all$label <- "bc_all"

results_bc_all_mac <- rbind(opt_mac_depth6_bc_all, opt_mac_depth8_bc_all, opt_mac_depth10_bc_all)

dir_gridsearch <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/1GridSearchResults/"
setDT(results_bc_all_mac)
fwrite(results_bc_all_mac, file = paste0(dir_gridsearch, "GridSearchResults_bc_all_mac.csv"))

#######
class(opt_aac_depth6_no2_all)
opt_aac_depth6_no2_all$depth <- 6
opt_aac_depth8_no2_all$depth <- 8
opt_aac_depth10_no2_all$depth <- 10

opt_aac_depth6_no2_all$col_rate <- rownames(opt_aac_depth6_no2_all)
opt_aac_depth8_no2_all$col_rate <- rownames(opt_aac_depth8_no2_all)
opt_aac_depth10_no2_all$col_rate <- rownames(opt_aac_depth10_no2_all)

opt_aac_depth6_no2_all$label <- "no2_all"
opt_aac_depth8_no2_all$label <- "no2_all"
opt_aac_depth10_no2_all$label <- "no2_all"

results_no2_all <- rbind(opt_aac_depth6_no2_all, opt_aac_depth8_no2_all, opt_aac_depth10_no2_all)

dir_gridsearch <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/1GridSearchResults/"
setDT(results_no2_all)
fwrite(results_no2_all, file = paste0(dir_gridsearch, "GridSearchResults_no2_all.csv"))

class(opt_mac_depth6_no2_all)
opt_mac_depth6_no2_all$depth <- 6
opt_mac_depth8_no2_all$depth <- 8
opt_mac_depth10_no2_all$depth <- 10

opt_mac_depth6_no2_all$col_rate <- rownames(opt_mac_depth6_no2_all)
opt_mac_depth8_no2_all$col_rate <- rownames(opt_mac_depth8_no2_all)
opt_mac_depth10_no2_all$col_rate <- rownames(opt_mac_depth10_no2_all)

opt_mac_depth6_no2_all$label <- "no2_all"
opt_mac_depth8_no2_all$label <- "no2_all"
opt_mac_depth10_no2_all$label <- "no2_all"

results_no2_all_mac <- rbind(opt_mac_depth6_no2_all, opt_mac_depth8_no2_all, opt_mac_depth10_no2_all)

dir_gridsearch <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/1GridSearchResults/"
setDT(results_no2_all_mac)
fwrite(results_no2_all_mac, file = paste0(dir_gridsearch, "GridSearchResults_no2_all_mac.csv"))
