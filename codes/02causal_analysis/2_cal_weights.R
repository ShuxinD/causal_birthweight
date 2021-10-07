#' Project: Causal black carbon and NO2 on birth weight in MA
#' Code: Estimate IPW with hyperparameter we set
#' Input: "MAbirth_for_analyses.csv" birth data
#' Input: grid search results
#' Output: weights dataset "MAbirth_ipw.csv" with ID and weights
#' Author: Shuxin Dong
#' First create date: 2021-01-27

## 0. setup ----
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

## set default parameters for H2O
min.rows <- 10
learn.rate <- 0.005

## read grid search results
gs_bc_all <- read.csv(paste0(dir_in_gridsearch, "GridSearchResults_bc_all.csv"))
setDT(gs_bc_all)
gs_bc_all[objective==min(objective),]
best_bc_all <- data.table(depth = 5,
                          col_rate = 0.9,
                          n.trees = 19960)

gs_no2_all <- read.csv(paste0(dir_in_gridsearch, "GridSearchResults_no2_all.csv"))
setDT(gs_no2_all)
gs_no2_all[objective==min(objective),]
best_no2_all <- data.table(depth = 8,
                          col_rate = 0.8,
                          n.trees = 13782)

gs_bc_tri <- read.csv(paste0(dir_in_gridsearch, "GridSearchResults_bc_tri.csv"))
setDT(gs_bc_tri)
gs_bc_tri[1:15, ][objective==min(objective),]
best_bc_30d <- data.table(depth = 4,
                          col_rate = 0.8,
                          n.trees = 15450)
gs_bc_tri[16:30, ][objective==min(objective),]
best_bc_3090d <- data.table(depth = 8,
                          col_rate = 1,
                          n.trees = 11329)
gs_bc_tri[31:45, ][objective==min(objective),]
best_bc_90280d <- data.table(depth = 8,
                            col_rate = 1,
                            n.trees = 17105)
# [1]  "n.trees"  "depth"    "col_rate"

# ## read bc_all and no2_all grid search results
# best <- data.frame(label = c("bc_all", "bc_all_mac", "no2_all", "no2_all_mac"),
#            bestAC = c(0.095428055, 0.988065078, 0.092660732, 0.947540067),
#            n.trees = c(8283, 6101, 9156, 14513),
#            depth = c(6, 6, 10, 6),
#            col_rate = c(0.8, 0.9, 0.9, 0.9))
# setDT(best)

########################## 1. Load birth data #################################
## load data
birth <- fread(paste0(dir_in_birth, "MAbirth_for_analyses.csv"))
names(birth)
# [1] "uniqueid_yr"    "year"           "sex"            "married"        "mage"           "mrace"         
# [7] "m_edu"          "cigdpp"         "cigddp"         "clinega"        "kotck"          "pncgov"        
# [13] "bdob"           "bwg"            "rf_db_gest"     "rf_db_other"    "rf_hbp_chronic" "rf_hbp_pregn"  
# [19] "rf_cervix"      "rf_prev_4kg"    "rf_prev_sga"    "mhincome"       "mhvalue"        "percentPoverty"
# [25] "bc_30d"         "bc_3090d"       "bc_90280d"      "no2_30d"        "no2_3090d"      "no2_90280d"    
# [31] "bc_all"         "no2_all"        "lbw"            "firstborn"      "m_wg_cat"       "smoker_ddp"    
# [37] "smoker_dpp"     "mrace_1"        "mrace_2"        "mrace_3"        "mrace_4"        "log_mhincome"  
# [43] "log_mhvalue"    "b_spring"       "b_summer"       "b_autumn"       "b_winter" 
# summary(birth)
birth$year <- as.factor(birth$year)
birth$m_edu <- as.factor(birth$m_edu)
birth$kotck <- as.factor(birth$kotck)
birth$m_wg_cat <- as.factor(birth$m_wg_cat)
# birth <- sample_n(birth, floor(0.001*dim(birth)[1])) # sample as an example

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

############################# 2. Fit GBM ######################################
h2o.init(nthreads = n_cores, min_mem_size = "200G", port = 54345)

## 2.011 bc_all ----
## construct data
dt <- birth[ , var, with = F]
dt[, T := bc_all] # change
dt[, bc_all := NULL] # change

independent <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
                 "clinega", "kotck","pncgov", "rf_db_gest","rf_db_other",
                 "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
                 "rf_prev_sga", "percentPoverty",
                 # "bc_all",
                 "no2_all",
                 "firstborn","m_wg_cat", "smoker_ddp", "smoker_dpp",
                 "mrace_1", "mrace_2", "mrace_3", "mrace_4",
                 "log_mhincome", "log_mhvalue") # change

model_best <- best[label=="bc_all",] #change

## fit GBM
birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
cat("fitting..bc_all...")
gbm <- h2o.gbm(y = "T",
               x = independent,
               training_frame = birth.hex,
               ntrees = model_best$n.trees,
               max_depth = model_best$depth,
               col_sample_rate = model_best$col_rate,
               min_rows = min.rows,
               learn_rate = learn.rate, 
               distribution = "gaussian")
cat("predicting...bc_all...")
pred_gbm <- h2o.predict(object = gbm, newdata = birth.hex)
GBM.fitted_bc_all <- as.vector(pred_gbm) # GBM.fitted_30d
h2o.removeAll()
h2o:::.h2o.garbageCollect()
h2o:::.h2o.garbageCollect()
h2o:::.h2o.garbageCollect()

## 2.012 no2_all ----
## construct data
dt <- birth[ , var, with = F]
dt[, T := no2_all] # change
dt[, no2_all := NULL] # change

independent <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
                 "clinega", "kotck","pncgov", "rf_db_gest","rf_db_other",
                 "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
                 "rf_prev_sga", "percentPoverty",
                 "bc_all",
                 # "no2_all",
                 "firstborn","m_wg_cat", "smoker_ddp", "smoker_dpp",
                 "mrace_1", "mrace_2", "mrace_3", "mrace_4",
                 "log_mhincome", "log_mhvalue") # change

model_best <- best[label=="no2_all",] #change

## fit GBM
birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
cat("fitting..no2_all...")
gbm <- h2o.gbm(y = "T",
               x = independent,
               training_frame = birth.hex,
               ntrees = model_best$n.trees,
               max_depth = model_best$depth,
               col_sample_rate = model_best$col_rate,
               min_rows = min.rows,
               learn_rate = learn.rate, 
               distribution = "gaussian")
cat("predicting...no2_all...")
pred_gbm <- h2o.predict(object = gbm, newdata = birth.hex)
GBM.fitted_no2_all <- as.vector(pred_gbm) # GBM.fitted_30d
h2o.removeAll()
h2o:::.h2o.garbageCollect()
h2o:::.h2o.garbageCollect()
h2o:::.h2o.garbageCollect()

## 2.021 bc_all_mac ----
## construct data
dt <- birth[ , var, with = F]
dt[, T := bc_all] # change
dt[, bc_all := NULL] # change

independent <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
                 "clinega", "kotck","pncgov", "rf_db_gest","rf_db_other",
                 "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
                 "rf_prev_sga", "percentPoverty",
                 # "bc_all",
                 "no2_all",
                 "firstborn","m_wg_cat", "smoker_ddp", "smoker_dpp",
                 "mrace_1", "mrace_2", "mrace_3", "mrace_4",
                 "log_mhincome", "log_mhvalue") # change

model_best <- best[label=="bc_all_mac",] #change

## fit GBM
birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
cat("fitting..bc_all_mac...")
gbm <- h2o.gbm(y = "T",
               x = independent,
               training_frame = birth.hex,
               ntrees = model_best$n.trees,
               max_depth = model_best$depth,
               col_sample_rate = model_best$col_rate,
               min_rows = min.rows,
               learn_rate = learn.rate, 
               distribution = "gaussian")
cat("predicting...bc_all_mac...")
pred_gbm <- h2o.predict(object = gbm, newdata = birth.hex)
GBM.fitted_bc_all_mac <- as.vector(pred_gbm) # GBM.fitted_30d
h2o.removeAll()
h2o:::.h2o.garbageCollect()
h2o:::.h2o.garbageCollect()
h2o:::.h2o.garbageCollect()

## 2.012 no2_all_mac ----
## construct data
dt <- birth[ , var, with = F]
dt[, T := no2_all] # change
dt[, no2_all := NULL] # change

independent <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
                 "clinega", "kotck","pncgov", "rf_db_gest","rf_db_other",
                 "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
                 "rf_prev_sga", "percentPoverty",
                 "bc_all",
                 # "no2_all",
                 "firstborn","m_wg_cat", "smoker_ddp", "smoker_dpp",
                 "mrace_1", "mrace_2", "mrace_3", "mrace_4",
                 "log_mhincome", "log_mhvalue") # change

model_best <- best[label=="no2_all_mac",] #change

## fit GBM
birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
cat("fitting..no2_all_mac...")
gbm <- h2o.gbm(y = "T",
               x = independent,
               training_frame = birth.hex,
               ntrees = model_best$n.trees,
               max_depth = model_best$depth,
               col_sample_rate = model_best$col_rate,
               min_rows = min.rows,
               learn_rate = learn.rate, 
               distribution = "gaussian")
cat("predicting...no2_all_mac...")
pred_gbm <- h2o.predict(object = gbm, newdata = birth.hex)
GBM.fitted_no2_all_mac <- as.vector(pred_gbm) # GBM.fitted_30d
h2o.removeAll()
h2o:::.h2o.garbageCollect()
h2o:::.h2o.garbageCollect()
h2o:::.h2o.garbageCollect()

h2o.shutdown(prompt = FALSE)
############################## 2.1 bc_30d ######################################
## construct data
dt <- birth[ , c(ps_exposures,ps_vars), with = F]
dt[, T := bc_30d] # change
dt[, bc_30d := NULL] # change

independent <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
                 "clinega", "kotck","pncgov", "rf_db_gest","rf_db_other",
                 "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
                 "rf_prev_sga", "percentPoverty",
                 # "bc_30d",
                 "bc_3090d", "bc_90280d", 
                 "no2_30d", "no2_3090d", "no2_90280d",
                 "firstborn","m_wg_cat", "smoker_ddp", "smoker_dpp",
                 "mrace_1", "mrace_2", "mrace_3", "mrace_4",
                 "log_mhincome", "log_mhvalue") # change

## fit GBM
birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
cat("fitting..bc_30d...")
gbm <- h2o.gbm(y = "T",
               x = independent,
               training_frame = birth.hex,
               ntrees = best_bc_30d$n.trees,
               max_depth = best_bc_30d$depth,
               col_sample_rate = best_bc_30d$col_rate,
               min_rows = min.rows,
               learn_rate = learn.rate, 
               distribution = "gaussian")
cat("predicting...bc_30d...")
pred_gbm <- h2o.predict(object = gbm, newdata = birth.hex)
GBM.fitted_bc_30d <- as.vector(pred_gbm) # GBM.fitted_30d
h2o.removeAll()
h2o:::.h2o.garbageCollect()
h2o:::.h2o.garbageCollect()
h2o:::.h2o.garbageCollect()

############################## 2.2 bc_3090d ######################################
## construct data
dt <- birth[ , var, with = F]
dt[, T := bc_3090d] # change
dt[, bc_3090d := NULL] # change

independent <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
                 "clinega", "kotck","pncgov", "rf_db_gest","rf_db_other",
                 "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
                 "rf_prev_sga", "percentPoverty",
                 "bc_30d",
                 # "bc_3090d",
                 "bc_90280d", 
                 "no2_30d", "no2_3090d", "no2_90280d",
                 "firstborn","m_wg_cat", "smoker_ddp", "smoker_dpp",
                 "mrace_1", "mrace_2", "mrace_3", "mrace_4",
                 "log_mhincome", "log_mhvalue",
                 "b_spring", "b_summer", "b_autumn", "b_winter") # change

## fit GBM
birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
cat("fitting..bc_3090d...")
gbm <- h2o.gbm(y = "T",
               x = independent,
               training_frame = birth.hex,
               ntrees = best_bc_3090d$n.trees,
               max_depth = best_bc_3090d$depth,
               col_sample_rate = best_bc_3090d$col_rate,
               min_rows = min.rows,
               learn_rate = learn.rate, 
               distribution = "gaussian")
cat("predicting...bc_3090d...")
pred_gbm <- h2o.predict(object = gbm, newdata = birth.hex)
GBM.fitted_bc_3090d <- as.vector(pred_gbm) # GBM.fitted_30d
h2o.removeAll()
h2o:::.h2o.garbageCollect()
h2o:::.h2o.garbageCollect()
h2o:::.h2o.garbageCollect()

############################## 2.3 bc_90280d ######################################
## construct data
dt <- birth[ , var, with = F]
dt[, T := bc_90280d] # change
dt[, bc_90280d := NULL] # change

independent <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
                 "clinega", "kotck","pncgov", "rf_db_gest","rf_db_other",
                 "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
                 "rf_prev_sga", "percentPoverty",
                 "bc_30d",
                 "bc_3090d",
                 # "bc_90280d", 
                 "no2_30d", "no2_3090d", "no2_90280d",
                 "firstborn","m_wg_cat", "smoker_ddp", "smoker_dpp",
                 "mrace_1", "mrace_2", "mrace_3", "mrace_4",
                 "log_mhincome", "log_mhvalue", 
                 "b_spring", "b_summer", "b_autumn", "b_winter") # change


## fit GBM
birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
cat("fitting..bc_90280d...")
gbm <- h2o.gbm(y = "T",
               x = independent,
               training_frame = birth.hex,
               ntrees = best_bc_90280d$n.trees,
               max_depth = best_bc_90280d$depth,
               col_sample_rate = best_bc_90280d$col_rate,
               min_rows = min.rows,
               learn_rate = learn.rate, 
               distribution = "gaussian")
cat("predicting...bc_90280d...")
pred_gbm <- h2o.predict(object = gbm, newdata = birth.hex)
GBM.fitted_bc_90280d <- as.vector(pred_gbm) # GBM.fitted_30d
h2o.removeAll()
h2o:::.h2o.garbageCollect()
h2o:::.h2o.garbageCollect()
h2o:::.h2o.garbageCollect()

############################## 2.4 no2_30d ######################################
## construct data
dt <- birth[ , var, with = F]
dt[, T := no2_30d] # change
dt[, no2_30d := NULL] # change

independent <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
                 "clinega", "kotck","pncgov", "rf_db_gest","rf_db_other",
                 "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
                 "rf_prev_sga", "percentPoverty",
                 "bc_30d", "bc_3090d", "bc_90280d", 
                 # "no2_30d",
                 "no2_3090d", "no2_90280d",
                 "firstborn","m_wg_cat", "smoker_ddp", "smoker_dpp",
                 "mrace_1", "mrace_2", "mrace_3", "mrace_4",
                 "log_mhincome", "log_mhvalue") # change

model_best <- best[label=="no2_30d",] #change

## fit GBM
birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
cat("fitting..no2_30d...")
gbm <- h2o.gbm(y = "T",
               x = independent,
               training_frame = birth.hex,
               ntrees = model_best$n.trees,
               max_depth = model_best$depth,
               col_sample_rate = model_best$col_rate,
               min_rows = min.rows,
               learn_rate = learn.rate, 
               distribution = "gaussian")
cat("predicting...no2_30d...")
pred_gbm <- h2o.predict(object = gbm, newdata = birth.hex)
GBM.fitted_no2_30d <- as.vector(pred_gbm) # GBM.fitted_30d
h2o.removeAll()
h2o:::.h2o.garbageCollect()
h2o:::.h2o.garbageCollect()
h2o:::.h2o.garbageCollect()

############################## 2.5 no2_3090d ##################################
## construct data
dt <- birth[ , var, with = F]
dt[, T := no2_3090d] # change
dt[, no2_3090d := NULL] # change

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
                 "log_mhincome", "log_mhvalue") # change

model_best <- best[label=="no2_3090d",] #change

## fit GBM
birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
cat("fitting..no2_3090d...")
gbm <- h2o.gbm(y = "T",
               x = independent,
               training_frame = birth.hex,
               ntrees = model_best$n.trees,
               max_depth = model_best$depth,
               col_sample_rate = model_best$col_rate,
               min_rows = min.rows,
               learn_rate = learn.rate, 
               distribution = "gaussian")
cat("predicting...no2_3090d...")
pred_gbm <- h2o.predict(object = gbm, newdata = birth.hex)
GBM.fitted_no2_3090d <- as.vector(pred_gbm) # GBM.fitted_30d
h2o.removeAll()
h2o:::.h2o.garbageCollect()
h2o:::.h2o.garbageCollect()
h2o:::.h2o.garbageCollect()

############################## 2.6 no2_90280d ##################################
## construct data
dt <- birth[ , var, with = F]
dt[, T := no2_90280d] # change
dt[, no2_90280d := NULL] # change

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
                 "log_mhincome", "log_mhvalue") # change

model_best <- best[label=="no2_90280d",] #change

## fit GBM
birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
cat("fitting..no2_90280d...")
gbm <- h2o.gbm(y = "T",
               x = independent,
               training_frame = birth.hex,
               ntrees = model_best$n.trees,
               max_depth = model_best$depth,
               col_sample_rate = model_best$col_rate,
               min_rows = min.rows,
               learn_rate = learn.rate, 
               distribution = "gaussian")
cat("predicting...no2_90280d...")
pred_gbm <- h2o.predict(object = gbm, newdata = birth.hex)
GBM.fitted_no2_90280d <- as.vector(pred_gbm) # GBM.fitted_30d
h2o.removeAll()
h2o:::.h2o.garbageCollect()
h2o:::.h2o.garbageCollect()
h2o:::.h2o.garbageCollect()

h2o.shutdown(prompt = FALSE)
gc()

############################# 3. Calculate Weights ############################
## bc_all
model.num_bc_all = lm(bc_all ~ 1, data = birth) 
ps.num_bc_all <- dnorm((birth$bc_all - model.num_bc_all$fitted)/(summary(model.num_bc_all))$sigma,0,1)
# scale.res <- scale(residuals(model.num_bc_all))
# df <- approxfun(density(scale.res, kernel = "gaussian", n = dim(birth)[1]))
# ps.num_bc_all <- df(res)
ps.den_bc_all <- dnorm((birth$bc_all - GBM.fitted_bc_all)/sd(birth$bc_all - GBM.fitted_bc_all),0,1)

bc_all.wt <- ps.num_bc_all/ps.den_bc_all
summary(bc_all.wt)

max.wt <- quantile(bc_all.wt, 0.975)
min.wt <- quantile(bc_all.wt, 0.025)
bc_all.wt <- fifelse(bc_all.wt>max.wt, max.wt, bc_all.wt) # truncate
bc_all.wt <- fifelse(bc_all.wt<min.wt, min.wt, bc_all.wt) # truncate
summary(bc_all.wt)

## no2_all
model.num_no2_all = lm(no2_all ~ 1, data = birth) 
ps.num_no2_all <- dnorm((birth$no2_all - model.num_no2_all$fitted)/(summary(model.num_no2_all))$sigma,0,1)
# scale.res <- scale(residuals(model.num_no2_all))
# df <- approxfun(density(scale.res, kernel = "gaussian", n = dim(birth)[1]))
# ps.num_no2_all <- df(res)
ps.den_no2_all <- dnorm((birth$no2_all - GBM.fitted_no2_all)/sd(birth$no2_all - GBM.fitted_no2_all),0,1)

no2_all.wt <- ps.num_no2_all/ps.den_no2_all
summary(no2_all.wt)

max.wt <- quantile(no2_all.wt, 0.975)
min.wt <- quantile(no2_all.wt, 0.025)
no2_all.wt <- fifelse(no2_all.wt>max.wt, max.wt, no2_all.wt) # truncate
no2_all.wt <- fifelse(no2_all.wt<min.wt, min.wt, no2_all.wt) # truncate
summary(no2_all.wt)

## bc_all_mac
model.num_bc_all = lm(bc_all ~ 1, data = birth) 
ps.num_bc_all <- dnorm((birth$bc_all - model.num_bc_all$fitted)/(summary(model.num_bc_all))$sigma,0,1)
# scale.res <- scale(residuals(model.num_bc_all))
# df <- approxfun(density(scale.res, kernel = "gaussian", n = dim(birth)[1]))
# ps.num_bc_all <- df(res)
ps.den_bc_all_mac <- dnorm((birth$bc_all - GBM.fitted_bc_all_mac)/sd(birth$bc_all - GBM.fitted_bc_all_mac),0,1)

bc_all_mac.wt <- ps.num_bc_all/ps.den_bc_all_mac
summary(bc_all_mac.wt)

max.wt <- quantile(bc_all_mac.wt, 0.975)
min.wt <- quantile(bc_all_mac.wt, 0.025)
bc_all_mac.wt <- fifelse(bc_all_mac.wt>max.wt, max.wt, bc_all_mac.wt) # truncate
bc_all_mac.wt <- fifelse(bc_all_mac.wt<min.wt, min.wt, bc_all_mac.wt) # truncate
summary(bc_all_mac.wt)

## no2_all_mac
model.num_no2_all = lm(no2_all ~ 1, data = birth) 
ps.num_no2_all <- dnorm((birth$no2_all - model.num_no2_all$fitted)/(summary(model.num_no2_all))$sigma,0,1)
# scale.res <- scale(residuals(model.num_no2_all))
# df <- approxfun(density(scale.res, kernel = "gaussian", n = dim(birth)[1]))
# ps.num_no2_all <- df(res)
ps.den_no2_all_mac <- dnorm((birth$no2_all - GBM.fitted_no2_all_mac)/sd(birth$no2_all - GBM.fitted_no2_all_mac),0,1)

no2_all_mac.wt <- ps.num_no2_all/ps.den_no2_all_mac
summary(no2_all_mac.wt)

max.wt <- quantile(no2_all_mac.wt, 0.975)
min.wt <- quantile(no2_all_mac.wt, 0.025)
no2_all_mac.wt <- fifelse(no2_all_mac.wt>max.wt, max.wt, no2_all_mac.wt) # truncate
no2_all_mac.wt <- fifelse(no2_all_mac.wt<min.wt, min.wt, no2_all_mac.wt) # truncate
summary(no2_all_mac.wt)

## bc_30d
model.num_bc_30d = lm(bc_30d ~ 1, data = birth) 
ps.num_bc_30d <- dnorm((birth$bc_30d - model.num_bc_30d$fitted)/(summary(model.num_bc_30d))$sigma,0,1)
# scale.res <- scale(residuals(model.num_bc_30d))
# df <- approxfun(density(scale.res, kernel = "gaussian", n = dim(birth)[1]))
# ps.num_bc_30d <- df(res)
ps.den_bc_30d <- dnorm((birth$bc_30d - GBM.fitted_bc_30d)/sd(birth$bc_30d - GBM.fitted_bc_30d),0,1)

bc_30d.wt <- ps.num_bc_30d/ps.den_bc_30d
summary(bc_30d.wt)

max.wt <- quantile(bc_30d.wt, 0.975)
min.wt <- quantile(bc_30d.wt, 0.025)
bc_30d.wt <- fifelse(bc_30d.wt>max.wt, max.wt, bc_30d.wt) # truncate
bc_30d.wt <- fifelse(bc_30d.wt<min.wt, min.wt, bc_30d.wt) # truncate
summary(bc_30d.wt)

## bc_3090d
model.num_bc_3090d = lm(bc_3090d ~ 1, data = birth) 

ps.num_bc_3090d <- dnorm((birth$bc_3090d - model.num_bc_3090d$fitted)/(summary(model.num_bc_3090d))$sigma,0,1)
#df <- approxfun(density(birth[, bc_3090d], kernel = "gaussian", n = dim(birth)[1]))
#ps.num_bc_3090d <- df(birth[,bc_3090d])
ps.den_bc_3090d <- dnorm((birth$bc_3090d - GBM.fitted_bc_3090d)/sd(birth$bc_3090d - GBM.fitted_bc_3090d),0,1)

bc_3090d.wt <- ps.num_bc_3090d/ps.den_bc_3090d
summary(bc_3090d.wt)

max.wt <- quantile(bc_3090d.wt, 0.975)
min.wt <- quantile(bc_3090d.wt, 0.025)
bc_3090d.wt <- fifelse(bc_3090d.wt>max.wt, max.wt, bc_3090d.wt) # truncate
bc_3090d.wt <- fifelse(bc_3090d.wt<min.wt, min.wt, bc_3090d.wt) # truncate
summary(bc_3090d.wt)

## bc_90280d
model.num_bc_90280d = lm(bc_90280d ~ 1, data = birth) 

ps.num_bc_90280d <- dnorm((birth$bc_90280d - model.num_bc_90280d$fitted)/(summary(model.num_bc_90280d))$sigma,0,1)
# df <- approxfun(density(birth[, bc_90280d], kernel = "gaussian", n = dim(birth)[1]))
# ps.num_bc_90280d <- df(birth[,bc_90280d])
ps.den_bc_90280d <- dnorm((birth$bc_90280d - GBM.fitted_bc_90280d)/sd(birth$bc_90280d - GBM.fitted_bc_90280d),0,1)

bc_90280d.wt <- ps.num_bc_90280d/ps.den_bc_90280d
summary(bc_90280d.wt)

max.wt <- quantile(bc_90280d.wt, 0.975)
min.wt <- quantile(bc_90280d.wt, 0.025)
bc_90280d.wt <- fifelse(bc_90280d.wt>max.wt, max.wt, bc_90280d.wt) # truncate
bc_90280d.wt <- fifelse(bc_90280d.wt<min.wt, min.wt, bc_90280d.wt) # truncate
summary(bc_90280d.wt)
# sd(bc_90280d.wt)/sqrt(dim(birth)[1])

## no2_30d
model.num_no2_30d = lm(no2_30d ~ 1, data = birth) 
hist(birth[,no2_30d])
hist((birth$no2_30d - model.num_no2_30d$fitted)/(summary(model.num_no2_30d))$sigma)
ps.num_no2_30d <- dnorm((birth$no2_30d - model.num_no2_30d$fitted)/(summary(model.num_no2_30d))$sigma,0,1)
# df <- approxfun(density(birth[, no2_30d], kernel = "gaussian", n = dim(birth)[1]))
# ps.num_no2_30d <- df(birth[,no2_30d])
ps.den_no2_30d <- dnorm((birth$no2_30d - GBM.fitted_no2_30d)/sd(birth$no2_30d - GBM.fitted_no2_30d),0,1)

no2_30d.wt <- ps.num_no2_30d/ps.den_no2_30d
summary(no2_30d.wt)

max.wt <- quantile(no2_30d.wt, 0.975)
min.wt <- quantile(no2_30d.wt, 0.025)
no2_30d.wt <- fifelse(no2_30d.wt>max.wt, max.wt, no2_30d.wt) # truncate
no2_30d.wt <- fifelse(no2_30d.wt<min.wt, min.wt, no2_30d.wt) # truncate
summary(no2_30d.wt)

## no2_3090d
model.num_no2_3090d = lm(no2_3090d ~ 1, data = birth) 

ps.num_no2_3090d <- dnorm((birth$no2_3090d - model.num_no2_3090d$fitted)/(summary(model.num_no2_3090d))$sigma,0,1)
# df <- approxfun(density(birth[, no2_3090d], kernel = "gaussian", n = floor(sqrt(dim(birth)[1]))))
# ps.num_no2_3090d <- df(birth[,no2_3090d])
ps.den_no2_3090d <- dnorm((birth$no2_3090d - GBM.fitted_no2_3090d)/sd(birth$no2_3090d - GBM.fitted_no2_3090d),0,1)

no2_3090d.wt <- ps.num_no2_3090d/ps.den_no2_3090d
summary(no2_3090d.wt)

max.wt <- quantile(no2_3090d.wt, 0.975)
min.wt <- quantile(no2_3090d.wt, 0.025)
no2_3090d.wt <- fifelse(no2_3090d.wt>max.wt, max.wt, no2_3090d.wt) # truncate
no2_3090d.wt <- fifelse(no2_3090d.wt<min.wt, min.wt, no2_3090d.wt) # truncate
summary(no2_3090d.wt)

## no2_90280d
model.num_no2_90280d = lm(no2_90280d ~ 1, data = birth) 

ps.num_no2_90280d <- dnorm((birth$no2_90280d - model.num_no2_90280d$fitted)/(summary(model.num_no2_90280d))$sigma,0,1)
# df <- approxfun(density(birth[, no2_90280d], kernel = "gaussian", n = floor(sqrt(dim(birth)[1]))))
# ps.num_no2_90280d <- df(birth[,no2_90280d])
ps.den_no2_90280d <- dnorm((birth$no2_90280d - GBM.fitted_no2_90280d)/sd(birth$no2_90280d - GBM.fitted_no2_90280d),0,1)

no2_90280d.wt <- ps.num_no2_90280d/ps.den_no2_90280d
summary(no2_90280d.wt)

max.wt <- quantile(no2_90280d.wt, 0.975)
min.wt <- quantile(no2_90280d.wt, 0.025)
no2_90280d.wt <- fifelse(no2_90280d.wt>max.wt, max.wt, no2_90280d.wt) # truncate
no2_90280d.wt <- fifelse(no2_90280d.wt<min.wt, min.wt, no2_90280d.wt) # truncate
summary(no2_90280d.wt)

######################### check ipw distribution ##############################
# pdf(paste0(dir_output_ipwplots,"hist_stdBC.pdf"))
# par(mfrow = c(2,2))
# hist((birth$bc_30d - model.num_bc_30d$fitted)/(summary(model.num_bc_30d))$sigma, 
#      main = "Histogram of standardized bc_30d",
#      xlab = NULL)
# hist((birth$bc_3090d - model.num_bc_3090d$fitted)/(summary(model.num_bc_3090d))$sigma, 
#      main = "Histogram of standardized bc_3090d",
#      xlab = NULL)
# hist((birth$bc_90280d - model.num_bc_90280d$fitted)/(summary(model.num_bc_90280d))$sigma, 
#      main = "Histogram of standardized bc_90280d",
#      xlab = NULL)
# dev.off()
# 
# pdf(paste0(dir_output_ipwplots,"hist_stdNO2.pdf"))
# par(mfrow = c(2,2))
# hist((birth$no2_30d - model.num_no2_30d$fitted)/(summary(model.num_no2_30d))$sigma, 
#      main = "Histogram of standardized no2_30d",
#      xlab = NULL)
# hist((birth$no2_3090d - model.num_no2_3090d$fitted)/(summary(model.num_no2_3090d))$sigma, 
#      main = "Histogram of standardized no2_3090d",
#      xlab = NULL)
# hist((birth$no2_90280d - model.num_no2_90280d$fitted)/(summary(model.num_no2_90280d))$sigma, 
#      main = "Histogram of standardized no2_90280d",
#      xlab = NULL)
# dev.off()

pdf(paste0(dir_output_ipwplots, file="hist_stdfitBCres.pdf"))
par(mfrow = c(2,2))
hist((birth$bc_30d - GBM.fitted_bc_30d)/sd(birth$bc_30d - GBM.fitted_bc_30d),
     main = "Histogram of standardized \n fitted bc_30d residuals",
     xlab = NULL)
hist((birth$bc_3090d - GBM.fitted_bc_3090d)/sd(birth$bc_3090d - GBM.fitted_bc_3090d),
     main = "Histogram of standardized \n fitted bc_3090d residuals",
     xlab = NULL)
hist((birth$bc_90280d - GBM.fitted_bc_90280d)/sd(birth$bc_90280d - GBM.fitted_bc_90280d),
     main = "Histogram of standardized \n fitted bc_90280d residuals",
     xlab = NULL)
dev.off()

pdf(paste0(dir_output_ipwplots, file="hist_stdfitNO2res.pdf"))
par(mfrow = c(2,2))
hist((birth$no2_30d - GBM.fitted_no2_30d)/sd(birth$no2_30d - GBM.fitted_no2_30d),
     main = "Histogram of standardized \n fitted no2_30d residuals",
     xlab = NULL)
hist((birth$no2_3090d - GBM.fitted_no2_3090d)/sd(birth$no2_3090d - GBM.fitted_no2_3090d),
     main = "Histogram of standardized \n fitted no2_3090d residuals",
     xlab = NULL)
hist((birth$no2_90280d - GBM.fitted_no2_90280d)/sd(birth$no2_90280d - GBM.fitted_no2_90280d),
     main = "Histogram of standardized \n fitted no2_90280d residuals",
     xlab = NULL)
dev.off()

pdf(paste0(dir_output_ipwplots, "hist_ipw_trct.pdf"))
par(mfrow = c(2,3))
hist(bc_30d.wt, main = "Histogram of truncated IPWs for bc_30d")
hist(bc_3090d.wt, main = "Histogram of truncated IPWs for bc_3090d")
hist(bc_90280d.wt, main = "Histogram of truncated IPWs for bc_90280d")
hist(no2_30d.wt, main = "Histogram of truncated IPWs for no2_30d")
hist(no2_3090d.wt, main = "Histogram of truncated IPWs for no2_3090d")
hist(no2_90280d.wt, main = "Histogram of truncated IPWs for no2_90280d")
dev.off()

############################# output IPW data #################################
birth_ipw <- data.table(birth[,uniqueid_yr],
                        bc_30d.wt, bc_3090d.wt, bc_90280d.wt,
                        ps.num_bc_30d, ps.den_bc_30d,
                        ps.num_bc_3090d, ps.den_bc_3090d,
                        ps.num_bc_90280d, ps.den_bc_90280d,
                        no2_30d.wt, no2_3090d.wt, no2_90280d.wt,
                        ps.num_no2_30d, ps.den_no2_30d,
                        ps.num_no2_3090d, ps.den_no2_3090d,
                        ps.num_no2_90280d, ps.den_no2_90280d)
colnames(birth_ipw)[1] <- "uniqueid_yr"
head(birth_ipw)
fwrite(birth_ipw, paste0(dir_birthdata, "MAbirth_ipw.csv"))

## bc and no2 all-period average
birth_ipw_2 <- data.table(birth[,uniqueid_yr],
                        bc_30d.wt, bc_3090d.wt,
                        bc_90280d.wt)
colnames(birth_ipw_2)[1] <- "uniqueid_yr"
head(birth_ipw_2)
dir_birthdata <- "/media/qnap3/Shuxin/airPollution_MAbirth/"
fwrite(birth_ipw_2, paste0(dir_birthdata, "MAbirth_temp_bc.csv"))
