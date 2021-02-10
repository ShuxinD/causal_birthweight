###############################################################################
# Project: Causal black carbon and NO2 on birth weight in MA                  #
# Code: (sensAnalysis) Estimate IPW with hyperparameter we set                #
# Input: "MAbirth_for_analyses.csv" birth data                                #
# Input: "BestCombinations_gbm_sens.csv" grid search results                  #                                        
# Output: weights dataset "MAbirth_ipw_sens.csv" with ID and weights          #
# Author: Shuxin Dong                                                         #
# Date: 2021-02-03                                                            #
###############################################################################

############################# 0. Setup ########################################
rm(list = ls())
gc()

library(data.table)
library(polycor)
library(h2o)
library(parallel)
n_cores <- detectCores() - 1 

# dir_input <- "/Users/shuxind/Desktop/BC_birthweight_data/"
setwd("/media/gate/Shuxin")
dir_birthdata <- "/media/gate/Shuxin/MAbirth/"
dir_gridsearch <- "/media/gate/Shuxin/MAbirth/results/SensitivityAnalysis/1GridSearchResults/"
dir_output_ipwplots <- "/media/gate/Shuxin/MAbirth/results/SensitivityAnalysis/2ipwPlots/"

## set default parameters for H2O
min.rows <- 10
learn.rate <- 0.005

## read grid search results
best <- fread(paste0(dir_gridsearch, "BestCombinations_gbm_sens.csv"))
names(best)
colnames(best)[3] <- "n.trees"
names(best)
# [1] "label"    "bestAAC"  "n.trees"  "depth"    "col_rate"

########################## 1. Load birth data #################################
## load data
birth <- fread(paste0(dir_birthdata, "MAbirth_for_analyses.csv"))
names(birth)
# > names(birth)
# [1] "uniqueid_yr"    "year"           "sex"            "married"        "mage"          
# [6] "mrace"          "m_edu"          "cigdpp"         "cigddp"         "clinega"       
# [11] "kotck"          "pncgov"         "bwg"            "rf_db_gest"     "rf_db_other"   
# [16] "rf_hbp_chronic" "rf_hbp_pregn"   "rf_cervix"      "rf_prev_4kg"    "rf_prev_sga"   
# [21] "mhincome"       "mhvalue"        "percentPoverty" "bc_30d"         "bc_3090d"      
# [26] "bc_90280d"      "no2_30d"        "no2_3090d"      "no2_90280d"     "lbw"           
# [31] "firstborn"      "m_wg_cat"       "smoker_ddp"     "smoker_dpp"     "mrace_1"       
# [36] "mrace_2"        "mrace_3"        "mrace_4"        "log_mhincome"   "log_mhvalue"   
summary(birth)
birth$year <- as.factor(birth$year)
birth$m_edu <- as.factor(birth$m_edu)
birth$kotck <- as.factor(birth$kotck)
birth$m_wg_cat <- as.factor(birth$m_wg_cat)
# birth <- sample_n(birth, floor(0.001*dim(birth)[1])) # sample as an example

var <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
         "clinega", "kotck","pncgov", "rf_db_gest","rf_db_other",
         "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
         "rf_prev_sga", "percentPoverty",
         "bc_30d","bc_3090d", "bc_90280d", 
         "no2_30d", "no2_3090d", "no2_90280d",
         "firstborn","m_wg_cat", "smoker_ddp", "smoker_dpp",
         "mrace_1", "mrace_2", "mrace_3", "mrace_4",
         "log_mhincome", "log_mhvalue")

############################# 2. Fit GBM ######################################

h2o.init(nthreads = n_cores, min_mem_size = "200G", port = 54345)

############################## 2.1 bc_30d ######################################
## construct data
dt <- birth[ , var, with = F]
dt[, T := bc_30d] # change
dt[, bc_30d := NULL] # change

independent <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
                 "clinega", "kotck","pncgov", "rf_db_gest",
                 "rf_db_other",
                 "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
                 "rf_prev_sga", "percentPoverty",
                 # "bc_30d",
                 "bc_3090d", "bc_90280d", 
                 "no2_30d", "no2_3090d", "no2_90280d",
                 "firstborn","m_wg_cat", "smoker_ddp", "smoker_dpp",
                 "mrace_1", "mrace_2", "mrace_3", "mrace_4",
                 "log_mhincome", "log_mhvalue") # change
independent <- independent[(independent!="rf_db_gest") & (independent!="rf_hbp_pregn") & (independent!="m_wg_cat")]

model_best <- best[label=="bc_30d",] #change

## fit GBM
birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
cat("fitting..bc_30d...")
gbm <- h2o.gbm(y = "T",
               x = independent,
               training_frame = birth.hex,
               ntrees = model_best$n.trees,
               max_depth = model_best$depth,
               col_sample_rate = model_best$col_rate,
               min_rows = min.rows,
               learn_rate = learn.rate, 
               distribution = "gaussian")
cat("predicting...bc_30d...\n")
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
                 "log_mhincome", "log_mhvalue") # change
independent <- independent[(independent!="rf_db_gest") & (independent!="rf_hbp_pregn") & (independent!="m_wg_cat")]

model_best <- best[label=="bc_3090d",] #change

## fit GBM
birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
cat("fitting..bc_3090d...\n")
gbm <- h2o.gbm(y = "T",
               x = independent,
               training_frame = birth.hex,
               ntrees = model_best$n.trees,
               max_depth = model_best$depth,
               col_sample_rate = model_best$col_rate,
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
                 "log_mhincome", "log_mhvalue") # change
independent <- independent[(independent!="rf_db_gest") & (independent!="rf_hbp_pregn") & (independent!="m_wg_cat")]

model_best <- best[label=="bc_90280d",] #change

## fit GBM
birth.hex <- as.h2o(dt, destination_frame = "birth.hex")
cat("fitting..bc_90280d...")
gbm <- h2o.gbm(y = "T",
               x = independent,
               training_frame = birth.hex,
               ntrees = model_best$n.trees,
               max_depth = model_best$depth,
               col_sample_rate = model_best$col_rate,
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
independent <- independent[(independent!="rf_db_gest") & (independent!="rf_hbp_pregn") & (independent!="m_wg_cat")]

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
independent <- independent[(independent!="rf_db_gest") & (independent!="rf_hbp_pregn") & (independent!="m_wg_cat")]

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
independent <- independent[(independent!="rf_db_gest") & (independent!="rf_hbp_pregn") & (independent!="m_wg_cat")]

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
## bc_30d
model.num_bc_30d = lm(bc_30d ~ 1, data = birth) 

ps.num_bc_30d <- dnorm((birth$bc_30d - model.num_bc_30d$fitted)/(summary(model.num_bc_30d))$sigma,0,1)
ps.den_bc_30d <- dnorm((birth$bc_30d - GBM.fitted_bc_30d)/sd(birth$bc_30d - GBM.fitted_bc_30d),0,1)

bc_30d.wt <- ps.num_bc_30d/ps.den_bc_30d
summary(bc_30d.wt)

max.wt <- quantile(bc_30d.wt, 0.99)
min.wt <- quantile(bc_30d.wt, 0.01)
bc_30d.wt <- fifelse(bc_30d.wt>max.wt, max.wt, bc_30d.wt) # truncate
bc_30d.wt <- fifelse(bc_30d.wt<min.wt, min.wt, bc_30d.wt) # truncate
summary(bc_30d.wt)

## bc_3090d
model.num_bc_3090d = lm(bc_3090d ~ 1, data = birth) 

ps.num_bc_3090d <- dnorm((birth$bc_3090d - model.num_bc_3090d$fitted)/(summary(model.num_bc_3090d))$sigma,0,1)
ps.den_bc_3090d <- dnorm((birth$bc_3090d - GBM.fitted_bc_3090d)/sd(birth$bc_3090d - GBM.fitted_bc_3090d),0,1)

bc_3090d.wt <- ps.num_bc_3090d/ps.den_bc_3090d
summary(bc_3090d.wt)

max.wt <- quantile(bc_3090d.wt, 0.99)
min.wt <- quantile(bc_3090d.wt, 0.01)
bc_3090d.wt <- fifelse(bc_3090d.wt>max.wt, max.wt, bc_3090d.wt) # truncate
bc_3090d.wt <- fifelse(bc_3090d.wt<min.wt, min.wt, bc_3090d.wt) # truncate
summary(bc_3090d.wt)

## bc_90280d
model.num_bc_90280d = lm(bc_90280d ~ 1, data = birth) 

ps.num_bc_90280d <- dnorm((birth$bc_90280d - model.num_bc_90280d$fitted)/(summary(model.num_bc_90280d))$sigma,0,1)
ps.den_bc_90280d <- dnorm((birth$bc_90280d - GBM.fitted_bc_90280d)/sd(birth$bc_90280d - GBM.fitted_bc_90280d),0,1)

bc_90280d.wt <- ps.num_bc_90280d/ps.den_bc_90280d
summary(bc_90280d.wt)

max.wt <- quantile(bc_90280d.wt, 0.99)
min.wt <- quantile(bc_90280d.wt, 0.01)
bc_90280d.wt <- fifelse(bc_90280d.wt>max.wt, max.wt, bc_90280d.wt) # truncate
bc_90280d.wt <- fifelse(bc_90280d.wt<min.wt, min.wt, bc_90280d.wt) # truncate
summary(bc_90280d.wt)

## no2_30d
model.num_no2_30d = lm(no2_30d ~ 1, data = birth) 

ps.num_no2_30d <- dnorm((birth$no2_30d - model.num_no2_30d$fitted)/(summary(model.num_no2_30d))$sigma,0,1)
ps.den_no2_30d <- dnorm((birth$no2_30d - GBM.fitted_no2_30d)/sd(birth$no2_30d - GBM.fitted_no2_30d),0,1)

no2_30d.wt <- ps.num_no2_30d/ps.den_no2_30d
summary(no2_30d.wt)

max.wt <- quantile(no2_30d.wt, 0.99)
min.wt <- quantile(no2_30d.wt, 0.01)
no2_30d.wt <- fifelse(no2_30d.wt>max.wt, max.wt, no2_30d.wt) # truncate
no2_30d.wt <- fifelse(no2_30d.wt<min.wt, min.wt, no2_30d.wt) # truncate
summary(no2_30d.wt)

## no2_3090d
model.num_no2_3090d = lm(no2_3090d ~ 1, data = birth) 

ps.num_no2_3090d <- dnorm((birth$no2_3090d - model.num_no2_3090d$fitted)/(summary(model.num_no2_3090d))$sigma,0,1)
ps.den_no2_3090d <- dnorm((birth$no2_3090d - GBM.fitted_no2_3090d)/sd(birth$no2_3090d - GBM.fitted_no2_3090d),0,1)

no2_3090d.wt <- ps.num_no2_3090d/ps.den_no2_3090d
summary(no2_3090d.wt)

max.wt <- quantile(no2_3090d.wt, 0.99)
min.wt <- quantile(no2_3090d.wt, 0.01)
no2_3090d.wt <- fifelse(no2_3090d.wt>max.wt, max.wt, no2_3090d.wt) # truncate
no2_3090d.wt <- fifelse(no2_3090d.wt<min.wt, min.wt, no2_3090d.wt) # truncate
summary(no2_3090d.wt)

## no2_90280d
model.num_no2_90280d = lm(no2_90280d ~ 1, data = birth) 

ps.num_no2_90280d <- dnorm((birth$no2_90280d - model.num_no2_90280d$fitted)/(summary(model.num_no2_90280d))$sigma,0,1)
ps.den_no2_90280d <- dnorm((birth$no2_90280d - GBM.fitted_no2_90280d)/sd(birth$no2_90280d - GBM.fitted_no2_90280d),0,1)

no2_90280d.wt <- ps.num_no2_90280d/ps.den_no2_90280d
summary(no2_90280d.wt)

max.wt <- quantile(no2_90280d.wt, 0.99)
min.wt <- quantile(no2_90280d.wt, 0.01)
no2_90280d.wt <- fifelse(no2_90280d.wt>max.wt, max.wt, no2_90280d.wt) # truncate
no2_90280d.wt <- fifelse(no2_90280d.wt<min.wt, min.wt, no2_90280d.wt) # truncate
summary(no2_90280d.wt)

######################### check ipw distribution ##############################
pdf(paste0(dir_output_ipwplots,"hist_stdBC.pdf"))
par(mfrow = c(2,2))
hist((birth$bc_30d - model.num_bc_30d$fitted)/(summary(model.num_bc_30d))$sigma, 
     main = "Histogram of standardized bc_30d",
     xlab = NULL)
hist((birth$bc_3090d - model.num_bc_3090d$fitted)/(summary(model.num_bc_3090d))$sigma, 
     main = "Histogram of standardized bc_3090d",
     xlab = NULL)
hist((birth$bc_90280d - model.num_bc_90280d$fitted)/(summary(model.num_bc_90280d))$sigma, 
     main = "Histogram of standardized bc_90280d",
     xlab = NULL)
dev.off()

pdf(paste0(dir_output_ipwplots,"hist_stdNO2.pdf"))
par(mfrow = c(2,2))
hist((birth$no2_30d - model.num_no2_30d$fitted)/(summary(model.num_no2_30d))$sigma, 
     main = "Histogram of standardized no2_30d",
     xlab = NULL)
hist((birth$no2_3090d - model.num_no2_3090d$fitted)/(summary(model.num_no2_3090d))$sigma, 
     main = "Histogram of standardized no2_3090d",
     xlab = NULL)
hist((birth$no2_90280d - model.num_no2_90280d$fitted)/(summary(model.num_no2_90280d))$sigma, 
     main = "Histogram of standardized no2_90280d",
     xlab = NULL)
dev.off()

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
fwrite(birth_ipw, paste0(dir_birthdata, "MAbirth_ipw_sens.csv"))

