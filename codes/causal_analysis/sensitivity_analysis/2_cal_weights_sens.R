###############################################################################
# Project: Causal black carbon on birth weight in MA                          #
# Code: Estimate GPS with hyperparameter we set                               #
# Input: birth_final.csv, hyperparameters                                     #
# Output: weights dataset "birth_ipw.csv" with ID and weights                 #
# Author: Shuxin Dong                                                         #
# Date: Nov 16, 2020                                                          #
###############################################################################

############################# 0. Setup ########################################
rm(list = ls())
gc()

library(dplyr)
library(data.table)
library(polycor)
library(h2o)
library(parallel)
n_cores <- detectCores() - 1 

# dir_input <- "/Users/shuxind/Desktop/BC_birthweight_data/"
setwd("/media/qnap3/Shuxin")
dir_input <- "/media/qnap3/Shuxin/"
dir_output <- "/media/qnap3/Shuxin/"
## set default parameters for H2O
min.rows <- 10
learn.rate <- 0.005
## set the hyperparameter values after grid search
n.trees_bc30d <- 19691
max.depth_bc30d <- 6
col.sample.rate_bc30d <- 1.0
n.trees_bc3090d <- 15451
max.depth_bc3090d <- 6
col.sample.rate_bc3090d <- 1.0
n.trees_bc90280d <- 9675
max.depth_bc90280d <- 6
col.sample.rate_bc90280d <- 1.0
############################# 1. Fit GBM ######################################
############################# 1.1 bc_30d ######################################
## load data
birth <- fread(paste0(dir_input, "birth_final.csv"),
               drop = "V1")
birth$year <- as.factor(birth$year)
birth$m_edu <- as.factor(birth$m_edu)
birth$kotck <- as.factor(birth$kotck)
birth$m_wg_cat <- as.factor(birth$m_wg_cat)

var <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
         "clinega","kotck","pncgov", 
         # "rf_db_gest",
         "rf_db_other",
         "rf_hbp_chronic", 
         # "rf_hbp_pregn",
         "rf_cervix","rf_prev_4kg",
         "rf_prev_sga", "bc_30d","bc_3090d", "bc_90280d", "firstborn",
         # "m_wg_cat",
         "log_mhincome", "log_mhvalue", "percentPoverty",
         "mrace_1", "mrace_2", "mrace_3", "mrace_4")
birth <- birth[ , var, with = F]
birth[, T := bc_30d]
birth[, bc_30d := NULL]

h2o.init(nthreads = n_cores, min_mem_size = "200G", port = 54345)
independent <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
                 "clinega","kotck","pncgov", 
                 # "rf_db_gest",
                 "rf_db_other",
                 "rf_hbp_chronic", 
                 # "rf_hbp_pregn",
                 "rf_cervix","rf_prev_4kg",
                 "rf_prev_sga","firstborn",
                 # "m_wg_cat",
                 "log_mhincome", "log_mhvalue", "percentPoverty",
                 "mrace_1", "mrace_2", "mrace_3", "mrace_4",
                 "bc_3090d", "bc_90280d")
birth.hex <- as.h2o(birth, destination_frame = "birth.hex")
gbm_30d <- h2o.gbm(y = "T",
                   x = independent,
                   training_frame = birth.hex,
                   ntrees = n.trees_bc30d, 
                   max_depth = max.depth_bc30d,
                   min_rows = min.rows,
                   learn_rate = learn.rate, 
                   col_sample_rate = col.sample.rate_bc30d,
                   distribution = "gaussian")
pred.gbm_30d <- h2o.predict(object = gbm_30d, newdata = birth.hex)
GBM.fitted_30d <- as.vector(pred.gbm_30d)
h2o.shutdown(prompt = FALSE)

############################# 1.2 bc_3090d ####################################
birth <- fread(paste0(dir_input, "birth_final.csv"),
               drop = "V1")
birth$year <- as.factor(birth$year)
birth$m_edu <- as.factor(birth$m_edu)
birth$kotck <- as.factor(birth$kotck)
birth$m_wg_cat <- as.factor(birth$m_wg_cat)

var <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
         "clinega","kotck","pncgov", 
         # "rf_db_gest",
         "rf_db_other",
         "rf_hbp_chronic", 
         # "rf_hbp_pregn",
         "rf_cervix","rf_prev_4kg",
         "rf_prev_sga", "bc_30d","bc_3090d", "bc_90280d", "firstborn",
         # "m_wg_cat",
         "log_mhincome", "log_mhvalue", "percentPoverty",
         "mrace_1", "mrace_2", "mrace_3", "mrace_4")
birth <- birth[ , var, with = F]
birth[, T := bc_3090d]
birth[, bc_3090d := NULL]

h2o.init(nthreads = n_cores, min_mem_size = "200G", port = 54345)
independent <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
                 "clinega","kotck","pncgov", 
                 # "rf_db_gest",
                 "rf_db_other",
                 "rf_hbp_chronic", 
                 # "rf_hbp_pregn",
                 "rf_cervix","rf_prev_4kg",
                 "rf_prev_sga","firstborn",
                 # "m_wg_cat",
                 "log_mhincome", "log_mhvalue", "percentPoverty",
                 "mrace_1", "mrace_2", "mrace_3", "mrace_4",
                 "bc_30d", "bc_90280d")
birth.hex <- as.h2o(birth, destination_frame = "birth.hex")
gbm_3090d <- h2o.gbm(y = "T",
                     x = independent,
                     training_frame = birth.hex,
                     ntrees = n.trees_bc3090d, 
                     max_depth = max.depth_bc3090d, # change
                     min_rows = min.rows,
                     learn_rate = learn.rate, 
                     col_sample_rate = col.sample.rate_bc3090d,
                     distribution = "gaussian")
pred.gbm_3090d <- h2o.predict(object = gbm_3090d, newdata = birth.hex)
GBM.fitted_3090d <- as.vector(pred.gbm_3090d)
h2o.shutdown(prompt = FALSE)

############################# 1.3 bc_90280d ###################################
birth <- fread(paste0(dir_input, "birth_final.csv"),
               drop = "V1")
birth$year <- as.factor(birth$year)
birth$m_edu <- as.factor(birth$m_edu)
birth$kotck <- as.factor(birth$kotck)
birth$m_wg_cat <- as.factor(birth$m_wg_cat)

var <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
         "clinega","kotck","pncgov", 
         # "rf_db_gest",
         "rf_db_other",
         "rf_hbp_chronic", 
         # "rf_hbp_pregn",
         "rf_cervix","rf_prev_4kg",
         "rf_prev_sga", "bc_30d","bc_3090d", "bc_90280d", "firstborn",
         # "m_wg_cat",
         "log_mhincome", "log_mhvalue", "percentPoverty",
         "mrace_1", "mrace_2", "mrace_3", "mrace_4")
birth <- birth[ , var, with = F]
birth[, T := bc_90280d]
birth[, bc_90280d := NULL]

h2o.init(nthreads = n_cores, min_mem_size = "200G", port = 54345)
independent <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
                 "clinega","kotck","pncgov", 
                 # "rf_db_gest",
                 "rf_db_other",
                 "rf_hbp_chronic", 
                 # "rf_hbp_pregn",
                 "rf_cervix","rf_prev_4kg",
                 "rf_prev_sga","firstborn",
                 # "m_wg_cat",
                 "log_mhincome", "log_mhvalue", "percentPoverty",
                 "mrace_1", "mrace_2", "mrace_3", "mrace_4",
                 "bc_30d", "bc_3090d")
birth.hex <- as.h2o(birth, destination_frame = "birth.hex")
gbm_90280d <- h2o.gbm(y = "T",
                     x = independent,
                     training_frame = birth.hex,
                     ntrees = n.trees_bc90280d, 
                     max_depth = max.depth_bc90280d, # change
                     min_rows = min.rows,
                     learn_rate = learn.rate, 
                     col_sample_rate = col.sample.rate_bc90280d,
                     distribution = "gaussian")
pred.gbm_90280d <- h2o.predict(object = gbm_90280d, newdata = birth.hex)
GBM.fitted_90280d <- as.vector(pred.gbm_90280d)
h2o.shutdown(prompt = FALSE)
gc()
############################# 3. Calculate Weights ############################
## load data
birth <- fread(paste0(dir_input, "birth_final.csv"),
               drop = "V1")
birth$year <- as.factor(birth$year)
birth$m_edu <- as.factor(birth$m_edu)
birth$kotck <- as.factor(birth$kotck)
birth$m_wg_cat <- as.factor(birth$m_wg_cat)

############################# calculate weights ###############################
## bc_30d
model.num_30d = lm(bc_30d ~ 1, data = birth) 
ps.num_30d <- dnorm((birth$bc_30d - model.num_30d$fitted)/(summary(model.num_30d))$sigma,0,1)
ps.den_30d <- dnorm((birth$bc_30d - GBM.fitted_30d)/sd(birth$bc_30d - GBM.fitted_30d),0,1)
bc_30d.wt <- ps.num_30d/ps.den_30d
bc_30d.wt.t <- bc_30d.wt
bc_30d.wt.t <- fifelse(bc_30d.wt.t>quantile(bc_30d.wt, 0.99), 
                       quantile(bc_30d.wt, 0.99), bc_30d.wt.t)
bc_30d.wt.t <- fifelse(bc_30d.wt.t<quantile(bc_30d.wt, 0.01), 
                       quantile(bc_30d.wt, 0.01), bc_30d.wt.t)
summary(bc_30d.wt)
summary(bc_30d.wt.t)

## bc_3090d
model.num_3090d = lm(bc_3090d ~ 1, data = birth) 
ps.num_3090d <- dnorm((birth$bc_3090d - model.num_3090d$fitted)/(summary(model.num_3090d))$sigma,0,1)
ps.den_3090d <- dnorm((birth$bc_3090d - GBM.fitted_3090d)/sd(birth$bc_3090d - GBM.fitted_3090d),0,1)
bc_3090d.wt <- ps.num_3090d/ps.den_3090d
bc_3090d.wt.t <- bc_3090d.wt
bc_3090d.wt.t <- fifelse(bc_3090d.wt.t>quantile(bc_3090d.wt, 0.99), 
                       quantile(bc_3090d.wt, 0.99), bc_3090d.wt.t)
bc_3090d.wt.t <- fifelse(bc_3090d.wt.t<quantile(bc_3090d.wt, 0.01), 
                       quantile(bc_3090d.wt, 0.01), bc_3090d.wt.t)
summary(bc_3090d.wt)
summary(bc_3090d.wt.t)

## bc_90280d
model.num_90280d = lm(bc_90280d ~ 1, data = birth) 
ps.num_90280d <- dnorm((birth$bc_90280d - model.num_90280d$fitted)/(summary(model.num_90280d))$sigma,0,1)
ps.den_90280d <- dnorm((birth$bc_90280d - GBM.fitted_90280d)/sd(birth$bc_90280d - GBM.fitted_90280d),0,1)
bc_90280d.wt <- ps.num_90280d/ps.den_90280d
bc_90280d.wt.t <- bc_90280d.wt
bc_90280d.wt.t <- fifelse(bc_90280d.wt.t>quantile(bc_90280d.wt, 0.99), 
                       quantile(bc_90280d.wt, 0.99), bc_90280d.wt.t)
bc_90280d.wt.t <- fifelse(bc_90280d.wt.t<quantile(bc_90280d.wt, 0.01), 
                       quantile(bc_90280d.wt, 0.01), bc_90280d.wt.t)
summary(bc_90280d.wt)
summary(bc_90280d.wt.t)

############################# check distribution ##############################
pdf(file="hist_stdBCsens.pdf")
par(mfrow = c(2,2))
hist((birth$bc_30d - model.num_30d$fitted)/(summary(model.num_30d))$sigma, 
     main = "Histogram of standardized bc_30d",
     xlab = NULL)
hist((birth$bc_3090d - model.num_3090d$fitted)/(summary(model.num_3090d))$sigma, 
     main = "Histogram of standardized bc_3090d",
     xlab = NULL)
hist((birth$bc_90280d - model.num_90280d$fitted)/(summary(model.num_90280d))$sigma, 
     main = "Histogram of standardized bc_90280d",
     xlab = NULL)
dev.off()

pdf(file="hist_stdfitBCressens.pdf")
par(mfrow = c(2,2))
hist((birth$bc_30d - GBM.fitted_30d)/sd(birth$bc_30d - GBM.fitted_30d),
     main = "Histogram of standardized \n fitted bc_30d residuals",
     xlab = NULL)
hist((birth$bc_3090d - GBM.fitted_3090d)/sd(birth$bc_3090d - GBM.fitted_3090d),
     main = "Histogram of standardized \n fitted bc_3090d residuals",
     xlab = NULL)
hist((birth$bc_90280d - GBM.fitted_90280d)/sd(birth$bc_90280d - GBM.fitted_90280d),
     main = "Histogram of standardized \n fitted bc_90280d residuals",
     xlab = NULL)
dev.off()

pdf(file="boxplot_IPW_noTrctsens.pdf")
par(mfrow = c(2,2))
boxplot(bc_30d.wt, xlab = "IPWs for bc_30d w/o truncation")
boxplot(bc_3090d.wt, xlab = "IPWs for bc_3090d w/o truncation")
boxplot(bc_90280d.wt, xlab = "IPWs for bc_90280d w/o truncation")
dev.off()

pdf(file="hist_IPW_trctsens.pdf")
par(mfrow = c(2,2))
hist(bc_30d.wt.t, main = "Histogram of truncated IPWs for bc_30d")
hist(bc_3090d.wt.t, main = "Histogram of truncated IPWs for bc_3090d")
hist(bc_90280d.wt.t, main = "Histogram of truncated IPWs for bc_90280d")
dev.off()

############################# output IPW data #################################
birth_ipw <- data.table(birth$uniqueid_yr,
                        bc_30d.wt.t, bc_3090d.wt.t, bc_90280d.wt.t,
                        bc_30d.wt, bc_3090d.wt, bc_90280d.wt.t,
                        ps.num_30d, ps.den_30d,
                        ps.num_3090d, ps.den_3090d,
                        ps.num_90280d, ps.den_90280d)
head(birth_ipw)
fwrite(birth_ipw, file = paste0(dir_output, "birth_ipwsens.csv"))
