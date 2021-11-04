#' Project: causal_birthweight
#' Code: estimate stabilized ipw with grid search results and check balance
#' Input: "MAbirth_for_analyses.csv" birth data
#' Output: weights dataset "MAbirth_ipw.csv" with ID and weights
#' Author: Shuxin Dong
#' First create date: 2021-01-27

## setup ----
rm(list = ls())
gc()

library(data.table)
library(polycor)
library(h2o)
library(parallel)
library(wCorr)
library(ggplot2)
n_cores <- detectCores() - 1 

setwd("/media/gate/Shuxin/")
dir_in <- "/media/qnap3/Shuxin/airPollution_MAbirth/"
dir_in_gs <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/1GridSearchResults/"
dir_out_ipwraw <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/2ipw/"

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

generate_ipw <- function(model_fitted_value, birth_raw_data, exposure_of_interest){
  numerator <- dnorm((birth_raw_data[,get(exposure_of_interest)]-mean(birth_raw_data[,get(exposure_of_interest)]))/sd(birth_raw_data[,get(exposure_of_interest)]),0,1)
  gps <- dnorm((birth_raw_data[,get(exposure_of_interest)]-model_fitted_value)/sd(birth_raw_data[,get(exposure_of_interest)]-model_fitted_value),0,1)
  ipw <- numerator/gps
  return(ipw)
}

## load birth data ----
birth <- fread(paste0(dir_in, "MAbirth_for_analyses.csv"))
# names(birth)
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

## load grid search results ----
gs <- fread(paste0(dir_in_gs, "GridSearchResults1028.csv"))
gs

## model GPS generate IPW ----
#' set default parameters for H2O
min.rows <- 10 
learn.rate <- 0.01

# n.trees <- 800
# max.depth <- 8
# col.sample.rate <- 1

h2o.init(nthreads = n_cores, min_mem_size = "250G", port = 54345)
for (ps_exposures_i in ps_exposures) {
  cat("FITTING ", ps_exposures_i, "\n")
  response <- ps_exposures_i
  predictor <- c(ps_vars, ps_exposures[ps_exposures!=response])
  birth_h2o <- h2o.importFile(paste0(dir_in, "MAbirth_for_analyses.csv"))
  birth_h2o$year <- as.factor(birth_h2o$year)
  gbm <- h2o.gbm(y = response, # exposure of interest
                 x = predictor,
                 training_frame = birth_h2o,
                 ntrees = gs[exposure==response, ntree], 
                 max_depth = gs[exposure==response, depth] ,
                 min_rows = min.rows,
                 learn_rate = learn.rate, 
                 col_sample_rate = gs[exposure==response, rate],
                 distribution = "gaussian")
  pred.gbm <- h2o.predict(object = gbm, newdata = birth_h2o)
  ipw_raw <- generate_ipw(model_fitted_value = as.vector(pred.gbm), birth_raw_data = birth, exposure_of_interest = response)
  summary(ipw_raw)
  ipw_id <- data.table(uniqueid_yr = birth[,uniqueid_yr], ipw_raw = ipw_raw)
  fwrite(ipw_id, file = paste0(dir_out_ipwraw, ps_exposures_i, ".csv"))
}
h2o.shutdown(prompt = FALSE)

## truncate IPW based on balance results ----
#' change for each exposure
exposure_interest <- "no2_90280d"
upper_percentile <- 0.9999
lower_percentile <- 0.0001

#' bc_30d 0.9995 0.0005
#' bc_3090d 0.9998 0.0002
#' bc_90280d 0.9993 0.0007
#' no2_30d 0.9993 0.0007
#' no2_3090d 0.9967 0.0033
#' no2_90280d 0.9999 0.0001

response <- exposure_interest
predictor <- c(ps_vars, ps_exposures[ps_exposures!=response])
ipw_id <- fread(paste0(dir_out_ipwraw, exposure_interest, ".csv"), colClasses = c("ipw_raw"="numeric"))
summary(as.numeric(ipw_id[,ipw_raw]))
ipw <- truncate_ipw(as.numeric(ipw_id[,ipw_raw]), upper_percentile, lower_percentile)
summary(ipw)
ac <- matrix(NA,length(predictor),2)
for (j in 1:length(predictor)){ ## unweighted
  ac[j,1] <-  ifelse(!is.factor(birth[,..predictor][[j]]), 
                     cor(birth[,get(response)], birth[,..predictor][[j]], method = "pearson"),
                     polyserial(birth[,get(response)], birth[,..predictor][[j]]))
  ac[j,2] <-  ifelse(!is.factor(birth[,..predictor][[j]]), 
                     weightedCorr(birth[,get(response)], birth[,..predictor][[j]], method = "pearson", weights = ipw),
                     weightedCorr(birth[,get(response)], birth[,..predictor][[j]], method = "polyserial", weights = ipw))
  cat("finish column", j,"/", length(predictor), "\n")
}
balance <- data.frame(name.x=predictor,
                      corr=c(ac[,1],ac[,2]),
                      weighted = c(rep("unweighted", length(predictor)),rep("weighted", length(predictor))))
ggplot(balance, aes(x = name.x, y = corr, group = weighted)) +
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

assign(paste0("ipw_", exposure_interest), ipw)

IPWs <- data.table(uniqueid_yr = birth[,uniqueid_yr],
                   ipw_bc_30d,
                   ipw_bc_3090d,
                   ipw_bc_90280d,
                   ipw_no2_30d,
                   ipw_no2_3090d,
                   ipw_no2_90280d)
fwrite(IPWs, file = "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/2ipw/IPWs1028.csv")


