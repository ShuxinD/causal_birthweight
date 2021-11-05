#' Project: causal_birthweight
#' Code: estimate stabilized ipw with grid search results and check balance
#' Input: "MAbirth_for_analyses.csv" birth data
#' Input: grid search results
#' Output: weights dataset "MAbirth_all_ipw.csv" with ID and weights
#' Author: Shuxin Dong
#' First Ceate Date: 2021-01-21

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
dir_data <- "/media/qnap3/Shuxin/airPollution_MAbirth/"
dir_gridsearch <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/1GridSearchResults/"
dir_ipwraw <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/2ipw/"

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
birth <- fread(paste0(dir_data, "MAbirth_for_analyses.csv"))
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

## load grid search results ----
gs_bc_all <- fread(paste0(dir_gridsearch, "GridSearch_bc_all1104.csv"))
gs_no2_all <- fread(paste0(dir_gridsearch, "GridSearch_no2_all1104.csv"))

gs <- rbind(gs_bc_all[objective==min(objective),][,.(objective, rate, depth, ntree)][,exposure:="bc_all"],
      gs_no2_all[objective==min(objective),][,.(objective, rate, depth, ntree)][,exposure:="no2_all"])
print(gs)

## model GPS generate IPW ----
#' set default parameters for H2O
min.rows <- 10 
learn.rate <- 0.01

h2o.init(nthreads = n_cores, min_mem_size = "250G", port = 54345)
for (ps_exposures_i in ps_exposures) {
  cat("FITTING ", ps_exposures_i, "\n")
  response <- ps_exposures_i
  predictor <- c(ps_vars, ps_exposures[ps_exposures!=response])
  birth_h2o <- h2o.importFile(paste0(dir_data, "MAbirth_for_analyses.csv"))
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
  fwrite(ipw_id, file = paste0(dir_ipwraw, ps_exposures_i, "_raw_gbm.csv"))
}
h2o.shutdown(prompt = FALSE)

## simple linear regression for GPS ----
# response <- "no2_all" # change everytime 
# ps_vars <- c("year","sex","married","mage", # "cigdpp","cigddp",
#              "clinega","pncgov", 
#              "rf_db_gest","rf_db_other", "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg", "rf_prev_sga", 
#              "smoker_ddp", # "smoker_dpp",
#              "mrace_1", "mrace_2", "mrace_3", # "mrace_4",
#              "m_edu_1", "m_edu_2", "m_edu_3","m_edu_4", # "m_edu_5",
#              "kotck_1","kotck_2","kotck_3","kotck_4",
#              "m_wg_cat_1","m_wg_cat_2","m_wg_cat_3","m_wg_cat_4",# "m_wg_cat_5",
#              # "log_mhvalue", "log_mhincome",
#              "percentPoverty", "firstborn")
for (ps_exposures_i in ps_exposures) {
  glm_exposure <- glm(get(ps_exposures_i) ~ year + sex + married + mage + 
                        clinega + pncgov +
                        rf_db_gest + rf_db_other + rf_hbp_chronic + rf_hbp_pregn + rf_cervix + rf_prev_4kg + rf_prev_sga + 
                        smoker_ddp + 
                        mrace_1 + mrace_2 + mrace_3 + mrace_4 +
                        m_edu_1 + m_edu_2 + m_edu_3 + 
                        kotck_1 + kotck_2 + kotck_3 + kotck_4 + 
                        m_wg_cat_1 + m_wg_cat_2 + m_wg_cat_3 + m_wg_cat_4 + 
                        percentPoverty + firstborn,
                      data = birth)
  pred.glm <- glm_exposure$fitted
  ipw_raw <- generate_ipw(model_fitted_value = as.vector(pred.glm), birth_raw_data = birth, exposure_of_interest = ps_exposures_i)
  summary(ipw_raw)
  ipw_id <- data.table(uniqueid_yr = birth[,uniqueid_yr], ipw_raw = ipw_raw)
  fwrite(ipw_id, file = paste0(dir_ipwraw, ps_exposures_i, "_raw_glm.csv"))
}

## truncate IPW based on balance results ----
#' change for each exposure
exposure_interest <- "bc_all"
upper_percentile <- 0.9995
lower_percentile <- 0.0005

response <- exposure_interest
predictor <- c(ps_vars, ps_exposures[ps_exposures!=response])

ipw_id_gbm <- fread(paste0(dir_ipwraw, exposure_interest, ".csv"), colClasses = c("ipw_raw"="numeric"))
summary(as.numeric(ipw_id_gbm[,ipw_raw]))
ipw_gbm <- truncate_ipw(as.numeric(ipw_id_gbm[,ipw_raw]), upper_percentile, lower_percentile)
summary(ipw_gbm)

ipw_id_glm <- fread(paste0(dir_ipwraw, exposure_interest, "_glm.csv"), colClasses = c("ipw_raw"="numeric"))
summary(as.numeric(ipw_id_glm[,ipw_raw]))
ipw_glm <- truncate_ipw(as.numeric(ipw_id_glm[,ipw_raw]), upper_percentile, lower_percentile)
summary(ipw_glm)

ac <- matrix(NA,length(predictor),3)
for (j in 1:length(predictor)){ 
  ac[j,1] <-  ifelse(!is.factor(birth[,..predictor][[j]]), 
                     cor(birth[,get(response)], birth[,..predictor][[j]], method = "pearson"),
                     polyserial(birth[,get(response)], birth[,..predictor][[j]]))
  ac[j,2] <-  ifelse(!is.factor(birth[,..predictor][[j]]), 
                     weightedCorr(birth[,get(response)], birth[,..predictor][[j]], method = "pearson", weights = ipw_gbm),
                     weightedCorr(birth[,get(response)], birth[,..predictor][[j]], method = "polyserial", weights = ipw_gbm))
  ac[j,3] <-  ifelse(!is.factor(birth[,..predictor][[j]]), 
                     weightedCorr(birth[,get(response)], birth[,..predictor][[j]], method = "pearson", weights = ipw_glm),
                     weightedCorr(birth[,get(response)], birth[,..predictor][[j]], method = "polyserial", weights = ipw_glm))
  cat("finish column", j,"/", length(predictor), "\n")
}

balance <- data.frame(name.x=predictor,
                      corr=c(ac[,1],ac[,2],ac[,3]),
                      weighted = c(rep("unweighted", length(predictor)),rep("weighted_gbm", length(predictor)), rep("weighted_glm", length(predictor))))

ggplot(balance, aes(x = name.x, y = corr, group = weighted)) +
  # geom_line(aes(colour = weighted), size = 0.5, linetype = "dashed")+
  geom_point(aes(colour = weighted), size = 2) +
  geom_line(aes(colour = weighted), size = 0.5) +
  geom_hline(yintercept=0, size=0.2) +
  geom_hline(yintercept=-0.1, size=0.2, linetype = "dashed") +
  geom_hline(yintercept=0.1, size=0.2, linetype = "dashed") +
  # geom_hline(yintercept=-0.2, size=0.1, linetype = "dashed") +
  # geom_hline(yintercept=0.2, size=0.1, linetype = "dashed") +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  coord_flip() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        text = element_text(size=16))

## save the ipw
assign(paste0("ipw_gbm_", exposure_interest), ipw_gbm)
assign(paste0("ipw_glm_", exposure_interest), ipw_glm)

IPWs <- data.table(uniqueid_yr = birth[,uniqueid_yr],
                   ipw_bc_all,
                   ipw_no2_all)
fwrite(IPWs, file = "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/2ipw/IPWs_all_1104.csv")
