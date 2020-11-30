###############################################################################
# Project: Causal black carbon on birth weight in MA                          #
# Code: GBM to estimate GPS with balance stopping rule                        #
# Code : from three non-overlapping periods                                   #
# Input: clean birth data, BC averages                                        #
# Output: best hyperpamater combination for three BC periods                  #
# Author: Shuxin Dong                                                         #
# Date: Oct 13, 2020                                                          #
###############################################################################

############################# 0. Setup ########################################
rm(list = ls())
gc()

library(dplyr)
library(h2o)
library(data.table)
library(polycor)
library(funique)
library(parallel)
n_cores <- detectCores() - 1 

# dir_input <- "/Users/shuxind/Desktop/BC_birthweight_data/"
setwd("/media/gate/Shuxin")
dir_input <- "/media/gate/Shuxin/"

## set parameters for h2o.gbm model
min.rows <- 10
learn.rate <- 0.005
rep.bootstrap <- 50

## set hyperparameters range for h2o.gbm model
n.trees <- 25000
max.depth <- c(6, 8, 10, 12)
col.sample.rate <- c(0.8, 0.9, 1.0)

########################### functions to compute AAC ##########################
## the following codes incorporate Zhu's code with large dataset
## change gbm package to h2o.gbm packge
## use data.table to manipulate the data instead of default
## use subsample of bootstrapping instead of whole sample 
F.aac.iter <- function(i, data, data.hex, ps.num, ps.model.pred, rep) {
  # i: number of iterations (number of trees) 
  # data: dataset containing the treatment and the covariates not in h2o structure.
  # data.hex: dataset containing the treatment and the covariates in h2o env.
  # ps.model.pred: the staged prediction results of boosting model to estimate (p(T_iX_i)) 
  # ps.num: the estimated p(T_i) 
  # rep: number of replications in bootstrap 
  GBM.fitted <- as.vector(ps.model.pred[,floor(i)])
  ps.den <- dnorm((data$T - GBM.fitted)/sd(data$T - GBM.fitted),0,1)
  wtnt <- ps.num/ps.den
  wt <- wtnt
  wt <- fifelse(wt>quantile(wtnt, 0.99), quantile(wtnt, 0.99), wt)
  wt <- fifelse(wt<quantile(wtnt, 0.01), quantile(wtnt, 0.01), wt)
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

############################# 2. data manipulation ############################
## load data
# birth <- fread(paste0(dir_input, "birth_final.csv"),
#                drop = "V1")
# birth$year <- as.factor(birth$year)
# birth$m_edu <- as.factor(birth$m_edu)
# birth$kotck <- as.factor(birth$kotck)
# birth$m_wg_cat <- as.factor(birth$m_wg_cat)
# var <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
#          "clinega","kotck","pncgov", "bwg", "rf_db_gest","rf_db_other",
#          "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
#          "rf_prev_sga", "bc_30d","bc_3090d", "bc_90280d", "firstborn","m_wg_cat",
#          "log_mhincome", "log_mhvalue", "percentPoverty",
#          "mrace_1", "mrace_2", "mrace_3", "mrace_4")
# birth <- birth[ , var, with = F]
# birth <- sample_n(birth, floor(0.001*dim(birth)[1])) # sample as an example

############################## 3. do h2o.gbm ###################################

############################## 3.1 bc_30d ######################################
## load data
birth <- fread(paste0(dir_input, "birth_final.csv"),
               drop = "V1")
birth$year <- as.factor(birth$year)
birth$m_edu <- as.factor(birth$m_edu)
birth$kotck <- as.factor(birth$kotck)
birth$m_wg_cat <- as.factor(birth$m_wg_cat)
var <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
         "clinega","kotck","pncgov", "rf_db_gest","rf_db_other",
         "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
         "rf_prev_sga", "bc_30d","bc_3090d", "bc_90280d", "firstborn","m_wg_cat",
         "log_mhincome", "log_mhvalue", "percentPoverty",
         "mrace_1", "mrace_2", "mrace_3", "mrace_4")
birth <- birth[ , var, with = F]

birth[, T := bc_30d]
birth[, bc_30d := NULL]
model.num = lm(T~1, data = birth) 
ps.num <- dnorm((birth$T-model.num$fitted)/(summary(model.num))$sigma,0,1)
independent <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
                 "clinega","kotck","pncgov", "rf_db_gest","rf_db_other",
                 "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
                 "rf_prev_sga","firstborn","m_wg_cat",
                 "log_mhincome", "log_mhvalue", "percentPoverty",
                 "mrace_1", "mrace_2", "mrace_3", "mrace_4",
                 "bc_3090d", "bc_90280d")

h2o.init(nthreads = n_cores, min_mem_size = "460G", port = 54345)
##### depth 6 ######
opt_aac_depth6_30d <- NULL
for (rate in col.sample.rate){
  birth.hex <- as.h2o(birth, destination_frame = "birth.hex")
  gbm_30d <- h2o.gbm(y = "T",
                     x = independent,
                     training_frame = birth.hex,
                     ntrees = n.trees, 
                     max_depth = max.depth[1], # change
                     min_rows = min.rows,
                     learn_rate = learn.rate, 
                     col_sample_rate = rate,
                     distribution = "gaussian")
  pred.gbm_30d <- h2o.staged_predict_proba(object = gbm_30d, newdata = birth.hex)
  opt_aac_i  <- optimize(F.aac.iter, interval = c(1, n.trees), data = birth,
                         data.hex = birth.hex,
                         ps.num = ps.num, 
                         ps.model.pred = pred.gbm_30d, # change
                         rep = rep.bootstrap)
  opt_aac_i <- as.data.frame(opt_aac_i)
  opt_aac_depth6_30d <- rbind(opt_aac_depth6_30d, opt_aac_i)
  rm(opt_aac_i)
  h2o.removeAll()
  gc()
}
row.names(opt_aac_depth6_30d) <- col.sample.rate

##### depth 8 ######
opt_aac_depth8_30d <- NULL
# h2o.init(nthreads = n_cores, min_mem_size = "460G")
for (rate in col.sample.rate){
  birth.hex <- as.h2o(birth, destination_frame = "birth.hex")
  gbm_30d <- h2o.gbm(y = "T",
                     x = independent,
                     training_frame = birth.hex,
                     ntrees = n.trees, 
                     max_depth = max.depth[2], # change
                     min_rows = min.rows,
                     learn_rate = learn.rate, 
                     col_sample_rate = rate,
                     distribution = "gaussian")
  pred.gbm_30d <- h2o.staged_predict_proba(object = gbm_30d, newdata = birth.hex)
  opt_aac_i  <- optimize(F.aac.iter, interval = c(1, n.trees), data = birth,
                         data.hex = birth.hex,
                         ps.num = ps.num, 
                         ps.model.pred = pred.gbm_30d, # change
                         rep = rep.bootstrap)
  opt_aac_i <- as.data.frame(opt_aac_i)
  opt_aac_depth8_30d <- rbind(opt_aac_depth8_30d, opt_aac_i)
  rm(opt_aac_i)
  h2o.removeAll()
  gc()
}
row.names(opt_aac_depth8_30d) <- col.sample.rate

###### depth 10 #######
opt_aac_depth10_30d <- NULL
h2o.init(nthreads = n_cores, min_mem_size = "460G")
for (rate in col.sample.rate){
  birth.hex <- as.h2o(birth, destination_frame = "birth.hex")
  gbm_30d <- h2o.gbm(y = "T",
                     x = independent,
                     training_frame = birth.hex,
                     ntrees = n.trees, 
                     max_depth = max.depth[3], # change
                     min_rows = min.rows,
                     learn_rate = learn.rate, 
                     col_sample_rate = rate,
                     distribution = "gaussian")
  pred.gbm_30d <- h2o.staged_predict_proba(object = gbm_30d, newdata = birth.hex)
  opt_aac_i  <- optimize(F.aac.iter, interval = c(1, n.trees), data = birth,
                         data.hex = birth.hex,
                         ps.num = ps.num, 
                         ps.model.pred = pred.gbm_30d, # change
                         rep = rep.bootstrap)
  opt_aac_i <- as.data.frame(opt_aac_i)
  opt_aac_depth10_30d <- rbind(opt_aac_depth10_30d, opt_aac_i)
  rm(opt_aac_i)
  h2o.removeAll()
  gc()
}
row.names(opt_aac_depth10_30d) <- col.sample.rate
h2o.shutdown(prompt = FALSE)

############################## 3.2 bc_3090d ######################################
## load data
birth <- fread(paste0(dir_input, "birth_final.csv"),
               drop = "V1")
birth$year <- as.factor(birth$year)
birth$m_edu <- as.factor(birth$m_edu)
birth$kotck <- as.factor(birth$kotck)
birth$m_wg_cat <- as.factor(birth$m_wg_cat)
var <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
         "clinega","kotck","pncgov", "rf_db_gest","rf_db_other",
         "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
         "rf_prev_sga", "bc_30d","bc_3090d", "bc_90280d", "firstborn","m_wg_cat",
         "log_mhincome", "log_mhvalue", "percentPoverty",
         "mrace_1", "mrace_2", "mrace_3", "mrace_4")
birth <- birth[ , var, with = F]

birth[, T := bc_3090d]
birth[, bc_3090d := NULL]
model.num = lm(T~1, data = birth) 
ps.num <- dnorm((birth$T-model.num$fitted)/(summary(model.num))$sigma,0,1)
independent <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
                 "clinega","kotck","pncgov", "rf_db_gest","rf_db_other",
                 "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
                 "rf_prev_sga","firstborn","m_wg_cat",
                 "log_mhincome", "log_mhvalue", "percentPoverty",
                 "mrace_1", "mrace_2", "mrace_3", "mrace_4",
                 "bc_30d", "bc_90280d")

h2o.init(nthreads = n_cores, min_mem_size = "460G", port = 54345)

##### depth 6 ######
opt_aac_depth6_3090d <- NULL
for (rate in col.sample.rate){
  birth.hex <- as.h2o(birth, destination_frame = "birth.hex")
  gbm_30d <- h2o.gbm(y = "T",
                     x = independent,
                     training_frame = birth.hex,
                     ntrees = n.trees, 
                     max_depth = max.depth[1], # change
                     min_rows = min.rows,
                     learn_rate = learn.rate, 
                     col_sample_rate = rate,
                     distribution = "gaussian")
  pred.gbm_30d <- h2o.staged_predict_proba(object = gbm_30d, newdata = birth.hex)
  opt_aac_i  <- optimize(F.aac.iter, interval = c(1, n.trees), data = birth,
                         data.hex = birth.hex,
                         ps.num = ps.num, 
                         ps.model.pred = pred.gbm_30d, # change
                         rep = rep.bootstrap)
  opt_aac_i <- as.data.frame(opt_aac_i)
  opt_aac_depth6_3090d <- rbind(opt_aac_depth6_3090d, opt_aac_i)
  rm(opt_aac_i)
  h2o.removeAll()
  gc()
}
row.names(opt_aac_depth6_3090d) <- col.sample.rate

##### depth 8 ######
opt_aac_depth8_3090d <- NULL
for (rate in col.sample.rate){
  birth.hex <- as.h2o(birth, destination_frame = "birth.hex")
  gbm_30d <- h2o.gbm(y = "T",
                     x = independent,
                     training_frame = birth.hex,
                     ntrees = n.trees, 
                     max_depth = max.depth[2], # change
                     min_rows = min.rows,
                     learn_rate = learn.rate, 
                     col_sample_rate = rate,
                     distribution = "gaussian")
  pred.gbm_30d <- h2o.staged_predict_proba(object = gbm_30d, newdata = birth.hex)
  opt_aac_i  <- optimize(F.aac.iter, interval = c(1, n.trees), data = birth,
                         data.hex = birth.hex,
                         ps.num = ps.num, 
                         ps.model.pred = pred.gbm_30d, # change
                         rep = rep.bootstrap)
  opt_aac_i <- as.data.frame(opt_aac_i)
  opt_aac_depth8_3090d <- rbind(opt_aac_depth8_3090d, opt_aac_i)
  rm(opt_aac_i)
  h2o.removeAll()
  gc()
}
row.names(opt_aac_depth8_3090d) <- col.sample.rate

##### depth 10 ######
opt_aac_depth10_3090d <- NULL
for (rate in col.sample.rate){
  birth.hex <- as.h2o(birth, destination_frame = "birth.hex")
  gbm_30d <- h2o.gbm(y = "T",
                     x = independent,
                     training_frame = birth.hex,
                     ntrees = n.trees, 
                     max_depth = max.depth[3], # change
                     min_rows = min.rows,
                     learn_rate = learn.rate, 
                     col_sample_rate = rate,
                     distribution = "gaussian")
  pred.gbm_30d <- h2o.staged_predict_proba(object = gbm_30d, newdata = birth.hex)
  opt_aac_i  <- optimize(F.aac.iter, interval = c(1, n.trees), data = birth,
                         data.hex = birth.hex,
                         ps.num = ps.num, 
                         ps.model.pred = pred.gbm_30d, # change
                         rep = rep.bootstrap)
  opt_aac_i <- as.data.frame(opt_aac_i)
  opt_aac_depth10_3090d <- rbind(opt_aac_depth10_3090d, opt_aac_i)
  rm(opt_aac_i)
  h2o.removeAll()
  gc()
}
row.names(opt_aac_depth10_3090d) <- col.sample.rate
h2o.shutdown(prompt = FALSE)

############################## 3.3 bc_90280d ##################################
## load data
birth <- fread(paste0(dir_input, "birth_final.csv"),
               drop = "V1")
birth$year <- as.factor(birth$year)
birth$m_edu <- as.factor(birth$m_edu)
birth$kotck <- as.factor(birth$kotck)
birth$m_wg_cat <- as.factor(birth$m_wg_cat)
var <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
         "clinega","kotck","pncgov", "rf_db_gest","rf_db_other",
         "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
         "rf_prev_sga", "bc_30d","bc_3090d", "bc_90280d", "firstborn","m_wg_cat",
         "log_mhincome", "log_mhvalue", "percentPoverty",
         "mrace_1", "mrace_2", "mrace_3", "mrace_4")
birth <- birth[ , var, with = F]

birth[, T := bc_90280d]
birth[, bc_90280d := NULL]
model.num = lm(T~1, data = birth) 
ps.num <- dnorm((birth$T-model.num$fitted)/(summary(model.num))$sigma,0,1)
independent <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
                 "clinega","kotck","pncgov", "rf_db_gest","rf_db_other",
                 "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
                 "rf_prev_sga","firstborn","m_wg_cat",
                 "log_mhincome", "log_mhvalue", "percentPoverty",
                 "mrace_1", "mrace_2", "mrace_3", "mrace_4",
                 "bc_30d", "bc_3090d")

h2o.init(nthreads = n_cores, min_mem_size = "460G", port = 54345)
##### depth 6 ######
opt_aac_depth6_90280d <- NULL
for (rate in col.sample.rate){
  birth.hex <- as.h2o(birth, destination_frame = "birth.hex")
  gbm_30d <- h2o.gbm(y = "T",
                     x = independent,
                     training_frame = birth.hex,
                     ntrees = n.trees, 
                     max_depth = max.depth[1], # change
                     min_rows = min.rows,
                     learn_rate = learn.rate, 
                     col_sample_rate = rate,
                     distribution = "gaussian")
  pred.gbm_30d <- h2o.staged_predict_proba(object = gbm_30d, newdata = birth.hex)
  opt_aac_i  <- optimize(F.aac.iter, interval = c(1, n.trees), data = birth,
                         data.hex = birth.hex,
                         ps.num = ps.num, 
                         ps.model.pred = pred.gbm_30d, # change
                         rep = rep.bootstrap)
  opt_aac_i <- as.data.frame(opt_aac_i)
  opt_aac_depth6_90280d <- rbind(opt_aac_depth6_90280d, opt_aac_i)
  rm(opt_aac_i)
  h2o.removeAll()
  gc()
}
row.names(opt_aac_depth6_90280d) <- col.sample.rate

##### depth 8 ######
opt_aac_depth8_90280d <- NULL
for (rate in col.sample.rate){
  birth.hex <- as.h2o(birth, destination_frame = "birth.hex")
  gbm_30d <- h2o.gbm(y = "T",
                     x = independent,
                     training_frame = birth.hex,
                     ntrees = n.trees, 
                     max_depth = max.depth[2], # change
                     min_rows = min.rows,
                     learn_rate = learn.rate, 
                     col_sample_rate = rate,
                     distribution = "gaussian")
  pred.gbm_30d <- h2o.staged_predict_proba(object = gbm_30d, newdata = birth.hex)
  opt_aac_i  <- optimize(F.aac.iter, interval = c(1, n.trees), data = birth,
                         data.hex = birth.hex,
                         ps.num = ps.num, 
                         ps.model.pred = pred.gbm_30d, # change
                         rep = rep.bootstrap)
  opt_aac_i <- as.data.frame(opt_aac_i)
  opt_aac_depth8_90280d <- rbind(opt_aac_depth8_90280d, opt_aac_i)
  rm(opt_aac_i)
  h2o.removeAll()
  gc()
}
row.names(opt_aac_depth8_90280d) <- col.sample.rate

##### depth 10 ######
opt_aac_depth10_90280d <- NULL
for (rate in col.sample.rate){
  birth.hex <- as.h2o(birth, destination_frame = "birth.hex")
  gbm_30d <- h2o.gbm(y = "T",
                     x = independent,
                     training_frame = birth.hex,
                     ntrees = n.trees, 
                     max_depth = max.depth[3], # change
                     min_rows = min.rows,
                     learn_rate = learn.rate, 
                     col_sample_rate = rate,
                     distribution = "gaussian")
  pred.gbm_30d <- h2o.staged_predict_proba(object = gbm_30d, newdata = birth.hex)
  opt_aac_i  <- optimize(F.aac.iter, interval = c(1, n.trees), data = birth,
                         data.hex = birth.hex,
                         ps.num = ps.num, 
                         ps.model.pred = pred.gbm_30d, # change
                         rep = rep.bootstrap)
  opt_aac_i <- as.data.frame(opt_aac_i)
  opt_aac_depth10_90280d <- rbind(opt_aac_depth10_90280d, opt_aac_i)
  rm(opt_aac_i)
  h2o.removeAll()
  gc()
}
row.names(opt_aac_depth10_90280d) <- col.sample.rate
h2o.shutdown(prompt = FALSE)

############################## 4. summarize results ###########################
result_30d <-
  rbind(
    opt_aac_depth6_30d %>% mutate(depth = 6,col_rate = col.sample.rate),
    opt_aac_depth8_30d %>% mutate(depth = 8,col_rate = col.sample.rate),
    opt_aac_depth10_30d %>% mutate(depth = 10,col_rate = col.sample.rate))
colnames(result_30d) <- c("optimal_ntrees", "AAC", "depth", "col_rate")

result_3090d <-
  rbind(
    opt_aac_depth6_3090d %>% mutate(depth = 6,col_rate = col.sample.rate),
    opt_aac_depth8_3090d %>% mutate(depth = 8,col_rate = col.sample.rate),
    opt_aac_depth10_3090d %>% mutate(depth = 10,col_rate = col.sample.rate))
colnames(result_3090d) <- c("optimal_ntrees", "AAC", "depth", "col_rate")

result_90280d <-
  rbind(
    opt_aac_depth6_90280d %>% mutate(depth = 6,col_rate = col.sample.rate),
    opt_aac_depth8_90280d %>% mutate(depth = 8,col_rate = col.sample.rate),
    opt_aac_depth10_90280d %>% mutate(depth = 10,col_rate = col.sample.rate))
colnames(result_90280d) <- c("optimal_ntrees", "AAC", "depth", "col_rate")

write.csv(result_30d, "gridSearch_30d.csv")
write.csv(result_3090d, "gridSearch_3090d.csv")
write.csv(result_90280d, "gridSearch_90280d.csv")

gc()
