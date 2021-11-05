#' Project: causal_birthweight
#' Code: outcome model and check balance again
#' Input: "MAbirth_for_analyses.csv" "IPWs_all_1104.csv"
#' Output: ...
#' Author: Shuxin Dong
#' First create date: 2021-01-27

## setup ----
rm(list = ls())
gc()

library(data.table)
library(mgcv)
library(sandwich)

setwd("/media/gate/Shuxin/")
dir_birth <- "/media/qnap3/Shuxin/airPollution_MAbirth/"
dir_ipw <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/2ipw/"

## load birth data ----
birth <- fread(paste0(dir_birth, "MAbirth_for_analyses.csv"))
# names(birth)
birth$year <- as.factor(birth$year)

## load IPWs ----
IPWs <- fread(paste0(dir_ipw, "IPWs_all_1104.csv"))

## merge ----
dt <- merge(birth, IPWs, by="uniqueid_yr")
names(dt)

## IQR ----
IQRs <- data.table(bc_all=IQR(dt[,bc_all]),
                   no2_all=IQR(dt[,no2_all]))

## bwg ~ spline ----
dir_splines <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/4mainEffects/"
exposures <- c("bc_all", "no2_all")

pdf(file = paste0(dir_splines,"splines.pdf"))
par(mfrow = c(2,1))
for (exposures_i in exposures){
  splines.c <- gam(bwg ~ s(get(exposures_i), bs = "cr", fx=F), family = gaussian, 
                   data = dt, weights = get(paste0("ipw_gbm_", exposures_i)))
  plot(splines.c, rug = T)
  cat("plot ", exposures_i, "\n")
}
dev.off()

## bwg ~ exposure ----
dir_bwg <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/4mainEffects/"
exposures <- c("bc_all", "no2_all")

results_glm.bwg <- NULL
for (exposures_i in exposures){
  results_i <- NULL
  cat("model ", exposures_i, "\n")
  glm.bwg <- glm(bwg ~ get(exposures_i), 
                 family = gaussian(link = "identity"), data = dt, 
                 weights = get(paste0("ipw_gbm_", exposures_i)))
  results_i <- c(coef(glm.bwg)[2],
                 sqrt(vcovHC(glm.bwg)[2,2]),
                 coef(glm.bwg)[2]*IQRs[,get(exposures_i)],
                 (coef(glm.bwg)[2]-qnorm(0.975)*sqrt(vcovHC(glm.bwg)[2,2]))*IQRs[,get(exposures_i)],
                 (coef(glm.bwg)[2]+qnorm(0.975)*sqrt(vcovHC(glm.bwg)[2,2]))*IQRs[,get(exposures_i)])
  results_glm.bwg <- rbind(results_glm.bwg, results_i)
}
rownames(results_glm.bwg) <- exposures
colnames(results_glm.bwg) <- c("coef", "robust.se", "diff", "diff_lci", "diff_uci")
print(results_glm.bwg)
write.csv(results_glm.bwg, file = paste0(dir_bwg, "results_all_bwg.csv"))

## lbw ~ exposure ----
dir_lbw <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/4mainEffects/"
exposures <- c("bc_all", "no2_all")

results_glm.lbw <- NULL
for (exposures_i in exposures){
  results_i <- NULL
  cat("model ", exposures_i, "\n")
  glm.lbw <- glm(lbw ~ get(exposures_i), family = quasibinomial(link = "logit"), data = dt, weights = get(paste0("ipw_gbm_", exposures_i)))
  results_i <- c(coef(glm.lbw)[2],
                 sqrt(vcovHC(glm.lbw)[2,2]),
                 exp(coef(glm.lbw)[2]*IQRs[,get(exposures_i)]),
                 exp((coef(glm.lbw)[2]-qnorm(0.975)*sqrt(vcovHC(glm.lbw)[2,2]))*IQRs[,get(exposures_i)]),
                 exp((coef(glm.lbw)[2]+qnorm(0.975)*sqrt(vcovHC(glm.lbw)[2,2]))*IQRs[,get(exposures_i)]))
  results_glm.lbw <- rbind(results_glm.lbw, results_i)
}
rownames(results_glm.lbw) <- exposures
colnames(results_glm.lbw) <- c("coef", "robust.se", "OR", "OR_lci", "OR_uci")
print(results_glm.lbw)
write.csv(results_glm.lbw, file = paste0(dir_bwg, "results_glm_lbw.csv"))

## export balance ----
dir_balance <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/3balancePlots/"
ps_exposures <- c("bc_all","no2_all")
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

response <- "no2_all" # change everytime
predictor <- c(ps_vars, ps_exposures[ps_exposures!=response])

ac <- matrix(NA,length(predictor),3)
for (j in 1:length(predictor)){ 
  ac[j,1] <-  ifelse(!is.factor(birth[,..predictor][[j]]), 
                     cor(birth[,get(response)], birth[,..predictor][[j]], method = "pearson"),
                     polyserial(birth[,get(response)], birth[,..predictor][[j]]))
  ac[j,2] <-  ifelse(!is.factor(birth[,..predictor][[j]]), 
                     weightedCorr(birth[,get(response)], birth[,..predictor][[j]], method = "pearson", weights = dt[,get(paste0("ipw_gbm_",response))]),
                     weightedCorr(birth[,get(response)], birth[,..predictor][[j]], method = "polyserial", weights = dt[,get(paste0("ipw_gbm_",response))]))
  ac[j,3] <-  ifelse(!is.factor(birth[,..predictor][[j]]), 
                     weightedCorr(birth[,get(response)], birth[,..predictor][[j]], method = "pearson", weights = dt[,get(paste0("ipw_glm_",response))]),
                     weightedCorr(birth[,get(response)], birth[,..predictor][[j]], method = "polyserial", weights = dt[,get(paste0("ipw_glm_",response))]))
  cat("finish column", j,"/", length(predictor), "\n")
}

balance <- data.frame(name.x=predictor,
                      corr=c(ac[,1],ac[,2],ac[,3]),
                      weighted = c(rep("unweighted", length(predictor)),rep("weighted_gbm", length(predictor)), rep("weighted_glm", length(predictor))))

pdf(file = paste0(dir_balance, response, "_balance.pdf"))
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
dev.off()

