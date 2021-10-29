#' Project: causal_birthweight
#' Code: outcome model and check balance again
#' Input: "MAbirth_for_analyses.csv" "IPWs1028.csv"
#' Output: weights dataset "MAbirth_ipw.csv" with ID and weights
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

## load birth registry ----
birth <- fread(paste0(dir_birth, "MAbirth_for_analyses.csv"))
# names(birth)
birth$year <- as.factor(birth$year)

## load IPWs ----
IPWs <- fread(paste0(dir_ipw, "IPWs1028.csv"))

## merge ----
dt <- merge(birth, IPWs, by="uniqueid_yr")
names(dt)

## IQR ----
IQRs <- data.table(bc_30d=IQR(dt[,bc_30d]),
                   bc_3090d=IQR(dt[,bc_3090d]),
                   bc_90280d=IQR(dt[,bc_90280d]),
                   no2_30d=IQR(dt[,no2_30d]),
                   no2_3090d=IQR(dt[,no2_3090d]),
                   no2_90280d=IQR(dt[,no2_90280d]))

## bwg ~ spline ----
dir_splines <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/4mainEffects/"
exposures <- c("bc_30d", "bc_3090d", "bc_90280d","no2_30d", "no2_3090d","no2_90280d")

pdf(file = paste0(dir_splines,"splines.pdf"))
par(mfrow = c(2,3))
for (exposures_i in exposures){
  splines.c <- gam(bwg ~ s(get(exposures_i), bs = "cr", fx=F) + as.factor(year), family = gaussian, 
                   data = dt, weights = get(paste0("ipw_", exposures_i)))
  plot(splines.c)
  cat("plot ", exposures_i, "\n")
}
dev.off()

## bwg ~ exposure ----
dir_bwg <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/4mainEffects/"
exposures <- c("bc_30d", "bc_3090d", "bc_90280d","no2_30d", "no2_3090d","no2_90280d")

results_glm.bwg <- NULL
for (exposures_i in exposures){
  results_i <- NULL
  cat("model ", exposures_i, "\n")
  glm.bwg <- glm(bwg ~ get(exposures_i) + as.factor(year), family = gaussian(link = "identity"), data = dt, weights = get(paste0("ipw_", exposures_i)))
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
write.csv(results_glm.bwg, file = paste0(dir_bwg, "results_glm_bwg.csv"))

## lbw ~ exposure ----
dir_lbw <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/4mainEffects/"
exposures <- c("bc_30d", "bc_3090d", "bc_90280d","no2_30d", "no2_3090d","no2_90280d")

results_glm.lbw <- NULL
for (exposures_i in exposures){
  results_i <- NULL
  cat("model ", exposures_i, "\n")
  glm.lbw <- glm(lbw ~ get(exposures_i) + as.factor(year), family = quasibinomial(link = "logit"), data = dt, weights = get(paste0("ipw_", exposures_i)))
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
ps_exposures <- c("bc_30d","bc_3090d", "bc_90280d", 
                  "no2_30d", "no2_3090d", "no2_90280d")
ps_vars <- c(# "year",
             "sex","married","mage", "cigdpp","cigddp",
             "clinega","pncgov", 
             # "rf_db_gest","rf_db_other", "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg", "rf_prev_sga", 
             # "smoker_ddp", "smoker_dpp",
             "mrace_1", "mrace_2", "mrace_3", "mrace_4",
             "m_edu_1", "m_edu_2", "m_edu_3","m_edu_4", "m_edu_5",
             "kotck_1","kotck_2","kotck_3","kotck_4",
             "m_wg_cat_1","m_wg_cat_2","m_wg_cat_3","m_wg_cat_4","m_wg_cat_5",
             # "log_mhvalue", "log_mhincome",
             "percentPoverty", "firstborn")

ps_exposures_i <- "no2_90280d"

response <- ps_exposures_i
predictor <- c(ps_vars, ps_exposures[ps_exposures!=response])
ipw <- dt[, get(paste0("ipw_", response))]
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
cat("plot", response, "\n")
pdf(file = paste0(dir_balance, response, "_balance.pdf"))
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
dev.off()
  
## bootstrapping ----
coef_bootstrap <- NULL
for (rep_i in 1:100) {
  bootstrap_sample <- NULL
  cat("bootstrapping ", rep_i, "\n")
  set.seed(rep_i)
  bo <- sample(dim(dt)[1], floor(dim(dt)[1]^0.7))
  bootstrap_sample <- dt[bo,]
  glm.bwg <- glm(bwg ~ bc_3090d + as.factor(year), family = gaussian(link = "identity"), data = bootstrap_sample, weights = ipw_bc_3090d)
  coef_bootstrap<- c(coef_bootstrap,coef(glm.bwg)[2])
}
par(mfrow=c(1,1))
hist(coef_bootstrap)
qnorm(0.975)*sd(coef_bootstrap)*sqrt(floor(dim(dt)[1]^0.7)/dim(dt)[1])



glm.bwg <- glm(bwg ~ get(exposures_i) + as.factor(year), family = gaussian(link = "identity"), data = dt, weights = get(paste0("ipw_", exposures_i)))
results_i <- c(coef(glm.bwg)[2],
               sqrt(vcovHC(glm.bwg)[2,2]),
               coef(glm.bwg)[2]*IQRs[,get(exposures_i)],
               (coef(glm.bwg)[2]-qnorm(0.975)*sqrt(vcovHC(glm.bwg)[2,2]))*IQRs[,get(exposures_i)],
               (coef(glm.bwg)[2]+qnorm(0.975)*sqrt(vcovHC(glm.bwg)[2,2]))*IQRs[,get(exposures_i)])
results_glm.bwg <- rbind(results_glm.bwg, results_i)
