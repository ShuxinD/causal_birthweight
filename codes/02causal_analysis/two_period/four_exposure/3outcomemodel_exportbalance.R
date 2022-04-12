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
dir_ipw <- "/media/qnap3/Shuxin/airPollution_MAbirth/data/ipw_two_period/four_exposure/"

## load birth data ----
birth <- fread(paste0(dir_birth, "MAbirth_for_analyses.csv"))
# names(birth)
birth$year <- as.factor(birth$year)

## load IPWs ----
IPWs <- fread(paste0(dir_ipw, "IPWs_all_truncated.csv"))

## merge ----
dt <- merge(birth, IPWs, by="uniqueid_yr")
names(dt)

## IQR ----
IQRs <- data.table(bc_30d=IQR(dt[,bc_30d]),
                   bc_30280d=IQR(dt[,bc_30280d]),
                   no2_30d=IQR(dt[,no2_30d]),
                   no2_30280d=IQR(dt[,no2_30280d]))
MEDIANs_overall <- data.table(bc_30d=median(dt[,bc_30d]),
                      bc_30280d=median(dt[,bc_30280d]),
                      no2_30d=median(dt[,no2_30d]),
                      no2_30280d=median(dt[,no2_30280d]))
MEDIANs_normal <- data.table(bc_30d=median(dt[lbw==0,bc_30d]),
                              bc_30280d=median(dt[lbw==0,bc_30280d]),
                              no2_30d=median(dt[lbw==0,no2_30d]),
                              no2_30280d=median(dt[lbw==0,no2_30280d]))
MEDIANs_low <- data.table(bc_30d=median(dt[lbw==1,bc_30d]),
                             bc_30280d=median(dt[lbw==1,bc_30280d]),
                             no2_30d=median(dt[lbw==1,no2_30d]),
                             no2_30280d=median(dt[lbw==1,no2_30280d]))
quantile(dt[,bc_30280d])
quantile(dt[lbw==0, bc_30280d])
quantile(dt[lbw==1, bc_30280d])
quantile(dt[,no2_30d])
quantile(dt[lbw==0, no2_30d])
quantile(dt[lbw==1, no2_30d])
## bwg ~ spline ----
dir_splines <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/2mainEffects/two_period/four_exposure/"
exposures <- c("bc_30d", "bc_30280d", "no2_30d", "no2_30280d")

pdf(file = paste0(dir_splines,"splines.pdf"))
par(mfrow = c(2,2))
for (exposures_i in exposures){
  splines.c <- gam(bwg ~ s(get(exposures_i), bs = "cr", fx=F), family = gaussian, 
                   data = dt, weights = dt[,get(paste0("ipw_gbm_", exposures_i))])
  plot(splines.c, rug = T)
  cat("plot ", exposures_i, "\n")
}
dev.off()

## bwg ~ exposure ----
dir_bwg <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/2mainEffects/two_period/four_exposure/"
exposures <- c("bc_30d", "bc_30280d", "no2_30d", "no2_30280d")

results_bwg <- NULL
for (exposures_i in exposures){
  results_i <- NULL
  cat("model ", exposures_i, "\n")
  mod <- glm(bwg ~ get(exposures_i), 
             family = gaussian(link = "identity"), data = dt, 
             weights = get(paste0("ipw_gbm_", exposures_i)))
  results_i <- c(coef(mod)[2],
                 sqrt(vcovHC(mod)[2,2]),
                 coef(mod)[2]*IQRs[,get(exposures_i)],
                 (coef(mod)[2]-qnorm(0.975)*sqrt(vcovHC(mod)[2,2]))*IQRs[,get(exposures_i)],
                 (coef(mod)[2]+qnorm(0.975)*sqrt(vcovHC(mod)[2,2]))*IQRs[,get(exposures_i)])
  results_bwg <- rbind(results_bwg, results_i)
}
rownames(results_bwg) <- exposures
colnames(results_bwg) <- c("coef", "robust.se", "diff", "diff_lci", "diff_uci")
print(results_bwg)
write.csv(results_bwg, file = paste0(dir_bwg, "bwg_gbm.csv"))

results_bwg <- NULL
for (exposures_i in exposures){
  results_i <- NULL
  cat("model ", exposures_i, "\n")
  mod <- glm(bwg ~ get(exposures_i), 
             family = gaussian(link = "identity"), data = dt, 
             weights = get(paste0("ipw_glm_", exposures_i)))
  results_i <- c(coef(mod)[2],
                 sqrt(vcovHC(mod)[2,2]),
                 coef(mod)[2]*IQRs[,get(exposures_i)],
                 (coef(mod)[2]-qnorm(0.975)*sqrt(vcovHC(mod)[2,2]))*IQRs[,get(exposures_i)],
                 (coef(mod)[2]+qnorm(0.975)*sqrt(vcovHC(mod)[2,2]))*IQRs[,get(exposures_i)])
  results_bwg <- rbind(results_bwg, results_i)
}
rownames(results_bwg) <- exposures
colnames(results_bwg) <- c("coef", "robust.se", "diff", "diff_lci", "diff_uci")
print(results_bwg)
write.csv(results_bwg, file = paste0(dir_bwg, "bwg_glm.csv"))

## lbw ~ exposure ----
dir_lbw <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/2mainEffects/two_period/four_exposure/"
exposures <- c("bc_30d", "bc_30280d", "no2_30d", "no2_30280d")

results_lbw <- NULL
for (exposures_i in exposures){
  results_i <- NULL
  cat("model ", exposures_i, "\n")
  mod <- glm(lbw ~ get(exposures_i), family = quasibinomial(link = "logit"), data = dt, 
             weights = get(paste0("ipw_gbm_", exposures_i)))
  results_i <- c(coef(mod)[2],
                 sqrt(vcovHC(mod)[2,2]),
                 exp(coef(mod)[2]*IQRs[,get(exposures_i)]),
                 exp((coef(mod)[2]-qnorm(0.975)*sqrt(vcovHC(mod)[2,2]))*IQRs[,get(exposures_i)]),
                 exp((coef(mod)[2]+qnorm(0.975)*sqrt(vcovHC(mod)[2,2]))*IQRs[,get(exposures_i)]))
  results_lbw <- rbind(results_lbw, results_i)
}
rownames(results_lbw) <- exposures
colnames(results_lbw) <- c("coef", "robust.se", "OR", "OR_lci", "OR_uci")
print(results_lbw)
write.csv(results_lbw, file = paste0(dir_lbw, "lbw_gbm.csv"))

results_lbw <- NULL
for (exposures_i in exposures){
  results_i <- NULL
  cat("model ", exposures_i, "\n")
  mod <- glm(lbw ~ get(exposures_i), family = quasibinomial(link = "logit"), data = dt, 
             weights = get(paste0("ipw_glm_", exposures_i)))
  results_i <- c(coef(mod)[2],
                 sqrt(vcovHC(mod)[2,2]),
                 exp(coef(mod)[2]*IQRs[,get(exposures_i)]),
                 exp((coef(mod)[2]-qnorm(0.975)*sqrt(vcovHC(mod)[2,2]))*IQRs[,get(exposures_i)]),
                 exp((coef(mod)[2]+qnorm(0.975)*sqrt(vcovHC(mod)[2,2]))*IQRs[,get(exposures_i)]))
  results_lbw <- rbind(results_lbw, results_i)
}
rownames(results_lbw) <- exposures
colnames(results_lbw) <- c("coef", "robust.se", "OR", "OR_lci", "OR_uci")
print(results_lbw)
write.csv(results_lbw, file = paste0(dir_lbw, "lbw_glm.csv"))

## export balance ----
library(ggplot2)
library(stringr)
dir_balance <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/3balancePlots/two_period/four_exposure/"
ps_exposures <-  c("bc_30d", "bc_30280d", "no2_30d", "no2_30280d")
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

response <- "no2_30280d" # change everytime
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

balance <- data.table(name.x=predictor,
                      corr=c(ac[,1],ac[,2],ac[,3]),
                      weighted = c(rep("unweighted", length(predictor)),rep("weighted_gbm", length(predictor)), rep("weighted_glm", length(predictor))))

balance[,.(AC=mean(abs(corr))), by = weighted]
wrap_label <- function(x){str_wrap(x, width = 10)}
balance_plot <- ggplot(balance[weighted!="weighted_glm",], aes(x = name.x, y = corr, group = weighted, colour = weighted)) +
  # geom_line(aes(colour = weighted), size = 0.5, linetype = "dashed")+
  geom_point(size = 2) +
  geom_line(size = 0.5) +
  geom_hline(yintercept=0, size=0.2) +
  geom_hline(yintercept=-0.1, size=0.2, linetype = "dashed") +
  geom_hline(yintercept=0.1, size=0.2, linetype = "dashed") +
  # geom_hline(yintercept=-0.2, size=0.1, linetype = "dashed") +
  # geom_hline(yintercept=0.2, size=0.1, linetype = "dashed") +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("year" = "Birth year (cat.)",
                              "sex" = "Female baby (b.)",
                              "married" = "Married mother (b.)",
                              "mage" = "Mother's age (cont.)",
                              "clinega" = "Clinical gestational age (cont.)",
                              "pncgov" = "Goverment support for prenatal care (b.)",
                              "rf_db_gest" =  "Gestational diabetes (b.)",
                              "rf_db_other" = "Chronic diabetes (b.)",
                              "rf_hbp_chronic" = "Chronic hypertension (b.)",
                              "rf_hbp_pregn" = "Gestational hypertention (b.)",
                              "rf_cervix" = "Incompetent cervix (b.)",
                              "rf_prev_4kg" = "With previous >4kg infant (b.)",
                              "rf_prev_sga" = "With previous SGA infant (b.)",
                              "smoker_ddp" = "Smoke during pregnancy (b.)",
                              "mrace_1" = "White mother (b.)",
                              "mrace_2" = "Black mother (b.)",
                              "mrace_3" = "Asian/Pacific Islander mother (b.)",
                              "m_edu_1" = "< high school (b.)",
                              "m_edu_2" = "high school/G.E.D. (b.)",
                              "m_edu_3" = "some college (b.)",
                              "m_edu_4" = "bachelor (b.)",
                              "kotck_1" = "APNCU index = 1 (b.)",
                              "kotck_2" = "APNCU index = 2 (b.)",
                              "kotck_3" = "APNCU index = 3 (b.)",
                              "kotck_4" = "APNCU index = 4 (b.)",
                              "m_wg_cat_1" = "Maternal weight change <0 lb. (b.)",
                              "m_wg_cat_2" = "Maternal weight change [0,15) lb. (b.)",
                              "m_wg_cat_3" = "Maternal weight change [15-25) lb. (b.)",
                              "m_wg_cat_4" = "Maternal weight change [25-36) lb. (b.)",
                              "percentPoverty" = "% Poverty at residence address (cont.)",
                              "firstborn" = "First-born infant (b.)",
                              "no2_30d" = "Average of daily NO[2] exposure over 30d* (cont.)",
                              "bc_30d" = "Average of daily BC exposure over 30d* (cont.)",
                              "no2_30280d" = "Average of daily NO[2] exposure over 30-280d* (cont.)",
                              "bc_30280d" = "Average of daily BC exposure over30-280d* (cont.)")) + 
  scale_colour_discrete(breaks=c("unweighted", "weighted_gbm"),
                        labels=c("Unweighted correlation", "Weighted correlation")) +
  coord_flip() + 
  theme(legend.position = "bottom", legend.title = element_blank(),
        text = element_text(size=14))
balance_plot

pdf(file = paste0(dir_balance, response, ".pdf"))
balance_plot
dev.off()
