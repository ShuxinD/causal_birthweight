###############################################################################
# Project: Causal black carbon and NO2 on birth weight in MA                  #
# Code: Outcome regression model with weights                                 #
# Input: "MAbirth_for_analyses.csv" birth data                                #
# Input: "MAbirth_ipw.csv" inverse-probability weights                        #
# Output: BC effects and plots                                                #
# Author: Shuxin Dong                                                         #
# Date: 2021-01-27                                                            #
###############################################################################

############################# 0. Setup ########################################
rm(list = ls())
gc()

library(data.table)
library(mgcv)
library(sandwich)

setwd("/media/gate/Shuxin")
dir_data <- "/media/qnap3/Shuxin/airPollution_MAbirth/"
dir_output_mod <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/4mainEffects/"

############################ 1. Load data #####################################
## birth data
birth <- fread(paste0(dir_data, "MAbirth_for_analyses.csv"))
birth[,':='(year = as.factor(year),
            m_edu = as.factor(m_edu),
            kotck = as.factor(kotck),
            m_wg_cat = as.factor(m_wg_cat))]

var <- c("uniqueid_yr", "bwg", "lbw",
         "bc_all", "no2_all")
         # "bc_30d","bc_3090d", "bc_90280d", 
         # "no2_30d", "no2_3090d", "no2_90280d")
sub_birth <- birth[ , var, with = F]
head(sub_birth)

## inverse-prob weights
# ipw <- fread(paste0(dir_data, "MAbirth_ipw.csv"),
#              select = c("uniqueid_yr",
#                         "bc_30d.wt", "bc_3090d.wt", "bc_90280d.wt",
#                         "no2_30d.wt", "no2_3090d.wt", "no2_90280d.wt"))
# head(ipw)
ipw <- fread(paste0(dir_data, "MAbirth_ipw_2.csv"),
             select = c("uniqueid_yr",
                        "bc_all.wt", "bc_all_mac.wt",
                        "no2_all.wt", "no2_all_mac.wt"))
head(ipw)

## merge together
dt <- merge(sub_birth, ipw, by = "uniqueid_yr")
head(dt)

rm(sub_birth)
rm(ipw)
gc()
gc()

############################# 2. splines ######################################
pdf(file = paste0(dir_output_mod,"splines.pdf"))
par(mfrow = c(2,3))

splines.c <- gam(bwg ~ s(bc_30d, bs = "cr"), family = gaussian, 
                    data = dt, weights = bc_30d.wt)
plot(splines.c)

splines.c <- gam(bwg ~ s(bc_3090d, bs = "cr"), family = gaussian, # change
                 data = dt, weights = bc_3090d.wt) # change
plot(splines.c)

splines.c <- gam(bwg ~ s(bc_90280d, bs = "cr"), family = gaussian, # change
                 data = dt, weights = bc_90280d.wt) # change
plot(splines.c)

splines.c <- gam(bwg ~ s(no2_30d, bs = "cr"), family = gaussian, # change
                 data = dt, weights = no2_30d.wt) # change
plot(splines.c)

splines.c <- gam(bwg ~ s(no2_3090d, bs = "cr"), family = gaussian, # change
                 data = dt, weights = no2_3090d.wt) # change
plot(splines.c)

splines.c <- gam(bwg ~ s(no2_90280d, bs = "cr"), family = gaussian, # change
                 data = dt, weights = no2_90280d.wt) # change
plot(splines.c)

dev.off()

######################## 3. outcome regression mod ############################
bc_280d <- (dt$bc_30d*30 + dt$bc_3090d*60 + dt$bc_90280d*190)/280
no2_280d <- (dt$no2_30d*30 + dt$no2_3090d*60 + dt$no2_90280d*190)/280
IQRs <- data.table(IQR(dt$bc_30d), IQR(dt$bc_3090d), IQR(dt$bc_90280d), IQR(bc_280d),
              IQR(dt$no2_30d), IQR(dt$no2_3090d), IQR(dt$no2_90280d), IQR(no2_280d))
colnames(IQRs) <- c("bc_30d", "bc_3090d", "bc_90280d", "bc_280d",
                    "no2_30d", "no2_3090d", "no2_90280d", "no2_280d")
IQRs
fwrite(IQRs, paste0(dir_output_mod, "IQRs.csv"))

# library(survey)
# ####################### 3.1 outcome as cont. birth weight ####################
# design_bc_30d <- svydesign(ids=~1,weights=~bc_30d.wt, data=dt)
# bwg_bc_30d <- svyglm(formula = bwg ~ bc_30d, design = design_bc_30d,
#                      family = gaussian)
# summary(bwg_bc_30d)
# 
# design_bc_3090d <- svydesign(ids=~1,weights=~bc_3090d.wt, data=dt)
# bwg_bc_3090d <- svyglm(formula = bwg ~ bc_3090d, design = design_bc_3090d, 
#                      family = gaussian)
# summary(bwg_bc_3090d)
# 
# design_bc_90280d <- svydesign(ids=~1,weights=~bc_90280d.wt, data=dt)
# bwg_bc_90280d <- svyglm(formula = bwg ~ bc_90280d, design = design_bc_90280d, 
#                      family = gaussian)
# summary(bwg_bc_90280d)
# 
# design_no2_30d <- svydesign(ids=~1,weights=~no2_30d.wt, data=dt)
# bwg_no2_30d <- svyglm(formula = bwg ~ no2_30d, design = design_no2_30d, 
#                      family = gaussian)
# summary(bwg_no2_30d)
# 
# design_no2_3090d <- svydesign(ids=~1,weights=~no2_3090d.wt, data=dt)
# bwg_no2_3090d <- svyglm(formula = bwg ~ no2_3090d, design = design_no2_3090d, 
#                        family = gaussian)
# summary(bwg_no2_3090d)
# 
# design_no2_90280d <- svydesign(ids=~1,weights=~no2_90280d.wt, data=dt)
# bwg_no2_90280d <- svyglm(formula = bwg ~ no2_90280d, design = design_no2_90280d, 
#                         family = gaussian)
# summary(bwg_no2_90280d)
bwg_bc_30d <- glm(bwg ~ bc_30d, family = gaussian, data = dt, 
                  weights = bc_30d.wt)
summary(bwg_bc_30d)
bwg_bc_3090d <- glm(bwg ~ bc_3090d, family = gaussian, data = dt,
                    weights = bc_3090d.wt)
summary(bwg_bc_3090d)
bwg_bc_90280d <- glm(bwg ~ bc_90280d, family = gaussian, data = dt, 
                     weights = bc_90280d.wt)
summary(bwg_bc_90280d)

bwg_no2_30d <- glm(bwg ~ no2_30d, family = gaussian, data = dt, 
                   weights = no2_30d.wt)
summary(bwg_no2_30d)
bwg_no2_3090d <- glm(bwg ~ no2_3090d, family = gaussian, data = dt,
                     weights = no2_3090d.wt)
summary(bwg_no2_3090d)
bwg_no2_90280d <- glm(bwg ~ no2_90280d, family = gaussian, data = dt, 
                      weights = no2_90280d.wt)
summary(bwg_no2_90280d)

#####
bwg_bc_all <- glm(bwg ~ bc_all, family = gaussian, data = dt, 
                      weights = bc_all.wt)
summary(bwg_bc_all)
bwg_bc_all_mac <- glm(bwg ~ bc_all, family = gaussian, data = dt, 
                  weights = bc_all_mac.wt)
summary(bwg_bc_all_mac)
bwg_no2_all <- glm(bwg ~ no2_all, family = gaussian, data = dt, 
                  weights = bc_all.wt)
summary(bwg_no2_all)
bwg_no2_all_mac <- glm(bwg ~ no2_all, family = gaussian, data = dt, 
                  weights = no2_all_mac.wt)
summary(bwg_no2_all_mac)

#################### 3.2 outcome as binary low birth weight ###########
# lbw_bc_30d <- svyglm(formula = lbw ~ bc_30d, design = design_bc_30d,
#                      family = quasibinomial(link = "logit"))
# summary(lbw_bc_30d)
# 
# lbw_bc_3090d <- svyglm(formula = lbw ~ bc_3090d, design = design_bc_3090d,
#                      family = quasibinomial(link = "logit"))
# summary(lbw_bc_3090d)
# 
# lbw_bc_90280d <- svyglm(formula = lbw ~ bc_90280d, design = design_bc_90280d,
#                        family = quasibinomial(link = "logit"))
# summary(lbw_bc_90280d)
# 
# lbw_no2_30d <- svyglm(formula = lbw ~ no2_30d, design = design_no2_30d,
#                      family = quasibinomial(link = "logit"))
# summary(lbw_no2_30d)
# 
# lbw_no2_3090d <- svyglm(formula = lbw ~ no2_3090d, design = design_no2_3090d,
#                        family = quasibinomial(link = "logit"))
# summary(lbw_no2_3090d)
# 
# lbw_no2_90280d <- svyglm(formula = lbw ~ no2_90280d, design = design_no2_90280d,
#                         family = quasibinomial(link = "logit"))
# summary(lbw_no2_90280d)
lbw_bc_30d <- glm(lbw ~ bc_30d, family = quasibinomial(link = "logit"), data = dt, 
                  weights = bc_30d.wt)
summary(lbw_bc_30d)
lbw_bc_3090d <- glm(lbw ~ bc_3090d, family = quasibinomial(link = "logit"), data = dt, 
                    weights = bc_3090d.wt)
summary(lbw_bc_3090d)
lbw_bc_90280d <- glm(lbw ~ bc_90280d, family = quasibinomial(link = "logit"), data = dt, 
                     weights = bc_90280d.wt)
summary(lbw_bc_90280d)

lbw_no2_30d <- glm(lbw ~ no2_30d, family = quasibinomial(link = "logit"), data = dt, 
                   weights = no2_30d.wt)
summary(lbw_no2_30d)
lbw_no2_3090d <- glm(lbw ~ no2_3090d, family = quasibinomial(link = "logit"), data = dt, 
                     weights = no2_3090d.wt)
summary(lbw_no2_3090d)
lbw_no2_90280d <- glm(lbw ~ no2_90280d, family = quasibinomial(link = "logit"), data = dt, 
                      weights = no2_90280d.wt)
summary(lbw_no2_90280d)

#####
lbw_bc_all <- glm(lbw ~ bc_all, family = quasibinomial(link = "logit"), data = dt, 
                  weights = bc_all.wt)
summary(lbw_bc_all)
lbw_bc_all_mac <- glm(lbw ~ bc_all, family = quasibinomial(link = "logit"), data = dt, 
                      weights = bc_all_mac.wt)
summary(lbw_bc_all_mac)
lbw_no2_all <- glm(lbw ~ no2_all, family = quasibinomial(link = "logit"), data = dt, 
                   weights = bc_all.wt)
summary(lbw_no2_all)
lbw_no2_all_mac <- glm(lbw ~ no2_all, family = quasibinomial(link = "logit"), data = dt, 
                   weights = no2_all_mac.wt)
summary(lbw_no2_all_mac)

############################# 4. results & plots ##############################
############################# 4.1 results tables ##############################
################################# bwg as outcome ##############################
## bc_all -----
effect <- rbind(coef(bwg_bc_all)[2]*IQR(dt[,bc_all]),
                coef(bwg_bc_all_mac)[2]*IQR(dt[,bc_all]))
halfCI <- rbind(qnorm(0.975)*sqrt(vcovHC(bwg_bc_all)[2,2])*IQR(dt[,bc_all]),
                qnorm(0.975)*sqrt(vcovHC(bwg_bc_all_mac)[2,2])*IQR(dt[,bc_all]))
result_bwg_bc <- cbind(effect,effect-halfCI,effect+halfCI)
colnames(result_bwg_bc) <- c("effect for IQR", "lower CI", "upper CI")
rownames(result_bwg_bc) <- c("bc_all", "bc_all_mac")
print(result_bwg_bc)
write.csv(result_bwg_bc, paste0(dir_output_mod, "result_bwg_bc_all.csv"))
## no2_all -----
effect <- rbind(coef(bwg_no2_all)[2]*IQR(dt[,no2_all]),
                coef(bwg_no2_all_mac)[2]*IQR(dt[,no2_all]))
halfCI <- rbind(qnorm(0.975)*sqrt(vcovHC(bwg_no2_all)[2,2])*IQR(dt[,no2_all]),
                qnorm(0.975)*sqrt(vcovHC(bwg_no2_all_mac)[2,2])*IQR(dt[,no2_all]))
result_bwg_no2 <- cbind(effect,effect-halfCI,effect+halfCI)
colnames(result_bwg_no2) <- c("effect for IQR", "lower CI", "upper CI")
rownames(result_bwg_no2) <- c("no2_all", "no2_all_mac")
print(result_bwg_no2)
write.csv(result_bwg_no2, paste0(dir_output_mod, "result_bwg_no2_all.csv"))

## bc_all_mac -----
effect <- rbind(coef(bwg_bc_all_mac)[2]*IQR(dt[,bc_all]),
                coef(bwg_bc_all_mac)[2]*IQR(dt[,bc_all]))
halfCI <- rbind(qnorm(0.975)*sqrt(vcovHC(bwg_bc_all)[2,2])*IQR(dt[,bc_all]),
                qnorm(0.975)*sqrt(vcovHC(bwg_bc_all_mac)[2,2])*IQR(dt[,bc_all]))
result_bwg_bc <- cbind(effect,effect-halfCI,effect+halfCI)
colnames(result_bwg_bc) <- c("effect for IQR", "lower CI", "upper CI")
rownames(result_bwg_bc) <- c("bc_all", "bc_all_mac")
print(result_bwg_bc)
write.csv(result_bwg_bc, paste0(dir_output_mod, "result_bwg_bc_all.csv"))
## no2_all -----
effect <- rbind(coef(bwg_no2_all)[2]*IQR(dt[,no2_all]),
                coef(bwg_no2_all_mac)[2]*IQR(dt[,no2_all]))
halfCI <- rbind(qnorm(0.975)*sqrt(vcovHC(bwg_no2_all)[2,2])*IQR(dt[,no2_all]),
                qnorm(0.975)*sqrt(vcovHC(bwg_no2_all_mac)[2,2])*IQR(dt[,no2_all]))
result_bwg_no2 <- cbind(effect,effect-halfCI,effect+halfCI)
colnames(result_bwg_no2) <- c("effect for IQR", "lower CI", "upper CI")
rownames(result_bwg_no2) <- c("no2_all", "no2_all_mac")
print(result_bwg_no2)
write.csv(result_bwg_no2, paste0(dir_output_mod, "result_bwg_no2_all.csv"))

## bc
effect <- rbind(coef(bwg_bc_30d)[2]*IQRs$bc_30d,
                coef(bwg_bc_3090d)[2]*IQRs$bc_3090d,
                coef(bwg_bc_90280d)[2]*IQRs$bc_90280d,
                (coef(bwg_bc_30d)[2]+coef(bwg_bc_3090d)[2]+coef(bwg_bc_90280d)[2])*IQRs$bc_280d)
halfCI <- rbind(qnorm(0.975)*sqrt(vcovHC(bwg_bc_30d)[2,2])*IQRs$bc_30d,
                qnorm(0.975)*sqrt(vcovHC(bwg_bc_3090d)[2,2])*IQRs$bc_3090d,
                qnorm(0.975)*sqrt(vcovHC(bwg_bc_90280d)[2,2])*IQRs$bc_90280d,
                qnorm(0.975)*sqrt(vcovHC(bwg_bc_30d)[2,2]+vcovHC(bwg_bc_3090d)[2,2]+vcovHC(bwg_bc_90280d)[2,2])*IQRs$bc_280d)

result_bwg_bc <- cbind(effect,effect-halfCI,effect+halfCI)
colnames(result_bwg_bc) <- c("effect for IQR", "lower CI", "upper CI")
rownames(result_bwg_bc) <- c("bc_30d", "bc_3090d", "bc_90280d", "bc_280d")
print(result_bwg_bc)
write.csv(result_bwg_bc, paste0(dir_output_mod, "result_bwg_bc.csv"))

## no2
effect <- rbind(coef(bwg_no2_30d)[2]*IQRs$no2_30d,
                coef(bwg_no2_3090d)[2]*IQRs$no2_3090d,
                coef(bwg_no2_90280d)[2]*IQRs$no2_90280d,
                (coef(bwg_no2_30d)[2]+coef(bwg_no2_3090d)[2]+coef(bwg_no2_90280d)[2])*IQRs$no2_280d)
halfCI <- rbind(qnorm(0.975)*sqrt(vcovHC(bwg_no2_30d)[2,2])*IQRs$no2_30d,
                qnorm(0.975)*sqrt(vcovHC(bwg_no2_3090d)[2,2])*IQRs$no2_3090d,
                qnorm(0.975)*sqrt(vcovHC(bwg_no2_90280d)[2,2])*IQRs$no2_90280d,
                qnorm(0.975)*sqrt(vcovHC(bwg_no2_30d)[2,2]+vcovHC(bwg_no2_3090d)[2,2]+vcovHC(bwg_no2_90280d)[2,2])*IQRs$no2_280d)

result_bwg_no2 <- cbind(effect,effect-halfCI,effect+halfCI)
colnames(result_bwg_no2) <- c("effect for IQR", "lower CI", "upper CI")
rownames(result_bwg_no2) <- c("no2_30d", "no2_3090d", "no2_90280d", "no2_280d")
print(result_bwg_no2)
write.csv(result_bwg_no2, paste0(dir_output_mod, "result_bwg_no2.csv"))

################################# lbw as outcome ##############################

## bc_all -----
effect <- rbind(coef(lbw_bc_all)[2]*IQR(dt[,bc_all]),
                coef(lbw_bc_all_mac)[2]*IQR(dt[,bc_all]))
halfCI <- rbind(qnorm(0.975)*sqrt(vcovHC(lbw_bc_all)[2,2])*IQR(dt[,bc_all]),
                qnorm(0.975)*sqrt(vcovHC(lbw_bc_all_mac)[2,2])*IQR(dt[,bc_all]))
result_lbw_bc <- exp(cbind(effect,effect-halfCI,effect+halfCI))
colnames(result_lbw_bc) <- c("effect for IQR", "lower CI", "upper CI")
rownames(result_lbw_bc) <- c("bc_all", "bc_all_mac")
print(result_lbw_bc)
write.csv(result_lbw_bc, paste0(dir_output_mod, "result_lbw_bc_all.csv"))
## no2_all -----
effect <- rbind(coef(lbw_no2_all)[2]*IQR(dt[,no2_all]),
                coef(lbw_no2_all_mac)[2]*IQR(dt[,no2_all]))
halfCI <- rbind(qnorm(0.975)*sqrt(vcovHC(lbw_no2_all)[2,2])*IQR(dt[,no2_all]),
                qnorm(0.975)*sqrt(vcovHC(lbw_no2_all_mac)[2,2])*IQR(dt[,no2_all]))
result_lbw_no2 <- exp(cbind(effect,effect-halfCI,effect+halfCI))
colnames(result_lbw_no2) <- c("effect for IQR", "lower CI", "upper CI")
rownames(result_lbw_no2) <- c("no2_all", "no2_all_mac")
print(result_lbw_no2)
write.csv(result_lbw_no2, paste0(dir_output_mod, "result_lbw_no2_all.csv"))

## bc_all_mac -----
effect <- rbind(coef(lbw_bc_all_mac)[2]*IQR(dt[,bc_all]),
                coef(lbw_bc_all_mac)[2]*IQR(dt[,bc_all]))
halfCI <- rbind(qnorm(0.975)*sqrt(vcovHC(lbw_bc_all)[2,2])*IQR(dt[,bc_all]),
                qnorm(0.975)*sqrt(vcovHC(lbw_bc_all_mac)[2,2])*IQR(dt[,bc_all]))
result_lbw_bc <- exp(cbind(effect,effect-halfCI,effect+halfCI))
colnames(result_lbw_bc) <- c("effect for IQR", "lower CI", "upper CI")
rownames(result_lbw_bc) <- c("bc_all", "bc_all_mac")
print(result_lbw_bc)
write.csv(result_lbw_bc, paste0(dir_output_mod, "result_lbw_bc_all.csv"))
## no2_all -----
effect <- rbind(coef(lbw_no2_all)[2]*IQR(dt[,no2_all]),
                coef(lbw_no2_all_mac)[2]*IQR(dt[,no2_all]))
halfCI <- rbind(qnorm(0.975)*sqrt(vcovHC(lbw_no2_all)[2,2])*IQR(dt[,no2_all]),
                qnorm(0.975)*sqrt(vcovHC(lbw_no2_all_mac)[2,2])*IQR(dt[,no2_all]))
result_lbw_no2 <- exp(cbind(effect,effect-halfCI,effect+halfCI))
colnames(result_lbw_no2) <- c("effect for IQR", "lower CI", "upper CI")
rownames(result_lbw_no2) <- c("no2_all", "no2_all_mac")
print(result_lbw_no2)
write.csv(result_lbw_no2, paste0(dir_output_mod, "result_lbw_no2_all.csv"))


## bc
effect <- rbind(coef(lbw_bc_30d)[2]*IQRs$bc_30d,
                coef(lbw_bc_3090d)[2]*IQRs$bc_3090d,
                coef(lbw_bc_90280d)[2]*IQRs$bc_90280d,
                (coef(lbw_bc_30d)[2]+coef(lbw_bc_3090d)[2]+coef(lbw_bc_90280d)[2])*IQRs$bc_280d)
halfCI <- rbind(qnorm(0.975)*sqrt(vcovHC(lbw_bc_30d)[2,2])*IQRs$bc_30d,
                qnorm(0.975)*sqrt(vcovHC(lbw_bc_3090d)[2,2])*IQRs$bc_3090d,
                qnorm(0.975)*sqrt(vcovHC(lbw_bc_90280d)[2,2])*IQRs$bc_90280d,
                qnorm(0.975)*sqrt(vcovHC(lbw_bc_30d)[2,2]+vcovHC(lbw_bc_3090d)[2,2]+vcovHC(lbw_bc_90280d)[2,2])*IQRs$bc_280d)

result_lbw_bc <- cbind(exp(effect), exp(effect-halfCI), exp(effect+halfCI))
colnames(result_lbw_bc) <- c("OR for IQR", "lower CI", "upper CI")
rownames(result_lbw_bc) <- c("bc_30d", "bc_3090d", "bc_90280d", "bc_280d")
print(result_lbw_bc)
write.csv(result_lbw_bc, paste0(dir_output_mod, "result_lbw_bc.csv"))

## no2
effect <- rbind(coef(lbw_no2_30d)[2]*IQRs$no2_30d,
                coef(lbw_no2_3090d)[2]*IQRs$no2_3090d,
                coef(lbw_no2_90280d)[2]*IQRs$no2_90280d,
                (coef(lbw_no2_30d)[2]+coef(lbw_no2_3090d)[2]+coef(lbw_no2_90280d)[2])*IQRs$no2_280d)
halfCI <- rbind(qnorm(0.975)*sqrt(vcovHC(lbw_no2_30d)[2,2])*IQRs$no2_30d,
                qnorm(0.975)*sqrt(vcovHC(lbw_no2_3090d)[2,2])*IQRs$no2_3090d,
                qnorm(0.975)*sqrt(vcovHC(lbw_no2_90280d)[2,2])*IQRs$no2_90280d,
                qnorm(0.975)*sqrt(vcovHC(lbw_no2_30d)[2,2]+vcovHC(lbw_no2_3090d)[2,2]+vcovHC(lbw_no2_90280d)[2,2])*IQRs$no2_280d)

result_lbw_no2 <- cbind(exp(effect), exp(effect-halfCI), exp(effect+halfCI))
colnames(result_lbw_no2) <- c("OR for IQR", "lower CI", "upper CI")
rownames(result_lbw_no2) <- c("no2_30d", "no2_3090d", "no2_90280d", "no2_280d")
print(result_lbw_no2)
write.csv(result_lbw_no2, paste0(dir_output_mod, "result_lbw_no2.csv"))

############################ 4.2 results plots ################################
################################# bwg as outcome ##############################
library(ggplot2)
effect_bc <- read.csv(paste0(dir_output_mod,"result_bwg_bc.csv"), row.names = 1)[-4,]
effect_no2 <- read.csv(paste0(dir_output_mod,"result_bwg_no2.csv"), row.names = 1)[-4,]

bcdays <- c(1:3)
bcdays <- factor(bcdays,
                 levels = c(1,2,3),
                 labels = c("0-30d", "31-90d", "91-280d"))
effect_bc <- cbind(bcdays, Exposure = rep("BC",3), effect_bc)
effect_no2 <- cbind(bcdays, Exposure = rep("NO[2]",3), effect_no2)

effect_bwg_all <- rbind(effect_bc, effect_no2)

## independent plot
p_bwg_bc <- 
  ggplot(effect_bc, aes(x = bcdays, y = effect.for.IQR)) + 
  geom_point(size = 0.01) +
  ylim(5,-20) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "red", size=1) +
  # geom_ribbon(aes(ymin = lower.CI, ymax = upper.CI)) +
  geom_pointrange(aes(ymin = lower.CI, ymax = upper.CI)) +
  xlab("Moving averages of BC exposure \n prior to the delivery date") +
  ylab("Change of birth weight with 95% CI (g)\n for one IQR increase in BC")
p_bwg_bc

p_bwg_no2 <- 
  ggplot(effect_no2, aes(x = bcdays, y = effect.for.IQR)) + 
  geom_point(size = 0.01) +
  ylim(5,-20) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "red", size=1) +
  geom_pointrange(aes(ymin = lower.CI, ymax = upper.CI)) +
  xlab("Moving averages of NO[2] exposure \n prior to the delivery date") +
  ylab("Change of birth weight with 95% CI (g)\n for one IQR increase in NO[2]")
p_bwg_no2

## two pollutants together
p_bwg_all <- 
  ggplot(effect_bwg_all, aes(x = bcdays, y = effect.for.IQR, colour = Exposure)) + 
  geom_errorbar(aes(ymin = lower.CI, ymax = upper.CI, colour = Exposure, width = 0.05)) +
  geom_point(size = 1) +
  ylim(5,-20) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = 1, size = 0.5) +
  xlab("Time windows prior to the delivery date \n of BC and NO[2] moving averages") +
  ylab("Change of birth weight with 95% CI (g)\n for one IQR increase in moving averages") +
  theme(legend.position="top")
p_bwg_all

################################# lbw as outcome ##############################
library(ggplot2)
effect_bc <- read.csv(paste0(dir_output_mod,"result_lbw_bc.csv"), row.names = 1)[-4,]
effect_no2 <- read.csv(paste0(dir_output_mod,"result_lbw_no2.csv"), row.names = 1)[-4,]

bcdays <- c(1:3)
bcdays <- factor(bcdays,
                 levels = c(1,2,3),
                 labels = c("0-30d", "31-90d", "91-280d"))
effect_bc <- cbind(bcdays, Exposure = rep("BC",3), effect_bc)
effect_no2 <- cbind(bcdays, Exposure = rep("NO[2]",3), effect_no2)

effect_lbw_all <- rbind(effect_bc, effect_no2)

## independent plots
p_lbw_bc <- 
  ggplot(effect_bc, aes(x = bcdays, y = OR.for.IQR)) + 
  geom_point(size = 1) +
  ylim(0.90,1.20) +
  geom_hline(yintercept=1, linetype="dashed", 
             color = "red", size=1) +
  geom_pointrange(aes(ymin = lower.CI, ymax = upper.CI)) +
  xlab("Moving averages of BC exposure \n prior to the delivery date") +
  ylab("OR of LBW \n for one IQR increase in BC with 95% CI")
p_lbw_bc

p_lbw_no2 <- 
  ggplot(effect_no2, aes(x = bcdays, y = OR.for.IQR)) + 
  geom_point(size = 1) +
  ylim(0.90,1.20) +
  geom_hline(yintercept=1, linetype="dashed", 
             color = "red", size=1) +
  geom_pointrange(aes(ymin = lower.CI, ymax = upper.CI)) +
  xlab("Moving averages of NO[2] exposure \n prior to the delivery date") +
  ylab("OR of LBW \n for one IQR increase in NO[2] with 95% CI")
p_lbw_no2

## two pollutants together
p_lbw_all <- 
  ggplot(effect_lbw_all, aes(x = bcdays, y = OR.for.IQR, colour = Exposure)) + 
  geom_errorbar(aes(ymin = lower.CI, ymax = upper.CI, colour = Exposure, width = 0.05)) +
  geom_point(size = 1) +
  ylim(0.90,1.20) +
  geom_hline(yintercept=1, linetype="dashed", 
             color = 1, size = 0.5) +
  xlab("Time windows prior to the delivery date \n of BC and NO[2] moving averages") +
  ylab("OR of LBW with 95% CI \n for one IQR increase in moving averages") +
  theme(legend.position="top")
p_lbw_all

# p_legend <- get_legend(p_lbw_all)
# plot_grid(p_bwg_all, p_lbw_all, labels = "AUTO")

library(lemon)
pdf(paste0(dir_output_mod, "results_all.pdf"))
grid_arrange_shared_legend(p_bwg_all, p_lbw_all, position = "top")
dev.off()

