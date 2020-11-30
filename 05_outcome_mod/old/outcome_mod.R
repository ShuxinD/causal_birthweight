###############################################################################
# Project: Causal black carbon on birth weight in MA                          #
# Code: Outcome regression model with weighting                               #
# Input: clean birth data, bc averages and weights                            #
# Output: BC effects and plots                                                #
# Author: Shuxin Dong                                                         #
# Date: Sep 17, 2020                                                          #
###############################################################################

############################# 0. Setup ########################################
rm(list = ls())
gc()

library(dplyr)

dir_input1 <- "/Users/shuxind/Desktop/BC_birthweight_data/" # change
dir_input2 <- "/Users/shuxind/Desktop/BC_birthweight_analysis/04_GPS_calc_weights/for_time_period/result/"
dir_output <- "/Users/shuxind/Desktop/BC_birthweight_analysis/06_outcome_mod/result/" # change

new_plot_name <- "effect_timeperiods.png"

## load data
birth <- readRDS(paste0(dir_input1, "birth_tfm.rds"))
birth <- birth %>% select(-bwg, -bc_30d)
bc <- readRDS(paste0(dir_input2, "gbm_bc_wt.rds"))
birth_wt <- left_join(bc, birth, by = "uniqueid_yr")
remove(bc)
remove(birth)
gc()


############################# 1. fit outcome regression model #################
# y as continuous birth weight
gbm_bc30d_cont <- lm(bwg ~ bc_30d, data = birth_wt, weights = wt_bc_30d)
summary(gbm_bc30d_cont)
gbm_bc3090d_cont <- lm(bwg ~ bc_3090d, data = birth_wt, weights = wt_bc_3090d)
summary(gbm_bc3090d_cont)
gbm_bc90280d_cont <- lm(bwg ~ bc_90280d, data = birth_wt, weights = wt_bc_90280d)
summary(gbm_bc90280d_cont)
# y as binary low birth weight
gbm_bc30d_binr <- glm(lbw ~ bc_30d, data = birth_wt, weights = wt_bc_30d, 
                      family = binomial(link = "logit"))
summary(gbm_bc30d_binr)
gbm_bc3090d_binr <- glm(lbw ~ bc_3090d, data = birth_wt, weights = wt_bc_3090d, 
                      family = binomial(link = "logit"))
summary(gbm_bc3090d_binr)
gbm_bc90280d_binr <- glm(lbw ~ bc_90280d, data = birth_wt, weights = wt_bc_90280d, 
                       family = binomial(link = "logit"))
summary(gbm_bc90280d_binr)

############################# 2. calculate the results ########################
iqr_bc30d <- IQR(birth_wt$bc_30d)
iqr_bc3090d <- IQR(birth_wt$bc_3090d)
iqr_bc90280d <- IQR(birth_wt$bc_90280d)

iqr_bc30d
iqr_bc3090d
iqr_bc90280d

effect <- rbind(gbm_bc30d_cont$coefficients[2]*iqr_bc30d,
                gbm_bc3090d_cont$coefficients[2]*iqr_bc3090d,
                gbm_bc90280d_cont$coefficients[2]*iqr_bc90280d)
halfCI <- rbind(qnorm(0.975)*summary(gbm_bc30d_cont)$coef[2,2]*iqr_bc30d,
                qnorm(0.975)*summary(gbm_bc3090d_cont)$coef[2,2]*iqr_bc3090d,
                qnorm(0.975)*summary(gbm_bc90280d_cont)$coef[2,2]*iqr_bc90280d)
result1 <- cbind(effect,effect-halfCI,effect+halfCI)
colnames(result1) <- c("effect for IQR", "lower CI", "upper CI")
rownames(result1) <- c("bc_30d", "bc_3090d", "bc_90280d")
print(result1)
# effect for IQR  lower CI  upper CI
# bc_30d         -12.75474 -14.09636 -11.41313
# bc_3090d       -14.53860 -15.94141 -13.13578
# bc_90280d      -11.55651 -12.94948 -10.16354

ORs <- rbind(exp(gbm_bc30d_binr$coefficients[2]*iqr_bc30d),
             exp(gbm_bc3090d_binr$coefficients[2]*iqr_bc3090d),
             exp(gbm_bc90280d_binr$coefficients[2]*iqr_bc90280d))
uppCI <- rbind(exp((gbm_bc30d_binr$coefficients[2]+qnorm(0.975)*summary(gbm_bc30d_binr)$coef[2,2])*iqr_bc30d),
               exp((gbm_bc3090d_binr$coefficients[2]+qnorm(0.975)*summary(gbm_bc3090d_binr)$coef[2,2])*iqr_bc3090d),
               exp((gbm_bc90280d_binr$coefficients[2]+qnorm(0.975)*summary(gbm_bc90280d_binr)$coef[2,2])*iqr_bc90280d))
lowCI <- rbind(exp((gbm_bc30d_binr$coefficients[2]-qnorm(0.975)*summary(gbm_bc30d_binr)$coef[2,2])*iqr_bc30d),
               exp((gbm_bc3090d_binr$coefficients[2]-qnorm(0.975)*summary(gbm_bc3090d_binr)$coef[2,2])*iqr_bc3090d),
               exp((gbm_bc90280d_binr$coefficients[2]-qnorm(0.975)*summary(gbm_bc90280d_binr)$coef[2,2])*iqr_bc90280d))
result2 <- cbind(ORs, lowCI, uppCI)
colnames(result2) <- c("OR for IQR", "lower CI", "upper CI")
rownames(result2) <- c("bc_30d", "bc_3090d", "bc_90280d")
print(result2)
# OR for IQR lower CI upper CI
# bc_30d      1.068878 1.044129 1.094213
# bc_3090d    1.091048 1.062901 1.119940
# bc_90280d   1.086250 1.059747 1.113415

############################# 3. draw the plots ##############################
bcdays <- c(1:3)
diff_lowCI <- result2[,1]-result2[,2]
diff_uppCI <- result2[,3]-result2[,1]

par(mfrow=c(1,2), mar=c(4,4,4,4))
plot(bcdays, effect,
     xlim = range(c(0.5,3.5)),
     ylim = range(c(-20, 0.2)),
     xaxt = "n",
     pch=19, xlab="Moving averages of BC exposure \n (a)", ylab="Change of birth weight with 95% CI (g)"
)
arrows(bcdays, effect-halfCI, bcdays, effect+halfCI, length=0.05, angle=90, code=3)
abline(h=0,col = 1, lty = 3)
axis(1, at=1:3, labels=c("30d","30-90d","90-280d"))

plot(bcdays, ORs,
     xlim = range(c(0.5,3.5)),
     ylim=range(c(0.99, 1.2)),
     xaxt = "n",
     pch=19, xlab="Moving averages of BC exposure \n (b)", ylab="OR of LBW with 95% CI"
)
arrows(bcdays, ORs-diff_lowCI, bcdays, ORs+diff_uppCI, length=0.05, angle=90, code=3)
abline(h=1,col = 1, lty = 3)
axis(1, at=1:3, labels=c("30d","30-90d","90-280d"))
