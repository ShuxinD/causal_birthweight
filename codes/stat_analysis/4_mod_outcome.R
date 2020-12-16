###############################################################################
# Project: Causal black carbon on birth weight in MA                          #
# Code: Outcome regression model with weights                                 #
# Input: "birth_all.csv"                                                      #
# Output: BC effects and plots                                                #
# Author: Shuxin Dong                                                         #
# Date: Nov 18, 2020                                                          #
###############################################################################

############################# 0. Setup ########################################
rm(list = ls())
gc()

library(dplyr)
library(data.table)
library(mgcv)

# dir_input <- "/Users/shuxind/Desktop/BC_birthweight_data/"
# setwd("/media/gate/Shuxin")
dir_input <- "/media/gate/Shuxin/"
dir_output <- "/media/gate/Shuxin/"

## load data
birth <- fread(paste0(dir_input, "birth_all.csv"),
                   select = c("bwg", 
                              "bc_30d", "bc_3090d", "bc_90280d",
                              "bc_30d.wt.t", "bc_3090d.wt.t", "bc_90280d.wt.t"))
birth$lbw <- 0
birth$lbw[birth$bwg<2500] <- 1

############################# 1. splines ######################################
############################# 1.1 bc_30d ######################################
spline_30d.c <- gam(bwg ~ s(bc_30d), family = gaussian, 
                    data = birth, weights = bc_30d.wt.t)
spline_30d.l <- gam(lbw ~ s(bc_30d), family = binomial, 
                    data = birth, weights = bc_30d.wt.t)

pdf(file = paste0(dir_output, "spline_30d_c.pdf"))
plot(spline_30d.c)
dev.off()

pdf(file = paste0(dir_output, "spline_30d_l.pdf"))
plot(spline_30d.l)
dev.off()

############################# 1.2 bc_3090d ####################################
spline_3090d.c <- gam(bwg ~ s(bc_3090d), family = gaussian, 
                    data = birth, weights = bc_3090d.wt.t)
spline_3090d.l <- gam(lbw ~ s(bc_3090d), family = binomial, 
                    data = birth, weights = bc_3090d.wt.t)

pdf(file = paste0(dir_output, "spline_3090d_c.pdf"))
plot(spline_3090d.c)
dev.off()

pdf(file = paste0(dir_output, "spline_3090d_l.pdf"))
plot(spline_3090d.l)
dev.off()

############################# 1.3 bc_90280d ###################################
spline_90280d.c <- gam(bwg ~ s(bc_90280d), family = gaussian, 
                      data = birth, weights = bc_90280d.wt.t)
spline_90280d.l <- gam(lbw ~ s(bc_90280d), family = binomial, 
                      data = birth, weights = bc_90280d.wt.t)

pdf(file = paste0(dir_output, "spline_90280d_c.pdf"))
plot(spline_90280d.c)
dev.off()

pdf(file = paste0(dir_output, "spline_90280d_l.pdf"))
plot(spline_90280d.l)
dev.off()

############################# 2. calculate results ############################
iqr_30d <- IQR(birth$bc_30d)
iqr_3090d <- IQR(birth$bc_3090d)
iqr_90280d <- IQR(birth$bc_90280d)
print(iqr_30d)
print(iqr_3090d)
print(iqr_90280d)
############################# 2.1 outcome as birth weight cont.################
bwg_30d <- lm(bwg ~ bc_30d, data = birth, weights = bc_30d.wt.t)
summary(bwg_30d)
bwg_3090d <- lm(bwg ~ bc_3090d, data = birth, weights = bc_3090d.wt.t)
summary(bwg_3090d)
bwg_90280d <- lm(bwg ~ bc_90280d, data = birth, weights = bc_90280d.wt.t)
summary(bwg_90280d)
############################# 2.2 outcome as low birth weight binary###########
lbw_30d <- glm(lbw ~ bc_30d, data = birth, weights = bc_30d.wt.t, 
                      family = binomial(link = "logit"))
summary(lbw_30d)
lbw_3090d <- glm(lbw ~ bc_3090d, data = birth, weights = bc_3090d.wt.t, 
               family = binomial(link = "logit"))
summary(lbw_3090d)
lbw_90280d <- glm(lbw ~ bc_90280d, data = birth, weights = bc_90280d.wt.t, 
               family = binomial(link = "logit"))
summary(lbw_90280d)
############################# 2.3 results table ##############################
effect <- rbind(bwg_30d$coefficients[2]*iqr_30d,
                bwg_3090d$coefficients[2]*iqr_3090d,
                bwg_90280d$coefficients[2]*iqr_90280d)
halfCI <- rbind(qnorm(0.975)*summary(bwg_30d)$coef[2,2]*iqr_30d,
                qnorm(0.975)*summary(bwg_3090d)$coef[2,2]*iqr_3090d,
                qnorm(0.975)*summary(bwg_90280d)$coef[2,2]*iqr_90280d)
result1 <- cbind(effect,effect-halfCI,effect+halfCI)
colnames(result1) <- c("effect for IQR", "lower CI", "upper CI")
rownames(result1) <- c("bc_30d", "bc_3090d", "bc_90280d")
print(result1)
write.csv(result1, file = paste0(dir_output, "result1_effect.csv"))

ORs <- rbind(exp(lbw_30d$coefficients[2]*iqr_30d),
             exp(lbw_3090d$coefficients[2]*iqr_3090d),
             exp(lbw_90280d$coefficients[2]*iqr_90280d))
uppCI <- rbind(exp((lbw_30d$coefficients[2]+qnorm(0.975)*summary(lbw_30d)$coef[2,2])*iqr_30d),
               exp((lbw_3090d$coefficients[2]+qnorm(0.975)*summary(lbw_3090d)$coef[2,2])*iqr_3090d),
               exp((lbw_90280d$coefficients[2]+qnorm(0.975)*summary(lbw_90280d)$coef[2,2])*iqr_90280d))
lowCI <- rbind(exp((lbw_30d$coefficients[2]-qnorm(0.975)*summary(lbw_30d)$coef[2,2])*iqr_30d),
               exp((lbw_3090d$coefficients[2]-qnorm(0.975)*summary(lbw_3090d)$coef[2,2])*iqr_3090d),
               exp((lbw_90280d$coefficients[2]-qnorm(0.975)*summary(lbw_90280d)$coef[2,2])*iqr_90280d))
result2 <- cbind(ORs, lowCI, uppCI)
colnames(result2) <- c("OR for IQR", "lower CI", "upper CI")
rownames(result2) <- c("bc_30d", "bc_3090d", "bc_90280d")
print(result2)
write.csv(result2, file = paste0(dir_output, "result2_ORs.csv"))

############################# 3. plot results #################################
library(ggplot2)
library(cowplot)

effect <- read.csv("/Users/shuxind/Documents/GitHub/causal_BC_birthweight/05_outcome_mod/result1_effect.csv",
                   row.names = 1)
ORs <- read.csv("/Users/shuxind/Documents/GitHub/causal_BC_birthweight/05_outcome_mod/result2_ORs.csv",
                row.names = 1)
bcdays <- c(1:3)
bcdays <- factor(bcdays,
                 levels = c(1,2,3),
                 labels = c("0-30d", "31-90d", "91-280d"))
effect <- cbind(bcdays, effect)
ORs <- cbind(bcdays, ORs)

## plot
p_e <- 
  ggplot(effect, aes(x = bcdays, y = effect.for.IQR)) + 
  ylim(0,-12) +
  geom_pointrange(aes(ymin = lower.CI, ymax = upper.CI)) +
  xlab("Moving averages of BC exposure \n prior to the delivery date") +
  ylab("Change of birth weight \n for one IQR increase in BC with 95% CI (g)")

p_OR <- 
ggplot(ORs, aes(x = bcdays, y = OR.for.IQR)) + 
  geom_pointrange(aes(ymin = lower.CI, ymax = upper.CI)) +
  geom_hline(yintercept=1, linetype="dashed", 
             color = "red", size=1) + 
  xlab("Moving averages of BC exposure \n prior to the delivery date") +
  ylab("OR of LBW \n for one IQR increase in BC with 95% CI")

plot_grid(p_e, p_OR, labels = "AUTO")

pdf(file = "/Users/shuxind/Documents/GitHub/causal_BC_birthweight/results/resultPlot.pdf",
    width = 10.5)
plot_grid(p_e, p_OR, labels = "AUTO")
dev.off()

