###############################################################################
# Project: Causal black carbon and NO2 on birth weight in MA                  #
# Code: balance checking                                                      #
# Input: "MAbirth_for_analyses.csv" birth data                                #
# Input: "MAbirth_ipw.csv" inverse-probability weights                        #
# Output: balance plots                                                       #
# Author: Shuxin Dong                                                         #
# Date: 2021-01-27                                                            #
###############################################################################

############################# 0. Setup ########################################
rm(list = ls())
gc()

library(data.table)
library(wCorr)
library(ggplot2)

setwd("/media/gate/Shuxin/")
dir_data <- "/media/qnap3/Shuxin/airPollution_MAbirth/"
dir_balancePlots <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/3balancePlots/"

############################# 1. load data ####################################
## birth data
rawbirth <- fread(paste0(dir_data, "MAbirth_for_analyses.csv"))
names(rawbirth)
# [1] "uniqueid_yr"    "year"           "sex"            "married"        "mage"          
# [6] "mrace"          "m_edu"          "cigdpp"         "cigddp"         "clinega"       
# [11] "kotck"          "pncgov"         "bwg"            "rf_db_gest"     "rf_db_other"   
# [16] "rf_hbp_chronic" "rf_hbp_pregn"   "rf_cervix"      "rf_prev_4kg"    "rf_prev_sga"   
# [21] "mhincome"       "mhvalue"        "percentPoverty" "bc_30d"         "bc_3090d"      
# [26] "bc_90280d"      "no2_30d"        "no2_3090d"      "no2_90280d"     "bc_all"        
# [31] "no2_all"        "lbw"            "firstborn"      "m_wg_cat"       "smoker_ddp"    
# [36] "smoker_dpp"     "mrace_1"        "mrace_2"        "mrace_3"        "mrace_4"       
# [41] "log_mhincome"   "log_mhvalue"   
summary(rawbirth)
rawbirth$year <- as.factor(rawbirth$year)
rawbirth$m_edu <- as.factor(rawbirth$m_edu)
rawbirth$kotck <- as.factor(rawbirth$kotck)
rawbirth$m_wg_cat <- as.factor(rawbirth$m_wg_cat)

var <- c("uniqueid_yr",
         "year","sex","married","mage","m_edu", "cigdpp","cigddp",
         "clinega", "kotck","pncgov", "rf_db_gest","rf_db_other",
         "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
         "rf_prev_sga", "percentPoverty",
         "bc_all", "no2_all",
         "bc_30d","bc_3090d", "bc_90280d", 
         "no2_30d", "no2_3090d", "no2_90280d",
         "firstborn","m_wg_cat", "smoker_ddp", "smoker_dpp",
         "mrace_1", "mrace_2", "mrace_3", "mrace_4",
         "log_mhincome", "log_mhvalue")
birth <- rawbirth[ , var, with = F]

## inverse-prob weights
ipw <- fread(paste0(dir_data, "MAbirth_temp_bc.csv"))
head(ipw)

## merge together
birth_all <- merge(birth, ipw, by = "uniqueid_yr")
names(birth_all)

rm(birth)
rm(ipw)
gc()
gc()

############################# 2. subset data ##################################
## only binary and continuous variables could be checked, ordinal cannot
## select columns including binary and continuous, and i-p weights
var <- c("sex", "married", "mage",
         "mrace_1", "mrace_2", "mrace_3", "mrace_4",
         "smoker_dpp","smoker_ddp", "cigdpp", "cigddp",
         "clinega", "firstborn",
         "rf_db_gest","rf_db_other",
         "rf_hbp_pregn","rf_hbp_chronic", "rf_cervix","rf_prev_4kg",
         "rf_prev_sga",
         "pncgov",
         "log_mhincome", "log_mhvalue", "percentPoverty",
         "bc_30d","bc_3090d", "bc_90280d",
         "no2_30d", "no2_3090d", "no2_90280d",
         "bc_30d.wt", "bc_3090d.wt", "bc_90280d.wt")
# var <- c("sex", "married", "mage", 
#          "mrace_1", "mrace_2", "mrace_3", "mrace_4",
#          "smoker_dpp","smoker_ddp", "cigdpp", "cigddp",
#          "clinega", "firstborn", 
#          "rf_db_gest","rf_db_other",
#          "rf_hbp_pregn","rf_hbp_chronic", "rf_cervix","rf_prev_4kg",
#          "rf_prev_sga", 
#          "pncgov", 
#          "log_mhincome", "log_mhvalue", "percentPoverty",
#          "bc_all", "no2_all",
#          "bc_all.wt", "bc_all_mac.wt",
#          "no2_all.wt", "no2_all_mac.wt")
balanceALL <- birth_all[ , var, with = F]

description <- c("New-born sex = Female","Married","Mother's age",
                 "Maternal race = white", "Maternal race = black", 
                 "Maternal race = Asian/Pacific Islander", 
                 "Maternal race = others",
                 "Smoking before pregnancy",
                 "Smoking during pregnancy",
                 "# of daily cigarettes smoking before pregnancy",
                 "# of daily cigarettes smoking during pregnancy",
                 "Clinical gestational age",
                 "First-born child",
                 "Gestational diabetes",
                 "Maternal diabetes",
                 "Gestational hypertention",
                 "Maternal chronic hypertension", 
                 "Maternal incompetent cervix",
                 "Mother with previous infant over 4000 g",
                 "Mother with previous small-for-gestational-age infant", 
                 "Government support for prenatal care",
                 "Log. median household income (census-tract level)", 
                 "Log. median value of house (census-tract level)", 
                 "% Poverty (census-tract level)",
                 "Average of daily BC exposure levels over 0-280 days prior to the delivery date",
                 "Average of daily NO[2] exposure levels over 0-280 days prior to the delivery date", 
                 "bc_all.wt", "bc_all_mac.wt" , 
                 "no2_all.wt", "no2_all_mac.wt")
describe_x <- data.table(description, var)
describe_x

###################### 3. get balance results #################################
library(dplyr)
## 3.01 bc_all ----
x <- balanceALL %>% select(-bc_all, # change
                           -bc_all.wt, -bc_all_mac.wt, 
                           -no2_all.wt, -no2_all_mac.wt)
name.x <- names(x)
unweightedCor <- matrix(NA,dim(x)[2],1)
weightedCor <- matrix(NA,dim(x)[2],1)

for (j in 1:dim(x)[2]){
  unweightedCor[j,1] = stats::cor(balanceALL$bc_all, x[[j]], method = "pearson") # change
}
for (j in 1:dim(x)[2]){
  weightedCor[j,1] = weightedCorr(balanceALL$bc_all, x[[j]],method = "pearson", # change
                                  weights = balanceALL$bc_all.wt) # change
}

unwtCorr <- data.table(name.x, unweightedCor, weighted = "unweighted")
names(unwtCorr)
colnames(unwtCorr)[2] <- "corr"

wtCorr <- data.table(name.x, weightedCor, weighted="weighted")
names(wtCorr)
colnames(wtCorr)[2] <- "corr"

balance <- rbind(unwtCorr,wtCorr)
balance <- merge(balance, describe_x, by.x = "name.x", by.y = "var")
head(balance)

fwrite(balance, paste0(dir_balancePlots, "balance_bc_all.csv")) # change

## 3.02 no2_all ----
x <- balanceALL %>% select(-no2_all, # change
                           -bc_all.wt, -bc_all_mac.wt, 
                           -no2_all.wt, -no2_all_mac.wt)
name.x <- names(x)
unweightedCor <- matrix(NA,dim(x)[2],1)
weightedCor <- matrix(NA,dim(x)[2],1)

for (j in 1:dim(x)[2]){
  unweightedCor[j,1] = stats::cor(balanceALL$no2_all, x[[j]], method = "pearson") # change
}
for (j in 1:dim(x)[2]){
  weightedCor[j,1] = weightedCorr(balanceALL$no2_all, x[[j]],method = "pearson", # change
                                  weights = balanceALL$no2_all.wt) # change
}

unwtCorr <- data.table(name.x, unweightedCor, weighted = "unweighted")
names(unwtCorr)
colnames(unwtCorr)[2] <- "corr"

wtCorr <- data.table(name.x, weightedCor, weighted="weighted")
names(wtCorr)
colnames(wtCorr)[2] <- "corr"

balance <- rbind(unwtCorr,wtCorr)
balance <- merge(balance, describe_x, by.x = "name.x", by.y = "var")
head(balance)

fwrite(balance, paste0(dir_balancePlots, "balance_no2_all.csv")) # change

## 3.011 bc_all_mac ----
x <- balanceALL %>% select(-bc_all, # change
                           -bc_all.wt, -bc_all_mac.wt, 
                           -no2_all.wt, -no2_all_mac.wt)
name.x <- names(x)
unweightedCor <- matrix(NA,dim(x)[2],1)
weightedCor <- matrix(NA,dim(x)[2],1)

for (j in 1:dim(x)[2]){
  unweightedCor[j,1] = stats::cor(balanceALL$bc_all, x[[j]], method = "pearson") # change
}
for (j in 1:dim(x)[2]){
  weightedCor[j,1] = weightedCorr(balanceALL$bc_all, x[[j]],method = "pearson", # change
                                  weights = balanceALL$bc_all_mac.wt) # change
}

unwtCorr <- data.table(name.x, unweightedCor, weighted = "unweighted")
names(unwtCorr)
colnames(unwtCorr)[2] <- "corr"

wtCorr <- data.table(name.x, weightedCor, weighted="weighted")
names(wtCorr)
colnames(wtCorr)[2] <- "corr"

balance <- rbind(unwtCorr,wtCorr)
balance <- merge(balance, describe_x, by.x = "name.x", by.y = "var")
head(balance)

fwrite(balance, paste0(dir_balancePlots, "balance_bc_all_mac.csv")) # change

## 3.021 no2_all_mac ----
x <- balanceALL %>% select(-no2_all, # change
                           -bc_all.wt, -bc_all_mac.wt, 
                           -no2_all.wt, -no2_all_mac.wt)
name.x <- names(x)
unweightedCor <- matrix(NA,dim(x)[2],1)
weightedCor <- matrix(NA,dim(x)[2],1)

for (j in 1:dim(x)[2]){
  unweightedCor[j,1] = stats::cor(balanceALL$no2_all, x[[j]], method = "pearson") # change
}
for (j in 1:dim(x)[2]){
  weightedCor[j,1] = weightedCorr(balanceALL$no2_all, x[[j]],method = "pearson", # change
                                  weights = balanceALL$no2_all_mac.wt) # change
}

unwtCorr <- data.table(name.x, unweightedCor, weighted = "unweighted")
names(unwtCorr)
colnames(unwtCorr)[2] <- "corr"

wtCorr <- data.table(name.x, weightedCor, weighted="weighted")
names(wtCorr)
colnames(wtCorr)[2] <- "corr"

balance <- rbind(unwtCorr,wtCorr)
balance <- merge(balance, describe_x, by.x = "name.x", by.y = "var")
head(balance)

fwrite(balance, paste0(dir_balancePlots, "balance_no2_all_mac.csv")) # change

############################# 3.1 bc_30d ######################################
x <- balanceALL %>% select(-bc_30d, # change
                           -bc_30d.wt, -bc_3090d.wt, -bc_90280d.wt)
name.x <- names(x)
unweightedCor <- matrix(NA,dim(x)[2],1)
weightedCor <- matrix(NA,dim(x)[2],1)

for (j in 1:dim(x)[2]){
  unweightedCor[j,1] = stats::cor(balanceALL$bc_30d, x[[j]], method = "pearson") # change
}
for (j in 1:dim(x)[2]){
  weightedCor[j,1] = weightedCorr(balanceALL$bc_30d, x[[j]],method = "pearson", # change
                                weights = balanceALL$bc_30d.wt) # change
}

unwtCorr <- data.table(name.x, unweightedCor, weighted = "unweighted")
names(unwtCorr)
colnames(unwtCorr)[2] <- "corr"

wtCorr <- data.table(name.x, weightedCor, weighted="weighted")
names(wtCorr)
colnames(wtCorr)[2] <- "corr"

balance <- rbind(unwtCorr,wtCorr)
balance <- merge(balance, describe_x, by.x = "name.x", by.y = "var")
head(balance)

fwrite(balance, paste0(dir_balancePlots, "balance_bc_30d.csv")) # change

############################# 3.2 bc_3090d #####################################
x <- balanceALL %>% select(-bc_3090d, # change
                           -bc_30d.wt, -bc_3090d.wt, -bc_90280d.wt)
name.x <- names(x)
unweightedCor <- matrix(NA,dim(x)[2],1)
weightedCor <- matrix(NA,dim(x)[2],1)

for (j in 1:dim(x)[2]){
  unweightedCor[j,1] = stats::cor(balanceALL$bc_3090d, x[[j]],method = "pearson") # change
}
for (j in 1:dim(x)[2]){
  weightedCor[j,1] = weightedCorr(balanceALL$bc_3090d, x[[j]],method = "pearson", # change
                                  weights = balanceALL$bc_3090d.wt) # change
}

unwtCorr <- data.table(name.x, unweightedCor, weighted = "unweighted")
names(unwtCorr)
colnames(unwtCorr)[2] <- "corr"

wtCorr <- data.table(name.x, weightedCor, weighted="weighted")
names(wtCorr)
colnames(wtCorr)[2] <- "corr"

balance <- rbind(unwtCorr,wtCorr)
balance <- merge(balance, describe_x, by.x = "name.x", by.y = "var")
head(balance)

fwrite(balance, paste0(dir_balancePlots, "balance_bc_3090d.csv")) # change

############################# 3.3 bc_90280d ###################################
x <- balanceALL %>% select(-bc_90280d, # change
                           -bc_30d.wt, -bc_3090d.wt, -bc_90280d.wt)
name.x <- names(x)
unweightedCor <- matrix(NA,dim(x)[2],1)
weightedCor <- matrix(NA,dim(x)[2],1)

for (j in 1:dim(x)[2]){
  unweightedCor[j,1] = stats::cor(balanceALL$bc_90280d, x[[j]],method = "pearson") # change
}
for (j in 1:dim(x)[2]){
  weightedCor[j,1] = weightedCorr(balanceALL$bc_90280d, x[[j]],method = "pearson", # change
                                  weights = balanceALL$bc_90280d.wt) # change
}

unwtCorr <- data.table(name.x, unweightedCor, weighted = "unweighted")
names(unwtCorr)
colnames(unwtCorr)[2] <- "corr"

wtCorr <- data.table(name.x, weightedCor, weighted="weighted")
names(wtCorr)
colnames(wtCorr)[2] <- "corr"

balance <- rbind(unwtCorr,wtCorr)
balance <- merge(balance, describe_x, by.x = "name.x", by.y = "var")
head(balance)

fwrite(balance, paste0(dir_balancePlots, "balance_bc_90280d.csv")) # change

############################# 3.4 no2_30d ######################################
x <- balanceALL %>% select(-no2_30d, # change
                           -bc_30d.wt, -bc_3090d.wt, -bc_90280d.wt,
                           -no2_30d.wt, -no2_3090d.wt, -no2_90280d.wt)
name.x <- names(x)
unweightedCor <- matrix(NA,dim(x)[2],1)
weightedCor <- matrix(NA,dim(x)[2],1)

for (j in 1:dim(x)[2]){
  unweightedCor[j,1] = stats::cor(balanceALL$no2_30d, x[[j]],method = "pearson") # change
}
for (j in 1:dim(x)[2]){
  weightedCor[j,1] = weightedCorr(balanceALL$no2_30d, x[[j]],method = "pearson", # change
                                  weights = balanceALL$no2_30d.wt) # change
}

unwtCorr <- data.table(name.x, unweightedCor, weighted = "unweighted")
names(unwtCorr)
colnames(unwtCorr)[2] <- "corr"

wtCorr <- data.table(name.x, weightedCor, weighted="weighted")
names(wtCorr)
colnames(wtCorr)[2] <- "corr"

balance <- rbind(unwtCorr,wtCorr)
balance <- merge(balance, describe_x, by.x = "name.x", by.y = "var")
head(balance)

fwrite(balance, paste0(dir_balancePlots, "balance_no2_30d.csv")) # change

############################# 3.5 no2_3090d #####################################
x <- balanceALL %>% select(-no2_3090d, # change
                           -bc_30d.wt, -bc_3090d.wt, -bc_90280d.wt,
                           -no2_30d.wt, -no2_3090d.wt, -no2_90280d.wt)
name.x <- names(x)
unweightedCor <- matrix(NA,dim(x)[2],1)
weightedCor <- matrix(NA,dim(x)[2],1)

for (j in 1:dim(x)[2]){
  unweightedCor[j,1] = stats::cor(balanceALL$no2_3090d, x[[j]],method = "pearson") # change
}
for (j in 1:dim(x)[2]){
  weightedCor[j,1] = weightedCorr(balanceALL$no2_3090d, x[[j]],method = "pearson", # change
                                  weights = balanceALL$no2_3090d.wt) # change
}

unwtCorr <- data.table(name.x, unweightedCor, weighted = "unweighted")
names(unwtCorr)
colnames(unwtCorr)[2] <- "corr"

wtCorr <- data.table(name.x, weightedCor, weighted="weighted")
names(wtCorr)
colnames(wtCorr)[2] <- "corr"

balance <- rbind(unwtCorr,wtCorr)
balance <- merge(balance, describe_x, by.x = "name.x", by.y = "var")
head(balance)

fwrite(balance, paste0(dir_balancePlots, "balance_no2_3090d.csv")) # change

############################# 3.6 no2_90280d ###################################
x <- balanceALL %>% select(-no2_90280d, # change
                           -bc_30d.wt, -bc_3090d.wt, -bc_90280d.wt,
                           -no2_30d.wt, -no2_3090d.wt, -no2_90280d.wt)
name.x <- names(x)
unweightedCor <- matrix(NA,dim(x)[2],1)
weightedCor <- matrix(NA,dim(x)[2],1)

for (j in 1:dim(x)[2]){
  unweightedCor[j,1] = stats::cor(balanceALL$no2_90280d, x[[j]],method = "pearson") # change
}
for (j in 1:dim(x)[2]){
  weightedCor[j,1] = weightedCorr(balanceALL$no2_90280d, x[[j]],method = "pearson", # change
                                  weights = balanceALL$no2_90280d.wt) # change
}

unwtCorr <- data.table(name.x, unweightedCor, weighted = "unweighted")
names(unwtCorr)
colnames(unwtCorr)[2] <- "corr"

wtCorr <- data.table(name.x, weightedCor, weighted="weighted")
names(wtCorr)
colnames(wtCorr)[2] <- "corr"

balance <- rbind(unwtCorr,wtCorr)
balance <- merge(balance, describe_x, by.x = "name.x", by.y = "var")
head(balance)

fwrite(balance, paste0(dir_balancePlots, "balance_no2_90280d.csv")) # change


############################# 4. plot balance #################################
rm(list = ls())
dir_data <- "/media/gate/Shuxin/MAbirth/"
dir_balancePlots <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/3balancePlots/"
gc()

####################### codes for complete plot ###############################
# ggplot(balance, aes(x = description, y = corr)) +
#   geom_point(aes(colour = weighted), size = 1) +
#   geom_hline(yintercept=0, size=0.2) +
#   geom_hline(yintercept=-0.1, size=0.1, linetype = "dashed") +
#   geom_hline(yintercept=0.1, size=0.1, linetype = "dashed") +
#   theme(plot.title = element_blank(),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank()) +
#   coord_flip() +
#   theme(legend.position = "top", legend.title = element_blank())

## bc_all -----
balance <- fread(paste0(dir_balancePlots, "balance_bc_all.csv"))
# balance_display <- balance[ !name.x %in% c("no2_30d", "no2_3090d"),] # change?
balance <- balance[-c(47:50),]
balance$description[1:4] <- c(rep("Smoking during pregnancy",2), rep("Smoking before pregnancy",2))
balance$description[9:12] <- c(rep("Median household income",2), rep("Median value of house",2))
balance$description[25:26] <- rep("The other exposure",2)
balance$description[27:28] <- rep("% poverty",2)
balance$description[43:44] <- rep("Previous SGA infant")
balance$description[41:42] <- rep("Previous infant over 4 kg")
bc_balance <- balance

bc_plot <- ggplot(balance, aes(x = name.x, y = corr, group = weighted)) +
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
bc_plot

balance <- fread(paste0(dir_balancePlots, "balance_no2_all.csv"))
# balance_display <- balance[ !name.x %in% c("no2_30d", "no2_3090d"),] # change?
View(balance)
balance <- balance[-c(47:50),]
balance$description[3:6] <- c(rep("Smoking during pregnancy",2), rep("Smoking before pregnancy",2))
balance$description[11:14] <- c(rep("Median household income",2), rep("Median value of house",2))
balance$description[1:2] <- rep("The other exposure",2)
balance$description[27:28] <- rep("% poverty",2)
balance$description[43:44] <- rep("Previous SGA infant")
balance$description[41:42] <- rep("Previous infant over 4 kg")
no2_balance <- balance

no2_plot <- ggplot(no2_balance, aes(x = description, y = corr, group = weighted)) +
  geom_point(aes(colour = weighted), size = 2) +
  coord_flip() +
  geom_hline(yintercept=0, size=0.2) +
  geom_hline(yintercept=-0.1, size=0.2, linetype = "dashed") +
  geom_hline(yintercept=0.1, size=0.2, linetype = "dashed") +
  geom_hline(yintercept=-0.2, size=0.1, linetype = "dashed") +
  geom_hline(yintercept=0.2, size=0.1, linetype = "dashed") +
  theme(plot.title = element_blank(),
        # axis.text.y = element_blank(),
        axis.title.x = e