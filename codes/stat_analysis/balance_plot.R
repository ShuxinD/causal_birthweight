###############################################################################
# Project: Causal black carbon on birth weight in MA                          #
# Code: balance checking                                                      #
# Input: clean birth data "birth_final.csv", ipweights "birth_ipw.csv"        #
# Output: balance plot                                                        #
# Author: Shuxin Dong                                                         #
# Date: Nov 16, 2020                                                          #
###############################################################################

############################# 0. Setup ########################################
rm(list = ls())
gc()

library(dplyr)
library(data.table)
library(wCorr)
library(ggplot2)

dir_input <- "/Users/shuxind/Desktop/BC_birthweight_data/"
dir_output <- "/Users/shuxind/Documents/GitHub/causal_BC_birthweight/04_balance_checking/"
# setwd("/media/gate/Shuxin")
# dir_input <- "/media/gate/Shuxin/"
# dir_output <- "/media/gate/Shuxin/"

############################# 1. data manipulation ############################
# load data
birth <- fread(paste0(dir_input, "birth_final.csv"),
               drop = "V1")
birth$year <- as.factor(birth$year)
birth$m_edu <- as.factor(birth$m_edu)
birth$kotck <- as.factor(birth$kotck)
birth$m_wg_cat <- as.factor(birth$m_wg_cat)
var <- c("uniqueid_yr",
         "year","sex","married","mage","m_edu", "cigdpp","cigddp",
         "clinega","kotck","pncgov", "bwg", "rf_db_gest","rf_db_other",
         "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
         "rf_prev_sga", "bc_30d","bc_3090d", "bc_90280d", "firstborn","m_wg_cat",
         "log_mhincome", "log_mhvalue", "percentPoverty",
         "mrace_1", "mrace_2", "mrace_3", "mrace_4")
birth <- birth[ , var, with = F]

ipw <- fread(paste0(dir_input, "birth_ipw.csv"),
             select = c("birth.uniqueid_yr","bc_30d.wt.t", "bc_3090d.wt.t", "bc_90280d.wt.t"))

# merge IPW into birth data
birth_all <- left_join(birth, ipw, by=c("uniqueid_yr" = "birth.uniqueid_yr"))
write.csv(birth_all, file = paste0(dir_input, "birth_all.csv"))
rm(birth)
rm(ipw)
gc()

############################# 2. check balance ###############################
## load data
birth_all <- fread(paste0(dir_input, "birth_all.csv"),
                   drop = "V1")
## restrict the birth_all data to 
## the 0/1 and continuous variables; 
## birth weight
## IPWs
var <- c("sex","married","mage","cigdpp","cigddp",
         "clinega","pncgov", "rf_db_gest","rf_db_other",
         "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
         "rf_prev_sga", "firstborn",
         "log_mhincome", "log_mhvalue", "percentPoverty",
         "mrace_1", "mrace_2", "mrace_3", "mrace_4",
         "bc_30d","bc_3090d", "bc_90280d",
         "bc_30d.wt.t", "bc_3090d.wt.t", "bc_90280d.wt.t")
balanceALL <- birth_all[ , var, with = F]

description <- c("Newborn sex = female","Married","Mother's age",
           "# of daily cigarettes smoking before pregnancy",
           "# of daily cigarettes smoking during pregnancy",
           "Clinical gestational age (weeks)",
           "Has government support for prenatal care",
           "Gestational diabetes","Maternal diabetes",
           "Maternal chronic high blood pressure", 
           "Maternal high blood pressure during pregnancy",
           "Maternal incompetent cervix",
           "Mother with previous infant over 4000 g",
           "Mother with previous small-for-gestational-age infant", 
           "First born child",
           "Log. median household income (census-tract level)", 
           "Log. median value of house (census-tract level)", 
           "% Poverty (census-tract level)",
           "Maternal race = white", "Maternal race = black", 
           "Maternal race = Asian/Pacific Islander", 
           "Maternal race = others",
           "Average of daily BC exposure over 0-30 days prior to the delivery date",
           "Average of daily BC exposure over 31-90 days prior to the delivery date", 
           "Average of daily BC exposure over 91-280 days prior to the delivery date",
           "bc_30d.wt.t", "bc_3090d.wt.t", "bc_90280d.wt.t")
describe_x <- data.frame(description, var)

############################# 2.1 bc_30d #####################################
x <- balanceALL %>% select(-bc_30d, -bc_30d.wt.t, -bc_3090d.wt.t,
                                    -bc_90280d.wt.t)
name.x <- names(x)
# j.drop <- match(c("T"), names(data))
# j.drop <- j.drop[!is.na(j.drop)]
x <- as.data.frame(x)
unweightedCor <- matrix(NA,dim(x)[2],1)
weightedCor <- matrix(NA,dim(x)[2],1)
for (j in 1:dim(x)[2]){
  unweightedCor[j] = stats::cor(balanceALL$bc_30d, x[,j],method = "pearson")
}
for (j in 1:dim(x)[2]){
  weightedCor[j] = weightedCorr(balanceALL$bc_30d, x[,j],method = "pearson",
                                weights = balanceALL$bc_30d.wt.t)
}
unwtCorr <- data.frame(name.x, unweightedCor, weighted="unweighted")
colnames(unwtCorr)[2] <- "corr"
wtCorr <- data.frame(name.x, weightedCor, weighted="weighted")
colnames(wtCorr)[2] <- "corr"
balance_30d <- rbind(unwtCorr,wtCorr)
balance_30d <- left_join(balance_30d, describe_x, by = c("name.x" = "var"))
head(balance_30d)

############################# 2.2 bc_3090d #####################################
x <- balanceALL %>% select(-bc_3090d, -bc_30d.wt.t, -bc_3090d.wt.t,
                           -bc_90280d.wt.t)
name.x <- names(x)
# j.drop <- match(c("T"), names(data))
# j.drop <- j.drop[!is.na(j.drop)]
x <- as.data.frame(x)
unweightedCor <- matrix(NA,dim(x)[2],1)
weightedCor <- matrix(NA,dim(x)[2],1)
for (j in 1:dim(x)[2]){
  unweightedCor[j] = stats::cor(balanceALL$bc_3090d, x[,j],method = "pearson")
}
for (j in 1:dim(x)[2]){
  weightedCor[j] = weightedCorr(balanceALL$bc_3090d, x[,j],method = "pearson",
                                weights = balanceALL$bc_3090d.wt.t)
}
unwtCorr <- data.frame(name.x, unweightedCor, weighted="unweighted")
colnames(unwtCorr)[2] <- "corr"
wtCorr <- data.frame(name.x, weightedCor, weighted="weighted")
colnames(wtCorr)[2] <- "corr"
balance_3090d <- rbind(unwtCorr,wtCorr)
balance_3090d <- left_join(balance_3090d, describe_x, by = c("name.x" = "var"))
head(balance_3090d)

############################# 2.3 bc_90280d ###################################
x <- balanceALL %>% select(-bc_90280d, -bc_30d.wt.t, -bc_3090d.wt.t,
                           -bc_90280d.wt.t)
name.x <- names(x)
# j.drop <- match(c("T"), names(data))
# j.drop <- j.drop[!is.na(j.drop)]
x <- as.data.frame(x)
unweightedCor <- matrix(NA,dim(x)[2],1)
weightedCor <- matrix(NA,dim(x)[2],1)
for (j in 1:dim(x)[2]){
  unweightedCor[j] = stats::cor(balanceALL$bc_90280d, x[,j],method = "pearson")
}
for (j in 1:dim(x)[2]){
  weightedCor[j] = weightedCorr(balanceALL$bc_90280d, x[,j],method = "pearson",
                                weights = balanceALL$bc_90280d.wt.t)
}
unwtCorr <- data.frame(name.x, unweightedCor, weighted="unweighted")
colnames(unwtCorr)[2] <- "corr"
wtCorr <- data.frame(name.x, weightedCor, weighted="weighted")
colnames(wtCorr)[2] <- "corr"
balance_90280d <- rbind(unwtCorr,wtCorr)
balance_90280d <- left_join(balance_90280d, describe_x, by = c("name.x" = "var"))
head(balance_90280d)

############################# 3. plot balance #################################
# rm(balanceALL)
# rm(birth_all)
gc()

library(gridExtra)

# balance_30d <- read.csv(file = paste0(dir_output, "balance_30d.csv"),
#                         row.names = 1)
balance_3090d <- read.csv(file = paste0(dir_output, "balance_3090d.csv"),
                        row.names = 1)
balance_3090d <- balance_3090d %>% filter(!name.x %in% c("bc_30d", "bc_90280d"))
balance_3090d$description[balance_3090d$description=="Newborn sex = female"] <-
  "New-born sex = Female"
balance_3090d$description[balance_3090d$description=="Has government support for prenatal care"] <-
  "Government support for prenatal care"
balance_3090d$description[balance_3090d$description=="First born child"] <-
  "First-born child"
balance_3090d$description[balance_3090d$description=="Clinical gestational age (weeks)"] <-
  "Clinical gestational age"
# balance_90280d <- read.csv(file = paste0(dir_output, "balance_90280d.csv"),
#                         row.names = 1)

# pdf(file = paste0(dir_output,"balance_30d.pdf"))
# bc30d <- 
#   ggplot(balance_30d, aes(x = description, y = corr)) +
#   geom_point(aes(colour = weighted), size = 1) +
#   geom_hline(yintercept=0, size=0.2) +
#   geom_hline(yintercept=-0.1, size=0.1, linetype = "dashed") +
#   geom_hline(yintercept=0.1, size=0.1, linetype = "dashed") +
#   theme(plot.title = element_blank(),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank()) +
#   coord_flip() +
#   theme(legend.position = "top", legend.title = element_blank())
# dev.off()

pdf(file = paste0(dir_output,"balance_3090d-new.pdf"))
# bc3090d <- 
ggplot(balance_3090d, aes(x = description, y = corr)) +
  geom_point(aes(colour = weighted), size = 1) +
  geom_hline(yintercept=0, size=0.2) +
  geom_hline(yintercept=-0.1, size=0.1, linetype = "dashed") +
  geom_hline(yintercept=0.1, size=0.1, linetype = "dashed") +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  coord_flip() +
  theme(legend.position = "bottom", legend.title = element_blank())
dev.off()

# pdf(file = paste0(dir_output,"balance_90280d.pdf"))
# bc90280d <-
# ggplot(balance_90280d, aes(x = description, y = corr)) +
#   geom_point(aes(colour = weighted), size = 1) +
#   geom_hline(yintercept=0, size=0.2) +
#   geom_hline(yintercept=-0.1, size=0.1, linetype = "dashed") +
#   geom_hline(yintercept=0.1, size=0.1, linetype = "dashed") +
#   theme(plot.title = element_blank(),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank()) +
#   coord_flip() 
# dev.off()

# three <- plot_grid(bc30d + theme(legend.position = "none"), 
#           bc3090d + theme(legend.position = "none"), 
#           bc90280d + theme(legend.position = "none"), 
#           nrow = 3, ncol = 1, labels = "AUTO")

# legend <- get_legend(bc30d)
# pdf(file = "/Users/shuxind/Documents/GitHub/causal_BC_birthweight/results/balancePlot.pdf",
#     width = 10, height = 15)
# grid.arrange(three, legend, ncol=1, nrow=2,
#              widths = 10, heights = c(15, 0.2))
# dev.off()
