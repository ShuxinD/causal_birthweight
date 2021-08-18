###############################################################################
# Project: Causal black carbon and NO2 on birth weight in MA                  #
# Code: Table one and normality                                               #
# Input: "MAbirth_merged.csv"                                                 #
# Output: "MAbirth_for_analyses.csv" for future analyses                      #
# Author: Shuxin Dong                                                         #
# Date: 2021-01-21.                                                           #
###############################################################################

############################# 0. Setup ########################################
rm(list = ls())
gc()

library(data.table)
library(fastDummies)
library(tableone)
library(rtf)

setwd("/media/gate/Shuxin/")
dir_input_birth <- "/media/gate/Shuxin/MAbirth/"
dir_descrp <- "/media/gate/Shuxin/MAbirth/results/0TableOne_checkVarDistribution/"

birth <- fread(paste0(dir_input_birth, "MAbirth_merged.csv"))
birth[, lbw:=0][]
birth$lbw[birth$bwg<2500] <- 1
# > names(birth)
# [1] "uniqueid_yr"    "year"           "sex"            "married"        "mage"          
# [6] "mrace"          "m_edu"          "cigdpp"         "cigddp"         "parit"         
# [11] "clinega"        "kotck"          "pncgov"         "bwg"            "rf_db_gest"    
# [16] "rf_db_other"    "rf_hbp_chronic" "rf_hbp_pregn"   "rf_cervix"      "rf_prev_4kg"   
# [21] "rf_prev_sga"    "m_wg"           "mhincome"       "mhvalue"        "percentPoverty"
# [26] "bc_30d"         "bc_3090d"       "bc_90280d"      "no2_30d"        "no2_3090d"     
# [31] "no2_90280d"     "lbw"

######################## 1. Prepare for table one #############################
## Transform continuous to categorical
## `parity` to `firstborn`, because of some extreme large unexplained value.
birth[, firstborn := parit][]
birth$firstborn[birth$parit>1] <- 0
birth[, parit := NULL]
## `m_wg` to `m_wg_cat` based on reference. Extreme large unexplained values exist.
plot(density(birth$m_wg, bw=1))
hist(birth$m_wg)
birth[, m_wg_cat := 0][]
birth$m_wg_cat[birth$m_wg<0] <- 1
birth$m_wg_cat[birth$m_wg>=0 & birth$m_wg<15] <- 2
birth$m_wg_cat[birth$m_wg>=15 & birth$m_wg<25] <- 3
birth$m_wg_cat[birth$m_wg>=25 & birth$m_wg<36] <- 4
birth$m_wg_cat[birth$m_wg>=36] <- 5
birth[, m_wg :=NULL]
## create non-smoker variable as binary variable
birth$smoker_ddp <- 0
birth$smoker_ddp[birth$cigddp!=0] <- 1
birth$smoker_dpp <- 0
birth$smoker_dpp[birth$cigdpp!=0] <- 1
summary(as.factor(birth$smoker_ddp))
summary(as.factor(birth$smoker_dpp))

## Get categorical and dummy variables
birth$mrace[birth$mrace==5] <- 4
birth$year <- as.factor(birth$year)
birth$sex <- birth$sex - 1 # 1 is girl; 0 is boy
birth$m_edu <- as.factor(birth$m_edu)
birth$kotck <- as.factor(birth$kotck)
birth$m_wg_cat <- as.factor(birth$m_wg_cat)
birth <- fastDummies::dummy_cols(birth, select_columns = "mrace")

######################## 2. Check variables ###################################
## check categorical variables with multiple levels
## check normality on continuous variables: take log for those skewed variables
## check extreme values
attach(birth)

pdf(paste0(dir_descrp,"check1.pdf"))
par(mfrow=c(2,3))
plot(year, main="year")
plot(as.factor(mrace), main="mrace")
plot(m_edu, main="m_edu")
plot(kotck, main="kotck")
plot(m_wg_cat, main="maternal weight change")
dev.off()

pdf(paste0(dir_descrp,"check2.pdf"))
par(mfrow = c(1,1))
hist(bwg)
dev.off()

pdf(paste0(dir_descrp,"check3.pdf"))
par(mfrow=c(3,2))
hist(bc_30d)
hist(bc_3090d)
hist(bc_90280d)
hist(no2_30d)
hist(no2_3090d)
hist(no2_90280d)
dev.off()

pdf(paste0(dir_descrp,"check4.pdf"))
par(mfrow=c(1,2))
plot(density(mage,  bw=2), main = "mage")
hist(clinega)
dev.off()

# cigdpp %>% sort(decreasing = TRUE) %>% head(30)
pdf(paste0(dir_descrp,"check5.pdf"))
par(mfrow=c(2,2))
plot(density(cigdpp,  bw=0.5))
hist(log(cigdpp))
plot(density(log(cigdpp),  bw=0.5))
dev.off()
# max_cigdpp <- quantile(cigdpp, 0.99995)
summary(cigdpp)

pdf(paste0(dir_descrp,"check6.pdf"))
par(mfrow=c(2,2))
plot(density(cigddp,  bw=0.5))
hist(log(cigddp))
plot(density(log(cigddp),  bw=0.5))
dev.off()
# max_cigddp <- quantile(cigddp, 0.99995)
summary(cigddp)

# med_hs_inc %>% sort(decreasing = TRUE) %>% head(30)
pdf(paste0(dir_descrp,"check7.pdf"))
par(mfrow=c(2,2))
plot(density(mhincome, bw=1))
hist(log(mhincome))
plot(density(log(mhincome),  bw=0.5))
birth[, log_mhincome:= log(mhincome)][]
dev.off()
summary(mhincome)

pdf(paste0(dir_descrp,"check8.pdf"))
par(mfrow=c(2,2))
plot(density(mhvalue, bw=1))
hist(log(mhvalue))
plot(density(log(mhvalue),  bw=0.5))
dev.off()
birth[, log_mhvalue:= log(mhvalue)][]
summary(mhvalue)

## truncate the extreme values for smoking
max_cigdpp <- quantile(cigdpp, 0.99995)
max_cigddp <- quantile(cigddp, 0.99995)
birth$cigdpp <- fifelse(cigdpp > max_cigdpp, max_cigdpp, cigdpp)
birth$cigddp <- fifelse(cigddp > max_cigddp, max_cigddp, cigddp)
summary(birth$cigdpp)
summary(birth$cigddp)

detach(birth)

names(birth)
# > names(birth)
# [1] "uniqueid_yr"    "year"           "sex"            "married"        "mage"          
# [6] "mrace"          "m_edu"          "cigdpp"         "cigddp"         "clinega"       
# [11] "kotck"          "pncgov"         "bwg"            "rf_db_gest"     "rf_db_other"   
# [16] "rf_hbp_chronic" "rf_hbp_pregn"   "rf_cervix"      "rf_prev_4kg"    "rf_prev_sga"   
# [21] "mhincome"       "mhvalue"        "percentPoverty" "bc_30d"         "bc_3090d"      
# [26] "bc_90280d"      "no2_30d"        "no2_3090d"      "no2_90280d"     "lbw"           
# [31] "firstborn"      "m_wg_cat"       "smoker_ddp"     "smoker_dpp"     "mrace_1"       
# [36] "mrace_2"        "mrace_3"        "mrace_4"        "log_mhincome"   "log_mhvalue" 
summary(birth)
######################## 3. Export table 1 ####################################
listVars <- c("bwg",
              "bc_30d","bc_3090d", "bc_90280d",
              "no2_30d", "no2_3090d", "no2_90280d",
              "year",
              "mage","sex","married", "mrace","m_edu", "pncgov",
              "mhincome", "mhvalue", "percentPoverty",
              "smoker_dpp", "smoker_ddp", "cigdpp","cigddp",
              "clinega", "firstborn", "kotck", "m_wg_cat",
              "rf_db_gest", "rf_db_other", "rf_hbp_pregn", "rf_hbp_chronic", 
              "rf_cervix", "rf_prev_4kg", "rf_prev_sga")
catVars <- c("year","sex","married","mrace","m_edu","pncgov",
             "smoker_dpp", "smoker_ddp",
             "firstborn", "kotck", "m_wg_cat", 
             "rf_db_gest", "rf_db_other", "rf_hbp_pregn", "rf_hbp_chronic", 
             "rf_cervix", "rf_prev_4kg", "rf_prev_sga")
bcVars <- c("bc_30d", "bc_3090d", "bc_90280d")
no2Vars <- c("no2_30d", "no2_3090d", "no2_90280d")

rawtable1 <- CreateTableOne(vars = listVars, factorVars = catVars, 
                            strata = "lbw", addOverall = TRUE, test = FALSE,
                            data = birth)
print(rawtable1, nonnormal = c(bcVars, no2Vars),
      formatOption = list(decimal.mark = ".",  big.mark = ",", scientific = FALSE),
      contDigits = 2)
table1 <- print(rawtable1, nonnormal = c(bcVars, no2Vars),
                formatOption = list(decimal.mark = ".",  big.mark = ",", scientific = FALSE),
                contDigits = 2)

rtffile <- RTF(file = paste0(dir_descrp, "table1.doc"))  # this can be an .rtf or a .doc
addParagraph(rtffile, "Table")
addTable(rtffile, cbind(rownames(table1), table1))
done(rtffile)

## revise the output of smoking variable
cigddp_smoker <- birth$cigddp[birth$cigddp!=0] # during pregnancy
cigdpp_smoker <- birth$cigdpp[birth$cigdpp!=0] # prior to pregnancy
cigddp_smoker_nbw <- birth$cigddp[birth$cigddp!=0 & birth$lbw==0]
cigdpp_smoker_nbw <- birth$cigdpp[birth$cigdpp!=0 & birth$lbw==0]
cigddp_smoker_lbw <- birth$cigddp[birth$cigddp!=0 & birth$lbw==1]
cigdpp_smoker_lbw <- birth$cigdpp[birth$cigdpp!=0 & birth$lbw==1]

num_smoker <- c(length(cigddp_smoker), length(cigdpp_smoker),
                length(cigddp_smoker_nbw), length(cigdpp_smoker_nbw),
                length(cigddp_smoker_lbw), length(cigdpp_smoker_lbw))
mean_cig <- c(mean(cigddp_smoker), mean(cigdpp_smoker),
              mean(cigddp_smoker_nbw), mean(cigdpp_smoker_nbw),
              mean(cigddp_smoker_lbw), mean(cigdpp_smoker_lbw))
sd_cig <- c(sd(cigddp_smoker), sd(cigdpp_smoker),
            sd(cigddp_smoker_nbw), sd(cigdpp_smoker_nbw),
            sd(cigddp_smoker_lbw), sd(cigdpp_smoker_lbw))
smokerInfo <- cbind(num_smoker,mean_cig,sd_cig)
rownames(smokerInfo) <- c("during pregnancy", "before pregnancy",
                          "during pregnancy with NBW", "before preganncy with NBW",
                          "during pregnancy with LBW", "before preganncy with LBW")
smokerInfo
# > smokerInfo
# num_smoker  mean_cig   sd_cig
# during pregnancy               61703  7.587373 5.525075
# before pregnancy              113588 12.193271 8.190175
# during pregnancy with NBW      58565  7.539480 5.501169
# before preganncy with NBW     109511 12.133674 8.174378
# during pregnancy with LBW       3138  8.481198 5.883473
# before preganncy with LBW       4077 13.794076 8.448792
write.csv(smokerInfo, paste0(dir_descrp, "table1_smokerInfo.csv"))

mean_bc <- rbind(c(mean(birth[,bc_30d]), mean(birth[lbw==0, bc_30d]),mean(birth[lbw==1, bc_30d])),
             c(mean(birth[,bc_3090d]), mean(birth[lbw==0, bc_3090d]),mean(birth[lbw==1, bc_3090d])),
             c(mean(birth[,bc_90280d]), mean(birth[lbw==0, bc_90280d]),mean(birth[lbw==1, bc_90280d])))
rownames(mean_bc) <- c("bc_30d", "bc_3090d", "bc_90280d")
colnames(mean_bc) <- c("overall", "nbw", "lbw")
mean_bc
# > mean_bc
# overall       nbw       lbw
# bc_30d    0.4770028 0.4767743 0.4875047
# bc_3090d  0.4775469 0.4773306 0.4874869
# bc_90280d 0.4798614 0.4796608 0.4890775

IQR_bc <- rbind(quantile(birth[,bc_30d]), quantile(birth[lbw==0, bc_30d]),quantile(birth[lbw==1, bc_30d]),
                 quantile(birth[,bc_3090d]), quantile(birth[lbw==0, bc_3090d]),quantile(birth[lbw==1, bc_3090d]),
                 quantile(birth[,bc_90280d]), quantile(birth[lbw==0, bc_90280d]),quantile(birth[lbw==1, bc_90280d]))
IQR_bc
# 0%       25%       50%       75%     100%
# [1,] 0.1630848 0.3807746 0.4480162 0.5445652 1.826903
# [2,] 0.1630848 0.3806807 0.4478211 0.5442752 1.826903
# [3,] 0.1866778 0.3853661 0.4580873 0.5599305 1.473415
# [4,] 0.1788135 0.3840173 0.4498046 0.5424744 1.776729
# [5,] 0.1788135 0.3839377 0.4496089 0.5421738 1.776729
# [6,] 0.1994498 0.3879015 0.4602769 0.5570112 1.354540
# [7,] 0.1968862 0.3917479 0.4540661 0.5386503 1.642583
# [8,] 0.1968862 0.3917036 0.4539207 0.5383319 1.642583
# [9,] 0.2221745 0.3937925 0.4616684 0.5540890 1.401050

mean_no2 <- rbind(c(mean(birth[,no2_30d]), mean(birth[lbw==0, no2_30d]),mean(birth[lbw==1, no2_30d])),
                 c(mean(birth[,no2_3090d]), mean(birth[lbw==0, no2_3090d]),mean(birth[lbw==1, no2_3090d])),
                 c(mean(birth[,no2_90280d]), mean(birth[lbw==0, no2_90280d]),mean(birth[lbw==1, no2_90280d])))
rownames(mean_no2) <- c("no2_30d", "no2_3090d", "no2_90280d")
colnames(mean_no2) <- c("overall", "nbw", "lbw")
mean_no2
# > mean_no2
# overall      nbw      lbw
# no2_30d    22.74911 22.73950 23.19087
# no2_3090d  22.93494 22.92671 23.31304
# no2_90280d 23.06401 23.05510 23.47348

IQR_no2 <- rbind(quantile(birth[,no2_30d]), quantile(birth[lbw==0, no2_30d]),quantile(birth[lbw==1, no2_30d]),
                quantile(birth[,no2_3090d]), quantile(birth[lbw==0, no2_3090d]),quantile(birth[lbw==1, no2_3090d]),
                quantile(birth[,no2_90280d]), quantile(birth[lbw==0, no2_90280d]),quantile(birth[lbw==1, no2_90280d]))
IQR_no2
# > IQR_no2
# 0%      25%      50%      75%     100%
# [1,] 0.8735659 16.09391 22.58374 28.88830 75.22628
# [2,] 0.8735659 16.08221 22.56944 28.88265 75.22628
# [3,] 1.4735506 16.68021 23.22411 29.16093 62.31848
# [4,] 1.0452284 16.53250 22.92165 28.97721 83.84932
# [5,] 1.0842349 16.52219 22.91188 28.96856 83.84932
# [6,] 1.0452284 16.99337 23.35559 29.29263 62.30313
# [7,] 1.4542685 17.66789 23.02038 28.41093 64.61834
# [8,] 1.4542685 17.65838 23.01136 28.40323 64.61834
# [9,] 3.4550755 18.14295 23.43051 28.75076 63.98316

summary(birth)
# > summary(birth)
# uniqueid_yr             year             sex            married           mage           mrace      
# Length:844554      2001   : 69464   Min.   :0.0000   Min.   :0.000   Min.   :12.00   Min.   :1.000  
# Class :character   2002   : 68659   1st Qu.:0.0000   1st Qu.:0.000   1st Qu.:26.00   1st Qu.:1.000  
# Mode  :character   2003   : 67866   Median :0.0000   Median :1.000   Median :30.58   Median :1.000  
# 2004   : 64118   Mean   :0.4901   Mean   :0.687   Mean   :30.10   Mean   :1.571  
# 2006   : 64025   3rd Qu.:1.0000   3rd Qu.:1.000   3rd Qu.:34.33   3rd Qu.:2.000  
# 2008   : 63703   Max.   :1.0000   Max.   :1.000   Max.   :58.08   Max.   :4.000  
# (Other):446719                                                                   
# m_edu          cigdpp           cigddp           clinega      kotck          pncgov            bwg      
# 1: 88944   Min.   :  0.00   Min.   : 0.0000   Min.   :37.00   1: 72056   Min.   :0.0000   Min.   : 505  
# 2:197554   1st Qu.:  0.00   1st Qu.: 0.0000   1st Qu.:39.00   2: 66124   1st Qu.:0.0000   1st Qu.:3120  
# 3:192470   Median :  0.00   Median : 0.0000   Median :39.00   3:412321   Median :0.0000   Median :3430  
# 4:220839   Mean   :  1.64   Mean   : 0.5543   Mean   :39.32   4:294053   Mean   :0.3324   Mean   :3441  
# 5:144747   3rd Qu.:  0.00   3rd Qu.: 0.0000   3rd Qu.:40.00              3rd Qu.:1.0000   3rd Qu.:3742  
# Max.   :155.45   Max.   :45.0000   Max.   :42.00              Max.   :1.0000   Max.   :5982  
# 
# rf_db_gest       rf_db_other       rf_hbp_chronic    rf_hbp_pregn       rf_cervix         rf_prev_4kg      
# Min.   :0.00000   Min.   :0.000000   Min.   :0.0000   Min.   :0.00000   Min.   :0.000000   Min.   :0.000000  
# 1st Qu.:0.00000   1st Qu.:0.000000   1st Qu.:0.0000   1st Qu.:0.00000   1st Qu.:0.000000   1st Qu.:0.000000  
# Median :0.00000   Median :0.000000   Median :0.0000   Median :0.00000   Median :0.000000   Median :0.000000  
# Mean   :0.04156   Mean   :0.008305   Mean   :0.0119   Mean   :0.03368   Mean   :0.003594   Mean   :0.006714  
# 3rd Qu.:0.00000   3rd Qu.:0.000000   3rd Qu.:0.0000   3rd Qu.:0.00000   3rd Qu.:0.000000   3rd Qu.:0.000000  
# Max.   :1.00000   Max.   :1.000000   Max.   :1.0000   Max.   :1.00000   Max.   :1.000000   Max.   :1.000000  
# 
# rf_prev_sga         mhincome         mhvalue        percentPoverty       bc_30d          bc_3090d     
# Min.   :0.00000   Min.   : 10706   Min.   :  51840   Min.   : 0.000   Min.   :0.1631   Min.   :0.1788  
# 1st Qu.:0.00000   1st Qu.: 43653   1st Qu.: 204480   1st Qu.: 4.125   1st Qu.:0.3808   1st Qu.:0.3840  
# Median :0.00000   Median : 59452   Median : 269600   Median : 7.599   Median :0.4480   Median :0.4498  
# Mean   :0.00783   Mean   : 61956   Mean   : 293825   Mean   :11.610   Mean   :0.4770   Mean   :0.4775  
# 3rd Qu.:0.00000   3rd Qu.: 76389   3rd Qu.: 349350   3rd Qu.:15.963   3rd Qu.:0.5446   3rd Qu.:0.5425  
# Max.   :1.00000   Max.   :217583   Max.   :1294600   Max.   :70.176   Max.   :1.8269   Max.   :1.7767  
# 
# bc_90280d         no2_30d          no2_3090d        no2_90280d          lbw           firstborn     
# Min.   :0.1969   Min.   : 0.8736   Min.   : 1.045   Min.   : 1.454   Min.   :0.0000   Min.   :0.0000  
# 1st Qu.:0.3917   1st Qu.:16.0939   1st Qu.:16.532   1st Qu.:17.668   1st Qu.:0.0000   1st Qu.:0.0000  
# Median :0.4541   Median :22.5837   Median :22.922   Median :23.020   Median :0.0000   Median :0.0000  
# Mean   :0.4799   Mean   :22.7491   Mean   :22.935   Mean   :23.064   Mean   :0.0213   Mean   :0.4526  
# 3rd Qu.:0.5387   3rd Qu.:28.8883   3rd Qu.:28.977   3rd Qu.:28.411   3rd Qu.:0.0000   3rd Qu.:1.0000  
# Max.   :1.6426   Max.   :75.2263   Max.   :83.849   Max.   :64.618   Max.   :1.0000   Max.   :1.0000  
# 
# m_wg_cat     smoker_ddp        smoker_dpp        mrace_1          mrace_2          mrace_3       
# 1:  2766   Min.   :0.00000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000   Min.   :0.00000  
# 2: 70462   1st Qu.:0.00000   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.00000  
# 3:185464   Median :0.00000   Median :0.0000   Median :1.0000   Median :0.0000   Median :0.00000  
# 4:359019   Mean   :0.07306   Mean   :0.1345   Mean   :0.7269   Mean   :0.0856   Mean   :0.07721  
# 5:226843   3rd Qu.:0.00000   3rd Qu.:0.0000   3rd Qu.:1.0000   3rd Qu.:0.0000   3rd Qu.:0.00000  
# Max.   :1.00000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.00000  
# 
# mrace_4        log_mhincome     log_mhvalue   
# Min.   :0.0000   Min.   : 9.279   Min.   :10.86  
# 1st Qu.:0.0000   1st Qu.:10.684   1st Qu.:12.23  
# Median :0.0000   Median :10.993   Median :12.50  
# Mean   :0.1103   Mean   :10.942   Mean   :12.50  
# 3rd Qu.:0.0000   3rd Qu.:11.244   3rd Qu.:12.76  
# Max.   :1.0000   Max.   :12.290   Max.   :14.07  

fwrite(birth, paste0(dir_input_birth, "MAbirth_for_analyses.csv"))
