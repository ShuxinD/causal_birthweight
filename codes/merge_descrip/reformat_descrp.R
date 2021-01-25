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
dir_output_table1 <- "/media/gate/Shuxin/MAbirth/results/"

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

par(mfrow=c(2,3))
plot(year, main="year")
plot(as.factor(mrace), main="mrace")
plot(m_edu, main="m_edu")
plot(kotck, main="kotck")
plot(m_wg_cat, main="maternal weight change")

par(mfrow = c(1,1))
hist(bwg)

par(mfrow=c(3,2))
hist(bc_30d)
hist(bc_3090d)
hist(bc_90280d)
hist(no2_30d)
hist(no2_3090d)
hist(no2_90280d)

par(mfrow=c(1,2))
plot(density(mage,  bw=2), main = "mage")
hist(clinega)

# cigdpp %>% sort(decreasing = TRUE) %>% head(30)
par(mfrow=c(2,2))
plot(density(cigdpp,  bw=0.5))
hist(log(cigdpp))
plot(density(log(cigdpp),  bw=0.5))
# max_cigdpp <- quantile(cigdpp, 0.99995)
summary(cigdpp)

par(mfrow=c(2,2))
plot(density(cigddp,  bw=0.5))
hist(log(cigddp))
plot(density(log(cigddp),  bw=0.5))
# max_cigddp <- quantile(cigddp, 0.99995)
summary(cigddp)

# med_hs_inc %>% sort(decreasing = TRUE) %>% head(30)
par(mfrow=c(2,2))
plot(density(mhincome, bw=1))
hist(log(mhincome))
plot(density(log(mhincome),  bw=0.5))
birth[, log_mhincome:= log(mhincome)][]
summary(mhincome)

par(mfrow=c(2,2))
plot(density(mhvalue, bw=1))
hist(log(mhvalue))
plot(density(log(mhvalue),  bw=0.5))
birth[, log_mhvalue:= log(mhvalue)][]
summary(mhvalue)

detach(birth)

## truncate the extreme values for smoking
max_cigdpp <- quantile(cigdpp, 0.99995)
max_cigddp <- quantile(cigddp, 0.99995)
birth$cigdpp <- fifelse(cigdpp > max_cigdpp, max_cigdpp, cigdpp)
birth$cigddp <- fifelse(cigddp > max_cigddp, max_cigddp, cigddp)
summary(birth$cigdpp)
summary(birth$cigddp)

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

rtffile <- RTF(file = paste0(dir_output_table1, "table1.doc"))  # this can be an .rtf or a .doc
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
write.csv(smokerInfo, paste0(dir_output_table1, "table1_smokerInfo.csv"))

fwrite(birth, paste0(dir_input_birth, "MAbirth_for_analyses.csv"))
