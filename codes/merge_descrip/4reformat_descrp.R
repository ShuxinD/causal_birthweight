###############################################################################
# Project: Causal black carbon on birth weight in MA                          #
# Code: Table one and normality                                               #
# Input: "mabirths_02NOV18.csv" original birth data                           #
# Input: "birth_SES.csv" birth data with merged SES variables                 #
# Output: "birth_final.csv" for future analyses                               #
# Author: Shuxin Dong                                                         #
# Date: Oct 28, 2020                                                          #
###############################################################################

############################# 0. Setup ########################################
rm(list = ls())
gc()

library(dplyr)
library(data.table)
library(fastDummies)
library(tableone)
library(rtf)

dir_input_birth <- "/Users/shuxind/Desktop/BC_birthweight_analysis/"
dir_output_table1 <- "/Users/shuxind/Desktop/BC_birthweight_analysis/"

birth <- fread(paste0(dir_input_birth, "birth_complt.csv"), drop = "V1")
birth$lbw <- 0
birth$lbw[birth$bwg<2500] <- 1
names(birth)
# > names(birth)
# [1] "uniqueid_yr"    "year"           "sex"            "married"        "mage"          
# [6] "mrace"          "m_edu"          "cigdpp"         "cigddp"         "parit"         
# [11] "clinega"        "kotck"          "pncgov"         "bwg"            "rf_db_gest"    
# [16] "rf_db_other"    "rf_hbp_chronic" "rf_hbp_pregn"   "rf_cervix"      "rf_prev_4kg"   
# [21] "rf_prev_sga"    "m_wg"           "mhincome"       "mhvalue"        "percentPoverty"
# [26] "bc_30d"         "bc_90d"         "bc_280d"        "lbw"  

######################## 1. Prepare for table one #############################
## Transform continuous to categorical
## `parity` to `firstborn`, because of some extreme large unexplained value.
## `m_wg` to `m_wg_cat` based on reference. Extreme large unexplained values exist.
birth$firstborn <- birth$parit
birth$firstborn[birth$parit>1 & !is.na(birth$parit)] <- 0
plot(density(birth$m_wg, bw=1))
hist(birth$m_wg)
birth$m_wg_cat <- 0
birth$m_wg_cat[birth$m_wg<0] <- 1
birth$m_wg_cat[birth$m_wg>=0 & birth$m_wg<15] <- 2
birth$m_wg_cat[birth$m_wg>=15 & birth$m_wg<25] <- 3
birth$m_wg_cat[birth$m_wg>=25 & birth$m_wg<36] <- 4
birth$m_wg_cat[birth$m_wg>=36] <- 5
birth <- birth %>% select(-parit, -m_wg)
summary(birth)
## get categorical and dummy variables
data <- birth
data$mrace[data$mrace==5] <- 4
data$year <- as.factor(data$year)
data$sex <- data$sex - 1 # 1 is girl; 0 is boy
data$m_edu <- as.factor(data$m_edu)
data$kotck <- as.factor(data$kotck)
data$m_wg_cat <- as.factor(data$m_wg_cat)
data <- fastDummies::dummy_cols(data, select_columns = "mrace")
data$bc_3090d <- with(data, (bc_90d*90-bc_30d*30)/(90-30))
data$bc_90280d <- with(data, (bc_280d*280-bc_90d*90)/(280-90))

######################## 2. Export table 1 ####################################
listVars <- c("bwg","bc_30d",
              "bc_90d","bc_280d",
              "bc_3090d", "bc_90280d",
              "year",
              "mage","sex","married", "mrace","m_edu", "pncgov",
              "mhincome", "mhvalue", "percentPoverty",
              "cigdpp","cigddp","clinega", "firstborn", "kotck", "m_wg_cat",
              "rf_db_gest", "rf_db_other", "rf_hbp_pregn", "rf_hbp_chronic", 
              "rf_cervix", "rf_prev_4kg", "rf_prev_sga")
catVars <- c("year","sex","married","mrace","m_edu","pncgov","firstborn",
             "kotck", "m_wg_cat", "rf_db_gest", "rf_db_other", "rf_hbp_pregn", "rf_hbp_chronic", 
             "rf_cervix", "rf_prev_4kg", "rf_prev_sga")
bcVars <- c("bc_30d",
            "bc_90d","bc_280d",
            "bc_3090d", "bc_90280d")

rawtable1 <- CreateTableOne(vars = listVars, factorVars = catVars, 
                            strata = "lbw", addOverall = TRUE, test = FALSE,
                            data = data)
table1 <- print(rawtable1, nonnormal = bcVars,
                formatOption = list(decimal.mark = ".",  big.mark = ",", scientific = FALSE), contDigits = 3)

rtffile <- RTF(file = paste0(dir_output_table1, "table1.doc"))  # this can be an .rtf or a .doc
addParagraph(rtffile, "Table")
addTable(rtffile, cbind(rownames(table1), table1))
done(rtffile)

######################## 3. Check variables ###################################
## check categorical variables with multiple levels
## check normality on continuous variables: take log for those skewed variables
attach(data)

par(mfrow=c(2,3))
plot(year, main="year")
plot(mrace, main="mrace")
plot(m_edu, main="m_edu")
plot(kotck, main="kotck")
plot(m_wg_cat, main="maternal weight change")

par(mfrow=c(3,2))
hist(bwg)
hist(bc_30d)
hist(bc_90d)
hist(bc_280d)
hist(bc_3090d)
hist(bc_90280d)

par(mfrow=c(1,2))
plot(density(mage,  bw=2), main = "mage")
hist(clinega)

# cigdpp %>% sort(decreasing = TRUE) %>% head(30)
par(mfrow=c(2,2))
plot(density(cigdpp,  bw=0.5))
hist(log(cigdpp))
plot(density(log(cigdpp),  bw=0.5))
summary(cigdpp)

par(mfrow=c(2,2))
plot(density(cigddp,  bw=0.5))
hist(log(cigddp))
plot(density(log(cigddp),  bw=0.5))
summary(cigddp)

# med_hs_inc %>% sort(decreasing = TRUE) %>% head(30)
par(mfrow=c(2,2))
plot(density(mhincome, bw=1))
hist(log(mhincome))
plot(density(log(mhincome),  bw=0.5))
data$log_mhincome <- log(birth$mhincome)
summary(mhincome)

par(mfrow=c(2,2))
plot(density(mhvalue, bw=1))
hist(log(mhvalue))
plot(density(log(mhvalue),  bw=0.5))
data$log_mhvalue <- log(birth$mhvalue)
summary(mhvalue)

detach(data)

write.csv(data, "/Users/shuxind/Desktop/BC_birthweight_data/birth_final.csv")