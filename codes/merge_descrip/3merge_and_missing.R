###############################################################################
# Project: Causal black carbon on birth weight in MA                          #
# Code: merge datasets and detect missingness                                 #
# Input: "mabirths_02NOV18.csv" original birth data                           #
# Input: "birth_SES.csv" birth data with merged SES variables                 #
# Output: "birth_complt.csv" with complete info                               #
# Author: Shuxin Dong                                                         #
# Date: Oct 28, 2020                                                          #
###############################################################################

############################# 0. Setup ########################################
rm(list = ls())
gc()

library(data.table)
library(dplyr)

dir_input_birth <- "/Users/shuxind/Desktop/BC_birthweight_data/MAbirth_02Nov18/"
dir_input_SES <- "/Users/shuxind/Desktop/BC_birthweight_data/MAbirth_w_SES/"

birth0 <- fread(paste0(dir_input_birth, "mabirths_02NOV18.csv"), header = TRUE)
dim(birth0) # dimension of original birth data
# > dim(birth0)
# [1] 1199288      55

############################# 1. Merge and clean ##############################
## only include those born between 2001 and 2015
birth0 <- birth0 %>% filter(!year==2000)
dim(birth0) # dimension after omitting 2000 born babies
# > dim(birth0)
# [1] 1119011      55

## merge SES variables get "birthSES"
SES <- fread(paste0(dir_input_SES, "birth_SES.csv"), drop = c("V1","year"))
dim(SES) # SES variable dimension
# > dim(SES) # SES variable dimension
# [1] 1174107       4
birthSES <- left_join(birth0, SES, by = "uniqueid_yr")
remove(SES)
dim(birthSES) # dimension of birth data with SES info
# > dim(birthSES) # dimension of birth data with SES info
# [1] 1119011      58
remove(birth0)

## merge BC exposure data into birthSES and get "birth"
bc <- readRDS("/Users/shuxind/Desktop/BC_birthweight_data/bc_mean.rds")
birthBC <- left_join(birthSES, bc, by = "uniqueid_yr")
remove(bc)
remove(birthSES)
birth <- birthBC # create a copy
dim(birth)[1] # study population
# > dim(birth)[1] # study population
# [1] 1119011

## Re-encode to incorporate missingness (kotck & mrace)
birth$mrace[birth$mrace==6] <- NA # mrace==6 coded as missing
birth$kotck[birth$kotck==0] <- NA # kotck==0 coded as missing

####################### 2. Restrict: complete exposure history ################
## exclude those with missing average exposure
birth <- birth %>% filter(!is.na(bc_30d) & !is.na(bc_90d) & !is.na(bc_280d)) 
print(paste("1st: Number of observation with complete exposure history", dim(birth)[1]))
# [1] "1st: Number of observation with complete exposure history 1093420"

###################### 3. Inclusion and Exclusion criteria ####################
birth <- birth %>% filter(plurality==1 & baby_alive==1) # singleton and live birth (from 2001 to 2015)
birth <- birth %>% filter(bwg>=500 & bwg<=6000) # exclude those with birthweight less than 500g and larger than 6000g
birth <- birth %>% filter(clinega>=37 & clinega<=42) # clinical gestational age between 37 and 42 weeks
print(paste("2st: Number of observation after inclusion and exclusion", dim(birth)[1]))
# [1] "2st: Number of observation after inclusion and exclusion 927592"

############################# 4. Variable Selection ###########################
## PS model include all the variables could predict the outcomes 
## as long as they are not colliders
names(birth)
# > names(birth)
# [1] "uniqueid_yr"      "year"             "sex"              "married"          "mage"            
# [6] "fage"             "mrace"            "frace"            "m_edu"            "cigdpp"          
# [11] "cigddp"           "modvag"           "modfor"           "modvac"           "modpcs"          
# [16] "modrcs"           "modvbac"          "parit"            "gacalc"           "clinega"         
# [21] "kotck"            "kess"             "pncgov"           "bdob"             "date1prenat_vis" 
# [26] "bwg"              "tract"            "long"             "lat"              "gravid"          
# [31] "apgar1"           "apgar5"           "jaund"            "lbnl"             "lbnd"            
# [36] "date_last_menses" "rf_lung"          "rf_anem"          "rf_card"          "rf_db_gest"      
# [41] "rf_db_other"      "rf_hbp_chronic"   "rf_hbp_pregn"     "rf_cervix"        "rf_prev_bd"      
# [46] "rf_prev_4kg"      "rf_prev_sga"      "rf_renal"         "rf_rh_sens"       "rf_sickle"       
# [51] "rf_ut_bld"        "baby_alive"       "plurality"        "m_wg"             "BLOCK10"         
# [56] "mhincome"         "mhvalue"          "percentPoverty"   "bc_30d"           "bc_90d"          
# [61] "bc_280d"  

## drop father related variables | precise location info | date of first prenatal care visit
## drop the geographic `BLOCK10`, `tract`
birth <- birth %>% select(-fage, - frace, -long, -lat, -date1prenat_vis, -BLOCK10, -tract)

## Drop the downstream variables of outcome
birth <- birth %>% select(-apgar1, -apgar5, -jaund, -lbnl, -lbnd) # the current status of live birth (alive/dead), Agnar score, and Jaundice

## Drop the variables with duplicated info:
## drop the variables used in the criteria
## choose Kotelchuck over Kessner index
## drop baby's date of birth and date of last menses for gestational age
## drop other gestational age variables
birth <- birth %>% select(-plurality, -baby_alive)
birth <- birth %>% select(-kess, -bdob, -date_last_menses, -gacalc)

## Drop uncommon predictors:
## drop some risk factors according to Joel
## drop the mode of delivery
birth <- birth %>% select(-rf_lung, -rf_anem, -rf_card, -rf_prev_bd, -rf_renal,-rf_rh_sens, -rf_sickle, -rf_ut_bld, -modvag, -modfor, -modvac, -modpcs, -modrcs, -modvbac)
birth <- birth %>% select(-gravid)

############################# 5. Summarize the current ########################
summary(birth)
# > summary(birth)
# uniqueid_yr             year           sex          married            mage           mrace      
# Length:927592      Min.   :2001   Min.   :1.00   Min.   :0.0000   Min.   :12.00   Min.   :1.000  
# Class :character   1st Qu.:2004   1st Qu.:1.00   1st Qu.:0.0000   1st Qu.:26.00   1st Qu.:1.000  
# Mode  :character   Median :2007   Median :1.00   Median :1.0000   Median :30.66   Median :1.000  
# Mean   :2008   Mean   :1.49   Mean   :0.6853   Mean   :30.12   Mean   :1.673  
# 3rd Qu.:2011   3rd Qu.:2.00   3rd Qu.:1.0000   3rd Qu.:34.33   3rd Qu.:2.000  
# Max.   :2015   Max.   :2.00   Max.   :1.0000   Max.   :58.08   Max.   :5.000  
# NA's   :2        NA's   :1       NA's   :3066   
#      m_edu           cigdpp             cigddp             parit            clinega     
#  Min.   :1.000   Min.   :   0.000   Min.   :  0.0000   Min.   :  0.000   Min.   :37.00  
#  1st Qu.:2.000   1st Qu.:   0.000   1st Qu.:  0.0000   1st Qu.:  1.000   1st Qu.:39.00  
#  Median :3.000   Median :   0.000   Median :  0.0000   Median :  2.000   Median :39.00  
#  Mean   :3.171   Mean   :   1.595   Mean   :  0.5404   Mean   :  1.852   Mean   :39.32  
#  3rd Qu.:4.000   3rd Qu.:   0.000   3rd Qu.:  0.0000   3rd Qu.:  2.000   3rd Qu.:40.00  
#  Max.   :5.000   Max.   :1200.000   Max.   :466.6667   Max.   :199.000   Max.   :42.00  
#  NA's   :3564    NA's   :707        NA's   :709        NA's   :975                      
#      kotck           pncgov            bwg         rf_db_gest      rf_db_other     rf_hbp_chronic  
#  Min.   :1.000   Min.   :0.0000   Min.   : 505   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000  
#  1st Qu.:3.000   1st Qu.:0.0000   1st Qu.:3118   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000  
#  Median :3.000   Median :0.0000   Median :3430   Median :0.0000   Median :0.0000   Median :0.0000  
#  Mean   :3.095   Mean   :0.3362   Mean   :3439   Mean   :0.0425   Mean   :0.0085   Mean   :0.0121  
#  3rd Qu.:4.000   3rd Qu.:1.0000   3rd Qu.:3742   3rd Qu.:0.0000   3rd Qu.:0.0000   3rd Qu.:0.0000  
#  Max.   :4.000   Max.   :1.0000   Max.   :5982   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000  
#  NA's   :13299   NA's   :1324                    NA's   :1505     NA's   :1505     NA's   :1505    
# rf_hbp_pregn      rf_cervix       rf_prev_4kg      rf_prev_sga          m_wg          mhincome     
# Min.   :0.0000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000   Min.   :-88.0   Min.   :  5833  
# 1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.: 22.0   1st Qu.: 43571  
# Median :0.0000   Median :0.0000   Median :0.0000   Median :0.0000   Median : 30.0   Median : 59787  
# Mean   :0.0347   Mean   :0.0037   Mean   :0.0066   Mean   :0.0079   Mean   : 29.9   Mean   : 62435  
# 3rd Qu.:0.0000   3rd Qu.:0.0000   3rd Qu.:0.0000   3rd Qu.:0.0000   3rd Qu.: 37.0   3rd Qu.: 77061  
# Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :601.0   Max.   :226181  
# NA's   :1505     NA's   :1505     NA's   :1505     NA's   :1505     NA's   :18586   NA's   :129     
# mhvalue        percentPoverty       bc_30d           bc_90d          bc_280d      
# Min.   :  51840   Min.   : 0.000   Min.   :0.1545   Min.   :0.1853   Min.   :0.1997  
# 1st Qu.: 206570   1st Qu.: 4.200   1st Qu.:0.3779   1st Qu.:0.3830   1st Qu.:0.3893  
# Median : 272710   Median : 7.776   Median :0.4453   Median :0.4472   Median :0.4514  
# Mean   : 297317   Mean   :11.928   Mean   :0.4757   Mean   :0.4761   Mean   :0.4780  
# 3rd Qu.: 352880   3rd Qu.:16.400   3rd Qu.:0.5431   3rd Qu.:0.5391   3rd Qu.:0.5352  
# Max.   :1333300   Max.   :83.000   Max.   :1.8477   Max.   :1.7935   Max.   :1.6123  
# NA's   :12873     NA's   :110 
mod1 <- lm(bwg ~ ., data = subset(birth, select = -c(uniqueid_yr, bc_30d, bc_90d, bc_280d)))
summary(mod1)
# > summary(mod1)
# 
# Call:
#   lm(formula = bwg ~ ., data = subset(birth, select = -c(uniqueid_yr, 
#                                                          bc_30d, bc_90d, bc_280d)))
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -5640.3  -281.0   -14.7   267.7  5446.6 
# 
# Coefficients:
#   Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)     2.927e+03  2.381e+02   12.294  < 2e-16 ***
#   year           -2.292e+00  1.184e-01  -19.362  < 2e-16 ***
#   sex            -1.296e+02  9.056e-01 -143.154  < 2e-16 ***
#   married         3.526e+01  1.276e+00   27.632  < 2e-16 ***
#   mage            2.000e+00  1.015e-01   19.707  < 2e-16 ***
#   mrace          -2.344e+01  3.901e-01  -60.083  < 2e-16 ***
#   m_edu           7.362e+00  5.001e-01   14.723  < 2e-16 ***
#   cigdpp         -8.257e-01  1.063e-01   -7.770 7.87e-15 ***
#   cigddp         -1.341e+01  2.285e-01  -58.678  < 2e-16 ***
#   parit           5.007e+01  4.785e-01  104.642  < 2e-16 ***
#   clinega         1.276e+02  4.078e-01  312.931  < 2e-16 ***
#   kotck           4.285e+00  5.367e-01    7.984 1.42e-15 ***
#   pncgov         -2.147e+01  1.289e+00  -16.663  < 2e-16 ***
#   rf_db_gest      9.760e+01  2.285e+00   42.712  < 2e-16 ***
#   rf_db_other     1.690e+02  5.037e+00   33.553  < 2e-16 ***
#   rf_hbp_chronic -5.261e+01  4.177e+00  -12.595  < 2e-16 ***
#   rf_hbp_pregn   -4.570e+01  2.497e+00  -18.304  < 2e-16 ***
#   rf_cervix      -1.057e+02  7.639e+00  -13.841  < 2e-16 ***
#   rf_prev_4kg     4.372e+02  5.588e+00   78.251  < 2e-16 ***
#   rf_prev_sga    -1.241e+02  5.123e+00  -24.223  < 2e-16 ***
#   m_wg            5.154e+00  3.536e-02  145.768  < 2e-16 ***
#   mhincome        2.168e-04  3.572e-05    6.069 1.28e-09 ***
#   mhvalue        -7.872e-05  5.234e-06  -15.041  < 2e-16 ***
#   percentPoverty -1.705e+00  7.047e-02  -24.195  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 424.9 on 881642 degrees of freedom
# (45926 observations deleted due to missingness)
# Multiple R-squared:  0.1858,	Adjusted R-squared:  0.1858 
# F-statistic:  8746 on 23 and 881642 DF,  p-value: < 2.2e-16
library(questionr)
freq.na(birth[,-1])
# > freq.na(birth[,-1])
# missing %
# m_wg             18586 2
# kotck            13299 1
# mhvalue          12873 1
# m_edu             3564 0
# mrace             3066 0
# rf_db_gest        1505 0
# rf_db_other       1505 0
# rf_hbp_chronic    1505 0
# rf_hbp_pregn      1505 0
# rf_cervix         1505 0
# rf_prev_4kg       1505 0
# rf_prev_sga       1505 0
# pncgov            1324 0
# parit              975 0
# cigddp             709 0
# cigdpp             707 0
# mhincome           129 0
# percentPoverty     110 0
# married              2 0
# mage                 1 0
# year                 0 0
# sex                  0 0
# clinega              0 0
# bwg                  0 0
# bc_30d               0 0
# bc_90d               0 0
# bc_280d              0 0

birth_complt <- birth %>% na.omit()
write.csv(birth_complt, file = "/Users/shuxind/Desktop/BC_birthweight_data/birth_complt.csv")
print(paste("Final Sample Size is", dim(birth_complt)[1]))
# [1] "Final Sample Size is 881666"