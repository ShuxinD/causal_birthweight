###############################################################################
# Project: Causal black carbon and NO2 on birth weight in MA                  #
# Code: merge datasets and detect missingness                                 #
# Input: "mabirths_02NOV18.csv" original birth data                           #
# Input: "birth_SES.csv", "no2_birth_average.csv", "bc_birth_average.csv"     #
# Output: "MAbirth_merged.csv"                                                #
# Author: Shuxin Dong                                                         #
# Date: 2021-01-19                                                            #
###############################################################################

############################# 0. Setup ########################################
rm(list = ls())
gc()

library(data.table)

setwd("/media/gate/Shuxin/")
dir_input <- "/media/gate/Shuxin/MAbirth/"
dir_input_SES <- "/media/gate/Shuxin/MAbirth/MAbirth_SES/"

birth0 <- fread(paste0(dir_input, "mabirths_02NOV18.csv"), header = TRUE)
setDT(birth0)

## only include those born between 2001 and 2015
birth0 <- birth0[year!=2000,]
dim(birth0) # dimension after omitting 2000 born babies
# > dim(birth0)
# [1] 1119011      55
######################## 1. read in bc, no2, SES ##############################
## read SES data
SES <- fread(paste0(dir_input_SES, "birth_SES.csv"), drop = "year")
dim(SES) # SES variable dimension
# > dim(SES) # SES variable dimension
# [1] 1174107       4

## read BC data
bc <- fread(paste0(dir_input, "bc_birth_average.csv"), drop = c("bc_90d", "bc_280d"))
dim(bc)
# > dim(bc)
# [1] 1093528       4

## read NO2 data
no2 <- fread(paste0(dir_input, "no2_birth_average.csv"))
dim(no2)
# > dim(no2)
# [1] 1094520       4
######################## 2. merge together ####################################
birth <- merge(birth0, SES, by = "uniqueid_yr", all.x = TRUE)
birth <- merge(birth, bc, by = "uniqueid_yr", all.x = TRUE)
birth <- merge(birth, no2, by = "uniqueid_yr", all.x = TRUE)
remove(birth0)
gc()
################# 3. re-encoding to incorporate missingness ###################
birth$mrace[birth$mrace==6] <- NA # mrace==6 coded as missing
birth$kotck[birth$kotck==0] <- NA # kotck==0 coded as missing
##################### 4. Restrict to complete exposure history ################
## exclude those with missing average exposures
birth_compltE <- birth[!is.na(bc_30d) & !is.na(bc_3090d) & !is.na(bc_90280d) & !is.na(no2_30d)
      & !is.na(no2_3090d) & !is.na(no2_90280d),]
print(paste("1st: Number of observation with complete exposure history is", dim(birth_compltE)[1]))
# [1] "1st: Number of observation with complete exposure history is 1040055"
print(paste("Exclude", dim(birth)[1] - dim(birth_compltE)[1], "from", dim(birth)[1], 
            "as",(dim(birth)[1] - dim(birth_compltE)[1])/dim(birth)[1]))
# [1] "Exclude 78956 from 1119011 as 0.0705587344539062"

###################### 5. Inclusion and Exclusion criteria ####################
## include singleton and live birth
## exclude those with birthweight less than 500g and larger than 6000g
## include clinical gestational age between 37 and 42 weeks
birth_inout <- 
  birth_compltE[plurality==1 & baby_alive==1 & bwg>=500 & bwg<=6000 & clinega>=37 & clinega<=42,]
print(paste("2nd: Number of observation after inclusion and exclusion is", dim(birth_inout)[1]))
# [1] "2nd: Number of observation after inclusion and exclusion is 886893"
print(paste("Exclude", dim(birth_compltE)[1] - dim(birth_inout)[1], "from", dim(birth_compltE)[1], 
            "as",(dim(birth_compltE)[1] - dim(birth_inout)[1])/dim(birth_compltE)[1]))
# [1] "Exclude 153162 from 1040055 as 0.147263365879689"

############################# 6. Variable Selection ###########################
## PS model include all the variables could predict the outcomes 
## as long as they are not colliders
names(birth_inout)
# > names(birth_inout)
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
# [56] "mhincome"         "mhvalue"          "percentPoverty"   "bc_30d"           "bc_3090d"        
# [61] "bc_90280d"        "no2_30d"          "no2_3090d"        "no2_90280d" 
## drop father related variables | location info | date of first prenatal care visit
birth_inout[, ':=' (fage = NULL, frace = NULL, 
                     long = NULL, lat = NULL, BLOCK10 = NULL, tract = NULL,
                     date1prenat_vis = NULL)][]
## drop the downstream variables of outcome
birth_inout[, ':=' (apgar1 = NULL, apgar5 = NULL, jaund = NULL, lbnl = NULL, lbnd = NULL)][]

## drop the variables with duplicated info:
birth_inout[, ':=' (plurality = NULL, baby_alive = NULL, # drop the variables used in the criteria
                     kess = NULL, # choose Kotelchuck over Kessner index
                     bdob = NULL, date_last_menses = NULL, # drop baby's date of birth and date of last menses for gestational age
                     gacalc = NULL)] # drop other gestational age variables

## drop uncommon predictors:
birth_inout[, ':='( # drop some risk factors according to Joel
                    rf_lung = NULL, rf_anem = NULL, rf_card = NULL, rf_prev_bd = NULL, 
                    rf_renal = NULL, rf_rh_sens = NULL, rf_sickle = NULL, 
                    rf_ut_bld = NULL, modvag = NULL, modfor = NULL, modvac = NULL, 
                    modpcs = NULL, modrcs = NULL, modvbac = NULL,
                    # drop the mode of delivery
                    gravid = NULL)]
dim(birth_inout)
# > dim(birth_inout)
# [1] 886893     31
############################# 7. Summarize the current ########################
summary(birth_inout)
# > summary(birth_inout)
# uniqueid_yr             year           sex          married            mage           mrace      
# Length:886893      Min.   :2001   Min.   :1.00   Min.   :0.0000   Min.   :12.00   Min.   :1.000  
# Class :character   1st Qu.:2004   1st Qu.:1.00   1st Qu.:0.0000   1st Qu.:26.00   1st Qu.:1.000  
# Mode  :character   Median :2007   Median :1.00   Median :1.0000   Median :30.59   Median :1.000  
# Mean   :2007   Mean   :1.49   Mean   :0.6857   Mean   :30.10   Mean   :1.685  
# 3rd Qu.:2011   3rd Qu.:2.00   3rd Qu.:1.0000   3rd Qu.:34.33   3rd Qu.:2.000  
# Max.   :2015   Max.   :2.00   Max.   :1.0000   Max.   :58.08   Max.   :5.000  
# NA's   :2        NA's   :1       NA's   :2468   
#      m_edu          cigdpp            cigddp             parit            clinega     
#  Min.   :1.00   Min.   :   0.00   Min.   :  0.0000   Min.   :  0.000   Min.   :37.00  
#  1st Qu.:2.00   1st Qu.:   0.00   1st Qu.:  0.0000   1st Qu.:  1.000   1st Qu.:39.00  
#  Median :3.00   Median :   0.00   Median :  0.0000   Median :  2.000   Median :39.00  
#  Mean   :3.16   Mean   :   1.62   Mean   :  0.5484   Mean   :  1.851   Mean   :39.32  
#  3rd Qu.:4.00   3rd Qu.:   0.00   3rd Qu.:  0.0000   3rd Qu.:  2.000   3rd Qu.:40.00  
#  Max.   :5.00   Max.   :1200.00   Max.   :466.6667   Max.   :199.000   Max.   :42.00  
#  NA's   :2803   NA's   :707       NA's   :709        NA's   :952                      
#      kotck           pncgov            bwg         rf_db_gest      rf_db_other     rf_hbp_chronic 
#  Min.   :1.000   Min.   :0.0000   Min.   : 505   Min.   :0.0000   Min.   :0.0000   Min.   :0.000  
#  1st Qu.:3.000   1st Qu.:0.0000   1st Qu.:3118   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.000  
#  Median :3.000   Median :0.0000   Median :3430   Median :0.0000   Median :0.0000   Median :0.000  
#  Mean   :3.095   Mean   :0.3352   Mean   :3439   Mean   :0.0418   Mean   :0.0086   Mean   :0.012  
#  3rd Qu.:4.000   3rd Qu.:1.0000   3rd Qu.:3742   3rd Qu.:0.0000   3rd Qu.:0.0000   3rd Qu.:0.000  
#  Max.   :4.000   Max.   :1.0000   Max.   :5982   Max.   :1.0000   Max.   :1.0000   Max.   :1.000  
#  NA's   :12564   NA's   :1131                    NA's   :1493     NA's   :1493     NA's   :1493   
# rf_hbp_pregn     rf_cervix       rf_prev_4kg      rf_prev_sga          m_wg       
# Min.   :0.000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000   Min.   :-88.00  
# 1st Qu.:0.000   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.: 22.00  
# Median :0.000   Median :0.0000   Median :0.0000   Median :0.0000   Median : 30.00  
# Mean   :0.034   Mean   :0.0038   Mean   :0.0066   Mean   :0.0078   Mean   : 29.82  
# 3rd Qu.:0.000   3rd Qu.:0.0000   3rd Qu.:0.0000   3rd Qu.:0.0000   3rd Qu.: 36.00  
# Max.   :1.000   Max.   :1.0000   Max.   :1.0000   Max.   :1.0000   Max.   :601.00  
# NA's   :1493    NA's   :1493     NA's   :1493     NA's   :1493     NA's   :16527   
#     mhincome         mhvalue        percentPoverty       bc_30d          bc_3090d     
#  Min.   :  5833   Min.   :  51840   Min.   : 0.000   Min.   :0.1631   Min.   :0.1788  
#  1st Qu.: 43398   1st Qu.: 205500   1st Qu.: 4.177   1st Qu.:0.3806   1st Qu.:0.3838  
#  Median : 59382   Median : 271490   Median : 7.700   Median :0.4483   Median :0.4501  
#  Mean   : 62042   Mean   : 295287   Mean   :11.867   Mean   :0.4785   Mean   :0.4790  
#  3rd Qu.: 76570   3rd Qu.: 350650   3rd Qu.:16.288   3rd Qu.:0.5462   3rd Qu.:0.5440  
#  Max.   :226181   Max.   :1294600   Max.   :83.000   Max.   :1.8477   Max.   :1.7767  
#  NA's   :123      NA's   :12600     NA's   :105                                       
# bc_90280d         no2_30d          no2_3090d        no2_90280d    
# Min.   :0.1969   Min.   : 0.8736   Min.   : 1.045   Min.   : 1.454  
# 1st Qu.:0.3915   1st Qu.:16.1290   1st Qu.:16.568   1st Qu.:17.709  
# Median :0.4542   Median :22.6335   Median :22.969   Median :23.081  
# Mean   :0.4813   Mean   :22.7929   Mean   :22.985   Mean   :23.127  
# 3rd Qu.:0.5401   3rd Qu.:28.9280   3rd Qu.:29.037   3rd Qu.:28.484  
# Max.   :1.6426   Max.   :78.8389   Max.   :83.849   Max.   :70.033 
mod <- lm(bwg ~ ., data = subset(birth_inout, select = -c(uniqueid_yr, 
                                                          bc_30d, bc_3090d, bc_90280d,
                                                          no2_30d, no2_3090d, no2_90280d)))
summary(mod)
# > summary(mod)
# 
# Call:
#   lm(formula = bwg ~ ., data = subset(birth_inout, select = -c(uniqueid_yr, 
#                                                                bc_30d, bc_3090d, bc_90280d, no2_30d, no2_3090d, no2_90280d)))
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -5630.0  -281.1   -14.7   267.8  5440.5 
# 
# Coefficients:
#   Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)     2.539e+03  2.558e+02    9.925  < 2e-16 ***
#   year           -2.116e+00  1.272e-01  -16.639  < 2e-16 ***
#   sex            -1.292e+02  9.254e-01 -139.611  < 2e-16 ***
#   married         3.559e+01  1.307e+00   27.222  < 2e-16 ***
#   mage            1.998e+00  1.036e-01   19.292  < 2e-16 ***
#   mrace          -2.299e+01  3.957e-01  -58.093  < 2e-16 ***
#   m_edu           7.376e+00  5.110e-01   14.436  < 2e-16 ***
#   cigdpp         -7.965e-01  1.074e-01   -7.414 1.23e-13 ***
#   cigddp         -1.342e+01  2.312e-01  -58.062  < 2e-16 ***
#   parit           4.999e+01  4.891e-01  102.212  < 2e-16 ***
#   clinega         1.285e+02  4.182e-01  307.306  < 2e-16 ***
#   kotck           3.689e+00  5.491e-01    6.718 1.84e-11 ***
#   pncgov         -2.154e+01  1.321e+00  -16.312  < 2e-16 ***
#   rf_db_gest      9.801e+01  2.355e+00   41.617  < 2e-16 ***
#   rf_db_other     1.653e+02  5.139e+00   32.161  < 2e-16 ***
#   rf_hbp_chronic -5.043e+01  4.289e+00  -11.758  < 2e-16 ***
#   rf_hbp_pregn   -4.508e+01  2.579e+00  -17.483  < 2e-16 ***
#   rf_cervix      -1.062e+02  7.743e+00  -13.712  < 2e-16 ***
#   rf_prev_4kg     4.398e+02  5.684e+00   77.385  < 2e-16 ***
#   rf_prev_sga    -1.297e+02  5.272e+00  -24.606  < 2e-16 ***
#   m_wg            5.147e+00  3.633e-02  141.670  < 2e-16 ***
#   mhincome        2.491e-04  3.684e-05    6.762 1.36e-11 ***
#   mhvalue        -8.132e-05  5.406e-06  -15.043  < 2e-16 ***
#   percentPoverty -1.660e+00  7.235e-02  -22.942  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 424.9 on 844530 degrees of freedom
# (42339 observations deleted due to missingness)
# Multiple R-squared:  0.1863,	Adjusted R-squared:  0.1863 
# F-statistic:  8405 on 23 and 844530 DF,  p-value: < 2.2e-16

melt(birth_inout[, lapply(.SD, function(x) sum(is.na(x)))])[order(-value)][, freq := value/dim(birth_inout)[1]][]
# variable value         freq
# 1:           m_wg 16527 1.863472e-02
# 2:        mhvalue 12600 1.420690e-02
# 3:          kotck 12564 1.416631e-02
# 4:          m_edu  2803 3.160471e-03
# 5:          mrace  2468 2.782748e-03
# 6:     rf_db_gest  1493 1.683405e-03
# 7:    rf_db_other  1493 1.683405e-03
# 8: rf_hbp_chronic  1493 1.683405e-03
# 9:   rf_hbp_pregn  1493 1.683405e-03
# 10:      rf_cervix  1493 1.683405e-03
# 11:    rf_prev_4kg  1493 1.683405e-03
# 12:    rf_prev_sga  1493 1.683405e-03
# 13:         pncgov  1131 1.275238e-03
# 14:          parit   952 1.073410e-03
# 15:         cigddp   709 7.994200e-04
# 16:         cigdpp   707 7.971649e-04
# 17:       mhincome   123 1.386864e-04
# 18: percentPoverty   105 1.183908e-04
# 19:        married     2 2.255063e-06
# 20:           mage     1 1.127532e-06
# 21:    uniqueid_yr     0 0.000000e+00
# 22:           year     0 0.000000e+00
# 23:            sex     0 0.000000e+00
# 24:        clinega     0 0.000000e+00
# 25:            bwg     0 0.000000e+00
# 26:         bc_30d     0 0.000000e+00
# 27:       bc_3090d     0 0.000000e+00
# 28:      bc_90280d     0 0.000000e+00
# 29:        no2_30d     0 0.000000e+00
# 30:      no2_3090d     0 0.000000e+00
# 31:     no2_90280d     0 0.000000e+00
# variable value         freq

birth_for_analysis <- na.omit(birth_inout)
fwrite(birth_for_analysis, paste0(dir_input, "MAbirth_merged.csv"))
print(paste("Final Sample Size is", dim(birth_for_analysis)[1]))
# [1] "Final Sample Size is 844554"