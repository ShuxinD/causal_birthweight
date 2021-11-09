#' Project: Causal black carbon and NO2 on birth weight in MA
#' Code: Table one and normality
#' Input: "MAbirth_merged.csv"
#' Output: "MAbirth_for_analyses.csv" for future analyses
#' Author: Shuxin Dong
#' First create date: 2021-01-21.

## 0. Setup -----
rm(list = ls())
gc()

library(data.table)
library(fastDummies)
library(tableone)
library(rtf)

setwd("/media/gate/Shuxin/")
dir_in <- "/media/qnap3/Shuxin/airPollution_MAbirth/"
dir_out <- "/media/qnap3/Shuxin/airPollution_MAbirth/causal_birthweight/results/0TableOne_checkVarDistribution/"

birth <- fread(paste0(dir_in, "MAbirth_merged.csv"))
birth[, lbw:=0][]
birth$lbw[birth$bwg<2500] <- 1
names(birth)
# [1] "uniqueid_yr"    "year"           "sex"            "married"        "mage"           "mrace"         
# [7] "m_edu"          "cigdpp"         "cigddp"         "parit"          "clinega"        "kotck"         
# [13] "pncgov"         "bdob"           "bwg"            "rf_db_gest"     "rf_db_other"    "rf_hbp_chronic"
# [19] "rf_hbp_pregn"   "rf_cervix"      "rf_prev_4kg"    "rf_prev_sga"    "m_wg"           "mhincome"      
# [25] "mhvalue"        "percentPoverty" "bc_30d"         "bc_3090d"       "bc_90280d"      "no2_30d"       
# [31] "no2_3090d"      "no2_90280d"     "bc_all"         "no2_all"        "lbw"       

## 1. Prepare for table one ----
#' Transform continuous to categorical
#'`parity` to `firstborn`, because of some extreme large unexplained value.
birth[, firstborn := parit][]
birth$firstborn[birth$parit>1] <- 0
birth[, parit := NULL]
#' `m_wg` to `m_wg_cat` based on reference. Extreme large unexplained values exist.
plot(density(birth$m_wg, bw=1))
hist(birth$m_wg)
birth[, m_wg_cat := 0][]
birth$m_wg_cat[birth$m_wg<0] <- 1
birth$m_wg_cat[birth$m_wg>=0 & birth$m_wg<15] <- 2
birth$m_wg_cat[birth$m_wg>=15 & birth$m_wg<25] <- 3
birth$m_wg_cat[birth$m_wg>=25 & birth$m_wg<36] <- 4
birth$m_wg_cat[birth$m_wg>=36] <- 5
birth[, m_wg :=NULL]
#' create non-smoker variable as binary variable
birth$smoker_ddp <- 0
birth$smoker_ddp[birth$cigddp!=0] <- 1
birth$smoker_dpp <- 0
birth$smoker_dpp[birth$cigdpp!=0] <- 1
summary(as.factor(birth$smoker_ddp))
summary(as.factor(birth$smoker_dpp))

#' Get categorical and dummy variables
birth$mrace[birth$mrace==5] <- 4
birth$year <- as.factor(birth$year)
birth$sex <- birth$sex - 1 # 1 is girl; 0 is boy
birth$m_edu <- as.factor(birth$m_edu)
birth$kotck <- as.factor(birth$kotck)
birth$m_wg_cat <- as.factor(birth$m_wg_cat)
birth <- fastDummies::dummy_cols(birth, select_columns = "mrace")

## 2. Check variables ----
#' check normality on continuous variables: take log for those skewed variables
#'  check extreme values
attach(birth)

pdf(paste0(dir_out,"check1.pdf"))
par(mfrow=c(2,3))
plot(year, main="year")
plot(as.factor(mrace), main="mrace")
plot(m_edu, main="m_edu")
plot(kotck, main="kotck")
plot(m_wg_cat, main="maternal weight change")
dev.off()

pdf(paste0(dir_out,"check2.pdf"))
par(mfrow = c(1,1))
hist(bwg)
dev.off()

pdf(paste0(dir_out,"check3.pdf"))
par(mfrow=c(3,2))
hist(bc_30d)
hist(bc_3090d)
hist(bc_90280d)
hist(no2_30d)
hist(no2_3090d)
hist(no2_90280d)
dev.off()

pdf(paste0(dir_out,"check4.pdf"))
par(mfrow=c(1,2))
plot(density(mage,  bw=2), main = "mage")
hist(clinega)
dev.off()

# cigdpp %>% sort(decreasing = TRUE) %>% head(30)
pdf(paste0(dir_out,"check5.pdf"))
par(mfrow=c(2,2))
plot(density(cigdpp,  bw=0.5))
hist(log(cigdpp))
plot(density(log(cigdpp),  bw=0.5))
dev.off()
# max_cigdpp <- quantile(cigdpp, 0.99995)
summary(cigdpp)

pdf(paste0(dir_out,"check6.pdf"))
par(mfrow=c(2,2))
plot(density(cigddp,  bw=0.5))
hist(log(cigddp))
plot(density(log(cigddp),  bw=0.5))
dev.off()
# max_cigddp <- quantile(cigddp, 0.99995)
summary(cigddp)

# med_hs_inc %>% sort(decreasing = TRUE) %>% head(30)
pdf(paste0(dir_out,"check7.pdf"))
par(mfrow=c(2,2))
plot(density(mhincome, bw=1))
hist(log(mhincome))
plot(density(log(mhincome),  bw=0.5))
dev.off()

summary(mhincome)
birth[, log_mhincome:= log(mhincome)][]

pdf(paste0(dir_out,"check8.pdf"))
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
# [1] "uniqueid_yr"    "year"           "sex"            "married"        "mage"          
# [6] "mrace"          "m_edu"          "cigdpp"         "cigddp"         "clinega"       
# [11] "kotck"          "pncgov"         "bwg"            "rf_db_gest"     "rf_db_other"   
# [16] "rf_hbp_chronic" "rf_hbp_pregn"   "rf_cervix"      "rf_prev_4kg"    "rf_prev_sga"   
# [21] "mhincome"       "mhvalue"        "percentPoverty" "bc_30d"         "bc_3090d"      
# [26] "bc_90280d"      "no2_30d"        "no2_3090d"      "no2_90280d"     "bc_all"        
# [31] "no2_all"        "lbw"            "firstborn"      "m_wg_cat"       "smoker_ddp"    
# [36] "smoker_dpp"     "mrace_1"        "mrace_2"        "mrace_3"        "mrace_4"       
# [41] "log_mhincome"   "log_mhvalue"  
# summary(birth)

## 3. add the born season variable ----
birth[, bdob:=as.Date(bdob, "%m/%d/%Y")]
birth[,`:=`(b_spring = (month(bdob) %in% c(3:5)),
            b_summer = (month(bdob) %in% c(6:8)),
            b_autumn = (month(bdob) %in% c(9:11)),
            b_winter = (month(bdob) %in% c(1,2,12)))]

## 4. Export table 1 ------
listVars <- c("bwg",
              "bc_all","bc_30d","bc_3090d", "bc_90280d",
              "no2_all","no2_30d", "no2_3090d", "no2_90280d",
              "year", "b_spring", "b_summer", "b_autumn", "b_winter",
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
bcVars <- c("bc_all", "bc_30d", "bc_3090d", "bc_90280d")
no2Vars <- c("no2_all", "no2_30d", "no2_3090d", "no2_90280d")

rawtable1 <- CreateTableOne(vars = listVars, factorVars = catVars, 
                            strata = "lbw", addOverall = TRUE, test = FALSE,
                            data = birth)
print(rawtable1, nonnormal = c(bcVars, no2Vars),
      formatOption = list(decimal.mark = ".",  big.mark = ",", scientific = FALSE),
      contDigits = 2)
table1 <- print(rawtable1, nonnormal = c(bcVars, no2Vars),
                formatOption = list(decimal.mark = ".",  big.mark = ",", scientific = FALSE),
                contDigits = 2)

rtffile <- RTF(file = paste0(dir_out, "table1.doc"))  # this can be an .rtf or a .doc
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
write.csv(smokerInfo, paste0(dir_out, "table1_smokerInfo.csv"))

mean_bc <- rbind(c(mean(birth[,bc_all]), mean(birth[lbw==0, bc_all]),mean(birth[lbw==1, bc_all])),
                 c(mean(birth[,bc_30d]), mean(birth[lbw==0, bc_30d]),mean(birth[lbw==1, bc_30d])),
                 c(mean(birth[,bc_3090d]), mean(birth[lbw==0, bc_3090d]),mean(birth[lbw==1, bc_3090d])),
                 c(mean(birth[,bc_90280d]), mean(birth[lbw==0, bc_90280d]),mean(birth[lbw==1, bc_90280d])))
rownames(mean_bc) <- c("bc_all","bc_30d", "bc_3090d", "bc_90280d")
colnames(mean_bc) <- c("overall", "nbw", "lbw")
mean_bc
#           overall       nbw       lbw
# bc_all    0.4790591 0.4788522 0.4885682
# bc_30d    0.4770028 0.4767743 0.4875047
# bc_3090d  0.4775469 0.4773306 0.4874869
# bc_90280d 0.4798614 0.4796608 0.4890775

IQR_bc <- rbind(quantile(birth[,bc_all]), quantile(birth[lbw==0, bc_all]),quantile(birth[lbw==1, bc_all]),
                quantile(birth[,bc_30d]), quantile(birth[lbw==0, bc_30d]),quantile(birth[lbw==1, bc_30d]),
                 quantile(birth[,bc_3090d]), quantile(birth[lbw==0, bc_3090d]),quantile(birth[lbw==1, bc_3090d]),
                 quantile(birth[,bc_90280d]), quantile(birth[lbw==0, bc_90280d]),quantile(birth[lbw==1, bc_90280d]))
IQR_bc
rownames(IQR_bc) <- paste0(rep(c("overall", "nbw", "lbw"),4), 
                           c(rep("_bc_all",3), rep("_bc_30d",3), rep("_bc_3090d",3), rep("_bc_90280d",3)))
IQR_bc
#                         0%       25%       50%       75%     100%
# overall_bc_all    0.1997006 0.3918021 0.4538842 0.5364093 1.612341
# nbw_bc_all        0.1997006 0.3917655 0.4536912 0.5360798 1.612341
# lbw_bc_all        0.2231169 0.3934087 0.4631700 0.5516666 1.351354
# overall_bc_30d    0.1630848 0.3807746 0.4480162 0.5445652 1.826903
# nbw_bc_30d        0.1630848 0.3806807 0.4478211 0.5442752 1.826903
# lbw_bc_30d        0.1866778 0.3853661 0.4580873 0.5599305 1.473415
# overall_bc_3090d  0.1788135 0.3840173 0.4498046 0.5424744 1.776729
# nbw_bc_3090d      0.1788135 0.3839377 0.4496089 0.5421738 1.776729
# lbw_bc_3090d      0.1994498 0.3879015 0.4602769 0.5570112 1.354540
# overall_bc_90280d 0.1968862 0.3917479 0.4540661 0.5386503 1.642583
# nbw_bc_90280d     0.1968862 0.3917036 0.4539207 0.5383319 1.642583
# lbw_bc_90280d     0.2221745 0.3937925 0.4616684 0.5540890 1.401050

mean_no2 <- rbind(c(mean(birth[,no2_all]), mean(birth[lbw==0, no2_all]),mean(birth[lbw==1, no2_all])),
                 c(mean(birth[,no2_30d]), mean(birth[lbw==0, no2_30d]),mean(birth[lbw==1, no2_30d])),
                 c(mean(birth[,no2_3090d]), mean(birth[lbw==0, no2_3090d]),mean(birth[lbw==1, no2_3090d])),
                 c(mean(birth[,no2_90280d]), mean(birth[lbw==0, no2_90280d]),mean(birth[lbw==1, no2_90280d])))
rownames(mean_no2) <- c("no2_all","no2_30d", "no2_3090d", "no2_90280d")
colnames(mean_no2) <- c("overall", "nbw", "lbw")
mean_no2
#             overall      nbw      lbw
# no2_all    23.00261 22.99377 23.40882
# no2_30d    22.74911 22.73950 23.19087
# no2_3090d  22.93494 22.92671 23.31304
# no2_90280d 23.06401 23.05510 23.47348

IQR_no2 <- rbind(quantile(birth[,no2_all]), quantile(birth[lbw==0, no2_all]),quantile(birth[lbw==1, no2_all]),
                quantile(birth[,no2_30d]), quantile(birth[lbw==0, no2_30d]),quantile(birth[lbw==1, no2_30d]),
                quantile(birth[,no2_3090d]), quantile(birth[lbw==0, no2_3090d]),quantile(birth[lbw==1, no2_3090d]),
                quantile(birth[,no2_90280d]), quantile(birth[lbw==0, no2_90280d]),quantile(birth[lbw==1, no2_90280d]))
IQR_no2
rownames(IQR_no2) <- paste0(rep(c("overall", "nbw", "lbw"),4), 
                           c(rep("_no2_all",3), rep("_no2_30d",3), rep("_no2_3090d",3), rep("_no2_90280d",3)))
IQR_no2
#                         0%      25%      50%      75%     100%
# overall_no2_all    2.1801237 18.02544 22.86641 27.82676 62.59741
# nbw_no2_all        2.1801237 18.01744 22.85781 27.81661 62.59741
# lbw_no2_all        3.3845054 18.40268 23.28793 28.29211 58.75638
# overall_no2_30d    0.8735659 16.09391 22.58374 28.88830 75.22628
# nbw_no2_30d        0.8735659 16.08221 22.56944 28.88265 75.22628
# lbw_no2_30d        1.4735506 16.68021 23.22411 29.16093 62.31848
# overall_no2_3090d  1.0452284 16.53250 22.92165 28.97721 83.84932
# nbw_no2_3090d      1.0842349 16.52219 22.91188 28.96856 83.84932
# lbw_no2_3090d      1.0452284 16.99337 23.35559 29.29263 62.30313
# overall_no2_90280d 1.4542685 17.66789 23.02038 28.41093 64.61834
# nbw_no2_90280d     1.4542685 17.65838 23.01136 28.40323 64.61834
# lbw_no2_90280d     3.4550755 18.14295 23.43051 28.75076 63.98316

## export data ready for analyses ----
summary(birth)
dir_outdata <- "/media/qnap3/Shuxin/airPollution_MAbirth/"
fwrite(birth, paste0(dir_outdata, "MAbirth_for_analyses.csv"))

birth <- fread(paste0(dir_outdata, "MAbirth_for_analyses.csv"))
# names(birth)
birth$year <- as.factor(birth$year)
birth$m_edu <- as.factor(birth$m_edu)
birth$kotck <- as.factor(birth$kotck)
birth$m_wg_cat <- as.factor(birth$m_wg_cat)
birth <- fastDummies::dummy_cols(birth, select_columns = c("m_edu", "kotck","m_wg_cat"))
fwrite(birth, paste0(dir_outdata, "MAbirth_for_analyses.csv"))

birth <- fread(paste0(dir_outdata, "MAbirth_for_analyses.csv"))
birth$bc_30280d <- (birth$bc_3090d*60 + birth$bc_90280d*190)/(60+190)
birth$no2_30280d <- (birth$no2_3090d*60 + birth$no2_90280d*190)/(60+190)
fwrite(birth, paste0(dir_outdata, "MAbirth_for_analyses.csv"))

## 5. correlation plot ----
birth <- fread(paste0(dir_in, "MAbirth_for_analyses.csv"))
#' for continuous variable
M <- cor(na.omit(birth)[,.(sex, married, mage, mrace_1, mrace_2, mrace_3, mrace_4, smoker_ddp, smoker_dpp, cigdpp, cigddp, clinega, pncgov, log_mhincome, log_mhvalue, firstborn, rf_db_gest, rf_db_other,  rf_hbp_chronic, rf_hbp_pregn, rf_cervix, rf_prev_4kg, rf_prev_sga, percentPoverty, bc_30d, bc_3090d, bc_90280d, no2_30d, no2_3090d, no2_90280d, bc_all, no2_all)])

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
p.mat <- cor.mtest(na.omit(birth)[,.(sex, married, mage, mrace_1, mrace_2, mrace_3, mrace_4, smoker_ddp, smoker_dpp, cigdpp, cigddp, clinega, pncgov, log_mhincome, log_mhvalue, firstborn, rf_db_gest, rf_db_other,  rf_hbp_chronic, rf_hbp_pregn, rf_cervix, rf_prev_4kg, rf_prev_sga, percentPoverty, bc_30d, bc_3090d, bc_90280d, no2_30d, no2_3090d, no2_90280d, bc_all, no2_all)])
# colnames(M) <- c("Particle radiation", "PM[2.5]", "Summer average temperature", "Winter average temperature", "Median household income", "Median value of house", "Percentage of Hispanic", "Percentage of black", "Pencentage below poverty line", "Percentage without high school diploma", "Population density", "Mean BMI", "Smoking rate")
# rownames(M) <- c("Particle radiation", "PM[2.5]", "Summer average temperature", "Winter average temperature", "Median household income", "Median value of house", "Percentage of Hispanic", "Percentage of black", "Pencentage below poverty line", "Percentage without high school diploma", "Population density", "Mean BMI", "Smoking rate")
par(mfrow=c(1,1))
corrplot(M, method="number", type = "lower", p.mat = p.mat, sig.level = 0.05)

pdf(paste0(dir_out,"corrTable.pdf"), width = 16, height = 16)
corrplot(M, method="number", type = "lower", p.mat = p.mat, sig.level = 0.05)
dev.off()


