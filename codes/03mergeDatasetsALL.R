###############################################################################
# Project: beta radiation and CVD death in MA
# Code: merge all datasets together for analyses
# Input: all datasets
# Output: "finalDT.rds"
# Author: Shuxin Dong                                                         
###############################################################################

## 0. set up ------------------------------------------------------------------
rm(list = ls())
gc()

setwd("/media/qnap3/Shuxin/ParticalRadiation_MAdeath/")

library(data.table)
## except for census data,data are all from Joel's server

## 0. load zcta and year combination ------------------------------------------
zcta_year <- readRDS("zcta_year.rds")

## 1.1 annual death on ZCTA ---------------------------------------------------
## 1.1.1 all age ---------
dir_death_count <- "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/"
count_month <- readRDS(paste0(dir_death_count, "MAdeath_count_ZIP_ageall.rds")) # this is by year and month
count_year <- count_month[, .(CVD = sum(CVD),
                              MI = sum(MI),
                              CHF = sum(CHF),
                              stroke = sum(stroke),
                              TOT = sum(TOT)), by = .(year, ZCTA5CE10)] # aggregate on year
summary(count_year[,year])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2001    2004    2008    2008    2012    2015 
sum(count_year[,TOT])
# [1] 730134

## 1.1.2. age 18-65 ----
dir_death_count <- "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/"
count_month1865 <- readRDS(paste0(dir_death_count, "MAdeath_count_ZIP_age1865.rds")) # this is by year and month
count_year1865 <- count_month1865[, .(CVD = sum(CVD),
                              MI = sum(MI),
                              CHF = sum(CHF),
                              stroke = sum(stroke),
                              TOT = sum(TOT)), by = .(year, ZCTA5CE10)] # aggregate on year
summary(count_year1865[,year])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2001    2004    2008    2008    2012    2015 
sum(count_year1865[,TOT])
names(count_year1865)[3:7] <- paste0(names(count_year1865)[3:7], "1865")

## 1.1.3. age 65-85 ----
dir_death_count <- "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/"
count_month6585 <- readRDS(paste0(dir_death_count, "MAdeath_count_ZIP_age6585.rds")) # this is by year and month
count_year6585 <- count_month6585[, .(CVD = sum(CVD),
                                    MI = sum(MI),
                                    CHF = sum(CHF),
                                    stroke = sum(stroke),
                                    TOT = sum(TOT)), by = .(year, ZCTA5CE10)] # aggregate on year
summary(count_year6585[,year])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2001    2004    2008    2008    2012    2015
sum(count_year6585[,TOT])
# [1] 323645
names(count_year6585)[3:7] <- paste0(names(count_year6585)[3:7], "6585")

## 1.1.4. age85+ ----
dir_death_count <- "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/"
count_month85 <- readRDS(paste0(dir_death_count, "MAdeath_count_ZIP_age85.rds")) # this is by year and month
count_year85 <- count_month85[, .(CVD = sum(CVD),
                                    MI = sum(MI),
                                    CHF = sum(CHF),
                                    stroke = sum(stroke),
                                    TOT = sum(TOT)), by = .(year, ZCTA5CE10)] # aggregate on year
summary(count_year85[,year])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2001    2004    2008    2008    2012    2015 
sum(count_year85[,TOT])
# [1] 275387
names(count_year85)[3:7] <- paste0(names(count_year85)[3:7], "85")

remove(count_month)
remove(count_month1865)
remove(count_month6585)
remove(count_month85)
gc()

## 1.2 annual beta radiation on ZCTA ------------------------------------------
load("MA_ZIP_Annual_Beta_V3.RData") # radiation data from Longxiang
setDT(ma_annual)
betaR <- ma_annual[,.(ZIP, Beta, Year)]
rm(ma_annual)
gc()
summary(betaR[,Year])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2001    2005    2009    2009    2013    2017 

## ** all following data should be 2001-2015 ---------------------------------
## 2.1 seasonal temperature on ZCTA ------------------------------------------
seasonaltemp <- fread("seasonal_tmean_2001_2015.csv", colClasses = c("zcta" = "character"))

## 2.2 annual PM2.5 on ZIP code -----------------------------------------------
dir_pm25 <- "/media/qnap3/Exposure modeling/3 National scale/USA/4 PM2.5 v2.2000_16/9 Zipcode average/Annual/"
pm25_files <- list.files(path = dir_pm25)[2:16]
pm25 <- NULL
for (i in 1:length(pm25_files)) {
  pm25_temp <- readRDS(paste0(dir_pm25, pm25_files[i]))
  setDT(pm25_temp)
  pm25_temp[, year:= i+2000]
  pm25 <- rbind(pm25, pm25_temp)
  cat("finish", 2000+i, "\n")
}
rm(pm25_temp)

## 2.3 population count by ZCTA -----------------------------------------
dir_pcount <- "/media/qnap3/Covariates/Census/PopulationCounts_zipcode_US_2000_2016/"
pcount_wide <- readRDS(paste0(dir_pcount, "pop_2000_2016_agr.rds"))
setDT(pcount_wide)
pcount_wide <- pcount_wide[(RACE=="ALL")&(SEX=="Both")&(AGELEVEL=="all")&(GEOID %in% count_year[,ZCTA5CE10]),]
head(pcount_wide)
pcount_wide[,`:=`(RACE=NULL, SEX=NULL, AGEGP=NULL, AGELEVEL=NULL)]
pcount_long <- melt(pcount_wide, id.vars = c("GEOID"), measure.vars = patterns("^count"), variable.name = "year", value.name = c("pcount"))
pcount_long[, year:= substr(year, 6, 9)][]
pcount_long[, year:= as.integer(year)]
names(pcount_long)[1] <- "ZCTA5"

## 2.4 other census variables -----------------------------------------------------
census_ses <- readRDS("census0115.rds")
census_ses <- census_ses[(zip %in% count_year[, ZCTA5CE10]),]

## 3. merge on zcta and year combination ---------------------------------------
## 3.1 merge count and exposure --------
## merge count
dt <- merge(zcta_year, count_year, by = c("year", "ZCTA5CE10"), all.x = T)
dt <- dt[,lapply(.SD,function(x){ifelse(is.na(x),0,x)})]
## merge beta exposure
dt <- merge(dt, betaR, by.x = c("year", "ZCTA5CE10"), by.y = c("Year", "ZIP"), all.x = TRUE) # X and Y from 2001 to 2015
dim(dt)
# [1] 8070    8
sum(dt[is.na(dt[,Beta]), TOT])
# [1] 10532 # number of individuals without Beta info
sum(dt[, TOT])
# [1] 738896
sum(dt[is.na(dt[,Beta]), TOT])/sum(dt[, TOT])
# [1] 0.01442475

sum(is.na(dt[,Beta]))
# [1] 765
dt <- na.omit(dt)
sum(dt[, TOT])
# [1] 719602

## 3.2 merge covariates -----
dt <- merge(dt, seasonaltemp, by.x = c("year", "ZCTA5CE10"), by.y = c("year", "zcta"), all.x = TRUE)
dt <- merge(dt, pm25, by.x = c("year", "ZCTA5CE10"), by.y = c("year", "ZIP"), all.x = TRUE)
dt <- merge(dt, pcount_long, by.x = c("year", "ZCTA5CE10"), by.y = c("year", "ZCTA5"), all.x = TRUE)
dt <- merge(dt, census_ses, by.x = c("year", "ZCTA5CE10"), by.y = c("year", "zip"), all.x = TRUE)

sum(dt[!complete.cases(dt), TOT]) # 2949 people without covariates info
sum(dt[!complete.cases(dt), TOT])/sum(dt[, TOT]) # [1] 0.004098099

final <- na.omit(dt)
sum(final[,TOT])
# [1] 716653
dim(final)
# [1] 7054   21
#saveRDS(final, "finalDT.rds")

#final_allage <- readRDS("finalDT.rds")
names(final)
final <- merge(final, count_year1865, by = c("year", "ZCTA5CE10"), all.x = T)
final <- merge(final, count_year6585, by = c("year", "ZCTA5CE10"), all.x = T)
final <- merge(final, count_year85, by = c("year", "ZCTA5CE10"), all.x = T)
final
anyNA(final)
final <- final[,lapply(.SD,function(x){ifelse(is.na(x),0,x)})] # convert NA to 0
anyNA(final)
final
saveRDS(final, "finalDT.rds")
