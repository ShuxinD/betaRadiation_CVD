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

## 1.1 annual death on ZCTA ---------------------------------------------------
dir_death_count <- "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/"
count_month <- readRDS(paste0(dir_death_count, "MAdeath_count_ZIP.rds")) # this is by year and month
count_year <- count_month[, .(CVD = sum(CVD),
                              MI = sum(MI),
                              CHF = sum(CHF),
                              stroke = sum(stroke),
                              TOT = sum(TOT)), by = .(year, ZCTA5CE10)] # aggregate on year
summary(count_year[,year])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2000    2003    2008    2008    2012    2015 
sum(count_year[,TOT])
# [1] 738913
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

## 2.4 other census variables -------------------------------------------------------
census_ses <- readRDS("census0115.rds")
census_ses <- census_ses[(zip %in% count_year[, ZCTA5CE10]),]

## 3. merge all together --------
## 3.1 merge exposure and health outcome -----
dt <- merge(count_year, betaR, by.x = c("year", "ZCTA5CE10"), by.y = c("Year", "ZIP"), all.x = TRUE) # X and Y from 2001 to 2015
dim(dt)
# [1] 7703    8
sum(dt[is.na(dt[,Beta]), TOT])
# [1] 10652 # number of individuals without Beta info
sum(dt[, TOT])
# [1] 738896
sum(dt[is.na(dt[,Beta]), TOT])/sum(dt[, TOT])
# [1] 0.0144161

sum(is.na(dt[,Beta]))
# [1] 528
dt <- na.omit(dt)
sum(dt[, TOT])
# [1] 728244

## 3.2 merge covariates -----
dt <- merge(dt, seasonaltemp, by.x = c("year", "ZCTA5CE10"), by.y = c("year", "zcta"), all.x = TRUE)
dt <- merge(dt, pm25, by.x = c("year", "ZCTA5CE10"), by.y = c("year", "ZIP"), all.x = TRUE)
dt <- merge(dt, pcount_long, by.x = c("year", "ZCTA5CE10"), by.y = c("year", "ZCTA5"), all.x = TRUE)
dt <- merge(dt, census_ses, by.x = c("year", "ZCTA5CE10"), by.y = c("year", "zip"), all.x = TRUE)

sum(dt[!complete.cases(dt), TOT]) # 2975 people without covariates info
sum(dt[!complete.cases(dt), TOT])/sum(dt[, TOT]) # [1] 0.004085169

final <- na.omit(dt)
sum(final[,TOT])
# [1] 725269
dim(final)
# [1] 6985   21
saveRDS(final, "finalDT.rds")
