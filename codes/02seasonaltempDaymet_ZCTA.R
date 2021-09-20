#' Project: beta radiation and CVD death in MA
#' Code: convert daily temperature to seasonal average for each year
#' Input: Daymet daily temperature files on qnap4 prepared by Edgar
#' Output: seasonal_tmean_2001_2015.csv
#' Author: Shuxin Dong                                                         

## 0. set up ------------------------------------------------------------------
rm(list = ls())
gc()

setwd("/media/qnap3/Shuxin/ParticalRadiation_MAdeath/")

library(data.table)
library(sf)

## the average temperature is defined as the average of the daily maximum and minimum temperatures for the season.
## 1. load temperatures  ---------------------------------------------------------
dir_temperature <- "/media/qnap4/Daymet/geoid-issues/census-aggregates/" # currently only with miscoding zcta
crosswalk <- fread("/media/qnap4/Daymet/geoid-issues/zcta_crosswalk.csv", colClasses = c("zcta" = "character"))
zcta_MA <- setDT(st_transform(sf::st_read("tl_2010_25_zcta510.shp"), 4326))[,.(STATEFP10, ZCTA5CE10)]
crosswalk_MA <- crosswalk[zcta %in% zcta_MA[,ZCTA5CE10],]

seasonaltemp <- NULL
for (year_i in 2001:2015){
  dailytmax <- fread(paste0(dir_temperature, "zcta_daily_mean_tmax_", year_i, ".csv.gz"), colClasses = c("geoid" = "character"))
  dailytmax <- dailytmax[geoid %in% zcta_MA[,ZCTA5CE10],]
  cat("load", year_i, "daily tmax\n")
  dailytmin <- fread(paste0(dir_temperature, "zcta_daily_mean_tmin_", year_i, ".csv.gz"))
  dailytmin <- merge(crosswalk_MA, dailytmin, by.x = "factor_level", by.y = "geoid")
  dailytmin[, factor_level:=NULL]
  colnames(dailytmin)[1] <- "geoid"
  dailytmin <- dailytmin[geoid %in% zcta_MA[,ZCTA5CE10],]
  cat("load", year_i, "daily tmin\n")
  
  if (year_i %in% c(2004, 2008, 2012)){ # leap year but all the data without 12/31 so ---
    dailytmax[, ":=" (summer_tmax = rowMeans(dailytmax[,(153+1):(244+1)], na.rm = TRUE),
                      winter_tmax = rowMeans(dailytmax[,c(2:(60+1), (336+1):366)], na.rm = TRUE))]
    dailytmin[, ":=" (summer_tmin = rowMeans(dailytmin[,(153+1):(244+1)], na.rm = TRUE),
                      winter_tmin = rowMeans(dailytmin[,c(2:(60+1), (336+1):366)], na.rm = TRUE))]
  } else {
    dailytmax[, ":=" (summer_tmax = rowMeans(dailytmax[,153:244], na.rm = TRUE),
                      winter_tmax = rowMeans(dailytmax[,c(2:60, 336:366)], na.rm = TRUE))]
    dailytmin[, ":=" (summer_tmin = rowMeans(dailytmin[,153:244], na.rm = TRUE),
                      winter_tmin = rowMeans(dailytmin[,c(2:60, 336:366)], na.rm = TRUE))]
  }
  seasonaltmax <- dailytmax[,.(geoid, summer_tmax, winter_tmax)]
  seasonaltmin <- dailytmin[,.(geoid, summer_tmin, winter_tmin)]
  
  seasonal_c <- merge(seasonaltmax, seasonaltmin, by = "geoid")
  seasonal_c[, ":=" (summer_tmean = (summer_tmax+summer_tmin)/2,
                     winter_tmean = (winter_tmax+winter_tmin)/2,
                     year = year_i)]
  seasonaltemp <- rbind(seasonaltemp, seasonal_c[,.(geoid, summer_tmean, winter_tmean, year)])
}

rm(dailytmax)
rm(dailytmin)
rm(seasonaltmax)
rm(seasonaltmin)
rm(seasonal_c)
gc()

head(seasonaltemp)
# geoid summer_tmean winter_tmean year
# 1: 01001     21.52846    -1.465979 2001
# 2: 01002     20.34408    -2.874618 2001
# 3: 01003     20.85149    -2.610668 2001
# 4: 01005     19.80180    -3.174362 2001
# 5: 01007     20.54024    -2.574793 2001
# 6: 01008     18.83868    -3.190020 2001
colnames(seasonaltemp)[1] <- "zcta"
dir <- getwd()
fwrite(seasonaltemp, paste0(dir, "/seasonal_tmean_2001_2015.csv"))

