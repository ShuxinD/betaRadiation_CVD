##################################################################
# Project: beta radiation and CVD death in MA
# Code: convert MA death record from one-row-per-person to one-row-per-zcta
# Input: "ma_2000_2015.sas7bdat"
# Output: "zcta_year.rds" as combination of zcta and year
# Output: "MAdeath_count_ZIP_age'group'.rds"
# Author: Shuxin Dong
##################################################################

## 0. set up --------------------------------------------------------------
rm(list = ls())
gc()

setwd("/media/qnap3/Shuxin/ParticalRadiation_MAdeath/")

library(haven)
library(data.table)
library(tigris)
library(sf)
# library(icd)
# library(sjlabelled)

## 1.ICD codes ---------------------------------------------------------
# - total cardiovascular disease (CVD): I00-I51
# - congestive heart failure (CHF): I50-I51
# - myocardial infarction (MI): I21-I22
# - stroke: I60-I69
# - all-cause mortality (TOT): A00-R99

## 2. load MA death records  -------------------------------------
origin <- read_sas("/media/qnap2/MortalityStates/GeocodedMortality/ma_2000_2015.sas7bdat")
# get_label(origin)
# year 
# "Year of death" 
# date 
# "Date of death" 
# age 
# "Age in years" 
# sex 
# "Sex (1 = Male, 2 = Female)" 
# icd 
# "TRX_CAUSE_ACME" 
# lat 
# "Latitude" 
# long 
# "Longitude" 
# edu3 
# "Education recode (1 = Less than HS, 2 = HS, 3 = More than HS)" 
# race3 
# "Race recode (1 = White, 2 = Black, 3 = Other)" 
# fips 
# "FIPs County code, 5-digit number in text format" 
# ind 
# "Industry (text or number code of unknown nature)" 
# occ 
# "Occupation (text or number code of unknown nature)" 
# mar 
# "Marital Status (1 = never, 2 = married/separated, 3 = widowed, 4 = divorced)" 
# hisp 
# "Hispanic status (1 = yes, 0 = no)" 
# state 
# "Reporting State, as well as State of Residence, and State of Occurence" 
setDT(origin)

## add IDs
origin[, ID:=1:.N][]
dim(origin)
# [1] 805563     16

## only 2001-2015
origin <- origin[year(date) %in% 2001:2015,]
dim(origin)
# [1] 752805     16

## only 18+
origin <- origin[age>=18,]

## remove missing
cat("Number of individual missing location info:", sum(is.na(origin[,lat])), "/", dim(origin)[1], sum(is.na(origin[,lat]))/dim(origin)[1])
# Number of individual missing location info: 13477 / 743873 0.01811734
cat("Number of individual missing location info:", sum(is.na(origin[,long])), "/", dim(origin)[1], sum(is.na(origin[,long]))/dim(origin)[1])
# Number of individual missing location info: 13477 / 743873 0.01811734
cat("Number of individual missing ICD info:", sum(origin[,icd]==""), "/", dim(origin)[1], sum(origin[,icd]=="")/dim(origin)[1])
# Number of individual missing ICD info: 259 / 743873 0.0003481777
cat("Number of individual missing age info:", sum(is.na(origin[,age])), "/", dim(origin)[1], sum(is.na(origin[,age]))/dim(origin)[1])
# Number of individual missing age info: 0 / 743873 0

## missingness for location and ICD
cat("Number of individual missing location info or ICD or:", sum(is.na(origin[,lat])| origin[,icd]=="" | is.na(origin[,age])), "/", dim(origin)[1], sum(is.na(origin[,lat])| origin[,icd]=="" | is.na(origin[,age]))/dim(origin)[1])
# Number of individual missing location info or ICD or: 13722 / 743873 0.0184467

mort <- na.omit(origin[,.(year, date, age, icd, lat, long, ID)])
mort <- mort[icd != "",]
dim(mort)
# [1] 730151      7
head(mort)
# year       date age  icd      lat      long    ID
# 1: 2001 2001-04-01  73 J449 42.08182 -72.61861 52758
# 2: 2001 2001-04-01  80 J449 42.05481 -72.62351 52760
# 3: 2001 2001-08-01  88 I219 42.08514 -72.62851 52761
# 4: 2001 2001-12-01  98  F03 42.08185 -72.61878 52762
# 5: 2001 2001-02-01  66  I64 42.08078 -72.59399 52763
# 6: 2001 2001-07-01  76  F03 42.07990 -72.61476 52764

## check age ---------------------
summary(mort[,age])
hist(mort[,age]) #  age 0
## check the death cause for age 0
tb <- sort(table(mort[age==0, icd]), decreasing = T)
write.csv(tb, file = "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/age0_icd.csv")

mort[, `:=` (icd_str1 = substr(icd,1,1),
           icd_str2 = substr(icd,1,2),
           icd_str3 = substr(icd,1,3))]
mort[, `:=` (CVD_TF = (icd_str2 %in% c("I0", "I1", "I2", "I3", "I4") | icd_str3 %in% c("I50", "I51")),
           MI_TF = icd_str3 %in% c("I21", "I22"),
           CHF_TF = icd_str3 %in% c("I50", "I51"),
           stroke_TF = icd_str3 %in% c("I60", "I61", "I62", "I63", "I64", "I65", "I66", "I67", "I68", "I69"),
           TOT_TF = icd_str1 %in% LETTERS[1:18])]

pdf(file = "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/TOT_age_distribution.pdf", width = 9, height = 15)
par(mfrow=c(5,3))
for (year_i in 2001:2015){
  hist(mort[year==year_i,age], main = paste("Total death age distribution in", year_i))
  print(year_i)
}
dev.off()

pdf(file = "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/CVD_age_distribution.pdf", width = 9, height = 15)
par(mfrow=c(5,3))
for (year_i in 2001:2015){
  hist(mort[CVD_TF & year==year_i,age], main = paste("CVD age distribution in", year_i)) 
  print(year_i)
}
dev.off()

pdf(file = "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/MI_age_distribution.pdf", width = 9, height = 15)
par(mfrow=c(5,3))
for (year_i in 2001:2015){
  hist(mort[MI_TF & year==year_i,age], main = paste("MI age distribution in", year_i)) 
  print(year_i)
}
dev.off()

pdf(file = "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/CHF_age_distribution.pdf", width = 9, height = 15)
par(mfrow=c(5,3))
for (year_i in 2001:2015){
  hist(mort[CHF_TF & year==year_i,age], main = paste("CHF age distribution in", year_i)) 
  print(year_i)
}
dev.off()

pdf(file = "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/stroke_age_distribution.pdf", width = 9, height = 15)
par(mfrow=c(5,3))
for (year_i in 2001:2015){
  hist(mort[stroke_TF & year==year_i,age], main = paste("stroke age distribution in", year_i)) 
  print(year_i)
}
dev.off()

## 3. load ZCTA shape file - the 2010 version ---------------------------------
point_dt <- st_as_sf(mort, coords = c("long", "lat"), crs = 4326) # spatial data
# download.file("ftp://ftp2.census.gov/geo/tiger/TIGER2010/ZCTA5/2010/tl_2010_25_zcta510.zip", "tl_2010_25_zcta510.zip")
# unzip("tl_2010_25_zcta510.zip")
zcta_shp10 <- st_transform(sf::st_read("tl_2010_25_zcta510.shp"), 4326)
zcta10 <- st_join(point_dt, zcta_shp10, join = st_intersects)
setDT(zcta10)

zcta_list <- setDT(zcta_shp10[, "ZCTA5CE10"])[,ZCTA5CE10]
zcta_list <- rep(zcta_list, 15)
year_list <- NULL
for (year_i in 2001:2015) {
  year_list <- c(year_list, rep(year_i, 538))
}
zcta_year <- setDT(data.frame(ZCTA5CE10 = zcta_list, 
                              year = year_list))
saveRDS(zcta_year, "zcta_year.rds")

## 4. aggregate MA death on ZCTA per year -------------------------------
## assign ZCTA to each row of death record
mort_loc <- merge(mort, zcta10[, .(ID, ZCTA5CE10)], by = "ID")
## create age label for each group
mort_loc[,`:=` (age18_65 = age>=18&age<65, 
                age65_85 = age>=65&age<85, 
                age85 = age>=85)]
## aggregate death records on ZCTA for year and month, and by age
dt <- mort_loc[, .(ID, year, date, icd, ZCTA5CE10, age18_65, age65_85, age85)]
dt[, month := month(date)]
dt[, `:=` (icd_str1 = substr(icd,1,1),
           icd_str2 = substr(icd,1,2),
           icd_str3 = substr(icd,1,3))]
dt[, `:=` (CVD_TF = (icd_str2 %in% c("I0", "I1", "I2", "I3", "I4") | icd_str3 %in% c("I50", "I51")),
           MI_TF = icd_str3 %in% c("I21", "I22"),
           CHF_TF = icd_str3 %in% c("I50", "I51"),
           stroke_TF = icd_str3 %in% c("I60", "I61", "I62", "I63", "I64", "I65", "I66", "I67", "I68", "I69"),
           TOT_TF = icd_str1 %in% LETTERS[1:18])]

setorder(dt, year, month, ZCTA5CE10)[]
count_all <- dt[, .(CVD = sum(CVD_TF),
                MI = sum(MI_TF),
                CHF = sum(CHF_TF),
                stroke = sum(stroke_TF),
                TOT = sum(TOT_TF)), by = .(year, month, ZCTA5CE10)]
count18_65 <- dt[(age18_65), .(CVD = sum(CVD_TF),
                             MI = sum(MI_TF),
                             CHF = sum(CHF_TF),
                             stroke = sum(stroke_TF),
                             TOT = sum(TOT_TF)), by = .(year, month, ZCTA5CE10)]
count65_85 <- dt[(age65_85), .(CVD = sum(CVD_TF),
                               MI = sum(MI_TF),
                               CHF = sum(CHF_TF),
                               stroke = sum(stroke_TF),
                               TOT = sum(TOT_TF)), by = .(year, month, ZCTA5CE10)]
count85 <- dt[(age85), .(CVD = sum(CVD_TF),
                         MI = sum(MI_TF),
                         CHF = sum(CHF_TF),
                         stroke = sum(stroke_TF),
                         TOT = sum(TOT_TF)), by = .(year, month, ZCTA5CE10)]
# sum <- c(sum(count_all[,TOT]), sum(count0_65[,TOT]), sum(count65_85[,TOT]), sum(count85[,TOT]))
# sum(sum[2:4])
# [1] 737383
sum(count_all[is.na(ZCTA5CE10),TOT])
# [1] 17 # number of individuals without assigned ZCTA
final_all <- na.omit(count_all)
final18_65 <- na.omit(count18_65)
final65_85 <- na.omit(count65_85)
final85 <- na.omit(count85)

sum(final_all[,TOT])
# [1] 730134
sum(final18_65[,TOT])
# [1] 131102
sum(final65_85[,TOT])
# [1] 323645
sum(final85[,TOT])
# [1] 275387

saveRDS(final_all, "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/MAdeath_count_ZIP_ageall.rds")
saveRDS(final18_65, "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/MAdeath_count_ZIP_age1865.rds")
saveRDS(final65_85, "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/MAdeath_count_ZIP_age6585.rds")
saveRDS(final85, "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/MAdeath_count_ZIP_age85.rds")
