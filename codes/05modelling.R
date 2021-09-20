#' Project: beta radiation and CVD death in MA
#' Code: Statistical Modelling
#' Input: "finalDT.rds"
#' Output: "/modresults_details/"
#' Output: "RRiqr_beta_beta.csv"
#' Output: "RRiqr_beta_betaPM.csv"
#' Output: "RRiqr_pm25_betaPM.csv"
#' Author: Shuxin Dong                                                         

## 0. set up ------------------------------------------------------------------
rm(list = ls())
gc()

setwd("/media/qnap3/Shuxin/ParticalRadiation_MAdeath/")
library(mgcv)

dt <- readRDS("finalDT.rds")
dt <- na.omit(dt)
dt_pcount <- dt[,pcount]

## 1. models ------------------------------------------------------------------
dir_modresults <- "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/betaRadiation_CVD/results/modresults_details/"
deathc <- c("TOT", "CVD", "MI", "stroke", 
            paste0(c("TOT", "CVD", "MI", "stroke"), 1865), 
            paste0(c("TOT", "CVD", "MI", "stroke"), 6585), 
            paste0(c("TOT", "CVD", "MI", "stroke"), 85))

## 1.1 DID mod0 ---------------------------------------------------------------
## beta + PM25 setA -----
for (deathc_i in deathc) {
  mod <- gam(get(deathc_i) ~ Beta + pm25 + summer_tmean + winter_tmean + as.factor(year) + as.factor(ZCTA5CE10) + offset(log(dt_pcount)),
               data = dt,
               family = quasipoisson(link = "log"))
  tb <- summary(mod)$p.table
  write.table(tb, file = paste0(dir_modresults, "mod0A_", deathc_i, ".csv"))
  cat(paste0("finish mod0A ", deathc_i, "\n"))
}

## only beta setB-----
for (deathc_i in deathc) {
  mod <- gam(get(deathc_i) ~ Beta + summer_tmean + winter_tmean + as.factor(year) + as.factor(ZCTA5CE10) + offset(log(dt_pcount)),
               data = dt,
               family = quasipoisson(link = "log"))
  tb <- summary(mod)$p.table
  write.table(tb, file = paste0(dir_modresults, "mod0B_", deathc_i, ".csv"))
  cat(paste0("finish mod0B ", deathc_i, "\n"))
}

## 1.2 mixed-effect model mod1 ------------------------------------------------
## beta and PM25 setA----
for (deathc_i in deathc) {
  mod <- gamm(get(deathc_i) ~ Beta + pm25 + summer_tmean + winter_tmean + medhouseholdincome + medianhousevalue + hispanic + pct_blk + poverty + education + popdensity + mean_bmi + smoke_rate + as.factor(year) + offset(log(dt_pcount)),
              data = dt,
              random = list(ZCTA5CE10=~1),
              family = quasipoisson(link = "log"))
  tb <- summary(mod$gam)$p.table
  write.table(tb, file = paste0(dir_modresults, "mod1A_", deathc_i, ".csv"))
  cat(paste0("finish mod1A ", deathc_i, "\n"))
}

## only beta setB-----
for (deathc_i in deathc) {
  mod <- gamm(get(deathc_i) ~ Beta + summer_tmean + winter_tmean + medhouseholdincome + medianhousevalue + hispanic + pct_blk + poverty + education + popdensity + mean_bmi + smoke_rate + as.factor(year) + offset(log(dt_pcount)),
              data = dt,
              random = list(ZCTA5CE10=~1),
              family = quasipoisson(link = "log"))
  tb <- summary(mod$gam)$p.table
  write.table(tb, file = paste0(dir_modresults, "mod1B_", deathc_i, ".csv"))
  cat(paste0("finish mod1B ", deathc_i, "\n"))
}

## 2. Results tables ----------------------------------------------------------
dir_results.table <- "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/betaRadiation_CVD/results/"
deathc <- c("TOT", "CVD", "MI", "stroke", 
            paste0(c("TOT", "CVD", "MI", "stroke"), 1865), 
            paste0(c("TOT", "CVD", "MI", "stroke"), 6585), 
            paste0(c("TOT", "CVD", "MI", "stroke"), 85))
## IQR ----
IQRs <- setDT(data.frame(beta = IQR(dt[,Beta]),
                         pm25 = IQR(dt[,pm25])))
print(IQRs)
# beta     pm25
# 1: 0.05530013 2.819183

## setA (beta + PM25) RR of beta for an IQR increase in two models ----
## model 0A
RR0A_tb <- NULL
rownames_list <- NULL

for (deathc_i in deathc) {
  tb <- read.table(paste0(dir_modresults, "mod0A_", deathc_i, ".csv"))
  estimates <- exp(c(tb[2,1], tb[2,1]-1.96*tb[2,2], tb[2,1]+1.96*tb[2,2])*IQRs[,beta])
  RR0A_tb <- rbind(RR0A_tb, estimates)
  rownames_list <- c(rownames_list, paste0(deathc_i, "age"))
}
rownames(RR0A_tb) <- rownames_list
colnames(RR0A_tb) <- c("0ARR", "0ARR_lci", "0ARR_uci")
print(RR0A_tb)

## model 1A
RR1A_tb <- NULL
rownames_list <- NULL
for (deathc_i in deathc) {
  tb <- read.table(paste0(dir_modresults, "mod1A_", deathc_i, ".csv"))
  estimates <- exp(c(tb[2,1], tb[2,1]-1.96*tb[2,2], tb[2,1]+1.96*tb[2,2])*IQRs[,beta])
  RR1A_tb <- rbind(RR1A_tb, estimates)
  rownames_list <- c(rownames_list, paste0(deathc_i, "age"))
}
rownames(RR1A_tb) <- rownames_list
colnames(RR1A_tb) <- c("1ARR", "1ARR_lci", "1ARR_uci")
print(RR1A_tb)

results_A <- cbind(RR0A_tb, RR1A_tb)
print(results_A)

write.table(results_A, file = paste0(dir_results.table, "RRiqr_beta_betaPM.csv"))
# kable(results_A, caption = "Estimated RR for one IQR increase in beta radiation, with Exposure Set A (radiation + PM2.5)")

## setA (beta + PM25) RR of PM25 for one IQR increase in two models ----
## model 0A
RR0A_tb <- NULL
rownames_list <- NULL

for (deathc_i in deathc) {
  tb <- read.table(paste0(dir_modresults, "mod0A_", deathc_i, ".csv"))
  estimates <- exp(c(tb[3,1], tb[3,1]-1.96*tb[3,2], tb[3,1]+1.96*tb[3,2])*IQRs[,pm25])
  RR0A_tb <- rbind(RR0A_tb, estimates)
  rownames_list <- c(rownames_list, paste0(deathc_i, "age"))
}

rownames(RR0A_tb) <- rownames_list
colnames(RR0A_tb) <- c("0ARR", "0ARR_lci", "0ARR_uci")
print(RR0A_tb)

## model 1A
RR1A_tb <- NULL
rownames_list <- NULL
for (deathc_i in deathc) {
  tb <- read.table(paste0(dir_modresults, "mod1A_", deathc_i, ".csv"))
  estimates <- exp(c(tb[3,1], tb[3,1]-1.96*tb[3,2], tb[3,1]+1.96*tb[3,2])*IQRs[,pm25])
  RR1A_tb <- rbind(RR1A_tb, estimates)
  rownames_list <- c(rownames_list, paste0(deathc_i, "age"))
}
rownames(RR1A_tb) <- rownames_list
colnames(RR1A_tb) <- c("1ARR", "1ARR_lci", "1ARR_uci")
print(RR1A_tb)

results_A <- cbind(RR0A_tb, RR1A_tb)
print(results_A)

write.table(results_A, file = paste0(dir_results.table, "RRiqr_pm25_betaPM.csv"))

## setB (only beta): RR of beta for one IQR increase in two models ------
## model 0B
RR0B_tb <- NULL
rownames_list <- NULL
for (deathc_i in deathc) {
  tb <- read.table(paste0(dir_modresults, "mod0B_", deathc_i, ".csv"))
  estimates <- exp(c(tb[2,1], tb[2,1]-1.96*tb[2,2], tb[2,1]+1.96*tb[2,2])*IQRs[,beta])
  RR0B_tb <- rbind(RR0B_tb, estimates)
  rownames_list <- c(rownames_list, paste0(deathc_i, "age"))
}
rownames(RR0B_tb) <- rownames_list
colnames(RR0B_tb) <- c("0BRR", "0BRR_lci", "0BRR_uci")
print(RR0B_tb)

## model 1B
RR1B_tb <- NULL
rownames_list <- NULL
for (deathc_i in deathc) {
  tb <- read.table(paste0(dir_modresults, "mod1B_", deathc_i, ".csv"))
  estimates <- exp(c(tb[2,1], tb[2,1]-1.96*tb[2,2], tb[2,1]+1.96*tb[2,2])*IQRs[,beta])
  RR1B_tb <- rbind(RR1B_tb, estimates)
  rownames_list <- c(rownames_list, paste0(deathc_i, "age"))
}
rownames(RR1B_tb) <- rownames_list
colnames(RR1B_tb) <- c("1BRR", "1BRR_lci", "1BRR_uci")
print(RR1B_tb)

results_B <- cbind(RR0B_tb, RR1B_tb)
print(results_B)

write.table(results_B, file = paste0(dir_results.table, "RRiqr_beta_beta.csv"))
# kable(results_B, caption = "Estimated RR for one IQR increase in beta radiation, with Exposure Set B (only radiation)")
