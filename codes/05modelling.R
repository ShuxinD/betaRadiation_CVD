###############################################################################
# Project: beta radiation and CVD death in MA
# Code: Statistical Modelling
# Input: "finalDT.rds"
# Output: 
# Author: Shuxin Dong                                                         
###############################################################################

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
## 1.1 DID --------------------------------------------------------------------
## beta + PM25-----
deathc <- c("TOT", "CVD", "MI", "stroke", 
            paste0(c("TOT", "CVD", "MI", "stroke"), 1865), 
            paste0(c("TOT", "CVD", "MI", "stroke"), 6585), 
            paste0(c("TOT", "CVD", "MI", "stroke"), 85))
for (deathc_i in deathc) {
  mod <- gam(get(deathc_i) ~ Beta + pm25 + summer_tmean + winter_tmean + as.factor(year) + as.factor(ZCTA5CE10) + offset(log(dt_pcount)),
               data = dt,
               family = quasipoisson(link = "log"))
  tb <- summary(mod)$p.table
  write.table(tb, file = paste0(dir_modresults, "mod0A_", deathc_i, ".csv"))
  cat(paste0("finish mod0A ", deathc_i, "\n"))
}

## only beta -----
deathc <- c("TOT", "CVD", "MI", "stroke", 
            paste0(c("TOT", "CVD", "MI", "stroke"), 1865), 
            paste0(c("TOT", "CVD", "MI", "stroke"), 6585), 
            paste0(c("TOT", "CVD", "MI", "stroke"), 85))
for (deathc_i in deathc) {
  mod <- gam(get(deathc_i) ~ Beta + summer_tmean + winter_tmean + as.factor(year) + as.factor(ZCTA5CE10) + offset(log(dt_pcount)),
               data = dt,
               family = quasipoisson(link = "log"))
  tb <- summary(mod)$p.table
  write.table(tb, file = paste0(dir_modresults, "mod0B_", deathc_i, ".csv"))
  cat(paste0("finish mod0B ", deathc_i, "\n"))
}

## 1.2 mixed-effect model -----------------------------------------------------
## beta and PM25 ----
deathc <- c("TOT", "CVD", "MI", "stroke", 
            paste0(c("TOT", "CVD", "MI", "stroke"), 1865), 
            paste0(c("TOT", "CVD", "MI", "stroke"), 6585), 
            paste0(c("TOT", "CVD", "MI", "stroke"), 85))
for (deathc_i in deathc) {
  mod <- gamm(get(deathc_i) ~ Beta + pm25 + summer_tmean + winter_tmean + medhouseholdincome + medianhousevalue + hispanic + pct_blk + poverty + education + popdensity + mean_bmi + smoke_rate + as.factor(year) + offset(log(dt_pcount)),
              data = dt,
              random = list(ZCTA5CE10=~1),
              family = quasipoisson(link = "log"))
  tb <- summary(mod$gam)$p.table
  write.table(tb, file = paste0(dir_modresults, "mod1A_", deathc_i, ".csv"))
  cat(paste0("finish mod1A ", deathc_i, "\n"))
}

## only beta -----
deathc <- c("TOT", "CVD", "MI", "stroke", 
            paste0(c("TOT", "CVD", "MI", "stroke"), 1865), 
            paste0(c("TOT", "CVD", "MI", "stroke"), 6585), 
            paste0(c("TOT", "CVD", "MI", "stroke"), 85))
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
dir_results <- "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/betaRadiation_CVD/results/"
## IQR ----
IQRs <- setDT(data.frame(beta=IQR(dt[,Beta]),
                         pm25=IQR(dt[,pm25])))
print(IQRs)
# beta     pm25
# 1: 0.05530013 2.819183

## Exposure A (beta + PM25) RR of beta for one IQR increase in two models ----
## model 0A
deathc <- c("TOT", "CVD", "MI", "stroke", 
            paste0(c("TOT", "CVD", "MI", "stroke"), 1865), 
            paste0(c("TOT", "CVD", "MI", "stroke"), 6585), 
            paste0(c("TOT", "CVD", "MI", "stroke"), 85))
RR0A_tb <- NULL
rownames_list <- NULL

for (deathc_i in deathc) {
  tb <- read.table(paste0(dir_modresults, "mod0A_", deathc_i, ".csv"))
  estimates <- exp(c(tb[2,1], tb[2,1]-1.96*tb[2,2], tb[2,1]-1.96*tb[2,2])*IQRs[,beta])
  RR0A_tb <- rbind(RR0A_tb, estimates)
  rownames_list <- c(rownames_list, paste0(deathc_i, "age"))
}
rownames(RR0A_tb) <- rownames_list
colnames(RR0A_tb) <- c("0ARR", "0ARR_lci", "0ARR_uci")
print(RR0A_tb)

## model 1A
deathc <- c("TOT", "CVD", "MI", "stroke", 
            paste0(c("TOT", "CVD", "MI", "stroke"), 1865), 
            paste0(c("TOT", "CVD", "MI", "stroke"), 6585), 
            paste0(c("TOT", "CVD", "MI", "stroke"), 85))
RR1A_tb <- NULL
rownames_list <- NULL
for (deathc_i in deathc) {
  tb <- read.table(paste0(dir_modresults, "mod1A_", deathc_i, ".csv"))
  estimates <- exp(c(tb[2,1], tb[2,1]-1.96*tb[2,2], tb[2,1]-1.96*tb[2,2])*IQRs[,beta])
  RR1A_tb <- rbind(RR1A_tb, estimates)
  rownames_list <- c(rownames_list, paste0(deathc_i, "age"))
}
rownames(RR1A_tb) <- rownames_list
colnames(RR1A_tb) <- c("1ARR", "1ARR_lci", "1ARR_uci")
print(RR1A_tb)

results_A <- cbind(RR0A_tb, RR1A_tb)
print(results_A)

write.table(results_A, file = paste0(dir_results, "RRiqr_beta_betaPM.csv"))
# kable(results_A, caption = "Estimated RR for one IQR increase in beta radiation, with Exposure Set A (radiation + PM2.5)")

## Exposure A (beta + PM25) RR of PM25 for one IQR increase in two models ----
## model 0A
deathc <- c("TOT", "CVD", "MI", "stroke", 
            paste0(c("TOT", "CVD", "MI", "stroke"), 1865), 
            paste0(c("TOT", "CVD", "MI", "stroke"), 6585), 
            paste0(c("TOT", "CVD", "MI", "stroke"), 85))
RR0A_tb <- NULL
rownames_list <- NULL

for (deathc_i in deathc) {
  tb <- read.table(paste0(dir_modresults, "mod0A_", deathc_i, ".csv"))
  estimates <- exp(c(tb[3,1], tb[3,1]-1.96*tb[3,2], tb[3,1]-1.96*tb[3,2])*IQRs[,pm25])
  RR0A_tb <- rbind(RR0A_tb, estimates)
  rownames_list <- c(rownames_list, paste0(deathc_i, "age"))
}

rownames(RR0A_tb) <- rownames_list
colnames(RR0A_tb) <- c("0ARR", "0ARR_lci", "0ARR_uci")
print(RR0A_tb)

## model 1A
deathc <- c("TOT", "CVD", "MI", "stroke", 
            paste0(c("TOT", "CVD", "MI", "stroke"), 1865), 
            paste0(c("TOT", "CVD", "MI", "stroke"), 6585), 
            paste0(c("TOT", "CVD", "MI", "stroke"), 85))
RR1A_tb <- NULL
rownames_list <- NULL
for (deathc_i in deathc) {
  tb <- read.table(paste0(dir_modresults, "mod1A_", deathc_i, ".csv"))
  estimates <- exp(c(tb[3,1], tb[3,1]-1.96*tb[3,2], tb[3,1]-1.96*tb[3,2])*IQRs[,beta])
  RR1A_tb <- rbind(RR1A_tb, estimates)
  rownames_list <- c(rownames_list, paste0(deathc_i, "age"))
}
rownames(RR1A_tb) <- rownames_list
colnames(RR1A_tb) <- c("1ARR", "1ARR_lci", "1ARR_uci")
print(RR1A_tb)

results_A <- cbind(RR0A_tb, RR1A_tb)
print(results_A)

write.table(results_A, file = paste0(dir_results, "RRiqr_pm25_betaPM.csv"))

## Exposure B (only beta): RR of beta for one IQR increase in two models ------
## model 0B
deathc <- c("TOT", "CVD", "MI", "stroke", 
            paste0(c("TOT", "CVD", "MI", "stroke"), 1865), 
            paste0(c("TOT", "CVD", "MI", "stroke"), 6585), 
            paste0(c("TOT", "CVD", "MI", "stroke"), 85))
RR0B_tb <- NULL
rownames_list <- NULL
for (deathc_i in deathc) {
  tb <- read.table(paste0(dir_modresults, "mod0B_", deathc_i, ".csv"))
  estimates <- exp(c(tb[2,1], tb[2,1]-1.96*tb[2,2], tb[2,1]-1.96*tb[2,2])*IQRs[,beta])
  RR0B_tb <- rbind(RR0B_tb, estimates)
  rownames_list <- c(rownames_list, paste0(deathc_i, "age"))
}
rownames(RR0B_tb) <- rownames_list
colnames(RR0B_tb) <- c("0BRR", "0BRR_lci", "0BRR_uci")
print(RR0B_tb)

## model 1B
deathc <- c("TOT", "CVD", "MI", "stroke", 
            paste0(c("TOT", "CVD", "MI", "stroke"), 1865), 
            paste0(c("TOT", "CVD", "MI", "stroke"), 6585), 
            paste0(c("TOT", "CVD", "MI", "stroke"), 85))
RR1B_tb <- NULL
rownames_list <- NULL
for (deathc_i in deathc) {
  tb <- read.table(paste0(dir_modresults, "mod1B_", deathc_i, ".csv"))
  estimates <- exp(c(tb[2,1], tb[2,1]-1.96*tb[2,2], tb[2,1]-1.96*tb[2,2])*IQRs[,beta])
  RR1B_tb <- rbind(RR1B_tb, estimates)
  rownames_list <- c(rownames_list, paste0(deathc_i, "age"))
}
rownames(RR1B_tb) <- rownames_list
colnames(RR1B_tb) <- c("1BRR", "1BRR_lci", "1BRR_uci")
print(RR1B_tb)

results_B <- cbind(RR0B_tb, RR1B_tb)
print(results_B)

write.table(results_B, file = paste0(dir_results, "RRiqr_beta_beta.csv"))
# kable(results_B, caption = "Estimated RR for one IQR increase in beta radiation, with Exposure Set B (only radiation)")

## 3. Results Plots ---------------------------------------------------------
library(ggplot2)
dir_results <- "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/betaRadiation_CVD/results/"
## 3.1 Plot for beta radiation --------------
results_A <- read.table(file = paste0(dir_results, "RRiqr_beta_betaPM.csv"))
results_B <- read.table(file = paste0(dir_results, "RRiqr_beta_beta.csv"))

RR_A <- c(results_A[,1], results_A[,4])
lowCI_A <- c(results_A[,2], results_A[,5])
highCI_A <- c(results_A[,3], results_A[,6])

RR_B <- c(results_B[,1], results_B[,4])
lowCI_B <- c(results_B[,2], results_B[,5])
highCI_B <- c(results_B[,3], results_B[,6])

n_model <- length(RR_A)/2
results_beta_all <- rbind(cbind(RR_A, lowCI_A, highCI_A), cbind(RR_B, lowCI_B, highCI_B))
print(results_beta_all)
results_beta_all <- data.frame(results_beta_all)
View(results_beta_all)
dim(results_beta_all)

results_beta_all$mod <- rep(c(rep("Differences in differences", n_model), rep("Mixed-effect", n_model)),2)
results_beta_all$cause <- rep(c("TOT", "CVD", "MI", "stroke"), dim(results_beta_all)[1]/4)
results_beta_all$exposures <-  c(rep("beta radiation + PM[2.5]", n_model*2), rep("only beta radiation", n_model*2))
results_beta_all$age_group <- rep(c(rep("18+", 4), rep("18-65", 4),rep("65-85", 4),rep("85+", 4)),2)
colnames(results_beta_all) <- c("RR", "lowCI", "highCI", "mod", "cause", "exposures", "age_group")
View(results_beta_all)

write.csv(results_beta_all, paste0(dir_results, "results_beta_all_plot.csv"))

# ggplot(results_beta_all, aes(x = cause, y = RR)) +
#   geom_errorbar(aes(ymin = lowCI, ymax = highCI, color = mod, position = position_dodge(0.3), width = 0.3)) +
#   geom_point(aes(shape = age_group, color = mod, linetype = exposures), position = position_dodge(0.3)) +
#   geom_hline(yintercept = 1, linetype="dashed", color = 1, size = 0.5) +
#   xlab("Death Cause") + ylab("Rate Ratio for an IQR increase \n with 95% Confidence Interval") +
#   theme(legend.position="bottom") +
#   labs(color = "Models") +
#   labs(linetype = "Exposure sets") +
#   guides(color=guide_legend(nrow=2, byrow=TRUE), linetype=guide_legend(nrow=2, byrow=TRUE), shape=guide_legend(nrow=2, byrow=TRUE))
# 
# pdf(paste0(dir_results, "RRiqr_beta.pdf"), height = 5)
# ggplot(results_beta_all, aes(x = cause, y = RR)) +
#   geom_errorbar(aes(ymin = lowCI, ymax = highCI, color = mod, linetype = exposures),position = position_dodge(0.3), width = 0.3) +
#   geom_point(aes(color = mod, linetype = exposures), position = position_dodge(0.3)) + 
#   geom_hline(yintercept = 1, linetype="dashed", color = 1, size = 0.5) + 
#   xlab("Death Cause") + ylab("Rate Ratio for an IQR increase \n with 95% Confidence Interval") +
#   theme(legend.position="bottom") +
#   labs(color = "Models") + 
#   labs(linetype = "Exposure sets") +
#   guides(color=guide_legend(nrow=2, byrow=F), linetype=guide_legend(nrow=2, byrow=F))
# dev.off()

# 3.2 Plot for PM25 --------------
results_A <- read.table(file = paste0(dir_results, "RRiqr_pm25_betaPM.csv"))

RR_A <- c(results_A[,1], results_A[,4])
lowCI_A <- c(results_A[,2], results_A[,5])
highCI_A <- c(results_A[,3], results_A[,6])

results_pm_all <- cbind(RR_A, lowCI_A, highCI_A)
results_pm_all <- data.frame(results_pm_all)
View(results_pm_all)

results_pm_all$mod <- c(rep("Differences in differences",16), rep("Mixed-effect",16))
results_pm_all$cause <- rep(c("TOT", "CVD", "MI", "stroke"),8)
colnames(results_pm_all) <- c("RR", "lowCI", "highCI", "mod", "cause")

write.csv(results_pm_all, paste0(dir_results, "results_pm_all_plot.csv"))
# 
# ggplot(results_pm_all, aes(x = cause, y = RR)) +
#   geom_errorbar(aes(ymin = lowCI, ymax = highCI, color = mod),position = position_dodge(0.3), width = 0.3) +
#   geom_point(aes(color = mod), position = position_dodge(0.3)) + 
#   geom_hline(yintercept = 1, linetype="dashed", color = 1, size = 0.5) + 
#   xlab("Death Cause") + ylab("Rate Ratio for an IQR increase \n with 95% Confidence Interval") +
#   theme(legend.position="bottom") +
#   labs(color = "Models")
# # ggtitle("Rate ratio for an IQR increase in PM[2.5] in different models")
# 
# pdf(paste0(dir_results, "RRiqr_pm25.pdf"), height = 5)
# ggplot(results_pm_all, aes(x = cause, y = RR)) +
#   geom_errorbar(aes(ymin = lowCI, ymax = highCI, color = mod),position = position_dodge(0.3), width = 0.3) +
#   geom_point(aes(color = mod), position = position_dodge(0.3)) + 
#   geom_hline(yintercept = 1, linetype="dashed", color = 1, size = 0.5) + 
#   xlab("Death Cause") + ylab("Rate Ratio for an IQR increase \n with 95% Confidence Interval") +
#   theme(legend.position="bottom") +
#   labs(color = "Models")
# dev.off()