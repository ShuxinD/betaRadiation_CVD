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
dir_results <- "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/betaRadiation_CVD/results/"

library(mgcv)

dt <- readRDS("finalDT.rds")
dt <- na.omit(dt)
dt_pcount <- dt[,pcount]

## 1. models ------------------------------------------------------------------
## 1.1 DID ------------------------------------------------------------------
## beta + PM25
mod0A_TOT <- gam(TOT ~ Beta + pm25 + summer_tmean + winter_tmean + as.factor(year) + as.factor(ZCTA5CE10) + offset(log(dt_pcount)),
                 data = dt,
                 family = quasipoisson(link = "log"))
mod0A_CVD <- gam(CVD ~ Beta + pm25 + summer_tmean + winter_tmean + as.factor(year) + as.factor(ZCTA5CE10) + offset(log(dt_pcount)),
                 data = dt,
                 family = quasipoisson(link = "log"))
mod0A_MI <- gam(MI ~ Beta + pm25 + summer_tmean + winter_tmean + as.factor(year) + as.factor(ZCTA5CE10) + offset(log(dt_pcount)),
                data = dt,
                family = quasipoisson(link = "log"))
mod0A_stroke <- gam(stroke ~ Beta + pm25 + summer_tmean + winter_tmean + as.factor(year) + as.factor(ZCTA5CE10) + offset(log(dt_pcount)),
                    data = dt,
                    family = quasipoisson(link = "log"))
## only beta
mod0B_TOT <- gam(TOT ~ Beta + summer_tmean + winter_tmean + as.factor(year) + as.factor(ZCTA5CE10) + offset(log(dt_pcount)),
                 data = dt,
                 family = quasipoisson(link = "log"))
mod0B_CVD <- gam(CVD ~ Beta + summer_tmean + winter_tmean + as.factor(year) + as.factor(ZCTA5CE10) + offset(log(dt_pcount)),
                 data = dt,
                 family = quasipoisson(link = "log"))
mod0B_MI <- gam(MI ~ Beta + summer_tmean + winter_tmean + as.factor(year) + as.factor(ZCTA5CE10) + offset(log(dt_pcount)),
                data = dt,
                family = quasipoisson(link = "log"))
mod0B_stroke <- gam(stroke ~ Beta + summer_tmean + winter_tmean + as.factor(year) + as.factor(ZCTA5CE10) + offset(log(dt_pcount)),
                    data = dt,
                    family = quasipoisson(link = "log"))
## 1.2 mixed-effect model -----------------------------------------------------
## beta and PM25
mod1A_TOT <- gamm(TOT ~ Beta + pm25 + summer_tmean + winter_tmean + medhouseholdincome + medianhousevalue + hispanic + pct_blk + poverty + education + popdensity + mean_bmi + smoke_rate + as.factor(year) + offset(log(dt_pcount)),
                  data = dt,
                  random = list(ZCTA5CE10=~1),
                  family = quasipoisson(link = "log"))
mod1A_CVD <- gamm(CVD ~ Beta + pm25 + summer_tmean + winter_tmean + medhouseholdincome + medianhousevalue + hispanic + pct_blk + poverty + education + popdensity + mean_bmi + smoke_rate + as.factor(year) + offset(log(dt_pcount)),
                  data = dt,
                  random = list(ZCTA5CE10=~1),
                  family = quasipoisson(link = "log"))
mod1A_MI <- gamm(MI ~ Beta + pm25 + summer_tmean + winter_tmean + medhouseholdincome + medianhousevalue + hispanic + pct_blk + poverty + education + popdensity + mean_bmi + smoke_rate + as.factor(year) + offset(log(dt_pcount)),
                 data = dt,
                 random = list(ZCTA5CE10=~1),
                 family = quasipoisson(link = "log"))
mod1A_stroke <- gamm(stroke ~ Beta + pm25 + summer_tmean + winter_tmean + medhouseholdincome + medianhousevalue + hispanic + pct_blk + poverty + education + popdensity + mean_bmi + smoke_rate + as.factor(year) + offset(log(dt_pcount)),
                     data = dt,
                     random = list(ZCTA5CE10=~1),
                     family = quasipoisson(link = "log"))
## only beta
mod1B_TOT <- gamm(TOT ~ Beta + summer_tmean + winter_tmean + medhouseholdincome + medianhousevalue + hispanic + pct_blk + poverty + education + popdensity + mean_bmi + smoke_rate + as.factor(year) + offset(log(dt_pcount)),
                  data = dt,
                  random = list(ZCTA5CE10=~1),
                  family = quasipoisson(link = "log"))
mod1B_CVD <- gamm(CVD ~ Beta + summer_tmean + winter_tmean + medhouseholdincome + medianhousevalue + hispanic + pct_blk + poverty + education + popdensity + mean_bmi + smoke_rate + as.factor(year) + offset(log(dt_pcount)),
                  data = dt,
                  random = list(ZCTA5CE10=~1),
                  family = quasipoisson(link = "log"))
mod1B_MI <- gamm(MI ~ Beta + summer_tmean + winter_tmean + medhouseholdincome + medianhousevalue + hispanic + pct_blk + poverty + education + popdensity + mean_bmi + smoke_rate + as.factor(year) + offset(log(dt_pcount)),
                 data = dt,
                 random = list(ZCTA5CE10=~1),
                 family = quasipoisson(link = "log"))
mod1B_stroke <- gamm(stroke ~ Beta + summer_tmean + winter_tmean + medhouseholdincome + medianhousevalue + hispanic + pct_blk + poverty + education + popdensity + mean_bmi + smoke_rate + as.factor(year) + offset(log(dt_pcount)),
                     data = dt,
                     random = list(ZCTA5CE10=~1),
                     family = quasipoisson(link = "log"))

## 2. Results tables ---------------------------------------------------------
## IQR ----
iqr_beta <- IQR(dt[,Beta])
iqr_pm25 <- IQR(dt[,pm25])
cat("IQR of annual beta is:", iqr_beta, "\n")
# IQR of annual beta is: 0.05529197 
cat("IQR of annual pm25 is:", iqr_pm25)
# IQR of annual pm25 is: 2.81567

## RR for one IQR increase in two models with Exposure set A (beta + PM25) ----
TOT_A <- c(exp(summary(mod0A_TOT)$p.table[2,1]*iqr_beta), exp((summary(mod0A_TOT)$p.table[2,1]-1.96*summary(mod0A_TOT)$p.table[2,2])*iqr_beta), exp((summary(mod0A_TOT)$p.table[2,1]+1.96*summary(mod0A_TOT)$p.table[2,2])*iqr_beta), exp(summary(mod1A_TOT$gam)$p.table[2,1]*iqr_beta), exp((summary(mod1A_TOT$gam)$p.table[2,1]-1.96*summary(mod1A_TOT$gam)$p.table[2,2])*iqr_beta), exp((summary(mod1A_TOT$gam)$p.table[2,1]+1.96*summary(mod1A_TOT$gam)$p.table[2,2])*iqr_beta))
CVD_A <- c(exp(summary(mod0A_CVD)$p.table[2,1]*iqr_beta), exp((summary(mod0A_CVD)$p.table[2,1]-1.96*summary(mod0A_CVD)$p.table[2,2])*iqr_beta), exp((summary(mod0A_CVD)$p.table[2,1]+1.96*summary(mod0A_CVD)$p.table[2,2])*iqr_beta), exp(summary(mod1A_CVD$gam)$p.table[2,1]*iqr_beta), exp((summary(mod1A_CVD$gam)$p.table[2,1]-1.96*summary(mod1A_CVD$gam)$p.table[2,2])*iqr_beta), exp((summary(mod1A_CVD$gam)$p.table[2,1]+1.96*summary(mod1A_CVD$gam)$p.table[2,2])*iqr_beta))
MI_A <- c(exp(summary(mod0A_MI)$p.table[2,1]*iqr_beta), exp((summary(mod0A_MI)$p.table[2,1]-1.96*summary(mod0A_MI)$p.table[2,2])*iqr_beta), exp((summary(mod0A_MI)$p.table[2,1]+1.96*summary(mod0A_MI)$p.table[2,2])*iqr_beta), exp(summary(mod1A_MI$gam)$p.table[2,1]*iqr_beta), exp((summary(mod1A_MI$gam)$p.table[2,1]-1.96*summary(mod1A_MI$gam)$p.table[2,2])*iqr_beta), exp((summary(mod1A_MI$gam)$p.table[2,1]+1.96*summary(mod1A_MI$gam)$p.table[2,2])*iqr_beta))
stroke_A <- c(exp(summary(mod0A_stroke)$p.table[2,1]*iqr_beta), exp((summary(mod0A_stroke)$p.table[2,1]-1.96*summary(mod0A_stroke)$p.table[2,2])*iqr_beta), exp((summary(mod0A_stroke)$p.table[2,1]+1.96*summary(mod0A_stroke)$p.table[2,2])*iqr_beta), exp(summary(mod1A_stroke$gam)$p.table[2,1]*iqr_beta), exp((summary(mod1A_stroke$gam)$p.table[2,1]-1.96*summary(mod1A_stroke$gam)$p.table[2,2])*iqr_beta), exp((summary(mod1A_stroke$gam)$p.table[2,1]+1.96*summary(mod1A_stroke$gam)$p.table[2,2])*iqr_beta))
results_A <- rbind(TOT_A, CVD_A, MI_A, stroke_A)
results_A <- data.frame(results_A)
colnames(results_A) <- c("Model 0 RR", "Model 0 95% CI (lowest bound)", "Model 0 95% CI (highest bound)", "Model 1 RR", "Model 1 95% CI (lowest bound)", "Model 1 95% CI (highest bound)")
print(results_A)
write.csv(results_A, file = paste0(dir_results, "RRiqr_beta_betaPM.csv"))
# kable(results_A, caption = "Estimated RR for one IQR increase in beta radiation, with Exposure Set A (radiation + PM2.5)")

## RR for one IQR increase in two models with Exposure set B (only beta) ------
TOT_B <- c(exp(summary(mod0B_TOT)$p.table[2,1]*iqr_beta), exp((summary(mod0B_TOT)$p.table[2,1]-1.96*summary(mod0B_TOT)$p.table[2,2])*iqr_beta), exp((summary(mod0B_TOT)$p.table[2,1]+1.96*summary(mod0B_TOT)$p.table[2,2])*iqr_beta), exp(summary(mod1B_TOT$gam)$p.table[2,1]*iqr_beta), exp((summary(mod1B_TOT$gam)$p.table[2,1]-1.96*summary(mod1B_TOT$gam)$p.table[2,2])*iqr_beta), exp((summary(mod1B_TOT$gam)$p.table[2,1]+1.96*summary(mod1B_TOT$gam)$p.table[2,2])*iqr_beta))
CVD_B <- c(exp(summary(mod0B_CVD)$p.table[2,1]*iqr_beta), exp((summary(mod0B_CVD)$p.table[2,1]-1.96*summary(mod0B_CVD)$p.table[2,2])*iqr_beta), exp((summary(mod0B_CVD)$p.table[2,1]+1.96*summary(mod0B_CVD)$p.table[2,2])*iqr_beta), exp(summary(mod1B_CVD$gam)$p.table[2,1]*iqr_beta), exp((summary(mod1B_CVD$gam)$p.table[2,1]-1.96*summary(mod1B_CVD$gam)$p.table[2,2])*iqr_beta), exp((summary(mod1B_CVD$gam)$p.table[2,1]+1.96*summary(mod1B_CVD$gam)$p.table[2,2])*iqr_beta))
MI_B <- c(exp(summary(mod0B_MI)$p.table[2,1]*iqr_beta), exp((summary(mod0B_MI)$p.table[2,1]-1.96*summary(mod0B_MI)$p.table[2,2])*iqr_beta), exp((summary(mod0B_MI)$p.table[2,1]+1.96*summary(mod0B_MI)$p.table[2,2])*iqr_beta), exp(summary(mod1B_MI$gam)$p.table[2,1]*iqr_beta), exp((summary(mod1B_MI$gam)$p.table[2,1]-1.96*summary(mod1B_MI$gam)$p.table[2,2])*iqr_beta), exp((summary(mod1B_MI$gam)$p.table[2,1]+1.96*summary(mod1B_MI$gam)$p.table[2,2])*iqr_beta))
stroke_B <- c(exp(summary(mod0B_stroke)$p.table[2,1]*iqr_beta), exp((summary(mod0B_stroke)$p.table[2,1]-1.96*summary(mod0B_stroke)$p.table[2,2])*iqr_beta), exp((summary(mod0B_stroke)$p.table[2,1]+1.96*summary(mod0B_stroke)$p.table[2,2])*iqr_beta), exp(summary(mod1B_stroke$gam)$p.table[2,1]*iqr_beta), exp((summary(mod1B_stroke$gam)$p.table[2,1]-1.96*summary(mod1B_stroke$gam)$p.table[2,2])*iqr_beta), exp((summary(mod1B_stroke$gam)$p.table[2,1]+1.96*summary(mod1B_stroke$gam)$p.table[2,2])*iqr_beta))
results_B <- rbind(TOT_B, CVD_B, MI_B, stroke_B)
results_B <- data.frame(results_B)
colnames(results_B) <- c("Model 0 RR", "Model 0 95% CI (lowest bound)", "Model 0 95% CI (highest bound)", "Model 1 RR", "Model 1 95% CI (lowest bound)", "Model 1 95% CI (highest bound)")
print(results_B)
write.csv(results_B, file = paste0(dir_results, "RRiqr_beta_beta.csv"))
# kable(results_B, caption = "Estimated RR for one IQR increase in beta radiation, with Exposure Set B (only radiation)")

## 3. Results Plots ---------------------------------------------------------
## 3.1 Plot for beta radiation --------------
RR_A <- c(results_A[,1], results_A[,4])
lowCI_A <- c(results_A[,2], results_A[,5])
highCI_A <- c(results_A[,3], results_A[,6])

RR_B <- c(results_B[,1], results_B[,4])
lowCI_B <- c(results_B[,2], results_B[,5])
highCI_B <- c(results_B[,3], results_B[,6])

results_beta_all <- rbind(cbind(RR_A, lowCI_A, highCI_A), cbind(RR_B, lowCI_B, highCI_B))
results_beta_all <- data.frame(results_beta_all)

results_beta_all$mod <- rep(c(rep("Differences in differences",4), rep("Mixed-effect",4)),2)
results_beta_all$cause <- rep(c("TOT", "CVD", "MI", "stroke"),4)
results_beta_all$exposures <-  c(rep("beta radiation + PM[2.5]",8), rep("only beta radiation",8))
colnames(results_beta_all) <- c("RR", "lowCI", "highCI", "mod", "cause", "exposures")

ggplot(results_beta_all, aes(x = cause, y = RR)) +
  geom_errorbar(aes(ymin = lowCI, ymax = highCI, color = mod, linetype = exposures),position = position_dodge(0.3), width = 0.3) +
  geom_point(aes(color = mod, linetype = exposures), position = position_dodge(0.3)) + 
  geom_hline(yintercept = 1, linetype="dashed", color = 1, size = 0.5) + 
  xlab("Death Cause") + ylab("Rate Ratio for an IQR increase \n with 95% Confidence Interval") +
  theme(legend.position="bottom") +
  labs(color = "Models") + 
  labs(linetype = "Exposure sets") +
  guides(color=guide_legend(nrow=2, byrow=TRUE), linetype=guide_legend(nrow=2, byrow=TRUE))

pdf(paste0(dir_results, "RRiqr_beta.pdf"), height = 5)
ggplot(results_beta_all, aes(x = cause, y = RR)) +
  geom_errorbar(aes(ymin = lowCI, ymax = highCI, color = mod, linetype = exposures),position = position_dodge(0.3), width = 0.3) +
  geom_point(aes(color = mod, linetype = exposures), position = position_dodge(0.3)) + 
  geom_hline(yintercept = 1, linetype="dashed", color = 1, size = 0.5) + 
  xlab("Death Cause") + ylab("Rate Ratio for an IQR increase \n with 95% Confidence Interval") +
  theme(legend.position="bottom") +
  labs(color = "Models") + 
  labs(linetype = "Exposure sets") +
  guides(color=guide_legend(nrow=2, byrow=F), linetype=guide_legend(nrow=2, byrow=F))
dev.off()

## 3.2 Plot for PM25 --------------
RR_A <- c(results_A[,1], results_A[,4])
lowCI_A <- c(results_A[,2], results_A[,5])
highCI_A <- c(results_A[,3], results_A[,6])

results_pm_all <- cbind(RR_A, lowCI_A, highCI_A)
results_pm_all <- data.frame(results_pm_all)
results_pm_all$mod <- c(rep("Differences in differences",4), rep("Mixed-effect",4))
results_pm_all$cause <- rep(c("TOT", "CVD", "MI", "stroke"),2)
colnames(results_pm_all) <- c("RR", "lowCI", "highCI", "mod", "cause")

ggplot(results_pm_all, aes(x = cause, y = RR)) +
  geom_errorbar(aes(ymin = lowCI, ymax = highCI, color = mod),position = position_dodge(0.3), width = 0.3) +
  geom_point(aes(color = mod), position = position_dodge(0.3)) + 
  geom_hline(yintercept = 1, linetype="dashed", color = 1, size = 0.5) + 
  xlab("Death Cause") + ylab("Rate Ratio for an IQR increase \n with 95% Confidence Interval") +
  theme(legend.position="bottom") +
  labs(color = "Models")
# ggtitle("Rate ratio for an IQR increase in PM[2.5] in different models")

pdf(paste0(dir_results, "RRiqr_pm25.pdf"), height = 5)
ggplot(results_pm_all, aes(x = cause, y = RR)) +
  geom_errorbar(aes(ymin = lowCI, ymax = highCI, color = mod),position = position_dodge(0.3), width = 0.3) +
  geom_point(aes(color = mod), position = position_dodge(0.3)) + 
  geom_hline(yintercept = 1, linetype="dashed", color = 1, size = 0.5) + 
  xlab("Death Cause") + ylab("Rate Ratio for an IQR increase \n with 95% Confidence Interval") +
  theme(legend.position="bottom") +
  labs(color = "Models")
dev.off()