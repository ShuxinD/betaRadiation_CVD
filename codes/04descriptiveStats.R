#' Project: beta radiation and CVD death in MA
#' Code: Output descriptive stats
#' Input: "finalDT.rds"
#' Output: "table1.doc" "corrTable.pdf"
#' Author: Shuxin Dong                                                       

## 0. set up ---------------------------------------------------------------
rm(list = ls())
gc()

setwd("/media/qnap3/Shuxin/ParticalRadiation_MAdeath/")
dir_results <- "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/betaRadiation_CVD/results/"

library(data.table)
library(ggplot2)
library(tableone)
library(rtf)
library(corrplot)

dt <- readRDS("finalDT.rds")

IQR(dt$Beta)
# 0.05530013

## 1. table one  -----------------------------------------------------------
summary(dt)
length(unique(dt[,ZCTA5CE10]))

listVars <- c("CVD", "MI", "stroke", "TOT",
              "CVD1865", "MI1865", "stroke1865", "TOT1865",
              "CVD6585", "MI6585", "stroke6585", "TOT6585",
              "CVD85", "MI85", "stroke85", "TOT85",
              "Beta", "pm25", "summer_tmean", "winter_tmean", "medhouseholdincome", "medianhousevalue","hispanic", "pct_blk", "poverty", "education", "popdensity", "mean_bmi", "smoke_rate", "pcount")
rawtable1 <- CreateTableOne(vars = listVars, data = dt)
table1 <- print(rawtable1, 
                formatOption = list(decimal.mark = ".",  big.mark = ",", scientific = FALSE),
                contDigits = 3)
rtffile <- RTF(file = paste0(dir_results, "table1.doc"))  # this can be an .rtf or a .doc
addParagraph(rtffile, "Table")
addTable(rtffile, cbind(rownames(table1), table1))
done(rtffile)

dt[,`:=`(CHF=NULL,CHF1865=NULL,CHF6585=NULL,CHF85=NULL,
         ZCTA5CE10=NULL)]
dt[,`:=`(CVDrate=CVD/pcount,MIrate=MI/pcount,strokerate=stroke/pcount,TOTrate=TOT/pcount)][]
summary <- rbind(dt[, lapply(.SD, mean)], dt[, lapply(.SD, sd)],
      dt[, lapply(.SD, quantile)])
rownames(summary) <- c("mean", "sd", "quantile0", "quantile25", "quantile50","quantile75", "quantile1")
write.csv(summary, paste0(dir_results, "tableone_spatial.csv"))

dt[,.(CVD_mean=mean(CVD)), by=year]

## 2. correlation table ----------------------------------------------------
M <- cor(na.omit(dt)[,.(Beta, pm25, summer_tmean, winter_tmean, medhouseholdincome, medianhousevalue, hispanic, pct_blk, poverty, education, popdensity, mean_bmi, smoke_rate)])
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
p.mat <- cor.mtest(na.omit(dt)[,.(Beta, pm25, summer_tmean, winter_tmean, medhouseholdincome, medianhousevalue, hispanic, pct_blk, poverty, education, popdensity, mean_bmi, smoke_rate)])
colnames(M) <- c("Particle radiation", "PM[2.5]", "Summer average temperature", "Winter average temperature", "Median household income", "Median value of house", "Percentage of Hispanic", "Percentage of black", "Pencentage below poverty line", "Percentage without high school diploma", "Population density", "Mean BMI", "Smoking rate")
rownames(M) <- c("Particle radiation", "PM[2.5]", "Summer average temperature", "Winter average temperature", "Median household income", "Median value of house", "Percentage of Hispanic", "Percentage of black", "Pencentage below poverty line", "Percentage without high school diploma", "Population density", "Mean BMI", "Smoking rate")
par(mfrow=c(1,1))
corrplot(M, method="number", type = "lower", p.mat = p.mat, sig.level = 0.05)

pdf(paste0(dir_results,"corrTable.pdf"), width = 10, height = 10)
corrplot(M, method="number", type = "lower", p.mat = p.mat, sig.level = 0.05)
dev.off()
