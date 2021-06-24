###############################################################################
# Project: beta radiation and CVD death in MA
# Code: Output descriptive stats
# Input: "finalDT.rds"
# Output: "table1.doc" "corrTable.pdf"
# Author: Shuxin Dong                                                         
###############################################################################

## 0. set up ------------------------------------------------------------------
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

## 1. table one  ------------------------------------------------------------------
summary(dt)
listVars <- c("CVD", "MI", "stroke", "TOT", "Beta", "pm25", "summer_tmean", "winter_tmean", "medhouseholdincome", "medianhousevalue","hispanic", "pct_blk", "poverty", "education", "popdensity", "mean_bmi", "smoke_rate", "pcount")
rawtable1 <- CreateTableOne(vars = listVars, data = dt)
table1 <- print(rawtable1, 
                formatOption = list(decimal.mark = ".",  big.mark = ",", scientific = FALSE),
                contDigits = 3)
rtffile <- RTF(file = paste0(dir_results, "table1.doc"))  # this can be an .rtf or a .doc
addParagraph(rtffile, "Table")
addTable(rtffile, cbind(rownames(table1), table1))
done(rtffile)

## 2. correlation table -------------------------------------------------------
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
corrplot(M, method="number", type = "lower", p.mat = p.mat, sig.level = 0.01)

pdf(paste0(dir_results,"corrTable.pdf"), width = 10, height = 10)
corrplot(M, method="number", type = "lower", p.mat = p.mat, sig.level = 0.01)
dev.off()
