#' Project: beta radiation and CVD death in MA
#' Code: plot results
#' Input: "RRiqr_pm25_PM.csvf"
#' Input: "RRiqr_beta_beta.csv"
#' Input: "RRiqr_beta_betaPM.csv"
#' Input: "RRiqr_pm25_betaPM.csv"
#' Output: ...
#' Author: Shuxin Dong                                                         

## 0. set up ----
rm(list = ls())
gc()

setwd("/media/qnap3/Shuxin/ParticalRadiation_MAdeath/")
library(ggplot2)
library(ggpubr)
library(data.table)

## 1. prepare the dataset for plotting -----
dir_results.table <- "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/github_repo/results/main_PRPM/supplements/"
dir_results.forplot <- "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/github_repo/results/main_PRPM/"

## 1.1 beta radiation -------
results_A <- read.table(file = paste0(dir_results.table, "RRiqr_beta_betaPM.csv"))
results_B <- read.table(file = paste0(dir_results.table, "RRiqr_beta_beta.csv"))

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
# View(results_beta_all)
dim(results_beta_all)

results_beta_all$mod <- rep(c(rep("Difference in differences", n_model), rep("Generalized linear mixed model", n_model)),2)
results_beta_all$cause <- rep(c("TOT", "CVD", "MI", "stroke"), dim(results_beta_all)[1]/4)
results_beta_all$exposures <-  c(rep("gross beta-activity + PM[2.5]", n_model*2), rep("only gross beta-activity", n_model*2))
results_beta_all$age_group <- rep(c(rep("18+", 4), rep("18-65", 4),rep("65-85", 4),rep("85+", 4)),2)
colnames(results_beta_all) <- c("RR", "lowCI", "highCI", "mod", "cause", "exposures", "age_group")
# View(results_beta_all)

write.csv(results_beta_all, paste0(dir_results.forplot, "results_beta_all_plot.csv"))

## 1.2 PM2.5 -------
results_A <- read.table(file = paste0(dir_results.table, "RRiqr_pm25_betaPM.csv"))
results_B <- read.table(file = paste0(dir_results.table, "RRiqr_pm25_PM.csv"))

RR_A <- c(results_A[,1], results_A[,4])
lowCI_A <- c(results_A[,2], results_A[,5])
highCI_A <- c(results_A[,3], results_A[,6])

RR_B <- c(results_B[,1], results_B[,4])
lowCI_B <- c(results_B[,2], results_B[,5])
highCI_B <- c(results_B[,3], results_B[,6])

n_model <- length(RR_A)/2
results_PM_all <- rbind(cbind(RR_A, lowCI_A, highCI_A), cbind(RR_B, lowCI_B, highCI_B))
print(results_PM_all)
results_PM_all <- data.frame(results_PM_all)
# View(results_PM_all)
dim(results_PM_all)

results_PM_all$mod <- rep(c(rep("Difference in differences", n_model), rep("Generalized linear mixed model", n_model)),2)
results_PM_all$cause <- rep(c("TOT", "CVD", "MI", "stroke"), dim(results_PM_all)[1]/4)
results_PM_all$exposures <-  c(rep("gross beta-activity + PM[2.5]", n_model*2), rep("only PM[2.5]", n_model*2))
results_PM_all$age_group <- rep(c(rep("18+", 4), rep("18-65", 4),rep("65-85", 4),rep("85+", 4)),2)
colnames(results_PM_all) <- c("RR", "lowCI", "highCI", "mod", "cause", "exposures", "age_group")
# View(results_PM_all)

write.csv(results_PM_all, paste0(dir_results.forplot, "results_PM_all_plot.csv"))

## 2. plot beta radiation -----------------------------------------------------
dir_plot <- "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/github_repo/results/main_PRPM/"
## 2.1 beta radiation with age group ----
## all exposure sets
plotDT <- results_beta_all
setDT(plotDT)
plotDT[, age_group:= factor(age_group, levels = c("18+", "18-65","65-85", "85+"))]
plotbeta <- ggplot(plotDT, aes(x = cause, y = RR)) +
  geom_pointrange(size=0.3, aes(ymin = lowCI, ymax = highCI, shape = age_group, linetype = exposures, color = mod), position = position_dodge(0.8)) +
  geom_hline(yintercept = 1, linetype="dashed", color = 1, size = 0.2) +
  ylab("Rate ratio for an IQR increase\nwith 95% confidence interval") + xlab("Death cause") +
  labs(color = "Models") +
  labs(linetype = "Exposure sets") +
  labs(shape = "Age groups") +
  guides(color=guide_legend(nrow=2, override.aes=list(shape=c(NA,NA))), linetype=guide_legend(nrow=2, override.aes=list(shape=c(NA,NA))), shape=guide_legend(nrow=2, override.aes=list(linetype=c(0,0)))) +
  scale_x_discrete(labels=c("CVD" = "Cardiovascular\ndisease", "MI" = "Myocardial\ninfarction", "TOT" = "All-causes")) +
  theme_minimal()
plotbeta

cairo_pdf(paste0(dir_plot, "RRiqr_beta_age.pdf"), height = 3.5)
plotbeta
dev.off()

## exposure set as beta+pm25
plotDT <- results_beta_all
setDT(plotDT)
plotDT <- plotDT[exposures=="gross beta-activity + PM[2.5]"]
plotDT[, age_group:= factor(age_group, levels = c("18+", "18-65","65-85", "85+"))]
plotbeta <- ggplot(plotDT, aes(x = cause, y = RR)) +
  geom_pointrange(size=0.5, aes(ymin = lowCI, ymax = highCI, shape = age_group, color = mod), position = position_dodge(0.8)) +
  geom_hline(yintercept = 1, linetype="dashed", color = 1, size = 0.2) +
  ylab("Rate ratio") + xlab("Death cause") +
  labs(color = "Models") +
  labs(shape = "Age groups") +
  guides(color=guide_legend(nrow=2, override.aes=list(shape=c(NA,NA))), shape=guide_legend(nrow=2, override.aes=list(linetype=c(0,0)))) +
  scale_x_discrete(labels=c("CVD" = "Cardiovascular\ndisease", "MI" = "Myocardial\ninfarction", "TOT" = "Non-accidental\n all causes")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_rect(colour = "black", fill=NA, size=1))
plotbeta

cairo_pdf(paste0(dir_plot, "RRiqr_beta_age_betaPM.pdf"), height = 3.5)
plotbeta
dev.off()

## 2.2 beta radiation only all ages ----
plotDT <- results_beta_all
setDT(plotDT)
plotDT <- plotDT[age_group=="18+",]
plotbeta_main <- ggplot(plotDT, aes(x = cause, y = RR)) +
  geom_pointrange(size=0.5, aes(ymin = lowCI, ymax = highCI,linetype = exposures, color = mod), position = position_dodge(0.8)) +
  geom_hline(yintercept = 1, linetype="dashed", color = 1, size = 0.2) +
  ylab("Rate ratio") + xlab("Death cause") +
  labs(color = "Models") +
  labs(linetype = "Exposure sets") +
  guides(color=guide_legend(nrow=2, override.aes=list(shape=c(NA,NA))), shape=guide_legend(nrow=2, override.aes=list(linetype=c(0,0)))) +
  scale_x_discrete(labels=c("CVD" = "Cardiovascular\ndisease", "MI" = "Myocardial\ninfarction", "stroke" = "Stroke","TOT" = "Non-accidental\nall causes")) +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_rect(colour = "black", fill=NA, size=1))
plotbeta_main

cairo_pdf(paste0(dir_plot, "RRiqr_beta_main.pdf"), height = 3.5)
plotbeta_main
dev.off()

# 3. plot PM2.5 --------------
dir_plot <- "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/github_repo/results/main_PRPM/"
## 3.1 pm25 with age group ----
## all exposure sets
plotDT <- results_PM_all
setDT(plotDT)
plotDT[, age_group:= factor(age_group, levels = c("18+", "18-65","65-85", "85+"))]
plotpm <- ggplot(plotDT, aes(x = cause, y = RR)) +
  geom_pointrange(size=0.3, aes(ymin = lowCI, ymax = highCI, shape = age_group, linetype = exposures, color = mod), position = position_dodge(0.8)) +
  geom_hline(yintercept = 1, linetype="dashed", color = 1, size = 0.2) +
  ylab("Rate ratio for an IQR increase\nwith 95% confidence interval") + xlab("Death cause") +
  labs(color = "Models") +
  labs(linetype = "Exposure sets") +
  labs(shape = "Age groups") +
  guides(color=guide_legend(nrow=2, override.aes=list(shape=c(NA,NA))), linetype=guide_legend(nrow=2, override.aes=list(shape=c(NA,NA))), shape=guide_legend(nrow=2, override.aes=list(linetype=c(0,0)))) +
  scale_x_discrete(labels=c("CVD" = "Cardiovascular\ndisease", "MI" = "Myocardial\ninfarction", "TOT" = "All-causes")) +
  theme_minimal()
plotpm

cairo_pdf(paste0(dir_plot, "RRiqr_pm_age.pdf"), height = 3.5)
plotpm
dev.off()

## exposure set as beta+pm25
plotDT <- results_PM_all
setDT(plotDT)
plotDT <- plotDT[exposures=="gross beta-activity + PM[2.5]"]
plotDT[, age_group:= factor(age_group, levels = c("18+", "18-65","65-85", "85+"))]
plotpm <- ggplot(plotDT, aes(x = cause, y = RR)) +
  geom_pointrange(size=0.5, aes(ymin = lowCI, ymax = highCI, shape = age_group, color = mod), position = position_dodge(0.8)) +
  geom_hline(yintercept = 1, linetype="dashed", color = 1, size = 0.2) +
  ylab("Rate ratio") + xlab("Death cause") +
  labs(color = "Models") +
  labs(shape = "Age groups") +
  guides(color=guide_legend(nrow=2, override.aes=list(shape=c(NA,NA))), shape=guide_legend(nrow=2, override.aes=list(linetype=c(0,0)))) +
  scale_x_discrete(labels=c("CVD" = "Cardiovascular\ndisease", "MI" = "Myocardial\ninfarction", "TOT" = "Non-accidental\n all causes")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_rect(colour = "black", fill=NA, size=1))
plotpm

cairo_pdf(paste0(dir_plot, "RRiqr_pm_age_betaPM.pdf"), height = 3.5)
plotpm
dev.off()

## 3.2 pm25 only all ages ----
plotDT <- results_PM_all
setDT(plotDT)
plotDT <- plotDT[age_group=="18+",]
plotpm_main <- ggplot(plotDT, aes(x = cause, y = RR)) +
  geom_pointrange(size=0.5, aes(ymin = lowCI, ymax = highCI,linetype = exposures, color = mod), position = position_dodge(0.8)) +
  geom_hline(yintercept = 1, linetype="dashed", color = 1, size = 0.2) +
  ylab("Rate ratio") + xlab("Death cause") +
  labs(color = "Models") +
  labs(linetype = "Exposure sets") +
  guides(color=guide_legend(nrow=2, override.aes=list(shape=c(NA,NA))), shape=guide_legend(nrow=2, override.aes=list(linetype=c(0,0)))) +
  scale_x_discrete(labels=c("CVD" = "Cardiovascular\ndisease", "MI" = "Myocardial\ninfarction", "stroke" = "Stroke","TOT" = "Non-accidental\nall causes")) +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_rect(colour = "black", fill=NA, size=1))
plotpm_main

cairo_pdf(paste0(dir_plot, "RRiqr_pm_main.pdf"), height = 3.5)
plotpm_main
dev.off()

## 4. configure main figures for each pollutant together ----
plotDT <- results_beta_all
setDT(plotDT)
plotDT <- plotDT[age_group=="18+",]
plotDT[exposures=="only gross beta-activity",]$exposures <- "only gross \u03B2-activity / PM[2.5]"
plotDT[exposures=="gross beta-activity + PM[2.5]",]$exposures <- "gross \u03B2-activity + PM[2.5]"
plotbeta_main <- ggplot(plotDT, aes(x = cause, y = RR)) +
  geom_pointrange(size=0.5, aes(ymin = lowCI, ymax = highCI, linetype = exposures, color = mod), position = position_dodge(0.8)) +
  geom_hline(yintercept = 1, linetype="dashed", color = 1, size = 0.2) +
  ylab("Rate ratio") + xlab("Death cause") +
  labs(color = "Models", linetype = "Exposure sets") +
  guides(color = guide_legend(nrow=2, override.aes=list(shape=c(NA,NA))), 
         linetype = guide_legend(nrow=2, override.aes=list(shape=c(NA,NA)))) +
  scale_x_discrete(labels=c("CVD" = "Cardiovascular\ndisease", "MI" = "Myocardial\ninfarction", "stroke" = "Stroke","TOT" = "Non-accidental\nall causes")) +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_rect(colour = "black", fill=NA, size=1))
plotbeta_main

plotDT <- results_PM_all
setDT(plotDT)
plotDT <- plotDT[age_group=="18+",]
plotDT[exposures=="only PM[2.5]",]$exposures <- "only gross \u03B2-activity / PM[2.5]"
plotDT[exposures=="gross beta-activity + PM[2.5]",]$exposures <- "gross \u03B2-activity + PM[2.5]"
plotpm_main <- ggplot(plotDT, aes(x = cause, y = RR)) +
  geom_pointrange(size=0.5, aes(ymin = lowCI, ymax = highCI,linetype = exposures, color = mod), position = position_dodge(0.8)) +
  geom_hline(yintercept = 1, linetype="dashed", color = 1, size = 0.2) +
  ylab("Rate ratio") + xlab("Death cause") +
  labs(color = "Models",linetype = "Exposure sets") +
  guides(color = guide_legend(nrow=2, override.aes=list(shape=c(NA,NA))), 
         linetype = guide_legend(nrow=2, override.aes=list(shape=c(NA,NA)))) +
  scale_x_discrete(labels=c("CVD" = "Cardiovascular\ndisease", "MI" = "Myocardial\ninfarction", "stroke" = "Stroke","TOT" = "Non-accidental\nall causes")) +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_rect(colour = "black", fill=NA, size=1))
plotpm_main

dir_plot <- "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/github_repo/results/main_PRPM/"
figure <- ggarrange(plotbeta_main, plotpm_main,
                    labels = c("A", "B"),
                    ncol = 1, nrow = 2,
                    legend = "top",
                    common.legend = TRUE)
figure

cairo_pdf(paste0(dir_plot, "RRiqr_main_configured.pdf"), height = 8)
figure
dev.off()
