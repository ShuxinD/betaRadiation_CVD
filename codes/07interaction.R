#' Project: beta radiation and CVD death in MA
#' Code: interaction analysis
#' Input: 
#' Output: ...
#' Author: Shuxin Dong


# 0. setup------
rm(list=ls())
gc()

library(mgcv)
library(data.table)
library(ggplot2)

dir_in <- "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/"
dir_out_detail <- "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/betaRadiation_CVD/results/interactionmodresults_details/"
dir_out_tb <- "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/betaRadiation_CVD/results/"

dt <- readRDS(paste0(dir_in, "finalDT.rds"))
dt <- na.omit(dt)
dt_pcount <- dt[,pcount]
names(dt)

# 1. model ----
quartiles <- quantile(dt[,Beta])
deathc <- c("TOT", "CVD", "MI", "stroke")

# 1.1 DID ----
for (deathc_i in deathc) {
  mod <- gam(get(deathc_i) ~ pm25 + Beta + I(pm25*Beta) + 
               summer_tmean + winter_tmean + as.factor(year) + as.factor(ZCTA5CE10) + offset(log(dt_pcount)),
             data = dt,
             family = quasipoisson(link = "log"))
  tb <- summary(mod)$p.table
  vcovm <- vcov(mod)
  coef_pm <- matrix(NA, nrow=3, ncol = 2)
  #' 1st quartile
  # quartiles[2]
  coef_pm[1,1] <- tb[2,1] + tb[4,1]*quartiles[2] # point estimate
  coef_pm[1,2] <- sqrt(vcovm[2,2] + (quartiles[2]^2)*vcovm[4,4]+ 2*quartiles[2]*vcovm[2,4])
  #' median
  # quartiles[3]
  coef_pm[2,1] <- tb[2,1] + tb[4,1]*quartiles[3] # point estimate
  coef_pm[2,2] <- sqrt(vcovm[2,2] + (quartiles[3]^2)*vcovm[4,4]+ 2*quartiles[3]*vcovm[2,4])
  #' 3rd quartile
  # quartiles[4]
  coef_pm[3,1] <- tb[2,1] + tb[4,1]*quartiles[4] # point estimate
  coef_pm[3,2] <- sqrt(vcovm[2,2] + (quartiles[4]^2)*vcovm[4,4]+ 2*quartiles[4]*vcovm[2,4])
  coef_pm
  rownames(coef_pm) <- c("pm25_beta25","pm25_beta50","pm25_beta75")
  colnames(coef_pm) <- c("estimate", "se")
  write.csv(coef_pm, file = paste0(dir_out_detail, "toxicity_DID_X_", deathc_i, ".csv"))
  cat(paste0("finish DID ", deathc_i, "\n"))
}

toxDID_tb <- NULL
rownames_list <- NULL
for (deathc_i in deathc) {
  RR_i <- NULL
  tb <- read.csv(paste0(dir_out_detail, "toxicity_DID_X_", deathc_i, ".csv"), row.names = 1)
  RR_i <- cbind(tb[,1], tb[,1]-1.96*tb[,2], tb[,1]+1.96*tb[,2])
  RR_i <- exp(RR_i*IQR(dt[,pm25]))
  toxDID_tb <- rbind(toxDID_tb, RR_i)
  rownames_list <- c(rownames_list, rep(deathc_i,3))
}
rownames(toxDID_tb) <- rownames_list
colnames(toxDID_tb) <- c("RR", "RR_lci", "RR_uci")
print(toxDID_tb)
write.csv(toxDID_tb, paste0(dir_out_tb, "toxDID_X_tb.csv"))

# 1.2 GLMM ----
for (deathc_i in deathc) {
  mod <- gamm(get(deathc_i) ~pm25 + Beta + I(pm25*Beta) + 
                summer_tmean + winter_tmean + medhouseholdincome + medianhousevalue + hispanic + pct_blk + poverty + education + popdensity + mean_bmi + smoke_rate + as.factor(year) + offset(log(dt_pcount)),
              data = dt,
              random = list(ZCTA5CE10=~1),
              family = quasipoisson(link = "log"))
  tb <- summary(mod$gam)$p.table
  vcovm <- vcov(mod$gam)
  coef_pm <- matrix(NA, nrow=3, ncol = 2)
  #' 1st quartile
  # quartiles[2]
  coef_pm[1,1] <- tb[2,1] + tb[4,1]*quartiles[2] # point estimate
  coef_pm[1,2] <- sqrt(vcovm[2,2] + (quartiles[2]^2)*vcovm[4,4]+ 2*quartiles[2]*vcovm[2,4])
  #' median
  # quartiles[3]
  coef_pm[2,1] <- tb[2,1] + tb[4,1]*quartiles[3] # point estimate
  coef_pm[2,2] <- sqrt(vcovm[2,2] + (quartiles[3]^2)*vcovm[4,4]+ 2*quartiles[3]*vcovm[2,4])
  #' 3rd quartile
  # quartiles[4]
  coef_pm[3,1] <- tb[2,1] + tb[4,1]*quartiles[4] # point estimate
  coef_pm[3,2] <- sqrt(vcovm[2,2] + (quartiles[4]^2)*vcovm[4,4]+ 2*quartiles[4]*vcovm[2,4])
  coef_pm
  rownames(coef_pm) <- c("pm25_beta25","pm25_beta50","pm25_beta75")
  colnames(coef_pm) <- c("estimate", "se")
  write.csv(coef_pm, file = paste0(dir_out_detail, "toxicity_GLMM_X_", deathc_i, ".csv"))
  cat(paste0("finish GLMM ", deathc_i, "\n"))
}

toxGLMM_tb <- NULL
rownames_list <- NULL
for (deathc_i in deathc) {
  RR_i <- NULL
  tb <- read.csv(paste0(dir_out_detail, "toxicity_GLMM_X_", deathc_i, ".csv"), row.names = 1)
  RR_i <- cbind(tb[,1], tb[,1]-1.96*tb[,2], tb[,1]+1.96*tb[,2])
  RR_i <- exp(RR_i*IQR(dt[,pm25]))
  toxGLMM_tb <- rbind(toxGLMM_tb, RR_i)
  rownames_list <- c(rownames_list, rep(deathc_i,3))
}
rownames(toxGLMM_tb) <- rownames_list
colnames(toxGLMM_tb) <- c("RR", "RR_lci", "RR_uci")
print(toxGLMM_tb)
write.csv(toxGLMM_tb, paste0(dir_out_tb, "toxGLMM_X_tb.csv"))

# 2. plots ----
dir_plot <- "/media/qnap3/Shuxin/ParticalRadiation_MAdeath/betaRadiation_CVD/results/"
## prepare dataset
toxDID_tb <- read.csv(paste0(dir_out_tb, "toxDID_X_tb.csv"))
toxGLMM_tb <- read.csv(paste0(dir_out_tb, "toxGLMM_X_tb.csv"))
plot_toxDID <- as.data.table(toxDID_tb)
plot_toxDID[,`:=`(cause=X,
                  beta_group=rep(c("beta25","beta50","beta75"), 4))]
plot_toxGLMM <- as.data.table(toxGLMM_tb)
plot_toxGLMM[,`:=`(cause=X,
                   beta_group=rep(c("beta25","beta50","beta75"), 4))]
plot_toxDID[,mod:="DID"]
plot_toxGLMM[,mod:="GLMM"]

## all models together
plotDT <- rbind(plot_toxDID, plot_toxGLMM)
plotDT[,beta_group:=as.factor(beta_group)]

plotpm <- ggplot(plotDT, aes(x = cause, y = RR)) +
  geom_pointrange(size=0.5, aes(ymin = RR_lci, ymax = RR_uci, color = beta_group, shape = mod), position = position_dodge(0.8)) + 
  geom_hline(yintercept = 1, linetype="dashed", color = 1, size = 0.2) + 
  ylab("Rate ratio") + xlab("Death cause") +
  labs(color = "Gross beta activity quartile", shape = "Model") + 
  # guides(color=guide_legend(nrow=2, override.aes=list(shape=c(NA,NA))), shape=guide_legend(nrow=2, override.aes=list(linetype=c(0,0)))) +
  scale_x_discrete(labels=c("CVD" = "Cardiovascular\n disease", "MI" = "Myocardial\n infarction","stroke"="Stroke", "TOT" = "Non-accidental\nall-causes")) + 
  scale_color_discrete(labels=c("1st quartile", "2nd quartile", "3rd quartile")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), panel.border = element_rect(colour = "black", fill=NA, size=1))
plotpm
  
pdf(paste0(dir_plot, "PMxBeta_PM.pdf"), height = 3.5)
plotpm
dev.off()
