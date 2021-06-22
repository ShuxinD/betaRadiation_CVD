###############################################################################
# Project: beta radiation and CVD death in MA
# Code: get ICD codes to identify deaths
# Input: CVD related ICD codes from the grant
# Output: "icd_mort.RData"
# Author: Shuxin Dong                                                         
###############################################################################

############################# 0. Setup ########################################
rm(list = ls())
gc()

library(data.table)
# library(icd)
library(haven)
setDTthreads(threads = 0)

dir_data <- "/media/gate/Shuxin/betaRadiation_CVD/"
dir_output <- "/media/gate/Shuxin/betaRadiation_CVD/"

########################### 1.  ICD codes #####################################
# total cardiovascular disease (CVD, ICD-9: 390-429, ICD‑10: I00 through I519)
# congestive heart failure (CHF, ICD-9: 428; ICD-10: I50-I51)
# myocardial infarction (MI, ICD-9: 410; ICD-10: I21-I22)
# ischemic stroke (ICD-9: 434,436)
# hemorrhagic stroke (ICD-9: 431)
# stroke (ICD-10: I60 through I69)
# all-cause mortality (TOT; ICD‑10: A00 through R99)

CVD <- expand_range("I00", "I519")
CHF <- expand_range("I50", "I51")
MI <- expand_range("I21", "I22")
stroke <- expand_range("I60", "I69")
allcause <- expand_range("A00", "R99")
save.image("icd_mort.RData")
#load("~/Documents/GitHub/betaRadiation_CVD/icd_mort.RData")
