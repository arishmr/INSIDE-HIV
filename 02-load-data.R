rm(list = ls())

## Load datasets
alldata <- read.csv("Data/Full_Dataset.csv")
colnames(alldata)

# convert variables to correct type

alldata <- alldata %>% mutate_at(c(
  "PID", "PHQ", "age", "edu",
  "CD4",   ## clinical
  colnames(alldata[,18:26]),   ## PHQ9 scores
  colnames(alldata[,27:44]),   ## MRS data
  colnames(alldata[,46:50]),   ## DW-MRS and DCE-MRI data
  colnames(alldata[,51:70])    ## blood proteins
), as.numeric)

alldata <- alldata %>% mutate_at(c(
  "Group", "sex", "gender", "ethnicity", "smoking", "alcohol", "recdrugs",
  "curr.antidep", "ever.antidep", "curr.therapy", "ever.therapy", "ViralLoad"
), as.factor)

alldata$Group <- relevel(alldata$Group, ref = "Low Depression Severity")
alldata$scandate <- as.Date(alldata$scandate)


## Extract different datasets
DCE <- alldata %>%
  dplyr::select(PID, K.trans_WB, K.trans_thal, Group)

DW.MRS <- alldata %>%
  dplyr::select(PID, scandate, ADC.NAA, ADC.Cr, ADC.Cho, Group)

MRS <- alldata %>%
  dplyr::select(PID, Asp, Cr, Glc, Glu.Gln, Ins, Lac, NAA, Cho,
                Asp.Cr, Glc.Cr, Glu.Gln.Cr, Ins.Cr, Lac.Cr, NAA.Cr, Cho.Cr,
                fGM, fWM, fCSF,
                Group)

blooddata <- alldata %>% dplyr::select(PID, BDNF:YKL40, Group, PHQ)
  