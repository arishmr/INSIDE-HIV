####### GOAL: compare means of biomarkers across the two depression groups using Wilcoxon rank sum tests

library(plyr, include.only = c("rbind.fill"))

############
## DW-MRS
############

## Create output dataframe 
output <- data.frame(
  W = numeric(0),
  p = numeric(0),
  stringsAsFactors = F
)

## Select relevant data
data <- alldata %>%
  dplyr::select(Group,
                ADC.NAA:ADC.Cho) %>%
  gather(key = "biomarker",
         value = "ADC",
         ADC.NAA:ADC.Cho)
data$ADC <- log(data$ADC)

data <- data %>%
  group_by(biomarker) %>%
  group_split

output <- lapply(data, FUN = function(x) {
  test <- wilcox.test(ADC ~ Group, data = x)
  W <- as.numeric(test$statistic)
  p <- as.numeric(test$p.value)
  output[nrow(output)+1,] <- c(W, p)
  output
})
output <- rbind.fill(output)

output$biomarker <- rep(unlist(lapply(data, function (x) {c(unique(x$biomarker))})), each = 1)
output$p.fdr <- p.adjust(output$p, method = "fdr")
output <- output[c(
  "biomarker",
  "W", "p", "p.fdr"
)]

write.csv(output, "Results/Group Comparisons Non-Parametric - DW-MRS.csv", row.names = F)



############
## MRS
############

## Create output dataframe 
output <- data.frame(
  W = numeric(0),
  p = numeric(0),
  stringsAsFactors = F
)

## Select relevant data
data <- alldata %>%
  dplyr::select(Group,
                Asp:Cho, Asp.Cr:Cho.Cr) %>%
  gather(key = "biomarker",
         value = "conc",
         Asp:Cho.Cr)
data$conc <- log(data$conc)

data$referencing <- ifelse(str_detect(data$biomarker, '.Cr'), "Creatine Referenced", "Water Referenced")
data$referencing <- as.factor(data$referencing)
data <- data %>%
  dplyr::filter(referencing == "Creatine Referenced") %>%
  mutate(biomarker = recode(biomarker,
                            Asp.Cr = "Aspartate", Cho.Cr = "Choline", Glc.Cr = "Glucose", Glu.Gln.Cr = "Glutamate + Glutamine", Ins.Cr = "Myo-Inositol", Lac.Cr = "Lactate", NAA.Cr = "NAA"))

data <- data %>%
  group_by(referencing, biomarker) %>%
  group_split

output <- lapply(data, FUN = function(x) {
  test <- wilcox.test(conc ~ Group, data = x)
  W <- as.numeric(test$statistic)
  p <- test$p.value
  output[nrow(output)+1,] <- c(W, p)
  output
})
output <- rbind.fill(output)

output$referencing <- rep(unlist(lapply(data, function (x) {c(unique(x$referencing))})), each = 1)
output$biomarker <- rep(unlist(lapply(data, function (x) {c(unique(x$biomarker))})), each = 1)
output$p.fdr <- p.adjust(output$p, method = "fdr")
output <- output[c(
  "referencing", "biomarker",
  "W", "p", "p.fdr"
)]

write.csv(output, "Results/Group Comparisons Non-Parametric - MRS (Creatine).csv", row.names = F)

## Create output dataframe 
output <- data.frame(
  W = numeric(0),
  p = numeric(0),
  stringsAsFactors = F
)

## Select relevant data
data <- alldata %>%
  dplyr::select(Group,
                Asp:Cho, Asp.Cr:Cho.Cr) %>%
  gather(key = "biomarker",
         value = "conc",
         Asp:Cho.Cr)
data$conc <- log(data$conc)

data$referencing <- ifelse(str_detect(data$biomarker, '.Cr'), "Creatine Referenced", "Water Referenced")
data$referencing <- as.factor(data$referencing)
data <- data %>%
  dplyr::filter(referencing == "Water Referenced") %>%
  mutate(biomarker = recode(biomarker,
                            Asp = "Aspartate", Cr = "Creatine", Cho = "Choline", Glc = "Glucose", Glu.Gln = "Glutamate + Glutamine", Ins = "Myo-Inositol", Lac = "Lactate", NAA = "NAA"))

data <- data %>%
  group_by(referencing, biomarker) %>%
  group_split

output <- lapply(data, FUN = function(x) {
  test <- wilcox.test(conc ~ Group, data = x)
  W <- as.numeric(test$statistic)
  p <- test$p.value
  output[nrow(output)+1,] <- c(W, p)
  output
})
output <- rbind.fill(output)

output$referencing <- rep(unlist(lapply(data, function (x) {c(unique(x$referencing))})), each = 1)
output$biomarker <- rep(unlist(lapply(data, function (x) {c(unique(x$biomarker))})), each = 1)
output$p.fdr <- p.adjust(output$p, method = "fdr")
output <- output[c(
  "referencing", "biomarker",
  "W", "p", "p.fdr"
)]

write.csv(output, "Results/Group Comparisons Non-Parametric - MRS (Water).csv", row.names = F)



############
## DCE-MRI
############

## Create output dataframe 
output <- data.frame(
  W = numeric(0),
  p = numeric(0),
  stringsAsFactors = F
)

## Select relevant data
data <- alldata %>%
  dplyr::select(Group,
                K.trans_WB, K.trans_thal) %>%
  gather(key = "biomarker",
         value = "K.trans",
         K.trans_WB:K.trans_thal)
data$K.trans <- log(data$K.trans)

data <- data %>%
  group_by(biomarker) %>%
  group_split

output <- lapply(data, FUN = function(x) {
  test <- wilcox.test(K.trans ~ Group, data = x)
  W <- as.numeric(test$statistic)
  p <- as.numeric(test$p.value)
  output[nrow(output)+1,] <- c(W, p)
  output
})
output <- rbind.fill(output)

output$biomarker <- rep(unlist(lapply(data, function (x) {c(unique(x$biomarker))})), each = 1)
output$p.fdr <- p.adjust(output$p, method = "fdr")
output <- output[c(
  "biomarker",
  "W", "p", "p.fdr"
)]

write.csv(output, "Results/Group Comparisons Non-Parametric - DCE-MRI.csv", row.names = F)



rm(data, output)
