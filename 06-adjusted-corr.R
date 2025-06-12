library(plyr, include.only = c("rbind.fill"))


################################################################
## Correlations between DW-MRS biomarkers and PHQ-9 adjusted for years of education / CD4 count
################################################################

##############################
#### DW-MRS adjusted for edu

set.seed(123)

data <- alldata %>%
  dplyr::select(PID, Group, PHQ,
                ADC.NAA:ADC.Cho,
                age, edu) %>%
  gather(key = "biomarker",
         value = "conc",
         ADC.NAA:ADC.Cho)
data$conc <- log(data$conc)
data <- na.omit(data)

splitdata <- data %>%
  group_by(biomarker) %>%
  group_split

cor.test(splitdata[[2]]$edu, splitdata[[2]]$conc)

output <- data.frame(
  r = numeric(0),
  p = numeric(0),
  n = numeric(0)
)

output <- lapply(splitdata, FUN = function (x) {
  fit <- pcor.test(x$PHQ, x$conc, x$edu, method = "spearman")
  fit
  r <- as.numeric(fit$estimate)
  p <- as.numeric(fit$p.value)
  n <- as.numeric(fit$n)
  output[nrow(output)+1,] <- c(r, p, n)
  output
})
output <- rbind.fill(output)

output$p.fdr <- p.adjust(output$p, method = "fdr")
output$biomarker <- unlist(lapply(splitdata, function (x) {c(unique(x$biomarker))}))
output <- output[c(5,3,1,2,4)]

write.csv(output, "Results/Biomarker x PHQ9 Correlations - DW-MRS (Edu Adjusted).csv", row.names = F)

rm(splitdata, data, output)

##############################
#### DW-MRS adjusted for CD4

set.seed(123)

data <- alldata %>%
  dplyr::select(PID, Group, PHQ,
                ADC.NAA:ADC.Cho,
                age, CD4) %>%
  gather(key = "biomarker",
         value = "conc",
         ADC.NAA:ADC.Cho)
data$conc <- log(data$conc)
data <- na.omit(data)

splitdata <- data %>%
  group_by(biomarker) %>%
  group_split

output <- data.frame(
  r = numeric(0),
  p = numeric(0),
  n = numeric(0)
)

output <- lapply(splitdata, FUN = function (x) {
  fit <- pcor.test(x$PHQ, x$conc, x$CD4, method = "spearman")
  fit
  r <- as.numeric(fit$estimate)
  p <- as.numeric(fit$p.value)
  n <- as.numeric(fit$n)
  output[nrow(output)+1,] <- c(r, p, n)
  output
})
output <- rbind.fill(output)

output$p.fdr <- p.adjust(output$p, method = "fdr")
output$biomarker <- unlist(lapply(splitdata, function (x) {c(unique(x$biomarker))}))
output <- output[c(5,3,1,2,4)]

write.csv(output, "Results/Biomarker x PHQ9 Correlations - DW-MRS (CD4 Adjusted).csv", row.names = F)

rm(splitdata, data, output)



################################################################
## Correlations between DCE-MRI biomarkers and PHQ-9 adjusted for years of education / CD4 count
################################################################


##############################
#### DCE-MRI adjusted for edu

set.seed(123)

data <- alldata %>%
  dplyr::select(PID, Group, PHQ,
                K.trans_WB, K.trans_thal,
                age, edu) %>%
  gather(key = "biomarker",
         value = "conc",
         K.trans_WB, K.trans_thal)
data$conc <- log(data$conc)
data <- na.omit(data)

splitdata <- data %>%
  group_by(biomarker) %>%
  group_split

output <- data.frame(
  r = numeric(0),
  p = numeric(0),
  n = numeric(0)
)

output <- lapply(splitdata, FUN = function (x) {
  fit <- pcor.test(x$PHQ, x$conc, x$edu, method = "spearman")
  fit
  r <- as.numeric(fit$estimate)
  p <- as.numeric(fit$p.value)
  n <- as.numeric(fit$n)
  output[nrow(output)+1,] <- c(r, p, n)
  output
})
output <- rbind.fill(output)

output$p.fdr <- p.adjust(output$p, method = "fdr")
output$biomarker <- unlist(lapply(splitdata, function (x) {c(unique(x$biomarker))}))
output <- output[c(5,3,1,2,4)]

write.csv(output, "Results/Biomarker x PHQ9 Correlations - DCE-MRI (Edu Adjusted).csv", row.names = F)

rm(splitdata, data, output)

##############################
#### DCE-MRI adjusted for CD4

set.seed(123)

data <- alldata %>%
  dplyr::select(PID, Group, PHQ,
                K.trans_WB, K.trans_thal,
                age, CD4) %>%
  gather(key = "biomarker",
         value = "conc",
         K.trans_WB, K.trans_thal)
data$conc <- log(data$conc)
data <- na.omit(data)

splitdata <- data %>%
  group_by(biomarker) %>%
  group_split

output <- data.frame(
  r = numeric(0),
  p = numeric(0),
  n = numeric(0)
)

output <- lapply(splitdata, FUN = function (x) {
  fit <- pcor.test(x$PHQ, x$conc, x$CD4, method = "spearman")
  fit
  r <- as.numeric(fit$estimate)
  p <- as.numeric(fit$p.value)
  n <- as.numeric(fit$n)
  output[nrow(output)+1,] <- c(r, p, n)
  output
})
output <- rbind.fill(output)

output$p.fdr <- p.adjust(output$p, method = "fdr")
output$biomarker <- unlist(lapply(splitdata, function (x) {c(unique(x$biomarker))}))
output <- output[c(5,3,1,2,4)]

write.csv(output, "Results/Biomarker x PHQ9 Correlations - DCE-MRI (CD4 Adjusted).csv", row.names = F)

rm(splitdata, data, output)




################################################################
## Correlations between MRS biomarkers and PHQ-9 adjusted for years of education / CD4 count + fWM
################################################################



##############################
#### MRS adjusted for edu

set.seed(123)

data <- alldata %>%
  dplyr::select(PID, Group, PHQ,
                Asp:Cho.Cr,
                fWM, age, edu, CD4) %>%
  gather(key = "biomarker",
         value = "conc",
         Asp:Cho.Cr)
data$conc <- log(data$conc)

data$referencing <- ifelse(str_detect(data$biomarker, '.Cr'), "Creatine Referenced", "Water Referenced")
data$referencing <- as.factor(data$referencing)
data <- data  %>%
  dplyr::filter(referencing == "Water Referenced") %>%
  mutate(biomarker = recode(biomarker,
                            Asp = "Aspartate", Cr = "Creatine", Cho = "Choline", Glc = "Glucose", Glu.Gln = "Glutamate + Glutamine", Ins = "Myo-Inositol", Lac = "Lactate", NAA = "NAA"))

data <- na.omit(data)

splitdata <- data %>%
  group_by(referencing, biomarker) %>%
  group_split

output <- data.frame(
  r = numeric(0),
  p = numeric(0),
  n = numeric(0)
)

output <- lapply(splitdata, FUN = function (x) {
  fit <- pcor.test(x$PHQ, x$conc, list(x$fWM, x$edu), method = "spearman")
  fit
  r <- as.numeric(fit$estimate)
  p <- as.numeric(fit$p.value)
  n <- as.numeric(fit$n)
  output[nrow(output)+1,] <- c(r, p, n)
  output
})
output <- rbind.fill(output)

output$p.fdr <- p.adjust(output$p, method = "fdr")
output$referencing <- unlist(lapply(splitdata, function (x) {c(unique(x$referencing))}))
output$biomarker <- unlist(lapply(splitdata, function (x) {c(unique(x$biomarker))}))
output <- output[c(5:6,3,1,2,4)]

write.csv(output, "Results/Biomarker x PHQ9 Correlations - MRS (Water, Edu & fWM Adjusted).csv", row.names = F)



set.seed(123)

data <- alldata %>%
  dplyr::select(PID, Group, PHQ,
                Asp:Cho.Cr,
                fWM, age, edu, CD4) %>%
  gather(key = "biomarker",
         value = "conc",
         Asp:Cho.Cr)
data$conc <- log(data$conc)

data$referencing <- ifelse(str_detect(data$biomarker, '.Cr'), "Creatine Referenced", "Water Referenced")
data$referencing <- as.factor(data$referencing)
data <- data  %>%
  dplyr::filter(referencing == "Creatine Referenced") %>%
  mutate(biomarker = recode(biomarker,
                            Asp.Cr = "Aspartate", Cho.Cr = "Choline", Glc.Cr = "Glucose", Glu.Gln.Cr = "Glutamate + Glutamine", Ins.Cr = "Myo-Inositol", Lac.Cr = "Lactate", NAA.Cr = "NAA"))

data <- na.omit(data)

splitdata <- data %>%
  group_by(referencing, biomarker) %>%
  group_split

output <- data.frame(
  r = numeric(0),
  p = numeric(0),
  n = numeric(0)
)

output <- lapply(splitdata, FUN = function (x) {
  fit <- pcor.test(x$PHQ, x$conc, list(x$fWM, x$edu), method = "spearman")
  fit
  r <- as.numeric(fit$estimate)
  p <- as.numeric(fit$p.value)
  n <- as.numeric(fit$n)
  output[nrow(output)+1,] <- c(r, p, n)
  output
})
output <- rbind.fill(output)

output$p.fdr <- p.adjust(output$p, method = "fdr")
output$referencing <- unlist(lapply(splitdata, function (x) {c(unique(x$referencing))}))
output$biomarker <- unlist(lapply(splitdata, function (x) {c(unique(x$biomarker))}))
output <- output[c(5:6,3,1,2,4)]

write.csv(output, "Results/Biomarker x PHQ9 Correlations - MRS (Creatine, Edu & fWM Adjusted).csv", row.names = F)

rm(splitdata, data, output)


##############################
#### MRS adjusted for CD4

set.seed(123)

data <- alldata %>%
  dplyr::select(PID, Group, PHQ,
                Asp:Cho.Cr,
                fWM, age, edu, CD4) %>%
  gather(key = "biomarker",
         value = "conc",
         Asp:Cho.Cr)
data$conc <- log(data$conc)

data$referencing <- ifelse(str_detect(data$biomarker, '.Cr'), "Creatine Referenced", "Water Referenced")
data$referencing <- as.factor(data$referencing)
data <- data  %>%
  dplyr::filter(referencing == "Water Referenced") %>%
  mutate(biomarker = recode(biomarker,
                            Asp = "Aspartate", Cr = "Creatine", Cho = "Choline", Glc = "Glucose", Glu.Gln = "Glutamate + Glutamine", Ins = "Myo-Inositol", Lac = "Lactate", NAA = "NAA"))

data <- na.omit(data)

splitdata <- data %>%
  group_by(referencing, biomarker) %>%
  group_split

output <- data.frame(
  r = numeric(0),
  p = numeric(0),
  n = numeric(0)
)

output <- lapply(splitdata, FUN = function (x) {
  fit <- pcor.test(x$PHQ, x$conc, list(x$fWM, x$CD4), method = "spearman")
  fit
  r <- as.numeric(fit$estimate)
  p <- as.numeric(fit$p.value)
  n <- as.numeric(fit$n)
  output[nrow(output)+1,] <- c(r, p, n)
  output
})
output <- rbind.fill(output)

output$p.fdr <- p.adjust(output$p, method = "fdr")
output$referencing <- unlist(lapply(splitdata, function (x) {c(unique(x$referencing))}))
output$biomarker <- unlist(lapply(splitdata, function (x) {c(unique(x$biomarker))}))
output <- output[c(5:6,3,1,2,4)]

write.csv(output, "Results/Biomarker x PHQ9 Correlations - MRS (Water, CD4 & fWM Adjusted).csv", row.names = F)



set.seed(123)

data <- alldata %>%
  dplyr::select(PID, Group, PHQ,
                Asp:Cho.Cr,
                fWM, age, edu, CD4) %>%
  gather(key = "biomarker",
         value = "conc",
         Asp:Cho.Cr)
data$conc <- log(data$conc)

data$referencing <- ifelse(str_detect(data$biomarker, '.Cr'), "Creatine Referenced", "Water Referenced")
data$referencing <- as.factor(data$referencing)
data <- data  %>%
  dplyr::filter(referencing == "Creatine Referenced") %>%
  mutate(biomarker = recode(biomarker,
                            Asp.Cr = "Aspartate", Cho.Cr = "Choline", Glc.Cr = "Glucose", Glu.Gln.Cr = "Glutamate + Glutamine", Ins.Cr = "Myo-Inositol", Lac.Cr = "Lactate", NAA.Cr = "NAA"))

data <- na.omit(data)

splitdata <- data %>%
  group_by(referencing, biomarker) %>%
  group_split

output <- data.frame(
  r = numeric(0),
  p = numeric(0),
  n = numeric(0)
)

output <- lapply(splitdata, FUN = function (x) {
  fit <- pcor.test(x$PHQ, x$conc, list(x$fWM, x$CD4), method = "spearman")
  fit
  r <- as.numeric(fit$estimate)
  p <- as.numeric(fit$p.value)
  n <- as.numeric(fit$n)
  output[nrow(output)+1,] <- c(r, p, n)
  output
})
output <- rbind.fill(output)

output$p.fdr <- p.adjust(output$p, method = "fdr")
output$referencing <- unlist(lapply(splitdata, function (x) {c(unique(x$referencing))}))
output$biomarker <- unlist(lapply(splitdata, function (x) {c(unique(x$biomarker))}))
output <- output[c(5:6,3,1,2,4)]

write.csv(output, "Results/Biomarker x PHQ9 Correlations - MRS (Creatine, CD4 & fWM Adjusted).csv", row.names = F)

rm(splitdata, data, output)

