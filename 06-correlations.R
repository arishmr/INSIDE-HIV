library(plyr, include.only = c("rbind.fill"))

################################################################
## Correlations between DW-MRS / DCE-MRI biomarkers and PHQ-9
################################################################

set.seed(123)

data <- alldata %>%
  dplyr::select(PID, Group, PHQ,
         ADC.NAA:ADC.Cho) %>%
  gather(key = "biomarker",
         value = "conc",
         ADC.NAA:ADC.Cho)
data$conc <- log(data$conc)

splitdata <- data %>%
  group_by(biomarker) %>%
  group_split

output <- data.frame(
  r = numeric(0),
  p.value = numeric(0)
)

output <- lapply(splitdata, FUN = function (x) {
  fit <- cor.test(x$PHQ, x$conc, method = "spearman")
  fit
  r <- as.numeric(fit$estimate)
  p <- as.numeric(fit$p.value)
  output[nrow(output)+1,] <- c(r, p)
  output
})
output <- rbind.fill(output)

output$p.fdr <- p.adjust(output$p, method = "fdr")
output$biomarker <- unlist(lapply(splitdata, function (x) {c(unique(x$biomarker))}))
output <- output[c(4,1:3)]

write.csv(output, "Results/Biomarker x PHQ9 Correlations - DW-MRS.csv", row.names = F)



set.seed(123)

data <- alldata %>%
  dplyr::select(PID, Group, PHQ,
                K.trans_WB, K.trans_thal) %>%
  gather(key = "biomarker",
         value = "conc",
         K.trans_WB:K.trans_thal)
data$conc <- log(data$conc)

splitdata <- data %>%
  group_by(biomarker) %>%
  group_split

output <- data.frame(
  r = numeric(0),
  p.value = numeric(0)
)

output <- lapply(splitdata, FUN = function (x) {
  fit <- cor.test(x$PHQ, x$conc, method = "spearman")
  fit
  r <- as.numeric(fit$estimate)
  p <- as.numeric(fit$p.value)
  output[nrow(output)+1,] <- c(r, p)
  output
})
output <- rbind.fill(output)

output$p.fdr <- p.adjust(output$p, method = "fdr")
output$biomarker <- unlist(lapply(splitdata, function (x) {c(unique(x$biomarker))}))
output <- output[c(4,1:3)]

write.csv(output, "Results/Biomarker x PHQ9 Correlations - DCE-MRI.csv", row.names = F)

rm(splitdata, data, output)



################################################################
## Correlations between MRS biomarkers and PHQ-9
################################################################

set.seed(123)

data <- alldata %>%
  dplyr::select(PID, Group, PHQ,
                Asp:Cho.Cr) %>%
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

splitdata <- data %>%
  group_by(referencing, biomarker) %>%
  group_split

output <- data.frame(
  r = numeric(0),
  p.value = numeric(0)
)

output <- lapply(splitdata, FUN = function (x) {
  fit <- cor.test(x$PHQ, x$conc, method = "spearman")
  fit
  r <- as.numeric(fit$estimate)
  p <- as.numeric(fit$p.value)
  output[nrow(output)+1,] <- c(r, p)
  output
})
output <- rbind.fill(output)

output$p.fdr <- p.adjust(output$p, method = "fdr")
output$referencing <- unlist(lapply(splitdata, function (x) {c(unique(x$referencing))}))
output$biomarker <- unlist(lapply(splitdata, function (x) {c(unique(x$biomarker))}))
output <- output[c(4:5,1:3)]

write.csv(output, "Results/Biomarker x PHQ9 Correlations - MRS (Water).csv", row.names = F)

set.seed(123)

data <- alldata %>%
  dplyr::select(PID, Group, PHQ,
                Asp:Cho.Cr) %>%
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

splitdata <- data %>%
  group_by(referencing, biomarker) %>%
  group_split

output <- data.frame(
  r = numeric(0),
  p.value = numeric(0)
)

output <- lapply(splitdata, FUN = function (x) {
  fit <- cor.test(x$PHQ, x$conc, method = "spearman")
  fit
  r <- as.numeric(fit$estimate)
  p <- as.numeric(fit$p.value)
  output[nrow(output)+1,] <- c(r, p)
  output
})
output <- rbind.fill(output)

output$p.fdr <- p.adjust(output$p, method = "fdr")
output$referencing <- unlist(lapply(splitdata, function (x) {c(unique(x$referencing))}))
output$biomarker <- unlist(lapply(splitdata, function (x) {c(unique(x$biomarker))}))
output <- output[c(4:5,1:3)]

write.csv(output, "Results/Biomarker x PHQ9 Correlations - MRS (Creatine).csv", row.names = F)


rm(splitdata, data, output)


################################################################
## Correlations between MRS biomarkers and PHQ-9 adjusted for fWM
################################################################

set.seed(123)

data <- alldata %>%
  dplyr::select(PID, Group, PHQ,
                Asp:Cho.Cr,
                fWM) %>%
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
  fit <- pcor.test(x$PHQ, x$conc, x$fWM, method = "spearman")
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

write.csv(output, "Results/Biomarker x PHQ9 Correlations - MRS (Water, fWM Adjusted).csv", row.names = F)



set.seed(123)

data <- alldata %>%
  dplyr::select(PID, Group, PHQ,
                Asp:Cho.Cr,
                fWM) %>%
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
  fit <- pcor.test(x$PHQ, x$conc, x$fWM, method = "spearman")
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

write.csv(output, "Results/Biomarker x PHQ9 Correlations - MRS (Creatine, fWM Adjusted).csv", row.names = F)

rm(splitdata, data, output)

