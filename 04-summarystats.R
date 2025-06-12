## Looking at depression scores for both questionnaires
hist(alldata$PHQ)
hist(log(alldata$PHQ+1))

median(alldata$age)
quantile(alldata$age, 0.25)
quantile(alldata$age, 0.75)

colnames(alldata)

## Looking at biomarker values
plotdata <- alldata %>%
  dplyr::select(PID, Group,
                PHQ,
                Asp:Cho, Asp.Cr:Cho.Cr,
                ADC.NAA:ADC.Cho,
                K.trans_WB, K.trans_thal) %>%
  gather(key = "biomarker",
         value = "conc",
         Asp:Cho, Asp.Cr:Cho.Cr,
         ADC.NAA:ADC.Cho,
         K.trans_WB, K.trans_thal,
         factor_key = T)

#plotdata$conc <- log(plotdata$conc)

ggplot(plotdata, aes(conc)) +
  stat_density() +
  facet_wrap(~ biomarker, scale = "free") +
  theme_minimal() +
  labs(x = "conc") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Figures/Biomarker Distributions.png", scale = 1.5, bg = "white")



################# SUMMARISE MEAN (SD) DEPRESSION SCORES BY DEPRESSION SEVERITY STATUS, AND OVERALL

summarydata <- alldata %>%
  dplyr::select(PID, Group, PHQ) %>%
  group_by(Group) %>%
  summarise(n = length(PHQ[!is.na(PHQ)]),
            median = median(PHQ, na.rm = T),
            Q1 = quantile(PHQ, 0.25, na.rm = TRUE),
            Q3 = quantile(PHQ, 0.75, na.rm = TRUE))

write.csv(summarydata, "Results/Summary Stats - Depression Score by Group.csv", row.names = F)

summarydata <- alldata %>%
  dplyr::select(PID, Group, PHQ) %>%
  summarise(n = length(PHQ[!is.na(PHQ)]),
            median = median(PHQ, na.rm = T),
            Q1 = quantile(PHQ, 0.25, na.rm = TRUE),
            Q3 = quantile(PHQ, 0.75, na.rm = TRUE))

write.csv(summarydata, "Results/Summary Stats - Depression Score Overall.csv", row.names = F)

#### WILCOXON RANK SUM TESTS FOR GROUP DIFFERENCES IN PHQ-9 SCORE
wilcox.test(alldata$PHQ ~ alldata$Group)


################# SUMMARISE MEAN (SD) BIOMARKER CONCENTRATIONS BY DEPRESSION SEVERITY STATUS, AND OVERALL


##################
#### MRS
##################

summarydata <- MRS %>%
  gather(key = "biomarker",
         value = "conc",
         Asp:Cho,
         factor_key = T)

summarydata$referencing <- rep(c("Water Referenced"))
summarydata$referencing <- as.factor(summarydata$referencing)

summarydata <- summarydata %>%
  group_by(referencing, biomarker, Group) %>%
  summarise(n = length(conc[!is.na(conc)]),
            median = median(conc, na.rm = T),
            Q1 = quantile(conc, 0.25, na.rm = TRUE),
            Q3 = quantile(conc, 0.75, na.rm = TRUE))

write.csv(summarydata, "Results/Summary Stats - MRS (Water) by Depression.csv", row.names = F)

summarydata <- MRS %>%
  gather(key = "biomarker",
         value = "conc",
         Asp:Cho,
         factor_key = T)

summarydata$referencing <- rep(c("Water Referenced"))
summarydata$referencing <- as.factor(summarydata$referencing)

summarydata <- summarydata %>%
  group_by(referencing, biomarker) %>%
  summarise(n = length(conc[!is.na(conc)]),
            median = median(conc, na.rm = T),
            Q1 = quantile(conc, 0.25, na.rm = TRUE),
            Q3 = quantile(conc, 0.75, na.rm = TRUE))

write.csv(summarydata, "Results/Summary Stats - MRS (Water) Overall.csv", row.names = F)


summarydata <- MRS %>%
  gather(key = "biomarker",
         value = "conc",
         Asp.Cr:Cho.Cr,
         factor_key = T)

summarydata$referencing <- rep(c("Creatine Referenced"))
summarydata$referencing <- as.factor(summarydata$referencing)

summarydata <- summarydata %>%
  group_by(referencing, biomarker, Group) %>%
  summarise(n = length(conc[!is.na(conc)]),
            median = median(conc, na.rm = T),
            Q1 = quantile(conc, 0.25, na.rm = TRUE),
            Q3 = quantile(conc, 0.75, na.rm = TRUE))

write.csv(summarydata, "Results/Summary Stats - MRS (Creatine) by Depression.csv", row.names = F)

summarydata <- MRS %>%
  gather(key = "biomarker",
         value = "conc",
         Asp.Cr:Cho.Cr,
         factor_key = T)

summarydata$referencing <- rep(c("Creatine Referenced"))
summarydata$referencing <- as.factor(summarydata$referencing)

summarydata <- summarydata %>%
  group_by(referencing, biomarker) %>%
  summarise(n = length(conc[!is.na(conc)]),
            median = median(conc, na.rm = T),
            Q1 = quantile(conc, 0.25, na.rm = TRUE),
            Q3 = quantile(conc, 0.75, na.rm = TRUE))

write.csv(summarydata, "Results/Summary Stats - MRS (Creatine) Overall.csv", row.names = F)


##################
#### DW-MRS
##################

summarydata <- DW.MRS %>%
  gather(key = "biomarker",
         value = "conc",
         ADC.NAA:ADC.Cho,
         factor_key = T) %>%
  group_by(biomarker, Group) %>%
  summarise(n = length(conc[!is.na(conc)]),
            median = median(conc, na.rm = T),
            Q1 = quantile(conc, 0.25, na.rm = TRUE),
            Q3 = quantile(conc, 0.75, na.rm = TRUE))

write.csv(summarydata, "Results/Summary Stats - DW-MRS by Depression.csv", row.names = F)

summarydata <- DW.MRS %>%
  gather(key = "biomarker",
         value = "conc",
         ADC.NAA:ADC.Cho,
         factor_key = T) %>%
  group_by(biomarker) %>%
  summarise(n = length(conc[!is.na(conc)]),
            median = median(conc, na.rm = T),
            Q1 = quantile(conc, 0.25, na.rm = TRUE),
            Q3 = quantile(conc, 0.75, na.rm = TRUE))

write.csv(summarydata, "Results/Summary Stats - DW-MRS Overall.csv", row.names = F)

##################
#### DCE-MRI
##################

summarydata <- DCE %>%
  group_by(Group) %>%
  summarise(n = length(K.trans_WB[!is.na(K.trans_WB)]),
            median = median(K.trans_WB, na.rm = T),
            Q1 = quantile(K.trans_WB, 0.25, na.rm = TRUE),
            Q3 = quantile(K.trans_WB, 0.75, na.rm = TRUE))

write.csv(summarydata, "Results/Summary Stats - DCE-MRI by Depression.csv", row.names = F)

summarydata <- DCE %>%
  summarise(n = length(K.trans_WB[!is.na(K.trans_WB)]),
            median = median(K.trans_WB, na.rm = T),
            Q1 = quantile(K.trans_WB, 0.25, na.rm = TRUE),
            Q3 = quantile(K.trans_WB, 0.75, na.rm = TRUE))

write.csv(summarydata, "Results/Summary Stats - DCE-MRI Overall.csv", row.names = F)


##################
#### BLOOD PROTEINS
##################

summarydata <- blooddata %>%
  gather(key = "biomarker",
         value = "conc",
         BDNF:YKL40,
         factor_key = T) %>%
  group_by(biomarker, Group) %>%
  summarise(n = length(conc[!is.na(conc)]),
            median = median(conc, na.rm = T),
            Q1 = quantile(conc, 0.25, na.rm = TRUE),
            Q3 = quantile(conc, 0.75, na.rm = TRUE))

write.csv(summarydata, "Results/Summary Stats - Blood Proteins by Depression.csv", row.names = F)

summarydata <- blooddata %>%
  gather(key = "biomarker",
         value = "conc",
         BDNF:YKL40,
         factor_key = T) %>%
  group_by(biomarker) %>%
  summarise(n = length(conc[!is.na(conc)]),
            median = median(conc, na.rm = T),
            Q1 = quantile(conc, 0.25, na.rm = TRUE),
            Q3 = quantile(conc, 0.75, na.rm = TRUE))

write.csv(summarydata, "Results/Summary Stats - Blood Proteins Overall.csv", row.names = F)


rm(summarydata, plotdata)
