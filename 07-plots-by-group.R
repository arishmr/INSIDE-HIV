################### GOAL: Plot neuroimaging variables by group

#### MRS DATA

plotdata <- alldata %>%
  dplyr::select(PID, Group, PHQ,
                Asp:Cho, Asp.Cr:Cho.Cr)
levels(plotdata$Group) <- c("Low", "High")

plotdata <- plotdata %>%
  gather(key = "biomarker",
         value = "conc",
         Asp:Cho, Asp.Cr:Cho.Cr,
         factor_key = T)

plotdata$referencing <- ifelse(str_detect(plotdata$biomarker, '.Cr'), "Creatine Referenced", "Water Referenced")
plotdata$referencing <- as.factor(plotdata$referencing)

plotdata <- plotdata %>%
  mutate(biomarker = recode(biomarker,
                            Asp = "Aspartate", Cr = "Creatine", Cho = "Choline (PCho + GPC)", Glc = "Glucose", Glu.Gln = "Glutamate + Glutamine", Ins = "Myo-Inositol", Lac = "Lactate", NAA = "NAA + NAAG", Asp.Cr = "Aspartate", Cho.Cr = "Choline (PCho + GPC)", Glc.Cr = "Glucose", Glu.Gln.Cr = "Glutamate + Glutamine", Ins.Cr = "Myo-Inositol", Lac.Cr = "Lactate", NAA.Cr = "NAA + NAAG"))

plotdataCr <- plotdata %>%
  filter(referencing == "Creatine Referenced")

gplot <- ggplot(plotdataCr, aes(Group, conc, fill = Group, color = Group)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center') +
  #stat_summary(geom = "errorbar", fun.data = mean_se, width = 0.6, alpha = 0.1) +
  #stat_summary(geom = "errorbar", fun.min = function(z) { quantile(z,0.25) },
  #             fun = median, fun.max = function(z) { quantile(z,0.75) }, width = 0.6, alpha = 0.1) +
  stat_summary(geom = "errorbar", fun.min = median, fun = median, fun.max = median, width = 0.6, alpha = 0.3) +
  facet_wrap(factor(plotdataCr$biomarker), scale="free_y") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Depressive Symptom Severity", y = "Concentration (Ratio to Creatine)") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(title = element_text(face = "bold"), axis.title = element_text(face="bold"))
gplot
ggsave("Figures/Group Comparisons - MRS (Creatine).png", gplot, width = 6, height = 6)


plotdataW <- plotdata %>%
  filter(referencing == "Water Referenced")

gplot <- ggplot(plotdataW, aes(Group, conc, fill = Group, color = Group)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center') +
  #stat_summary(geom = "errorbar", fun.data = mean_se, width = 0.6, alpha = 0.1) +
  #stat_summary(geom = "errorbar", fun.min = function(z) { quantile(z,0.25) },
  #             fun = median, fun.max = function(z) { quantile(z,0.75) }, width = 0.6, alpha = 0.1) +
  stat_summary(geom = "errorbar", fun.min = median, fun = median, fun.max = median, width = 0.6, alpha = 0.3) +
  facet_wrap(factor(plotdataW$biomarker), scale="free_y") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Depressive Symptom Severity", y = "Concentration (Ratio to Water)") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(title = element_text(face = "bold"), axis.title = element_text(face="bold"))
gplot
ggsave("Figures/Group Comparisons - MRS (Water).png", gplot, width = 6, height = 6)




#### DW-MRS DATA

plotdata <- alldata %>%
  dplyr::select(PID, Group, PHQ,
                ADC.NAA:ADC.Cho)
levels(plotdata$Group) <- c("Low", "High")

plotdata <- plotdata %>%
  gather(key = "biomarker",
         value = "conc",
         ADC.NAA:ADC.Cho,
         factor_key = T)

plotdata <- plotdata %>%
  mutate(biomarker = recode(biomarker,
                            ADC.NAA = "NAA + NAAG",
                            ADC.Cr = "Creatine (Cr + PCr)", ADC.Cho = "Choline (PCho + GPC)"))

gplot <- ggplot(plotdata, aes(Group, conc, fill = Group, color = Group)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center') +
  #stat_summary(geom = "errorbar", fun.data = mean_se, width = 0.6, alpha = 0.1) +
  #stat_summary(geom = "errorbar", fun.min = function(z) { quantile(z,0.25) },
  #             fun = median, fun.max = function(z) { quantile(z,0.75) }, width = 0.6, alpha = 0.1) +
  stat_summary(geom = "errorbar", fun.min = median, fun = median, fun.max = median, width = 0.6, alpha = 0.3) +
  facet_wrap(factor(plotdata$biomarker), scale="free_y") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Depressive Symptom Severity", y = "Apparent Diffusion Coefficient (ADC) (\u00b5m\u00b3/ms)") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(title = element_text(face = "bold"), axis.title = element_text(face="bold"))
gplot
ggsave("Figures/Group Comparisons - DW-MRS.png", gplot, width = 10, height = 4)




#### DCE-MRI DATA

plotdata <- alldata %>%
  dplyr::select(PID, Group, PHQ,
                K.trans_WB)
levels(plotdata$Group) <- c("Low", "High")

gplot <- ggplot(plotdata, aes(Group, K.trans_WB, fill = Group, color = Group)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center') +
  #stat_summary(geom = "errorbar", fun.data = mean_se, width = 0.6, alpha = 0.1) +
  #stat_summary(geom = "errorbar", fun.min = function(z) { quantile(z,0.25) },
  #             fun = median, fun.max = function(z) { quantile(z,0.75) }, width = 0.6, alpha = 0.1) +
  stat_summary(geom = "errorbar", fun.min = median, fun = median, fun.max = median, width = 0.6, alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Depressive Symptom Severity", y = "Volume Transfer Constant (K.trans) (/min)") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(title = element_text(face = "bold"), axis.title = element_text(face="bold"))
gplot
ggsave("Figures/Group Comparisons - DCE-MRI (Whole-Brain).png", gplot, width = 4, height = 4)



plotdata <- alldata %>%
  dplyr::select(PID, Group, PHQ,
                K.trans_thal)
levels(plotdata$Group) <- c("Low", "High")

gplot <- ggplot(plotdata, aes(Group, K.trans_thal, fill = Group, color = Group)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center') +
  #stat_summary(geom = "errorbar", fun.data = mean_se, width = 0.6, alpha = 0.1) +
  #stat_summary(geom = "errorbar", fun.min = function(z) { quantile(z,0.25) },
  #             fun = median, fun.max = function(z) { quantile(z,0.75) }, width = 0.6, alpha = 0.1) +
  stat_summary(geom = "errorbar", fun.min = median, fun = median, fun.max = median, width = 0.6, alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Depressive Symptom Severity", y = "Volume Transfer Constant (K.trans) (/min)") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(title = element_text(face = "bold"), axis.title = element_text(face="bold"))
gplot
ggsave("Figures/Group Comparisons - DCE-MRI (Thalamus).png", gplot, width = 4, height = 4)


rm(plotdataCr, plotdataW, plotdata, gplot)


