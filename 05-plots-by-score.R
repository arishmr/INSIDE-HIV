## GOAL: GENERATE SCATTER PLOTS OF BIOMARKERS BY DEPRESSION SCORES, AND ASSESS CORRELATIONS BETWEEN THESE VARIABLES

## Histogram of PHQ-9 scores

hist(alldata$PHQ)

plotdata <- alldata
levels(plotdata$Group)[levels(plotdata$Group)=="Low Depression Severity"] <- "Low"
levels(plotdata$Group)[levels(plotdata$Group)=="High Depression Severity"] <- "High"

ggplot(plotdata) +
  geom_dotplot(aes(x = PHQ))

ggplot(plotdata, aes(x = Group, y = PHQ, colour = Group)) +
  geom_dotplot(aes(fill = Group, colour = Group), binaxis = 'y',
               stackdir = 'center',
               binwidth = 0.5) +
  stat_summary(geom = "errorbar", fun.min = median, fun = median, fun.max = median, width = 0.6, alpha = 0.3) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Depressive Symptom Severity", y = "Total PHQ-9 Score") +
  theme(title = element_text(face = "bold"), axis.title = element_text(face="bold"))
ggsave("Figures/PHQ-9 Score by Group.png", width = 4, height = 4)


## Scatterplots of each biomarker vs depression score


#### MRS DATA

plotdata <- alldata %>%
  dplyr::select(PID, Group, PHQ,
         Asp:Cho, Asp.Cr:Cho.Cr)

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

gplot <- ggplot(plotdataCr, aes(PHQ, conc)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", se = T) +
  facet_wrap( ~ factor(plotdataCr$biomarker), scale="free") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Total PHQ-9 Score", y = "Concentration (Ratio to Creatine)") +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(panel.spacing.x = unit(1, "lines"),
        panel.spacing.y = unit(2, "lines")) +
  theme(title = element_text(face = "bold"), axis.title = element_text(face="bold"))
gplot

ggsave("Figures/PHQ-9 Score and MRS (Creatine).png", gplot, width = 6, height = 6)


plotdataW <- plotdata %>%
  filter(referencing == "Water Referenced")

gplot <- ggplot(plotdataW, aes(PHQ, conc)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", se = T) +
  facet_wrap( ~ factor(plotdataW$biomarker), scale="free") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Total PHQ-9 Score", y = "Concentration (Ratio to Water)") +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(panel.spacing.x = unit(1, "lines"),
        panel.spacing.y = unit(2, "lines")) +
  theme(title = element_text(face = "bold"), axis.title = element_text(face="bold"))
gplot

ggsave("Figures/PHQ-9 Score and MRS (Water).png", gplot, width = 6, height = 6)



#### DW-MRS DATA

plotdata <- alldata %>%
  dplyr::select(PID, Group, PHQ,
         ADC.NAA:ADC.Cho)

plotdata <- plotdata %>%
  gather(key = "biomarker",
         value = "conc",
         ADC.NAA:ADC.Cho,
         factor_key = T)

plotdata <- plotdata %>%
  mutate(biomarker = recode(biomarker,
                            ADC.NAA = "NAA + NAAG",
                            ADC.Cr = "Creatine (Cr + PCr)", ADC.Cho = "Choline (PCho + GPC)"))

gplot <- ggplot(plotdata, aes(PHQ, conc)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", se = T) +
  facet_wrap(factor(plotdata$biomarker), scale="free_y") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Total PHQ-9 Score", y = "Apparent Diffusion Coefficient (ADC) (\u00b5m\u00b3/ms)") +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(panel.spacing.x = unit(1, "lines"),
        panel.spacing.y = unit(2, "lines")) +
  theme(title = element_text(face = "bold"), axis.title = element_text(face="bold"))
gplot

ggsave("Figures/PHQ-9 Score and DW-MRS.png", gplot, width = 10, height = 4)



#### DCE-MRI DATA

plotdata <- alldata %>%
  dplyr::select(PID, Group, PHQ,
         K.trans_WB)

gplot <- ggplot(plotdata, aes(PHQ, K.trans_WB)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", se = T) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Total PHQ-9 Score", y = "Volume Transfer Constant (K.trans) (/min)") +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(panel.spacing.x = unit(1, "lines"),
        panel.spacing.y = unit(2, "lines")) +
  theme(title = element_text(face = "bold"), axis.title = element_text(face="bold"))
gplot

ggsave("Figures/PHQ-9 Score and DCE-MRI (Whole Brain).png", gplot, width = 4, height = 4)


plotdata <- alldata %>%
  dplyr::select(PID, Group, PHQ,
                K.trans_thal)

gplot <- ggplot(plotdata, aes(PHQ, K.trans_thal)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", se = T) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Total PHQ-9 Score", y = "Volume Transfer Constant (K.trans) (/min)") +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(panel.spacing.x = unit(1, "lines"),
        panel.spacing.y = unit(2, "lines")) +
  theme(title = element_text(face = "bold"), axis.title = element_text(face="bold"))
gplot

ggsave("Figures/PHQ-9 Score and DCE-MRI (Thalamus).png", gplot, width = 4, height = 4)

#cor.test(plotdata$PHQ, plotdata$K.trans_thal, method = "spearman")
#p.adjust(c(0.53, 0.13, 0.20), method = "fdr")

rm(gplot, plotdata, plotdataCr, plotdataW)
