## GOAL: GENERATE SCATTER PLOTS OF BLOOD BIOMARKERS BY DEPRESSION SCORES, AND ASSESS CORRELATIONS BETWEEN THESE VARIABLES

## Scatterplots of each biomarker vs depression score

plotdata <- blooddata
colnames(plotdata)
plotdata[,c(2:21,23)] <- apply(plotdata[,c(2:21,23)], 2, scale)

plotdata <- plotdata %>%
  gather(key = "biomarker",
         value = "conc",
         BDNF:YKL40,
         factor_key = T)

gplot <- ggplot(plotdata, aes(PHQ, conc)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", se = T) +
  facet_wrap( ~ factor(plotdata$biomarker), scale="free") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Total PHQ-9 Score", y = "Concentration (z-score)") +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(panel.spacing.x = unit(1, "lines"),
        panel.spacing.y = unit(2, "lines")) +
  theme(title = element_text(face = "bold"), axis.title = element_text(face="bold"))
gplot

ggsave("Figures/PHQ-9 Score and Blood Biomarkers.png", gplot, width = 9, height = 8)

rm(gplot, plotdata)
