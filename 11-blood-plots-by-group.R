################### GOAL: Plot blood biomarker variables by group

plotdata <- blooddata
colnames(plotdata)
plotdata[,c(2:21,23)] <- apply(plotdata[,c(2:21,23)], 2, scale)

plotdata$Group <- factor(plotdata$Group)
plotdata$Group <- relevel(plotdata$Group, ref = "Low Depression Severity")
levels(plotdata$Group) <- c("Low", "High")

plotdata <- plotdata %>%
  gather(key = "biomarker",
         value = "conc",
         BDNF:YKL40,
         factor_key = T)

gplot <- ggplot(plotdata, aes(Group, conc, fill = Group, color = Group)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center') +
  stat_summary(geom = "errorbar", fun.min = median, fun = median, fun.max = median, width = 0.6, alpha = 0.3) +
  facet_wrap(factor(plotdata$biomarker), scale="free_y") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Depressive Symptom Severity", y = "Concentration (z-score)") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(title = element_text(face = "bold"), axis.title = element_text(face="bold"))
gplot
ggsave("Figures/Group Comparisons - Blood Biomarkers.png", gplot, width = 8, height = 8)

rm(gplot, plotdata)