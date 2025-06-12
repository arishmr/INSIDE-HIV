## Goal: calculate correlations for all blood proteins with all neuroimaging parameters
## Then plot these in a multi-panel heatmap with a facet for neuroimaging category

data <- alldata %>%
  dplyr::select(Asp:Cho.Cr,
                ADC.NAA:K.trans_thal,
                BDNF:YKL40)

output <- rcorr(as.matrix(data), type = "spearman")
#View(output$r)

output$r <- output$r %>% as.data.frame %>%
  dplyr::select(BDNF:YKL40)
output$r <- output$r[1:25,]

output$P <- output$P %>% as.data.frame %>%
  dplyr::select(BDNF:YKL40)
output$P <- output$P[1:25,]

output$n <- output$n %>% as.data.frame %>%
  dplyr::select(BDNF:YKL40)
output$n <- output$n[1:25,]

r <- as.data.frame(output$r)
r$var1 <- row.names(r)
r <- gather(r,
               key = "var2",
               value = "r",
               BDNF:YKL40)

p <- as.data.frame(output$P)
p$var1 <- row.names(p)
p <- gather(p,
            key = "var2",
            value = "p",
            BDNF:YKL40)


n <- as.data.frame(output$n)
n$var1 <- row.names(n)
n <- gather(n,
            key = "var2",
            value = "n",
            BDNF:YKL40)

output <- merge(r, p)
output <- merge(output, n)
rm(r, p, n, data)

output[,3] <- round(output[,3], 2)
#output[,4] <- round(output[,4], 3)

output$category <- ifelse(str_detect(output$var1, 'ADC'), "DW-MRS",
                          ifelse(str_detect(output$var1, '.Cr'), "MRS (Ratio to Creatine)",
                                ifelse(str_detect(output$var1, 'K.trans'), "DCE-MRI",
                                               "MRS (Ratio to Water)")))
output$category <- as.factor(output$category)

splitdata <- output %>%
  group_by(category) %>%
  group_split()
splitdata <- lapply(splitdata, FUN = function (x) {x <- x %>% dplyr::mutate(p.fdr = p.adjust(x$p, method = "fdr"))})
output <- dplyr::bind_rows(splitdata)

write.csv(output, "Results/Blood x Neuroimaging Correlations.csv", row.names = F)

output <- output %>%
  mutate(var1 = recode(var1,
                       Asp = "Aspartate", Cr = "Creatine", Cho = "Choline", Glc = "Glucose", Glu.Gln = "Glutamate\n + Glutamine", Ins = "Myo-Inositol", Lac = "Lactate", NAA = "NAA",
                       Asp.Cr = "Aspartate", Cho.Cr = "Choline", Glc.Cr = "Glucose", Glu.Gln.Cr = "Glutamate\n + Glutamine", Ins.Cr = "Myo-Inositol", Lac.Cr = "Lactate", NAA.Cr = "NAA",
                       ADC.Cho = "Choline ADC", ADC.Cr = "Creatine ADC", ADC.NAA = "NAA ADC",
                       K.trans_WB = "Whole Brain K.trans", K.trans_thal = "Thalamus K.trans"))

output <- output %>%
  filter(category != "MRS (Ratio to Water)")

plotdata <- output

ggplot(plotdata, aes(x = var1, y = var2, fill = r)) +
  geom_tile() +
  geom_text(aes(label=r), size = 3) +
  facet_grid(~ factor(category), scale="free", space = "free") +
  scale_fill_gradient(name = "Spearman's \U03C1", low = "darkgoldenrod1", high = "blueviolet", na.value = "white") +
  labs(x = "Neuroimaging Parameters",
       y = "Blood Proteins") +
  theme_classic() +
  theme(panel.spacing.x = unit(0.5, "lines")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title = element_text(face="bold")) +
  theme(legend.position = "top")
ggsave("Figures/Blood-Neuroimaging Correlations.png", width = 10, height = 10)


plotdata <- output %>%
  mutate(r = ifelse(p<0.05, r, NA))

ggplot(plotdata, aes(x = var1, y = var2, fill = r)) +
  geom_tile() +
  geom_text(aes(label=r), size = 3) +
  facet_grid(~ factor(category), scale="free", space = "free") +
  scale_fill_gradient(name = "Spearman's \U03C1", low = "darkgoldenrod1", high = "blueviolet", na.value = "white") +
  labs(x = "Neuroimaging Parameters",
       y = "Blood Proteins") +
  theme_classic() +
  theme(panel.spacing.x = unit(0.5, "lines")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title = element_text(face="bold")) +
  theme(legend.position = "top")
ggsave("Figures/Blood-Neuroimaging Correlations (Significant Only).png", width = 10, height = 10)


#plotdata <- output %>%
#  mutate(r = ifelse(p.fdr<0.05, r, NA))

#ggplot(plotdata, aes(x = var1, y = var2, fill = r)) +
#  geom_tile() +
#  geom_text(aes(label=r), size = 3) +
#  facet_grid(~ factor(category), scale="free", space = "free") +
#  scale_fill_gradient(name = "Spearman's \U03C1", low = "darkgoldenrod1", high = "blueviolet", na.value = "white") +
#  labs(x = "Neuroimaging Parameters",
#       y = "Blood Proteins") +
#  theme_classic() +
#  theme(panel.spacing.x = unit(0.5, "lines")) +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#  theme(axis.title = element_text(face="bold")) +
#  theme(legend.position = "top")
#ggsave("Figures/Blood-Neuroimaging Correlations (Significant Only, FDR-Corrected).png", width = 10, height = 10)


rm(output, plotdata, splitdata)
