library(plyr, include.only = c("rbind.fill"))

################################################################
## Correlations between blood biomarkers and PHQ-9
################################################################

set.seed(123)

data <- blooddata %>%
  gather(key = "biomarker",
         value = "conc",
         BDNF:YKL40)

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

write.csv(output, "Results/Biomarker x PHQ9 Correlations - Blood.csv", row.names = F)

rm(splitdata, data, output)