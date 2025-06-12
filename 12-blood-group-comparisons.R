####### GOAL: compare means of biomarkers across the two depression groups using Wilcoxon rank sum tests

library(plyr, include.only = c("rbind.fill"))

############
## BLOOD BIOMARKERS
############

## Create output dataframe 
output <- data.frame(
  W = numeric(0),
  p = numeric(0),
  stringsAsFactors = F
)

data <- blooddata %>%
  gather(key = "biomarker",
         value = "conc",
         BDNF:YKL40)

data <- data %>%
  group_by(biomarker) %>%
  group_split

output <- lapply(data, FUN = function(x) {
  test <- wilcox.test(conc ~ Group, data = x)
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

write.csv(output, "Results/Group Comparisons Non-Parametric - Blood Biomarkers.csv", row.names = F)

rm(data, output)
