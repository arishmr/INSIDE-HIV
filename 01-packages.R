# Define vector of package names
packages <- c("ggplot2", "RColorBrewer", "tidyverse", "psych", "summarytools", "htmlTable", "officer", "dplyr", "stringr", "lubridate", "gtsummary", "flextable", "reshape2", "lsr", "corrplot", "ppcor", "EnvStats", "Hmisc")

# Load packages
for (package in packages) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  } else {
    library(package, character.only = TRUE)
  }
}

rm(packages, package)

