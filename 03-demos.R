demodata <- alldata %>%
  dplyr::select(c(Group,
    age, sex, gender, ethnicity, edu, smoking, alcohol, recdrugs,
    curr.antidep, ever.antidep, curr.therapy, ever.therapy, ViralLoad, CD4)
    )

cor.test(alldata$edu, alldata$PHQ, method = "spearman")
cor.test(alldata$CD4, alldata$PHQ, method = "spearman")


## Create gt_summary table by depression status

tab1 <- tbl_summary(demodata,
                    by = Group,
                    type = list(age ~ "continuous",
                                sex ~ "categorical",
                                gender ~ "categorical",
                                #sexual.orient ~ "categorical",
                                ethnicity ~ "categorical",
                                edu ~ "continuous",
                                smoking ~ "categorical",
                                alcohol ~ "categorical",
                                recdrugs ~ "categorical",
                                curr.antidep ~ "categorical",
                                ever.antidep ~ "categorical",
                                curr.therapy ~ "categorical",
                                ever.therapy ~ "categorical",
                                ViralLoad ~ "categorical",
                                CD4 ~ "continuous"
                                #ART ~ "categorical",
                                #Antidepressant ~ "categorical"
                                ),
                    statistic = list(
                      all_continuous() ~ "{median} ({p25}, {p75})",
                      all_categorical() ~ "{n} ({p}%)"
                    ),
                    digits = all_continuous() ~ 1,
                    label = list(
                      age ~ "Age (years)",
                      sex ~ "Sex",
                      gender ~ "Gender",
                      #sexual.orient ~ "Sexual Orientation",
                      ethnicity ~ "Ethnicity",
                      edu ~ "Education (years)",
                      smoking ~ "Smoking",
                      alcohol ~ "Alcohol Use",
                      recdrugs ~ "Recreational Drug Use",
                      curr.antidep ~ "Currently Taking Antidepressants",
                      ever.antidep ~ "Ever Taken Antidepressants",
                      curr.therapy ~ "Currently Receiving Psychotherapy",
                      ever.therapy ~ "Ever Received Psychotherapy",
                      ViralLoad ~ "Viral Load",
                      CD4 ~ "CD4 Count"
                      #ART ~ "ART Regimen",
                      #Antidepressant ~ "Antidepressant Regimen"
                    ),
                    missing_text = "(Missing)"
) %>%
  add_p() %>%
  add_overall() %>%
  add_stat_label() %>%
  add_significance_stars()

## Save gt_summary table to Word doc

tab1.print <- tab1 %>% as_flex_table()
save_as_docx(tab1.print, 
             path = "Results/Demographics by Depression Status.docx",
             pr_section = prop_section(page_size(orient = "landscape")))

rm(tab1, tab1.print, demodata)
