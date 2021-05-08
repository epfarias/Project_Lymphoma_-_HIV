#' ---
#' title: "A correlation analysis of clinical variables of TCGA-DLBCLpatients"
#' output: 
#'   html_document: 
#'     default
#'   github_document: 
#'     df_print: paged
#'     html_preview: FALSE
#'     keep_html: TRUE
#'   pdf_document:
#'     latex_engine: xelatex
#' knit: (function(inputFile, encoding) {
#'   rmarkdown::render(inputFile, encoding = encoding, output_format = "all") })     
#' ---
#' 
#' This project contains a pipeline of clinical analysis of the Cancer Genome Atlas Diffuse Large B-Cell Lymphoma (TCGA-DLBCL) data of patients, from [Genomic Data Commons Data Portal]().
#' 
#' Previously, we presented [an exploratory preprocessing analysis](1.preprocessing.md). In this section, Chi-squared test is applied to compare two or more proportions of categorical variables and T-student test to compare the means of numeric ones regardind the levels of 'Overall_Survival_Status'. The Hypoteis test is performed and p-value indicates the strength of evidence in supportting the null hypothesis.
#' 
#' 

#' 
#' 
## ----message=FALSE, warning=FALSE, paged.print=FALSE--------------------------
# Set the packages of interest
packages = c("tidyverse","skimr","finalfit","rstatix", "ggpubr","GGally", "plotly")

# if a package is installed, it will be loaded
# otherwise, the missing package(s) will be installed and loaded
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

suppressMessages(library("tidyverse"))
rm(packages)
setwd(".")

#' 
#' ## 1. Importing data
#' 
## ----message=FALSE, warning=FALSE, paged.print=FALSE--------------------------
dlbcl_clin <- read_csv("DBLC_clinical/dlbcl_clin.csv")

#' 
#' ## 2. Taming data 
#' 
## -----------------------------------------------------------------------------
dlbcl_clin <- dlbcl_clin %>%
  mutate_if(is.character, as.factor) %>%
  mutate(patient_id = as.character(patient_id),
         age = as.integer(age),
         year_of_diagnosis = as.integer(year_of_diagnosis))

# check 
glimpse(dlbcl_clin)

#' 
#' ## 3.The dependent variable
#' 
## -----------------------------------------------------------------------------
# Check the number of levels. If greater than 2, it thas to run a ordinal logistic regression presents only two levels (otherwise, it has to run a ordinal logistic regression)
table(dlbcl_clin$treatment_type, exclude = NULL)

#' 
#' ## 4. Numeric variables vs. treatment type
#' 
#' Correlation matrix - graphic visualization
#' 
## ----message=FALSE, warning=FALSE---------------------------------------------
cols_numeric <- dlbcl_clin %>% 
  select_if(is.numeric) %>%
  names

dlbcl_clin_numeric <- dlbcl_clin %>%
                      select(one_of(c(cols_numeric, "treatment_type")))  

levels(dlbcl_clin_numeric$treatment_type) <- c("Pharma","Radi")

ggpairs(dlbcl_clin_numeric, columns = cols_numeric, 
        title="Correlation matrix",               
        mapping= aes(colour = treatment_type), 
        upper = list(combo = wrap("box_no_facet", alpha=0.1), 
                     continuous = wrap("cor", size = 2, alignPercent = 0.8)),
        lower = list(continuous = wrap("smooth", alpha = 0.3, size=0.2) )) +
        theme(panel.background = element_rect(color = "black", size=0.5, fill="white"),
          panel.grid.major = element_blank()) 
        

#' 
#' Run multiple T-tests on treatment type
#' 
#' Transform the data into long format
#' 
## -----------------------------------------------------------------------------
# Put all variables in the same column except `over_surv_stt`, the grouping variable

levels(dlbcl_clin_numeric$treatment_type) <- c("Pharma","Radia")

# Convert to Tidyverse
dlbcl_clin_numeric.long <- dlbcl_clin_numeric %>%
  pivot_longer(-treatment_type, names_to = "variables", values_to = "value")
dlbcl_clin_numeric.long <- dlbcl_clin_numeric.long[!is.na(dlbcl_clin_numeric.long$value), ]
dlbcl_clin_numeric.long$value.log <- log2(dlbcl_clin_numeric.long$value+1)

# OR
# kirc_clin_numeric.long <- kirc_clin_numeric %>% 
#   gather(key = 'variables', value = 'value', -over_surv_stt, na.rm = TRUE) %>%
#     mutate(value.log = log2(kirc_clin_numeric.long$value+1))

dlbcl_clin_numeric.long %>% sample_n(6)

#' 
#' Group the data by variables and compare treatment type groups
#' 
#' Adjust the p-values and add significance levels
#' 
## -----------------------------------------------------------------------------
stat.test <- dlbcl_clin_numeric.long %>%
  group_by(variables) %>%
  t_test(value ~ treatment_type) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test

#' 
#' 
## -----------------------------------------------------------------------------
# Create the plot on logscale
myplot <- ggboxplot(
  dlbcl_clin_numeric.long, x = "treatment_type", y = "value.log",
  fill = "treatment_type", palette = "npg", legend = "none", 
  ggtheme = theme_pubr(border = TRUE)
  ) +
  facet_wrap(~variables)

# Add statistical test p-values
# OBS: different p-values over vaule vs. log.value!! 
stat.test <- stat.test %>% add_xy_position(x = "treatment_type")
myplot + stat_pvalue_manual(stat.test, label = "p.adj.signif")

#' 
#' 
## -----------------------------------------------------------------------------
# Group the data by variables and do a graph for each variable
graphs <- dlbcl_clin_numeric.long %>%
  group_by(variables) %>%
  doo(
    ~ggboxplot(
      data =., x = "treatment_type", y = "value",
      fill = "treatment_type", palette = "npg", legend = "none",
      ggtheme = theme_pubr()
      )  +
      geom_jitter(width = 0.05, alpha = 0.2, color = "orange"), 
    result = "plots"
  )
graphs

#' 
#' 
## -----------------------------------------------------------------------------
 #"Add statitistical tests to each corresponding plot"
 #variables <- graphs$variables
 #for(i in 1:length(variables)){
  # graph.i <- graphs$plots[[i]] + 
  #    labs(title = variables[i]) +
    #stat_pvalue_manual(stat.test[i, ], label = "p.adj.signif")
   #print(graphs.i)
 #}
# Error in print(graph.i) : objeto 'graph.i' não encontrado

#' 
#' 
## -----------------------------------------------------------------------------
# ggplot(kirc_clin, aes(age, fill= over_surv_stt)) +
#   geom_histogram(bins = 15, position = "dodge")
# t.test(kirc_clin$age ~ kirc_clin$over_surv_stt) 
# 
# ggplot(kirc_clin, aes(year_diagnose, fill= over_surv_stt)) +
#   geom_histogram(bins = 15, position = "dodge")
# t.test(kirc_clin$year_diagnose ~ kirc_clin$over_surv_stt) 
# 
# ggplot(kirc_clin, aes(x=over_surv_stt, y=disease_free_mth)) +
#   geom_boxplot(width = .5) +
#   geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
# t.test(kirc_clin$disease_free_mth ~ kirc_clin$over_surv_stt)
# 
# ggplot(kirc_clin, aes(x=over_surv_stt, y=frac_genome_alter)) +
#   geom_boxplot(width = .5) +
#   geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
# t.test(kirc_clin$frac_genome_alter ~ kirc_clin$over_surv_stt)
# 
# ggplot(kirc_clin, aes(x=over_surv_stt, y=long_dim)) +
#   geom_boxplot(width = .5) +
#   geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
# t.test(kirc_clin$long_dim ~ kirc_clin$over_surv_stt)
# 
# ggplot(kirc_clin, aes(x=over_surv_stt, y=mutation_cnt)) +
#   geom_boxplot(width = .5) +
#   geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
# t.test(kirc_clin$mutation_cnt ~ kirc_clin$over_surv_stt)
# 
# ggplot(kirc_clin, aes(x=over_surv_stt, y=over_surv_mth)) +
#   geom_boxplot(width = .5) +
#   geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
# t.test(kirc_clin$over_surv_mth ~ kirc_clin$over_surv_stt)
# 
# ggplot(kirc_clin, aes(x=over_surv_stt, y=short_dim)) +
#   geom_boxplot(width = .5) +
#   geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
# t.test(kirc_clin$short_dim ~ kirc_clin$over_surv_stt)
# 
# ggplot(kirc_clin, aes(x=over_surv_stt, y=second_long_dim)) +
#   geom_boxplot(width = .5) +
#   geom_jitter(width = 0.05, alpha = 0.2, color = "orange")
# t.test(kirc_clin$second_long_dim ~ kirc_clin$over_surv_stt)

#' 
#' Summary for continuous explanatory variables 
#' use a parametric or non-parametric test?? 
#' 
## -----------------------------------------------------------------------------

explanatory_num <- dlbcl_clin %>%
  select(-treatment_type) %>%
  select_if(is.numeric) %>%
  names
dependent <- 'treatment_type'

table_num <- dlbcl_clin %>%
  summary_factorlist(dependent, explanatory_num, p=TRUE, 
                     add_dependent_label=TRUE,  na_include = TRUE)

knitr::kable(table_num, row.names=FALSE, align=c("l", "l", "r", "r", "r"))


#' 
## -----------------------------------------------------------------------------
# Correlation Matrix
# Pearson's (normal distribution) or Spearman (not-normal) correlations
corr_num <- dlbcl_clin %>%
     select_if(is.numeric) %>%
     drop_na()

# Check the correlation between variables to exclude the higly correlated
cor_matrix <- cor(corr_num, method = "spearman")
cor_matrix <- round(cor_matrix, 2)
cor_matrix

#' 
#' ## 5. Categorical variables vs. treatment_type
#' 
#' Tabulation and chi-square test
#' 
## ----warning=FALSE------------------------------------------------------------
 t_metas_stg <- table(dlbcl_clin$aa_clin_stg, dlbcl_clin$treatment_type, exclude = NULL)
 t_metas_stg <- addmargins(round(100*prop.table(t_metas_stg)))
 t_metas_stg
 chisq.test(x = dlbcl_clin$aa_clin_stg, y = dlbcl_clin$treatment_type) 

#' 
#' Summary for chategorical explanatory variables
#' Chi-squared warnings will be generated when the expected count in any cell is less than 5.
#' 
## -----------------------------------------------------------------------------

summary(dlbcl_clin)

explanatory_char <- dlbcl_clin %>%
  select(-treatment_type) %>%
  select_if(is.factor) %>%
  names

dependent <-  'treatment_type'

table_char <- dlbcl_clin %>%
  summary_factorlist(dependent, explanatory_char, p=TRUE, 
                     add_dependent_label=TRUE,  na_include = TRUE)

knitr::kable(table_char, row.names=FALSE, align=c("l", "l", "r", "r", "r"))


# Droping levels with narrow distributions -> check warnings ()
# Group some levels or drop one (NULL = 'level') when grouping is not possible 

dlbcl_clin2 <- dlbcl_clin %>%
     mutate(ann_arbor_clinical_stage = fct_collapse(ann_arbor_clinical_stage, 'Stage I-II' = c('Stage I','Stage II'), 'Stage III-IV' = c('Stage III','Stage IV')))
                                                       
dlbcl_clin2 <- dlbcl_clin2 %>%
     mutate(race = fct_recode(race, NULL = 'Asian'))


table_char2 <- dlbcl_clin2 %>%
  summary_factorlist(dependent, explanatory_char, p=TRUE, 
                     add_dependent_label=TRUE,  na_include = TRUE)

knitr::kable(table_char2, row.names=FALSE, align=c("l", "l", "r", "r", "r"))

summary(dlbcl_clin2)

#' 
#' ## 6. saving dataset for regression model
#' 
## -----------------------------------------------------------------------------
dlbcl_glm <- dlbcl_clin2 
write_csv(dlbcl_glm, path = "data/dlbcl_glm.csv")

#' 
#' ## Further analysis
#' 
#' - [A logistic regression analysis](3.logistic_regression.md) of each clinical variable weight.
#' 
## -----------------------------------------------------------------------------
sessionInfo()

#' 
