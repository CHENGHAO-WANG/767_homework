# Homework 1 - R Programming

library(dplyr)
library(tidyr)

# set working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# load data
data <- read.csv("../data/camp_teach.csv")
head(data)

# Q2: Dataset Description

# count the number of subjects
length(unique(data$id))
# 695

# number and timing of measurement occasions
data.wide <- data %>%
  pivot_wider(id_cols = id, names_from = visitc,
     values_from = fdays, names_prefix = "fdays_")

head(data.wide)

ncol(data.wide) - 1  # number of measurement occasions
# 20

sort(unique(data$visitc))
# [1]   0   2   4  12  16  24  28  36  40  44  48  52  56  60  64  72  84  96 108
# [20] 120

# Q4: Descriptive Statistics

# rename variables for clarity
data.tidy <- data %>%
  rename(
    PostBD_FEV = POSFEV,
    gender = GENDER,
    ethnic = ETHNIC)

# create discrete outcome variable
data.tidy <- data.tidy %>%
  mutate(PreBD_FFB = ifelse(PREFF < 70, 1, 0))

# use labels for categorical variables
data.tidy <- data.tidy %>%
  mutate(
    PreBD_FFB = recode(as.character(PreBD_FFB),
      "1" = "obstruction",
      "0" = "normal",
      .default = as.character(PreBD_FFB)
    ),
    TG = recode(as.character(TG),
      "A" = "budesonide",
      "B" = "nedocromil",
      "C" = "placebo",
      .default = as.character(TG)
    ),
    gender = recode(as.character(gender),
      "m" = "male",
      "f" = "female",
      .default = as.character(gender)
    ),
    ethnic = recode(as.character(ethnic),
      "w" = "white",
      "b" = "black",
      "h" = "hispanic",
      "o" = "other",
      .default = as.character(ethnic)
    ),
    anypet = recode(as.character(anypet),
      "1" = "yes",
      "2" = "no",
      .default = as.character(anypet)
    ),
    woodstove = recode(as.character(woodstove),
      "1" = "yes",
      "2" = "no",
      .default = as.character(woodstove)
    ),
    dehumid = recode(as.character(dehumid),
      "1" = "yes",
      "2" = "no",
      "3" = "unknown",
      .default = as.character(dehumid)
    ),
    parent_smokes = recode(as.character(parent_smokes),
      "1" = "yes",
      "2" = "no",
      .default = as.character(parent_smokes)
    ),
    any_smokes = recode(as.character(any_smokes),
      "1" = "yes",
      "2" = "no",
      .default = as.character(any_smokes)
    )
  )


# filter tidy data for baseline visit
data.baseline <- data.tidy %>%
  filter(visitc == 0)

continuous.variables <- c("age_rz", "hemog", "wbc", "agehome", "PostBD_FEV")
discrete.variables <- c("TG", "gender", "ethnic", "anypet", "woodstove", "dehumid",
                        "parent_smokes", "any_smokes", "PreBD_FFB" )

# compute descriptive statistics (auto)
summarize_continuous <- function(df, vars) {
  df %>%
    summarise(across(all_of(vars), list(
      n_non_missing = ~sum(!is.na(.)),
      mean = ~mean(., na.rm = TRUE),
      Std = ~sd(., na.rm = TRUE),
      median = ~median(., na.rm = TRUE),
      IQR = ~IQR(., na.rm = TRUE),
      pct_missing = ~mean(is.na(.)) * 100
    ), .names = "{.col}__{.fn}")) %>%
    pivot_longer(
      cols = everything(),
      names_to = c("variable", "stat"),
      names_sep = "__"
    ) %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    mutate(across(where(is.numeric), ~round(., 2)))
}

summarize_discrete <- function(df, vars) {
  out <- lapply(vars, function(v) {
    x <- df[[v]]
    n_total <- length(x)
    x2 <- addNA(x, ifany = FALSE)

    if (n_total == 0) {
      return(tibble(
        variable = v,
        category = factor(levels(x2)),
        n = integer(0),
        pct = numeric(0)
      ))
    }

    tab <- table(x2, useNA = "no")
    tibble(
      variable = v,
      category = names(tab),
      n = as.integer(tab),
      pct = as.numeric(tab) / n_total * 100
    )
  })
  bind_rows(out) %>%
    mutate(across(where(is.numeric), ~round(., 2)))
}

continuous.summary <- summarize_continuous(data.baseline, continuous.variables)
discrete.summary <- summarize_discrete(data.baseline, discrete.variables)

continuous.labels <- c(
  age_rz = "Age at randomization (years)",
  hemog = "Hemoglobin (g/dL)",
  wbc = "White blood cell count (1000 cells/\\textmu{} L)",
  agehome = "Age of current home (years)",
  PostBD_FEV = "PostBD FEV1 (L)"
)

discrete.labels <- c(
  TG = "Treatment group",
  gender = "Gender",
  ethnic = "Ethnicity",
  anypet = "Any pets in home",
  woodstove = "Use a wood stove for heating/cooking",
  dehumid = "Use a dehumidifier",
  parent_smokes = "Either parent/partner smokes in home",
  any_smokes = "Anyone smokes in home",
  PreBD_FFB = "PreBD FFB (obstruction if PreBD FEV1/FVC ratio \\textless{} 70\\%)"
)

continuous.summary <- continuous.summary %>%
  mutate(variable = recode(variable, !!!continuous.labels, .default = variable))

discrete.summary <- discrete.summary %>%
  mutate(variable = recode(variable, !!!discrete.labels, .default = variable))

write.csv(continuous.summary, "./output/continuous_summary.csv", row.names = FALSE)
write.csv(discrete.summary, "./output/discrete_summary.csv", row.names = FALSE)

# save data.tidy for future use
write.csv(data.tidy, "../data/data_tidy.csv", row.names = FALSE)

# Q5: Limitations and Challenges

# missingess example
sum(is.na(data.tidy$hemog))/nrow(data.tidy) * 100
# [1] 86.96089