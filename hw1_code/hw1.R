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
      "A" = "bud",
      "B" = "ned",
      "C" = "plbo",
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
    pivot_wider(names_from = stat, values_from = value)
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
  bind_rows(out)
}

continuous.summary <- summarize_continuous(data.baseline, continuous.variables)
discrete.summary <- summarize_discrete(data.baseline, discrete.variables)

write.csv(continuous.summary, "./output/continuous_summary.csv", row.names = FALSE)
write.csv(discrete.summary, "./output/discrete_summary.csv", row.names = FALSE)

