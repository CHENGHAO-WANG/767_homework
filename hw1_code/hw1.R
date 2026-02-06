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
  mutate(PreBD_FFB = ifelse(PREFF < 70, "obstruction", "normal"))    

# filter tidy data for baseline visit
data.baseline <- data.tidy %>%
  filter(visitc == 0)

continuous.variables <- c("age_rz", "hemog", "wbc", "age_home", "PostBD_FEV")
discrete.variables <- c("TG", "gender", "ethnic", "anypet", "woodstove", "dehumid",
                        "parent_smokes", "any_smokes", "PreBD_FFB" )

# compute descriptive statistics
