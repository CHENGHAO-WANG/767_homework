################################
# Part 1: Descriptive Analysis 
################################

### 1. Descriptive statistics (10 points): Create a table presenting the mean and standard devia-
### tion of your continuous outcome at each time point, separately by group.

library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

dat <- read.csv("../data/data_tidy.csv", na.strings = c("", "NA"))

# group variable: TG
# time points: visitc
# continous outcome variable: PostBD_FEV

# First, turn missing rows into explicit NA
all_visit <- sort(unique(dat$visitc))
dat.full <- dat %>%
  group_by(id, TG) %>%
  complete(visitc = all_visit) %>%
  ungroup()

dat.continuous_summary <- dat.full %>%
  group_by(TG, visitc) %>%
  summarise(
    mean_PostBD_FEV = mean(PostBD_FEV, na.rm = TRUE),
    sd_PostBD_FEV = sd(PostBD_FEV, na.rm = TRUE),
    missing_pct = sum(is.na(PostBD_FEV)) / n() * 100
  )

View(dat.continuous_summary)

# too many missing values
# check the missingness at each time point
dat.continuous_missing <- dat.full %>%
  group_by(visitc) %>%
  summarise(missing_pct = sum(is.na(PostBD_FEV)) / n() * 100)

View(dat.continuous_missing)
#   visitc missing_pct
#    <int>       <dbl>
# 1      0       0.144
# 2      2       4.60
# 3      4       3.74
# 4     12       3.17
# 5     16       4.46
# 6     24       4.89
# 7     28       8.06
# 8     36       6.33
# 9     40       9.50
#10     44      98.6  
#11     48       9.93
#12     52      39.4
#13     56      80
#14     60      35.0
#15     64      97.3
#16     72      15.5
#17     84      17.6
#18     96      17.4
#19    108      28.8
#20    120      89.8

# Since RM ANOVA and MANOVA require complete data,
# we will filter out time points with too many missing values
visitc_to_include <- dat.continuous_missing %>%
  filter(missing_pct < 20) %>%
  pull(visitc)
dat.full.filtered <- dat.full %>%
  filter(visitc %in% visitc_to_include)

# then check how many subjects have no missing data across the included time points
dat.full.filtered.complete_subj <- dat.full.filtered %>%
  group_by(id) %>%
  summarise(missing = ifelse(any(is.na(PostBD_FEV)), 1, 0)) %>%
  filter(missing == 0)

nrow(dat.full.filtered.complete_subj)
# 471 subjects (out of 695)
# We will only include these 471 subjects in the analysis.

sort(visitc_to_include)
# [1]  0  2  4 12 16 24 28 36 40 48 72 84 96
# These time points will be included in the analysis.

length(visitc_to_include)
# 13 time points (out of 20)

id_to_include <- dat.full.filtered.complete_subj %>%
  pull(id)

# how many subjects in each treatment group are kept
dat %>%
  #filter(id %in% id_to_include) %>%
  mutate(id_included = ifelse(id %in% id_to_include, id, NA)) %>%
  group_by(TG) %>%
  summarise(n = n_distinct(id), n_included = n_distinct(id_included, na.rm =TRUE)) %>%
  mutate(include_pct = n_included/n * 100)

# filtering for complete cases
dat.filtered <- dat.full %>%
  filter(id %in% id_to_include, visitc %in% visitc_to_include)

# save the filtered data
write.csv(dat.filtered, "../data/complete_case.csv", row.names = FALSE)

# re-compute the mean and sd for the included subjects and time points
dat.filtered.continuous_summary <- dat.filtered %>%
  group_by(TG, visitc) %>%
  summarise(
    mean = mean(PostBD_FEV, na.rm = TRUE),
    sd = sd(PostBD_FEV, na.rm = TRUE)
  ) %>%
  ungroup()

# long to wide
dat.continuous_summary.wide <- dat.filtered.continuous_summary %>%
  pivot_wider(
    id_cols = visitc,
    names_from = TG,
    values_from = c(mean, sd),
    names_sep = "_"
  )

View(dat.continuous_summary.wide)

# arrange columns
group_suffix <- sort(unique(dat$TG))

ordered_cols <- unlist(lapply(group_suffix, function(g) paste0(c("mean_", "sd_"), g)))

dat.continuous_summary.wide <- dat.continuous_summary.wide %>%
  select(visitc, any_of(ordered_cols))


write.csv(dat.continuous_summary.wide, "./output/descriptive_summary.csv", row.names = FALSE)

### 2. Interaction plot (5 points): Create a plot showing the mean response over time for each
### group. Based on your plot, comment on whether the profiles appear parallel and whether the
### groups appear to differ in their overall level.

library(ggplot2)

ggplot(
  dat.filtered.continuous_summary,
  aes(
    x = factor(visitc, levels = sort(unique(visitc))),
    y = mean,
    group = factor(TG),
    color = factor(TG)
  )
) +
  geom_line() +
  geom_point() +
  labs(
    x = "Visit (Months)",
    y = "Mean PostBD FEV1",
    color = "Treatment Group"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "orange", "green"))

ggsave("./output/interaction_plot.png", width = 8, height = 5)
