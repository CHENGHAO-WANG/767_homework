
### Part 1: Model Specification
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

dat <- read.csv("../data/data_tidy.csv", na.strings = c("", "NA"))

View(dat)

# Missingness percentage by column
missingness_pct <- data.frame(
  column = names(dat),
  missing_pct = colMeans(is.na(dat)) * 100,
  row.names = NULL
)
missingness_pct <- missingness_pct[order(-missingness_pct$missing_pct), ]

print(missingness_pct)
View(missingness_pct)

#          column missing_pct
# 20           wbc  87.4535036
# 7          hemog  86.9608927
# 21       agehome  68.6438122
# 24       dehumid  68.2215743
# 23     woodstove  68.1512014
# 25 parent_smokes  68.1310948
# 26    any_smokes  68.1310948
# 22        anypet  68.1210415
# 15         POSPF   1.0053282
# 11         PREPF   0.8042626
# 18      POSFEVPP   0.4523977
# 19      POSFVCPP   0.4523977
# 16      PREFEVPP   0.3920780
# 17      PREFVCPP   0.3920780
# 12    PostBD_FEV   0.3719714
# 13        POSFVC   0.3719714
# 14         POSFF   0.3719714
# 8         PREFEV   0.3116518
# 9         PREFVC   0.3116518
# 10         PREFF   0.3116518
# 29     PreBD_FFB   0.3116518
# 1             TX   0.0000000
# 2             TG   0.0000000
# 3             id   0.0000000
# 4         age_rz   0.0000000
# 5         gender   0.0000000
# 6         ethnic   0.0000000
# 27        visitc   0.0000000
# 28         fdays   0.0000000

# group variable: TG
# time points: visitc
# continous outcome variable: PostBD_FEV
# other covariates: age_rz, gender, ethnic


# Plot deviation from planned measurement time
library(ggplot2)

days_per_month <- 365.25 / 12
timing_dev <- subset(dat, !is.na(visitc) & !is.na(fdays))
timing_dev$planned_days <- timing_dev$visitc * days_per_month
timing_dev$deviation_days <- timing_dev$fdays - timing_dev$planned_days

# Shared plot theme with larger text for report readability
plot_theme_large_text <- theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13),
    plot.title = element_text(size = 17, face = "bold")
  )

deviation_plot <- ggplot(
  timing_dev,
  aes(x = factor(visitc, levels = sort(unique(visitc))), y = deviation_days)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_boxplot(fill = "lightblue", color = "gray30", outlier.alpha = 0.25) +
  labs(
    # title = "Deviation from Planned Measurement Time",
    x = "Planned visit time (month)",
    y = "Actual visit time (day)"
  ) +
  plot_theme_large_text

print(deviation_plot)

if (!dir.exists("./output")) {
  dir.create("./output", recursive = TRUE)
}

ggsave(
  "./output/measurement_timing_deviation.png",
  deviation_plot,
  width = 7,
  height = 4.5,
  dpi = 300
)

# check the missingness at each time point
all_visit <- sort(unique(dat$visitc))
dat.full <- dat %>%
  group_by(id, TG) %>%
  complete(visitc = all_visit) %>%
  ungroup()

dat.continuous_missing <- dat.full %>%
  group_by(visitc) %>%
  summarise(missing_pct = sum(is.na(PostBD_FEV)) / n() * 100)

print(dat.continuous_missing)
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

### Visualization
# Mean trajectories by treatment group over visit month (visitc)
dat.mean.trajectory <- dat %>%
  filter(!is.na(TG), !is.na(visitc), !is.na(PostBD_FEV)) %>%
  group_by(TG, visitc) %>%
  summarise(mean_PostBD_FEV = mean(PostBD_FEV), .groups = "drop")

mean_trajectory_plot <- ggplot(
  dat.mean.trajectory,
  aes(
    x = visitc,
    y = mean_PostBD_FEV,
    color = factor(TG),
    group = factor(TG)
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = sort(unique(dat.mean.trajectory$visitc))) +
  labs(
    x = "Visit month",
    y = "Mean PostBD FEV1",
    color = "Treatment\ngroup"
  ) +
  plot_theme_large_text +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(mean_trajectory_plot)
ggsave(
  "./output/mean_trajectories_by_TG_over_visitc.png",
  mean_trajectory_plot,
  width = 7,
  height = 4.5,
  dpi = 300
)

visitc_to_include <- dat.continuous_missing %>%
  filter(missing_pct < 75) %>%
  pull(visitc)
dat.filtered <- dat %>%
  filter(visitc %in% visitc_to_include)

dat.mean.trajectory2 <- dat.filtered  %>%
  filter(!is.na(TG), !is.na(visitc), !is.na(PostBD_FEV)) %>%
  group_by(TG, visitc) %>%
  summarise(mean_PostBD_FEV = mean(PostBD_FEV), .groups = "drop")

mean_trajectory_plot2 <- ggplot(
  dat.mean.trajectory2,
  aes(
    x = visitc,
    y = mean_PostBD_FEV,
    color = factor(TG),
    group = factor(TG)
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = sort(unique(dat.mean.trajectory2$visitc))) +
  labs(
    x = "Visit month",
    y = "Mean PostBD FEV1",
    color = "Treatment\ngroup"
  ) +
  plot_theme_large_text +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(mean_trajectory_plot2)
ggsave(
  "./output/filtered_mean_trajectories_by_TG_over_visitc.png",
  mean_trajectory_plot2,
  width = 7,
  height = 4.5,
  dpi = 300
)
