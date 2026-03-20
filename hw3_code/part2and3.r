
library(tidyverse)
library(nlme)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("./hw3_code")

dat <- read.csv("../data/data_tidy.csv", na.strings = c("", "NA"))

# fit glm using data where missingness < 75% at each time points
all_visitc <- dat %>%
    filter(!(visitc %in% c(44, 56, 64, 120))) %>%
    distinct(visitc) %>%
    arrange(visitc) %>%
    pull(visitc)

dat2 <- dat %>%
    filter(!(visitc %in% c(44, 56, 64, 120))) %>%
    mutate(visitc = as.numeric(visitc)) %>%
    arrange(id, visitc) %>%
    group_by(id) %>%
    complete(visitc = all_visitc) %>%
    fill(age_rz, gender, ethnic, TG, .direction = "downup") %>%
    ungroup() %>%
    mutate(
        id = factor(id),
        TG = factor(TG),
        gender = relevel(factor(gender), ref = "female"),
        ethnic = relevel(factor(ethnic), ref = "white"),
        visitc_f = factor(visitc)
    ) %>%
    arrange(id, visitc, TG)

###################################################
# warning: the runtime is too long! Don't run it.
fit.reml <- gls(
    PostBD_FEV ~ 0 + TG + TG:visitc + age_rz + gender + ethnic,
    data = dat2,
    method = "REML",
    correlation = corSymm(form = ~ 1 | id),
    weights = varIdent(form = ~ 1),
    na.action = na.omit
)
###################################################

# We take a subset of time points to roughly view the
# correlation structure.
visitc_keep <- c(0, 12, 24, 36, 48, 60, 72, 84, 96)

dat3 <- dat2 %>%
    filter(visitc %in% visitc_keep)

fit.reml.rough <- gls(
    PostBD_FEV ~ 0 + TG + TG:visitc + age_rz + gender + ethnic,
    data = dat3,
    method = "REML",
    correlation = corSymm(form = ~ 1 | id),
    weights = varIdent(form = ~ 1),
    na.action = na.omit
)

(corr.rough <- fit.reml.rough$modelStruct$corStruct)

# Export estimated response correlation matrix from fit.reml.rough
if (!dir.exists("./output")) {
    dir.create("./output", recursive = TRUE)
}

write_matrix_outputs <- function(
    mat, csv_path, tex_path,
    row_name_col = "visitc",
    header_label = "Visit") {
    mat <- as.matrix(mat)
    mat <- matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    if (is.null(rownames(mat))) {
        rownames(mat) <- as.character(seq_len(nrow(mat)))
    }
    if (is.null(colnames(mat))) {
        colnames(mat) <- as.character(seq_len(ncol(mat)))
    }

    mat_df <- as.data.frame(round(mat, 3)) %>%
        rownames_to_column(var = row_name_col)

    write.csv(
        mat_df,
        csv_path,
        row.names = FALSE,
        quote = FALSE
    )

    latex_lines <- c(
        "\\small",
        paste0("\\begin{tabular}{l", paste(rep("r", ncol(mat)), collapse = ""), "}"),
        "\\hline",
        paste0(header_label, " & ", paste(colnames(mat), collapse = " & "), " \\\\"),
        "\\hline"
    )

    for (i in seq_len(nrow(mat))) {
        row_text <- paste(
            c(
                rownames(mat)[i],
                formatC(mat[i, ], format = "f", digits = 3)
            ),
            collapse = " & "
        )
        latex_lines <- c(latex_lines, paste0(row_text, " \\\\"))
    }

    latex_lines <- c(latex_lines, "\\hline", "\\end{tabular}")

    writeLines(
        latex_lines,
        tex_path
    )
}

corr_mat_rough <- corMatrix(corr.rough)[[1]]
dimnames(corr_mat_rough) <- list(visitc_keep, visitc_keep)
write_matrix_outputs(
    corr_mat_rough,
    "./output/fit_reml_rough_correlation_matrix.csv",
    "./output/fit_reml_rough_correlation_matrix.tex"
)

# The correlation structure is close to Toeplitz, AR(1), Exponential
sort(unique(dat2$visitc))
# Since the times points are not equally spaced
# in `dat2`, the "more complete" data, we will assume
# Exponentail correlation structure.

# Nevertheless, the correlation with fixed lag
# slightly decreases as time increases.
# Want to check variance heterogeneity
# using a subset of data with less time points,
# as we have more parameters to estimate.

# visitc_keep <- c(0, 2, 12, 48, 72)
visitc_keep <- c(0, 24, 48, 72)

dat4 <- dat2 %>%
    filter(visitc %in% visitc_keep)

fit.reml.rough2 <- gls(
    PostBD_FEV ~ 0 + TG + TG:visitc + age_rz + gender + ethnic,
    data = dat4,
    method = "REML",
    correlation = corSymm(form = ~ 1 | id),
    weights = varIdent(form = ~ 1 | visitc_f),
    na.action = na.omit
)

(corr.rough2 <- fit.reml.rough2$modelStruct$corStruct)
(var.rough2 <- fit.reml.rough2$modelStruct$varStruct)

getVarCov(fit.reml.rough2, individual = "1")

# Export estimated response correlation matrix from fit.reml.rough2
corr_mat_rough2 <- corMatrix(corr.rough2)[[1]]
dimnames(corr_mat_rough2) <- list(visitc_keep, visitc_keep)
write_matrix_outputs(
    corr_mat_rough2,
    "./output/fit_reml_rough2_correlation_matrix.csv",
    "./output/fit_reml_rough2_correlation_matrix.tex"
)

# Export estimated response covariance matrix from fit.reml.rough2 (individual 1)
cov_mat_rough2 <- as.matrix(getVarCov(fit.reml.rough2, individual = "1"))
if (nrow(cov_mat_rough2) == length(visitc_keep)) {
    dimnames(cov_mat_rough2) <- list(visitc_keep, visitc_keep)
} else {
    dimnames(cov_mat_rough2) <- list(seq_len(nrow(cov_mat_rough2)), seq_len(ncol(cov_mat_rough2)))
}
write_matrix_outputs(
    cov_mat_rough2,
    "./output/fit_reml_rough2_covariance_matrix.csv",
    "./output/fit_reml_rough2_covariance_matrix.tex"
)


### model fit
fit.reml.exp <- gls(
    PostBD_FEV ~ 0 + TG + TG:visitc + age_rz + gender + ethnic,
    data = dat2,
    method = "REML",
    correlation = corExp(form = ~ visitc | id),
    weights = varIdent(form = ~ 1 | visitc_f),
    na.action = na.omit
)

(corr <- fit.reml.exp$modelStruct$corStruct)

(V_1 <- getVarCov(fit.reml.exp, individual = "1"))
(R_1 <- cov2cor(V_1))

# to fit the model with parsimonious variance structure,
# create a new time variable
# This is necessary, because baseline time is 0,
# which makes optimization numerically difficult
# for variance structure including power of time.
time_shift <- 0.1
dat2$visitc_shifted <- dat2$visitc + time_shift

fit.reml.exp.varconst <- gls(
    PostBD_FEV ~ 0 + TG + TG:visitc + age_rz + gender + ethnic,
    data = dat2,
    method = "REML",
    correlation = corExp(form = ~ visitc | id),
    weights = varIdent(form = ~ 1),
    na.action = na.omit
)

fit.reml.exp.varpower <- gls(
    PostBD_FEV ~ 0 + TG + TG:visitc + age_rz + gender + ethnic,
    data = dat2,
    method = "REML",
    correlation = corExp(form = ~ visitc | id),
    weights = varPower(form = ~ visitc_shifted),
    na.action = na.omit
)

fit.reml.exp.varconstpower <- gls(
    PostBD_FEV ~ 0 + TG + TG:visitc + age_rz + gender + ethnic,
    data = dat2,
    method = "REML",
    correlation = corExp(form = ~ visitc | id),
    weights = varConstPower(form = ~ visitc_shifted),
    na.action = na.omit
)

fit.reml.exp.varexp <- gls(
    PostBD_FEV ~ 0 + TG + TG:visitc + age_rz + gender + ethnic,
    data = dat2,
    method = "REML",
    correlation = corExp(form = ~ visitc | id),
    weights = varExp(form = ~ visitc_shifted),
    na.action = na.omit
)

fit.reml.exp.varconstprop <- gls(
    PostBD_FEV ~ 0 + TG + TG:visitc + age_rz + gender + ethnic,
    data = dat2,
    method = "REML",
    correlation = corExp(form = ~ visitc | id),
    weights = varConstProp(form = ~ visitc_shifted),
    na.action = na.omit
)

aic_model_compare <- AIC(
    fit.reml.exp, fit.reml.exp.varconst,
    fit.reml.exp.varpower,
    fit.reml.exp.varconstpower,
    fit.reml.exp.varexp, fit.reml.exp.varconstprop
)
bic_model_compare <- BIC(
    fit.reml.exp, fit.reml.exp.varconst,
    fit.reml.exp.varpower,
    fit.reml.exp.varconstpower,
    fit.reml.exp.varexp, fit.reml.exp.varconstprop
)

print(aic_model_compare)
print(bic_model_compare)

aic_model_compare_df <- as.data.frame(aic_model_compare) %>%
    rownames_to_column(var = "model")
bic_model_compare_df <- as.data.frame(bic_model_compare) %>%
    rownames_to_column(var = "model")

model_compare_aic_bic_df <- aic_model_compare_df %>%
    select(model, df, AIC) %>%
    left_join(
        bic_model_compare_df %>% select(model, BIC),
        by = "model"
    )

write.csv(
    aic_model_compare_df,
    "./output/fit_reml_exp_model_compare_aic.csv",
    row.names = FALSE,
    quote = FALSE
)
write.csv(
    bic_model_compare_df,
    "./output/fit_reml_exp_model_compare_bic.csv",
    row.names = FALSE,
    quote = FALSE
)
write.csv(
    model_compare_aic_bic_df,
    "./output/fit_reml_exp_model_compare_aic_bic.csv",
    row.names = FALSE,
    quote = FALSE
)

anova(fit.reml.exp, fit.reml.exp.varconst,
fit.reml.exp.varpower,
fit.reml.exp.varconstpower,
fit.reml.exp.varexp, fit.reml.exp.varconstprop)
# getVarCov(fit.reml.exp, individual = "1")

### baseline intercept

fit.ml.exp <- gls(
    PostBD_FEV ~ 0 + TG + TG:visitc + age_rz + gender + ethnic,
    data = dat2,
    method = "ML",
    correlation = corExp(form = ~ visitc | id),
    weights = varIdent(form = ~ 1 | visitc_f),
    na.action = na.omit
)

fit.ml.c.exp <- gls(
    PostBD_FEV ~ 1 + TG:visitc + age_rz + gender + ethnic,
    data = dat2,
    method = "ML",
    correlation = corExp(form = ~ visitc | id),
    weights = varIdent(form = ~ 1 | visitc_f),
    na.action = na.omit
)

coef(summary(fit.ml.exp))
coef(summary(fit.ml.c.exp))

# Likelihood ratio test: reduced model (common intercept) vs full model
# This also shows AIC and BIC
anova(fit.ml.c.exp, fit.ml.exp)

# Conclude that a parsimonious model that constrains
# basline mean to be equal is reasonable.

dat2 <- dat2 %>% arrange(id, visitc, TG)

fit.reml.c.exp <- gls(
    PostBD_FEV ~ 1 + TG:visitc + age_rz + gender + ethnic,
    data = dat2,
    method = "REML",
    correlation = corExp(form = ~ visitc | id),
    weights = varIdent(form = ~ 1 | visitc_f),
    na.action = na.omit
)

# Whether to model unequal covariance matrices across treatment group?
fit.reml.c.exp.unequal_corr <- gls(
    PostBD_FEV ~ 1 + TG:visitc + age_rz + gender + ethnic,
    data = dat2,
    method = "REML",
    correlation = corExp(form = ~ visitc | TG/id),
    weights = varIdent(form = ~ 1 | visitc_f),
    na.action = na.omit
)

fit.reml.c.exp$modelStruct$corStruct
fit.reml.c.exp.unequal_corr$modelStruct$corStruct

getVarCov(fit.reml.c.exp.unequal_corr, individual = "1")
getVarCov(fit.reml.c.exp.unequal_corr, individual = "9")
getVarCov(fit.reml.c.exp, individual = "1")

anova(fit.reml.c.exp, fit.reml.c.exp.unequal_corr)
coef(summary(fit.reml.c.exp))
coef(summary(fit.reml.c.exp.unequal_corr))

# Unfortunately, different correlation across groups 
# is not supported by gls()

# different weights (variances)
fit.reml.c.exp.unequal_var <- gls(
    PostBD_FEV ~ 1 + TG:visitc + age_rz + gender + ethnic,
    data = dat2,
    method = "REML",
    correlation = corExp(form = ~ visitc | id),
    weights = varIdent(form = ~ 1 | visitc_f * TG),
    na.action = na.omit
)

anova(fit.reml.c.exp, fit.reml.c.exp.unequal_var)

# AIC are close, whereas LRT p-val is significant.
# Now check the variance in each group

fit.reml.c.exp$modelStruct$varStruct
fit.reml.c.exp.unequal_var$modelStruct$varStruct

# getVarCov(fit.reml.c.exp, individual = "1")
# this doesn't give you covariance matrix over all time points.
# This only gives the covariance matrix of subject with id = 1.
# Because there are missing values.
# We have 16 time points, while subject 1 only has 15 time points,
# and then its covariance matrix is 15*15.

# equal variance across groups
getVarCov(fit.reml.c.exp, individual = "1")
fit.reml.c.exp$modelStruct$varStruct
fit.reml.c.exp$modelStruct

# unequal variance across groups
getVarCov(fit.reml.c.exp.unequal_var, individual = "1")
fit.reml.c.exp.unequal_var$modelStruct$varStruct

# sigma() estimates standard deviations
# we can thereby construct the covariance matrices.
sigma(fit.reml.c.exp)

sigma(fit.reml.c.exp.unequal_var)

################
# diagonal variances of fit.reml.c.exp
# 1) Baseline standard deviation from sigma()
baseline_sd <- sigma(fit.reml.c.exp)

# 2) Standard deviation multipliers by time from varStruct
visit_levels <- levels(dat2$visitc_f)
sd_multiplier <- setNames(rep(1, length(visit_levels)), visit_levels)

var_struct_coef <- coef(
    fit.reml.c.exp$modelStruct$varStruct,
    unconstrained = FALSE
)

coef_names <- names(var_struct_coef)
sd_multiplier[coef_names] <- unname(var_struct_coef)

# 3) Square the standard deviations to get variances
# Build a data.frame(time, variance)
df_var.fit.reml.c.exp <- data.frame(
    time = as.numeric(visit_levels),
    variance = (baseline_sd * as.numeric(sd_multiplier))^2
) %>%
    arrange(time)

df_var.fit.reml.c.exp

################
# diagonal variances of fit.reml.c.exp.unequal_var
# 1) Baseline standard deviation from sigma()
baseline_sd_unequal_var <- sigma(fit.reml.c.exp.unequal_var)

# 2) Standard deviation multipliers by (time, TG) from varStruct
var_struct_coef_unequal_var <- coef(
    fit.reml.c.exp.unequal_var$modelStruct$varStruct,
    unconstrained = FALSE
)

tg_time_grid <- dat2 %>%
    distinct(visitc_f, TG) %>%
    mutate(
        TG = as.character(TG),
        time = as.numeric(as.character(visitc_f)),
        key = paste0(as.character(visitc_f), "*", TG)
    )

sd_multiplier_unequal_var <- setNames(rep(1, nrow(tg_time_grid)), tg_time_grid$key)

coef_names <- names(var_struct_coef_unequal_var)
sd_multiplier_unequal_var[coef_names] <- var_struct_coef_unequal_var

# 3) Square the standard deviations to get variances
# Build a data.frame(TG, time, variance)
df_var.fit.reml.c.exp.unequal_var <- tg_time_grid %>%
    mutate(
        variance = (baseline_sd_unequal_var * as.numeric(sd_multiplier_unequal_var))^2
    ) %>%
    select(TG, time, variance) %>%
    arrange(TG, time)

df_var.fit.reml.c.exp.unequal_var

################
# compare variance estimates between fit.reml.c.exp and fit.reml.c.exp.unequal_var

df_var.fit.reml.c.exp <- df_var.fit.reml.c.exp %>%
    mutate(TG = "average") %>%
    mutate(model = "homogeneous variances")

df_var.fit.reml.c.exp.unequal_var_labeled <- df_var.fit.reml.c.exp.unequal_var %>%
    mutate(model = "group-specific variances")

df_var.compare <- bind_rows(
    df_var.fit.reml.c.exp,
    df_var.fit.reml.c.exp.unequal_var_labeled
) %>% arrange(TG, time, model)

df_var.compare$model <- as.factor(df_var.compare$model)
df_var.compare$TG <- as.factor(df_var.compare$TG)

variance_compare_plot <- ggplot(
    df_var.compare,
    aes(
        x = time,
        y = variance,
        color = TG,
        linetype = model,
        group = interaction(TG, model)
    )
) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    labs(
        x = "Time",
        y = "Variance",
        color = "Treatment\ngroup",
        linetype = "Model"#,
        #title = "Variance Estimates by Model",
        #subtitle = "Single-plot comparison across TG and model"
    ) +
    theme_bw()

print(variance_compare_plot)
ggsave(
  "./output/variance_compare.png",
  variance_compare_plot,
  width = 7,
  height = 4.5,
  dpi = 300
)

##############################
# The final model fit:
fit.reml.c.exp.unequal_var

coef(summary(fit.reml.c.exp.unequal_var))

saveRDS(fit.reml.c.exp.unequal_var, file = "model_fit.rds")
fit.reml.c.exp.unequal_var <- readRDS("model_fit.rds")

### fitted mean trajectories vs observed mean trajectories

# Baseline summaries from dat2 (which includes NA values)
baseline_dat <- dat2 %>%
    filter(visitc == 0, !is.na(TG))

baseline_age_by_tg <- baseline_dat %>%
    group_by(TG) %>%
    summarise(
        mean_age_rz = mean(age_rz, na.rm = TRUE),
        .groups = "drop"
    )

baseline_gender_prop_by_tg <- baseline_dat %>%
    filter(!is.na(gender)) %>%
    count(TG, gender, name = "n_gender") %>%
    group_by(TG) %>%
    complete(gender = levels(dat2$gender), fill = list(n_gender = 0)) %>%
    mutate(prop_gender = n_gender / sum(n_gender)) %>%
    ungroup()

baseline_ethnic_prop_by_tg <- baseline_dat %>%
    filter(!is.na(ethnic)) %>%
    count(TG, ethnic, name = "n_ethnic") %>%
    group_by(TG) %>%
    complete(ethnic = levels(dat2$ethnic), fill = list(n_ethnic = 0)) %>%
    mutate(prop_ethnic = n_ethnic / sum(n_ethnic)) %>%
    ungroup()

baseline_age_by_tg
baseline_gender_prop_by_tg
baseline_ethnic_prop_by_tg

# Fitted mean at each visitc for each TG:
# use baseline mean(age_rz | TG) and baseline marginal proportions
# of gender and ethnic within each TG as weights.
pred_grid <- expand_grid(
    TG = levels(dat2$TG),
    visitc = sort(unique(dat2$visitc)),
    gender = levels(dat2$gender),
    ethnic = levels(dat2$ethnic)
) %>%
    left_join(baseline_age_by_tg, by = "TG") %>%
    left_join(
        baseline_gender_prop_by_tg %>% select(TG, gender, prop_gender),
        by = c("TG", "gender")
    ) %>%
    left_join(
        baseline_ethnic_prop_by_tg %>% select(TG, ethnic, prop_ethnic),
        by = c("TG", "ethnic")
    ) %>%
    mutate(
        TG = factor(TG, levels = levels(dat2$TG)),
        gender = factor(gender, levels = levels(dat2$gender)),
        ethnic = factor(ethnic, levels = levels(dat2$ethnic)),
        age_rz = mean_age_rz,
        weight = coalesce(prop_gender, 0) * coalesce(prop_ethnic, 0)
    )

pred_grid$pred <- predict(fit.reml.c.exp.unequal_var, newdata = pred_grid)
# weight is used to marginalize model predictions over
# baseline gender/ethnic composition within each TG
fitted_mean_by_tg_visit <- pred_grid %>%
    group_by(TG, visitc) %>%
    summarise(
        fitted_mean = sum(pred * weight, na.rm = TRUE),
        weight_sum = sum(weight, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    arrange(TG, visitc)

if (any(abs(fitted_mean_by_tg_visit$weight_sum - 1) > 1e-8)) {
    warning("Weighted prediction rows do not sum to 1 for some TG/visitc combinations.")
}

fitted_mean_by_tg_visit

# Observed mean at each visitc in each TG
observed_mean_by_tg_visit <- dat2 %>%
    group_by(TG, visitc) %>%
    summarise(
        observed_mean = mean(PostBD_FEV, na.rm = TRUE),
        n_obs = sum(!is.na(PostBD_FEV)),
        .groups = "drop"
    ) %>%
    arrange(TG, visitc)

observed_mean_by_tg_visit

# Combined fitted vs observed means in one plot
mean_compare_by_tg_visit <- bind_rows(
    fitted_mean_by_tg_visit %>%
        transmute(TG, visitc, mean_value = fitted_mean, series = "Fitted"),
    observed_mean_by_tg_visit %>%
        transmute(TG, visitc, mean_value = observed_mean, series = "Observed")
) %>%
    mutate(series = factor(series, levels = c("Observed", "Fitted")))

mean_compare_plot <- ggplot(
    mean_compare_by_tg_visit,
    aes(
        x = visitc,
        y = mean_value,
        color = TG,
        linetype = series,
        shape = series,
        group = interaction(TG, series)
    )
) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    geom_point(size = 1.5, na.rm = TRUE) +
    labs(
        x = "Time (month)",
        y = "Mean PostBD_FEV",
        color = "Treatment\ngroup",
        linetype = "Trajectory",
        shape = "Trajectory"
    ) +
    theme_bw()

print(mean_compare_plot)
ggsave(
    "./output/fitted_observed_mean_compare.png",
    mean_compare_plot,
    width = 7,
    height = 4.5,
    dpi = 300
)

### Covariance matrices

# Construct estimated covariance matrices by treatment group using:
# covariance(time_i, time_j | TG) = sd(time_i, TG) * sd(time_j, TG) * corr(time_i, time_j)
baseline_sd_cov <- sigma(fit.reml.c.exp.unequal_var)
var_struct_coef_cov <- coef(
    fit.reml.c.exp.unequal_var$modelStruct$varStruct,
    unconstrained = FALSE
)

tg_time_grid_cov <- dat2 %>%
    distinct(visitc_f, TG) %>%
    mutate(
        TG = as.character(TG),
        time = as.numeric(as.character(visitc_f)),
        key = paste0(as.character(visitc_f), "*", TG)
    ) %>%
    arrange(TG, time)

sd_multiplier_cov <- setNames(rep(1, nrow(tg_time_grid_cov)), tg_time_grid_cov$key)
sd_multiplier_cov[names(var_struct_coef_cov)] <- var_struct_coef_cov

tg_time_grid_cov <- tg_time_grid_cov %>%
    mutate(
        sd = baseline_sd_cov * as.numeric(sd_multiplier_cov)
    ) %>%
    arrange(TG, time)

tg_time_grid_cov

visitc_cov_levels <- sort(unique(tg_time_grid_cov$time))

cor_value_cov <- coef(
    fit.reml.c.exp.unequal_var$modelStruct$corStruct,
    unconstrained = FALSE
)
cor_template_cov <- corExp(
    value = cor_value_cov,
    form = ~ visitc | id,
    nugget = FALSE,
    fixed = TRUE
)
cor_template_cov <- Initialize(
    cor_template_cov,
    data = data.frame(
        visitc = visitc_cov_levels,
        id = factor(rep("template_id", length(visitc_cov_levels)))
    )
)

(corr_mat_cov <- corMatrix(cor_template_cov))
# verify the correlation matrix we computed
cov2cor(getVarCov(fit.reml.c.exp.unequal_var, individual = "1"))


dimnames(corr_mat_cov) <- list(visitc_cov_levels, visitc_cov_levels)

for (tg_name in unique(tg_time_grid_cov$TG)) {
    sd_vec_tg <- tg_time_grid_cov %>%
        filter(TG == tg_name) %>%
        arrange(time) %>%
        pull(sd)

    cov_mat_tg <- (sd_vec_tg %o% sd_vec_tg) * corr_mat_cov
    dimnames(cov_mat_tg) <- list(visitc_cov_levels, visitc_cov_levels)

    write_matrix_outputs(
        cov_mat_tg,
        paste0("./output/fit_reml_c_exp_unequal_var_covariance_matrix_", tg_name, ".csv"),
        paste0("./output/fit_reml_c_exp_unequal_var_covariance_matrix_", tg_name, ".tex")
    )
}



#################
# parameter estimation

# Default output by R
model_based_coef_table <- coef(summary(fit.reml.c.exp.unequal_var)) %>%
    as.data.frame() %>%
    rownames_to_column(var = "term") %>%
    rename(
        estimate = Value,
        std_error = `Std.Error`,
        t_value = `t-value`,
        p_value = `p-value`
    )

# The following function computes SE that matches SAS
robust.cov <- function(u) {
  form <- formula(u)
  mf <- model.frame(form, getData(u))
  Xmat <- model.matrix(form, mf)
  resid_vec <- residuals(u, type = "response")

  # Use character IDs to avoid factor-indexing mismatch in getVarCov().
  group_vec <- as.character(getGroups(u))
  ids <- unique(group_vec)

  p <- ncol(Xmat)
  XtVinvX <- matrix(0, nrow = p, ncol = p)
  meat <- matrix(0, nrow = p, ncol = p)

  # Accumulate model-based and sandwich "meat" terms subject-by-subject.
  for (id_i in ids) {
    idx <- which(group_vec == id_i)
    Xi <- Xmat[idx, , drop = FALSE]
    ei <- resid_vec[idx]
    Vi <- as.matrix(getVarCov(u, individual = id_i, type = "marginal"))
    Vinv_i <- solve(Vi)

    XtVinvX <- XtVinvX + t(Xi) %*% Vinv_i %*% Xi
    meat <- meat + t(Xi) %*% Vinv_i %*% (ei %o% ei) %*% Vinv_i %*% Xi
  }

  Sig.model <- solve(XtVinvX)
  Sig.robust <- Sig.model %*% meat %*% Sig.model
  se.robust <- sqrt(diag(Sig.robust))
  se.model <- sqrt(diag(Sig.model))

  list(
    Sig.model = Sig.model,
    se.model = se.model,
    Sig.robust = Sig.robust,
    se.robust = se.robust
  )
}


model_based_coef_table

### we use Between-Within method for dof estimation
nsub       <- length(unique(dat2$id))
nobs       <- sum(!is.na(dat2$PostBD_FEV))
# number of between-subject effects
nbetween   <- 6
# number of within-subject effects
nwithin    <- 3

beta.model <- coef(fit.reml.c.exp.unequal_var)
vcov.model <- robust.cov(fit.reml.c.exp.unequal_var)
saveRDS(vcov.model, file = "beta_vcov.rds")
vcov.model <- readRDS("beta_vcov.rds")
se.model   <- vcov.model$se.model
t.model    <- beta.model/se.model
ddfm.model <- c(rep(nsub - nbetween, nbetween), rep(nobs - nsub - nwithin, nwithin))
p.model    <- pt(-abs(t.model), ddfm.model)

coef.model_based <- data.frame(
    Estimate = beta.model,
    StdError = se.model,
    t = t.model,
    ddfm = ddfm.model,
    p.val = p.model
) %>%
 tibble::rownames_to_column(var = "term") %>%
 mutate(
     p.val = if_else(p.val < 0.001, "<0.001", sprintf("%.3f", p.val))
 )
coef.model_based



write.csv(coef.model_based, "./output/coef_model_based.csv",
    row.names = FALSE, quote = FALSE)


######################
# sandwich estimates
library(clubSandwich)

dat2 %>%
  group_by(id) %>%
  summarise(n_obs = sum(!is.na(PostBD_FEV)), .groups = "drop") %>%
  filter(n_obs == 1)

# This verifies that there are singleton subjects
# (only 1 usable observations after removing na)

################################
# The following will trigger an error due to singleton subjects.
# Don't run it!
# Without singleton subjects, it should have given sandwich estimates.

V_sand <- vcovCR(fit.reml.c.exp.unequal_var, 
    cluster = getGroups(fit.reml.c.exp.unequal_var), type = "CR0")

sandwich_coef_raw <- coef_test(
    fit.reml.c.exp.unequal_var,
    vcov = V_sand,
    test = "Satterthwaite"
) %>%
    as.data.frame()

################################

beta.robust <- coef(fit.reml.c.exp.unequal_var)
vcov.robust <- vcov.model
se.robust   <- vcov.robust$se.robust
t.robust    <- beta.robust / se.robust
ddfm.robust <- c(rep(nsub - nbetween, nbetween), rep(nobs - nsub - nwithin, nwithin))
p.robust    <- pt(-abs(t.robust), ddfm.robust)

coef.robust <- data.frame(
    Estimate = beta.robust,
    StdError = se.robust,
    t = t.robust,
    ddfm = ddfm.robust,
    p.val = p.robust
) %>%
 tibble::rownames_to_column(var = "term") %>%
 mutate(
     p.val = if_else(p.val < 0.001, "<0.001", sprintf("%.3f", p.val))
 )
coef.robust

write.csv(coef.robust, "./output/coef_robust.csv",
    row.names = FALSE, quote = FALSE)



######################
# export combined coefficient table for LaTeX

escape_latex_term <- function(x) {
    gsub("_", "\\_", x, fixed = TRUE)
}

format_coef_rows <- function(df) {
    term <- escape_latex_term(df$term)
    estimate <- sprintf("%.4f", as.numeric(df$Estimate))
    std_error <- sprintf("%.4f", as.numeric(df$StdError))
    t_value <- sprintf("%.3f", as.numeric(df$t))
    ddfm <- as.character(as.integer(round(as.numeric(df$ddfm))))
    p_value <- as.character(df$p.val)

    paste0(
        term, " & ", estimate, " & ", std_error, " & ",
        t_value, " & ", ddfm, " & ", p_value, " \\\\"
    )
}

coef_model_vs_robust_tex <- c(
    "\\begin{tabular}{lrrrrl}",
    "\\hline",
    "\\multicolumn{6}{l}{Model-based} \\\\",
    "\\hline",
    "Term & Estimate & Std. Error & t & df & p-value \\\\",
    "\\hline",
    format_coef_rows(coef.model_based),
    "\\hline",
    "\\multicolumn{6}{l}{Robust} \\\\",
    "\\hline",
    format_coef_rows(coef.robust),
    "\\hline",
    "\\end{tabular}"
)

writeLines(coef_model_vs_robust_tex, "./output/coef_model_vs_robust.tex")

# how many individuals were included in the analysis
nlevels(getGroups(fit.reml.c.exp.unequal_var))
