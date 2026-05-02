library(lme4)
library(ggplot2)

setwd("./hw7_code")

output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

dat <- read.csv(
    "../data/data_tidy.csv",
    na.strings = c("", "NA")
)

missingness_var <- "PreBD_FFB"
missingness_grid <- expand.grid(
    id = sort(unique(dat$id)),
    visitc = sort(unique(dat$visitc))
)
missingness_observed <- dat[c("id", "visitc", missingness_var)]
missingness_complete <- merge(
    missingness_grid,
    missingness_observed,
    by = c("id", "visitc"),
    all.x = TRUE,
    sort = TRUE
)
missingness_visitc <- sort(unique(missingness_complete$visitc))
missingness_by_visitc <- do.call(
    rbind,
    lapply(missingness_visitc, function(visit_value) {
        visit_values <- missingness_complete[
            missingness_complete$visitc == visit_value,
            missingness_var
        ]
        data.frame(
            visitc = visit_value,
            total_n = length(visit_values),
            missing_n = sum(is.na(visit_values)),
            missing_pct = 100 * mean(is.na(visit_values))
        )
    })
)

included_visitc <- missingness_by_visitc$visitc[
    missingness_by_visitc$missing_pct <= 75
]

required_vars <- c(
    "TG", "id", "age_rz", "gender", "ethnic", "visitc", "PreBD_FFB"
)

dat_model <- dat[dat$visitc %in% included_visitc, required_vars]
dat_model <- dat_model[complete.cases(dat_model), ]
dat_model <- dat_model[order(dat_model$id, dat_model$visitc), ]

dat_model$PreBD_FFB <- factor(dat_model$PreBD_FFB)
if (!all(c("normal", "obstruction") %in% levels(dat_model$PreBD_FFB))) {
    stop("PreBD_FFB must contain both 'normal' and 'obstruction'.", call. = FALSE)
}

dat_model$PreBD_FFB_binary <- as.integer(dat_model$PreBD_FFB == "obstruction")
dat_model$id <- factor(dat_model$id)
dat_model$TG <- relevel(factor(dat_model$TG), ref = "placebo")
dat_model$gender <- relevel(factor(dat_model$gender), ref = "female")
dat_model$ethnic <- relevel(factor(dat_model$ethnic), ref = "white")
dat_model$visitc <- as.numeric(dat_model$visitc)
dat_model$visitc_f <- factor(dat_model$visitc)

fixed_effects_formula <- PreBD_FFB_binary ~
    TG * visitc_f + age_rz + gender + ethnic
glmm_formula <- update(fixed_effects_formula, . ~ . + (visitc | id))
glmm_random_slope_formula <- update(
    fixed_effects_formula,
    . ~ . + (0 + visitc | id)
)
glmm_random_intercept_formula <- update(
    fixed_effects_formula,
    . ~ . + (1 | id)
)

fit_glmm <- lme4::glmer(
    formula = glmm_formula,
    data = dat_model,
    family = binomial(link = "logit"),
    nAGQ = 0,
    control = lme4::glmerControl(optimizer = "bobyqa")
)

fit_glmm_random_slope <- lme4::glmer(
    formula = glmm_random_slope_formula,
    data = dat_model,
    family = binomial(link = "logit"),
    nAGQ = 0,
    control = lme4::glmerControl(optimizer = "bobyqa")
)

fit_glmm_random_intercept <- lme4::glmer(
    formula = glmm_random_intercept_formula,
    data = dat_model,
    family = binomial(link = "logit"),
    nAGQ = 0,
    control = lme4::glmerControl(optimizer = "bobyqa")
)

fit_fixed_effects <- glm(
    formula = fixed_effects_formula,
    data = dat_model,
    family = binomial(link = "logit")
)

aic_comparison <- data.frame(
    model = c(
        "random_intercept_and_slope",
        "random_slope_only",
        "random_intercept_only",
        "fixed_effects_only"
    ),
    random_effects = c(
        "(visitc | id)",
        "(0 + visitc | id)",
        "(1 | id)",
        "none"
    ),
    df = c(
        attr(logLik(fit_glmm), "df"),
        attr(logLik(fit_glmm_random_slope), "df"),
        attr(logLik(fit_glmm_random_intercept), "df"),
        attr(logLik(fit_fixed_effects), "df")
    ),
    AIC = c(
        AIC(fit_glmm),
        AIC(fit_glmm_random_slope),
        AIC(fit_glmm_random_intercept),
        AIC(fit_fixed_effects)
    ),
    row.names = NULL
)
aic_comparison <- aic_comparison[order(aic_comparison$AIC), ]

write.csv(
    aic_comparison,
    file.path(output_dir, "glmm_aic_comparison.csv"),
    row.names = FALSE,
    quote = TRUE
)

coef_table <- as.data.frame(summary(fit_glmm)$coefficients)
p_values <- coef_table[["Pr(>|z|)"]]
parameter_estimates <- data.frame(
    parameter = row.names(coef_table),
    estimate = coef_table[["Estimate"]],
    std_error = coef_table[["Std. Error"]],
    statistic = coef_table[["z value"]],
    p_value = ifelse(
        p_values < 0.001,
        "<0.001",
        formatC(p_values, format = "f", digits = 3)
    ),
    row.names = NULL,
    check.names = FALSE
)

write.csv(
    parameter_estimates,
    file.path(output_dir, "glmm_parameter_estimates.csv"),
    row.names = FALSE,
    quote = TRUE
)

interaction_terms <- grep(
    ":visitc_f",
    names(lme4::fixef(fit_glmm)),
    value = TRUE
)
coef_covariance <- as.matrix(vcov(fit_glmm))
interaction_estimates <- lme4::fixef(fit_glmm)[interaction_terms]
interaction_covariance <- coef_covariance[interaction_terms, interaction_terms]
parallel_wald_chisq <- as.numeric(
    t(interaction_estimates) %*%
        solve(interaction_covariance, interaction_estimates)
)
parallel_p_value <- pchisq(
    parallel_wald_chisq,
    df = length(interaction_terms),
    lower.tail = FALSE
)
parallel_trajectory_test <- data.frame(
    hypothesis = "All TG-by-visit interaction coefficients are zero",
    df = length(interaction_terms),
    wald_chisq = parallel_wald_chisq,
    p_value = ifelse(
        parallel_p_value < 0.001,
        "<0.001",
        formatC(parallel_p_value, format = "f", digits = 3)
    ),
    row.names = NULL
)

write.csv(
    parallel_trajectory_test,
    file.path(output_dir, "glmm_parallel_trajectory_test.csv"),
    row.names = FALSE,
    quote = TRUE
)

variance_components <- as.matrix(lme4::VarCorr(fit_glmm)$id)

write.csv(
    variance_components,
    file.path(output_dir, "glmm_variance_components.csv"),
    row.names = TRUE,
    quote = TRUE
)

gee_parameter_estimates <- read.csv(
    file.path(
        "..",
        "hw6_code",
        "output",
        "gee_parameter_estimates.csv"
    )
)
glmm_estimates <- setNames(
    parameter_estimates$estimate,
    parameter_estimates$parameter
)
estimate_comparison <- data.frame(
    parameter = gee_parameter_estimates$parameter,
    GEE_estimate = gee_parameter_estimates$independence,
    GLMM_estimate = as.numeric(glmm_estimates[gee_parameter_estimates$parameter]),
    row.names = NULL
)

write.csv(
    estimate_comparison,
    file.path(output_dir, "glmm_gee_independence_estimate_comparison.csv"),
    row.names = FALSE,
    quote = TRUE
)

random_effect_predictions <- lme4::ranef(fit_glmm)$id
random_effect_predictions <- data.frame(
    id = row.names(random_effect_predictions),
    random_intercept = random_effect_predictions[["(Intercept)"]],
    random_slope = random_effect_predictions[["visitc"]],
    row.names = NULL
)

write.csv(
    random_effect_predictions,
    file.path(output_dir, "glmm_random_effect_predictions.csv"),
    row.names = FALSE,
    quote = TRUE
)

random_intercept_breaks <- hist(
    random_effect_predictions$random_intercept,
    breaks = "FD",
    plot = FALSE
)$breaks
random_intercept_plot <- ggplot2::ggplot(
    random_effect_predictions,
    ggplot2::aes(x = random_intercept)
) +
    ggplot2::geom_histogram(
        breaks = random_intercept_breaks,
        fill = "steelblue",
        color = "white"
    ) +
    ggplot2::labs(
        title = "Predicted Random Intercepts by Subject",
        x = "Predicted random intercept",
        y = "Count"
    ) +
    ggplot2::theme_bw()
ggplot2::ggsave(
    filename = file.path(output_dir, "glmm_random_intercept_histogram.png"),
    plot = random_intercept_plot,
    width = 9,
    height = 6,
    dpi = 100
)

random_slope_breaks <- hist(
    random_effect_predictions$random_slope,
    breaks = "FD",
    plot = FALSE
)$breaks
random_slope_plot <- ggplot2::ggplot(
    random_effect_predictions,
    ggplot2::aes(x = random_slope)
) +
    ggplot2::geom_histogram(
        breaks = random_slope_breaks,
        fill = "darkorange",
        color = "white"
    ) +
    ggplot2::labs(
        title = "Predicted Random Slopes by Subject",
        x = "Predicted random slope for visitc",
        y = "Count"
    ) +
    ggplot2::theme_bw()
ggplot2::ggsave(
    filename = file.path(output_dir, "glmm_random_slope_histogram.png"),
    plot = random_slope_plot,
    width = 9,
    height = 6,
    dpi = 100
)

set.seed(767)
subject_profiles <- dat_model[
    !duplicated(dat_model$id),
    c("id", "TG", "age_rz", "gender", "ethnic")
]
subject_profiles <- subject_profiles[order(subject_profiles$TG, subject_profiles$id), ]
# Keep seeded subject sampling independent of the TG reference level.
sample_tg <- sort(unique(as.character(subject_profiles$TG)))
selected_subject_ids <- unlist(
    lapply(sample_tg, function(tg_value) {
        subject_ids <- subject_profiles$id[subject_profiles$TG == tg_value]
        sample(subject_ids, size = 4)
    }),
    use.names = FALSE
)
selected_subject_profiles <- subject_profiles[
    subject_profiles$id %in% selected_subject_ids,
]
selected_subject_profiles <- selected_subject_profiles[
    order(selected_subject_profiles$TG, selected_subject_profiles$id),
]

selected_subject_probabilities <- do.call(
    rbind,
    lapply(seq_len(nrow(selected_subject_profiles)), function(row_index) {
        subject_profile <- selected_subject_profiles[row_index, ]
        data.frame(
            id = subject_profile$id,
            TG = subject_profile$TG,
            age_rz = subject_profile$age_rz,
            gender = subject_profile$gender,
            ethnic = subject_profile$ethnic,
            visitc = sort(unique(dat_model$visitc))
        )
    })
)
selected_subject_probabilities$id <- factor(
    selected_subject_probabilities$id,
    levels = levels(dat_model$id)
)
selected_subject_probabilities$TG <- factor(
    selected_subject_probabilities$TG,
    levels = levels(dat_model$TG)
)
selected_subject_probabilities$gender <- factor(
    selected_subject_probabilities$gender,
    levels = levels(dat_model$gender)
)
selected_subject_probabilities$ethnic <- factor(
    selected_subject_probabilities$ethnic,
    levels = levels(dat_model$ethnic)
)
selected_subject_probabilities$visitc_f <- factor(
    selected_subject_probabilities$visitc,
    levels = levels(dat_model$visitc_f)
)
selected_subject_probabilities$predicted_probability <- predict(
    fit_glmm,
    newdata = selected_subject_probabilities,
    type = "response"
)

write.csv(
    selected_subject_probabilities[
        c("id", "TG", "visitc", "predicted_probability")
    ],
    file.path(output_dir, "glmm_selected_subject_probabilities.csv"),
    row.names = FALSE,
    quote = TRUE
)

selected_subject_plot_data <- selected_subject_probabilities
selected_subject_plot_data$id <- droplevels(selected_subject_plot_data$id)
selected_subject_plot_data$TG <- droplevels(selected_subject_plot_data$TG)
selected_subject_plot <- ggplot2::ggplot(
    selected_subject_plot_data,
    ggplot2::aes(
        x = visitc,
        y = predicted_probability,
        color = id,
        group = id
    )
) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::facet_wrap(ggplot2::vars(TG), ncol = 1) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
        title = "Subject-Specific Predicted Probabilities by Treatment Group",
        x = "Visit",
        y = "Predicted probability",
        color = "Subject id"
    ) +
    ggplot2::theme_bw()
ggplot2::ggsave(
    filename = file.path(output_dir, "glmm_selected_subject_probabilities.png"),
    plot = selected_subject_plot,
    width = 9,
    height = 9,
    dpi = 100
)
