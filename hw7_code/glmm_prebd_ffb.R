if (!requireNamespace("lme4", quietly = TRUE)) {
    stop(
        "Package 'lme4' is required for GLMM models. ",
        "Install it with install.packages('lme4') and rerun this script.",
        call. = FALSE
    )
}

script_path <- normalizePath(
    sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]),
    winslash = "/",
    mustWork = TRUE
)
script_dir <- dirname(script_path)
project_dir <- dirname(script_dir)

output_dir <- file.path(script_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

dat <- read.csv(
    file.path(project_dir, "data", "data_tidy.csv"),
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
dat_model$TG <- factor(dat_model$TG)
dat_model$gender <- relevel(factor(dat_model$gender), ref = "female")
dat_model$ethnic <- relevel(factor(dat_model$ethnic), ref = "white")
dat_model$visitc <- as.numeric(dat_model$visitc)
dat_model$visitc_f <- factor(dat_model$visitc)

glmm_formula <- PreBD_FFB_binary ~
    TG * visitc_f + age_rz + gender + ethnic + (visitc | id)

fit_glmm <- lme4::glmer(
    formula = glmm_formula,
    data = dat_model,
    family = binomial(link = "logit"),
    nAGQ = 0,
    control = lme4::glmerControl(optimizer = "bobyqa")
)

coef_table <- as.data.frame(summary(fit_glmm)$coefficients)
parameter_estimates <- data.frame(
    parameter = row.names(coef_table),
    estimate = coef_table[["Estimate"]],
    std_error = coef_table[["Std. Error"]],
    statistic = coef_table[["z value"]],
    p_value = coef_table[["Pr(>|z|)"]],
    row.names = NULL,
    check.names = FALSE
)

write.csv(
    parameter_estimates,
    file.path(output_dir, "glmm_parameter_estimates.csv"),
    row.names = FALSE,
    quote = TRUE
)

variance_components <- as.data.frame(lme4::VarCorr(fit_glmm))

write.csv(
    variance_components,
    file.path(output_dir, "glmm_variance_components.csv"),
    row.names = FALSE,
    quote = TRUE
)
