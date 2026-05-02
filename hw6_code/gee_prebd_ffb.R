library(stats)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("./hw6_code")

if (!requireNamespace("geepack", quietly = TRUE)) {
    stop(
        "Package 'geepack' is required for GEE models. ",
        "Install it with install.packages('geepack') and rerun this script.",
        call. = FALSE
    )
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
        "Package 'ggplot2' is required for plotting. ",
        "Install it with install.packages('ggplot2') and rerun this script.",
        call. = FALSE
    )
}

output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

dat <- read.csv("../data/data_tidy.csv", na.strings = c("", "NA"))

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

write.csv(
    missingness_by_visitc,
    file.path(output_dir, "prebd_ffb_missingness_by_visitc.csv"),
    row.names = FALSE,
    quote = TRUE
)

logit_cap <- 10
subject_tg <- unique(dat[c("id", "TG")])
logit_grid <- merge(
    missingness_grid,
    subject_tg,
    by = "id",
    all.x = TRUE,
    sort = TRUE
)
logit_observed <- dat[c("id", "visitc", missingness_var)]
logit_complete <- merge(
    logit_grid,
    logit_observed,
    by = c("id", "visitc"),
    all.x = TRUE,
    sort = TRUE
)
logit_groups <- unique(logit_complete[c("TG", "visitc")])
logit_groups <- logit_groups[order(logit_groups$TG, logit_groups$visitc), ]
logit_by_tg_visitc <- do.call(
    rbind,
    lapply(seq_len(nrow(logit_groups)), function(row_index) {
        group_values <- logit_complete[
            logit_complete$TG == logit_groups$TG[row_index] &
                logit_complete$visitc == logit_groups$visitc[row_index],
            missingness_var
        ]
        observed_values <- group_values[!is.na(group_values)]
        obstruction_n <- sum(observed_values == "obstruction")
        observed_n <- length(observed_values)
        obstruction_prob <- obstruction_n / observed_n
        logit_obstruction_prob <- log(obstruction_prob / (1 - obstruction_prob))
        logit_obstruction_prob[is.infinite(logit_obstruction_prob)] <-
            sign(logit_obstruction_prob[is.infinite(logit_obstruction_prob)]) *
            logit_cap

        data.frame(
            TG = logit_groups$TG[row_index],
            visitc = logit_groups$visitc[row_index],
            observed_n = observed_n,
            obstruction_n = obstruction_n,
            obstruction_prob = obstruction_prob,
            logit_obstruction_prob = logit_obstruction_prob
        )
    })
)

write.csv(
    logit_by_tg_visitc,
    file.path(output_dir, "prebd_ffb_logit_by_tg_visitc.csv"),
    row.names = FALSE,
    quote = TRUE
)

logit_by_tg_visitc$TG <- relevel(factor(logit_by_tg_visitc$TG), ref = "placebo")
logit_plot <- ggplot2::ggplot(
    logit_by_tg_visitc,
    ggplot2::aes(
        x = visitc,
        y = logit_obstruction_prob,
        color = TG,
        group = TG
    )
) +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::geom_point(size = 2.5, na.rm = TRUE) +
    ggplot2::labs(
        title = "Observed Logit of Obstruction by Treatment Group and Visit",
        x = "Visit",
        y = "logit(P(PreBD_FFB = obstruction))",
        color = "TG"
    ) +
    ggplot2::theme_bw()
ggplot2::ggsave(
    filename = file.path(output_dir, "prebd_ffb_logit_by_tg_visitc.png"),
    plot = logit_plot,
    width = 9,
    height = 6,
    dpi = 100
)

included_visitc <- missingness_by_visitc$visitc[
    missingness_by_visitc$missing_pct <= 75
]
logit_by_tg_visitc_le75 <- logit_by_tg_visitc[
    logit_by_tg_visitc$visitc %in% included_visitc,
]

write.csv(
    logit_by_tg_visitc_le75,
    file.path(output_dir, "prebd_ffb_logit_by_tg_visitc_missingness_le75.csv"),
    row.names = FALSE,
    quote = TRUE
)

logit_le75_plot <- ggplot2::ggplot(
    logit_by_tg_visitc_le75,
    ggplot2::aes(
        x = visitc,
        y = logit_obstruction_prob,
        color = TG,
        group = TG
    )
) +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::geom_point(size = 2.5, na.rm = TRUE) +
    ggplot2::labs(
        title = "Observed Logit by Treatment Group and Visit, Missingness <= 75%",
        x = "Visit",
        y = "logit(P(PreBD_FFB = obstruction))",
        color = "TG"
    ) +
    ggplot2::theme_bw()
ggplot2::ggsave(
    filename = file.path(output_dir, "prebd_ffb_logit_by_tg_visitc_missingness_le75.png"),
    plot = logit_le75_plot,
    width = 9,
    height = 6,
    dpi = 100
)

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

gee_formula <- PreBD_FFB_binary ~ TG * visitc_f + age_rz + gender + ethnic

fit_independence <- geepack::geeglm(
    formula = gee_formula,
    id = id,
    data = dat_model,
    family = binomial(link = "logit"),
    corstr = "independence"
)

fit_ar1 <- geepack::geeglm(
    formula = gee_formula,
    id = id,
    waves = visitc,
    data = dat_model,
    family = binomial(link = "logit"),
    corstr = "ar1"
)

independence_se_comparison <- data.frame(
    parameter = names(coef(fit_independence)),
    `model-based SE` = sqrt(diag(fit_independence$geese$vbeta.naiv)),
    `sandwich SE` = sqrt(diag(fit_independence$geese$vbeta)),
    row.names = NULL,
    check.names = FALSE
)

write.csv(
    independence_se_comparison,
    file.path(output_dir, "gee_independence_se_comparison.csv"),
    row.names = FALSE,
    quote = TRUE
)

interaction_terms <- grep(
    ":visitc_f",
    names(coef(fit_independence)),
    value = TRUE
)
coef_cov_sandwich <- fit_independence$geese$vbeta
dimnames(coef_cov_sandwich) <- list(
    names(coef(fit_independence)),
    names(coef(fit_independence))
)
interaction_estimates <- coef(fit_independence)[interaction_terms]
interaction_cov <- coef_cov_sandwich[interaction_terms, interaction_terms]
parallel_wald_chisq <- as.numeric(
    t(interaction_estimates) %*%
        solve(interaction_cov, interaction_estimates)
)
parallel_test <- data.frame(
    working_correlation = "independence",
    hypothesis = "All TG-by-visit interaction coefficients are zero",
    df = length(interaction_terms),
    wald_chisq = parallel_wald_chisq,
    p_value = pchisq(
        parallel_wald_chisq,
        df = length(interaction_terms),
        lower.tail = FALSE
    )
)

write.csv(
    parallel_test,
    file.path(output_dir, "gee_parallel_trajectory_test.csv"),
    row.names = FALSE,
    quote = TRUE
)

qic_names <- c("QIC", "QICu", "Quasi Lik", "CIC", "params")
qic_row <- function(fit, working_correlation) {
    qic_values <- geepack::QIC(fit)
    qic_table <- data.frame(
        working_correlation = working_correlation,
        t(as.numeric(qic_values[qic_names])),
        row.names = NULL,
        check.names = FALSE
    )
    names(qic_table) <- c("working_correlation", qic_names)
    qic_table
}

qic_comparison <- rbind(
    qic_row(fit_independence, "independence"),
    qic_row(fit_ar1, "ar1")
)

write.csv(
    qic_comparison,
    file.path(output_dir, "gee_qic_comparison.csv"),
    row.names = FALSE,
    quote = TRUE
)

independence_estimates <- coef(fit_independence)
ar1_estimates <- coef(fit_ar1)
parameters <- union(names(independence_estimates), names(ar1_estimates))

estimate_comparison <- data.frame(
    parameter = parameters,
    independence = as.numeric(independence_estimates[parameters]),
    ar1 = as.numeric(ar1_estimates[parameters]),
    row.names = NULL
)

write.csv(
    estimate_comparison,
    file.path(output_dir, "gee_parameter_estimates.csv"),
    row.names = FALSE,
    quote = TRUE
)
