library(tidyverse)
library(mmrm)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("./hw4_code")

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

mmrm_opts <- mmrm:::h_get_optimizers()
mmrm_ctrl_group <- mmrm_control(
    optimizers = list(BFGS = mmrm_opts[["BFGS"]], nlminb = mmrm_opts[["nlminb"]]),
    method = "Between-Within"
)

# fit <- mmrm(
#    PostBD_FEV ~ 1 + TG:visitc + age_rz + gender + ethnic + adh(visitc_f | TG / id),
#    data = dat2,
#    control = mmrm_ctrl_group,
#    reml = TRUE
#)

dat2$visitc.sp <- pmax(dat2$visitc - 60, 0) 

fit <- mmrm(
    PostBD_FEV ~ 1 + TG:visitc + TG:visitc.sp + age_rz + gender + ethnic + adh(visitc_f | TG / id),
    data = dat2,
    control = mmrm_ctrl_group,
    reml = FALSE
)


# Export estimated covariance matrices by treatment group from fit
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

cov_mats_by_tg <- mmrm::component(fit, "varcor")

for (tg_name in names(cov_mats_by_tg)) {
    cov_mat_tg <- as.matrix(cov_mats_by_tg[[tg_name]])
    write_matrix_outputs(
        cov_mat_tg,
        paste0("./output/fit_covariance_matrix_", tg_name, ".csv"),
        paste0("./output/fit_covariance_matrix_", tg_name, ".tex")
    )
}

# Fixed-effects table from fit
effect_table <- coef(summary(fit)) %>%
    as.data.frame() %>%
    rownames_to_column(var = "effect_name") %>%
    mutate(
        t_critical = qt(0.975, df = df),
        ci_lower = Estimate - t_critical * `Std. Error`,
        ci_upper = Estimate + t_critical * `Std. Error`,
        ci_95 = paste0(
            "[",
            formatC(ci_lower, format = "f", digits = 2),
            ", ",
            formatC(ci_upper, format = "f", digits = 2),
            "]"
        )
    ) %>%
    transmute(
        effect_name = effect_name,
        estimate = Estimate,
        se = `Std. Error`,
        t_statistic = `t value`,
        dof = df,
        p_value = case_when(
            is.na(`Pr(>|t|)`) ~ NA_character_,
            `Pr(>|t|)` < 0.001 ~ "<0.001",
            TRUE ~ formatC(`Pr(>|t|)`, format = "f", digits = 3)
        ),
        ci_95 = ci_95
    )

write.csv(
    effect_table,
    "./output/effect_table.csv",
    row.names = FALSE,
    quote = TRUE
)

# Joint hypothesis test:
# H0:
# TGbudesonide:visitc = TGnedocromil:visitc = TGplacebo:visitc
# and
# TGbudesonide:visitc.sp = TGnedocromil:visitc.sp = TGplacebo:visitc.sp
beta_names <- names(mmrm::component(fit, "beta_est"))

required_terms <- c(
    "TGbudesonide:visitc",
    "TGnedocromil:visitc",
    "TGplacebo:visitc",
    "TGbudesonide:visitc.sp",
    "TGnedocromil:visitc.sp",
    "TGplacebo:visitc.sp"
)

contrast_joint <- matrix(
    0,
    nrow = 4,
    ncol = length(beta_names),
    dimnames = list(
        c(
            "visitc_budesonide_minus_nedocromil",
            "visitc_budesonide_minus_placebo",
            "visitc_sp_budesonide_minus_nedocromil",
            "visitc_sp_budesonide_minus_placebo"
        ),
        beta_names
    )
)

contrast_joint["visitc_budesonide_minus_nedocromil", "TGbudesonide:visitc"] <- 1
contrast_joint["visitc_budesonide_minus_nedocromil", "TGnedocromil:visitc"] <- -1

contrast_joint["visitc_budesonide_minus_placebo", "TGbudesonide:visitc"] <- 1
contrast_joint["visitc_budesonide_minus_placebo", "TGplacebo:visitc"] <- -1

contrast_joint["visitc_sp_budesonide_minus_nedocromil", "TGbudesonide:visitc.sp"] <- 1
contrast_joint["visitc_sp_budesonide_minus_nedocromil", "TGnedocromil:visitc.sp"] <- -1

contrast_joint["visitc_sp_budesonide_minus_placebo", "TGbudesonide:visitc.sp"] <- 1
contrast_joint["visitc_sp_budesonide_minus_placebo", "TGplacebo:visitc.sp"] <- -1

parallelism_test <- mmrm::df_md(fit, contrast_joint)

parallelism_test_table <- tibble(
    hypothesis = "Equal TG slopes for visitc and visitc.sp",
    num_df = parallelism_test$num_df,
    denom_df = parallelism_test$denom_df,
    f_stat = parallelism_test$f_stat,
    p_value = parallelism_test$p_val
)

parallelism_test_table

write.csv(
    parallelism_test_table,
    "./output/parallelism_test.csv",
    row.names = FALSE,
    quote = TRUE
)


# Fitted mean trajectories vs observed mean trajectories (linear spline model)
# Build covariate-standardized fitted means by TG/visit using baseline
# covariate distributions within each treatment group.
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
    mutate(prop_gender = if (sum(n_gender) > 0) n_gender / sum(n_gender) else rep(0, n())) %>%
    ungroup()

baseline_ethnic_prop_by_tg <- baseline_dat %>%
    filter(!is.na(ethnic)) %>%
    count(TG, ethnic, name = "n_ethnic") %>%
    group_by(TG) %>%
    complete(ethnic = levels(dat2$ethnic), fill = list(n_ethnic = 0)) %>%
    mutate(prop_ethnic = if (sum(n_ethnic) > 0) n_ethnic / sum(n_ethnic) else rep(0, n())) %>%
    ungroup()

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
        visitc.sp = pmax(visitc - 60, 0),
        visitc_f = factor(visitc, levels = levels(dat2$visitc_f)),
        id = factor(levels(dat2$id)[1], levels = levels(dat2$id)),
        weight = coalesce(prop_gender, 0) * coalesce(prop_ethnic, 0)
    ) %>%
    filter(!is.na(age_rz))

pred_grid$fitted <- as.numeric(predict(fit, newdata = pred_grid))

fitted_mean_by_tg_visit <- pred_grid %>%
    group_by(TG, visitc) %>%
    summarise(
        fitted_mean = sum(fitted * weight, na.rm = TRUE),
        weight_sum = sum(weight, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    arrange(TG, visitc)

if (any(abs(fitted_mean_by_tg_visit$weight_sum - 1) > 1e-8)) {
    warning("Weighted prediction rows do not sum to 1 for some TG/visitc combinations.")
}

observed_mean_by_tg_visit <- dat2 %>%
    group_by(TG, visitc) %>%
    summarise(
        observed_mean = mean(PostBD_FEV, na.rm = TRUE),
        n_obs = sum(!is.na(PostBD_FEV)),
        .groups = "drop"
    ) %>%
    arrange(TG, visitc)

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
    geom_vline(xintercept = 60, linetype = "dashed", color = "grey50", linewidth = 0.5) +
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
