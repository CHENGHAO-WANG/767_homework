
#####################
# Part 1
#####################

library(tidyverse)
library(nlme)
library(future)
library(future.apply)
library(progressr)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("./hw5_code")

output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

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

fit <- nlme::lme(
    PostBD_FEV ~ 1 + TG:visitc + age_rz + gender + ethnic,
    random = ~ visitc | id,
    correlation = NULL,
    weights = varIdent(form = ~ 1 | TG), # form = ~ 1 | visitc_f * TG gives a convergence error
    data = dat2,
    na.action = na.omit,
    method = "REML"
)

random_effect_cov <- unclass(getVarCov(fit)) %>%
    as.data.frame() %>%
    rownames_to_column(var = "random_effect")

write.csv(
    random_effect_cov,
    file.path(output_dir, "random_effect_covariance_matrix.csv"),
    row.names = FALSE,
    quote = TRUE
)

residual_sd_ratio <- coef(fit$modelStruct$varStruct, unconstrained = FALSE, allCoef = TRUE)
residual_variance_table <- tibble(
    TG = names(residual_sd_ratio),
    residual_variance = sigma(fit)^2 * residual_sd_ratio^2
)

write.csv(
    residual_variance_table,
    file.path(output_dir, "residual_variance_by_group.csv"),
    row.names = FALSE,
    quote = TRUE
)

effect_table <- summary(fit)$tTable %>%
    as.data.frame() %>%
    rownames_to_column(var = "effect_name") %>%
    mutate(
        effect_name = recode(
            effect_name,
            "(Intercept)" = "Intercept",
            "age_rz" = "Age at randomization",
            "gendermale" = "Male",
            "ethnicblack" = "Black",
            "ethnichispanic" = "Hispanic",
            "ethnicother" = "Other ethnicity",
            "TGbudesonide:visitc" = "Budesonide x visitc",
            "TGnedocromil:visitc" = "Nedocromil x visitc",
            "TGplacebo:visitc" = "Placebo x visitc",
            .default = effect_name
        )
    ) %>%
    transmute(
        effect_name = effect_name,
        estimate = Value, #round(Value, 3),
        se = Std.Error, #round(Std.Error, 3),
        test_statistic = round(`t-value`, 3),
        p_value = case_when(
            is.na(`p-value`) ~ NA_character_,
            `p-value` < 0.001 ~ "<0.001",
            TRUE ~ formatC(`p-value`, format = "f", digits = 3)
        )
    )

write.csv(
    effect_table,
    file.path(output_dir, "effect_table.csv"),
    row.names = FALSE,
    quote = TRUE
)



#####################
# Part 2
#####################

### Evaluate the need for random effects

fit.no_rslope <- nlme::lme(
    PostBD_FEV ~ 1 + TG:visitc + age_rz + gender + ethnic,
    random = ~ 1 | id,
    correlation = NULL,
    weights = varIdent(form = ~ 1 | TG),
    data = dat2,
    na.action = na.omit,
    method = "REML"
)
lrt.rslope <- as.numeric(2*(logLik(fit) - logLik(fit.no_rslope)))
lrt.rslope
# 3691.1

pval.rslope <- .5 * pchisq(lrt.rslope, df = 1, lower.tail = FALSE) + .5 * pchisq(lrt.rslope, df = 2, lower.tail = FALSE)
pval.rslope
# < 0.001

fit.no_rint <- nlme::lme(
    PostBD_FEV ~ 1 + TG:visitc + age_rz + gender + ethnic,
    random = ~ 0 + visitc | id,
    correlation = NULL,
    weights = varIdent(form = ~ 1 | TG),
    data = dat2,
    na.action = na.omit,
    method = "REML"
)

lrt.rint <- as.numeric(2*(logLik(fit) - logLik(fit.no_rint)))
lrt.rint
10092.1

pval.rint <- .5 * pchisq(lrt.rint, df = 1, lower.tail = FALSE) + .5 * pchisq(lrt.rint, df = 2, lower.tail = FALSE)
pval.rint
# < 0.001

fit.no_re <- nlme::gls(
    PostBD_FEV ~ 1 + TG:visitc + age_rz + gender + ethnic,
    correlation = NULL,
    weights = varIdent(form = ~ 1 | TG),
    data = dat2,
    na.action = na.omit,
    method = "REML"
)

aic_comparison <- AIC(fit, fit.no_rslope, fit.no_rint, fit.no_re) %>%
    as.data.frame() %>%
    mutate(
        model = c(
            "Mixed model",
            "No random slope",
            "No random intercept",
            "No random effects"
        )
    ) %>%
    select(model, everything())

write.csv(
    aic_comparison,
    file.path(output_dir, "aic_comparison.csv"),
    row.names = FALSE,
    quote = TRUE
)


#####################
# Part 3
#####################

eblup <- ranef(fit) %>%
    rownames_to_column("id")


hist_re_intercept <- ggplot(eblup, aes(x = `(Intercept)`)) +
    geom_histogram(bins = 30, color = "white", fill = "steelblue") +
    labs(
        title = "Histogram of Random Intercepts",
        x = "Random intercept eBLUP",
        y = "Count"
    ) +
    theme_minimal()

hist_re_slope <- ggplot(eblup, aes(x = visitc)) +
    geom_histogram(bins = 30, color = "white", fill = "darkorange") +
    labs(
        title = "Histogram of Random Slopes",
        x = "Random slope eBLUP",
        y = "Count"
    ) +
    theme_minimal()

ggsave(file.path(output_dir, "random_intercept_histogram.png"), hist_re_intercept, width = 7, height = 5, dpi = 300)
ggsave(file.path(output_dir, "random_slope_histogram.png"), hist_re_slope, width = 7, height = 5, dpi = 300)

set.seed(767)
n_subjects_per_tg <- 4

sampled_subjects <- dat2 %>%
    distinct(id, TG) %>%
    filter(!is.na(TG)) %>%
    group_by(TG) %>%
    reframe(id = sample(as.character(id), size = min(n_subjects_per_tg, dplyr::n()))) %>%
    ungroup()

subject_traj_plot_data <- dat2 %>%
    mutate(
        id = as.character(id),
        # level = 1 makes predict() use eBLUPs. With level = 0, only fixed effect estimates are used.
        # if id is from the fitted data, then that id's eBLUP can be used;
        # if id is a brand-new subject not seen in fitting, then eBLUPs cannot be used.
        subject_fitted = as.numeric(predict(fit, newdata = dat2, level = 1))
    ) %>%
    semi_join(sampled_subjects, by = c("id", "TG")) %>%
    arrange(TG, id, visitc)

fitted_traj_plot <- ggplot(subject_traj_plot_data, aes(x = visitc, y = PostBD_FEV, group = id, color = id)) +
    geom_line(
        data = ~ dplyr::filter(.x, !is.na(PostBD_FEV)),
        aes(linetype = "Observed"),
        linewidth = 0.5,
        alpha = 0.45
    ) +
    geom_point(
        data = ~ dplyr::filter(.x, !is.na(PostBD_FEV)),
        size = 1.5,
        alpha = 0.7
    ) +
    geom_line(
        aes(y = subject_fitted, linetype = "Fitted"),
        linewidth = 0.9
    ) +
    scale_linetype_manual(
        name = "Trajectory",
        values = c(Observed = "dashed", Fitted = "solid")
    ) +
    facet_wrap(~ TG) +
    labs(
        # title = "Observed and Subject-Specific Fitted Trajectories",
        # subtitle = paste("Random sample of", n_subjects_per_tg, "subjects per treatment group"),
        x = "Visit Time (months)",
        y = "PostBD FEV1",
        color = "Subject ID"
    ) +
    theme_minimal()

ggsave(
    file.path(output_dir, "subject_specific_trajectory_plot.png"),
    fitted_traj_plot,
    width = 10,
    height = 6,
    dpi = 300
)



#####################
# Part 4
#####################

### residual diagnostics
png(file.path(output_dir, "fit_diagnostic_plot.png"), width = 1200, height = 900, res = 150)
plot(fit)
dev.off()

pearson_residual_df <- tibble(
    observation = seq_along(resid(fit, type = "pearson")),
    id = names(resid(fit, type = "pearson")),
    pearson_residual = as.numeric(resid(fit, type = "pearson"))
) %>%
    arrange(desc(abs(pearson_residual)))

head(pearson_residual_df)
# id 1038 has the extremest residual

# qqnorm(resid(fit, type = "normalized"))

# plot(fit, resid(., type = "normalized") ~ fitted(.))

### Influence diagnostics

ids <- levels(droplevels(dat2$id))
beta_full <- fixef(fit)
vcov_beta_full <- fit$varFix

mf_fixed <- model.frame(formula(fit), data = getData(fit), na.action = na.omit)
rank_x <- qr(model.matrix(formula(fit), data = mf_fixed))$rank

loo_control <- nlme::lmeControl(
    maxIter = 10,
    msMaxIter = 10,
    msMaxEval = 50,
    returnObject = TRUE
)

workers <- max(1L, future::availableCores() - 1L)
old_plan <- future::plan()
on.exit(future::plan(old_plan), add = TRUE)
future::plan(future::multisession, workers = workers)

progressr::handlers(global = TRUE)
cook_results <- progressr::with_progress({
    p <- progressr::progressor(along = ids)

    future.apply::future_lapply(
        ids,
        future.packages = c("nlme", "dplyr"),
        future.seed = TRUE,
        FUN = function(subject_id) {
            p(sprintf("Refit without subject %s", subject_id))

            dat_minus_i <- dplyr::filter(dat2, id != subject_id) %>%
                droplevels()

            refit <- tryCatch(
                update(fit, data = dat_minus_i, control = loo_control),
                error = function(e) e
            )

            if (inherits(refit, "error")) {
                return(data.frame(
                    id = as.character(subject_id),
                    cooks_distance = NA_real_,
                    converged = FALSE,
                    message = conditionMessage(refit)
                ))
            }

            beta_minus_i <- fixef(refit)
            delta_beta <- beta_minus_i - beta_full
            cooks_distance <- as.numeric(t(delta_beta) %*% solve(vcov_beta_full, delta_beta) / rank_x)

            data.frame(
                id = as.character(subject_id),
                cooks_distance = cooks_distance,
                converged = TRUE,
                message = NA_character_
            )
        }
    )
})

cook_df <- bind_rows(cook_results) %>%
    mutate(id = factor(id, levels = id[order(cooks_distance, na.last = TRUE)]))

write.csv(cook_df, file.path(output_dir, "subject_level_cooks_distance.csv"), row.names = FALSE)

all(cook_df$converged)
# TRUE

cook_plot <- cook_df %>%
    arrange(desc(cooks_distance)) %>%
    mutate(plot_index = seq_len(n())) %>%
    ggplot(aes(x = plot_index, y = cooks_distance)) +
    geom_point(size = 2, color = "steelblue") +
    scale_x_continuous(
        breaks = NULL,
        expand = expansion(mult = c(0.03, 0.03))
    ) +
    labs(
        title = "Subject-Level Cook's Distance",
        x = "Subject ID",
        y = "Cook's distance"
    ) +
    theme_minimal() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )

ggsave(
    file.path(output_dir, "subject_level_cooks_distance.png"),
    cook_plot,
    width = 8,
    height = 5,
    dpi = 300
)

cook_df %>% arrange(desc(cooks_distance)) %>% head(10)

# top 2 cook's distances: id = 803, 345


#####################
# Part 5
#####################

slope_coef_names <- c(
    "TGbudesonide:visitc",
    "TGnedocromil:visitc",
    "TGplacebo:visitc"
)

beta_hat <- fixef(fit)
vcov_beta <- fit$varFix

L <- matrix(
    0,
    nrow = 2,
    ncol = length(beta_hat),
    dimnames = list(
        c(
            "budesonide_vs_nedocromil",
            "budesonide_vs_placebo"
        ),
        names(beta_hat)
    )
)
L["budesonide_vs_nedocromil", "TGbudesonide:visitc"] <- 1
L["budesonide_vs_nedocromil", "TGnedocromil:visitc"] <- -1
L["budesonide_vs_placebo", "TGbudesonide:visitc"] <- 1
L["budesonide_vs_placebo", "TGplacebo:visitc"] <- -1

wald_diff <- as.vector(L %*% beta_hat)
wald_cov <- L %*% vcov_beta %*% t(L)
wald_stat <- as.numeric(t(wald_diff) %*% solve(wald_cov, wald_diff))
wald_df <- nrow(L)
wald_p_value <- pchisq(wald_stat, df = wald_df, lower.tail = FALSE)

part5_wald_test <- tibble(
    hypothesis = "TGbudesonide:visitc = TGnedocromil:visitc = TGplacebo:visitc",
    wald_chisq = wald_stat,
    df = wald_df,
    p_value = wald_p_value
)

write.csv(
    part5_wald_test,
    file.path(output_dir, "part5_equal_slope_wald_test.csv"),
    row.names = FALSE,
    quote = TRUE
)

part5_wald_test

