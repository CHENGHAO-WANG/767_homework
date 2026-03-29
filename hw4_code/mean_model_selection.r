library(tidyverse)
library(mmrm)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("./hw4_code")

if (!dir.exists("./output")) {
    dir.create("./output", recursive = TRUE)
}

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
    method = "Satterthwaite"
)

### saturated
fit.saturated <- mmrm(
    PostBD_FEV ~ 0 + TG + TG:visitc_f + age_rz + gender + ethnic + adh(visitc_f | TG / id),
    data = dat2,
    control = mmrm_ctrl_group,
    reml = FALSE
)

### linear trend
fit.linear <- mmrm(
    PostBD_FEV ~ 1 + TG:visitc + age_rz + gender + ethnic + adh(visitc_f | TG / id),
    data = dat2,
    control = mmrm_ctrl_group,
    reml = FALSE
)

### linear spline knot at 60
dat2$visitc.sp <- pmax(dat2$visitc - 60, 0) 

fit.lsp <- mmrm(
    PostBD_FEV ~ 1 + TG:visitc + TG:visitc.sp + age_rz + gender + ethnic + adh(visitc_f | TG / id),
    data = dat2,
    control = mmrm_ctrl_group,
    reml = FALSE
)

### quadratic with centered time
(visitc.mean <- mean(unique(dat2$visitc)))
dat2$visitc.centered <- dat2$visitc - visitc.mean 

fit.quad <- mmrm(
    PostBD_FEV ~ 1 + TG:visitc.centered + TG:I(visitc.centered^2) + age_rz + gender + ethnic + adh(visitc_f | TG / id),
    data = dat2,
    control = mmrm_ctrl_group,
    reml = FALSE
) 

##########################
# model comparison

aic <- function(loglik, n_par) {
    ll <- as.numeric(loglik)
    -2 * ll + 2 * n_par
}

bic_subject <- function(loglik, n_par, n_subject) {
    ll <- as.numeric(loglik)
    -2 * ll + log(n_subject) * n_par
}

dof_mmrm <- function(fit) {
    n_fix <- length(coef(fit, complete = FALSE))
    n_cov <- length(mmrm::component(fit, "theta_est"))
    n_fix + n_cov
}

mean_model_list <- list(
    fit.saturated = fit.saturated,
    fit.linear = fit.linear,
    fit.lsp = fit.lsp,
    fit.quad = fit.quad
)

dof_values_mean <- vapply(mean_model_list, dof_mmrm, numeric(1))
loglik_values_mean <- vapply(mean_model_list, function(fit) as.numeric(logLik(fit)), numeric(1))
n_subject <- dplyr::n_distinct(dat2$id)
aic_values_mean <- vapply(
    names(mean_model_list),
    function(nm) aic(logLik(mean_model_list[[nm]]), dof_values_mean[[nm]]),
    numeric(1)
)

bic_values_mean <- vapply(
    names(mean_model_list),
    function(nm) bic_subject(logLik(mean_model_list[[nm]]), dof_values_mean[[nm]], n_subject),
    numeric(1)
)

mean_aic_compare <- data.frame(
    model = c("Saturated", "Linear trend", "Linear spline (knot = 60)", "Quadratic"),
    dof = as.numeric(dof_values_mean),
    logLik = as.numeric(loglik_values_mean),
    AIC = as.numeric(aic_values_mean),
    BIC = as.numeric(bic_values_mean),
    row.names = NULL
) %>%
    arrange(AIC)

mean_aic_compare

write.csv(
    mean_aic_compare,
    "./output/mean_aic_compare.csv",
    row.names = FALSE,
    quote = FALSE
)

##########################
# linear spline mean model with candidate covariance structures
# REML estimation, homogeneous variance

fmla_lsp <- PostBD_FEV ~ 1 + TG:visitc + TG:visitc.sp + age_rz + gender + ethnic

fit.lsp.cs <- nlme::gls(
    fmla_lsp,
    data = dat2,
    method = "REML",
    correlation = nlme::corCompSymm(form = ~ 1 | id),
    na.action = na.omit
)

fit.lsp.ar1 <- nlme::gls(
    fmla_lsp,
    data = dat2,
    method = "REML",
    correlation = nlme::corAR1(form = ~ visitc | id),
    na.action = na.omit
)

fit.lsp.exp <- nlme::gls(
    fmla_lsp,
    data = dat2,
    method = "REML",
    correlation = nlme::corExp(form = ~ visitc | id),
    na.action = na.omit
)

mmrm_ctrl_reml <- mmrm_control(
    optimizers = list(BFGS = mmrm_opts[["BFGS"]]),
    method = "Satterthwaite"
)

fit.lsp.toep <- mmrm(
    PostBD_FEV ~ 1 + TG:visitc + TG:visitc.sp + age_rz + gender + ethnic + toep(visitc_f | id),
    data = dat2,
    control = mmrm_ctrl_reml,
    reml = TRUE
)

fit.lsp.ante <- mmrm(
    PostBD_FEV ~ 1 + TG:visitc + TG:visitc.sp + age_rz + gender + ethnic + ad(visitc_f | id),
    data = dat2,
    control = mmrm_ctrl_reml,
    reml = TRUE
)

dof_gls_style <- function(fit) {
    if (inherits(fit, "gls")) {
        fixSig <- attr(fit[["modelStruct"]], "fixedSigma")
        fixSig <- !is.null(fixSig) && fixSig
        p <- fit$dims$p
        return(p + length(coef(fit[["modelStruct"]])) + as.integer(!fixSig))
    }

    if (inherits(fit, c("mmrm", "mmrm_fit", "mmrm_tmb"))) {
        n_fix <- length(coef(fit, complete = FALSE))
        n_cov <- length(mmrm::component(fit, "theta_est"))
        return(n_fix + n_cov)
    }

    stop("Unsupported model class for dof computation.")
}

lsp_cov_model_list <- list(
    fit.lsp.cs = fit.lsp.cs,
    fit.lsp.ar1 = fit.lsp.ar1,
    fit.lsp.exp = fit.lsp.exp,
    fit.lsp.toep = fit.lsp.toep,
    fit.lsp.ante = fit.lsp.ante
)

lsp_cov_dof <- vapply(lsp_cov_model_list, dof_gls_style, numeric(1))
lsp_cov_loglik <- vapply(lsp_cov_model_list, function(fit) as.numeric(logLik(fit)), numeric(1))
lsp_cov_aic <- vapply(
    names(lsp_cov_model_list),
    function(nm) aic(logLik(lsp_cov_model_list[[nm]]), lsp_cov_dof[[nm]]),
    numeric(1)
)

lsp_cov_aic_compare <- data.frame(
    model = c("Compound Symmetry", "AR(1)", "Exponential", "Toeplitz", "ANTE"),
    dof = as.numeric(lsp_cov_dof),
    logLik = as.numeric(lsp_cov_loglik),
    AIC = as.numeric(lsp_cov_aic),
    row.names = NULL
) %>%
    arrange(AIC)

lsp_cov_aic_compare

write.csv(
    lsp_cov_aic_compare,
    "./output/lsp_cov_aic_compare.csv",
    row.names = FALSE,
    quote = FALSE
)


###########################
# verify variance structure

fit.lsp.var.ad <- mmrm(
    PostBD_FEV ~ 1 + TG:visitc + TG:visitc.sp + age_rz + gender + ethnic + ad(visitc_f | id),
    data = dat2,
    control = mmrm_ctrl_group,
    reml = TRUE
)

fit.lsp.var.adh <- mmrm(
    PostBD_FEV ~ 1 + TG:visitc + TG:visitc.sp + age_rz + gender + ethnic + adh(visitc_f | id),
    data = dat2,
    control = mmrm_ctrl_group,
    reml = TRUE
)

fit.lsp.var.ad_group <- mmrm(
    PostBD_FEV ~ 1 + TG:visitc + TG:visitc.sp + age_rz + gender + ethnic + ad(visitc_f | TG / id),
    data = dat2,
    control = mmrm_ctrl_group,
    reml = TRUE
)

fit.lsp.var.adh_group <- mmrm(
    PostBD_FEV ~ 1 + TG:visitc + TG:visitc.sp + age_rz + gender + ethnic + adh(visitc_f | TG / id),
    data = dat2,
    control = mmrm_ctrl_group,
    reml = TRUE
)

lsp_var_model_list <- list(
    fit.lsp.var.ad = fit.lsp.var.ad,
    fit.lsp.var.adh = fit.lsp.var.adh,
    fit.lsp.var.ad_group = fit.lsp.var.ad_group,
    fit.lsp.var.adh_group = fit.lsp.var.adh_group
)

lsp_var_dof <- vapply(lsp_var_model_list, dof_gls_style, numeric(1))
lsp_var_loglik <- vapply(lsp_var_model_list, function(fit) as.numeric(logLik(fit)), numeric(1))
lsp_var_aic <- vapply(
    names(lsp_var_model_list),
    function(nm) aic(logLik(lsp_var_model_list[[nm]]), lsp_var_dof[[nm]]),
    numeric(1)
)

lsp_var_aic_compare <- data.frame(
    model = c(
        "Homogeneous variance",
        "Heterogeneous variance",
        "Homogeneous variance with group-specific covariance",
        "Heterogeneous variance with group-specific covariance"
    ),
    dof = as.numeric(lsp_var_dof),
    logLik = as.numeric(lsp_var_loglik),
    AIC = as.numeric(lsp_var_aic),
    row.names = NULL
) %>%
    arrange(AIC)

lsp_var_aic_compare

write.csv(
    lsp_var_aic_compare,
    "./output/lsp_var_aic_compare.csv",
    row.names = FALSE,
    quote = FALSE
)

