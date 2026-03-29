library(tidyverse)
library(nlme)

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

#########################

# fmla <- PostBD_FEV ~ 0 + TG:visitc_f + age_rz + gender + ethnic
# This design matrix is singular. Don't use it!
# Use this one:
fmla <- PostBD_FEV ~ 0 + TG + TG:visitc_f + age_rz + gender + ethnic

### CS
fit.cs <- gls(
    fmla,
    data = dat2,
    method = "REML",
    correlation = corCompSymm(form = ~ 1 | id), # `~ visitc | id` makes no difference
    weights = varIdent(form = ~ 1),
    na.action = na.omit
)

### AR1
fit.ar1 <- gls(
    fmla,
    data = dat2,
    method = "REML",
    correlation = corAR1(form = ~ visitc | id), # should not use visitc_f (factor)!
    weights = varIdent(form = ~ 1),
    na.action = na.omit
) 

### Exp
fit.exp <- gls(
    fmla,
    data = dat2,
    method = "REML",
    correlation = corExp(form = ~ visitc | id), # corCAR1 also works here
    weights = varIdent(form = ~ 1),
    na.action = na.omit
)

library(mmrm)
mmrm_opts <- mmrm:::h_get_optimizers()
mmrm_ctrl <- mmrm_control(optimizers = list(BFGS = mmrm_opts[["BFGS"]]), method = "Satterthwaite")

### Toeplitz
fit.toep <- mmrm(
    PostBD_FEV ~ 0 + TG + TG:visitc_f + age_rz + gender + ethnic + toep(visitc_f | id), # note that visitc_f must be used here (factor)
    data = dat2,
    control = mmrm_ctrl,
    reml = TRUE
)

### ANTE
fit.ante <- mmrm(
    PostBD_FEV ~ 0 + TG + TG:visitc_f + age_rz + gender + ethnic + ad(visitc_f | id),
    data = dat2,
    control = mmrm_ctrl,
    reml = TRUE
)

aic <- function(loglik, n_par) {
    ll <- as.numeric(loglik)
    -2 * ll + 2 * n_par
}

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

model_list <- list(
    fit.cs = fit.cs,
    fit.ar1 = fit.ar1,
    fit.exp = fit.exp,
    fit.toep = fit.toep,
    fit.ante = fit.ante
)

dof_values <- vapply(model_list, dof_gls_style, numeric(1))
loglik_values <- vapply(model_list, function(fit) as.numeric(logLik(fit)), numeric(1))
aic_values <- vapply(
    names(model_list),
    function(nm) aic(logLik(model_list[[nm]]), dof_values[[nm]]),
    numeric(1)
)

names(model_list)
model_names <- c("Compound Symmetry", "AR(1)", "Exponential", "Toeplitz", "ANTE")
model_compare <- data.frame(
    model = model_names,
    dof = as.numeric(dof_values),
    logLik = as.numeric(loglik_values),
    AIC = as.numeric(aic_values),
    row.names = NULL
)
model_compare

write.csv(
    model_compare,
    "./output/model_compare.csv",
    row.names = FALSE,
    quote = FALSE
)

############
# Heterogeneous variance
mmrm_ctrl_group <- mmrm_control(
    optimizers = list(BFGS = mmrm_opts[["BFGS"]], nlminb = mmrm_opts[["nlminb"]]),
    method = "Satterthwaite"
)

### ANTE (heterogeneous)
fit.adh <- mmrm(
    PostBD_FEV ~ 0 + TG + TG:visitc_f + age_rz + gender + ethnic + adh(visitc_f | id),
    data = dat2,
    control = mmrm_ctrl_group,
    reml = TRUE
)

### ANTE (homogeneous), group-specific covariance by TG
fit.ad_group <- mmrm(
    PostBD_FEV ~ 0 + TG + TG:visitc_f + age_rz + gender + ethnic + ad(visitc_f | TG / id),
    data = dat2,
    control = mmrm_ctrl_group,
    reml = TRUE
)

### ANTE (heterogeneous), group-specific covariance by TG
fit.adh_group <- mmrm(
    PostBD_FEV ~ 0 + TG + TG:visitc_f + age_rz + gender + ethnic + adh(visitc_f | TG / id),
    data = dat2,
    control = mmrm_ctrl_group,
    reml = TRUE
)


model_list_mmrm_var <- list(
    fit.ante = fit.ante,
    fit.adh = fit.adh,
    fit.ad_group = fit.ad_group,
    fit.adh_group = fit.adh_group
)

dof_values_mmrm_var <- vapply(model_list_mmrm_var, dof_gls_style, numeric(1))
loglik_values_mmrm_var <- vapply(model_list_mmrm_var, function(fit) as.numeric(logLik(fit)), numeric(1))
aic_values_mmrm_var <- vapply(
    names(model_list_mmrm_var),
    function(nm) aic(logLik(model_list_mmrm_var[[nm]]), dof_values_mmrm_var[[nm]]),
    numeric(1)
)

names(model_list_mmrm_var)
model_names_mmrm_var <- c("Homogeneous variance", 
"Heterogeneous variance", 
"Homogeneous variance with group-specific covariance",
"Heterogeneous variance with group-specific covariance")
aic_compare_mmrm_var <- data.frame(
    model = model_names_mmrm_var,
    dof = as.numeric(dof_values_mmrm_var),
    logLik = as.numeric(loglik_values_mmrm_var),
    AIC = as.numeric(aic_values_mmrm_var),
    row.names = NULL
)
aic_compare_mmrm_var

write.csv(
    aic_compare_mmrm_var,
    "./output/aic_compare_mmrm_var.csv",
    row.names = FALSE,
    quote = FALSE
)
