################################
# Part 2: Univariate Repeated Measures ANOVA
################################

library(car)
library(dplyr)
library(tidyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

dat <- read.csv("../data/complete_case.csv")

# Only consider two variables: TG, visitc
# response variable: PostBD_FEV

# long to wide
dat.wide <- dat %>%
  select(id, TG, visitc, PostBD_FEV) %>%
  arrange(id, visitc) %>%
  pivot_wider(
    id_cols = c(id, TG),
    names_from = visitc,
    values_from = PostBD_FEV,
    names_prefix = "visitc_"
  )

# Set coding scheme
dat.wide$TG <- factor(dat.wide$TG)
contrasts(dat.wide$TG) <- contr.sum(nlevels(dat.wide$TG))

response <- as.matrix(dat.wide[, -(1:2)])
fit <- lm(response ~ TG, data = dat.wide)

idata <- data.frame(visitc = factor(sort(unique(dat$visitc)), levels = sort(unique(dat$visitc))))
fit_anova <- car::Anova(fit, idata = idata, idesign = ~ visitc, type = 3)

a <- summary(fit_anova, multivariate = FALSE, univariate = TRUE)
a
writeLines(capture.output(a), "./output/rm_anova.txt")


################################
# Part 3: MANOVA
################################

b <- summary(fit_anova, multivariate = TRUE, univariate = FALSE)
b
writeLines(capture.output(b), "./output/manova.txt")

### Test of group effect

fit_manova <- manova(response ~ TG, data = dat.wide)
res.manova_group_effect <- summary(fit_manova, test = "Wilks")
res.manova_group_effect
writeLines(capture.output(res.manova_group_effect), "./output/manova_group_effect.txt")

### Test of time effect

fit0 <- lm(response ~ TG - 1, data = dat.wide)

n_time <- ncol(response)
M.rows <- vector(mode = "list", length = n_time)
for (i in 2:n_time) {
  v <- rep(0, n_time)
  v[1] <- 1
  v[i] <- -1
  M.rows[[i]] <- v
}
M <- do.call(rbind, M.rows)

# colnames of M (rownames of P) must match the colnames of coef(fit0)
colnames(M) <- colnames(coef(fit0))
# rownames of M (colnames of P) are just labels for the contrasts; they don't affect the test results.
rownames(M) <- paste(colnames(response)[1], colnames(response)[2:n_time], sep = "-")

# car::linearHypothesis uses C %*% B %*% P = 0, so P is transpose of SAS M.
P <- t(M)

C_group <- diag(3)
# colnames of C must match the rownames of coef(fit0)
colnames(C_group) <- rownames(coef(fit0))
# rownames of C are just labels for the contrasts; they don't affect the test results.
rownames(C_group) <- rownames(coef(fit0))

test_time <- car::linearHypothesis(
  fit0,
  hypothesis.matrix = C_group,
  P = P
)

print(test_time)

writeLines(capture.output(test_time), "./output/manova_time_effect.txt")
