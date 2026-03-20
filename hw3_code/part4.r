library(tidyverse)
library(nlme)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("./hw3_code")

fit <- readRDS("model_fit.rds")

vcov.model <- readRDS("beta_vcov.rds")

# Test of parallelism

beta.model <- coef(fit)
Sig.model  <- vcov.model$Sig.model
L          <- rbind(c(0, 0, 0, 0, 0, 0, 1, 0, -1),
                    c(0, 0, 0, 0, 0, 0, 0, 1, -1))
L

### we use Between-Within method for dof estimation
nsub       <- nlevels(getGroups(fit.reml.c.exp.unequal_var))
nobs       <- nobs(fit)
# number of between-subject effects
nbetween   <- 6
# number of within-subject effects
nwithin    <- 3

L.beta     <- L %*% beta.model
L.beta.cov <- L %*% Sig.model %*% t(L)
F.stat     <- t(L.beta) %*% solve(L.beta.cov) %*% L.beta
ddfm       <- nobs - nsub - nwithin
ndfm       <- nrow(L)
prob.F     <- pf(F.stat,df1 = ndfm,df2 = ddfm, lower.tail=FALSE)


parallel.reml <- data.frame(
            'F' = F.stat, 
            'P-value' = prob.F,
            'NDFM'    = ndfm,
            'DDFM'    = ddfm
            )

library(kableExtra)
round(parallel.reml , 4) %>% 
    mutate_if(is.numeric, format, digits=4) %>% 
    kbl(caption = "Test of Parallelism" ) %>% 
    kable_styling('hover', full_width = T)  



### Additional test

# bud vs plbo

L <- rbind(c(0, 0, 0, 0, 0, 0, 1, 0, -1))
L

L.beta     <- L %*% beta.model
L.beta.cov <- L %*% Sig.model %*% t(L)
F.stat     <- t(L.beta) %*% solve(L.beta.cov) %*% L.beta
ddfm       <- nobs - nsub - nwithin
ndfm       <- nrow(L)
prob.F     <- pf(F.stat,df1 = ndfm,df2 = ddfm, lower.tail=FALSE)


parallel.reml <- data.frame(
            'F' = F.stat, 
            'P-value' = prob.F,
            'NDFM'    = ndfm,
            'DDFM'    = ddfm
            )

library(kableExtra)
round(parallel.reml , 4) %>% 
    mutate_if(is.numeric, format, digits=4) %>% 
    kbl(caption = "bud vs plbo" ) %>% 
    kable_styling('hover', full_width = T)  

# ned vs plbo

L <- rbind(c(0, 0, 0, 0, 0, 0, 0, 1, -1))
L

L.beta     <- L %*% beta.model
L.beta.cov <- L %*% Sig.model %*% t(L)
F.stat     <- t(L.beta) %*% solve(L.beta.cov) %*% L.beta
ddfm       <- nobs - nsub - nwithin
ndfm       <- nrow(L)
prob.F     <- pf(F.stat,df1 = ndfm,df2 = ddfm, lower.tail=FALSE)


parallel.reml <- data.frame(
            'F' = F.stat, 
            'P-value' = prob.F,
            'NDFM'    = ndfm,
            'DDFM'    = ddfm
            )

library(kableExtra)
round(parallel.reml , 4) %>% 
    mutate_if(is.numeric, format, digits=4) %>% 
    kbl(caption = "ned vs plbo" ) %>% 
    kable_styling('hover', full_width = T)  
