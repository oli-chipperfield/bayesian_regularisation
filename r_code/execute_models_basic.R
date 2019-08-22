source("data_prep_basic.R")

model.dump.string <- getwd()

### Execute the standard versions of the models

### OLS for reference

## Basic

ols.model <- lm(y ~ X-1)

saveRDS(ols.model, paste0(model.dump.string, "ols_model_basic.RDS"))

### Frequentist ridge

## Basic

# Run cross-validation to get lambda

lambda.seq <- exp(seq(-10, 10, 0.01))

freq.ridge.cv <- cv.glmnet(X, y, alpha = 0, intercept = FALSE, lambda = lambda.seq)

#plot(freq.ridge.cv)

opt.ridge.lambda <- freq.ridge.cv$lambda.min

freq.ridge.model <- glmnet(X, y, alpha = 0, lambda = opt.ridge.lambda, intercept = FALSE)

saveRDS(freq.ridge.model, paste0(model.dump.string, "freq_ridge_model_basic.RDS"))

### Frequentist LASSO

## Basic

# Run cross-validation to get lambda

lambda.seq <- exp(seq(-10, 10, 0.01))

freq.lasso.cv <- cv.glmnet(X, y, alpha = 1, intercept = FALSE, lambda = lambda.seq)

#plot(freq.lasso.cv.full)

opt.lasso.lambda <- freq.lasso.cv$lambda.min

freq.lasso.model <- glmnet(X, y, alpha = 1, lambda = opt.lasso.lambda)

saveRDS(freq.lasso.model, paste0(model.dump.string, "freq_lasso_model_basic.RDS"))

### BLM model

source("blm_stan_functions.R")

blm.sims <- blm_stan_simulation(y, X, a = 0.001, b = 0.001)

saveRDS(blm.sims, paste0(model.dump.string, "blm_sims_basic.RDS"))

rm(blm.sims)

### Ridge model

source("ridge_stan_functions.R")

ridge.sims <- ridge_stan_simulation("ridge.stan", y, X, pars = list(a = 0.001, b = 0.001, r = 0.001, d = 0.001))

saveRDS(ridge.sims, paste0(model.dump.string, "ridge_sims_basic.RDS"))

rm(ridge.sims)

### Lasso model

source("lasso_stan_functions.R")

lasso.sims <- lasso_stan_simulation("lasso.stan", y, X, pars = list(a = 0.001, b = 0.001, r = 0.001, d = 0.001))

saveRDS(lasso.sims, paste0(model.dump.string, "lasso_sims_basic.RDS"))

rm(lasso.sims)

