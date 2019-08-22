source("data_prep_basic.R")

model.dump.string <- getwd()

### Do an OLS k-fold cross-validation

ols.k.sims <- list()

for (i in 1:k.number) {
  
  # Recentre X and Y
  
  X.means <- apply(X.trains[[i]], 2, mean)
  
  X.train <- sweep(X.trains[[i]], 2, X.means)
  
  X.tilde <- sweep(X.tests[[i]], 2, X.means)
  
  y.bar <- mean(y.trains[[i]])
  
  y.train <- y.trains[[i]] - y.bar
  
  y.test <- y.tests[[i]] - y.bar
  
  ols.k.mod <- lm(y.train ~ X.train - 1)

  ols.k.y.tilde <- X.tilde %*% ols.k.mod$coefficients

  ols.k.sims[[i]] <- list(y.tilde = ols.k.y.tilde,
                          y.test = y.test)
  
  }
  
saveRDS(ols.k.sims, paste0(model.dump.string, "ols_k_sims_basic.RDS"))

### Do an frequentist ridge k-fold cross-validation

# Basic

freq.ridge.k.sims <- list()

for (i in 1:k.number) {
  
  # Recentre X and Y
  
  X.means <- apply(X.trains[[i]], 2, mean)
  
  X.train <- sweep(X.trains[[i]], 2, X.means)
  
  X.tilde <- sweep(X.tests[[i]], 2, X.means)
  
  y.bar <- mean(y.trains[[i]])
  
  y.train <- y.trains[[i]] - y.bar
  
  y.test <- y.tests[[i]] - y.bar
  
  ### Run freq ridge
  
  lambda.seq <- exp(seq(-10, 10, 0.01))
  
  ridge.cv.lambda <- cv.glmnet(X.train, y.train, alpha = 0, intercept = FALSE, lambda = lambda.seq)
  
  ridge.cv.opt.lambda <- ridge.cv.lambda$lambda.min
  
  freq.ridge.model <- glmnet(X.train, y.train, alpha = 0, lambda = ridge.cv.opt.lambda, intercept = FALSE)

  freq.ridge.k.ytilde <- X.tilde %*% freq.ridge.model$beta
  
  freq.ridge.k.sims[[i]] <- list(y.tilde = freq.ridge.k.ytilde,
                                 y.test = y.test)
  
  }
  
saveRDS(freq.ridge.k.sims, paste0(model.dump.string, "freq_ridge_k_sims_basic.RDS"))

### Do an frequentist lasso k-fold cross-validation

# Basic

freq.lasso.k.sims <- list()

for (i in 1:k.number) {
  
  # Recentre X and Y
  
  X.means <- apply(X.trains[[i]], 2, mean)
  
  X.train <- sweep(X.trains[[i]], 2, X.means)
  
  X.tilde <- sweep(X.tests[[i]], 2, X.means)
  
  y.bar <- mean(y.trains[[i]])
  
  y.train <- y.trains[[i]] - y.bar
  
  y.test <- y.tests[[i]] - y.bar
  
  # Run frequentist LASSO
  
  lambda.seq <- exp(seq(-10, 10, 0.01))
  
  lasso.cv.lambda <- cv.glmnet(X.train, y.train, alpha = 1, intercept = FALSE, lambda = lambda.seq)
  
  lasso.cv.opt.lambda <- lasso.cv.lambda$lambda.min
  
  freq.lasso.model <- glmnet(X.train, y.train, alpha = 1, lambda = lasso.cv.opt.lambda, intercept = FALSE)
  
  freq.lasso.k.ytilde <- X.tilde %*% freq.lasso.model$beta
  
  freq.lasso.k.sims[[i]] <- list(y.tilde = freq.lasso.k.ytilde,
                                 y.test = y.test)
  
}

saveRDS(freq.lasso.k.sims, paste0(model.dump.string, "freq_lasso_k_sims_basic.RDS"))

### BLM

source("blm_stan_functions.R")

blm.k.sims <- list()

for (i in 1:k.number) {
  
  blm.k.sims[[i]] <- blm_stan_simulation(y.trains[[i]], X.trains[[i]], a = 0.001, b = 0.001, y.tests[[i]], X.tests[[i]])
  
}

saveRDS(blm.k.sims, paste0(model.dump.string, "blm_k_sims_basic.RDS"))

rm(blm.k.sims)

### Ridge

source("ridge_stan_functions.R")

ridge.k.sims <- list()

for (i in 1:k.number) {
  
  ridge.k.sims[[i]] <- ridge_stan_simulation("ridge.stan", y.trains[[i]], X.trains[[i]], 
                                                    y.tests[[i]], X.tests[[i]],
                                                    pars = list(a = 0.001, b = 0.001, r = 0.001, d = 0.001))
  
}

saveRDS(ridge.k.sims, paste0(model.dump.string, "ridge_k_sims_basic.RDS"))

rm(ridge.k.sims)

### Lasso

source("lasso_stan_functions.R")

lasso.k.sims <- list()

for (i in 1:k.number) {
  
  lasso.k.sims[[i]] <- lasso_stan_simulation("lasso.stan",
                                             y.trains[[i]], 
                                             X.trains[[i]],
                                             y.tests[[i]], X.tests[[i]],
                                             pars = list(a = 0.001, b = 0.001, r = 0.001, d = 0.001))
  
  
}

saveRDS(lasso.k.sims, paste0(model.dump.string, "lasso_k_sims_basic.RDS"))

rm(lasso.k.sims)

