source("data_prep_basic.R")

model.dump.string <- getwd()

### Define the lambda prior parameters

lambda.seq <- c(0.0001, 0.001, 0.01, 1, 2)

lambda.array <- data.table(expand.grid(lambda.seq, lambda.seq))
setnames(lambda.array, c("r", "d"))

lambda.array <- lambda.array[d < 1,]

### Ridge model

source("ridge_stan_functions.R")

for (i in 1:NROW(lambda.array)) {

  # Full X matrix
  
  ridge.lambda.sims <- ridge_stan_simulation("ridge.stan",y , X,
                                             pars = list(a = 0.001, b = 0.001,
                                                         r = lambda.array$r[i], 
                                                         d = lambda.array$d[i]))
  
  saveRDS(ridge.lambda.sims, paste0(model.dump.string, "ridge_lambda_sims_basic_",i,".RDS"))
  
  ridge.lambda.k.sims <- list()
   
  for (j in 1:l.k.number) {

    # Full X matrix

    ridge.lambda.k.sims[[j]] <- ridge_stan_simulation("ridge.stan",
                                                      ly.trains[[j]],
                                                      lX.trains[[j]],
                                                      ly.tests[[j]],
                                                      lX.tests[[j]],
                                                      pars = list(a = 0.001, b = 0.001,
                                                                  r = lambda.array$r[i],
                                                                  d = lambda.array$d[i]))
    }

    saveRDS(ridge.lambda.k.sims, paste0(model.dump.string, "ridge_lambda_k_sims_basic_",i,".RDS"))

  }

### Lasso model

source("lasso_stan_functions.R")

for (i in 1:NROW(lambda.array)) {
  
  # Full X matrix
  
  lasso.lambda.sims <- lasso_stan_simulation("lasso.stan",y , X,
                                             pars = list(a = 0.001, b = 0.001,
                                                         r = lambda.array$r[i], 
                                                         d = lambda.array$d[i]))
  
  saveRDS(lasso.lambda.sims, paste0(model.dump.string, "lasso_lambda_sims_basic_",i,".RDS"))
  
  lasso.lambda.k.sims <- list()
  
  for (j in 1:l.k.number) {

    # Full X matrix

    lasso.lambda.k.sims[[j]] <- lasso_stan_simulation("lasso.stan",
                                                      ly.trains[[j]],
                                                      lX.trains[[j]],
                                                      ly.tests[[j]],
                                                      lX.tests[[j]],
                                                      pars = list(a = 0.001, b = 0.001,
                                                                  r = lambda.array$r[i],
                                                                  d = lambda.array$d[i]))

  }

  saveRDS(lasso.lambda.k.sims, paste0(model.dump.string, "lasso_lambda_k_sims_basic_",i,".RDS"))

  }