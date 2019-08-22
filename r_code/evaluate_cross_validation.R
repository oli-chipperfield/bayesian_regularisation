source("data_prep_basic.R")
source("evaluate_model_parameters.R")

model.dump.string <- getwd()

#### Assess k-fold cross validation  ####

#### Calculate the NULL RMSPE

get_null_rmspe<- function(y.tests) {
  
  y.err <- rep(as.numeric(NA), length(y.tests))
  
  y.n <- rep(as.numeric(NA), length(y.tests))
  
  for (i in 1:length(y.tests)) {
    
    ### Calculate error
    
    yi.error <- (y.tests[[i]] - mean(y.tests[[i]])) ^2
    
    yi.sum.error <- sum(yi.error)
    
    y.err[i] <- yi.sum.error
    
    yi.n <- length(y.tests[[i]])
    
    y.n[i] <- yi.n
    
  }
  
  null.rmse <- sqrt(sum(y.err) / sum(y.n))
  
  return(list(mean = null.rmse, y.err = y.err, y.n = y.n))
  
}

null.rmspe <- get_null_rmspe(y.tests)

#### Calculate the frequentist RMSPE

## OLS

get_ols_rmspe <- function(ols.k.sims) {
  
  k.rmse <- list()
  k.n <- list()
  
  for (i in 1:length(ols.k.sims)) {
    
    sse <- sum((ols.k.sims[[i]]$y.tilde - ols.k.sims[[i]]$y.test)^2)
    mse <- sse / length(ols.k.sims[[i]]$y.tilde)
    rmse.i <- sqrt(mse)
    
    k.rmse[[i]] <- rmse.i
    k.n[[i]] <- length(ols.k.sims[[i]]$y.tilde)
    
  }
  
  rmse <- weighted.mean(unlist(k.rmse), unlist(k.n))
  
  return(list(rmse = rmse, k.rmse = k.rmse, k.n = k.n))
  
}

## Freq elastic nets

#freq.k.sims <- freq.ridge.k.sims 

get_freq_ea_rmspe <- function(freq.k.sims) {
  
  k.rmse <- list()
  k.n <- list()
  
  for (i in 1:length(freq.k.sims)) {
    
    sse <- sum((freq.k.sims[[i]]$y.tilde[, 1] - freq.k.sims[[i]]$y.test)^2)
    mse <- sse / length(freq.k.sims[[i]]$y.tilde)
    
    rmse.i <- sqrt(mse)
    
    
    k.rmse[[i]] <- rmse.i
    k.n[[i]] <- length(freq.k.sims[[i]]$y.tilde)
    
  }
  
  rmse <- weighted.mean(unlist(k.rmse), unlist(k.n))
  
  return(list(rmse = rmse, k.rmse = k.rmse, k.n = k.n))
  
}

#### Calculate RMSPE ####

## Define function that'll calculate the posterior predictive RMSPE

#k.sims <- lasso.k.sims

calculate_pp_k_rmspe <- function(k.sims) {
  
  k.error <- rep(as.numeric(NA), length(k.sims))
  
  k.n <- rep(as.numeric(NA), length(k.sims))
  
  for (i in 1:length(k.sims)) {
    
    k.sim <- k.sims[[i]]
    
    y.test <- k.sim$y.test
    
    # Get posterior predictives
    
    y.tilde.matrix <- k.sim$y.tilde.matrix
    
    y.tilde.error <- sweep(y.tilde.matrix, 2, y.test)
    y.tilde.error <- y.tilde.error^2
    
    y.tilde.error <- apply(y.tilde.error, 1, sum)
    y.tilde.error <- y.tilde.error / k.sim$n.tilde
    y.tilde.error <- sqrt(y.tilde.error)
    
    k.error[i] <- mean(y.tilde.error)
    
    k.n[i] <- k.sim$n.tilde
    
  }
  
  rmse <- weighted.mean(k.error, k.n)
  
  return(list(rmse = rmse, k.rmse = k.error, k.n = k.n))
}

#### Define function that'll calculate the posterior MAP predictive RMSPE

calculate_map_k_rmspe <- function(k.sims) {
  
  k.error <- rep(as.numeric(NA), length(k.sims))
  
  k.n <- rep(as.numeric(NA), length(k.sims))
  
  for (i in 1:length(k.sims)) {
    
    sim <- k.sims[[i]]
    
    ### Extract X tilde
    
    X.tilde <- sim$X.test
    
    ### Use estimate_parameters to get MAP estimates
    
    b.ests <- estimate_parameters(sim, lambda = FALSE, sigma_sq = FALSE)
    b.ests <- b.ests$par.summary$map
    
    ### Estimate y tilde
    
    y.tilde <- X.tilde %*% b.ests
    
    ### Compare to y.test
    
    sq.error <- (y.tilde - sim$y.test)^2 
    
    k.error[i] <- sqrt(sum(sq.error) / length(sim$y.test))
    
    k.n[i] <- length(sim$y.test)
    
  }
  
  map.rmse <- weighted.mean(k.error, k.n)
  
  return(list(rmse = map.rmse, k.rmse = k.error, k.n = k.n))
  
}

#### Load models

ols.k.sims <- readRDS(paste0(model.dump.string, "ols_k_sims.RDS"))

freq.ridge.k.sims <- readRDS(paste0(model.dump.string, "freq_ridge_k_sims.RDS"))

freq.lasso.k.sims <- readRDS(paste0(model.dump.string, "freq_lasso_k_sims.RDS"))

blm.k.sims <- readRDS(paste0(model.dump.string, "blm_k_sims.RDS"))

ridge.k.sims <- readRDS(paste0(model.dump.string, "ridge_k_sims.RDS"))

lasso.k.sims <- readRDS(paste0(model.dump.string, "lasso_k_sims.RDS"))

#### OLS ####

ols.rmspe <- get_ols_rmspe(ols.k.sims)

#ols.rmspe

#### Freq Ridge ####

freq.ridge.rmspe <- get_freq_ea_rmspe(freq.ridge.k.sims)

# freq.ridge.rmspe

#### Freq LASSO ####

freq.lasso.rmspe <- get_freq_ea_rmspe(freq.lasso.k.sims)

# freq.lasso.rmspe

#### BLM ####

blm.pp.rmspe <- calculate_pp_k_rmspe(blm.k.sims)

blm.map.rmspe <- calculate_map_k_rmspe(blm.k.sims)

# blm.pp.rmspe
# blm.map.rmspe

#### Ridge ####

ridge.pp.rmspe <- calculate_pp_k_rmspe(ridge.k.sims)

ridge.map.rmspe <- calculate_map_k_rmspe(ridge.k.sims)

# ridge.pp.rmspe
# ridge.map.rmspe

#### LASSO

lasso.pp.rmspe <- calculate_pp_k_rmspe(lasso.k.sims)

lasso.map.rmspe <- calculate_map_k_rmspe(lasso.k.sims)

# lasso.pp.rmspe
# lasso.map.rmspe

# Compare models

mod.string <- c("Null", "OLS", "BLM", "Freq ridge", "Bayesian Ridge", "Freq LASSO", "Bayesian LASSO")

map.rmspe <- c(null.rmspe$mean, 
               ols.rmspe$rmse, blm.map.rmspe$rmse, 
               freq.ridge.rmspe$rmse, ridge.map.rmspe$rmse, 
               freq.lasso.rmspe$rmse, lasso.map.rmspe$rmse)

pp.rmspe <- c(NA, NA, blm.pp.rmspe$rmse, NA, ridge.pp.rmspe$rmse, NA, lasso.pp.rmspe$rmse)

d <- data.table(model = mod.string,
                MAP_RMSPE = map.rmspe,
                PP_RMSPE = pp.rmspe)

d

# Build box-plot

build_box_plot <- function(null.rmspe, ols.rmspe, 
                           freq.ridge.rmspe, freq.lasso.rmspe,
                           blm.map.rmspe, ridge.map.rmspe, lasso.map.rmspe) {
  
  null.d <- data.table(rmspe = unlist(null.rmspe$y.err / null.rmspe$y.n))
  null.d[, model := "Null"]
  
  ols.d <- data.table(rmspe = unlist(ols.rmspe$k.rmse))
  ols.d[, model := "OLS"]
  
  f.ridge.d <- data.table(rmspe = unlist(freq.ridge.rmspe$k.rmse))
  f.ridge.d[, model := "Freq. Ridge"]
  
  f.lasso.d <- data.table(rmspe = unlist(freq.lasso.rmspe$k.rmse))
  f.lasso.d[, model := "Freq. LASSO"]
  
  blm.d <- data.table(rmspe = unlist(blm.map.rmspe$k.rmse))
  blm.d[, model := "BLM"]
  
  ridge.d <- data.table(rmspe = unlist(ridge.map.rmspe$k.rmse))
  ridge.d[, model := "Ridge"]
  
  lasso.d <- data.table(rmspe = unlist(lasso.map.rmspe$k.rmse))
  lasso.d[, model := "LASSO"]
  
  d <- rbind(null.d,
             ols.d,
             f.ridge.d,
             f.lasso.d,
             blm.d,
             ridge.d,
             lasso.d)
  
  d[, model := factor(model, 
                      levels = c("Null", "OLS", 
                                 "Freq. Ridge", "Freq. LASSO", 
                                 "BLM", "Ridge", "LASSO"))]
  
  d.summary <- data.table(model = c("Null", "OLS", 
                                    "Freq. Ridge", "Freq. LASSO", 
                                    "BLM", "Ridge", "LASSO"),
                          mean_rmspe = c(null.rmspe$mean,
                                         ols.rmspe$rmse,
                                         freq.ridge.rmspe$rmse,
                                         freq.lasso.rmspe$rmse,
                                         blm.map.rmspe$rmse,
                                         ridge.map.rmspe$rmse,
                                         lasso.map.rmspe$rmse))
  
  dp <- ggplot(d, aes(x = model, y = rmspe)) + 
    geom_boxplot() +
    geom_point(data = d.summary, aes(y = mean_rmspe), shape = 21, fill = "white", colour = "red", size = 3) +
    theme_bw(base_size = 12) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("") + ylab("RMSPE")
  
  return(list(dp = dp, d.summary = d.summary))
  
}

rmspe.box.plot <- build_box_plot(null.rmspe, ols.rmspe, 
                                 freq.ridge.rmspe, freq.lasso.rmspe,
                                 blm.map.rmspe, ridge.map.rmspe, lasso.map.rmspe)  




