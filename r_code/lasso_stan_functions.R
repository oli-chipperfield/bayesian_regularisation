#source("data_prep.R")

#### Define LASSO model string

lasso_model_string <- function() {
  
  string <- "
  
  data {
  
  int<lower=0> n;             // Observations 
  int<lower=0> p;             // Parameters

  int<lower=0> n_tilde;       // Prediction n
  
  vector[n] y;                // Y vector
  
  matrix[n, p] X;             // X (design) matrix

  matrix[n_tilde, p] X_tilde; // Testing X
  
  real<lower=0> a;            // Sigma sq alpha paramater
  real<lower=0> b;            // Sigma sq beta parameter
  
  real<lower=0> r;            // Lambda r parameter
  real<lower=0> d;            // Lambda d parameter
  
  }
  
  transformed data {
  
  vector[p] beta_zero;                     // Prior betas
  
  beta_zero = rep_vector (0, p);
  
  }
  
  parameters {
  
  real<lower=0> lambda;           // lambda parameter  
  
  real<lower=0> sigma_sq;         // sigma sq parameter  

  vector<lower=0>[p] tau_sq;               // tau sq parameters  

  vector[p] beta;                 // beta parameters
  
  }
  
  model {
  
  // Priors
  
  lambda ~ gamma (r, d);
  
  sigma_sq ~ inv_gamma (a, b);

  tau_sq ~ exponential (0.5 * lambda);

  beta ~ multi_normal (beta_zero, sigma_sq * diag_matrix(tau_sq));

  // Likelihood

  y ~ normal (X * beta, sigma_sq);
  
  }

  generated quantities {
  
    vector[n_tilde] y_tilde;        // Predicted y values
  
    y_tilde = X_tilde * beta;
  
  }
  
  "
  
  return(string)
  
}

lasso.model.string <- lasso_model_string()

#### Save string as .stan file

write(lasso.model.string, "lasso.stan")

#### Define function to run model

# y <- y.trains[[1]] - mean(y.trains[[1]])
# X <- X.trains[[1]]
# y.test <- y.tests[[1]] - mean(y.tests[[1]])
# X.test <- X.tests[[1]]
# a <- 1
# b <- 1
# r <- 1
# d <- 1

# file.string <- "lasso_fixed_lambda.stan"
# # 
# pars <- list(a = 0.001, b = 0.001, lambda = 25)

lasso_stan_simulation <- function(file.string, y, X, y.test = NULL, X.test = NULL, pars = NULL) {
  
  # Define constants
  
  n <- length(y)
  
  p <- ncol(X)
  
  y.test <- if(is.null(y.test)) {y[1:10]} else {y.test}  # If there's no testing just define
  
  n.tilde <- length(y.test)
  
  X.tilde <- if(is.null(X.test)) {X[1:10,]} else {X.test}
  
  # Recentre X and y if using k-fold validation
  
  # Recentre X
  
  X.means <- apply(X, 2, mean)
  
  X <- sweep(X, 2, X.means)
  
  X.tilde <- sweep(X.tilde, 2, X.means)
  
  # Recentre Y
  
  y.bar <- mean(y)
  
  y <- y - y.bar
  
  y.test <- y.test - y.bar
  
  lasso.data <- list(n = n,
                     p = p,
                     n_tilde = n.tilde,
                     y = y,
                     X = X,
                     X_tilde = X.tilde)
  
  lasso.data <- c(lasso.data, pars)

  lasso.fit <- stan(file = file.string, data = lasso.data, cores = 4, iter = 3000, warmup = 1000)

  # Compile beta parameters
  
  beta.list <- list()
  
  for (i in 1:p) {
    
    beta.list[[i]] <- rstan::extract(lasso.fit, paste0("beta[",i,"]"))[[1]]    
    
  }
  
  max.m <- max(unlist(lapply(beta.list, length)))
  
  beta.matrix <- matrix(rep(as.numeric(NA), max.m * p), max.m, p)
  
  for (i in 1:p) {
    
    beta.matrix[, i] <- beta.list[[i]]
    
  }
  
  names(beta.matrix) <- colnames(X)  
  
  # Extract sigma sq
  
  sigma.sq.draws <- rstan::extract(lasso.fit, "sigma_sq")[[1]]
  
  # Extract lambda
  
  lambda.draws <- tryCatch(rstan::extract(lasso.fit, "lambda")[[1]],
                           error = function(err) {NULL})
  
  # Extract posterior predictive y's
  
  y.tilde.list <- list()
  
  for (i in 1:n.tilde) {
    
    y.tilde.list[[i]] <- rstan::extract(lasso.fit, paste0("y_tilde[",i,"]"))[[1]]    
    
  }
  
  max.y.m <- max(unlist(lapply(y.tilde.list, length)))
  
  y.tilde.matrix <- matrix(rep(as.numeric(NA), max.y.m * n.tilde), max.y.m, n.tilde)
  
  for (i in 1:n.tilde) {
    
    y.tilde.matrix[, i] <- y.tilde.list[[i]]
    
  }
  
  r.list <- list(n = n,
                 p = p,
                 y = y,
                 X = X,
                 n.tilde = n.tilde,
                 y.test = y.test,
                 X.test = X.test,
                 beta.matrix = beta.matrix,
                 sigma.sq.draws = sigma.sq.draws,
                 lambda.draws = lambda.draws,
                 y.tilde.matrix = y.tilde.matrix,
                 X.test = X.test,
                 model.object = lasso.fit)
  
  r.list <- c(r.list, pars)
  
  return(r.list)
  
  }

 