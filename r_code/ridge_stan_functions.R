#source("data_prep.R")

#### Define Ridge model string

ridge_model_string <- function() {
  
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
      
      matrix[p, p] id_p_matrix;                // Identity matrix for Sigma matrix
  
      vector[p] beta_zero;                     // Prior betas

      id_p_matrix = diag_matrix (rep_vector (1, p));

      beta_zero = rep_vector (0, p);

    }

    parameters {
  
      real<lower=0> sigma_sq;         // Sigma^2 parameter
  
      real<lower=0> lambda;           // lambda parameter  
  
      vector[p] beta;                 // beta parameters
  
    }

    model {

      // Priors

      lambda ~ gamma (r, d);

      sigma_sq ~ inv_gamma (a, b);

      beta ~ multi_normal (beta_zero, (sigma_sq / lambda) * id_p_matrix); 

      // Likelhood

      y ~ normal (X * beta, sigma_sq);

    }

    generated quantities {
  
      vector[n_tilde] y_tilde;        // Predicted y values
  
      y_tilde = X_tilde * beta;
  
    }

    "
  return(string)
  
  }

ridge.model.string <- ridge_model_string()

#### Save string as .stan file

write(ridge.model.string, "ridge.stan")

#### Define function to run model
# 
# y <- y.trains[[1]] - mean(y.trains[[1]])
# X <- X.trains[[1]]
# y.test <- y.tests[[1]] - mean(y.tests[[1]])
# X.test <- X.tests[[1]]
# # a <- 1
# # b <- 1
# # r <- 1
# # d <- 1
# 
# file.string <- "ridge_fixed_lambda.stan"
# 
# pars <- list(a = 1, b = 1, lambda = 4)

ridge_stan_simulation <- function(file.string, y, X, y.test = NULL, X.test = NULL, pars = NULL) {
  
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
  
  ridge.data <- list(n = n,
                     p = p,
                     n_tilde = n.tilde,
                     y = y,
                     X = X,
                     X_tilde = X.tilde)
  
  ridge.data <- c(ridge.data, pars)
  
  ridge.fit <- stan(file = file.string, data = ridge.data, cores = 4, iter = 3000, warmup = 1000)

  # Compile beta parameters
  
  beta.list <- list()
  
  for (i in 1:p) {
    
    beta.list[[i]] <- rstan::extract(ridge.fit, paste0("beta[",i,"]"))[[1]]    
    
  }
  
  max.m <- max(unlist(lapply(beta.list, length)))
  
  beta.matrix <- matrix(rep(as.numeric(NA), max.m * p), max.m, p)
  
  for (i in 1:p) {
    
    beta.matrix[, i] <- beta.list[[i]]
    
  }
  
  names(beta.matrix) <- colnames(X)  
  
  # Extract sigma sq
  
  sigma.sq.draws <- rstan::extract(ridge.fit, "sigma_sq")[[1]]
  
  # Extract lambda
  
  lambda.draws <- tryCatch(rstan::extract(ridge.fit, "lambda")[[1]],
                           error = function(err) {NULL})
                           
  # Extract posterior predictive y's
  
  y.tilde.list <- list()
  
  for (i in 1:n.tilde) {
    
    y.tilde.list[[i]] <- rstan::extract(ridge.fit, paste0("y_tilde[",i,"]"))[[1]]    
    
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
                 model.object = ridge.fit)  
  
  r.list <- c(r.list, pars)
  
  return(r.list)
  
  }




