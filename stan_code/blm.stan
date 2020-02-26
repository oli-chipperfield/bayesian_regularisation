
  
  data {
  
    int<lower=0> n;             // Observations 
    int<lower=0> p;             // Parameters
  
    int<lower=0> n_tilde;       // Prediction n
  
    vector[n] y;                // Y vector
  
    matrix[n, p] X;             // X (design) matrix
  
    matrix[n_tilde, p] X_tilde; // Testing n

    real<lower=0> a;            // Sigma sq alpha paramater

    real<lower=0> b;            // Sigma sq beta parameter
  
  }
  
  parameters {
  
    vector[p] beta;                 // beta parameters
  
    real<lower=0> sigma_sq;         // Sigma^2 parameter
  
  }    
  
  model {
  
  // Priors
  
  sigma_sq ~ inv_gamma (a, b);
  
  // Likelihood
  
  y ~ normal (X * beta, sigma_sq);
  
  }
  
  generated quantities {
  
  vector[n_tilde] y_tilde;        // Predicted y values
  
  y_tilde = X_tilde * beta;
  
  }
  
  
