
  
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
  
  
