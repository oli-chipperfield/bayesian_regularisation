
  
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

    
