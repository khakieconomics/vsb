  // our Stan model, saved as vsb.stan
// first we define the function that takes data and parameters and returns predicted market shares
functions {
  // calculates shares for a given market
  vector shares(real alpha, vector beta, matrix bigX_plain, matrix bigX_rc, matrix Sigma, vector xi, matrix z) {
    matrix[rows(z), rows(xi)+1] utilities; // NS x J_t + 1
    matrix[rows(z), rows(xi)+1] probs;
    vector[rows(xi)+1] shares;
    // 1. Rather than pass in p and x separately, we'll pass in bigX = append_col(p, X)
    // 2. append alpha_shock, beta_shock
    {
      matrix[rows(z), rows(xi)] tmp; // NS x J_t
      
      tmp = rep_matrix((bigX_plain*append_row(alpha, beta) + xi)', rows(z));
      
      // replace the append_col wing single values (might speed things up)
      utilities[1:rows(z), 1:rows(xi)] = tmp + z * cholesky_decompose(Sigma)' * bigX_rc';
      utilities[1:rows(z),(rows(xi)+1)] = rep_vector(0.0, rows(z)); //outside good has 0 utility
      
      for(i in 1:rows(z)) {
         probs[i] = softmax(utilities[i]')'; // Doing some unnecessary transposes here.
      }
      
    }
    
    for(j in 1:cols(probs)) {
      shares[j] = mean(col(probs, j))'; // Second unnecessary transpose to match the first
    }
    
    return(shares);
  }
}
// next define our data
data {
  int NS; // number of individuals in integration
  int T; // number of markets
  int J_m; // total number of products (\sum_t {J_t})
  int P_rc; // number of random coefficient features
  int P_plain; // number of non-random coefficient features
  int P_iv; // number of IV features

  int mkt_size[T]; // vector of market sizes for each mkt
  int mkt_numProds[T]; // number of products in each market
  int mkt_id[J_m]; // index for which market a product is in; this is 'cdid' in other blp scripts
  
  
  vector[J_m] price; // price for each unit
  int sales[J_m + T]; // unit sales across T markets for J products (including outside good)
  matrix[J_m, P_rc] X_rc; // Xs stacked on top of each other.
  matrix[J_m, P_plain] X_plain; // Xs stacked on top of each other.
  matrix[J_m, P_iv] Z_iv; // Zs stacked on top of each other.

  real<lower=0,upper=1> run_estimation;
  
  matrix[NS, P_rc+1] z; // normal(0,1) draws of the shocks
  real nu;
  int P0; // prior number of nonzero coefficients in price model
}
// next join the product data together into single matrices
transformed data {
  matrix[J_m, P_rc+1] bigX_rc; // X with price as the first column
  matrix[J_m, P_plain+1] bigX_plain; // X with price as the first column

  int mkt_start_indices[T];
  int mkt_end_indices[T];
  int last_index_seen;
  
  bigX_rc = append_col(price, X_rc);
  bigX_plain = append_col(price, X_plain);

  last_index_seen = 0;
    for(t in 1:T) {
      mkt_start_indices[t] = last_index_seen + 1;
      mkt_end_indices[t] = mkt_start_indices[t] + mkt_numProds[t]; //includes outside good
      last_index_seen = mkt_end_indices[t];
    }
}
// define parameters
parameters {
  real alpha; 
  vector[P_plain] beta;
  vector[P_iv + P_plain] mu;
  vector<lower = 0>[P_iv + P_plain] lambda2;
  real<lower = 0> price_scale;
  vector[J_m] xi;
  vector<lower = 0>[P_rc+1] scales;
  corr_matrix[P_rc+1] Omega;
  real<lower = 0> lambda;
}

transformed parameters {
  cov_matrix[P_rc+1] Sigma;
  real tau;
  Sigma = quad_form_diag(Omega, scales);
  tau = P0*pow(P_plain + P_iv - P0, -1)*((price_scale)/sqrt(J_m));
}
// and the model
model {
  vector[J_m + T] pred_shares;
  // priors
  alpha ~ normal(0, 1);
  beta[1] ~ cauchy(-10, 2);
  beta[2:P_plain] ~ normal(0, 1);
  mu ~ normal(0,tau*lambda2);
  lambda2 ~ cauchy(0, 1);
  lambda ~ normal(.5, 1);
  to_vector(xi) ~ normal(0, 1);
  scales ~ normal(0, 1);
  Omega ~ lkj_corr(nu);
  price_scale ~ normal(0, .5);
  // model of price -- this helps pin down xi
  for (i in 1:rows(price)) {
    log(price[i]) ~ normal(append_col(X_plain[i], Z_iv[i]) * mu + xi[i], price_scale);
  }
  
  // model of sales 
  {
    for(t in 1:T) {
      pred_shares[mkt_start_indices[t] : mkt_end_indices[t]] = shares(alpha, 
                                                                      beta, 
                                                                      bigX_plain[(mkt_start_indices[t] - (t-1) ) : (mkt_end_indices[t] - t)],
                                                                      bigX_rc[(mkt_start_indices[t] - (t-1) ) : (mkt_end_indices[t] - t)], 
                                                                      Sigma, 
                                                                      lambda *xi[(mkt_start_indices[t] - (t-1) ) : (mkt_end_indices[t] - t)],
                                                                      z); //includes the outisde good

      //print("sales", sales[mkt_start_indices[t] : mkt_end_indices[t]]);
      //print("xi ", xi[(mkt_start_indices[t] - (t-1) ) : (mkt_end_indices[t] - t)]);
      //print("beta", beta);
      //print("pred_shares", pred_shares[mkt_start_indices[t] : mkt_end_indices[t]]);
      sales[mkt_start_indices[t] : mkt_end_indices[t]] ~ multinomial(pred_shares[mkt_start_indices[t] : mkt_end_indices[t]]);

    }
  }
}
