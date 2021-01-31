
// Benchmarking to Effect Size, optimistic analysis prior 
// Date: 01 Sep 2019 
// Edit: 20 Feb 2020 - changed to pressimistic, informative prior for ATE.

// 1. Data Block
data {
  int<lower=0> N;  // number of observations
  real y[N];  // outcome variable
  vector[N] treat;  // treatment indicator
}

// 2. Parameters
parameters {
  real alpha;
  real beta;
}

// 3. transformed parameters (not applicable)

// 4. Model --> specify priors and likelihood for outcomes 
model {
  // priors
  alpha ~ normal(-1,10);  // our prior for alpha
  beta ~ normal(-0.1,0.3); // our prior for beta
  // log-likelihood
  y ~ normal(alpha + beta*treat,1);
} 

// in 4, we still assume fixed variance, can/should be more flexible. 

// Blank End


