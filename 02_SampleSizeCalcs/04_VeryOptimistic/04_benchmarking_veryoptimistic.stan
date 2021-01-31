
// Benchmarking to Effect Size, very optimistic analysis prior 
// Date: 05 DEC 2020 

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
  beta ~ normal(0.3,0.1); // our prior for beta
  // log-likelihood
  y ~ normal(alpha + beta*treat,1);
} 

// in 4, we still assume fixed variance, can/should be more flexible. 

// Blank End

