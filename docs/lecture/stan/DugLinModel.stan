  data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
  transformed data { 
    vector[N] logage = log(x);
  } 
parameters {
  real beta0;
  real beta1;
  real<lower=0> sigma;
}
  transformed parameters { 
    vector[N] mu;
      mu =  beta0+ beta1*logage; 
  } 
model {
  y ~ normal(mu, sigma);
  sigma ~ cauchy(0,1000);
} 
