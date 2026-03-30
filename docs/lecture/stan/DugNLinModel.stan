  data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  real alpha;
  real<lower=0> beta;
  real<lower=0.5, upper=1> gamma;
  real<lower=0> sigma;
}
  transformed parameters { 
    vector[N] mu;
      mu =  alpha - beta * pow(gamma,x); 

  } 
model {
  y ~ normal(mu, sigma);
  alpha ~ cauchy(0,1000);
  beta ~ cauchy(0,1000);
  sigma ~ cauchy(0,1000);

} 
