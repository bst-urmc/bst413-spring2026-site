  data { 
    int<lower=0> N;   //flexible dimension
    real Y[N];        //observed data
    real<lower=0> tau; //known std dev
  } 

  parameters { 
    real mu; 
  } 
   
  model { 
    Y ~ normal(mu , tau); 
  } 
