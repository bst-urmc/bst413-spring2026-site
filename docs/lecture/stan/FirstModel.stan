  data { 
    real Y[10];
  } 
  transformed data { 
    vector[10] Ysq;
    for( i in 1:10) { 
      Ysq[i] = Y[i]^2; 
    } 
  } 
   
  parameters { 
    real mu; 
    real<lower=0> tau; 
  } 
   
  model { 
    Y ~ normal(mu , tau); 
  } 
