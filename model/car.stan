data {
   int N;
   vector[N] degraded ;
   vector[N] deforested ;
   vector[N] climate ;
   matrix<lower=0>[N,N] spatial_weights ;
   matrix<lower=0,upper=1>[N,N] identity;
}
parameters {
   real intact ;
   real degrad ;
   real deforest ;
   real<lower = 0> sigma;
   real<lower=-1,upper=1> lambda;
}
model {
  climate ~ multi_normal_prec(intact + degrad * degraded + deforest * deforested,  
                              crossprod(identity - lambda * spatial_weights)/(sigma*sigma));
}

