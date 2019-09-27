data {
  int<lower=1> K;
  int<lower=1> D;
  int<lower=1> W;
  int P;
  matrix<lower=0>[D,W] count_matrix;
  vector[P] x[D];
}
parameters {
  simplex[W] phi[K];    
  vector<lower=0>[W] eta;
  vector<lower=0, upper=1>[K] v;
  real<lower=0> sigma;
  real<lower=0> gamma;
  real<lower=0> xi;
  vector[P] lambda[K];
  vector<lower=0, upper=1>[K] pi_dash[D];
} 
transformed parameters {
  vector[K] alpha[D];
  simplex[K] pi[D];
  vector[K] beta;
  vector[D] sum_pi_dash;
  beta[1]=v[1];
  for (k in 2:K){
    beta[k]=v[k]*prod(1-(v[1:(k-1)]));
  }
  for (d in 1:D){
    for (k in 1:K){
      alpha[d,k] = sum(x[d,1:P].*lambda[k,1:P]);
      alpha[d,k] = exp(alpha[d,k]);
    }
  } 
  for (d in 1:D){
    sum_pi_dash[d]=sum(pi_dash[d,1:K]);
    for (k in 1:K){
      pi[d,k]=pi_dash[d,k]/sum_pi_dash[d];
    }
  }
}
model {
  for (k in 1:K){
    for (p in 1:P){
      lambda[k,p] ~ normal(0, sigma);
    }
  }
  for (k in 1:K){
    v[k]~beta(1,gamma);
  }
  for (k in 1:K){
    for (d in 1:D){
      pi_dash[d,k]~gamma(xi*beta[k],alpha[d,k]); 
    }
  }
  for (k in 1:K)
    phi[k] ~ dirichlet(eta);
  for (d in 1:D) {
    for (w in 1:W) {
      real L[K];
      for (k in 1:K)
        L[k] = count_matrix[d,w]*(log(pi[d,k]) + log(phi[k,w]) );
      target+=log_sum_exp(L);
    }
  }
}