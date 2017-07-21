data {
  int<lower=1> nvals;
  int<lower=1> L;
  int<lower=1,upper=nvals> N;
  int<lower=1,upper=nvals> G;
  int<lower=1,upper=N> nn[nvals];
  int<lower=1,upper=G> gg[nvals];
  vector[nvals] y;
  real<lower=0> sw;
  real<lower=0> sv_rate; //controls ARD, large value=strong shrinkage
  real ymn;
  //mnar part
  int<lower=0,upper=1> Z[N,G]; //indicator of observed data, transpose from Y indices
  vector<lower=0,upper=1>[N] Q; //detection rates
  real b1_mn; #.5 is reasonable
  real<lower=0> b1_sd; #.1 is reasonable
}
transformed data {
  row_vector[N] b0_mn;
  b0_mn = logit(Q') - b1_mn*ymn;
}
parameters {
  real y0_raw;
  matrix[L,N] U;
  matrix[G,L] V_raw;
  vector[G] w_raw;
  real<lower=0> sy;
  vector<lower=0>[L] sv;
  row_vector[N] b0_raw; #missingness intercept
  row_vector[N] b1_raw; #missingness slope
}
transformed parameters {
  real y0;
  vector[G] w;
  matrix[G,L] V;
  row_vector[N] b0;
  row_vector[N] b1;
  y0 = 5*y0_raw + ymn; #super flat prior for intercept term
  w = sw*w_raw;
  V = diag_post_multiply(V_raw,sv);
  b0 = b0_mn+b0_raw;
  b1 = b1_mn+b1_sd*b1_raw;
}
model {
  matrix[G,N] eta; //imputed estimates for 'latent expression'
  vector[nvals] y_imp;
  // prior
  w_raw~normal(0,1);
  y0_raw ~ cauchy(0,1);
  b0_raw ~ cauchy(0,1);
  b1_raw ~ normal(0,1);
  sy ~ cauchy(0,1); #half cauchy
  sv ~ gamma(2,sv_rate);
  for(l in 1:L){
    V_raw[:,l] ~ normal(0,1);
  }
  for(n in 1:N){
    U[:,n] ~ normal(0,1);
  }
  // intermediate variables
  eta = y0 + rep_matrix(w,N) + V*U;
  for (val in 1:nvals) {
    y_imp[val] = eta[gg[val],nn[val]];
  }
  // likelihood
  for(n in 1:N){
    Z[n,:] ~ bernoulli_logit(b0[n]+b1[n]*eta[:,n]);
  }
  y ~ normal(y_imp, sy);
}
