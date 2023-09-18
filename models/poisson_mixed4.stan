data {
  int<lower=0> N;
  int<lower=1> K; // number of fixed effect levels
  int<lower=1> J; // number of random effect levels
  array[N] int<lower=0> y;    // vector of data
  array[N] int<lower=1, upper=J> id;
  array[N] int<lower=1, upper=K> fixed;
}
parameters {
  vector[K] beta;
  vector[J] alpha;
  real<lower=0, upper = 10> kappa;
}

model {
  beta ~ normal(0, 10);
  alpha ~ normal(0, kappa);
  for (n in 1:N) {
    y[n] ~ poisson_log(alpha[id[n]] + beta[fixed[n]]);
  }
}

