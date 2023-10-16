data {
  int<lower=0> N;
  int<lower=1> K; // number of fixed effect levels
  array[N] int<lower=0> y;    // vector of data
  array[N] int<lower=1, upper=K> fixed;
}
parameters {
  vector[K] beta;
}

model {
  beta ~ normal(0, 10);
  for (n in 1:N) {
    y[n] ~ poisson_log(beta[fixed[n]]);
  }
}

