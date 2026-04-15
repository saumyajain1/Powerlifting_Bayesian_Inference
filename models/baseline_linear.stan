data {
  int<lower=1> N;
  vector[N] y;
  array[N] int<lower=0, upper=1> male;
  vector[N] bw_z;
  vector[N] age_z;

  int<lower=1> N_pred;
  array[N_pred] int<lower=0, upper=1> male_pred;
  vector[N_pred] bw_z_pred;
  vector[N_pred] age_z_pred;

  int<lower=1> N_ppc;
  array[N_ppc] int<lower=0, upper=1> male_ppc;
  vector[N_ppc] bw_z_ppc;
  vector[N_ppc] age_z_ppc;

  real alpha_mean;
  real<lower=0> alpha_sd;
  real<lower=0> beta_male_sd;
  real<lower=0> beta_bw_sd;
  real<lower=0> beta_age_sd;
  real<lower=0> sigma_rate;
}

parameters {
  real alpha;
  real beta_male;
  real beta_bw;
  real beta_age;
  real<lower=0> sigma;
}

model {
  vector[N] mu;

  mu = alpha + beta_male * to_vector(male) + beta_bw * bw_z + beta_age * age_z;

  alpha ~ normal(alpha_mean, alpha_sd);
  beta_male ~ normal(0, beta_male_sd);
  beta_bw ~ normal(0, beta_bw_sd);
  beta_age ~ normal(0, beta_age_sd);
  sigma ~ exponential(sigma_rate);

  y ~ normal(mu, sigma);
}

generated quantities {
  vector[N_pred] mu_pred;
  vector[N_pred] y_pred;
  vector[N_ppc] mu_ppc;
  vector[N_ppc] y_ppc;

  for (n in 1:N_pred) {
    mu_pred[n] = alpha + beta_male * male_pred[n] + beta_bw * bw_z_pred[n] + beta_age * age_z_pred[n];
    y_pred[n] = normal_rng(mu_pred[n], sigma);
  }

  for (n in 1:N_ppc) {
    mu_ppc[n] = alpha + beta_male * male_ppc[n] + beta_bw * bw_z_ppc[n] + beta_age * age_z_ppc[n];
    y_ppc[n] = normal_rng(mu_ppc[n], sigma);
  }
}
