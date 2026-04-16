data {
  int<lower=1> N;
  vector[N] y;
  array[N] int<lower=0, upper=1> male;
  vector[N] bw_z;
  vector[N] bw_z2;
  vector[N] age_z;
  vector[N] age_z2;

  int<lower=1> N_pred;
  array[N_pred] int<lower=0, upper=1> male_pred;
  vector[N_pred] bw_z_pred;
  vector[N_pred] bw_z2_pred;
  vector[N_pred] age_z_pred;
  vector[N_pred] age_z2_pred;

  int<lower=1> N_ppc;
  array[N_ppc] int<lower=0, upper=1> male_ppc;
  vector[N_ppc] bw_z_ppc;
  vector[N_ppc] bw_z2_ppc;
  vector[N_ppc] age_z_ppc;
  vector[N_ppc] age_z2_ppc;

  real alpha_mean;
  real<lower=0> alpha_sd;
  real<lower=0> beta_male_sd;
  real<lower=0> beta_bw_sd;
  real<lower=0> beta_bw2_sd;
  real<lower=0> beta_age_sd;
  real<lower=0> beta_age2_sd;
  real<lower=0> sigma_rate;
  real<lower=0> nu_offset;
  real<lower=0> nu_rate;
}

parameters {
  real alpha;
  real beta_male;
  real beta_bw;
  real beta_bw2;
  real beta_age;
  real beta_age2;
  real<lower=0> sigma;
  real<lower=0> nu_shift;
}

transformed parameters {
  real<lower=nu_offset> nu = nu_offset + nu_shift;
}

model {
  vector[N] mu;

  mu = alpha +
    beta_male * to_vector(male) +
    beta_bw * bw_z +
    beta_bw2 * bw_z2 +
    beta_age * age_z +
    beta_age2 * age_z2;

  alpha ~ normal(alpha_mean, alpha_sd);
  beta_male ~ normal(0, beta_male_sd);
  beta_bw ~ normal(0, beta_bw_sd);
  beta_bw2 ~ normal(0, beta_bw2_sd);
  beta_age ~ normal(0, beta_age_sd);
  beta_age2 ~ normal(0, beta_age2_sd);
  sigma ~ exponential(sigma_rate);
  nu_shift ~ exponential(nu_rate);

  y ~ student_t(nu, mu, sigma);
}

generated quantities {
  vector[N_pred] mu_pred;
  vector[N_pred] y_pred;
  vector[N_ppc] mu_ppc;
  vector[N_ppc] y_ppc;

  for (n in 1:N_pred) {
    mu_pred[n] = alpha +
      beta_male * male_pred[n] +
      beta_bw * bw_z_pred[n] +
      beta_bw2 * bw_z2_pred[n] +
      beta_age * age_z_pred[n] +
      beta_age2 * age_z2_pred[n];

    y_pred[n] = student_t_rng(nu, mu_pred[n], sigma);
  }

  for (n in 1:N_ppc) {
    mu_ppc[n] = alpha +
      beta_male * male_ppc[n] +
      beta_bw * bw_z_ppc[n] +
      beta_bw2 * bw_z2_ppc[n] +
      beta_age * age_z_ppc[n] +
      beta_age2 * age_z2_ppc[n];

    y_ppc[n] = student_t_rng(nu, mu_ppc[n], sigma);
  }
}
