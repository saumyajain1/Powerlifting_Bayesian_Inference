source("scripts/utils.R")

suppressPackageStartupMessages(library(cmdstanr))
suppressPackageStartupMessages(library(posterior))

stan_file <- "models/main_nonlinear_student_t.stan"
tables_dir <- "datasets/intermediate"
artifacts_dir <- "artifacts/main"

main_data <- readRDS(file.path(tables_dir, "stan_data_main.rds"))
model_df <- readRDS(file.path(tables_dir, "model_data.rds"))
prediction_grid <- readRDS(file.path(tables_dir, "prediction_grid.rds"))
ppc_df <- baseline_ppc_subset(model_df)
stan_data <- main_stan_data(main_data, prediction_grid, ppc_df)

compile_dir <- file.path(artifacts_dir, "cmdstan")
hmc_dir <- file.path(artifacts_dir, "hmc")
meanfield_dir <- file.path(artifacts_dir, "vi_meanfield")
fullrank_dir <- file.path(artifacts_dir, "vi_fullrank")
hmc_fit_path <- file.path(artifacts_dir, "main_hmc_fit.rds")
meanfield_fit_path <- file.path(artifacts_dir, "main_vi_meanfield_fit.rds")
fullrank_fit_path <- file.path(artifacts_dir, "main_vi_fullrank_fit.rds")

ensure_dir(tables_dir)
ensure_dir(compile_dir)
ensure_dir(hmc_dir)
ensure_dir(meanfield_dir)
ensure_dir(fullrank_dir)

model <- cmdstan_model(stan_file, compile = FALSE)
model$compile(dir = compile_dir, quiet = TRUE)

if (file.exists(hmc_fit_path)) {
  hmc_fit <- readRDS(hmc_fit_path)
} else {
  hmc_fit <- model$sample(
    data = stan_data,
    seed = 407,
    chains = 4,
    parallel_chains = min(4, parallel::detectCores()),
    iter_warmup = 500,
    iter_sampling = 500,
    adapt_delta = 0.95,
    refresh = 200,
    output_dir = hmc_dir,
    output_basename = "main_hmc"
  )
  saveRDS(hmc_fit, hmc_fit_path)
}

if (file.exists(meanfield_fit_path)) {
  meanfield_fit <- readRDS(meanfield_fit_path)
} else {
  meanfield_fit <- model$variational(
    data = stan_data,
    seed = 407,
    algorithm = "meanfield",
    iter = 10000,
    output_samples = 1000,
    refresh = 200,
    output_dir = meanfield_dir,
    output_basename = "main_vi_meanfield"
  )
  saveRDS(meanfield_fit, meanfield_fit_path)
}

if (file.exists(fullrank_fit_path)) {
  fullrank_fit <- readRDS(fullrank_fit_path)
} else {
  fullrank_fit <- model$variational(
    data = stan_data,
    seed = 408,
    algorithm = "fullrank",
    iter = 10000,
    output_samples = 1000,
    refresh = 200,
    output_dir = fullrank_dir,
    output_basename = "main_vi_fullrank"
  )
  saveRDS(fullrank_fit, fullrank_fit_path)
}

parameter_summary <- rbind(
  fit_parameter_summary(hmc_fit, main_parameter_vars, "HMC"),
  fit_parameter_summary(meanfield_fit, main_parameter_vars, "Mean-field VI"),
  fit_parameter_summary(fullrank_fit, main_parameter_vars, "Full-rank VI")
)

prediction_summary <- rbind(
  transform(
    summarize_generated(as_draws_matrix(hmc_fit$draws(variables = "mu_pred")), "mu_pred", prediction_grid),
    method = "HMC"
  ),
  transform(
    summarize_generated(as_draws_matrix(meanfield_fit$draws(variables = "mu_pred")), "mu_pred", prediction_grid),
    method = "Mean-field VI"
  ),
  transform(
    summarize_generated(as_draws_matrix(fullrank_fit$draws(variables = "mu_pred")), "mu_pred", prediction_grid),
    method = "Full-rank VI"
  )
)

write.csv(parameter_summary, file.path(tables_dir, "main_parameter_summary.csv"), row.names = FALSE)
write.csv(prediction_summary, file.path(tables_dir, "main_prediction_summary.csv"), row.names = FALSE)
