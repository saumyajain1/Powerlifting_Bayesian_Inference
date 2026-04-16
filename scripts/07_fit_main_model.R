source("scripts/utils.R")

suppressPackageStartupMessages(library(cmdstanr))

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

fits <- list(
  "HMC" = load_or_sample_fit(
    model,
    stan_data,
    hmc_fit_path,
    hmc_dir,
    "main_hmc",
    seed = 407,
    adapt_delta = 0.95
  ),
  "Mean-field VI" = load_or_vi_fit(
    model,
    stan_data,
    meanfield_fit_path,
    meanfield_dir,
    "main_vi_meanfield",
    seed = 407,
    algorithm = "meanfield"
  ),
  "Full-rank VI" = load_or_vi_fit(
    model,
    stan_data,
    fullrank_fit_path,
    fullrank_dir,
    "main_vi_fullrank",
    seed = 408,
    algorithm = "fullrank"
  )
)

parameter_summary <- fit_list_parameter_summary(fits, main_parameter_vars)
prediction_summary <- fit_list_generated_summary(fits, "mu_pred", "mu_pred", prediction_grid)

write.csv(parameter_summary, file.path(tables_dir, "main_parameter_summary.csv"), row.names = FALSE)
write.csv(prediction_summary, file.path(tables_dir, "main_prediction_summary.csv"), row.names = FALSE)
