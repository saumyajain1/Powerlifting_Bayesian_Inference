source("scripts/utils.R")

suppressPackageStartupMessages(library(cmdstanr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(posterior))
suppressPackageStartupMessages(library(tidyr))

tables_dir <- "datasets/intermediate"
figures_dir <- "figures/baseline_checks"
artifacts_dir <- "artifacts/baseline"
recovery_dir <- file.path(artifacts_dir, "recovery")
recovery_fit_path <- file.path(recovery_dir, "baseline_recovery_hmc_fit.rds")
stan_file <- "models/baseline_linear.stan"

ensure_dir(figures_dir)
ensure_dir(recovery_dir)

exact_baseline_draws <- function(data_df, n_draws = 12000) {
  X <- cbind(1, data_df$male, data_df$bw_z, data_df$age_z)
  y <- data_df$TotalKg
  prior <- baseline_prior_data()
  prior_mean <- c(prior$alpha_mean, 0, 0, 0)
  prior_prec <- diag(1 / c(prior$alpha_sd^2, prior$beta_male_sd^2, prior$beta_bw_sd^2, prior$beta_age_sd^2))

  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  y_sq <- sum(y^2)
  n <- nrow(X)
  sigma_grid <- exp(seq(log(40), log(200), length.out = 1500))
  sigma_step <- c(diff(sigma_grid), tail(diff(sigma_grid), 1))

  log_post <- vapply(sigma_grid, function(sigma) {
    A <- prior_prec + XtX / sigma^2
    V <- chol2inv(chol(A))
    m <- V %*% (prior_prec %*% prior_mean + Xty / sigma^2)
    quad <- y_sq / sigma^2 + drop(crossprod(prior_mean, prior_prec %*% prior_mean)) - drop(crossprod(m, A %*% m))

    -n * log(sigma) - 0.5 * determinant(A, logarithm = TRUE)$modulus - 0.5 * quad - prior$sigma_rate * sigma
  }, numeric(1))

  weights <- exp(log_post + log(sigma_step) - max(log_post + log(sigma_step)))
  weights <- weights / sum(weights)
  sigma_draws <- sample(sigma_grid, n_draws, replace = TRUE, prob = weights)
  beta_draws <- matrix(NA_real_, n_draws, 4)

  for (i in seq_len(n_draws)) {
    sigma <- sigma_draws[i]
    A <- prior_prec + XtX / sigma^2
    V <- chol2inv(chol(A))
    m <- V %*% (prior_prec %*% prior_mean + Xty / sigma^2)
    beta_draws[i, ] <- drop(m + t(chol(V)) %*% rnorm(4))
  }

  out <- as.data.frame(beta_draws)
  names(out) <- c("alpha", "beta_male", "beta_bw", "beta_age")
  out$sigma <- sigma_draws
  out
}

draw_summary <- function(draw_df, method) {
  out <- data.frame(
    variable = names(draw_df),
    mean = vapply(draw_df, mean, numeric(1)),
    median = vapply(draw_df, median, numeric(1)),
    sd = vapply(draw_df, sd, numeric(1)),
    q5 = vapply(draw_df, quantile, numeric(1), probs = 0.05),
    q95 = vapply(draw_df, quantile, numeric(1), probs = 0.95),
    row.names = NULL
  )

  out <- align_summary_columns(out)
  out$method <- method
  out
}

ppc_density_df <- function(draw_matrix, observed, method, n_reps = 25) {
  keep <- sample(seq_len(nrow(draw_matrix)), n_reps)

  rbind(
    bind_rows(lapply(seq_along(keep), function(i) {
      dens <- density(draw_matrix[keep[i], ])
      data.frame(x = dens$x, y = dens$y, method = method, source = "Posterior predictive", rep = i)
    })),
    {
      dens <- density(observed)
      data.frame(x = dens$x, y = dens$y, method = method, source = "Observed data", rep = 0)
    }
  )
}

ppc_coverage_df <- function(draw_matrix, observed, method, levels = seq(0.1, 0.95, by = 0.05)) {
  bind_rows(lapply(levels, function(level) {
    alpha <- (1 - level) / 2
    lower <- apply(draw_matrix, 2, quantile, probs = alpha)
    upper <- apply(draw_matrix, 2, quantile, probs = 1 - alpha)

    data.frame(method = method, nominal = level, actual = mean(observed >= lower & observed <= upper))
  }))
}

hmc_fit <- readRDS(file.path(artifacts_dir, "baseline_hmc_fit.rds"))
meanfield_fit <- readRDS(file.path(artifacts_dir, "baseline_vi_meanfield_fit.rds"))
fullrank_fit <- readRDS(file.path(artifacts_dir, "baseline_vi_fullrank_fit.rds"))
model_df <- readRDS(file.path(tables_dir, "model_data.rds"))
prediction_summary <- read.csv(file.path(tables_dir, "baseline_prediction_summary.csv"), stringsAsFactors = FALSE)
ppc_df <- baseline_ppc_subset(model_df)

set.seed(405)
parameter_summary <- rbind(
  draw_summary(exact_baseline_draws(model_df[, c("TotalKg", "male", "bw_z", "age_z")]), "Exact"),
  fit_parameter_summary(hmc_fit, baseline_parameter_vars, "HMC"),
  fit_parameter_summary(meanfield_fit, baseline_parameter_vars, "Mean-field VI"),
  fit_parameter_summary(fullrank_fit, baseline_parameter_vars, "Full-rank VI")
)

write.csv(parameter_summary, file.path(tables_dir, "baseline_parameter_summary.csv"), row.names = FALSE)

trace_df <- as.data.frame(as_draws_df(hmc_fit$draws(variables = baseline_parameter_vars))) |>
  select(all_of(baseline_parameter_vars), .chain, .iteration) |>
  pivot_longer(cols = all_of(baseline_parameter_vars), names_to = "parameter", values_to = "value")

save_plot(
  ggplot(trace_df, aes(.iteration, value, color = factor(.chain), group = interaction(.chain, parameter))) +
    geom_line(alpha = 0.65, linewidth = 0.3) +
    facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
    labs(title = "Baseline HMC Trace Plots", x = "Iteration", y = "Value", color = "Chain") +
    theme_minimal(),
  file.path(figures_dir, "baseline_traceplot.png"),
  9,
  7
)

param_plot_df <- parameter_summary |>
  select(variable, mean, q5, q95, method) |>
  mutate(
    variable = factor(variable, levels = baseline_parameter_vars),
    method = factor(method, levels = c("Exact", "HMC", "Mean-field VI", "Full-rank VI"))
  )

save_plot(
  ggplot(param_plot_df, aes(mean, variable, color = method)) +
    geom_errorbarh(aes(xmin = q5, xmax = q95), height = 0.18, position = position_dodge(width = 0.55)) +
    geom_point(position = position_dodge(width = 0.55), size = 2) +
    labs(
      title = "Baseline Posterior Comparison",
      subtitle = "Semi-analytic reference against HMC and both VI approximations",
      x = "Estimate",
      y = NULL,
      color = "Method"
    ) +
    theme_minimal(),
  file.path(figures_dir, "baseline_exact_vs_hmc_vi.png"),
  8.5,
  5.2
)

hmc_y_ppc <- as_draws_matrix(hmc_fit$draws(variables = "y_ppc"))
meanfield_y_ppc <- as_draws_matrix(meanfield_fit$draws(variables = "y_ppc"))
fullrank_y_ppc <- as_draws_matrix(fullrank_fit$draws(variables = "y_ppc"))

save_plot(
  ggplot(
    bind_rows(
      ppc_density_df(hmc_y_ppc, ppc_df$TotalKg, "HMC"),
      ppc_density_df(meanfield_y_ppc, ppc_df$TotalKg, "Mean-field VI"),
      ppc_density_df(fullrank_y_ppc, ppc_df$TotalKg, "Full-rank VI")
    )
  ) +
    geom_line(
      data = function(x) subset(x, source == "Posterior predictive"),
      aes(x, y, group = interaction(method, rep)),
      color = "#94a3b8",
      alpha = 0.35,
      linewidth = 0.4
    ) +
    geom_line(
      data = function(x) subset(x, source == "Observed data"),
      aes(x, y),
      color = "#111827",
      linewidth = 1
    ) +
    facet_wrap(~ method, ncol = 1) +
    labs(
      title = "Baseline Posterior Predictive Check",
      subtitle = "Observed subset density against posterior predictive replicated densities",
      x = "Total (kg)",
      y = "Density"
    ) +
    theme_minimal(),
  file.path(figures_dir, "baseline_ppc_density.png"),
  8,
  7
)

coverage_df <- bind_rows(
  ppc_coverage_df(hmc_y_ppc, ppc_df$TotalKg, "HMC"),
  ppc_coverage_df(meanfield_y_ppc, ppc_df$TotalKg, "Mean-field VI"),
  ppc_coverage_df(fullrank_y_ppc, ppc_df$TotalKg, "Full-rank VI")
)

save_plot(
  ggplot(coverage_df, aes(nominal, actual, color = method)) +
    geom_abline(slope = 1, intercept = 0, color = "#9ca3af", linetype = "dashed") +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = "Baseline PPC Coverage",
      subtitle = "Actual central-interval coverage against nominal coverage",
      x = "Nominal coverage",
      y = "Actual coverage",
      color = "Method"
    ) +
    theme_minimal(),
  file.path(figures_dir, "baseline_ppc_coverage.png"),
  6.5,
  6
)

save_plot(
  ggplot(
    transform(
      prediction_summary,
      method = factor(method, levels = c("HMC", "Mean-field VI", "Full-rank VI")),
      Sex = factor(Sex, levels = c("F", "M"))
    ),
    aes(BodyweightKg, mean, color = Sex, fill = Sex)
  ) +
    geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.18, color = NA) +
    geom_line(linewidth = 0.9) +
    facet_grid(method ~ Age) +
    labs(
      title = "Baseline Predicted Strength Curves",
      subtitle = "Posterior mean and 90% intervals across bodyweight",
      x = "Bodyweight (kg)",
      y = "Predicted Total (kg)"
    ) +
    theme_minimal(),
  file.path(figures_dir, "baseline_prediction_curves.png"),
  10,
  7
)

recovery_df <- model_df[sample(seq_len(nrow(model_df)), 1500), c("male", "bw_z", "age_z")]
recovery_data <- c(
  list(
    N = nrow(recovery_df),
    y = with(recovery_df, rnorm(nrow(recovery_df), 370 + 173 * male + 62 * bw_z - 14 * age_z, 97)),
    male = recovery_df$male,
    bw_z = recovery_df$bw_z,
    age_z = recovery_df$age_z,
    N_pred = 1,
    male_pred = 0L,
    bw_z_pred = 0,
    age_z_pred = 0,
    N_ppc = 1,
    male_ppc = 0L,
    bw_z_ppc = 0,
    age_z_ppc = 0
  ),
  baseline_prior_data()
)

model <- cmdstan_model(stan_file, compile = FALSE)
model$compile(dir = file.path(artifacts_dir, "cmdstan"), quiet = TRUE)

if (file.exists(recovery_fit_path)) {
  recovery_fit <- readRDS(recovery_fit_path)
} else {
  recovery_fit <- model$sample(
    data = recovery_data,
    seed = 406,
    chains = 4,
    parallel_chains = min(4, parallel::detectCores()),
    iter_warmup = 250,
    iter_sampling = 250,
    adapt_delta = 0.9,
    refresh = 100,
    output_dir = recovery_dir,
    output_basename = "baseline_recovery_hmc"
  )
  saveRDS(recovery_fit, recovery_fit_path)
}

recovery_summary <- as.data.frame(recovery_fit$summary(variables = baseline_parameter_vars)) |>
  select(variable, mean, q5, q95) |>
  mutate(
    truth = c(370, 173, 62, -14, 97),
    variable = factor(variable, levels = rev(baseline_parameter_vars))
  )

save_plot(
  ggplot(recovery_summary, aes(y = variable)) +
    geom_errorbarh(aes(xmin = q5, xmax = q95), height = 0.18, color = "#2563eb") +
    geom_point(aes(x = mean), color = "#2563eb", size = 2) +
    geom_point(aes(x = truth), color = "#dc2626", size = 2.5, shape = 4, stroke = 1.2) +
    labs(
      title = "Baseline Synthetic Recovery Check",
      subtitle = "Blue: posterior mean and 90% interval, red x: true value",
      x = "Value",
      y = NULL
    ) +
    theme_minimal(),
  file.path(figures_dir, "baseline_recovery.png"),
  8,
  4.8
)
