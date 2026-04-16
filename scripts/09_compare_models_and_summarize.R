source("scripts/utils.R")

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(posterior))
suppressPackageStartupMessages(library(tidyr))

tables_dir <- "datasets/intermediate"
figures_dir <- "figures/model_comparison"

ensure_dir(tables_dir)
ensure_dir(figures_dir)

ppc_density_df <- function(draw_matrix, observed, model, n_reps = 25) {
  keep <- sample(seq_len(nrow(draw_matrix)), n_reps)

  rbind(
    bind_rows(lapply(seq_along(keep), function(i) {
      dens <- density(draw_matrix[keep[i], ])
      data.frame(x = dens$x, y = dens$y, model = model, source = "Posterior predictive", rep = i)
    })),
    {
      dens <- density(observed)
      data.frame(x = dens$x, y = dens$y, model = model, source = "Observed data", rep = 0)
    }
  )
}

ppc_coverage_df <- function(draw_matrix, observed, model, levels = seq(0.1, 0.95, by = 0.05)) {
  bind_rows(lapply(levels, function(level) {
    alpha <- (1 - level) / 2
    lower <- apply(draw_matrix, 2, quantile, probs = alpha)
    upper <- apply(draw_matrix, 2, quantile, probs = 1 - alpha)

    data.frame(model = model, nominal = level, actual = mean(observed >= lower & observed <= upper))
  }))
}

ppc_metrics <- function(draw_matrix, observed, model, method) {
  lower <- apply(draw_matrix, 2, quantile, probs = 0.05)
  upper <- apply(draw_matrix, 2, quantile, probs = 0.95)
  median_pred <- apply(draw_matrix, 2, median)

  data.frame(
    model = model,
    method = method,
    mae = mean(abs(observed - median_pred)),
    coverage90 = mean(observed >= lower & observed <= upper),
    mean_interval90 = mean(upper - lower),
    row.names = NULL
  )
}

baseline_mu_draws <- function(draws, male, bw_z, age_z) {
  draws[, "alpha"] +
    draws[, "beta_male"] * male +
    draws[, "beta_bw"] * bw_z +
    draws[, "beta_age"] * age_z
}

main_mu_draws <- function(draws, male, bw_z, age_z) {
  draws[, "alpha"] +
    draws[, "beta_male"] * male +
    draws[, "beta_bw"] * bw_z +
    draws[, "beta_bw2"] * bw_z^2 +
    draws[, "beta_age"] * age_z +
    draws[, "beta_age2"] * age_z^2
}

model_df <- readRDS(file.path(tables_dir, "model_data.rds"))
prediction_grid <- readRDS(file.path(tables_dir, "prediction_grid.rds"))
ppc_df <- baseline_ppc_subset(model_df)

baseline_prediction <- read.csv(file.path(tables_dir, "baseline_prediction_summary.csv"), stringsAsFactors = FALSE)
main_prediction <- read.csv(file.path(tables_dir, "main_prediction_summary.csv"), stringsAsFactors = FALSE)

baseline_hmc <- readRDS("artifacts/baseline/baseline_hmc_fit.rds")
baseline_mf <- readRDS("artifacts/baseline/baseline_vi_meanfield_fit.rds")
baseline_fr <- readRDS("artifacts/baseline/baseline_vi_fullrank_fit.rds")
main_hmc <- readRDS("artifacts/main/main_hmc_fit.rds")
main_mf <- readRDS("artifacts/main/main_vi_meanfield_fit.rds")
main_fr <- readRDS("artifacts/main/main_vi_fullrank_fit.rds")

baseline_hmc_y_ppc <- as_draws_matrix(baseline_hmc$draws(variables = "y_ppc"))
baseline_mf_y_ppc <- as_draws_matrix(baseline_mf$draws(variables = "y_ppc"))
baseline_fr_y_ppc <- as_draws_matrix(baseline_fr$draws(variables = "y_ppc"))
main_hmc_y_ppc <- as_draws_matrix(main_hmc$draws(variables = "y_ppc"))
main_mf_y_ppc <- as_draws_matrix(main_mf$draws(variables = "y_ppc"))
main_fr_y_ppc <- as_draws_matrix(main_fr$draws(variables = "y_ppc"))

model_comparison_summary <- bind_rows(
  ppc_metrics(baseline_hmc_y_ppc, ppc_df$TotalKg, "Baseline", "HMC"),
  ppc_metrics(baseline_mf_y_ppc, ppc_df$TotalKg, "Baseline", "Mean-field VI"),
  ppc_metrics(baseline_fr_y_ppc, ppc_df$TotalKg, "Baseline", "Full-rank VI"),
  ppc_metrics(main_hmc_y_ppc, ppc_df$TotalKg, "Main", "HMC"),
  ppc_metrics(main_mf_y_ppc, ppc_df$TotalKg, "Main", "Mean-field VI"),
  ppc_metrics(main_fr_y_ppc, ppc_df$TotalKg, "Main", "Full-rank VI")
)

write.csv(model_comparison_summary, file.path(tables_dir, "model_comparison_summary.csv"), row.names = FALSE)

metrics_long <- model_comparison_summary |>
  pivot_longer(cols = c(mae, coverage90, mean_interval90), names_to = "metric", values_to = "value") |>
  mutate(
    metric = factor(metric, levels = c("mae", "coverage90", "mean_interval90"), labels = c("MAE", "90% coverage", "Mean 90% interval width")),
    method = factor(method, levels = c("HMC", "Mean-field VI", "Full-rank VI")),
    model = factor(model, levels = c("Baseline", "Main"))
  )

save_plot(
  ggplot(metrics_long, aes(method, value, color = model, group = model)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    labs(
      title = "Predictive Comparison Across Models and Methods",
      subtitle = "Lower MAE and narrower, well-calibrated intervals are better",
      x = NULL,
      y = NULL,
      color = "Model"
    ) +
    theme_minimal(),
  file.path(figures_dir, "model_comparison_metrics.png"),
  11,
  4.8
)

set.seed(405)
save_plot(
  ggplot(
    bind_rows(
      ppc_density_df(baseline_hmc_y_ppc, ppc_df$TotalKg, "Baseline"),
      ppc_density_df(main_hmc_y_ppc, ppc_df$TotalKg, "Main")
    )
  ) +
    geom_line(
      data = function(x) subset(x, source == "Posterior predictive"),
      aes(x, y, group = interaction(model, rep)),
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
    facet_wrap(~ model, ncol = 1) +
    labs(
      title = "Baseline vs Main HMC Posterior Predictive Check",
      subtitle = "Observed subset density against posterior predictive replicated densities",
      x = "Total (kg)",
      y = "Density"
    ) +
    theme_minimal(),
  file.path(figures_dir, "baseline_vs_main_ppc_density.png"),
  8,
  6.5
)

save_plot(
  ggplot(
    bind_rows(
      ppc_coverage_df(baseline_hmc_y_ppc, ppc_df$TotalKg, "Baseline"),
      ppc_coverage_df(main_hmc_y_ppc, ppc_df$TotalKg, "Main")
    ),
    aes(nominal, actual, color = model)
  ) +
    geom_abline(slope = 1, intercept = 0, color = "#9ca3af", linetype = "dashed") +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = "Baseline vs Main HMC Coverage",
      subtitle = "Actual central-interval coverage against nominal coverage",
      x = "Nominal coverage",
      y = "Actual coverage",
      color = "Model"
    ) +
    theme_minimal(),
  file.path(figures_dir, "baseline_vs_main_ppc_coverage.png"),
  6.5,
  6
)

curve_compare <- bind_rows(
  transform(subset(baseline_prediction, method == "HMC"), model = "Baseline"),
  transform(subset(main_prediction, method == "HMC"), model = "Main")
) |>
  mutate(
    Sex = factor(Sex, levels = c("F", "M")),
    model = factor(model, levels = c("Baseline", "Main"))
  )

save_plot(
  ggplot(curve_compare, aes(BodyweightKg, mean, color = model)) +
    geom_line(linewidth = 0.95) +
    facet_grid(Sex ~ Age) +
    labs(
      title = "Baseline vs Main HMC Curves",
      subtitle = "The main model adds nonlinear bodyweight and age effects",
      x = "Bodyweight (kg)",
      y = "Predicted Total (kg)",
      color = "Model"
    ) +
    theme_minimal(),
  file.path(figures_dir, "baseline_vs_main_hmc_curves.png"),
  10,
  6
)

baseline_mu_pred <- as_draws_matrix(baseline_hmc$draws(variables = "mu_pred"))
main_mu_pred <- as_draws_matrix(main_hmc$draws(variables = "mu_pred"))
n_draws <- min(nrow(baseline_mu_pred), nrow(main_mu_pred))
mu_diff <- main_mu_pred[seq_len(n_draws), , drop = FALSE] - baseline_mu_pred[seq_len(n_draws), , drop = FALSE]

difference_summary <- cbind(
  prediction_grid,
  data.frame(
    mean = colMeans(mu_diff),
    q05 = apply(mu_diff, 2, quantile, probs = 0.05),
    q95 = apply(mu_diff, 2, quantile, probs = 0.95),
    row.names = NULL
  )
)

save_plot(
  ggplot(transform(difference_summary, Sex = factor(Sex, levels = c("F", "M"))), aes(BodyweightKg, mean, color = Sex, fill = Sex)) +
    geom_hline(yintercept = 0, color = "#d1d5db") +
    geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.18, color = NA) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ Age, nrow = 1) +
    labs(
      title = "Main Minus Baseline HMC Prediction",
      subtitle = "Difference in posterior mean curves",
      x = "Bodyweight (kg)",
      y = "Predicted total difference (kg)"
    ) +
    theme_minimal(),
  file.path(figures_dir, "baseline_vs_main_difference_curves.png"),
  10,
  4.8
)

bw_center <- mean(model_df$BodyweightKg)
bw_scale <- sd(model_df$BodyweightKg)
fixed_delta_grid <- data.frame(
  BodyweightKg = seq(
    unname(quantile(model_df$BodyweightKg, 0.02)),
    unname(quantile(model_df$BodyweightKg, 0.98)) - 5,
    length.out = 120
  )
)

baseline_draws <- as_draws_matrix(baseline_hmc$draws(variables = c("alpha", "beta_male", "beta_bw", "beta_age")))
main_draws <- as_draws_matrix(main_hmc$draws(variables = c("alpha", "beta_male", "beta_bw", "beta_bw2", "beta_age", "beta_age2")))

fixed_delta_compare <- bind_rows(lapply(seq_len(nrow(fixed_delta_grid)), function(i) {
  bw1 <- fixed_delta_grid$BodyweightKg[i]
  bw2 <- bw1 + 5
  bw1_z <- (bw1 - bw_center) / bw_scale
  bw2_z <- (bw2 - bw_center) / bw_scale

  baseline_delta <- baseline_mu_draws(baseline_draws, 0, bw2_z, 0) - baseline_mu_draws(baseline_draws, 0, bw1_z, 0)
  main_delta <- main_mu_draws(main_draws, 0, bw2_z, 0) - main_mu_draws(main_draws, 0, bw1_z, 0)

  bind_rows(
    data.frame(
      model = "Baseline",
      BodyweightKg = bw1,
      mean = mean(baseline_delta),
      q05 = unname(quantile(baseline_delta, 0.05)),
      q95 = unname(quantile(baseline_delta, 0.95)),
      row.names = NULL
    ),
    data.frame(
      model = "Main",
      BodyweightKg = bw1,
      mean = mean(main_delta),
      q05 = unname(quantile(main_delta, 0.05)),
      q95 = unname(quantile(main_delta, 0.95)),
      row.names = NULL
    )
  )
}))

save_plot(
  ggplot(transform(fixed_delta_compare, model = factor(model, levels = c("Baseline", "Main"))), aes(BodyweightKg, mean, color = model, fill = model)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.16, color = NA) +
    geom_line(linewidth = 0.95) +
    labs(
      title = "Baseline vs Main Gain From a +5 kg Bodyweight Increase",
      subtitle = "Baseline implies a constant gain; the main model allows diminishing returns",
      x = "Starting bodyweight (kg)",
      y = "Predicted change in total (kg)",
      color = "Model",
      fill = "Model"
    ) +
    theme_minimal(),
  file.path(figures_dir, "baseline_vs_main_fixed_delta_5kg.png"),
  8.5,
  5
)
