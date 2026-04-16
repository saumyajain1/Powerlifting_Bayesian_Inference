source("scripts/utils.R")

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(posterior))
suppressPackageStartupMessages(library(tidyr))

tables_dir <- "datasets/intermediate"
figures_dir <- "figures/main_checks"
artifacts_dir <- "artifacts/main"

ensure_dir(figures_dir)

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

main_mu_draws <- function(draws, male, bw_z, age_z) {
  draws[, "alpha"] +
    draws[, "beta_male"] * male +
    draws[, "beta_bw"] * bw_z +
    draws[, "beta_bw2"] * bw_z^2 +
    draws[, "beta_age"] * age_z +
    draws[, "beta_age2"] * age_z^2
}

main_hmc <- readRDS(file.path(artifacts_dir, "main_hmc_fit.rds"))
main_mf <- readRDS(file.path(artifacts_dir, "main_vi_meanfield_fit.rds"))
main_fr <- readRDS(file.path(artifacts_dir, "main_vi_fullrank_fit.rds"))
model_df <- readRDS(file.path(tables_dir, "model_data.rds"))
prediction_summary <- read.csv(file.path(tables_dir, "main_prediction_summary.csv"), stringsAsFactors = FALSE)
parameter_summary <- read.csv(file.path(tables_dir, "main_parameter_summary.csv"), stringsAsFactors = FALSE)
ppc_df <- baseline_ppc_subset(model_df)

trace_df <- as.data.frame(as_draws_df(main_hmc$draws(variables = main_parameter_vars))) |>
  select(all_of(main_parameter_vars), .chain, .iteration) |>
  pivot_longer(cols = all_of(main_parameter_vars), names_to = "parameter", values_to = "value")

save_plot(
  ggplot(trace_df, aes(.iteration, value, color = factor(.chain), group = interaction(.chain, parameter))) +
    geom_line(alpha = 0.65, linewidth = 0.3) +
    facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
    labs(title = "Main HMC Trace Plots", x = "Iteration", y = "Value", color = "Chain") +
    theme_minimal(),
  file.path(figures_dir, "main_traceplot.png"),
  9,
  8
)

param_plot_df <- parameter_summary |>
  select(variable, mean, q5, q95, method) |>
  mutate(
    variable = factor(variable, levels = main_parameter_vars),
    method = factor(method, levels = c("HMC", "Mean-field VI", "Full-rank VI"))
  )

save_plot(
  ggplot(param_plot_df, aes(mean, variable, color = method)) +
    geom_errorbarh(aes(xmin = q5, xmax = q95), height = 0.18, position = position_dodge(width = 0.55)) +
    geom_point(position = position_dodge(width = 0.55), size = 2) +
    labs(
      title = "Main Model Posterior Comparison",
      subtitle = "HMC against mean-field and full-rank VI",
      x = "Estimate",
      y = NULL,
      color = "Method"
    ) +
    theme_minimal(),
  file.path(figures_dir, "main_hmc_vs_vi.png"),
  8.5,
  5.6
)

hmc_y_ppc <- as_draws_matrix(main_hmc$draws(variables = "y_ppc"))
mf_y_ppc <- as_draws_matrix(main_mf$draws(variables = "y_ppc"))
fr_y_ppc <- as_draws_matrix(main_fr$draws(variables = "y_ppc"))

set.seed(405)
save_plot(
  ggplot(
    bind_rows(
      ppc_density_df(hmc_y_ppc, ppc_df$TotalKg, "HMC"),
      ppc_density_df(mf_y_ppc, ppc_df$TotalKg, "Mean-field VI"),
      ppc_density_df(fr_y_ppc, ppc_df$TotalKg, "Full-rank VI")
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
      title = "Main Posterior Predictive Check",
      subtitle = "Observed subset density against posterior predictive replicated densities",
      x = "Total (kg)",
      y = "Density"
    ) +
    theme_minimal(),
  file.path(figures_dir, "main_ppc_density.png"),
  8,
  8.5
)

coverage_df <- bind_rows(
  ppc_coverage_df(hmc_y_ppc, ppc_df$TotalKg, "HMC"),
  ppc_coverage_df(mf_y_ppc, ppc_df$TotalKg, "Mean-field VI"),
  ppc_coverage_df(fr_y_ppc, ppc_df$TotalKg, "Full-rank VI")
)

save_plot(
  ggplot(coverage_df, aes(nominal, actual, color = method)) +
    geom_abline(slope = 1, intercept = 0, color = "#9ca3af", linetype = "dashed") +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = "Main PPC Coverage",
      subtitle = "Actual central-interval coverage against nominal coverage",
      x = "Nominal coverage",
      y = "Actual coverage",
      color = "Method"
    ) +
    theme_minimal(),
  file.path(figures_dir, "main_ppc_coverage.png"),
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
      title = "Main Predicted Strength Curves",
      subtitle = "Posterior mean and 90% intervals across bodyweight",
      x = "Bodyweight (kg)",
      y = "Predicted Total (kg)"
    ) +
    theme_minimal(),
  file.path(figures_dir, "main_prediction_curves.png"),
  10,
  7
)

main_draws <- as_draws_matrix(main_hmc$draws(variables = c("alpha", "beta_male", "beta_bw", "beta_bw2", "beta_age", "beta_age2")))
bw_center <- mean(model_df$BodyweightKg)
bw_scale <- sd(model_df$BodyweightKg)
age_center <- mean(model_df$Age)
age_scale <- sd(model_df$Age)

weight_class_summary <- model_df |>
  filter(WeightClassKg != "") |>
  mutate(
    class_label = trimws(WeightClassKg),
    plus_class = grepl("\\+", class_label),
    upper = suppressWarnings(as.numeric(sub("\\+.*$", "", class_label)))
  ) |>
  filter(!plus_class, !is.na(upper)) |>
  group_by(Sex, class_label, upper) |>
  summarise(n = n(), bodyweight = median(BodyweightKg), .groups = "drop") |>
  filter(n >= 100) |>
  arrange(Sex, upper) |>
  group_by(Sex) |>
  mutate(
    next_class = lead(class_label),
    next_bodyweight = lead(bodyweight)
  ) |>
  filter(!is.na(next_class)) |>
  ungroup()

age_values <- c(25, 40, 60)
delta_rows <- vector("list", nrow(weight_class_summary) * length(age_values))
row_id <- 1

for (i in seq_len(nrow(weight_class_summary))) {
  row <- weight_class_summary[i, ]
  male <- ifelse(row$Sex == "M", 1, 0)
  bw1_z <- (row$bodyweight - bw_center) / bw_scale
  bw2_z <- (row$next_bodyweight - bw_center) / bw_scale

  for (age in age_values) {
    age_z <- (age - age_center) / age_scale
    delta <- main_mu_draws(main_draws, male, bw2_z, age_z) - main_mu_draws(main_draws, male, bw1_z, age_z)

    delta_rows[[row_id]] <- data.frame(
      Sex = row$Sex,
      Age = age,
      from_class = row$class_label,
      to_class = row$next_class,
      mean = mean(delta),
      q05 = unname(quantile(delta, 0.05)),
      q95 = unname(quantile(delta, 0.95)),
      row.names = NULL
    )
    row_id <- row_id + 1
  }
}

weight_class_delta_summary <- bind_rows(delta_rows)
write.csv(weight_class_delta_summary, file.path(tables_dir, "weight_class_delta_summary.csv"), row.names = FALSE)

save_plot(
  ggplot(
    transform(weight_class_delta_summary, transition = factor(paste(from_class, "\u2192", to_class), levels = unique(paste(from_class, "\u2192", to_class)))),
    aes(transition, mean, color = Sex)
  ) +
    geom_hline(yintercept = 0, color = "#d1d5db") +
    geom_errorbar(aes(ymin = q05, ymax = q95), width = 0.15, position = position_dodge(width = 0.45)) +
    geom_point(position = position_dodge(width = 0.45), size = 2) +
    facet_wrap(~ Age, ncol = 1) +
    labs(
      title = "Estimated Gain From Moving Up One Weight Class",
      subtitle = "Main HMC posterior mean and 90% intervals",
      x = "Adjacent weight-class transition",
      y = "Predicted change in total (kg)",
      color = "Sex"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  file.path(figures_dir, "weight_class_deltas.png"),
  10,
  8
)

fixed_delta_grid <- data.frame(
  BodyweightKg = seq(
    unname(quantile(model_df$BodyweightKg, 0.02)),
    unname(quantile(model_df$BodyweightKg, 0.98)) - 5,
    length.out = 120
  )
)

fixed_delta_summary <- bind_rows(lapply(seq_len(nrow(fixed_delta_grid)), function(i) {
  bw1 <- fixed_delta_grid$BodyweightKg[i]
  bw2 <- bw1 + 5
  bw1_z <- (bw1 - bw_center) / bw_scale
  bw2_z <- (bw2 - bw_center) / bw_scale
  delta <- main_mu_draws(main_draws, 0, bw2_z, 0) - main_mu_draws(main_draws, 0, bw1_z, 0)

  data.frame(
    BodyweightKg = bw1,
    mean = mean(delta),
    q05 = unname(quantile(delta, 0.05)),
    q95 = unname(quantile(delta, 0.95)),
    row.names = NULL
  )
}))

save_plot(
  ggplot(fixed_delta_summary, aes(BodyweightKg, mean)) +
    geom_line(color = "#0891b2", linewidth = 1) +
    geom_line(aes(y = q05), color = "#67e8f9", linewidth = 0.8, linetype = "dashed") +
    geom_line(aes(y = q95), color = "#67e8f9", linewidth = 0.8, linetype = "dashed") +
    labs(
      title = "Estimated Gain From a +5 kg Bodyweight Increase",
      subtitle = "Main HMC posterior mean with 90% interval bounds; same for all ages and sexes in the current model",
      x = "Starting bodyweight (kg)",
      y = "Predicted change in total (kg)"
    ) +
    theme_minimal(),
  file.path(figures_dir, "fixed_bodyweight_delta_5kg.png"),
  8.5,
  5
)
