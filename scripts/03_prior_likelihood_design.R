source("scripts/utils.R")

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

data_path <- "datasets/processed/openpowerlifting_2025_sbd.rds"
eda_path <- "datasets/intermediate/eda_group_summaries.csv"
tables_dir <- "datasets/intermediate"
figures_dir <- "figures/model_design"
sex_colors <- c(F = "#f97373", M = "#22c55e")

ensure_dir(tables_dir)
ensure_dir(figures_dir)

skewness <- function(x) {
  mean((x - mean(x))^3) / sd(x)^3
}

simulate_direct_prior <- function(data, n_draws, priors, nonlinear = FALSE, student_t = FALSE) {
  idx <- sample(nrow(data), n_draws, replace = TRUE)
  rows <- data[idx, ]

  mu <- rnorm(n_draws, priors$alpha_mean, priors$alpha_sd) +
    rnorm(n_draws, 0, priors$sex_sd) * rows$male +
    rnorm(n_draws, 0, priors$bw_sd) * rows$bw_z +
    rnorm(n_draws, 0, priors$age_sd) * rows$age_z

  if (nonlinear) {
    mu <- mu +
      rnorm(n_draws, 0, priors$quad_sd) * rows$bw_z^2 +
      rnorm(n_draws, 0, priors$quad_sd) * rows$age_z^2
  }

  sigma <- rexp(n_draws, rate = priors$sigma_rate)

  if (student_t) {
    nu <- priors$nu_offset + rexp(n_draws, rate = priors$nu_rate)
    y <- mu + sigma * rt(n_draws, df = nu)
  } else {
    y <- rnorm(n_draws, mu, sigma)
  }

  y
}

prior_summary <- function(x, label) {
  data.frame(
    candidate = label,
    prop_negative = mean(x < 0),
    p01 = unname(quantile(x, 0.01)),
    p05 = unname(quantile(x, 0.05)),
    median = median(x),
    p95 = unname(quantile(x, 0.95)),
    p99 = unname(quantile(x, 0.99)),
    row.names = NULL
  )
}

curve_draws <- function(grid, priors, n_curves = 80) {
  curves <- vector("list", n_curves)

  for (i in seq_len(n_curves)) {
    mu <- rnorm(1, priors$alpha_mean, priors$alpha_sd) +
      rnorm(1, 0, priors$sex_sd) * grid$male +
      rnorm(1, 0, priors$bw_sd) * grid$bw_z +
      rnorm(1, 0, priors$age_sd) * grid$age_z +
      rnorm(1, 0, priors$quad_sd) * grid$bw_z^2 +
      rnorm(1, 0, priors$quad_sd) * grid$age_z^2

    curves[[i]] <- data.frame(
      curve_id = i,
      Sex = grid$Sex,
      Age = grid$Age,
      BodyweightKg = grid$BodyweightKg,
      prior_mean = mu
    )
  }

  bind_rows(curves)
}

raw_df <- readRDS(data_path)
eda_groups <- read.csv(eda_path, stringsAsFactors = FALSE)

main_df <- raw_df |>
  filter(Sex %in% c("F", "M")) |>
  mutate(
    male = ifelse(Sex == "M", 1, 0),
    bw_z = as.numeric(scale(BodyweightKg)),
    age_z = as.numeric(scale(Age))
  )

bw_center <- mean(main_df$BodyweightKg)
bw_scale <- sd(main_df$BodyweightKg)
age_center <- mean(main_df$Age)
age_scale <- sd(main_df$Age)
total_mean <- mean(main_df$TotalKg)
total_sd <- sd(main_df$TotalKg)

bodyweight_spread <- subset(eda_groups, summary_type == "bodyweight_spread")
age_class_summary <- subset(eda_groups, summary_type == "age_class_summary")

bodyweight_sd_ratio <- bodyweight_spread |>
  group_by(group_1) |>
  summarise(sd_ratio = max(sd_total) / min(sd_total), .groups = "drop")

age_peaks <- age_class_summary |>
  group_by(group_1) |>
  slice_max(mean_total, n = 1, with_ties = FALSE) |>
  ungroup()

total_skew <- skewness(main_df$TotalKg)
log_total_skew <- skewness(log(main_df$TotalKg))

linear_fit <- lm(TotalKg ~ Sex + bw_z + age_z, data = main_df)
quadratic_fit <- lm(TotalKg ~ Sex + bw_z + I(bw_z^2) + age_z + I(age_z^2), data = main_df)
quad_comparison <- anova(linear_fit, quadratic_fit)

priors <- list(
  alpha_mean = total_mean,
  alpha_sd = 0.75 * total_sd,
  sex_sd = 0.60 * total_sd,
  bw_sd = 0.35 * total_sd,
  age_sd = 0.20 * total_sd,
  quad_sd = 0.10 * total_sd,
  sigma_rate = 1 / (0.45 * total_sd),
  nu_offset = 4,
  nu_rate = 1 / 6
)

# Final model specification chosen in Step 3:
# - Standardized predictors:
#     bw_z_i = (BodyweightKg_i - 83.58) / 20.86
#     age_z_i = (Age_i - 28.71) / 11.93
#     male_i = 1 for M, 0 for F
#
# - Baseline model:
#     TotalKg_i | mu_i, sigma ~ Normal(mu_i, sigma)
#     mu_i = alpha + beta_male * male_i + beta_bw * bw_z_i + beta_age * age_z_i
#
# - Main model:
#     TotalKg_i | mu_i, sigma, nu ~ Student_t(nu, mu_i, sigma)
#     mu_i = alpha + beta_male * male_i + beta_bw * bw_z_i + beta_bw2 * bw_z_i^2 +
#            beta_age * age_z_i + beta_age2 * age_z_i^2
#
# - Final priors:
#     alpha ~ Normal(483.2, 117.7)
#     beta_male ~ Normal(0, 94.2)
#     beta_bw ~ Normal(0, 54.9)
#     beta_age ~ Normal(0, 31.4)
#     beta_bw2 ~ Normal(0, 15.7)
#     beta_age2 ~ Normal(0, 15.7)
#     sigma ~ Exponential(rate = 0.0142)
#     nu = 4 + Exponential(rate = 0.1667)

set.seed(405)
prior_draws <- bind_rows(
  data.frame(candidate = "Direct Gaussian, linear mean", y = simulate_direct_prior(main_df, 5000, priors, nonlinear = FALSE, student_t = FALSE)),
  data.frame(candidate = "Direct Gaussian, nonlinear mean", y = simulate_direct_prior(main_df, 5000, priors, nonlinear = TRUE, student_t = FALSE)),
  data.frame(candidate = "Direct Student-t, nonlinear mean", y = simulate_direct_prior(main_df, 5000, priors, nonlinear = TRUE, student_t = TRUE)),
  data.frame(candidate = "Observed data", y = sample(main_df$TotalKg, 5000, replace = TRUE))
)

prior_predictive_summary <- bind_rows(
  prior_summary(subset(prior_draws, candidate == "Direct Gaussian, linear mean")$y, "Direct Gaussian, linear mean"),
  prior_summary(subset(prior_draws, candidate == "Direct Gaussian, nonlinear mean")$y, "Direct Gaussian, nonlinear mean"),
  prior_summary(subset(prior_draws, candidate == "Direct Student-t, nonlinear mean")$y, "Direct Student-t, nonlinear mean"),
  prior_summary(sample(main_df$TotalKg, 5000, replace = TRUE), "Observed data")
)

decisions <- data.frame(
  item = c(
    "analysis_scope",
    "tested_scope",
    "response_scale",
    "likelihood",
    "bodyweight_term",
    "age_term",
    "predictor_standardization",
    "baseline_model",
    "main_model",
    "intercept_prior",
    "sex_prior",
    "bodyweight_linear_prior",
    "age_linear_prior",
    "quadratic_prior",
    "sigma_prior",
    "nu_prior"
  ),
  decision = c(
    "Main models use Sex in {F, M} only",
    "Keep Tested exploratory only",
    "Model TotalKg directly on the original kg scale",
    "Use Student-t likelihood for the main model",
    "Use a nonlinear bodyweight term",
    "Use a nonlinear age term",
    "Center and scale BodyweightKg and Age before modeling",
    "TotalKg ~ Sex + bw_z + age_z",
    "TotalKg ~ Sex + bw_z + bw_z^2 + age_z + age_z^2",
    sprintf("alpha ~ normal(%.1f, %.1f)", priors$alpha_mean, priors$alpha_sd),
    sprintf("beta_male ~ normal(0, %.1f)", priors$sex_sd),
    sprintf("beta_bw ~ normal(0, %.1f)", priors$bw_sd),
    sprintf("beta_age ~ normal(0, %.1f)", priors$age_sd),
    sprintf("beta_quadratic ~ normal(0, %.1f)", priors$quad_sd),
    sprintf("sigma ~ exponential(rate = %.4f)", priors$sigma_rate),
    sprintf("nu = 4 + exponential(rate = %.4f)", priors$nu_rate)
  ),
  evidence = c(
    sprintf("Mx has %d rows out of %d in the frozen subset.", sum(raw_df$Sex == 'Mx'), nrow(raw_df)),
    "Tested is encoded as Yes or blank and is not central to the project question.",
    sprintf("TotalKg skew = %.3f while log(TotalKg) skew = %.3f, so the direct scale is already more symmetric.", total_skew, log_total_skew),
    sprintf("Spread rises across bodyweight bins: SD ratio is %.2f for F and %.2f for M.", bodyweight_sd_ratio$sd_ratio[bodyweight_sd_ratio$group_1 == 'F'], bodyweight_sd_ratio$sd_ratio[bodyweight_sd_ratio$group_1 == 'M']),
    sprintf("The quadratic upgrade over a linear mean model is very strong (F = %.1f, p < 2e-16).", quad_comparison$F[2]),
    paste0("Peak mean total by age class is ", paste(paste(age_peaks$group_1, age_peaks$group_2), collapse = '; '), "."),
    sprintf("Use bw mean/sd = %.2f / %.2f and age mean/sd = %.2f / %.2f.", bw_center, bw_scale, age_center, age_scale),
    "Baseline stays simple and interpretable for HMC vs VI comparison.",
    "Main model adds the nonlinear terms suggested by the EDA without jumping to splines yet.",
    "Keeps prior mass near the observed outcome scale while remaining weakly informative.",
    "Allows a large male-female shift without forcing it.",
    "Allows meaningful bodyweight effects on the standardized scale.",
    "Age effects were visibly weaker than bodyweight effects in the EDA, so the prior is tighter.",
    "Keeps curvature plausible but avoids wild prior mean shapes.",
    "Prior mean for sigma is about 45% of the observed SD, which keeps the prior predictive spread broad without too much mass below zero.",
    "Centers the t degrees of freedom in a moderate heavy-tail range."
  ),
  row.names = NULL
)

write.csv(decisions, file.path(tables_dir, "model_design_decisions.csv"), row.names = FALSE)
write.csv(prior_predictive_summary, file.path(tables_dir, "prior_predictive_summary.csv"), row.names = FALSE)

save_plot(
  ggplot(prior_draws, aes(x = y, color = candidate)) +
    geom_density(linewidth = 1) +
    coord_cartesian(xlim = c(0, 1200)) +
    labs(
      title = "Prior Predictive Comparison",
      subtitle = "Candidate forms against the observed scale of TotalKg",
      x = "Total (kg)",
      y = "Density"
    ) +
    theme_minimal(),
  file.path(figures_dir, "prior_predictive_candidates.png"),
  8,
  5
)

curve_grid <- expand.grid(
  BodyweightKg = seq(40, 140, by = 2.5),
  Age = c(25, 40, 60),
  Sex = c("F", "M"),
  KEEP.OUT.ATTRS = FALSE
)
curve_grid$male <- ifelse(curve_grid$Sex == "M", 1, 0)
curve_grid$bw_z <- (curve_grid$BodyweightKg - bw_center) / bw_scale
curve_grid$age_z <- (curve_grid$Age - age_center) / age_scale

curve_df <- curve_draws(curve_grid, priors)

save_plot(
  ggplot(curve_df, aes(x = BodyweightKg, y = prior_mean, group = interaction(curve_id, Sex), color = Sex)) +
    geom_line(alpha = 0.12) +
    facet_wrap(~Age) +
    scale_color_manual(values = sex_colors) +
    labs(
      title = "Prior Mean Curves for the Chosen Main Model",
      subtitle = "Nonlinear mean structure before fitting any Bayesian model",
      x = "Bodyweight (kg)",
      y = "Prior mean total (kg)"
    ) +
    theme_minimal(),
  file.path(figures_dir, "prior_mean_curves.png"),
  8,
  5
)

print(decisions[, c("item", "decision")])
print(prior_predictive_summary)

# Step-3 design notes:
# - TotalKg is modeled directly because the raw kg scale is already reasonably symmetric after filtering to F/M.
# - Student-t is preferred over Gaussian because the EDA showed larger spread at heavier bodyweights and some extreme totals.
# - Both BodyweightKg and Age will be centered and scaled later, and both get quadratic terms in the richer model.
# - The baseline model remains linear in bw_z and age_z so we still have a simple benchmark for HMC vs VI.
