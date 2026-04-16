ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

save_csv_rds <- function(data, csv_path, rds_path) {
  ensure_dir(dirname(csv_path))
  ensure_dir(dirname(rds_path))
  write.csv(data, csv_path, row.names = FALSE)
  saveRDS(data, rds_path)
}

save_plot <- function(plot, path, width, height) {
  ensure_dir(dirname(path))
  ggplot2::ggsave(path, plot = plot, width = width, height = height, dpi = 300, bg = "white")
}

summary_block <- function(section, item, value) {
  data.frame(
    section = section,
    item = item,
    value = as.character(value),
    row.names = NULL
  )
}

count_summary <- function(x, name) {
  out <- as.data.frame(table(x), stringsAsFactors = FALSE)
  names(out) <- c(name, "n")
  out$prop <- out$n / sum(out$n)
  out
}

numeric_summary_table <- function(data, vars) {
  do.call(rbind, lapply(vars, function(var) {
    x <- data[[var]]
    data.frame(
      variable = var,
      mean = mean(x),
      sd = sd(x),
      min = min(x),
      p05 = unname(quantile(x, 0.05)),
      median = median(x),
      p95 = unname(quantile(x, 0.95)),
      max = max(x),
      row.names = NULL
    )
  }))
}

numeric_summary_rows <- function(summary_table) {
  stats <- c("mean", "sd", "min", "p05", "median", "p95", "max")

  do.call(rbind, lapply(seq_len(nrow(summary_table)), function(i) {
    summary_block(
      "numeric_summary",
      paste(summary_table$variable[i], stats, sep = "_"),
      unlist(summary_table[i, stats])
    )
  }))
}

quantile_bin <- function(x, probs = seq(0, 1, by = 0.1)) {
  cut(x, breaks = unique(quantile(x, probs, na.rm = TRUE)), include.lowest = TRUE)
}

baseline_parameter_vars <- c("alpha", "beta_male", "beta_bw", "beta_age", "sigma")
main_parameter_vars <- c("alpha", "beta_male", "beta_bw", "beta_bw2", "beta_age", "beta_age2", "sigma", "nu")

baseline_prior_data <- function() {
  list(
    alpha_mean = 483.2,
    alpha_sd = 117.7,
    beta_male_sd = 94.2,
    beta_bw_sd = 54.9,
    beta_age_sd = 31.4,
    sigma_rate = 0.0142
  )
}

baseline_ppc_subset <- function(model_df, n = 300, seed = 405) {
  set.seed(seed)
  idx <- sort(sample(seq_len(nrow(model_df)), n))
  out <- model_df[idx, c("Sex", "Age", "BodyweightKg", "TotalKg", "male", "bw_z", "bw_z2", "age_z", "age_z2")]
  out$ppc_id <- seq_len(nrow(out))
  out
}

baseline_stan_data <- function(baseline_data, prediction_grid, ppc_df) {
  c(
    baseline_data,
    list(
      N_pred = nrow(prediction_grid),
      male_pred = prediction_grid$male,
      bw_z_pred = prediction_grid$bw_z,
      age_z_pred = prediction_grid$age_z,
      N_ppc = nrow(ppc_df),
      male_ppc = ppc_df$male,
      bw_z_ppc = ppc_df$bw_z,
      age_z_ppc = ppc_df$age_z
    ),
    baseline_prior_data()
  )
}

summarize_generated <- function(draw_matrix, prefix, meta) {
  values <- draw_matrix[, grep(paste0("^", prefix, "\\["), colnames(draw_matrix)), drop = FALSE]

  cbind(
    meta,
    data.frame(
      mean = colMeans(values),
      sd = apply(values, 2, sd),
      q05 = apply(values, 2, quantile, probs = 0.05),
      median = apply(values, 2, median),
      q95 = apply(values, 2, quantile, probs = 0.95),
      row.names = NULL
    )
  )
}

align_summary_columns <- function(summary_df) {
  keep <- c("variable", "mean", "median", "sd", "mad", "q5", "q95", "rhat", "ess_bulk", "ess_tail")

  for (name in setdiff(keep, names(summary_df))) {
    summary_df[[name]] <- NA_real_
  }

  summary_df[, keep]
}

fit_parameter_summary <- function(fit, variables, method) {
  out <- align_summary_columns(as.data.frame(fit$summary(variables = variables)))
  out$method <- method
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

ppc_density_df <- function(draw_matrix, observed, label, n_reps = 25) {
  keep <- sample(seq_len(nrow(draw_matrix)), n_reps)

  rbind(
    dplyr::bind_rows(lapply(seq_along(keep), function(i) {
      dens <- density(draw_matrix[keep[i], ])
      data.frame(x = dens$x, y = dens$y, label = label, source = "Posterior predictive", rep = i)
    })),
    {
      dens <- density(observed)
      data.frame(x = dens$x, y = dens$y, label = label, source = "Observed data", rep = 0)
    }
  )
}

ppc_coverage_df <- function(draw_matrix, observed, label, levels = seq(0.1, 0.95, by = 0.05)) {
  dplyr::bind_rows(lapply(levels, function(level) {
    alpha <- (1 - level) / 2
    lower <- apply(draw_matrix, 2, quantile, probs = alpha)
    upper <- apply(draw_matrix, 2, quantile, probs = 1 - alpha)

    data.frame(label = label, nominal = level, actual = mean(observed >= lower & observed <= upper))
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

load_or_sample_fit <- function(model, data, fit_path, output_dir, output_basename, seed, adapt_delta, iter_warmup = 500, iter_sampling = 500, refresh = 200) {
  if (file.exists(fit_path)) {
    return(readRDS(fit_path))
  }

  fit <- model$sample(
    data = data,
    seed = seed,
    chains = 4,
    parallel_chains = min(4, parallel::detectCores()),
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    refresh = refresh,
    output_dir = output_dir,
    output_basename = output_basename
  )
  saveRDS(fit, fit_path)
  fit
}

load_or_vi_fit <- function(model, data, fit_path, output_dir, output_basename, seed, algorithm, iter = 10000, output_samples = 1000, refresh = 200) {
  if (file.exists(fit_path)) {
    return(readRDS(fit_path))
  }

  fit <- model$variational(
    data = data,
    seed = seed,
    algorithm = algorithm,
    iter = iter,
    output_samples = output_samples,
    refresh = refresh,
    output_dir = output_dir,
    output_basename = output_basename
  )
  saveRDS(fit, fit_path)
  fit
}

fit_list_parameter_summary <- function(fits, variables) {
  dplyr::bind_rows(lapply(names(fits), function(method) {
    fit_parameter_summary(fits[[method]], variables, method)
  }))
}

fit_list_generated_summary <- function(fits, variable, prefix, meta) {
  dplyr::bind_rows(lapply(names(fits), function(method) {
    transform(
      summarize_generated(posterior::as_draws_matrix(fits[[method]]$draws(variables = variable)), prefix, meta),
      method = method
    )
  }))
}

main_prior_data <- function() {
  list(
    alpha_mean = 483.2,
    alpha_sd = 117.7,
    beta_male_sd = 94.2,
    beta_bw_sd = 54.9,
    beta_bw2_sd = 15.7,
    beta_age_sd = 31.4,
    beta_age2_sd = 15.7,
    sigma_rate = 0.0142,
    nu_offset = 4,
    nu_rate = 0.1667
  )
}

main_stan_data <- function(main_data, prediction_grid, ppc_df) {
  c(
    main_data,
    list(
      N_pred = nrow(prediction_grid),
      male_pred = prediction_grid$male,
      bw_z_pred = prediction_grid$bw_z,
      bw_z2_pred = prediction_grid$bw_z2,
      age_z_pred = prediction_grid$age_z,
      age_z2_pred = prediction_grid$age_z2,
      N_ppc = nrow(ppc_df),
      male_ppc = ppc_df$male,
      bw_z_ppc = ppc_df$bw_z,
      bw_z2_ppc = ppc_df$bw_z2,
      age_z_ppc = ppc_df$age_z,
      age_z2_ppc = ppc_df$age_z2
    ),
    main_prior_data()
  )
}
