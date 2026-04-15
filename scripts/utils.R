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
  out <- model_df[idx, c("Sex", "Age", "BodyweightKg", "TotalKg", "male", "bw_z", "age_z")]
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
