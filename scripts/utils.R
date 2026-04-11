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
