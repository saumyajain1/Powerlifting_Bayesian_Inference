ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

save_csv_rds <- function(data, csv_path, rds_path) {
  ensure_dir(dirname(csv_path))
  ensure_dir(dirname(rds_path))
  write.csv(data, csv_path, row.names = FALSE)
  saveRDS(data, rds_path)
}
