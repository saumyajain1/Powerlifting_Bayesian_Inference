# Build the exact model-ready data objects used by the Stan scripts.
#
# This step keeps only the final F/M analysis population, creates the
# standardized predictors chosen in Step 3, and saves reusable inputs for
# later fitting and plotting.

source("scripts/utils.R")

suppressPackageStartupMessages(library(dplyr))

data_path <- "datasets/processed/openpowerlifting_2025_sbd.rds"
output_dir <- "datasets/intermediate"

model_csv <- file.path(output_dir, "model_data.csv")
model_rds <- file.path(output_dir, "model_data.rds")
grid_csv <- file.path(output_dir, "prediction_grid.csv")
grid_rds <- file.path(output_dir, "prediction_grid.rds")
baseline_rds <- file.path(output_dir, "stan_data_baseline.rds")
main_rds <- file.path(output_dir, "stan_data_main.rds")

ensure_dir(output_dir)

raw_df <- readRDS(data_path)

model_df <- raw_df |>
  filter(Sex %in% c("F", "M")) |>
  mutate(
    male = ifelse(Sex == "M", 1L, 0L),
    bw_z = as.numeric(scale(BodyweightKg)),
    age_z = as.numeric(scale(Age))
  ) |>
  mutate(
    bw_z2 = bw_z^2,
    age_z2 = age_z^2
  )

scaling <- list(
  bw_center = mean(model_df$BodyweightKg),
  bw_scale = sd(model_df$BodyweightKg),
  age_center = mean(model_df$Age),
  age_scale = sd(model_df$Age)
)

prediction_grid <- bind_rows(lapply(c("F", "M"), function(sex_value) {
  sex_df <- filter(model_df, Sex == sex_value)
  bw_seq <- seq(
    unname(quantile(sex_df$BodyweightKg, 0.01)),
    unname(quantile(sex_df$BodyweightKg, 0.99)),
    length.out = 80
  )

  expand.grid(
    BodyweightKg = bw_seq,
    Age = c(25, 40, 60),
    Sex = sex_value,
    KEEP.OUT.ATTRS = FALSE
  )
})) |>
  mutate(
    male = ifelse(Sex == "M", 1L, 0L),
    bw_z = (BodyweightKg - scaling$bw_center) / scaling$bw_scale,
    age_z = (Age - scaling$age_center) / scaling$age_scale,
    bw_z2 = bw_z^2,
    age_z2 = age_z^2
  )

attr(model_df, "scaling") <- scaling
attr(prediction_grid, "scaling") <- scaling

stan_data_baseline <- list(
  N = nrow(model_df),
  y = model_df$TotalKg,
  male = model_df$male,
  bw_z = model_df$bw_z,
  age_z = model_df$age_z
)

stan_data_main <- list(
  N = nrow(model_df),
  y = model_df$TotalKg,
  male = model_df$male,
  bw_z = model_df$bw_z,
  bw_z2 = model_df$bw_z2,
  age_z = model_df$age_z,
  age_z2 = model_df$age_z2
)

save_csv_rds(model_df, model_csv, model_rds)
save_csv_rds(prediction_grid, grid_csv, grid_rds)
saveRDS(stan_data_baseline, baseline_rds)
saveRDS(stan_data_main, main_rds)
