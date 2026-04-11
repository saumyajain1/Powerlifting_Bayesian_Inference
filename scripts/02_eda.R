source("scripts/utils.R")

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

data_path <- "datasets/processed/openpowerlifting_2025_sbd.rds"
tables_dir <- "datasets/intermediate"
figures_dir <- "figures/eda"
sex_colors <- c(F = "#f97373", M = "#22c55e")

ensure_dir(tables_dir)
ensure_dir(figures_dir)

raw_df <- readRDS(data_path)

power_df <- raw_df |>
  mutate(
    TestedGroup = ifelse(Tested == "", "Blank", Tested),
    WeightClassSimple = suppressWarnings(as.numeric(gsub("\\+", "", WeightClassKg)))
  )

main_sex_df <- power_df |>
  filter(Sex %in% c("F", "M")) |>
  mutate(
    BodyweightBin = quantile_bin(BodyweightKg),
    AgeBin = quantile_bin(Age)
  )

key_vars <- c("Sex", "Tested", "Age", "BodyweightKg", "TotalKg", "Event", "Equipment", "Place")
sex_counts <- count_summary(power_df$Sex, "Sex")
tested_counts <- count_summary(power_df$TestedGroup, "TestedGroup")
numeric_summary <- numeric_summary_table(power_df, c("Age", "BodyweightKg", "TotalKg"))

summary_rows <- bind_rows(
  summary_block(
    "overview",
    c("rows", "columns", "unique_federations", "unique_meets"),
    c(nrow(raw_df), ncol(raw_df), n_distinct(raw_df$Federation), n_distinct(raw_df$MeetName))
  ),
  summary_block("missing_core", key_vars, sapply(key_vars, function(x) sum(is.na(power_df[[x]])))),
  summary_block("sex_count", sex_counts$Sex, sex_counts$n),
  summary_block("sex_prop", sex_counts$Sex, sprintf("%.6f", sex_counts$prop)),
  summary_block("tested_count", tested_counts$TestedGroup, tested_counts$n),
  summary_block("tested_prop", tested_counts$TestedGroup, sprintf("%.6f", tested_counts$prop)),
  numeric_summary_rows(numeric_summary)
)
write.csv(summary_rows, file.path(tables_dir, "eda_summary.csv"), row.names = FALSE)

bodyweight_spread <- main_sex_df |>
  group_by(Sex, BodyweightBin) |>
  summarise(n = n(), mean_total = mean(TotalKg), sd_total = sd(TotalKg), .groups = "drop")

age_spread <- main_sex_df |>
  group_by(Sex, AgeBin) |>
  summarise(n = n(), mean_total = mean(TotalKg), sd_total = sd(TotalKg), .groups = "drop")

weight_class_summary <- main_sex_df |>
  filter(!is.na(WeightClassSimple)) |>
  group_by(Sex, WeightClassKg, WeightClassSimple) |>
  summarise(n = n(), mean_total = mean(TotalKg), sd_total = sd(TotalKg), .groups = "drop") |>
  arrange(Sex, WeightClassSimple)

age_levels <- main_sex_df |>
  filter(AgeClass != "") |>
  transmute(AgeClass, lower = as.numeric(sub("-.*", "", AgeClass))) |>
  distinct() |>
  arrange(lower) |>
  pull(AgeClass)

age_class_summary <- main_sex_df |>
  filter(AgeClass != "") |>
  group_by(Sex, AgeClass) |>
  summarise(n = n(), mean_total = mean(TotalKg), sd_total = sd(TotalKg), .groups = "drop") |>
  mutate(AgeClass = factor(AgeClass, levels = age_levels))

group_summaries <- bind_rows(
  bodyweight_spread |> transmute(summary_type = "bodyweight_spread", group_1 = Sex, group_2 = as.character(BodyweightBin), n, mean_total, sd_total),
  age_spread |> transmute(summary_type = "age_spread", group_1 = Sex, group_2 = as.character(AgeBin), n, mean_total, sd_total),
  weight_class_summary |> transmute(summary_type = "weight_class_summary", group_1 = Sex, group_2 = WeightClassKg, n, mean_total, sd_total),
  age_class_summary |> transmute(summary_type = "age_class_summary", group_1 = Sex, group_2 = as.character(AgeClass), n, mean_total, sd_total)
)
write.csv(group_summaries, file.path(tables_dir, "eda_group_summaries.csv"), row.names = FALSE)

correlation_vars <- c("Age", "BodyweightKg", "Best3SquatKg", "Best3BenchKg", "Best3DeadliftKg", "TotalKg")
correlation_df <- as.data.frame(as.table(cor(select(power_df, all_of(correlation_vars)), use = "complete.obs")), stringsAsFactors = FALSE)
names(correlation_df) <- c("var_x", "var_y", "correlation")
correlation_df$var_x <- factor(correlation_df$var_x, levels = correlation_vars)
correlation_df$var_y <- factor(correlation_df$var_y, levels = rev(correlation_vars))

set.seed(405)
plot_df <- main_sex_df[sample(nrow(main_sex_df), min(20000, nrow(main_sex_df))), ]

save_plot(
  ggplot(power_df, aes(x = Age)) +
    geom_histogram(bins = 40, fill = "steelblue", color = "white") +
    labs(title = "Age Distribution") +
    theme_minimal(),
  file.path(figures_dir, "age_distribution.png"),
  7,
  4
)

save_plot(
  ggplot(main_sex_df, aes(x = BodyweightKg, fill = Sex, color = Sex)) +
    geom_density(alpha = 0.25) +
    scale_fill_manual(values = sex_colors) +
    scale_color_manual(values = sex_colors) +
    labs(title = "Bodyweight Distribution by Sex", x = "Bodyweight (kg)") +
    theme_minimal(),
  file.path(figures_dir, "bodyweight_distribution.png"),
  7,
  4
)

save_plot(
  ggplot(main_sex_df, aes(x = TotalKg, fill = Sex, color = Sex)) +
    geom_density(alpha = 0.25) +
    scale_fill_manual(values = sex_colors) +
    scale_color_manual(values = sex_colors) +
    labs(title = "Total Distribution by Sex", subtitle = "Mx omitted because n = 43", x = "Total (kg)") +
    theme_minimal(),
  file.path(figures_dir, "total_distribution.png"),
  7,
  4
)

save_plot(
  ggplot(plot_df, aes(x = BodyweightKg, y = TotalKg, color = Sex)) +
    geom_point(alpha = 0.18, size = 0.6) +
    geom_smooth(se = FALSE, linewidth = 0.9) +
    scale_color_manual(values = sex_colors) +
    labs(title = "Total vs Bodyweight", subtitle = "20,000-point sample for readability", x = "Bodyweight (kg)", y = "Total (kg)") +
    theme_minimal(),
  file.path(figures_dir, "total_vs_bodyweight.png"),
  7,
  5
)

save_plot(
  ggplot(plot_df, aes(x = Age, y = TotalKg, color = Sex)) +
    geom_point(alpha = 0.18, size = 0.6) +
    geom_smooth(se = FALSE, linewidth = 0.9) +
    scale_color_manual(values = sex_colors) +
    labs(title = "Total vs Age", subtitle = "20,000-point sample for readability", x = "Age", y = "Total (kg)") +
    theme_minimal(),
  file.path(figures_dir, "total_vs_age.png"),
  7,
  5
)

save_plot(
  ggplot(plot_df, aes(x = Age, y = BodyweightKg, color = Sex)) +
    geom_point(alpha = 0.18, size = 0.6) +
    geom_smooth(se = FALSE, linewidth = 0.9) +
    scale_color_manual(values = sex_colors) +
    labs(title = "Bodyweight vs Age", subtitle = "20,000-point sample for readability", x = "Age", y = "Bodyweight (kg)") +
    theme_minimal(),
  file.path(figures_dir, "bodyweight_vs_age.png"),
  7,
  5
)

save_plot(
  ggplot(correlation_df, aes(x = var_x, y = var_y, fill = correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", correlation)), size = 3) +
    scale_fill_gradient2(low = "#3b82f6", mid = "white", high = "#ef4444", midpoint = 0, limits = c(-1, 1)) +
    labs(title = "Correlation Heatmap", x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  file.path(figures_dir, "correlation_heatmap.png"),
  7,
  5.5
)

save_plot(
  ggplot(main_sex_df, aes(x = TestedGroup, y = TotalKg, fill = TestedGroup)) +
    geom_boxplot(outlier.alpha = 0.03) +
    facet_wrap(~Sex) +
    scale_fill_manual(values = c("Blank" = "#94a3b8", "Yes" = "#60a5fa")) +
    labs(title = "Total by Tested Group and Sex", x = "Tested group", y = "Total (kg)") +
    theme_minimal() +
    theme(legend.position = "none"),
  file.path(figures_dir, "total_by_tested_and_sex.png"),
  7,
  5
)

save_plot(
  ggplot(bodyweight_spread, aes(x = BodyweightBin, y = sd_total, group = Sex, color = Sex)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = sex_colors) +
    labs(title = "Spread of Total Across Bodyweight Bins", x = "Bodyweight bin", y = "SD of total (kg)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  file.path(figures_dir, "total_sd_by_bodyweight_bin.png"),
  8,
  4.8
)

save_plot(
  ggplot(age_spread, aes(x = AgeBin, y = sd_total, group = Sex, color = Sex)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = sex_colors) +
    labs(title = "Spread of Total Across Age Bins", x = "Age bin", y = "SD of total (kg)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  file.path(figures_dir, "total_sd_by_age_bin.png"),
  8,
  4.8
)

save_plot(
  ggplot(weight_class_summary, aes(x = WeightClassSimple, y = mean_total, color = Sex)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = sex_colors) +
    labs(title = "Mean Total by Weight Class", x = "Weight class (kg)", y = "Mean total (kg)") +
    theme_minimal(),
  file.path(figures_dir, "mean_total_by_weight_class.png"),
  7,
  4.8
)

save_plot(
  ggplot(age_class_summary, aes(x = AgeClass, y = mean_total, color = Sex, group = Sex)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = sex_colors) +
    labs(title = "Mean Total by Age Class", x = "Age class", y = "Mean total (kg)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  file.path(figures_dir, "mean_total_by_age_class.png"),
  8,
  4.8
)

print(filter(summary_rows, section == "overview"))
print(sex_counts)
print(tested_counts)
print(numeric_summary)

# EDA notes for the frozen 2025 SBD raw subset:
# - The working dataset has 113,293 rows and the core modeling fields have no missing values after filtering.
# - Sex counts are F = 39,149, M = 74,101, and Mx = 43, so Mx is about 0.038% of the dataset.
# - Decision: exclude Mx from the main Bayesian models because the sample is too small for stable estimation.
# - Tested is recorded as Yes or blank, so any Tested analysis should stay descriptive / exploratory.
# - The curved lines in the scatterplots are smoothed trend lines from ggplot; they are just EDA guides, not fitted Bayesian models.
# - Total increases strongly with bodyweight, while the age relationship looks weaker and worth checking for nonlinearity later.
# - The bodyweight-binned and age-binned spread plots help check whether the variance in Total changes across the predictor range.
# - Weight-class and age-class summaries are saved because they connect directly to the project questions and later report tables.
# - Bodyweight and Total are moderately correlated overall; Total is, by construction, extremely correlated with the three lift components.
