# Powerlifting Dataset preliminary Filtering and Preview
#
# This script assumes that you have already downloaded a CSV file containing the dataset.
# The script then filters to full‑power raw competitions in 2025 (excluding rows with NA entries or if individuals were disqualified) 
# The script also displays the first few rows.

file_path <- "openpowerlifting.csv"

power_df <- read.csv(file_path, stringsAsFactors = FALSE)

str(power_df)
head(power_df)
nrow(power_df)

subset_df <- subset(power_df,
                      !is.na(Date) &
                      substr(Date, 1, 4) == "2025" &
                      Event == "SBD" &
                      Equipment == "Raw" &
                      !is.na(TotalKg) &
                      !is.na(BodyweightKg) &
                      !is.na(Sex) &
                      !is.na(Age) &
                      Place != "DQ" &
                      Place != "DD" &
                      Place != "NS")
head(subset_df)
nrow(subset_df)
