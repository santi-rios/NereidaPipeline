# Load the data
marine_data <- read.csv("./data/raw/acropora_data.csv")

# Perform analysis (e.g., summary statistics)
summary_stats <- marine_data |>
  dplyr::group_by(scientificName) |>
  dplyr::summarize(count = dplyr::n())

# dplyr::glimpse(summary_stats)

# Save summary statistics
write.csv(summary_stats, file.path("./data/processed/acropora_summary_stats.csv"), row.names = FALSE)
