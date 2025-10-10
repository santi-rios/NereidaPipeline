create_visualization <- function(data) {
  library(ggplot2)
  ggplot(data, aes(x = decimalLongitude, y = decimalLatitude)) +
    geom_point() +
    labs(title = "Distribució
    n de Acropora")
}

stats_summary_visualization <- function(summary_data) {
  library(ggplot2)
  ggplot(summary_data, aes(x = scientificName, y = count)) +
  theme_minimal() +
  labs(title = "Marine Species Occurrences", x = "Species", y = "Count")
}

create_visualization(Acropora)
stats_summary_visualization(summary_stats)
