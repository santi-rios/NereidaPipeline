# Enhanced Marine Biodiversity Research Pipeline Examples

This document provides comprehensive examples and workflows for using the enhanced marine biodiversity research pipeline, focusing on conservation and evolutionary biology applications.

## Overview

The enhanced pipeline integrates multiple R packages to provide a comprehensive workflow for marine biodiversity research:

- **biomaRt**: Genomic sequence retrieval and analysis
- **wikitaxa**: Taxonomic information from Wikipedia
- **PRISM**: Environmental and climate data
- **taxa**: Standardized taxonomic data management
- **geotargets**: Geospatial data processing with targets
- **CoordinateCleaner**: Data quality control for occurrence records
- **nodbi**: NoSQL database integration with DuckDB
- **myTAI**: Evolutionary transcriptomics analysis

## Quick Start Example

```r
# Load the enhanced pipeline
source("_targets.R")

# Run the complete pipeline for marine coral species
tar_make()

# Check pipeline status
tar_visnetwork()

# Load results
cleaned_data <- tar_read(cleaned_occurrences)
spatial_results <- tar_read(spatial_analysis)
database <- tar_read(marine_database)
```

## Example 1: Coral Conservation Assessment

This example demonstrates how to assess the conservation status of coral species using occurrence data, environmental factors, and spatial analysis.

```r
# Define coral species of conservation concern
coral_species <- c(
  "Acropora cervicornis",   # Staghorn coral (Critically Endangered)
  "Acropora palmata",       # Elkhorn coral (Critically Endangered)
  "Orbicella faveolata"     # Mountainous star coral (Endangered)
)

# 1. Collect comprehensive occurrence data
coral_data <- list()
for (species in coral_species) {
  cat("Collecting data for", species, "\n")
  
  # Get OBIS data
  obis_data <- get_obis_data(species, out_dir = "data/coral_conservation")
  
  # Get GBIF marine data
  gbif_data <- get_gbif_marine_data(species, marine_only = TRUE)
  
  # Get taxonomic information
  taxonomy <- get_taxonomic_data(species)
  
  # Combine data sources
  if (!is.null(obis_data) && !is.null(gbif_data)) {
    combined <- rbind(
      transform(obis_data, data_source = "OBIS"),
      transform(gbif_data, data_source = "GBIF")
    )
    coral_data[[species]] <- combined
  }
}

# 2. Clean and validate the data
coral_occurrences <- do.call(rbind, coral_data)

# Apply comprehensive cleaning
cleaning_results <- clean_marine_biodiversity_data(
  occurrence_data = coral_occurrences,
  coordinate_cleaning_level = "strict",  # Use strict cleaning for conservation
  dataset_cleaning = TRUE,
  output_dir = "data/coral_conservation/cleaned"
)

clean_coral_data <- cleaning_results$final_cleaned

# 3. Create spatial analysis for conservation assessment
# Define Caribbean study region (main coral habitat)
caribbean_extent <- c(-90, -60, 10, 30)  # lon_min, lon_max, lat_min, lat_max

# Get marine environmental layers
env_layers <- get_marine_environmental_layers(
  extent = caribbean_extent,
  variables = c("depth", "sst", "salinity", "chlorophyll"),
  output_dir = "data/coral_conservation/environmental"
)

# Perform spatial analysis
coral_spatial_analysis <- analyze_marine_spatial_patterns(
  occurrence_data = clean_coral_data,
  env_layers = env_layers,
  output_dir = "data/coral_conservation/spatial"
)

# 4. Create habitat suitability models for each species
coral_habitat_models <- list()

for (species in coral_species) {
  species_data <- clean_coral_data[clean_coral_data$scientificName == species, ]
  
  if (nrow(species_data) >= 20) {  # Minimum records for reliable modeling
    # Extract environmental values
    species_with_env <- extract_environmental_values(species_data, env_layers)
    
    # Create habitat model
    habitat_model <- create_habitat_suitability_model(
      occurrence_data = species_with_env,
      env_layers = env_layers,
      species_column = "scientificName"
    )
    
    coral_habitat_models[[species]] <- habitat_model
  }
}

# 5. Generate conservation assessment report
conservation_assessment <- list(
  assessment_date = Sys.Date(),
  species_assessed = coral_species,
  study_region = "Caribbean",
  
  data_summary = list(
    total_records = nrow(coral_occurrences),
    cleaned_records = nrow(clean_coral_data),
    data_quality_improvement = cleaning_results$cleaning_summary$quality_improvement
  ),
  
  spatial_metrics = coral_spatial_analysis$spatial_statistics,
  
  habitat_suitability = lapply(coral_habitat_models, function(model) {
    if (!is.null(model)) {
      list(
        model_performance = summary(model$model),
        suitable_area_km2 = sum(terra::values(model$prediction) > 0.5, na.rm = TRUE) * 100,  # Assuming 10x10km pixels
        environmental_preferences = model$env_variables
      )
    }
  }),
  
  conservation_recommendations = c(
    "Prioritize protection of areas with high habitat suitability scores",
    "Monitor temperature and depth preferences for climate change adaptation",
    "Establish marine protected areas in regions with highest occurrence density",
    "Implement coral restoration programs in historically suitable but now depleted areas"
  )
)

# Save conservation assessment
jsonlite::write_json(
  conservation_assessment,
  "data/coral_conservation/conservation_assessment_report.json",
  pretty = TRUE
)
```

## Example 2: Marine Invasive Species Analysis

This example shows how to analyze the spread and environmental preferences of marine invasive species.

```r
# Define invasive marine species
invasive_species <- c(
  "Caulerpa taxifolia",     # Killer algae
  "Carcinus maenas",        # European green crab
  "Undaria pinnatifida"     # Wakame seaweed
)

# 1. Collect temporal occurrence data to track spread
invasive_analysis <- list()

for (species in invasive_species) {
  cat("Analyzing invasive species:", species, "\n")
  
  # Get comprehensive occurrence data
  species_data <- collect_marine_species_data(
    species = species,
    include_sequences = TRUE,    # Get genomic data for population analysis
    include_climate = TRUE,      # Get environmental data
    out_dir = "data/invasive_species"
  )
  
  # Clean the occurrence data
  if (!is.null(species_data$obis) || !is.null(species_data$gbif)) {
    combined_occurrences <- rbind(
      if (!is.null(species_data$obis)) transform(species_data$obis, source = "OBIS"),
      if (!is.null(species_data$gbif)) transform(species_data$gbif, source = "GBIF")
    )
    
    # Apply coordinate cleaning
    cleaned_invasive <- clean_marine_coordinates(
      occurrence_data = combined_occurrences,
      clean_level = "standard"
    )
    
    # Temporal analysis of spread
    if ("year" %in% names(cleaned_invasive)) {
      temporal_spread <- aggregate(
        cbind(n_records = rep(1, nrow(cleaned_invasive)),
              lat_mean = cleaned_invasive$decimalLatitude,
              lon_mean = cleaned_invasive$decimalLongitude),
        by = list(year = cleaned_invasive$year),
        FUN = function(x) if (is.numeric(x)) mean(x, na.rm = TRUE) else sum(x, na.rm = TRUE)
      )
      
      # Calculate rate of geographic expansion
      temporal_spread$geographic_range <- NA
      for (i in seq_len(nrow(temporal_spread))) {
        year_data <- cleaned_invasive[cleaned_invasive$year == temporal_spread$year[i], ]
        if (nrow(year_data) > 1) {
          lat_range <- diff(range(year_data$decimalLatitude, na.rm = TRUE))
          lon_range <- diff(range(year_data$decimalLongitude, na.rm = TRUE))
          temporal_spread$geographic_range[i] <- sqrt(lat_range^2 + lon_range^2)
        }
      }
    }
    
    invasive_analysis[[species]] <- list(
      occurrence_data = cleaned_invasive,
      temporal_spread = if (exists("temporal_spread")) temporal_spread else NULL,
      genomic_data = species_data$sequences,
      climate_data = species_data$climate
    )
  }
}

# 2. Environmental niche modeling for invasive potential
invasion_risk_models <- list()

for (species in names(invasive_analysis)) {
  species_data <- invasive_analysis[[species]]$occurrence_data
  
  if (nrow(species_data) >= 30) {
    # Create global environmental layers for invasion risk assessment
    global_extent <- c(-180, 180, -70, 70)
    
    global_env_layers <- get_marine_environmental_layers(
      extent = global_extent,
      variables = c("sst", "depth", "salinity"),
      output_dir = paste0("data/invasive_species/", gsub(" ", "_", species))
    )
    
    # Model current distribution
    current_model <- create_habitat_suitability_model(
      occurrence_data = species_data,
      env_layers = global_env_layers
    )
    
    if (!is.null(current_model)) {
      # Identify areas at risk of invasion (high suitability, no current presence)
      suitability_raster <- current_model$prediction
      
      # Calculate invasion risk score
      high_suitability <- terra::values(suitability_raster) > 0.7
      invasion_risk_area <- sum(high_suitability, na.rm = TRUE) * 100  # km²
      
      invasion_risk_models[[species]] <- list(
        model = current_model,
        invasion_risk_area_km2 = invasion_risk_area,
        environmental_preferences = current_model$env_variables,
        management_priority = if (invasion_risk_area > 10000) "HIGH" else if (invasion_risk_area > 1000) "MEDIUM" else "LOW"
      )
    }
  }
}

# 3. Generate invasion risk assessment
invasion_assessment <- list(
  assessment_type = "Marine Invasive Species Risk Assessment",
  assessment_date = Sys.Date(),
  species_analyzed = invasive_species,
  
  temporal_analysis = lapply(invasive_analysis, function(x) x$temporal_spread),
  
  invasion_risk = invasion_risk_models,
  
  management_recommendations = list(
    high_risk_species = names(invasion_risk_models)[
      sapply(invasion_risk_models, function(x) x$management_priority == "HIGH")
    ],
    
    monitoring_priorities = c(
      "Establish early detection networks in high-risk areas",
      "Monitor environmental conditions that favor invasive establishment",
      "Implement rapid response protocols for new detections",
      "Assess genetic diversity in invasive populations for management strategies"
    ),
    
    prevention_strategies = c(
      "Strengthen ballast water management in high-risk ports",
      "Enhance aquaculture biosecurity protocols",
      "Improve public awareness and reporting systems",
      "Develop species-specific control methods"
    )
  )
)

# Save invasion risk assessment
jsonlite::write_json(
  invasion_assessment,
  "data/invasive_species/invasion_risk_assessment.json",
  pretty = TRUE
)
```

## Example 3: Evolutionary Analysis of Marine Adaptation

This example demonstrates phylostratigraphic analysis to understand evolutionary adaptations in marine organisms.

```r
# This example requires expression data - here we show the workflow
# In practice, you would obtain RNA-seq data from marine organisms

# 1. Prepare example expression data for marine coral
# (In real analysis, this would be loaded from RNA-seq results)

# Simulate coral expression data across developmental stages
set.seed(123)
n_genes <- 1000
n_stages <- 6

# Create example gene IDs and descriptions
gene_annotations <- data.frame(
  gene_id = paste0("gene_", sprintf("%04d", 1:n_genes)),
  description = sample(c(
    "ribosomal protein L5", "histone H3", "actin beta",
    "collagen alpha-1", "neural cell adhesion", "immune response protein",
    "species-specific protein", "heat shock protein", "calcium binding protein",
    "extracellular matrix protein", "cytoskeletal protein", "metabolic enzyme"
  ), n_genes, replace = TRUE),
  stringsAsFactors = FALSE
)

# Create expression matrix (genes x developmental stages)
expression_matrix <- matrix(
  rnorm(n_genes * n_stages, mean = 100, sd = 30),
  nrow = n_genes,
  ncol = n_stages
)
rownames(expression_matrix) <- gene_annotations$gene_id
colnames(expression_matrix) <- paste0("Stage_", 1:n_stages)

# Ensure positive expression values
expression_matrix[expression_matrix < 0] <- abs(expression_matrix[expression_matrix < 0])

# Convert to data frame
coral_expression <- as.data.frame(expression_matrix)

# 2. Create phylostratigraphic mapping
coral_phylostrata <- create_marine_phylostrata(
  gene_annotations = gene_annotations,
  organism_name = "Acropora_cervicornis",
  reference_phylogeny = "cnidaria"
)

# 3. Prepare phylostratigraphic expression set
coral_phylo_set <- prepare_marine_expression_set(
  expression_data = coral_expression,
  phylostrata_data = coral_phylostrata,
  gene_id_col = "gene_id"
)

# 4. Calculate transcriptome age index (TAI)
developmental_stages <- c(
  "Egg", "Early_Embryo", "Late_Embryo", 
  "Larva", "Settlement", "Juvenile"
)

coral_tai <- calculate_marine_tai(
  phylo_expression_set = coral_phylo_set,
  developmental_stages = developmental_stages
)

# 5. Perform comprehensive phylostratigraphic analysis
coral_evolution_analysis <- analyze_marine_phylostratigraphy(
  phylo_expression_set = coral_phylo_set,
  analysis_type = "development"
)

# 6. Create stress response analysis (simulate stress vs control)
# Create comparison groups
stress_groups <- c("control", "control", "control", "stress", "stress", "stress")

coral_stress_analysis <- analyze_marine_phylostratigraphy(
  phylo_expression_set = coral_phylo_set,
  comparison_groups = stress_groups,
  analysis_type = "stress"
)

# 7. Generate conservation insights
coral_conservation_insights <- generate_conservation_insights(
  phylo_analysis_results = coral_evolution_analysis,
  species_name = "Acropora cervicornis",
  conservation_context = "endangered"
)

coral_stress_insights <- generate_conservation_insights(
  phylo_analysis_results = coral_stress_analysis,
  species_name = "Acropora cervicornis",
  conservation_context = "endangered"
)

# 8. Comparative analysis across coral species
# Simulate expression data for multiple coral species
coral_species_names <- c("Acropora_cervicornis", "Acropora_palmata", "Porites_astreoides")
coral_expression_sets <- list()

for (species in coral_species_names) {
  # Create species-specific expression patterns
  species_expression <- coral_expression + 
    matrix(rnorm(n_genes * n_stages, mean = 0, sd = 10), nrow = n_genes)
  species_expression[species_expression < 0] <- abs(species_expression[species_expression < 0])
  
  # Create phylostratum mapping for this species
  species_phylostrata <- create_marine_phylostrata(
    gene_annotations = gene_annotations,
    organism_name = species,
    reference_phylogeny = "cnidaria"
  )
  
  # Create phyloexpression set
  species_phylo_set <- prepare_marine_expression_set(
    expression_data = as.data.frame(species_expression),
    phylostrata_data = species_phylostrata
  )
  
  coral_expression_sets[[species]] <- species_phylo_set
}

# Comparative evolutionary analysis
coral_comparative_analysis <- compare_marine_species_evolution(
  expression_sets = coral_expression_sets,
  species_names = coral_species_names
)

# 9. Generate comprehensive evolutionary report
evolutionary_report <- list(
  analysis_type = "Marine Coral Evolutionary Transcriptomics",
  species_analyzed = coral_species_names,
  analysis_date = Sys.Date(),
  
  developmental_analysis = list(
    tai_patterns = coral_tai,
    phylostratigraphy_results = coral_evolution_analysis,
    conservation_implications = coral_conservation_insights$conservation_interpretation
  ),
  
  stress_response_analysis = list(
    stress_phylostratigraphy = coral_stress_analysis,
    stress_conservation_insights = coral_stress_insights$stress_insights,
    adaptive_capacity = coral_stress_insights$stress_insights$adaptive_capacity
  ),
  
  comparative_evolution = coral_comparative_analysis,
  
  conservation_recommendations = c(
    coral_conservation_insights$conservation_recommendations,
    "Monitor expression of ancient gene families during bleaching events",
    "Assess transcriptomic diversity in coral restoration programs",
    "Use phylostratigraphic patterns to guide selective breeding for resilience"
  ),
  
  research_priorities = c(
    "Expand phylostratigraphic analysis to more coral species",
    "Integrate environmental stress experiments with evolutionary analysis",
    "Develop transcriptomic biomarkers for coral health assessment",
    "Study evolutionary responses to ocean acidification and warming"
  )
)

# Save evolutionary analysis report
jsonlite::write_json(
  evolutionary_report,
  "data/evolutionary_analysis/coral_evolution_report.json",
  pretty = TRUE
)
```

## Example 4: Integrated Database Workflow

This example demonstrates how to use the integrated DuckDB database for comprehensive marine biodiversity research.

```r
# 1. Initialize marine biodiversity database
marine_db <- init_marine_database(
  db_path = "data/comprehensive_marine_db.duckdb",
  create_schema = TRUE
)

# 2. Populate database with multiple data types
# Load example data from previous analyses
coral_data <- tar_read(cleaned_occurrences)
spatial_data <- tar_read(spatial_analysis)

# Store occurrence data
n_occurrences <- store_occurrence_data(
  src = marine_db,
  occurrence_data = coral_data,
  container_name = "coral_occurrences",
  metadata = list(
    project = "Coral Conservation Assessment",
    region = "Caribbean",
    data_quality = "cleaned_and_validated"
  )
)

# Store environmental data
if (!is.null(spatial_data$occurrence_environmental)) {
  n_environmental <- store_environmental_data(
    src = marine_db,
    environmental_data = spatial_data$occurrence_environmental,
    container_name = "environmental_measurements",
    spatial_metadata = list(
      resolution = "0.1_degrees",
      variables = c("depth", "sst", "salinity")
    )
  )
}

# Store taxonomic data
coral_taxonomy_data <- list(
  taxonomic_hierarchy = "Animalia > Cnidaria > Anthozoa > Scleractinia",
  conservation_status = "Critically Endangered",
  habitat_preferences = "Shallow reef environments",
  threat_factors = c("Ocean warming", "Ocean acidification", "Disease", "Pollution")
)

n_taxonomy <- store_taxonomic_data(
  src = marine_db,
  taxonomic_data = coral_taxonomy_data,
  container_name = "species_profiles"
)

# 3. Query database for research questions
# Query 1: Find all coral occurrences in a specific region
caribbean_corals <- query_occurrence_data(
  src = marine_db,
  container_name = "coral_occurrences",
  spatial_bounds = list(
    min_lat = 18.0, max_lat = 28.0,
    min_lon = -85.0, max_lon = -75.0
  )
)

# Query 2: Find occurrences for specific species
acropora_records <- query_occurrence_data(
  src = marine_db,
  container_name = "coral_occurrences",
  species = "Acropora"
)

# Query 3: Find recent observations (last 20 years)
recent_corals <- query_occurrence_data(
  src = marine_db,
  container_name = "coral_occurrences",
  date_range = c("2000-01-01", "2023-12-31")
)

# 4. Generate database summaries
species_summary <- create_data_summary(
  src = marine_db,
  summary_type = "species",
  container_name = "coral_occurrences"
)

spatial_summary <- create_data_summary(
  src = marine_db,
  summary_type = "spatial", 
  container_name = "coral_occurrences"
)

temporal_summary <- create_data_summary(
  src = marine_db,
  summary_type = "temporal",
  container_name = "coral_occurrences"
)

# 5. Export data for external analysis
# Export to CSV for GIS analysis
export_success_csv <- export_marine_data(
  src = marine_db,
  container_name = "coral_occurrences",
  output_format = "csv",
  output_path = "data/exports/coral_occurrences_for_gis.csv"
)

# Export to JSON for web applications
export_success_json <- export_marine_data(
  src = marine_db,
  container_name = "coral_occurrences",
  output_format = "json",
  output_path = "data/exports/coral_occurrences_web.json"
)

# 6. Create comprehensive research database report
database_report <- list(
  database_creation_date = Sys.Date(),
  database_path = "data/comprehensive_marine_db.duckdb",
  
  data_summary = list(
    occurrence_records = nrow(caribbean_corals),
    environmental_records = if (exists("n_environmental")) n_environmental else 0,
    taxonomic_profiles = n_taxonomy,
    species_analyzed = length(unique(caribbean_corals$scientificName))
  ),
  
  spatial_coverage = list(
    latitude_range = range(caribbean_corals$decimalLatitude, na.rm = TRUE),
    longitude_range = range(caribbean_corals$decimalLongitude, na.rm = TRUE),
    grid_cells_occupied = nrow(spatial_summary)
  ),
  
  temporal_coverage = list(
    year_range = range(temporal_summary$year, na.rm = TRUE),
    records_by_decade = temporal_summary
  ),
  
  research_applications = list(
    conservation_assessment = "Database supports coral conservation prioritization",
    climate_impact_studies = "Environmental data enables climate change research",
    biodiversity_monitoring = "Temporal data supports trend analysis",
    restoration_planning = "Spatial data guides restoration site selection"
  ),
  
  data_quality_metrics = list(
    coordinate_completeness = sum(!is.na(caribbean_corals$decimalLatitude)) / nrow(caribbean_corals),
    taxonomic_completeness = sum(!is.na(caribbean_corals$scientificName)) / nrow(caribbean_corals),
    temporal_completeness = sum(!is.na(caribbean_corals$year)) / nrow(caribbean_corals)
  )
)

# Save database report
jsonlite::write_json(
  database_report,
  "data/exports/marine_database_report.json",
  pretty = TRUE
)

# Close database connection
DBI::dbDisconnect(marine_db$con, shutdown = TRUE)
```

## Running the Complete Enhanced Pipeline

To run the complete enhanced pipeline with all new capabilities:

```r
# 1. Ensure all required packages are installed
required_packages <- c(
  "targets", "tarchetypes", "robis", "rgbif", "dplyr",
  "CoordinateCleaner", "taxa", "terra", "nodbi", "DBI",
  "jsonlite", "biomaRt", "wikitaxa", "prism", "myTAI"
)

# Install missing packages
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

# 2. Load the enhanced pipeline
source("_targets.R")

# 3. Visualize the pipeline
tar_visnetwork()

# 4. Run the complete pipeline
tar_make()

# 5. Check for any errors
tar_meta(fields = c("name", "error", "warnings"))

# 6. Access results
final_report <- tar_read(pipeline_report)
print(final_report$pipeline_status)

# 7. Load specific analysis results
cleaned_data <- tar_read(cleaned_occurrences)
spatial_results <- tar_read(spatial_analysis)
habitat_models <- tar_read(habitat_models)
database <- tar_read(marine_database)

cat("Pipeline completed successfully!\n")
cat("Cleaned records:", nrow(cleaned_data), "\n")
cat("Species modeled:", length(habitat_models), "\n")
cat("Analysis outputs saved to: data/processed/final/\n")
```

This enhanced pipeline provides a comprehensive framework for marine biodiversity research, integrating data acquisition, quality control, spatial analysis, evolutionary studies, and database management to support conservation and evolutionary biology research.