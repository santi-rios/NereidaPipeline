# Enhanced Marine Biodiversity Research Pipeline with targets
# Comprehensive workflow integrating multiple data sources and analysis methods

library(targets)
library(tarchetypes)

# Load enhanced packages for marine biodiversity research
suppressPackageStartupMessages({
  library(robis)           # Ocean Biodiversity Information System
  library(rgbif)           # Global Biodiversity Information Facility
  library(dplyr)           # Data manipulation
  library(CoordinateCleaner)  # Coordinate cleaning and validation
  library(taxa)            # Taxonomic data management
  library(terra)           # Spatial data processing
  library(nodbi)           # NoSQL database interface
  library(biomartr)         # For genomic sequence retrieval
})

# Conditional loading of optional packages
optional_packages <- c("biomartr", "wikitaxa", "prism", "myTAI", "geotargets", "fastqcr")
for (pkg in optional_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    message("Loaded optional package: ", pkg)
  } else {
    message("Optional package not available: ", pkg)
  }
}

# Source all enhanced R functions
source("R/data_acquisition.R")
source("R/taxonomic_management.R")
source("R/geospatial_processing.R")
source("R/data_cleaning.R")
source("R/database_integration.R")
source("R/evolutionary_analysis.R")
source("R/metagenomic_analysis.R")

# Configuration parameters for the marine biodiversity pipeline
MARINE_SPECIES <- c(
  "Acropora cervicornis",
  "Acropora palmata",
  "Porites astreoides"
  )
STUDY_REGION <- list(
  min_lat = 24.0, max_lat = 26.0,
  min_lon = -82.0, max_lon = -80.0
)

# Enhanced Marine Biodiversity Research Pipeline
list(

  # === MICROBIOME ANALYSIS PHASE ===

  # 1. Define file paths for microbiome data
  tar_target(otu_file, "data/otutab.txt", format = "file"),
  tar_target(tax_file, "data/taxonomy.txt", format = "file"),
  tar_target(meta_file, "data/metadata.tsv", format = "file"),
  tar_target(tree_file, "data/otus.tree", format = "file"),

  # 2. Create a phyloseq object
  tar_target(
    phyloseq_object,
    create_phyloseq_object(otu_file, tax_file, meta_file, tree_file)
  ),

  # === DATA ACQUISITION PHASE ===
  
  # 1. Collect occurrence data from multiple sources
  tar_target(
    obis_occurrences,
    {
      all_obis_data <- list()
      for (species in MARINE_SPECIES) {
        species_data <- get_obis_data(
          species = species,
          out_dir = "data/raw/obis",
          write = TRUE
        )
        if (!is.null(species_data) && nrow(species_data) > 0) {
          all_obis_data[[species]] <- species_data
        }
      }
      
      if (length(all_obis_data) > 0) {
        all_columns <- unique(unlist(lapply(all_obis_data, names)))
        standardized_data <- lapply(all_obis_data, function(df) {
          for (col in all_columns) {
            if (!col %in% names(df)) {
              df[[col]] <- NA
            }
          }
          df[, all_columns, drop = FALSE]
        })
        do.call(rbind, standardized_data)
      } else {
        data.frame()
      }
    }
  ),
  
  tar_target(
    gbif_occurrences,
    {
      all_gbif_data <- list()
      for (species in MARINE_SPECIES) {
        species_data <- get_gbif_marine_data(
          species = species,
          out_dir = "data/raw/gbif",
          marine_only = TRUE
        )
        if (!is.null(species_data) && nrow(species_data) > 0) {
          all_gbif_data[[species]] <- species_data
        }
      }
      
      if (length(all_gbif_data) > 0) {
        all_columns <- unique(unlist(lapply(all_gbif_data, names)))
        standardized_data <- lapply(all_gbif_data, function(df) {
          for (col in all_columns) {
            if (!col %in% names(df)) {
              df[[col]] <- NA
            }
          }
          df[, all_columns, drop = FALSE]
        })
        do.call(rbind, standardized_data)
      } else {
        data.frame()
      }
    }
  ),
  
  # 2. Collect taxonomic information
  tar_target(
    taxonomic_data,
    {
      all_taxonomy <- list()
      for (species in MARINE_SPECIES) {
        tax_data <- get_taxonomic_data(
          species = species,
          out_dir = "data/raw/taxonomy"
        )
        if (!is.null(tax_data)) {
          all_taxonomy[[species]] <- tax_data
        }
      }
      all_taxonomy
    }
  ),
  
  # 3. Combine occurrence datasets
  tar_target(
    combined_occurrences,
    {
      combined_data <- list()
      
      if (!is.null(obis_occurrences) && nrow(obis_occurrences) > 0) {
        combined_data[["obis"]] <- obis_occurrences
      }
      
      if (!is.null(gbif_occurrences) && nrow(gbif_occurrences) > 0) {
        combined_data[["gbif"]] <- gbif_occurrences
      }
      
      if (length(combined_data) > 0) {
        standardized_data <- lapply(combined_data, function(df) {
          required_cols <- c("scientificName", "decimalLatitude", "decimalLongitude")
          missing_cols <- setdiff(required_cols, names(df))
          for (col in missing_cols) {
            df[[col]] <- NA
          }
          return(df)
        })
        do.call(rbind, standardized_data)
      } else {
        data.frame()
      }
    }
  ),
  
  # === DATA CLEANING PHASE ===
  # 4. Clean and validate occurrence data
  tar_target(
    cleaned_occurrences,
    {
      if (nrow(combined_occurrences) == 0) {
        message("No combined occurrences to clean")
        return(data.frame())
      }
      cleaned <- clean_marine_data(combined_occurrences, cleaning_level = "moderate")
      result <- as.data.frame(cleaned)
      result
    }
  ),
  # 4a. Alternative: Basic cleaning (less strict) 
  tar_target(
    cleaned_occurrences_basic,
    {
      if (nrow(combined_occurrences) == 0) {
        message("No combined occurrences to clean")
        return(data.frame())
      }
      cleaned <- clean_marine_data(combined_occurrences, cleaning_level = "basic")
      as.data.frame(cleaned)
    }
  ),
  # 5. Create taxonomic data structure
  tar_target(
    marine_taxonomy,
    {
      if (nrow(cleaned_occurrences) > 0) {
        parse_marine_taxonomy(
          occurrence_data = cleaned_occurrences,
          taxonomy_data = taxonomic_data
        )
      } else {
        NULL
      }
    }
  ),
  
  # === SPATIAL ANALYSIS PHASE ===
  
  # 6. Create environmental layers
  tar_target(
    environmental_layers,
    {
      if (nrow(cleaned_occurrences) > 0) {
        coords <- cleaned_occurrences[
          !is.na(cleaned_occurrences$decimalLatitude) & 
          !is.na(cleaned_occurrences$decimalLongitude),
        ]
        
        if (nrow(coords) > 0) {
          extent_vec <- c(
            min(coords$decimalLongitude, na.rm = TRUE) - 0.5,
            max(coords$decimalLongitude, na.rm = TRUE) + 0.5,
            min(coords$decimalLatitude, na.rm = TRUE) - 0.5,
            max(coords$decimalLatitude, na.rm = TRUE) + 0.5
          )
          
          get_marine_environmental_layers(
            extent = extent_vec,
            variables = c("depth", "sst", "salinity"),
            output_dir = "data/raw/environmental"
          )
        } else {
          NULL
        }
      } else {
        NULL
      }
    },
    format = "file"
  ),
  
  # # 7. Spatial analysis and modeling
  tar_target(
    spatial_analysis,
    {
      env_layers_rasters <- terra::rast(environmental_layers)
      if (nrow(cleaned_occurrences) > 0 && !is.null(env_layers_rasters)) {
        analyze_marine_spatial_patterns(
          occurrence_data = cleaned_occurrences,
          env_layers = env_layers_rasters,
          output_dir = "data/processed/spatial"
        )
      } else {
        NULL
      }
    }
  ),
  
  # 8. Create habitat suitability models
  tar_target(
    habitat_models,
    {
      if (nrow(cleaned_occurrences) > 0 && !is.null(environmental_layers)) {
        env_layers_rasters <- terra::rast(environmental_layers)
        species_list <- unique(cleaned_occurrences$scientificName)
        species_list <- species_list[!is.na(species_list)]
        
        habitat_model_results <- list()
        
        for (species in species_list) {
          species_data <- cleaned_occurrences[
            cleaned_occurrences$scientificName == species &
            !is.na(cleaned_occurrences$scientificName),
          ]
          
          if (nrow(species_data) >= 10) {
            species_with_env <- extract_environmental_values(
              occurrence_data = species_data,
              env_layers = env_layers_rasters
            )
            
            model_result <- create_habitat_suitability_model(
              occurrence_data = species_with_env,
              env_layers = env_layers_rasters,
              species_column = "scientificName"
            )
            
            if (!is.null(model_result)) {
              habitat_model_results[[species]] <- model_result
            }
          }
        }
        
        habitat_model_results
      } else {
        NULL
      }
    }
  ),
  
  # === DATABASE INTEGRATION PHASE ===
# Instalar la extensión JSON de DuckDB
### DBI::dbExecute(duckdb::dbConnect(duckdb::duckdb()), 'INSTALL json;')
### DBI::dbExecute(duckdb::dbConnect(duckdb::duckdb()), 'LOAD json;')
  # 9. Integrate data into DuckDB database
  tar_target(
    marine_database,
    {
      # Install DuckDB extensions first
      tryCatch({
        temp_con <- duckdb::dbConnect(duckdb::duckdb())
        DBI::dbExecute(temp_con, 'INSTALL json;')
        DBI::dbExecute(temp_con, 'LOAD json;')
        DBI::dbDisconnect(temp_con)
        message("DuckDB JSON extension installed successfully")
      }, error = function(e) {
        warning("Failed to install DuckDB extensions: ", e$message)
      })
      
      # Create database directory
      db_dir <- "data"
      if (!dir.exists(db_dir)) {
        dir.create(db_dir, recursive = TRUE)
      }
      
      # Initialize nodbi connection
      db_src <- nodbi::src_duckdb(dbdir = "data/marine_biodiversity.duckdb")
      
      # Store occurrence data directly with nodbi
      if (!is.null(cleaned_occurrences) && nrow(cleaned_occurrences) > 0) {
        message("Storing occurrence data...")
        
        # Prepare occurrence data for storage
        occurrence_list <- cleaned_occurrences %>%
          dplyr::mutate(
            id = row_number(),
            stored_at = Sys.time()
          ) %>%
          dplyr::select(
            id, scientificName, decimalLatitude, decimalLongitude, 
            eventDate, datasetName, stored_at, everything()
          )
        
        # Convert to list format for nodbi
        occurrence_records <- split(occurrence_list, seq_len(nrow(occurrence_list)))
        
        tryCatch({
          nodbi::docdb_create(db_src, "occurrences", occurrence_records)
          message("Successfully stored ", length(occurrence_records), " occurrence records")
        }, error = function(e) {
          message("Error storing occurrence data: ", e$message)
        })
      }
      
      # Store taxonomic data directly
      if (!is.null(taxonomic_data) && length(taxonomic_data) > 0) {
        message("Storing taxonomic data...")
        
        # Convert taxonomic data to appropriate format
        taxonomic_list <- lapply(names(taxonomic_data), function(species) {
          list(
            species = species,
            taxonomic_info = taxonomic_data[[species]],
            stored_at = Sys.time()
          )
        })
        
        tryCatch({
          nodbi::docdb_create(db_src, "taxonomy", taxonomic_list)
          message("Successfully stored taxonomic data for ", length(taxonomic_list), " species")
        }, error = function(e) {
          message("Error storing taxonomic data: ", e$message)
        })
      }
      
      # Store spatial analysis results directly
      if (!is.null(spatial_analysis) && !is.null(spatial_analysis$spatial_statistics)) {
        message("Storing spatial analysis results...")
        
        spatial_record <- list(
          analysis_type = "spatial_summary",
          results = spatial_analysis$spatial_statistics,
          stored_at = Sys.time()
        )
        
        tryCatch({
          nodbi::docdb_create(db_src, "spatial_analysis", list(spatial_record))
          message("Successfully stored spatial analysis results")
        }, error = function(e) {
          message("Error storing spatial analysis: ", e$message)
        })
      }
      
      # Store environmental data if available
      if (!is.null(spatial_analysis) && !is.null(spatial_analysis$occurrence_environmental)) {
        message("Storing environmental data...")
        
        env_data <- spatial_analysis$occurrence_environmental %>%
          dplyr::mutate(
            id = row_number(),
            stored_at = Sys.time()
          )
        
        env_records <- split(env_data, seq_len(nrow(env_data)))
        
        tryCatch({
          nodbi::docdb_create(db_src, "environmental_data", env_records)
          message("Successfully stored ", length(env_records), " environmental records")
        }, error = function(e) {
          message("Error storing environmental data: ", e$message)
        })
      }
      
      message("Database integration completed successfully")
      db_src
    }
  ),
  # === SUMMARY AND REPORTING PHASE ===
  # 10. Generate comprehensive summaries
  tar_target(
    data_summaries,
    {
      summaries <- list()
      
      # Species summary (direct from cleaned data)
      if (!is.null(cleaned_occurrences) && nrow(cleaned_occurrences) > 0) {
        summaries$species <- cleaned_occurrences %>%
          dplyr::count(scientificName, sort = TRUE, name = "n_records") %>%
          dplyr::filter(!is.na(scientificName))
        
        message("Generated species summary with ", nrow(summaries$species), " species")
      }
      
      # Spatial summary (direct from cleaned data)
      if (!is.null(cleaned_occurrences) && nrow(cleaned_occurrences) > 0) {
        valid_coords <- cleaned_occurrences %>%
          dplyr::filter(
            !is.na(decimalLatitude) & !is.na(decimalLongitude),
            decimalLatitude >= -90 & decimalLatitude <= 90,
            decimalLongitude >= -180 & decimalLongitude <= 180
          )
        
        if (nrow(valid_coords) > 0) {
          summaries$spatial <- list(
            n_records_with_coords = nrow(valid_coords),
            lat_range = range(valid_coords$decimalLatitude),
            lon_range = range(valid_coords$decimalLongitude),
            n_unique_locations = nrow(dplyr::distinct(valid_coords, decimalLatitude, decimalLongitude))
          )
        }
        
        message("Generated spatial summary")
      }
      
      # Taxonomic summary (from direct data)
      if (!is.null(marine_taxonomy)) {
        tryCatch({
          summaries$taxonomy <- summarize_marine_taxonomy(marine_taxonomy)
          message("Generated taxonomic summary")
        }, error = function(e) {
          message("Taxonomic summary failed: ", e$message)
          summaries$taxonomy <- list(message = "Taxonomic summary not available")
        })
      }
      
      # Spatial analysis summary (from direct data)
      if (!is.null(spatial_analysis)) {
        summaries$spatial_analysis <- spatial_analysis$spatial_statistics
        message("Added spatial analysis summary")
      }
      
      # Habitat model summary (from direct data)
      if (!is.null(habitat_models) && length(habitat_models) > 0) {
        model_summary <- lapply(habitat_models, function(model) {
          list(
            n_presence_points = sum(model$model_data$presence == 1, na.rm = TRUE),
            n_background_points = sum(model$model_data$presence == 0, na.rm = TRUE),
            model_AIC = if (inherits(model$model, "glm")) AIC(model$model) else NA,
            env_variables = model$env_variables
          )
        })
        summaries$habitat_models <- model_summary
        message("Added habitat models summary for ", length(habitat_models), " models")
      }
      
      summaries
    }
  ),
  
  # 11. Export results to multiple formats
  tar_target(
    exported_results,
    {
      export_paths <- c()
      
      # Create output directory
      output_dir <- "data/processed/final"
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
      }
      
      # Export cleaned occurrence data directly (not from database)
      if (!is.null(cleaned_occurrences) && nrow(cleaned_occurrences) > 0) {
        # Export as CSV
        csv_path <- file.path(output_dir, "marine_occurrences.csv")
        write.csv(cleaned_occurrences, csv_path, row.names = FALSE)
        export_paths <- c(export_paths, "marine_occurrences.csv")
        
        # Export as JSON
        json_path <- file.path(output_dir, "marine_occurrences.json")
        jsonlite::write_json(cleaned_occurrences, json_path, pretty = TRUE)
        export_paths <- c(export_paths, "marine_occurrences.json")
        
        message("Exported occurrence data to CSV and JSON formats")
      }
      
      # Export taxonomic data
      if (!is.null(taxonomic_data) && length(taxonomic_data) > 0) {
        tax_path <- file.path(output_dir, "taxonomic_data.json")
        jsonlite::write_json(taxonomic_data, tax_path, pretty = TRUE)
        export_paths <- c(export_paths, "taxonomic_data.json")
        
        message("Exported taxonomic data")
      }
      
      # Export spatial analysis results (only serializable parts)
      if (!is.null(spatial_analysis)) {
        # Extract only serializable components from spatial_analysis
        spatial_export <- list()
        
        # Add spatial statistics if available
        if (!is.null(spatial_analysis$spatial_statistics)) {
          spatial_export$spatial_statistics <- spatial_analysis$spatial_statistics
        }
        
        # Add occurrence environmental data if available (as data.frame)
        if (!is.null(spatial_analysis$occurrence_environmental)) {
          spatial_export$occurrence_environmental <- as.data.frame(spatial_analysis$occurrence_environmental)
        }
        
        # Add any other non-SpatRaster components
        non_raster_components <- spatial_analysis[
          !sapply(spatial_analysis, function(x) inherits(x, "SpatRaster"))
        ]
        
        # Merge with our manual extractions
        spatial_export <- c(spatial_export, non_raster_components[
          !names(non_raster_components) %in% c("spatial_statistics", "occurrence_environmental")
        ])
        
        spatial_path <- file.path(output_dir, "spatial_analysis.json")
        jsonlite::write_json(spatial_export, spatial_path, pretty = TRUE)
        export_paths <- c(export_paths, "spatial_analysis.json")
        
        message("Exported spatial analysis results (serializable components)")
      }
      
      # Export habitat models results (only serializable parts)
      if (!is.null(habitat_models) && length(habitat_models) > 0) {
        # Extract serializable components from habitat models
        models_export <- lapply(names(habitat_models), function(species_name) {
          model <- habitat_models[[species_name]]
          serializable_model <- list(
            species = species_name
          )
          
          # Add model statistics if available
          if (!is.null(model$model) && inherits(model$model, c("lm", "glm"))) {
            serializable_model$model_summary <- list(
              coefficients = as.list(coef(model$model)),
              AIC = AIC(model$model),
              deviance = deviance(model$model),
              formula = as.character(formula(model$model))[1],
              residual_deviance = model$model$deviance,
              null_deviance = model$model$null.deviance,
              df_residual = model$model$df.residual,
              df_null = model$model$df.null
            )
            
            # Add model performance metrics if available
            if (!is.null(model$model$fitted.values)) {
              serializable_model$model_summary$n_fitted <- length(model$model$fitted.values)
            }
          }
          
          # Add model data (as data.frame, excluding complex objects)
          if (!is.null(model$model_data)) {
            model_data_clean <- as.data.frame(model$model_data)
            # Remove any columns that might contain complex objects
            complex_cols <- sapply(model_data_clean, function(x) {
              inherits(x, c("SpatRaster", "lm", "glm", "list")) && !is.vector(x)
            })
            model_data_clean <- model_data_clean[, !complex_cols, drop = FALSE]
            serializable_model$model_data <- model_data_clean
          }
          
          # Add environmental variables used
          if (!is.null(model$env_variables)) {
            serializable_model$env_variables <- as.character(model$env_variables)
          }
          
          # Add other simple components (excluding complex objects)
          other_components <- model[!names(model) %in% c("model", "model_data", "env_variables")]
          for (comp_name in names(other_components)) {
            comp_value <- other_components[[comp_name]]
            # Only add if it's a simple object (not SpatRaster, lm, glm, etc.)
            if (!inherits(comp_value, c("SpatRaster", "lm", "glm")) && 
                (is.atomic(comp_value) || is.list(comp_value))) {
              # For lists, check if they contain only simple objects
              if (is.list(comp_value)) {
                simple_list <- all(sapply(comp_value, function(x) {
                  is.atomic(x) || is.null(x)
                }))
                if (simple_list) {
                  serializable_model[[comp_name]] <- comp_value
                }
              } else {
                serializable_model[[comp_name]] <- comp_value
              }
            }
          }
          
          return(serializable_model)
        })
        
        names(models_export) <- names(habitat_models)
        
        models_path <- file.path(output_dir, "habitat_models.json")
        jsonlite::write_json(models_export, models_path, pretty = TRUE)
        export_paths <- c(export_paths, "habitat_models.json")
        
        message("Exported habitat models (serializable components)")
      }
      
      # Save summaries
      if (!is.null(data_summaries)) {
        summaries_path <- file.path(output_dir, "analysis_summaries.json")
        jsonlite::write_json(data_summaries, summaries_path, pretty = TRUE)
        export_paths <- c(export_paths, "analysis_summaries.json")
        
        message("Exported analysis summaries")
      }
      
      export_paths
    }
  ),
  
  # 12. Generate final report data
  tar_target(
    pipeline_report,
    {
      report_data <- list(
        pipeline_completion_date = Sys.time(),
        species_analyzed = MARINE_SPECIES,
        study_region = STUDY_REGION,
        
        data_acquisition = list(
          obis_records = if (!is.null(obis_occurrences)) nrow(obis_occurrences) else 0,
          gbif_records = if (!is.null(gbif_occurrences)) nrow(gbif_occurrences) else 0,
          combined_records = if (!is.null(combined_occurrences)) nrow(combined_occurrences) else 0,
          cleaned_records = if (!is.null(cleaned_occurrences)) nrow(cleaned_occurrences) else 0
        ),
        
        analysis_results = list(
          taxonomic_analysis_completed = !is.null(marine_taxonomy),
          spatial_analysis_completed = !is.null(spatial_analysis),
          habitat_models_created = if (!is.null(habitat_models)) length(habitat_models) else 0,
          database_integration_completed = !is.null(marine_database)
        ),
        
        output_files = exported_results,
        
        pipeline_status = "completed"
      )
      
      # Save pipeline report
      jsonlite::write_json(
        report_data,
        "data/processed/final/pipeline_report.json",
        pretty = TRUE
      )
      
      report_data
    }
  ),
  # === REPORTING PHASE ===
  # 13. Generate Quarto report
  tar_quarto(
    reporte_biodiversidad,
    path = "./reportes/analisis_biodiversidad_marina.qmd",
    quiet = FALSE
  )
)