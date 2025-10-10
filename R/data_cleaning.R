# Data Cleaning Pipeline using CoordinateCleaner
# For standardizing and validating marine occurrence records from biological collection databases

#' Comprehensive coordinate cleaning for marine occurrence data
#' @param occurrence_data Data frame with occurrence records
#' @param lon_col Character. Name of longitude column
#' @param lat_col Character. Name of latitude column  
#' @param species_col Character. Name of species column
#' @param country_col Character. Name of country column (optional)
#' @param clean_level Character. Level of cleaning ("basic", "standard", "strict")
#' @return Data frame with cleaned occurrence records and cleaning flags
clean_marine_coordinates <- function(occurrence_data,
                                   lon_col = "decimalLongitude",
                                   lat_col = "decimalLatitude", 
                                   species_col = "scientificName",
                                   country_col = "countryCode",
                                   clean_level = "standard") {
  
  if (!requireNamespace("CoordinateCleaner", quietly = TRUE)) {
    stop("Package 'CoordinateCleaner' is required for coordinate cleaning.")
  }
  
  # Ensure required columns exist
  required_cols <- c(lon_col, lat_col, species_col)
  missing_cols <- setdiff(required_cols, names(occurrence_data))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Standardize column names for CoordinateCleaner
  cleaned_data <- occurrence_data
  
  # Rename columns if necessary
  if (lon_col != "decimalLongitude") {
    names(cleaned_data)[names(cleaned_data) == lon_col] <- "decimalLongitude"
  }
  if (lat_col != "decimalLatitude") {
    names(cleaned_data)[names(cleaned_data) == lat_col] <- "decimalLatitude"
  }
  if (species_col != "species") {
    names(cleaned_data)[names(cleaned_data) == species_col] <- "species"
  }
  
  # Convert coordinates to numeric
  cleaned_data$decimalLongitude <- as.numeric(cleaned_data$decimalLongitude)
  cleaned_data$decimalLatitude <- as.numeric(cleaned_data$decimalLatitude)
  
  # Remove obviously invalid coordinates
  cleaned_data <- cleaned_data[
    !is.na(cleaned_data$decimalLongitude) & 
    !is.na(cleaned_data$decimalLatitude) &
    cleaned_data$decimalLatitude >= -90 & 
    cleaned_data$decimalLatitude <= 90 &
    cleaned_data$decimalLongitude >= -180 & 
    cleaned_data$decimalLongitude <= 180,
  ]
  
  if (nrow(cleaned_data) == 0) {
    warning("No valid coordinates remaining after initial cleaning")
    return(cleaned_data)
  }
  
  # Apply CoordinateCleaner tests based on cleaning level
  if (clean_level == "basic") {
    # Basic tests only
    cleaned_data <- tryCatch({
      temp_data <- CoordinateCleaner::cc_val(cleaned_data, verbose = FALSE)
      temp_data <- CoordinateCleaner::cc_equ(temp_data, verbose = FALSE)
      temp_data <- CoordinateCleaner::cc_zero(temp_data, verbose = FALSE)
      temp_data
    }, error = function(e) {
      warning("Basic coordinate cleaning failed: ", e$message)
      cleaned_data
    })
    
  } else if (clean_level == "standard") {
    # Standard cleaning protocol for marine data
    cleaned_data <- tryCatch({
      temp_data <- CoordinateCleaner::cc_val(cleaned_data, verbose = FALSE)
      temp_data <- CoordinateCleaner::cc_equ(temp_data, verbose = FALSE)
      temp_data <- CoordinateCleaner::cc_zero(temp_data, verbose = FALSE)
      temp_data <- CoordinateCleaner::cc_cap(temp_data, verbose = FALSE)
      temp_data <- CoordinateCleaner::cc_cen(temp_data, verbose = FALSE)
      temp_data <- CoordinateCleaner::cc_gbif(temp_data, verbose = FALSE)
      temp_data <- CoordinateCleaner::cc_inst(temp_data, verbose = FALSE)
      temp_data <- CoordinateCleaner::cc_sea(temp_data, verbose = FALSE)
      temp_data <- CoordinateCleaner::cc_outl(temp_data, verbose = FALSE)
      temp_data
    }, error = function(e) {
      warning("Standard coordinate cleaning failed: ", e$message)
      cleaned_data
    })
    
  } else if (clean_level == "strict") {
    # Strict cleaning including additional tests
    cleaned_data <- tryCatch({
      # First apply standard cleaning
      intermediate <- CoordinateCleaner::cc_val(cleaned_data, verbose = FALSE)
      intermediate <- CoordinateCleaner::cc_equ(intermediate, verbose = FALSE)
      intermediate <- CoordinateCleaner::cc_zero(intermediate, verbose = FALSE)
      intermediate <- CoordinateCleaner::cc_cap(intermediate, verbose = FALSE)
      intermediate <- CoordinateCleaner::cc_cen(intermediate, verbose = FALSE)
      intermediate <- CoordinateCleaner::cc_gbif(intermediate, verbose = FALSE)
      intermediate <- CoordinateCleaner::cc_inst(intermediate, verbose = FALSE)
      intermediate <- CoordinateCleaner::cc_sea(intermediate, verbose = FALSE)
      intermediate <- CoordinateCleaner::cc_outl(intermediate, verbose = FALSE)
      
      # Additional strict tests
      if (country_col %in% names(occurrence_data)) {
        # Country mismatch test
        intermediate <- CoordinateCleaner::cc_coun(
          intermediate, 
          iso3 = country_col,
          verbose = FALSE
        )
      }
      
      # Duplicate records
      intermediate <- CoordinateCleaner::cc_dupl(intermediate, verbose = FALSE)
      
      intermediate
      
    }, error = function(e) {
      warning("Strict coordinate cleaning failed: ", e$message)
      cleaned_data
    })
  }
  
  # Add cleaning metadata
  cleaned_data$cleaning_level <- clean_level
  cleaned_data$cleaning_date <- Sys.Date()
  cleaned_data$n_records_before_cleaning <- nrow(occurrence_data)
  cleaned_data$n_records_after_cleaning <- nrow(cleaned_data)
  
  return(cleaned_data)
}

#' Marine-specific data quality assessment
#' @param occurrence_data Data frame with marine occurrence records
#' @param depth_col Character. Name of depth column (optional)
#' @param habitat_col Character. Name of habitat column (optional)
#' @return List with data quality assessment results
assess_marine_data_quality <- function(occurrence_data,
                                     depth_col = "depth",
                                     habitat_col = "habitat") {
  
  assessment <- list()
  
  # Basic data completeness
  total_records <- nrow(occurrence_data)
  assessment$total_records <- total_records
  
  # Coordinate completeness
  valid_coords <- sum(
    !is.na(occurrence_data$decimalLatitude) & 
    !is.na(occurrence_data$decimalLongitude) &
    occurrence_data$decimalLatitude >= -90 & 
    occurrence_data$decimalLatitude <= 90 &
    occurrence_data$decimalLongitude >= -180 & 
    occurrence_data$decimalLongitude <= 180
  )
  
  assessment$coordinate_completeness <- valid_coords / total_records
  assessment$n_valid_coordinates <- valid_coords
  
  # Species name completeness
  if ("scientificName" %in% names(occurrence_data)) {
    valid_species <- sum(!is.na(occurrence_data$scientificName) & 
                        nzchar(occurrence_data$scientificName))
    assessment$species_completeness <- valid_species / total_records
    assessment$n_unique_species <- length(unique(occurrence_data$scientificName[
      !is.na(occurrence_data$scientificName)
    ]))
  }
  
  # Temporal completeness
  date_cols <- c("eventDate", "year", "month", "day")
  available_date_cols <- intersect(date_cols, names(occurrence_data))
  
  if (length(available_date_cols) > 0) {
    if ("year" %in% available_date_cols) {
      valid_years <- sum(!is.na(occurrence_data$year) & 
                        occurrence_data$year > 1700 & 
                        occurrence_data$year <= as.numeric(format(Sys.Date(), "%Y")))
      assessment$temporal_completeness <- valid_years / total_records
      assessment$year_range <- range(occurrence_data$year, na.rm = TRUE)
    }
  }
  
  # Marine-specific quality indicators
  if (depth_col %in% names(occurrence_data)) {
    depth_data <- as.numeric(occurrence_data[[depth_col]])
    marine_depth_records <- sum(!is.na(depth_data) & depth_data > 0)
    assessment$marine_depth_records <- marine_depth_records
    assessment$depth_completeness <- marine_depth_records / total_records
    
    if (marine_depth_records > 0) {
      assessment$depth_range <- range(depth_data[depth_data > 0], na.rm = TRUE)
      assessment$mean_depth <- mean(depth_data[depth_data > 0], na.rm = TRUE)
    }
  }
  
  # Habitat information
  if (habitat_col %in% names(occurrence_data)) {
    habitat_data <- occurrence_data[[habitat_col]]
    marine_matches <- grepl("marine|ocean|sea|reef|pelagic|benthic|coral", 
                           habitat_data, ignore.case = TRUE)
    marine_matches[is.na(marine_matches)] <- FALSE
    marine_habitats <- sum(marine_matches)
    assessment$marine_habitat_records <- marine_habitats
  }
  
  # Coordinate precision assessment
  if (valid_coords > 0) {
    coord_precision <- list()
    
    # Check for rounded coordinates (potential precision issues)
    lat_rounded <- sum(occurrence_data$decimalLatitude %% 1 == 0, na.rm = TRUE)
    lon_rounded <- sum(occurrence_data$decimalLongitude %% 1 == 0, na.rm = TRUE)
    
    coord_precision$rounded_lat <- lat_rounded / valid_coords
    coord_precision$rounded_lon <- lon_rounded / valid_coords
    
    # Check coordinate uncertainty if available
    if ("coordinateUncertaintyInMeters" %in% names(occurrence_data)) {
      uncertainty_data <- as.numeric(occurrence_data$coordinateUncertaintyInMeters)
      uncertainty_records <- sum(!is.na(uncertainty_data))
      
      if (uncertainty_records > 0) {
        coord_precision$uncertainty_completeness <- uncertainty_records / valid_coords
        coord_precision$mean_uncertainty_m <- mean(uncertainty_data, na.rm = TRUE)
        coord_precision$median_uncertainty_m <- median(uncertainty_data, na.rm = TRUE)
      }
    }
    
    assessment$coordinate_precision <- coord_precision
  }
  
  # Data source assessment
  if ("basisOfRecord" %in% names(occurrence_data)) {
    basis_counts <- table(occurrence_data$basisOfRecord, useNA = "ifany")
    assessment$basis_of_record <- as.list(basis_counts)
  }
  
  if ("institutionCode" %in% names(occurrence_data)) {
    institution_counts <- table(occurrence_data$institutionCode, useNA = "ifany")
    assessment$top_institutions <- head(sort(institution_counts, decreasing = TRUE), 10)
  }
  
  return(assessment)
}

#' Clean dataset-level issues for marine biodiversity data
#' @param occurrence_data Data frame with occurrence records
#' @param min_year Numeric. Minimum year for records
#' @param max_uncertainty_km Numeric. Maximum coordinate uncertainty in kilometers
#' @param exclude_fossils Logical. Whether to exclude fossil records
#' @return Data frame with cleaned occurrence records
clean_marine_dataset <- function(occurrence_data,
                                min_year = 1950,
                                max_uncertainty_km = 100,
                                exclude_fossils = TRUE) {
  
  cleaned_data <- occurrence_data
  n_original <- nrow(cleaned_data)
  
  # Remove records before minimum year
  if ("year" %in% names(cleaned_data)) {
    year_filter <- !is.na(cleaned_data$year) & cleaned_data$year >= min_year
    cleaned_data <- cleaned_data[year_filter, ]
    message("Removed ", sum(!year_filter, na.rm = TRUE), " records before ", min_year)
  }
  
  # Remove records with high coordinate uncertainty
  if ("coordinateUncertaintyInMeters" %in% names(cleaned_data)) {
    uncertainty_m <- as.numeric(cleaned_data$coordinateUncertaintyInMeters)
    uncertainty_filter <- is.na(uncertainty_m) | (uncertainty_m / 1000 <= max_uncertainty_km)
    cleaned_data <- cleaned_data[uncertainty_filter, ]
    message("Removed ", sum(!uncertainty_filter, na.rm = TRUE), 
           " records with uncertainty > ", max_uncertainty_km, " km")
  }
  
  # Filter by basis of record (exclude fossils if requested)
  if (exclude_fossils && "basisOfRecord" %in% names(cleaned_data)) {
    fossil_records <- grepl("fossil|preserved_specimen", 
                           cleaned_data$basisOfRecord, ignore.case = TRUE)
    
    # Keep human observations and living specimens
    suitable_basis <- cleaned_data$basisOfRecord %in% c(
      "HUMAN_OBSERVATION", "OBSERVATION", "LIVING_SPECIMEN",
      "MACHINE_OBSERVATION", "MATERIAL_SAMPLE"
    ) | !fossil_records
    
    cleaned_data <- cleaned_data[suitable_basis, ]
    message("Removed ", sum(!suitable_basis, na.rm = TRUE), " fossil/unsuitable records")
  }
  
  # Remove records with suspicious depth values (terrestrial depths in marine context)
  if ("depth" %in% names(cleaned_data)) {
    depth_values <- as.numeric(cleaned_data$depth)
    
    # For marine records, depth should be positive (below sea level)
    # Remove records with negative depth unless they're clearly intertidal
    suspicious_depth <- !is.na(depth_values) & depth_values < -10
    
    if (sum(suspicious_depth) > 0) {
      cleaned_data <- cleaned_data[!suspicious_depth, ]
      message("Removed ", sum(suspicious_depth), " records with suspicious depth values")
    }
  }
  
  # Remove duplicates based on species, coordinates, and date
  duplicate_cols <- c("scientificName", "decimalLatitude", "decimalLongitude")
  if ("eventDate" %in% names(cleaned_data)) {
    duplicate_cols <- c(duplicate_cols, "eventDate")
  }
  
  available_duplicate_cols <- intersect(duplicate_cols, names(cleaned_data))
  
  if (length(available_duplicate_cols) >= 3) {
    cleaned_data <- cleaned_data[!duplicated(cleaned_data[, available_duplicate_cols]), ]
    message("Removed ", nrow(occurrence_data) - nrow(cleaned_data), " duplicate records")
  }
  
  # Add cleaning summary
  n_final <- nrow(cleaned_data)
  
  message("Dataset cleaning summary:")
  message("  Original records: ", n_original)
  message("  Final records: ", n_final)
  message("  Records removed: ", n_original - n_final, " (", 
         round((n_original - n_final) / n_original * 100, 1), "%)")
  
  # Add metadata
  cleaned_data$dataset_cleaned <- TRUE
  cleaned_data$cleaning_date <- Sys.Date()
  cleaned_data$min_year_filter <- min_year
  cleaned_data$max_uncertainty_km_filter <- max_uncertainty_km
  
  return(cleaned_data)
}

#' Comprehensive marine data cleaning pipeline
#' @param occurrence_data Data frame with raw occurrence records
#' @param coordinate_cleaning_level Character. Level of coordinate cleaning
#' @param dataset_cleaning Logical. Whether to apply dataset-level cleaning
#' @param output_dir Character. Directory to save cleaning results
#' @return List with cleaned data and cleaning summary
clean_marine_biodiversity_data <- function(occurrence_data,
                                         coordinate_cleaning_level = "standard",
                                         dataset_cleaning = TRUE,
                                         output_dir = "./data/processed/cleaned") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  message("Starting comprehensive marine data cleaning pipeline...")
  
  cleaning_results <- list()
  
  # 1. Initial data quality assessment
  message("Assessing initial data quality...")
  initial_quality <- assess_marine_data_quality(occurrence_data)
  cleaning_results$initial_quality <- initial_quality
  
  # 2. Coordinate cleaning
  message("Cleaning coordinates...")
  coord_cleaned <- clean_marine_coordinates(
    occurrence_data,
    clean_level = coordinate_cleaning_level
  )
  cleaning_results$coordinate_cleaned <- coord_cleaned
  
  # 3. Dataset-level cleaning (optional)
  if (dataset_cleaning) {
    message("Applying dataset-level cleaning...")
    final_cleaned <- clean_marine_dataset(coord_cleaned)
    cleaning_results$final_cleaned <- final_cleaned
  } else {
    final_cleaned <- coord_cleaned
    cleaning_results$final_cleaned <- final_cleaned
  }
  
  # 4. Final data quality assessment
  message("Assessing final data quality...")
  final_quality <- assess_marine_data_quality(final_cleaned)
  cleaning_results$final_quality <- final_quality
  
  # 5. Create cleaning summary
  cleaning_summary <- list(
    original_records = nrow(occurrence_data),
    coordinate_cleaned_records = nrow(coord_cleaned),
    final_cleaned_records = nrow(final_cleaned),
    coordinate_cleaning_level = coordinate_cleaning_level,
    dataset_cleaning_applied = dataset_cleaning,
    cleaning_date = Sys.Date(),
    quality_improvement = list(
      coordinate_completeness = final_quality$coordinate_completeness - initial_quality$coordinate_completeness,
      species_completeness = final_quality$species_completeness - initial_quality$species_completeness
    )
  )
  
  cleaning_results$cleaning_summary <- cleaning_summary
  
  # 6. Save results
  # Save cleaned data
  utils::write.csv(
    final_cleaned,
    file.path(output_dir, "cleaned_marine_occurrences.csv"),
    row.names = FALSE
  )
  
  # Save cleaning summary
  jsonlite::write_json(
    cleaning_summary,
    file.path(output_dir, "cleaning_summary.json"),
    pretty = TRUE
  )
  
  # Save quality assessments
  jsonlite::write_json(
    list(initial = initial_quality, final = final_quality),
    file.path(output_dir, "quality_assessment.json"),
    pretty = TRUE
  )
  
  message("Data cleaning pipeline completed successfully!")
  message("Cleaned data saved to: ", file.path(output_dir, "cleaned_marine_occurrences.csv"))
  
  return(cleaning_results)
}

#' Enhanced marine data cleaning pipeline with configurable strictness
#' @param data Combined occurrence data from multiple sources
#' @param cleaning_level Character: "strict", "moderate", or "basic" 
#' @return Cleaned data.frame with quality assessment
clean_marine_data <- function(data, cleaning_level = "moderate") {
  message("Starting comprehensive marine data cleaning pipeline (", cleaning_level, " mode)...")
  
  # Ensure we have a proper data.frame
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }
  
  # Check if data is empty
  if (nrow(data) == 0) {
    warning("Input data is empty")
    return(data.frame())
  }
  
  message("Assessing initial data quality...")
  
  # Fix column name issues - remove duplicates and standardize
  message("Standardizing column names...")
  col_names <- names(data)
  
  # Remove duplicate columns (keep first occurrence)
  duplicate_cols <- duplicated(col_names)
  if (any(duplicate_cols)) {
    message("Removing ", sum(duplicate_cols), " duplicate columns")
    data <- data[, !duplicate_cols, drop = FALSE]
  }
  
  # Ensure required columns exist
  required_cols <- c("decimalLongitude", "decimalLatitude", "scientificName")
  missing_cols <- setdiff(required_cols, names(data))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Convert coordinates to numeric if they're not already
  data$decimalLongitude <- as.numeric(data$decimalLongitude)
  data$decimalLatitude <- as.numeric(data$decimalLatitude)
  
  # Remove rows with NA coordinates
  before_na <- nrow(data)
  data <- data[!is.na(data$decimalLongitude) & !is.na(data$decimalLatitude), ]
  after_na <- nrow(data)
  if (before_na != after_na) {
    message("Removed ", before_na - after_na, " records with missing coordinates")
  }
  
  if (nrow(data) == 0) {
    warning("No data remaining after removing NA coordinates")
    return(data.frame())
  }
  
  message("Cleaning coordinates...")
  
  # Apply cleaning based on level
  cleaned_data <- switch(cleaning_level,
    "basic" = clean_basic(data),
    "moderate" = clean_moderate(data),
    "strict" = clean_strict(data),
    clean_moderate(data) # default
  )
  
  message("Applying dataset-level cleaning...")
  
  # Dataset-level cleaning
  original_count <- nrow(cleaned_data)
  
  # Remove very old records (configurable threshold)
  year_threshold <- switch(cleaning_level,
    "basic" = 1900,
    "moderate" = 1950, 
    "strict" = 1980,
    1950
  )
  
  if ("year" %in% names(cleaned_data)) {
    old_records <- is.na(cleaned_data$year) | cleaned_data$year < year_threshold
    cleaned_data <- cleaned_data[!old_records, ]
    removed_old <- sum(old_records, na.rm = TRUE)
    message("Removed ", removed_old, " records before ", year_threshold)
  }
  
  # Remove records with high coordinate uncertainty (configurable)
  uncertainty_threshold <- switch(cleaning_level,
    "basic" = 500000,    # 500 km
    "moderate" = 100000, # 100 km  
    "strict" = 10000,    # 10 km
    100000
  )
  
  if ("coordinateUncertaintyInMeters" %in% names(cleaned_data)) {
    high_uncertainty <- !is.na(cleaned_data$coordinateUncertaintyInMeters) & 
                       cleaned_data$coordinateUncertaintyInMeters > uncertainty_threshold
    cleaned_data <- cleaned_data[!high_uncertainty, ]
    removed_uncertain <- sum(high_uncertainty, na.rm = TRUE)
    message("Removed ", removed_uncertain, " records with uncertainty > ", uncertainty_threshold/1000, " km")
  }
  
  # Remove fossil records and other unsuitable basis of record
  if ("basisOfRecord" %in% names(cleaned_data)) {
    unsuitable_basis <- cleaned_data$basisOfRecord %in% 
                       c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN", "UNKNOWN")
    cleaned_data <- cleaned_data[!unsuitable_basis, ]
    removed_fossil <- sum(unsuitable_basis, na.rm = TRUE)
    message("Removed ", removed_fossil, " fossil/unsuitable records")
  }
  
  # Remove duplicates based on coordinates and species
  before_dupl <- nrow(cleaned_data)
  if (nrow(cleaned_data) > 0) {
    # For basic cleaning, use broader coordinate rounding to catch near-duplicates
    if (cleaning_level == "basic") {
      cleaned_data$lon_rounded <- round(cleaned_data$decimalLongitude, 2)
      cleaned_data$lat_rounded <- round(cleaned_data$decimalLatitude, 2)
      cleaned_data <- cleaned_data[!duplicated(cleaned_data[, c("lon_rounded", "lat_rounded", "scientificName")]), ]
      cleaned_data$lon_rounded <- NULL
      cleaned_data$lat_rounded <- NULL
    } else {
      cleaned_data <- cleaned_data[!duplicated(cleaned_data[, c("decimalLongitude", 
                                                              "decimalLatitude", 
                                                              "scientificName")]), ]
    }
  }
  after_dupl <- nrow(cleaned_data)
  removed_dupl <- before_dupl - after_dupl
  message("Removed ", removed_dupl, " duplicate records")
  
  # Summary
  final_count <- nrow(cleaned_data)
  removed_total <- original_count - final_count
  removal_percent <- round((removed_total / original_count) * 100, 1)
  
  message("Dataset cleaning summary:")
  message("  Original records: ", original_count)
  message("  Final records: ", final_count)
  message("  Records removed: ", removed_total, " (", removal_percent, "%)")
  
  # Add quality flags for remaining data
  if (nrow(cleaned_data) > 0) {
    cleaned_data$quality_flag <- assess_record_quality(cleaned_data)
  }
  
  message("Assessing final data quality...")
  
  # Ensure we return a proper data.frame
  result <- as.data.frame(cleaned_data)
  
  # Add cleaning metadata as attributes
  attr(result, "cleaning_summary") <- list(
    cleaning_level = cleaning_level,
    original_records = nrow(data),
    final_records = final_count,
    records_removed = nrow(data) - final_count,
    removal_percentage = round(((nrow(data) - final_count) / nrow(data)) * 100, 1)
  )
  
  return(result)
}

# Helper functions for different cleaning levels
clean_basic <- function(data) {
  # Basic cleaning - only obvious errors
  data_clean <- data[data$decimalLongitude >= -180 & 
                    data$decimalLongitude <= 180 &
                    data$decimalLatitude >= -90 & 
                    data$decimalLatitude <= 90, ]
  
  # Only remove obvious zero coordinates
  data_clean <- data_clean[!(data_clean$decimalLongitude == 0 & data_clean$decimalLatitude == 0), ]
  
  return(data_clean)
}

clean_moderate <- function(data) {
  # Moderate cleaning - standard CoordinateCleaner tests
  tryCatch({
    data_clean <- data
    
    # Apply core CoordinateCleaner functions
    data_clean <- cc_val(data_clean, lon = "decimalLongitude", lat = "decimalLatitude")
    data_clean <- cc_equ(data_clean, lon = "decimalLongitude", lat = "decimalLatitude")  
    data_clean <- cc_zero(data_clean, lon = "decimalLongitude", lat = "decimalLatitude")
    data_clean <- cc_cen(data_clean, lon = "decimalLongitude", lat = "decimalLatitude")
    
    # For marine species, remove terrestrial points (but less aggressively)
    data_clean <- cc_sea(data_clean, lon = "decimalLongitude", lat = "decimalLatitude", 
                         ref = NULL, scale = 110) # Use coarser scale
    
    # Outlier detection per species
    if (length(unique(data_clean$scientificName)) > 1) {
      data_clean <- cc_outl(data_clean, lon = "decimalLongitude", lat = "decimalLatitude", 
                           species = "scientificName", method = "quantile", mltpl = 5)
    }
    
    return(data_clean)
    
  }, error = function(e) {
    warning("Moderate cleaning failed, falling back to basic: ", e$message)
    clean_basic(data)
  })
}

clean_strict <- function(data) {
  # Strict cleaning - all tests with conservative parameters
  tryCatch({
    data_clean <- clean_moderate(data)
    
    # Additional strict tests
    data_clean <- cc_cap(data_clean, lon = "decimalLongitude", lat = "decimalLatitude")
    data_clean <- cc_inst(data_clean, lon = "decimalLongitude", lat = "decimalLatitude")
    
    return(data_clean)
    
  }, error = function(e) {
    warning("Strict cleaning failed, falling back to moderate: ", e$message)
    clean_moderate(data)
  })
}

# Function to assess quality of remaining records
assess_record_quality <- function(data) {
  quality_flags <- rep("good", nrow(data))
  
  # Flag records with high uncertainty
  if ("coordinateUncertaintyInMeters" %in% names(data)) {
    high_uncertainty <- !is.na(data$coordinateUncertaintyInMeters) & 
                       data$coordinateUncertaintyInMeters > 50000
    quality_flags[high_uncertainty] <- "uncertain"
  }
  
  # Flag old records  
  if ("year" %in% names(data)) {
    old_records <- !is.na(data$year) & data$year < 2000
    quality_flags[old_records & quality_flags == "good"] <- "old"
  }
  
  return(quality_flags)
}