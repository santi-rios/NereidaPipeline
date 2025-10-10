# Geospatial Data Processing Functions using geotargets and terra
# For handling marine occurrence data with spatial coordinates and environmental layers

#' Create spatial raster data for marine environments
#' @param occurrence_data Data frame with decimalLatitude and decimalLongitude columns
#' @param buffer_km Numeric. Buffer distance in kilometers around occurrences
#' @param resolution Numeric. Raster resolution in degrees
#' @return SpatRaster object with marine occurrence density
create_occurrence_raster <- function(occurrence_data, 
                                   buffer_km = 50, 
                                   resolution = 0.1) {
  
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for spatial data processing.")
  }
  
  # Clean coordinates
  clean_coords <- occurrence_data[
    !is.na(occurrence_data$decimalLatitude) & 
    !is.na(occurrence_data$decimalLongitude) &
    occurrence_data$decimalLatitude >= -90 & 
    occurrence_data$decimalLatitude <= 90 &
    occurrence_data$decimalLongitude >= -180 & 
    occurrence_data$decimalLongitude <= 180,
  ]
  
  if (nrow(clean_coords) == 0) {
    stop("No valid coordinates found in occurrence data")
  }
  
  # Create spatial extent
  lon_range <- range(clean_coords$decimalLongitude, na.rm = TRUE)
  lat_range <- range(clean_coords$decimalLatitude, na.rm = TRUE)
  
  # Add buffer to extent
  buffer_deg <- buffer_km / 111  # Approximate km to degrees conversion
  extent_vec <- c(
    lon_range[1] - buffer_deg, lon_range[2] + buffer_deg,
    lat_range[1] - buffer_deg, lat_range[2] + buffer_deg
  )
  
  # Create raster template
  raster_template <- terra::rast(
    extent = terra::ext(extent_vec),
    resolution = resolution,
    crs = "EPSG:4326"
  )
  
  # Create points vector
  coords_vect <- terra::vect(
    clean_coords, 
    geom = c("decimalLongitude", "decimalLatitude"),
    crs = "EPSG:4326"
  )
  
  # Rasterize occurrence points (density)
  occurrence_raster <- terra::rasterize(
    coords_vect, 
    raster_template, 
    fun = "count",
    background = 0
  )
  
  names(occurrence_raster) <- "occurrence_density"
  
  return(occurrence_raster)
}

#' Download and process marine environmental layers
#' @param extent SpatExtent or numeric vector (xmin, xmax, ymin, ymax)
#' @param variables Character vector of environmental variables
#' @param output_dir Character. Directory to save environmental layers
#' @return SpatRasterCollection with environmental layers
get_marine_environmental_layers <- function(extent, 
                                          variables = c("sst", "salinity", "depth", "chlorophyll"),
                                          output_dir = "./data/raw/environmental") {
  
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for environmental data processing.")
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  env_layers <- list()
  saved_files <- c() # <-- AÑADIR ESTO
  
  # Bio-ORACLE or other marine environmental data sources
  # This is a simplified example - in practice you would download from specific sources
  
  tryCatch({
    for (var in variables) {
      if (var == "depth") {
        # Create example bathymetry layer (replace with actual data source)
        depth_raster <- terra::rast(
          extent = extent,
          resolution = 0.1,
          crs = "EPSG:4326"
        )
        
        # Simulate depth data (negative values for ocean depth)
        terra::values(depth_raster) <- runif(terra::ncell(depth_raster), -5000, 0)
        names(depth_raster) <- "depth_m"
        
        env_layers[["depth"]] <- depth_raster
        
        # Save layer
        depth_file <- file.path(output_dir, "depth.tif")
        terra::writeRaster(depth_raster, depth_file, overwrite = TRUE)
        saved_files <- c(saved_files, depth_file) # <-- AÑADIR ESTO
        
      } else if (var == "sst") {
        # Create example sea surface temperature layer
        sst_raster <- terra::rast(
          extent = extent,
          resolution = 0.1,
          crs = "EPSG:4326"
        )
        
        # Simulate SST data (in Celsius)
        terra::values(sst_raster) <- runif(terra::ncell(sst_raster), 15, 30)
        names(sst_raster) <- "sst_celsius"
        
        env_layers[["sst"]] <- sst_raster
        
        # Save layer
        sst_file <- file.path(output_dir, "sst.tif")
        terra::writeRaster(sst_raster, sst_file, overwrite = TRUE)
        saved_files <- c(saved_files, sst_file) # <-- AÑADIR ESTO
        
      } else if (var == "salinity") {
        # Create example salinity layer
        sal_raster <- terra::rast(
          extent = extent,
          resolution = 0.1,
          crs = "EPSG:4326"
        )
        
        # Simulate salinity data (PSU)
        terra::values(sal_raster) <- runif(terra::ncell(sal_raster), 32, 37)
        names(sal_raster) <- "salinity_psu"
        
        env_layers[["salinity"]] <- sal_raster
        
        # Save layer
        sal_file <- file.path(output_dir, "salinity.tif")
        terra::writeRaster(sal_raster, sal_file, overwrite = TRUE)
        saved_files <- c(saved_files, sal_file) # <-- AÑADIR ESTO
        
      } else if (var == "chlorophyll") {
        # Create example chlorophyll-a layer
        chl_raster <- terra::rast(
          extent = extent,
          resolution = 0.1,
          crs = "EPSG:4326"
        )
        
        # Simulate chlorophyll-a data (mg/m³)
        terra::values(chl_raster) <- runif(terra::ncell(chl_raster), 0.1, 10)
        names(chl_raster) <- "chlorophyll_mg_m3"
        
        env_layers[["chlorophyll"]] <- chl_raster
        
        # Save layer
        chl_file <- file.path(output_dir, "chlorophyll.tif")
        terra::writeRaster(chl_raster, chl_file, overwrite = TRUE)
        saved_files <- c(saved_files, chl_file) # <-- AÑADIR ESTO
      }
    }
    
    return(saved_files) # <-- CAMBIAR VALOR DE RETORNO
    
  }, error = function(e) {
    warning("Environmental layer processing failed: ", e$message)
    return(NULL)
  })
}

#' Extract environmental values at occurrence points
#' @param occurrence_data Data frame with coordinate columns
#' @param env_layers SpatRasterCollection or list of SpatRaster objects
#' @return Data frame with occurrence data and environmental values
extract_environmental_values <- function(occurrence_data, env_layers) {
  
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for environmental data extraction.")
  }
  
  # Clean coordinates
  clean_coords <- occurrence_data[
    !is.na(occurrence_data$decimalLatitude) & 
    !is.na(occurrence_data$decimalLongitude),
  ]
  
  if (nrow(clean_coords) == 0) {
    warning("No valid coordinates found")
    return(occurrence_data)
  }
  
  # Create coordinate matrix
  coords_matrix <- cbind(
    clean_coords$decimalLongitude,
    clean_coords$decimalLatitude
  )
  
  # Extract values from environmental layers
  env_values <- list()
  
  if (inherits(env_layers, "SpatRaster")) { # <-- CAMBIAR ESTA LÓGICA
    # Handle single multi-layer raster
    values <- terra::extract(env_layers, coords_matrix, method = "bilinear")
    env_values <- values[, -1, drop = FALSE]  # Remove ID column
  } else if (inherits(env_layers, "SpatRasterCollection")) {
    # Handle raster collection
    for (i in seq_along(env_layers)) {
      layer <- env_layers[[i]]
      values <- terra::extract(layer, coords_matrix, method = "bilinear")
      env_values <- c(env_values, values[, -1, drop = FALSE])  # Remove ID column
    }
  }
  
  # Combine with occurrence data
  if (length(env_values) > 0) {
    env_df <- as.data.frame(env_values)
    
    # Match indices
    result_data <- occurrence_data
    result_data[seq_len(nrow(clean_coords)), names(env_df)] <- env_df
    
    return(result_data)
  }
  
  return(occurrence_data)
}

#' Create marine habitat suitability model
#' @param occurrence_data Data frame with species occurrences and environmental data
#' @param env_layers SpatRasterCollection with environmental predictors
#' @param species_column Character. Name of species column
#' @return List with model object and prediction raster
create_habitat_suitability_model <- function(occurrence_data, 
                                           env_layers, 
                                           species_column = "scientificName") {
  
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required for habitat modeling.")
  }
  
  # Prepare data for modeling
  model_data <- occurrence_data[!is.na(occurrence_data[[species_column]]), ]
  
  # Get environmental variable names
  env_names <- names(env_layers)
  
  # Check if environmental data is available
  available_env <- intersect(env_names, names(model_data))
  
  if (length(available_env) == 0) {
    stop("No environmental variables found in occurrence data")
  }
  
  # Create presence/absence data (simplified example)
  model_data$presence <- 1  # All occurrences are presences
  
  # Generate background/pseudo-absence points
  extent_coords <- terra::ext(env_layers[[1]])
  
  # Sample random background points
  n_background <- min(nrow(model_data) * 2, 1000)
  background_coords <- data.frame(
    decimalLongitude = runif(n_background, extent_coords[1], extent_coords[2]),
    decimalLatitude = runif(n_background, extent_coords[3], extent_coords[4])
  )
  
  # Extract environmental values for background points
  background_env <- extract_environmental_values(background_coords, env_layers)
  background_env$presence <- 0
  background_env[[species_column]] <- "background"
  
  # Combine presence and background data
  full_model_data <- rbind(
    model_data[, c(species_column, "decimalLongitude", "decimalLatitude", 
                  available_env, "presence")],
    background_env[, c(species_column, "decimalLongitude", "decimalLatitude", 
                      available_env, "presence")]
  )
  
  # Remove rows with missing environmental data
  complete_cases <- complete.cases(full_model_data[, available_env])
  full_model_data <- full_model_data[complete_cases, ]
  
  if (nrow(full_model_data) < 20) {
    warning("Insufficient data for modeling (n = ", nrow(full_model_data), ")")
    return(NULL)
  }
  
  # Fit GLM model (simple example - could use MaxEnt, Random Forest, etc.)
  model_formula <- as.formula(paste("presence ~", paste(available_env, collapse = " + ")))
  
  habitat_model <- tryCatch({
    glm(model_formula, data = full_model_data, family = "binomial")
  }, error = function(e) {
    warning("Model fitting failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(habitat_model)) {
    return(NULL)
  }
  
  # Create prediction raster
  prediction_raster <- tryCatch({
    terra::predict(env_layers, habitat_model, type = "response")
  }, error = function(e) {
    warning("Prediction failed: ", e$message)
    return(NULL)
  })
  
  names(prediction_raster) <- paste0(gsub(" ", "_", species_column), "_habitat_suitability")
  
  return(list(
    model = habitat_model,
    prediction = prediction_raster,
    model_data = full_model_data,
    formula = model_formula,
    env_variables = available_env
  ))
}

#' Create spatial analysis summary for marine species
#' @param occurrence_data Data frame with occurrence records
#' @param env_layers SpatRasterCollection with environmental layers
#' @param output_dir Character. Output directory for spatial products
#' @return List with spatial analysis results
analyze_marine_spatial_patterns <- function(occurrence_data, 
                                          env_layers, 
                                          output_dir = "./data/processed/spatial") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  spatial_results <- list()
  
  # 1. Create occurrence density raster
  message("Creating occurrence density raster...")
  occurrence_raster <- create_occurrence_raster(occurrence_data)
  spatial_results$occurrence_density <- occurrence_raster
  
  # Save occurrence raster
  terra::writeRaster(
    occurrence_raster, 
    file.path(output_dir, "occurrence_density.tif"), 
    overwrite = TRUE
  )
  
  # 2. Extract environmental values
  message("Extracting environmental values...")
  occurrence_with_env <- extract_environmental_values(occurrence_data, env_layers)
  spatial_results$occurrence_environmental <- occurrence_with_env
  
  # Save environmental data
  utils::write.csv(
    occurrence_with_env, 
    file.path(output_dir, "occurrence_environmental_data.csv"), 
    row.names = FALSE
  )
  
  # 3. Calculate spatial statistics
  message("Calculating spatial statistics...")
  
  # Extent of occurrences
  valid_coords <- occurrence_data[
    !is.na(occurrence_data$decimalLatitude) & 
    !is.na(occurrence_data$decimalLongitude),
  ]
  
  if (nrow(valid_coords) > 0) {
    spatial_stats <- list(
      n_occurrences = nrow(valid_coords),
      lat_range = range(valid_coords$decimalLatitude),
      lon_range = range(valid_coords$decimalLongitude),
      centroid_lat = mean(valid_coords$decimalLatitude),
      centroid_lon = mean(valid_coords$decimalLongitude)
    )
    
    # Calculate area of occupancy (AOO) - simplified 2x2 km grid
    coords_matrix <- cbind(valid_coords$decimalLongitude, valid_coords$decimalLatitude)
    
    # Create 2x2 km grid cells
    grid_size <- 2 / 111  # 2 km in degrees (approximate)
    
    # Round coordinates to grid
    grid_coords <- data.frame(
      grid_lon = round(coords_matrix[, 1] / grid_size) * grid_size,
      grid_lat = round(coords_matrix[, 2] / grid_size) * grid_size
    )
    
    # Count unique grid cells
    unique_cells <- unique(grid_coords)
    aoo_km2 <- nrow(unique_cells) * 4  # 2x2 km cells
    
    spatial_stats$area_of_occupancy_km2 <- aoo_km2
    spatial_stats$n_grid_cells <- nrow(unique_cells)
    
    spatial_results$spatial_statistics <- spatial_stats
  }
  
  # 4. Environmental niche analysis
  if (!is.null(occurrence_with_env)) {
    env_names <- names(env_layers)
    available_env <- intersect(env_names, names(occurrence_with_env))
    
    if (length(available_env) > 0) {
      env_summary <- list()
      
      for (var in available_env) {
        var_data <- occurrence_with_env[[var]]
        var_data <- var_data[!is.na(var_data)]
        
        if (length(var_data) > 0) {
          env_summary[[var]] <- list(
            mean = mean(var_data),
            sd = sd(var_data),
            min = min(var_data),
            max = max(var_data),
            median = median(var_data),
            q25 = quantile(var_data, 0.25),
            q75 = quantile(var_data, 0.75)
          )
        }
      }
      
      spatial_results$environmental_niche <- env_summary
    }
  }
  
  # Save spatial analysis summary
  jsonlite::write_json(
    spatial_results[!names(spatial_results) %in% c("occurrence_density", "occurrence_environmental")],
    file.path(output_dir, "spatial_analysis_summary.json"),
    pretty = TRUE
  )
  
  message("Spatial analysis completed.")
  return(spatial_results)
}