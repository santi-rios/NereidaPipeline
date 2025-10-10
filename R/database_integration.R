# NoSQL Database Integration using nodbi with DuckDB backend
# For efficient storage and querying of heterogeneous marine biodiversity data

#' Initialize DuckDB database connection for marine biodiversity data
#' @param db_path Character. Path to DuckDB database file (use ":memory:" for in-memory)
#' @param create_schema Logical. Whether to create predefined schema for marine data
#' @return nodbi source object for DuckDB connection
init_marine_database <- function(db_path = "./data/marine_biodiversity.duckdb",
                                create_schema = TRUE) {
  
  if (!requireNamespace("nodbi", quietly = TRUE)) {
    stop("Package 'nodbi' is required for NoSQL database operations.")
  }
  
  if (!requireNamespace("DBI", quietly = TRUE)) {
    stop("Package 'DBI' is required for database operations.")
  }
  
  # Create database directory if it doesn't exist
  if (db_path != ":memory:") {
    db_dir <- dirname(db_path)
    if (!dir.exists(db_dir)) {
      dir.create(db_dir, recursive = TRUE)
    }
  }
  
  # Initialize DuckDB connection
  src <- tryCatch({
    nodbi::src_duckdb(dbdir = db_path)
  }, error = function(e) {
    stop("Failed to initialize DuckDB connection: ", e$message)
  })
  
  if (create_schema) {
    # Create predefined containers (tables) for marine data types
    marine_containers <- c(
      "occurrences", "taxonomic_data", "environmental_data",
      "genomic_metadata", "cleaning_logs", "spatial_analysis"
    )
    
    for (container in marine_containers) {
      if (!nodbi::docdb_exists(src, container)) {
        # Create empty container with sample structure
        sample_doc <- list(
          created_at = Sys.time(),
          data_type = container,
          description = paste("Container for", container, "in marine biodiversity research")
        )
        
        tryCatch({
          nodbi::docdb_create(src, container, sample_doc)
          message("Created container: ", container)
        }, error = function(e) {
          warning("Failed to create container ", container, ": ", e$message)
        })
      }
    }
  }
  
  return(src)
}

#' Store occurrence data in DuckDB
#' @param src nodbi source object
#' @param occurrence_data Data frame with occurrence records
#' @param container_name Character. Name of the container to store data
#' @param metadata List. Additional metadata to attach to records
#' @return Number of records stored
store_occurrence_data <- function(src, 
                                 occurrence_data, 
                                 container_name = "occurrences",
                                 metadata = NULL) {
  
  if (!inherits(src, "src_duckdb")) {
    stop("Source must be a DuckDB nodbi connection")
  }
  
  if (nrow(occurrence_data) == 0) {
    warning("No data to store")
    return(0)
  }
  
  # Prepare data for storage
  storage_data <- occurrence_data
  
  # Add metadata
  storage_data$stored_at <- Sys.time()
  storage_data$data_type <- "occurrence"
  
  if (!is.null(metadata)) {
    for (name in names(metadata)) {
      storage_data[[name]] <- metadata[[name]]
    }
  }
  
  # Convert to list of records for nodbi
  record_list <- split(storage_data, seq_len(nrow(storage_data)))
  names(record_list) <- NULL
  
  # Store in database
  n_stored <- 0
  
  for (record in record_list) {
    tryCatch({
      nodbi::docdb_create(src, container_name, record)
      n_stored <- n_stored + 1
    }, error = function(e) {
      warning("Failed to store record: ", e$message)
    })
  }
  
  message("Stored ", n_stored, " occurrence records in container: ", container_name)
  return(n_stored)
}

#' Query occurrence data from DuckDB
#' @param src nodbi source object
#' @param container_name Character. Name of container to query
#' @param species Character. Species name to filter (optional)
#' @param date_range Character vector. Date range [start, end] (optional)
#' @param spatial_bounds List. Spatial bounding box with lat/lon limits (optional)
#' @return Data frame with query results
query_occurrence_data <- function(src,
                                 container_name = "occurrences", 
                                 species = NULL,
                                 date_range = NULL,
                                 spatial_bounds = NULL) {
  
  if (!inherits(src, "src_duckdb")) {
    stop("Source must be a DuckDB nodbi connection")
  }
  
  # Build query
  query_conditions <- list()
  
  # Species filter
  if (!is.null(species)) {
    query_conditions$scientificName <- list("$regex" = species)
  }
  
  # Date range filter
  if (!is.null(date_range) && length(date_range) == 2) {
    query_conditions$eventDate <- list(
      "$gte" = date_range[1],
      "$lte" = date_range[2]
    )
  }
  
  # Spatial bounds filter
  if (!is.null(spatial_bounds)) {
    if (all(c("min_lat", "max_lat", "min_lon", "max_lon") %in% names(spatial_bounds))) {
      query_conditions$decimalLatitude <- list(
        "$gte" = spatial_bounds$min_lat,
        "$lte" = spatial_bounds$max_lat
      )
      query_conditions$decimalLongitude <- list(
        "$gte" = spatial_bounds$min_lon,
        "$lte" = spatial_bounds$max_lon
      )
    }
  }
  
  # Execute query
  if (length(query_conditions) > 0) {
    query_json <- jsonlite::toJSON(query_conditions, auto_unbox = TRUE)
    results <- nodbi::docdb_query(src, container_name, query = query_json)
  } else {
    # Get all records if no conditions
    results <- nodbi::docdb_get(src, container_name)
  }
  
  return(results)
}

#' Store taxonomic data in DuckDB
#' @param src nodbi source object
#' @param taxonomic_data List or data frame with taxonomic information
#' @param container_name Character. Name of container for taxonomic data
#' @return Number of taxonomic records stored
store_taxonomic_data <- function(src,
                                taxonomic_data,
                                container_name = "taxonomic_data") {
  
  if (!inherits(src, "src_duckdb")) {
    stop("Source must be a DuckDB nodbi connection")
  }
  
  # Handle different input types
  if (is.data.frame(taxonomic_data)) {
    # Convert data frame to list of records
    taxonomic_list <- split(taxonomic_data, seq_len(nrow(taxonomic_data)))
    names(taxonomic_list) <- NULL
  } else if (is.list(taxonomic_data)) {
    taxonomic_list <- list(taxonomic_data)
  } else {
    stop("taxonomic_data must be a data frame or list")
  }
  
  # Add metadata to each record
  for (i in seq_along(taxonomic_list)) {
    taxonomic_list[[i]]$stored_at <- Sys.time()
    taxonomic_list[[i]]$data_type <- "taxonomy"
  }
  
  # Store records
  n_stored <- 0
  
  for (record in taxonomic_list) {
    tryCatch({
      nodbi::docdb_create(src, container_name, record)
      n_stored <- n_stored + 1
    }, error = function(e) {
      warning("Failed to store taxonomic record: ", e$message)
    })
  }
  
  message("Stored ", n_stored, " taxonomic records in container: ", container_name)
  return(n_stored)
}

#' Store environmental data in DuckDB
#' @param src nodbi source object
#' @param environmental_data Data frame with environmental measurements
#' @param container_name Character. Name of container for environmental data
#' @param spatial_metadata List. Additional spatial metadata
#' @return Number of environmental records stored
store_environmental_data <- function(src,
                                   environmental_data,
                                   container_name = "environmental_data",
                                   spatial_metadata = NULL) {
  
  if (!inherits(src, "src_duckdb")) {
    stop("Source must be a DuckDB nodbi connection")
  }
  
  if (nrow(environmental_data) == 0) {
    warning("No environmental data to store")
    return(0)
  }
  
  # Prepare environmental data for storage
  env_storage_data <- environmental_data
  env_storage_data$stored_at <- Sys.time()
  env_storage_data$data_type <- "environmental"
  
  # Add spatial metadata if provided
  if (!is.null(spatial_metadata)) {
    for (name in names(spatial_metadata)) {
      env_storage_data[[paste0("spatial_", name)]] <- spatial_metadata[[name]]
    }
  }
  
  # Convert to list of records
  record_list <- split(env_storage_data, seq_len(nrow(env_storage_data)))
  names(record_list) <- NULL
  
  # Store records
  n_stored <- 0
  
  for (record in record_list) {
    tryCatch({
      nodbi::docdb_create(src, container_name, record)
      n_stored <- n_stored + 1
    }, error = function(e) {
      warning("Failed to store environmental record: ", e$message)
    })
  }
  
  message("Stored ", n_stored, " environmental records in container: ", container_name)
  return(n_stored)
}

#' Create aggregated summaries from stored data
#' @param src nodbi source object
#' @param summary_type Character. Type of summary ("species", "spatial", "temporal")
#' @param container_name Character. Container to summarize
#' @return Data frame with summary statistics
create_data_summary <- function(src,
                               summary_type = "species",
                               container_name = "occurrences") {
  
  if (!inherits(src, "src_duckdb")) {
    stop("Source must be a DuckDB nodbi connection")
  }
  
  # Get all data from container
  all_data <- nodbi::docdb_get(src, container_name)
  
  if (is.null(all_data) || nrow(all_data) == 0) {
    warning("No data found in container: ", container_name)
    return(data.frame())
  }
  
  summary_result <- switch(summary_type,
    "species" = {
      # Species-level summaries
      if ("scientificName" %in% names(all_data)) {
        # Create species summary using base R approach
        species_counts <- aggregate(
          cbind(n_records = rep(1, nrow(all_data)), 
                lat = all_data$decimalLatitude, 
                lon = all_data$decimalLongitude),
          by = list(species = all_data$scientificName),
          FUN = function(x) if(is.numeric(x)) sum(x, na.rm = TRUE) else length(x),
          na.rm = TRUE
        )
        species_summary <- data.frame(
          scientificName = species_counts$species,
          n_records = species_counts$n_records,
          stringsAsFactors = FALSE
        )
        species_summary
      } else {
        data.frame(message = "No species information available")
      }
    },
    
    "spatial" = {
      # Spatial summaries
      if (all(c("decimalLatitude", "decimalLongitude") %in% names(all_data))) {
        # Create spatial grid summary (1-degree cells)
        all_data$grid_lat <- floor(all_data$decimalLatitude)
        all_data$grid_lon <- floor(all_data$decimalLongitude)
        
        # Create spatial summary using base R
        spatial_counts <- aggregate(
          cbind(n_records = rep(1, nrow(all_data))),
          by = list(grid_lat = all_data$grid_lat, grid_lon = all_data$grid_lon),
          FUN = sum,
          na.rm = TRUE
        )
        spatial_summary <- spatial_counts
        spatial_summary
      } else {
        data.frame(message = "No spatial information available")
      }
    },
    
    "temporal" = {
      # Temporal summaries
      if ("year" %in% names(all_data)) {
        # Create temporal summary using base R
        temporal_counts <- aggregate(
          cbind(n_records = rep(1, nrow(all_data))),
          by = list(year = all_data$year),
          FUN = sum,
          na.rm = TRUE
        )
        temporal_summary <- temporal_counts[order(temporal_counts$year), ]
        temporal_summary
      } else {
        data.frame(message = "No temporal information available")
      }
    },
    
    # Default case
    data.frame(message = paste("Unknown summary type:", summary_type))
  )
  
  return(summary_result)
}

#' Export data from DuckDB to various formats
#' @param src nodbi source object
#' @param container_name Character. Container to export
#' @param output_format Character. Output format ("csv", "json", "parquet")
#' @param output_path Character. Output file path
#' @param query Character. Optional query to filter data
#' @return Logical indicating success
export_marine_data <- function(src,
                              container_name,
                              output_format = "csv",
                              output_path,
                              query = NULL) {
  
  if (!inherits(src, "src_duckdb")) {
    stop("Source must be a DuckDB nodbi connection")
  }
  
  # Get data from container
  if (!is.null(query)) {
    data_to_export <- nodbi::docdb_query(src, container_name, query = query)
  } else {
    data_to_export <- nodbi::docdb_get(src, container_name)
  }
  
  if (is.null(data_to_export) || nrow(data_to_export) == 0) {
    warning("No data to export from container: ", container_name)
    return(FALSE)
  }
  
  # Create output directory if needed
  output_dir <- dirname(output_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Export based on format
  success <- tryCatch({
    switch(output_format,
      "csv" = {
        utils::write.csv(data_to_export, output_path, row.names = FALSE)
        TRUE
      },
      
      "json" = {
        jsonlite::write_json(data_to_export, output_path, pretty = TRUE)
        TRUE
      },
      
      "parquet" = {
        if (requireNamespace("arrow", quietly = TRUE)) {
          arrow::write_parquet(data_to_export, output_path)
          TRUE
        } else {
          warning("Package 'arrow' required for Parquet export")
          FALSE
        }
      },
      
      {
        warning("Unsupported output format: ", output_format)
        FALSE
      }
    )
  }, error = function(e) {
    warning("Export failed: ", e$message)
    FALSE
  })
  
  if (success) {
    message("Successfully exported ", nrow(data_to_export), " records to: ", output_path)
  }
  
  return(success)
}

#' Initialize DuckDB connection with required extensions
#' @param db_path Character. Path to the DuckDB database file
#' @return DuckDB connection object
initialize_duckdb_connection <- function(db_path) {
  
  if (!requireNamespace("duckdb", quietly = TRUE)) {
    stop("Package 'duckdb' is required for database integration.")
  }
  
  if (!requireNamespace("DBI", quietly = TRUE)) {
    stop("Package 'DBI' is required for database integration.")
  }
  
  tryCatch({
    # Create connection
    con <- duckdb::dbConnect(duckdb::duckdb(), dbdir = db_path)
    
    # Install and load required extensions
    message("Installing DuckDB extensions...")
    DBI::dbExecute(con, "INSTALL json;")
    DBI::dbExecute(con, "LOAD json;")
    
    # Optionally install other useful extensions
    tryCatch({
      DBI::dbExecute(con, "INSTALL spatial;")
      DBI::dbExecute(con, "LOAD spatial;")
    }, error = function(e) {
      message("Spatial extension installation failed (optional): ", e$message)
    })
    
    message("DuckDB connection initialized successfully")
    return(con)
    
  }, error = function(e) {
    stop("Failed to initialize DuckDB connection: ", e$message)
  })
}

#' Integrate marine biodiversity data into DuckDB database
#' @param occurrence_data Data frame with occurrence records
#' @param taxonomic_data List with taxonomic information  
#' @param environmental_data Data frame with environmental data (optional)
#' @param db_path Character. Path to DuckDB database file
#' @return nodbi connection object (not raw DBI connection)
integrate_marine_data <- function(occurrence_data, 
                                 taxonomic_data = NULL,
                                 environmental_data = NULL,
                                 db_path = "data/marine_biodiversity.duckdb") {
  
  message("Starting marine data integration workflow...")
  
  # Ensure output directory exists
  db_dir <- dirname(db_path)
  if (!dir.exists(db_dir)) {
    dir.create(db_dir, recursive = TRUE)
  }
  
  # Initialize connection with extensions (DBI connection)
  dbi_con <- initialize_duckdb_connection(db_path)
  
  # Close the DBI connection as we'll create a nodbi connection
  DBI::dbDisconnect(dbi_con)
  
  # Create nodbi connection for data operations
  src <- nodbi::src_duckdb(dbdir = db_path)
  
  # Create or update occurrences table
  if (!is.null(occurrence_data) && nrow(occurrence_data) > 0) {
    message("Storing occurrence data...")
    
    # Convert data to JSON format for nodbi
    occurrence_json <- occurrence_data %>%
      dplyr::mutate(
        id = row_number(),
        latitude = as.numeric(decimalLatitude),
        longitude = as.numeric(decimalLongitude),
        species = scientificName,
        source = ifelse(is.na(datasetName), "unknown", datasetName),
        date_collected = eventDate,
        created_at = Sys.time()
      ) %>%
      dplyr::select(-decimalLatitude, -decimalLongitude, -scientificName) %>%
      jsonlite::toJSON(dataframe = "rows")
    
    # Store in nodbi
    tryCatch({
      nodbi::docdb_create(src, "occurrences", occurrence_json)
      message("Occurrence data stored successfully")
    }, error = function(e) {
      # If collection exists, update it
      tryCatch({
        nodbi::docdb_update(src, "occurrences", occurrence_json)
        message("Occurrence data updated successfully")
      }, error = function(e2) {
        warning("Failed to store occurrence data: ", e2$message)
      })
    })
  }
  
  # Store taxonomic data if available
  if (!is.null(taxonomic_data) && length(taxonomic_data) > 0) {
    message("Storing taxonomic data...")
    
    taxonomic_json <- jsonlite::toJSON(taxonomic_data, auto_unbox = TRUE)
    
    tryCatch({
      nodbi::docdb_create(src, "taxonomy", taxonomic_json)
      message("Taxonomic data stored successfully")
    }, error = function(e) {
      tryCatch({
        nodbi::docdb_update(src, "taxonomy", taxonomic_json)
        message("Taxonomic data updated successfully")
      }, error = function(e2) {
        warning("Failed to store taxonomic data: ", e2$message)
      })
    })
  }
  
  # Store environmental data if available
  if (!is.null(environmental_data) && nrow(environmental_data) > 0) {
    message("Storing environmental data...")
    
    env_json <- environmental_data %>%
      dplyr::mutate(
        id = row_number(),
        created_at = Sys.time()
      ) %>%
      jsonlite::toJSON(dataframe = "rows")
    
    tryCatch({
      nodbi::docdb_create(src, "environmental", env_json)
      message("Environmental data stored successfully")
    }, error = function(e) {
      tryCatch({
        nodbi::docdb_update(src, "environmental", env_json)
        message("Environmental data updated successfully")
      }, error = function(e2) {
        warning("Failed to store environmental data: ", e2$message)
      })
    })
  }
  
  message("Marine data integration completed")
  return(src) # <-- Devolver la conexión nodbi, no la DBI
}