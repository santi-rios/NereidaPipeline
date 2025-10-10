# Enhanced Data Acquisition Functions for Marine Biodiversity Research
# Integrating OBIS, GBIF, biomaRt, wikitaxa, and PRISM for comprehensive data collection

#' Get occurrence data from OBIS with enhanced error handling
#' @param species Character. Scientific name of the species
#' @param out_dir Character. Output directory path
#' @param suffix Character. File suffix for output
#' @param write Logical. Whether to write to file
#' @return Data frame with occurrence records
get_obis_data <- function(species,
              out_dir = "./data/raw",
              suffix = "_obis_data.csv",
              write = TRUE) {
  if (missing(species) || !nzchar(as.character(species)[1])) {
    stop("Please provide a non-empty species name (character string).")
  }
  species <- as.character(species)[1]

  # sanitize species for a filename: lowercase, keep alnum, replace others with _
  safe_name <- tolower(gsub("[^a-z0-9]+", "_", species))
  safe_name <- gsub("^_|_$", "", safe_name)  # trim leading/trailing underscores
  filename <- file.path(out_dir, paste0(safe_name, suffix))

  if (write && !dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  res <- tryCatch(
    robis::occurrence(scientificname = species),
    error = function(e) {
      warning("OBIS request failed: ", e$message)
      return(NULL)
    }
  )

  if (is.null(res) || nrow(res) == 0) {
    warning("No data returned; nothing written to disk.")
    return(NULL)
  }

  # Add metadata
  res$data_source <- "OBIS"
  res$query_date <- Sys.Date()
  res$query_species <- species

  if (write) {
    utils::write.csv(res, filename, row.names = FALSE)
  }

  invisible(res)
}

#' Get occurrence data from GBIF for marine species
#' @param species Character. Scientific name of the species
#' @param out_dir Character. Output directory path
#' @param limit Integer. Maximum number of records to retrieve
#' @param marine_only Logical. Filter for marine environments only
#' @return Data frame with occurrence records
get_gbif_marine_data <- function(species,
                                out_dir = "./data/raw",
                                limit = 10000,
                                marine_only = TRUE) {
  
  safe_name <- tolower(gsub("[^a-z0-9]+", "_", species))
  safe_name <- gsub("^_|_$", "", safe_name)
  filename <- file.path(out_dir, paste0(safe_name, "_gbif_data.csv"))
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  res <- tryCatch({
    # Get occurrence data with marine habitat filters
    data <- rgbif::occ_search(
      scientificName = species,
      limit = limit,
      hasCoordinate = TRUE,
      basisOfRecord = c("HUMAN_OBSERVATION", "OBSERVATION", "PRESERVED_SPECIMEN"),
      year = paste0(1950, ",", format(Sys.Date(), "%Y")) # Recent records only
    )
    
    if (marine_only && !is.null(data$data)) {
      # Filter for marine environments (depth > 0 or marine keywords)
      marine_data <- data$data[
        (!is.na(data$data$depth) & data$data$depth > 0) |
        grepl("marine|ocean|sea|reef|pelagic", 
              paste(data$data$habitat, data$data$locality, data$data$waterBody), 
              ignore.case = TRUE) |
        (!is.na(data$data$decimalLatitude) & !is.na(data$data$decimalLongitude)),
      ]
      data$data <- marine_data
    }
    
    data$data
  }, error = function(e) {
    warning("GBIF request failed: ", e$message)
    return(NULL)
  })
  
  if (!is.null(res) && nrow(res) > 0) {
    res$data_source <- "GBIF"
    res$query_date <- Sys.Date()
    res$query_species <- species
    
    utils::write.csv(res, filename, row.names = FALSE)
  }
  
  invisible(res)
}

#' Get genomic/sequence data using biomaRt
#' @param species Character. Scientific name of the species
#' @param db Character. Database to search ("refseq", "genbank", "ensembl")
#' @param type Character. Type of data ("genome", "proteome", "CDS", "RNA")
#' @param out_dir Character. Output directory
#' @return Character vector of downloaded file paths
get_sequence_data <- function(species,
                             db = "refseq",
                             type = "genome",
                             out_dir = "./data/raw/sequences") {
  
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package 'biomaRt' is required for sequence retrieval.")
  }
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  safe_name <- tolower(gsub("[^a-z0-9]+", "_", species))
  
  tryCatch({
    # Check if organism is available in the database
    available_organisms <- biomaRt::listGenomes(db = db, type = type)
    
    # Search for the species in available organisms
    species_match <- grep(species, available_organisms$organism_name, 
                         ignore.case = TRUE, value = FALSE)
    
    if (length(species_match) == 0) {
      # Try genus search if species not found
      genus <- strsplit(species, " ")[[1]][1]
      genus_match <- grep(genus, available_organisms$organism_name, 
                         ignore.case = TRUE, value = FALSE)
      
      if (length(genus_match) == 0) {
        warning("No genomic data found for ", species, " in ", db)
        return(NULL)
      }
      species_match <- genus_match[1]  # Take first genus match
    } else {
      species_match <- species_match[1]  # Take first species match
    }
    
    organism_name <- available_organisms$organism_name[species_match]
    
    # Download the data
    if (type == "genome") {
      file_path <- biomaRt::getGenome(
        db = db,
        organism = organism_name,
        path = out_dir
      )
    } else if (type == "proteome") {
      file_path <- biomaRt::getProteome(
        db = db,
        organism = organism_name,
        path = out_dir
      )
    } else if (type == "CDS") {
      file_path <- biomaRt::getCDS(
        db = db,
        organism = organism_name,
        path = out_dir
      )
    }
    
    return(file_path)
    
  }, error = function(e) {
    warning("Sequence retrieval failed for ", species, ": ", e$message)
    return(NULL)
  })
}

#' Get taxonomic information using wikitaxa
#' @param species Character. Scientific name of the species
#' @param out_dir Character. Output directory
#' @return List with taxonomic information
get_taxonomic_data <- function(species, out_dir = "./data/raw/taxonomy") {
  
  if (!requireNamespace("wikitaxa", quietly = TRUE)) {
    stop("Package 'wikitaxa' is required for taxonomic data retrieval.")
  }
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  safe_name <- tolower(gsub("[^a-z0-9]+", "_", species))
  filename <- file.path(out_dir, paste0(safe_name, "_taxonomy.json"))
  
  tryCatch({
    # Get taxonomic information from Wikipedia
    wiki_data <- wikitaxa::wt_wikispecies(species)
    
    # Get additional taxonomic classification
    wiki_summary <- wikitaxa::wt_wikipedia(species)
    
    # Combine results
    taxonomic_info <- list(
      species = species,
      query_date = Sys.Date(),
      wikispecies_data = wiki_data,
      wikipedia_summary = wiki_summary
    )
    
    # Save as JSON
    jsonlite::write_json(taxonomic_info, filename, pretty = TRUE)
    
    return(taxonomic_info)
    
  }, error = function(e) {
    warning("Taxonomic data retrieval failed for ", species, ": ", e$message)
    return(NULL)
  })
}

#' Get environmental/climate data using PRISM
#' @param coordinates Data frame with lat/lon coordinates from occurrence data
#' @param variables Character vector of climate variables to download
#' @param years Numeric vector of years
#' @param out_dir Character. Output directory
#' @return Data frame with environmental data
get_environmental_data <- function(coordinates,
                                  variables = c("tmean", "ppt", "tmax", "tmin"),
                                  years = 2010:2023,
                                  out_dir = "./data/raw/climate") {
  
  if (!requireNamespace("prism", quietly = TRUE)) {
    stop("Package 'prism' is required for climate data retrieval.")
  }
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Set PRISM download directory
  prism::prism_set_dl_dir(out_dir)
  
  tryCatch({
    env_data_list <- list()
    
    for (var in variables) {
      for (year in years) {
        # Download annual data for each variable
        prism::get_prism_annual(type = var, year = year, keepZip = FALSE)
        
        # Get the raster file
        prism_files <- prism::prism_archive_ls()
        current_file <- prism_files[grep(paste0(var, "_", year), prism_files)]
        
        if (length(current_file) > 0) {
          # Extract values for coordinates
          raster_path <- file.path(prism::prism_get_dl_dir(), current_file[1], 
                                  paste0(current_file[1], ".bil"))
          
          if (file.exists(raster_path)) {
            climate_raster <- terra::rast(raster_path)
            
            # Extract climate data for occurrence coordinates
            if (!is.null(coordinates) && nrow(coordinates) > 0) {
              coords_matrix <- cbind(coordinates$decimalLongitude, 
                                   coordinates$decimalLatitude)
              climate_values <- terra::extract(climate_raster, coords_matrix)
              
              env_data_list[[paste0(var, "_", year)]] <- climate_values
            }
          }
        }
      }
    }
    
    # Combine environmental data
    if (length(env_data_list) > 0) {
      env_df <- do.call(cbind, env_data_list)
      env_df <- cbind(coordinates[, c("decimalLongitude", "decimalLatitude")], env_df)
      
      # Save environmental data
      filename <- file.path(out_dir, "environmental_data.csv")
      utils::write.csv(env_df, filename, row.names = FALSE)
      
      return(env_df)
    }
    
    return(NULL)
    
  }, error = function(e) {
    warning("Environmental data retrieval failed: ", e$message)
    return(NULL)
  })
}

#' Comprehensive marine species data collection
#' @param species Character. Scientific name of the marine species
#' @param include_sequences Logical. Whether to download genomic data
#' @param include_climate Logical. Whether to download environmental data
#' @param out_dir Character. Base output directory
#' @return List with all collected data
collect_marine_species_data <- function(species,
                                       include_sequences = TRUE,
                                       include_climate = TRUE,
                                       out_dir = "./data/raw") {
  
  message("Collecting comprehensive data for ", species, "...")
  
  collected_data <- list()
  
  # 1. Get OBIS occurrence data
  message("Fetching OBIS occurrence data...")
  obis_data <- get_obis_data(species, out_dir = out_dir)
  collected_data$obis <- obis_data
  
  # 2. Get GBIF occurrence data
  message("Fetching GBIF occurrence data...")
  gbif_data <- get_gbif_marine_data(species, out_dir = out_dir)
  collected_data$gbif <- gbif_data
  
  # 3. Get taxonomic information
  message("Fetching taxonomic information...")
  taxonomy <- get_taxonomic_data(species, out_dir = file.path(out_dir, "taxonomy"))
  collected_data$taxonomy <- taxonomy
  
  # 4. Get genomic data (optional)
  if (include_sequences) {
    message("Fetching genomic data...")
    genome_data <- get_sequence_data(species, 
                                   out_dir = file.path(out_dir, "sequences"))
    collected_data$sequences <- genome_data
  }
  
  # 5. Get environmental data (optional)
  if (include_climate && !is.null(obis_data) && nrow(obis_data) > 0) {
    message("Fetching environmental data...")
    coords <- obis_data[!is.na(obis_data$decimalLongitude) & 
                       !is.na(obis_data$decimalLatitude), ]
    
    if (nrow(coords) > 0) {
      climate_data <- get_environmental_data(coords, 
                                           out_dir = file.path(out_dir, "climate"))
      collected_data$climate <- climate_data
    }
  }
  
  message("Data collection completed for ", species)
  return(collected_data)
}

# example usage
# get_obis_data("Acropora")    # writes ./data/raw/acropora_obis_data.csv
# get_gbif_marine_data("Acropora tenuis") # writes ./data/raw/acropora_tenuis_gbif_data.csv
# collect_marine_species_data("Acropora cervicornis") # comprehensive data collection
