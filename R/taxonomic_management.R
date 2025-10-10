# Taxonomic Data Management Functions using taxa package
# For standardized handling of marine biodiversity taxonomic information

#' Parse taxonomic data from multiple sources into taxa objects
#' @param occurrence_data Data frame with occurrence records from OBIS/GBIF
#' @param taxonomy_data List with taxonomic information from wikitaxa
#' @return taxa::taxmap object with parsed taxonomic hierarchy
parse_marine_taxonomy <- function(occurrence_data, taxonomy_data = NULL) {
  
  if (!requireNamespace("taxa", quietly = TRUE)) {
    stop("Package 'taxa' is required for taxonomic data management.")
  }
  
  # Extract taxonomic columns from occurrence data
  tax_cols <- c("kingdom", "phylum", "class", "order", "family", "genus", "species", 
                "scientificName", "taxonomicStatus", "taxonRank")
  
  available_cols <- intersect(tax_cols, names(occurrence_data))
  
  if (length(available_cols) == 0) {
    stop("No taxonomic columns found in occurrence data")
  }
  
  # Create taxmap object from occurrence data
  taxmap_obj <- tryCatch({
    taxa::parse_tax_data(
      tax_data = occurrence_data,
      class_cols = available_cols,
      named_by_rank = TRUE,
      class_sep = ";",
      class_key = c("kingdom" = "kingdom", "phylum" = "phylum", "class" = "class",
                    "order" = "order", "family" = "family", "genus" = "genus", 
                    "species" = "species")
    )
  }, error = function(e) {
    warning("Failed to parse taxonomic data with taxa: ", e$message)
    return(NULL)
  })
  
  return(taxmap_obj)
}

#' Create marine-specific taxonomic hierarchy
#' @param species_list Character vector of marine species names
#' @return taxa::taxonomy object with marine phylogenetic structure
create_marine_taxonomy <- function(species_list) {
  
  if (!requireNamespace("taxa", quietly = TRUE)) {
    stop("Package 'taxa' is required for taxonomic management.")
  }
  
  # Create marine taxonomic hierarchy focusing on common marine groups
  marine_hierarchy <- list()
  
  for (species in species_list) {
    # Basic taxonomic parsing for marine organisms
    if (grepl("Acropora|Porites|Montipora", species, ignore.case = TRUE)) {
      hierarchy <- c("Animalia", "Cnidaria", "Anthozoa", "Scleractinia", 
                    "Acroporidae", species)
      names(hierarchy) <- c("kingdom", "phylum", "class", "order", "family", "species")
    } else if (grepl("Thalassiosira|Skeletonema|Chaetoceros", species, ignore.case = TRUE)) {
      hierarchy <- c("Chromista", "Bacillariophyta", "Bacillariophyceae", 
                    "Centrales", "Thalassiosiraceae", species)
      names(hierarchy) <- c("kingdom", "phylum", "class", "order", "family", "species")
    } else if (grepl("Zostera|Posidonia|Halophila", species, ignore.case = TRUE)) {
      hierarchy <- c("Plantae", "Tracheophyta", "Liliopsida", "Alismatales", 
                    "Hydrocharitaceae", species)
      names(hierarchy) <- c("kingdom", "phylum", "class", "order", "family", "species")
    } else {
      # Generic marine classification
      hierarchy <- c("Marine_Organism", "Unknown_Phylum", "Unknown_Class", 
                    "Unknown_Order", "Unknown_Family", species)
      names(hierarchy) <- c("kingdom", "phylum", "class", "order", "family", "species")
    }
    
    marine_hierarchy[[species]] <- hierarchy
  }
  
  # Create taxon objects
  taxon_list <- list()
  for (i in seq_along(marine_hierarchy)) {
    species_name <- names(marine_hierarchy)[i]
    hierarchy <- marine_hierarchy[[i]]
    
    # Create taxon for each rank
    for (j in seq_along(hierarchy)) {
      rank_name <- names(hierarchy)[j]
      taxon_name <- hierarchy[j]
      
      taxon_list[[paste0(species_name, "_", rank_name)]] <- taxa::taxon(
        name = taxa::taxon_name(taxon_name),
        rank = taxa::taxon_rank(rank_name)
      )
    }
  }
  
  # Create taxonomy object
  marine_taxonomy <- taxa::taxonomy(taxon_list)
  
  return(marine_taxonomy)
}

#' Validate and standardize marine species names
#' @param species_names Character vector of species names to validate
#' @param use_gbif Logical. Whether to use GBIF for name validation
#' @return Data frame with original names, standardized names, and taxonomic status
validate_marine_names <- function(species_names, use_gbif = TRUE) {
  
  if (use_gbif && !requireNamespace("rgbif", quietly = TRUE)) {
    warning("Package 'rgbif' not available for name validation")
    use_gbif <- FALSE
  }
  
  validation_results <- data.frame(
    original_name = species_names,
    standardized_name = NA,
    match_type = NA,
    taxonomic_status = NA,
    accepted_name = NA,
    confidence = NA,
    stringsAsFactors = FALSE
  )
  
  if (use_gbif) {
    for (i in seq_along(species_names)) {
      tryCatch({
        # Use GBIF name backbone for validation
        name_lookup <- rgbif::name_backbone(name = species_names[i])
        
        if (!is.null(name_lookup) && nrow(name_lookup) > 0) {
          validation_results$standardized_name[i] <- name_lookup$canonicalName[1]
          validation_results$match_type[i] <- name_lookup$matchType[1]
          validation_results$taxonomic_status[i] <- name_lookup$status[1]
          validation_results$accepted_name[i] <- name_lookup$species[1]
          validation_results$confidence[i] <- name_lookup$confidence[1]
        }
      }, error = function(e) {
        warning("Name validation failed for ", species_names[i], ": ", e$message)
      })
    }
  } else {
    # Basic name standardization without external validation
    validation_results$standardized_name <- trimws(species_names)
    validation_results$match_type <- "LOCAL"
  }
  
  return(validation_results)
}

#' Filter taxonomic data for marine environments
#' @param taxmap_obj taxa::taxmap object
#' @param marine_keywords Character vector of keywords indicating marine environment
#' @return Filtered taxa::taxmap object with marine taxa only
filter_marine_taxa <- function(taxmap_obj, 
                              marine_keywords = c("marine", "ocean", "sea", "reef", 
                                                 "pelagic", "benthic", "coral")) {
  
  if (!inherits(taxmap_obj, "Taxmap")) {
    stop("Input must be a taxa::taxmap object")
  }
  
  # Filter taxa based on marine keywords in habitat or locality information
  if ("tax_data" %in% names(taxmap_obj$data)) {
    # Create marine indicator
    marine_indicator <- rep(FALSE, nrow(taxmap_obj$data$tax_data))
    
    # Check various columns for marine indicators
    check_cols <- c("habitat", "locality", "waterBody", "marine", "depth", 
                   "establishmentMeans", "occurrenceRemarks")
    
    for (col in check_cols) {
      if (col %in% names(taxmap_obj$data$tax_data)) {
        col_data <- taxmap_obj$data$tax_data[[col]]
        if (!is.null(col_data)) {
          marine_matches <- grepl(paste(marine_keywords, collapse = "|"), 
                                col_data, ignore.case = TRUE)
          marine_matches[is.na(marine_matches)] <- FALSE
          marine_indicator <- marine_indicator | marine_matches
        }
      }
    }
    
    # Also include records with depth > 0 (indicating marine environment)
    if ("depth" %in% names(taxmap_obj$data$tax_data)) {
      depth_data <- as.numeric(taxmap_obj$data$tax_data$depth)
      marine_depth <- !is.na(depth_data) & depth_data > 0
      marine_indicator <- marine_indicator | marine_depth
    }
    
    # Filter the taxmap object
    if (sum(marine_indicator) > 0) {
      filtered_taxmap <- taxa::filter_obs(taxmap_obj, "tax_data", marine_indicator)
      return(filtered_taxmap)
    } else {
      warning("No marine taxa found with the specified criteria")
      return(taxmap_obj)
    }
  }
  
  return(taxmap_obj)
}

#' Create taxonomic summary for marine biodiversity analysis
#' @param taxmap_obj taxa::taxmap object with marine taxonomic data
#' @return List with taxonomic summaries and diversity metrics
summarize_marine_taxonomy <- function(taxmap_obj) {
  
  if (!inherits(taxmap_obj, "Taxmap")) {
    stop("Input must be a taxa::taxmap object")
  }
  
  summary_list <- list()
  
  # Basic taxonomic counts
  summary_list$n_taxa <- taxa::n_taxa(taxmap_obj)
  summary_list$n_edges <- taxa::n_supertaxa(taxmap_obj)
  
  # Taxa by rank
  if (length(taxmap_obj$taxon_ranks()) > 0) {
    rank_counts <- table(taxmap_obj$taxon_ranks())
    summary_list$taxa_by_rank <- rank_counts
  }
  
  # Unique taxa names by rank
  if ("tax_data" %in% names(taxmap_obj$data)) {
    tax_data <- taxmap_obj$data$tax_data
    
    for (rank in c("kingdom", "phylum", "class", "order", "family", "genus", "species")) {
      if (rank %in% names(tax_data)) {
        unique_taxa <- length(unique(tax_data[[rank]][!is.na(tax_data[[rank]])]))
        summary_list[[paste0("unique_", rank)]] <- unique_taxa
      }
    }
  }
  
  # Taxonomic diversity metrics
  if ("tax_data" %in% names(taxmap_obj$data) && "species" %in% names(taxmap_obj$data$tax_data)) {
    species_data <- taxmap_obj$data$tax_data$species
    species_counts <- table(species_data[!is.na(species_data)])
    
    # Shannon diversity
    if (length(species_counts) > 0) {
      p <- species_counts / sum(species_counts)
      shannon_div <- -sum(p * log(p))
      summary_list$shannon_diversity <- shannon_div
      
      # Simpson diversity
      simpson_div <- 1 - sum(p^2)
      summary_list$simpson_diversity <- simpson_div
      
      # Species richness
      summary_list$species_richness <- length(species_counts)
    }
  }
  
  return(summary_list)
}

#' Export taxonomic data in various formats
#' @param taxmap_obj taxa::taxmap object
#' @param output_dir Character. Output directory path
#' @param formats Character vector. Output formats ("csv", "json", "newick")
#' @return Character vector of created file paths
export_taxonomic_data <- function(taxmap_obj, 
                                 output_dir = "./data/processed/taxonomy",
                                 formats = c("csv", "json")) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  created_files <- character()
  
  # Export as CSV
  if ("csv" %in% formats && "tax_data" %in% names(taxmap_obj$data)) {
    csv_file <- file.path(output_dir, "taxonomic_data.csv")
    utils::write.csv(taxmap_obj$data$tax_data, csv_file, row.names = FALSE)
    created_files <- c(created_files, csv_file)
  }
  
  # Export as JSON
  if ("json" %in% formats) {
    json_file <- file.path(output_dir, "taxonomic_hierarchy.json")
    
    # Create hierarchical structure for JSON export
    hierarchy_list <- list(
      taxa = taxmap_obj$taxon_names(),
      ranks = taxmap_obj$taxon_ranks(),
      data = taxmap_obj$data
    )
    
    jsonlite::write_json(hierarchy_list, json_file, pretty = TRUE)
    created_files <- c(created_files, json_file)
  }
  
  # Export as Newick format (if tree structure available)
  if ("newick" %in% formats) {
    tryCatch({
      newick_file <- file.path(output_dir, "taxonomic_tree.newick")
      
      # Simple tree export (requires ape package)
      if (requireNamespace("ape", quietly = TRUE)) {
        # Create a basic phylogenetic tree structure
        # This is a simplified example - real implementation would need proper tree construction
        taxa_names <- taxmap_obj$taxon_names()
        if (length(taxa_names) > 2) {
          # Create a simple tree structure
          tree_string <- paste0("(", paste(taxa_names, collapse = ","), ");")
          writeLines(tree_string, newick_file)
          created_files <- c(created_files, newick_file)
        }
      }
    }, error = function(e) {
      warning("Newick export failed: ", e$message)
    })
  }
  
  return(created_files)
}