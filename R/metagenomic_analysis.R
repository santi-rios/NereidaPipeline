#' Metagenomic Analysis Functions for Marine Species Conservation
#' 
#' Integrated pipeline for genomic sequence retrieval using biomartr package
#' for organism-centric data access

library(biomartr)
library(dplyr)
library(Biostrings)

#' Retrieve genome for marine species using biomartr
#' 
#' @param species_list Character vector of species names
#' @param db Database to query ("refseq", "genbank", "ensembl")
#' @param out_dir Output directory for downloaded sequences
#' @return List of genome retrieval results
#' @export
retrieve_marine_genomes <- function(species_list, 
                                   db = "refseq",
                                   out_dir = "data/raw/genomic") {
  
  # Create output directory
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Set timeout for large downloads
  options(timeout = 30000)
  
  results <- list()
  
  for (species in species_list) {
    message("Retrieving genome for ", species)
    
    tryCatch({
      # Check if genome is available
      is_available <- biomartr::is.genome.available(
        organism = species,
        db = db,
        details = TRUE
      )
      
      if (nrow(is_available) > 0) {
        message("Found ", nrow(is_available), " genome assembly/assemblies for ", species)
        
        # Get the genome
        genome_file <- biomartr::getGenome(
          db = db,
          organism = species,
          path = out_dir
        )
        
        results[[species]] <- list(
          status = "success",
          file = genome_file,
          database = db,
          assembly_info = is_available,
          download_date = Sys.time()
        )
        
        message("✓ Successfully retrieved genome for ", species)
      } else {
        results[[species]] <- list(
          status = "not_available",
          message = paste(species, "genome not found in", db)
        )
        message("✗ ", species, " genome not available in ", db)
      }
      
    }, error = function(e) {
      results[[species]] <- list(
        status = "error",
        message = e$message
      )
      warning("Failed to retrieve genome for ", species, ": ", e$message)
    })
  }
  
  return(results)
}

#' Retrieve proteome for marine species using biomartr
#' 
#' @param species_list Character vector of species names
#' @param db Database to query ("refseq", "genbank", "ensembl")
#' @param out_dir Output directory for downloaded sequences
#' @return List of proteome retrieval results
#' @export
retrieve_marine_proteomes <- function(species_list, 
                                     db = "refseq",
                                     out_dir = "data/raw/genomic/proteomes") {
  
  # Create output directory
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Set timeout for large downloads
  options(timeout = 30000)
  
  results <- list()
  
  for (species in species_list) {
    message("Retrieving proteome for ", species)
    
    tryCatch({
      # Get the proteome
      proteome_file <- biomartr::getProteome(
        db = db,
        organism = species,
        path = out_dir
      )
      
      results[[species]] <- list(
        status = "success",
        file = proteome_file,
        database = db,
        download_date = Sys.time()
      )
      
      message("✓ Successfully retrieved proteome for ", species)
      
    }, error = function(e) {
      results[[species]] <- list(
        status = "error",
        message = e$message
      )
      warning("Failed to retrieve proteome for ", species, ": ", e$message)
    })
  }
  
  return(results)
}

#' Retrieve CDS (coding sequences) for marine species using biomartr
#' 
#' @param species_list Character vector of species names
#' @param db Database to query ("refseq", "genbank", "ensembl")
#' @param out_dir Output directory for downloaded sequences
#' @return List of CDS retrieval results
#' @export
retrieve_marine_cds <- function(species_list, 
                               db = "refseq",
                               out_dir = "data/raw/genomic/cds") {
  
  # Create output directory
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Set timeout for large downloads
  options(timeout = 30000)
  
  results <- list()
  
  for (species in species_list) {
    message("Retrieving CDS for ", species)
    
    tryCatch({
      # Get the CDS
      cds_file <- biomartr::getCDS(
        db = db,
        organism = species,
        path = out_dir
      )
      
      results[[species]] <- list(
        status = "success",
        file = cds_file,
        database = db,
        download_date = Sys.time()
      )
      
      message("✓ Successfully retrieved CDS for ", species)
      
    }, error = function(e) {
      results[[species]] <- list(
        status = "error",
        message = e$message
      )
      warning("Failed to retrieve CDS for ", species, ": ", e$message)
    })
  }
  
  return(results)
}

#' Retrieve GFF annotation for marine species using biomartr
#' 
#' @param species_list Character vector of species names
#' @param db Database to query ("refseq", "genbank", "ensembl")
#' @param out_dir Output directory for downloaded annotations
#' @return List of GFF retrieval results
#' @export
retrieve_marine_gff <- function(species_list, 
                               db = "refseq",
                               out_dir = "data/raw/genomic/gff") {
  
  # Create output directory
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Set timeout for large downloads
  options(timeout = 30000)
  
  results <- list()
  
  for (species in species_list) {
    message("Retrieving GFF annotation for ", species)
    
    tryCatch({
      # Get the GFF
      gff_file <- biomartr::getGFF(
        db = db,
        organism = species,
        path = out_dir
      )
      
      results[[species]] <- list(
        status = "success",
        file = gff_file,
        database = db,
        download_date = Sys.time()
      )
      
      message("✓ Successfully retrieved GFF for ", species)
      
    }, error = function(e) {
      results[[species]] <- list(
        status = "error",
        message = e$message
      )
      warning("Failed to retrieve GFF for ", species, ": ", e$message)
    })
  }
  
  return(results)
}

#' Retrieve RNA sequences for marine species using biomartr
#' 
#' @param species_list Character vector of species names
#' @param db Database to query ("refseq", "genbank", "ensembl")
#' @param out_dir Output directory for downloaded sequences
#' @return List of RNA retrieval results
#' @export
retrieve_marine_rna <- function(species_list, 
                               db = "refseq",
                               out_dir = "data/raw/genomic/rna") {
  
  # Create output directory
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Set timeout for large downloads
  options(timeout = 30000)
  
  results <- list()
  
  for (species in species_list) {
    message("Retrieving RNA for ", species)
    
    tryCatch({
      # Get the RNA
      rna_file <- biomartr::getRNA(
        db = db,
        organism = species,
        path = out_dir
      )
      
      results[[species]] <- list(
        status = "success",
        file = rna_file,
        database = db,
        download_date = Sys.time()
      )
      
      message("✓ Successfully retrieved RNA for ", species)
      
    }, error = function(e) {
      results[[species]] <- list(
        status = "error",
        message = e$message
      )
      warning("Failed to retrieve RNA for ", species, ": ", e$message)
    })
  }
  
  return(results)
}

#' Get assembly statistics for marine species
#' 
#' @param species_list Character vector of species names
#' @param db Database to query ("refseq", "genbank", "ensembl")
#' @param out_dir Output directory
#' @return Data frame with assembly statistics
#' @export
get_assembly_stats <- function(species_list,
                              db = "refseq",
                              out_dir = "data/raw/genomic/assembly_stats") {
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  all_stats <- list()
  
  for (species in species_list) {
    message("Retrieving assembly stats for ", species)
    
    tryCatch({
      # Get assembly statistics
      assembly_info <- biomartr::getAssemblyStats(
        db = db,
        organism = species,
        path = out_dir
      )
      
      if (!is.null(assembly_info) && nrow(assembly_info) > 0) {
        # Add species name
        assembly_info$species <- species
        assembly_info$retrieval_date <- Sys.time()
        
        all_stats[[species]] <- assembly_info
        message("✓ Retrieved assembly stats for ", species)
      }
      
    }, error = function(e) {
      warning("Failed to get stats for ", species, ": ", e$message)
    })
  }
  
  if (length(all_stats) > 0) {
    combined_stats <- dplyr::bind_rows(all_stats)
    
    # Save statistics
    write.csv(
      combined_stats,
      file.path(out_dir, "assembly_statistics.csv"),
      row.names = FALSE
    )
    
    return(combined_stats)
  } else {
    return(NULL)
  }
}

#' Check genome availability for marine species
#' 
#' @param species_list Character vector of species names
#' @param db Database to query ("refseq", "genbank", "ensembl")
#' @return Data frame with availability information
#' @export
check_genome_availability <- function(species_list, db = "refseq") {
  
  availability_list <- list()
  
  for (species in species_list) {
    message("Checking genome availability for ", species, " in ", db)
    
    tryCatch({
      is_available <- biomartr::is.genome.available(
        organism = species,
        db = db,
        details = TRUE
      )
      
      if (nrow(is_available) > 0) {
        is_available$species <- species
        is_available$database <- db
        availability_list[[species]] <- is_available
        message("✓ Found ", nrow(is_available), " assembly/assemblies for ", species)
      } else {
        message("✗ No genome found for ", species, " in ", db)
      }
      
    }, error = function(e) {
      warning("Error checking availability for ", species, ": ", e$message)
    })
  }
  
  if (length(availability_list) > 0) {
    return(dplyr::bind_rows(availability_list))
  } else {
    return(NULL)
  }
}

#' Check genome availability across multiple databases
#' 
#' @param species_list Character vector of species names
#' @return Data frame with availability information across databases
#' @export
check_multi_database_availability <- function(species_list) {
  
  results <- list()
  databases <- c("refseq", "genbank", "ensembl")
  
  for (db in databases) {
    message("Checking availability in ", db, "...")
    
    for (species in species_list) {
      tryCatch({
        is_available <- biomartr::is.genome.available(
          organism = species,
          db = db,
          details = TRUE
        )
        
        if (!is.null(is_available) && nrow(is_available) > 0) {
          # Convert date columns to character to avoid type conflicts
          if ("seq_rel_date" %in% names(is_available)) {
            is_available$seq_rel_date <- as.character(is_available$seq_rel_date)
          }
          if ("submission_date" %in% names(is_available)) {
            is_available$submission_date <- as.character(is_available$submission_date)
          }
          if ("release_date" %in% names(is_available)) {
            is_available$release_date <- as.character(is_available$release_date)
          }
          
          is_available$species <- species
          is_available$database <- db
          results[[paste(species, db, sep = "_")]] <- is_available
        }
        
      }, error = function(e) {
        message("  - ", species, " not found in ", db)
      })
    }
  }
  
  if (length(results) > 0) {
    combined <- dplyr::bind_rows(results)
    return(combined)
  } else {
    return(NULL)
  }
}

#' Retrieve complete metagenomic dataset with fallback databases
#' 
#' @param species_list Character vector of species names
#' @param data_types Vector of data types to retrieve
#' @param out_dir Base output directory
#' @param try_databases Vector of databases to try in order
#' @return List with all retrieval results
#' @export
retrieve_complete_metagenomic_data_smart <- function(species_list,
                                                      data_types = c("proteome", "cds", "gff"),
                                                      out_dir = "data/raw/genomic",
                                                      try_databases = c("refseq", "genbank")) {
  
  message("Starting smart metagenomic data retrieval for ", length(species_list), " species")
  message("Will try databases in order: ", paste(try_databases, collapse = ", "))
  
  # First, check availability across all databases
  availability <- check_multi_database_availability(species_list)
  
  if (!is.null(availability)) {
    write.csv(
      availability,
      file.path(out_dir, "species_availability.csv"),
      row.names = FALSE
    )
    message("Saved availability report to: ", file.path(out_dir, "species_availability.csv"))
  }
  
  results <- list()
  results$availability <- availability
  
  # For each species, try databases in order
  for (species in species_list) {
    species_results <- list()
    
    # Find which databases have this species
    if (!is.null(availability)) {
      available_dbs <- availability %>%
        filter(species == !!species) %>%
        pull(database)
    } else {
      available_dbs <- character(0)
    }
    
    # Select database (prefer order in try_databases)
    selected_db <- NULL
    for (db in try_databases) {
      if (db %in% available_dbs) {
        selected_db <- db
        break
      }
    }
    
    if (is.null(selected_db)) {
      message("✗ ", species, " not available in any database")
      species_results$status <- "not_available"
      results[[species]] <- species_results
      next
    }
    
    message("\nRetrieving ", species, " from ", selected_db)
    species_results$database_used <- selected_db
    
    # Retrieve each data type
    for (dtype in data_types) {
      tryCatch({
        result <- switch(dtype,
          "genome" = biomartr::getGenome(
            db = selected_db,
            organism = species,
            path = file.path(out_dir, "genomes"),
            reference = FALSE
          ),
          "proteome" = biomartr::getProteome(
            db = selected_db,
            organism = species,
            path = file.path(out_dir, "proteomes"),
            reference = FALSE
          ),
          "cds" = biomartr::getCDS(
            db = selected_db,
            organism = species,
            path = file.path(out_dir, "cds"),
            reference = FALSE
          ),
          "gff" = biomartr::getGFF(
            db = selected_db,
            organism = species,
            path = file.path(out_dir, "gff"),
            reference = FALSE
          ),
          "rna" = biomartr::getRNA(
            db = selected_db,
            organism = species,
            path = file.path(out_dir, "rna"),
            reference = FALSE
          )
        )
        
        species_results[[dtype]] <- list(
          status = "success",
          file = result,
          download_date = Sys.time()
        )
        
      }, error = function(e) {
        species_results[[dtype]] <- list(
          status = "error",
          message = e$message
        )
        message("  ✗ Failed to retrieve ", dtype, ": ", e$message)
      })
    }
    
    results[[species]] <- species_results
  }
  
  # Save summary report
  saveRDS(results, file.path(out_dir, "retrieval_summary.rds"))
  message("\n✓ Complete metagenomic data retrieval finished")
  message("Summary saved to: ", file.path(out_dir, "retrieval_summary.rds"))
  
  return(results)
}

#' DNA Barcoding Analysis for Species Identification
#' 
#' @param reference_sequences DNAbin object with reference sequences
#' @param query_sequences DNAbin object with query sequences
#' @param method Identification method ("fuzzyId", "Bayesian", "bpNewTraining")
#' @param out_dir Output directory for results
#' @return List with identification results
#' @export
# perform_dna_barcoding <- function(reference_sequences,
#                                  query_sequences,
#                                  method = "fuzzyId",
#                                  out_dir = "data/processed/barcoding") {
  
#   if (!dir.exists(out_dir)) {
#     dir.create(out_dir, recursive = TRUE)
#   }
  
#   if (!requireNamespace("BarcodingR", quietly = TRUE)) {
#     stop("BarcodingR package is required but not installed")
#   }
  
#   message("Performing DNA barcoding with method: ", method)
  
#   tryCatch({
#     # Perform species identification
#     barcoding_result <- BarcodingR::barcoding.spe.identify(
#       ref = reference_sequences,
#       que = query_sequences,
#       method = method
#     )
    
#     # Calculate barcoding gap
#     barcode_gap <- BarcodingR::barcoding.gap(
#       ref = reference_sequences,
#       dist = "K80"
#     )
    
#     # Summarize reference dataset
#     ref_summary <- BarcodingR::summarize.ref(
#       ref = reference_sequences,
#       taxonStat = TRUE,
#       seqStat = TRUE,
#       barcodeStat = TRUE
#     )
    
#     results <- list(
#       identification = barcoding_result,
#       barcoding_gap = barcode_gap,
#       reference_summary = ref_summary,
#       method = method,
#       analysis_date = Sys.time()
#     )
    
#     # Save results
#     saveRDS(results, file.path(out_dir, "barcoding_results.rds"))
    
#     message("✓ DNA barcoding analysis completed")
#     return(results)
    
#   }, error = function(e) {
#     warning("DNA barcoding analysis failed: ", e$message)
#     return(NULL)
#   })
# }

#' Evolutionary Transcriptomics Analysis
#' 
#' @param expression_data Phylostratigraphic expression data
#' @param analysis_type Type of analysis ("TAI", "TDI", "both")
#' @param out_dir Output directory
#' @return List with evolutionary analysis results
#' @export
# analyze_evolutionary_transcriptomics <- function(expression_data,
#                                                 analysis_type = "both",
#                                                 out_dir = "data/processed/evolutionary") {
  
#   if (!dir.exists(out_dir)) {
#     dir.create(out_dir, recursive = TRUE)
#   }
  
#   if (!requireNamespace("myTAI", quietly = TRUE)) {
#     stop("myTAI package is required but not installed")
#   }
  
#   results <- list()
  
#   tryCatch({
#     # Transcriptome Age Index (TAI)
#     if (analysis_type %in% c("TAI", "both")) {
#       message("Calculating Transcriptome Age Index (TAI)...")
      
#       # Calculate TAI
#       tai_results <- myTAI::TAI(expression_data)
      
#       # Plot TAI signature
#       tai_plot <- myTAI::PlotSignature(
#         ExpressionSet = expression_data,
#         measure = "TAI",
#         TestStatistic = "FlatLineTest"
#       )
      
#       results$TAI <- list(
#         values = tai_results,
#         plot = tai_plot
#       )
#     }
    
#     # Transcriptome Divergence Index (TDI)
#     if (analysis_type %in% c("TDI", "both")) {
#       message("Calculating Transcriptome Divergence Index (TDI)...")
      
#       # Calculate TDI
#       tdi_results <- myTAI::TDI(expression_data)
      
#       results$TDI <- list(
#         values = tdi_results
#       )
#     }
    
#     # Gene age category expression
#     message("Analyzing gene age category expression...")
    
#     category_expr <- myTAI::PlotCategoryExpr(
#       ExpressionSet = expression_data,
#       legendName = "Phylostratum",
#       log.expr = TRUE
#     )
    
#     results$category_expression <- category_expr
    
#     # Relative expression analysis
#     message("Performing relative expression analysis...")
    
#     re_analysis <- myTAI::PlotRE(
#       ExpressionSet = expression_data,
#       Groups = list(old = 1:3, young = 4:12),
#       legendName = "PS"
#     )
    
#     results$relative_expression <- re_analysis
    
#     # Save results
#     saveRDS(results, file.path(out_dir, "evolutionary_transcriptomics.rds"))
    
#     message("✓ Evolutionary transcriptomics analysis completed")
#     return(results)
    
#   }, error = function(e) {
#     warning("Evolutionary transcriptomics analysis failed: ", e$message)
#     return(NULL)
#   })
# }

#' Create phylostratigraphic map for gene age analysis using Ensembl
#' 
#' @param species Scientific name of species
#' @param out_dir Output directory
#' @return Phylostratigraphic map data
#' @export
# create_phylostratigraphic_map <- function(species,
#                                          out_dir = "data/processed/phylostrat") {
  
#   if (!dir.exists(out_dir)) {
#     dir.create(out_dir, recursive = TRUE)
#   }
  
#   message("Creating phylostratigraphic map for ", species)
  
#   tryCatch({
#     ensembl <- useEnsembl(biomart = "genes")
#     organism_name <- gsub(" ", "_", tolower(species))
    
#     # Search for matching dataset
#     matching_datasets <- searchDatasets(ensembl, pattern = organism_name)
    
#     if (nrow(matching_datasets) > 0) {
#       dataset_name <- matching_datasets$dataset[1]
#       species_mart <- useDataset(dataset = dataset_name, mart = ensembl)
      
#       # Get gene information with homology data
#       genes <- getBM(
#         attributes = c(
#           "ensembl_gene_id",
#           "external_gene_name",
#           "gene_biotype",
#           "description"
#         ),
#         mart = species_mart
#       )
      
#       # Create a basic phylostratigraphic map
#       # This is simplified - real implementation would need orthology data
#       phylostrat_map <- genes %>%
#         mutate(
#           GeneID = ensembl_gene_id,
#           Phylostratum = NA_integer_,  # Would need orthology analysis
#           Gene_Age = "Unknown",
#           Gene_Name = external_gene_name,
#           Description = description
#         ) %>%
#         select(GeneID, Gene_Name, Phylostratum, Gene_Age, Description)
      
#       # Save map
#       write.csv(
#         phylostrat_map,
#         file.path(out_dir, paste0(gsub(" ", "_", species), "_phylostrat_map.csv")),
#         row.names = FALSE
#       )
      
#       message("✓ Phylostratigraphic map created with ", nrow(phylostrat_map), " genes")
#       return(phylostrat_map)
#     } else {
#       warning("Species not found in Ensembl")
#       return(NULL)
#     }
    
#   }, error = function(e) {
#     warning("Failed to create phylostratigraphic map: ", e$message)
#     return(NULL)
#   })
# }

#' Assess genome assembly quality based on gene statistics
#' 
#' @param assembly_stats Data frame with gene statistics
#' @param quality_threshold Threshold for quality assessment
#' @return Data frame with quality assessments
#' @export
# assess_assembly_quality <- function(assembly_stats,
#                                    quality_threshold = list(
#                                      min_genes = 10000,
#                                      min_protein_coding = 5000,
#                                      min_chromosomes = 1
#                                    )) {
  
#   if (is.null(assembly_stats) || nrow(assembly_stats) == 0) {
#     warning("No assembly statistics provided")
#     return(NULL)
#   }
  
#   message("Assessing genome assembly quality...")
  
#   # Calculate quality scores
#   quality_assessment <- assembly_stats %>%
#     mutate(
#       gene_count_quality = ifelse(
#         total_genes >= quality_threshold$min_genes,
#         "Pass", "Fail"
#       ),
#       protein_coding_quality = ifelse(
#         protein_coding >= quality_threshold$min_protein_coding,
#         "Pass", "Fail"
#       ),
#       chromosome_quality = ifelse(
#         chromosomes >= quality_threshold$min_chromosomes,
#         "Pass", "Fail"
#       ),
#       overall_quality = case_when(
#         gene_count_quality == "Pass" & 
#           protein_coding_quality == "Pass" & 
#           chromosome_quality == "Pass" ~ "High",
#         gene_count_quality == "Pass" | 
#           protein_coding_quality == "Pass" ~ "Medium",
#         TRUE ~ "Low"
#       )
#     )
  
#   message("✓ Quality assessment completed for ", nrow(quality_assessment), " species")
  
#   return(quality_assessment)
# }

#' Integrate metagenomic data with occurrence data
#' 
#' @param occurrence_data Data frame with occurrence records
#' @param genomic_data List with genomic analysis results
#' @param taxonomic_data Taxonomic hierarchy data
#' @return Integrated dataset
#' @export
# integrate_metagenomic_occurrences <- function(occurrence_data,
#                                              genomic_data,
#                                              taxonomic_data) {
  
#   message("Integrating metagenomic data with occurrence records...")
  
#   if (is.null(occurrence_data) || nrow(occurrence_data) == 0) {
#     warning("No occurrence data provided")
#     return(NULL)
#   }
  
#   # Add genomic data availability flags
#   integrated_data <- occurrence_data %>%
#     dplyr::mutate(
#       has_genome = scientificName %in% names(genomic_data),
#       genome_quality = NA_character_,
#       phylostrat_available = FALSE
#     )
  
#   # Add quality information where available
#   for (species in names(genomic_data)) {
#     if (!is.null(genomic_data[[species]]$quality)) {
#       integrated_data <- integrated_data %>%
#         dplyr::mutate(
#           genome_quality = ifelse(
#             scientificName == species,
#             genomic_data[[species]]$quality,
#             genome_quality
#           )
#         )
#     }
#   }
  
#   message("✓ Integration completed for ", nrow(integrated_data), " records")
  
#   return(integrated_data)
# }

#' Generate metagenomic conservation priorities
#' 
#' @param integrated_data Integrated occurrence and genomic data
#' @param criteria Priority criteria
#' @return Data frame with conservation priorities
#' @export
# generate_metagenomic_priorities <- function(integrated_data,
#                                            criteria = list(
#                                              genome_available = 10,
#                                              high_quality = 5,
#                                              rare_species = 15,
#                                              evolutionary_unique = 20
#                                            )) {
  
#   message("Generating conservation priorities based on metagenomic data...")
  
#   if (is.null(integrated_data) || nrow(integrated_data) == 0) {
#     warning("No integrated data provided")
#     return(NULL)
#   }
  
#   # Calculate priority scores
#   priorities <- integrated_data %>%
#     dplyr::group_by(scientificName) %>%
#     dplyr::summarise(
#       n_occurrences = n(),
#       has_genome = any(has_genome),
#       high_quality_genome = any(genome_quality == "High", na.rm = TRUE),
#       geographic_range = sqrt(
#         (max(decimalLatitude, na.rm = TRUE) - min(decimalLatitude, na.rm = TRUE))^2 +
#         (max(decimalLongitude, na.rm = TRUE) - min(decimalLongitude, na.rm = TRUE))^2
#       ),
#       .groups = "drop"
#     ) %>%
#     dplyr::mutate(
#       rarity_score = case_when(
#         n_occurrences < 50 ~ criteria$rare_species,
#         n_occurrences < 100 ~ criteria$rare_species * 0.5,
#         TRUE ~ 0
#       ),
#       genome_score = ifelse(has_genome, criteria$genome_available, 0),
#       quality_score = ifelse(high_quality_genome, criteria$high_quality, 0),
#       total_priority_score = rarity_score + genome_score + quality_score,
#       priority_category = case_when(
#         total_priority_score >= 25 ~ "CRITICAL",
#         total_priority_score >= 15 ~ "HIGH",
#         total_priority_score >= 10 ~ "MEDIUM",
#         TRUE ~ "LOW"
#       )
#     ) %>%
#     dplyr::arrange(desc(total_priority_score))
  
#   message("✓ Conservation priorities generated for ", nrow(priorities), " species")
  
#   return(priorities)
# }
#' Assess quality of retrieved genomic data
#' 
#' @param genomic_sequences Results from retrieve_complete_metagenomic_data_smart
#' @return Data frame with quality metrics
#' @export
assess_genomic_data_quality <- function(genomic_sequences) {
  
  quality_results <- list()
  
  for (species in names(genomic_sequences)) {
    if (species == "availability") next
    
    sp_data <- genomic_sequences[[species]]
    
    quality_metrics <- list(
      species = species,
      database = sp_data$database_used %||% NA,
      data_completeness = 0,
      proteome_available = FALSE,
      cds_available = FALSE,
      gff_available = FALSE,
      proteome_size = NA,
      cds_count = NA
    )
    
    # Check proteome
    if (!is.null(sp_data$proteome) && sp_data$proteome$status == "success") {
      quality_metrics$proteome_available <- TRUE
      quality_metrics$data_completeness <- quality_metrics$data_completeness + 1
      
      # Get file size
      if (file.exists(sp_data$proteome$file)) {
        quality_metrics$proteome_size <- file.size(sp_data$proteome$file) / (1024^2) # MB
      }
    }
    
    # Check CDS
    if (!is.null(sp_data$cds) && sp_data$cds$status == "success") {
      quality_metrics$cds_available <- TRUE
      quality_metrics$data_completeness <- quality_metrics$data_completeness + 1
      
      # Count sequences if file exists
      if (file.exists(sp_data$cds$file)) {
        tryCatch({
          seqs <- Biostrings::readDNAStringSet(sp_data$cds$file)
          quality_metrics$cds_count <- length(seqs)
        }, error = function(e) {
          message("Could not read CDS file for ", species)
        })
      }
    }
    
    # Check GFF
    if (!is.null(sp_data$gff) && sp_data$gff$status == "success") {
      quality_metrics$gff_available <- TRUE
      quality_metrics$data_completeness <- quality_metrics$data_completeness + 1
    }
    
    quality_metrics$data_completeness <- quality_metrics$data_completeness / 3 * 100
    
    quality_results[[species]] <- quality_metrics
  }
  
  dplyr::bind_rows(quality_results)
}

#' Generate visualization of data retrieval success
#' 
#' @param quality_data Output from assess_genomic_data_quality
#' @param out_dir Output directory for plot
#' @export
plot_retrieval_summary <- function(quality_data, out_dir = "data/processed/genomic") {
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  library(ggplot2)
  
  # Reshape data for plotting
  plot_data <- quality_data %>%
    tidyr::pivot_longer(
      cols = c(proteome_available, cds_available, gff_available),
      names_to = "data_type",
      values_to = "available"
    ) %>%
    mutate(
      data_type = gsub("_available", "", data_type),
      data_type = tools::toTitleCase(data_type)
    )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = species, y = data_type, fill = available)) +
    geom_tile(color = "white", size = 1) +
    scale_fill_manual(
      values = c("TRUE" = "#2ecc71", "FALSE" = "#e74c3c"),
      labels = c("TRUE" = "Available", "FALSE" = "Not Available")
    ) +
    labs(
      title = "Genomic Data Retrieval Summary",
      subtitle = paste("Retrieved from:", unique(quality_data$database)),
      x = "Species",
      y = "Data Type",
      fill = "Status"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )
  
  ggsave(
    file.path(out_dir, "retrieval_summary.png"),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  message("Plot saved to: ", file.path(out_dir, "retrieval_summary.png"))
  
  return(p)
}

#' @title Create Phyloseq Object
#' @description Loads microbiome data from flat files and creates a phyloseq object.
#' @param otu_file Path to the OTU table.
#' @param tax_file Path to the taxonomy table.
#' @param meta_file Path to the sample metadata.
#' @param tree_file Path to the phylogenetic tree.
#' @return A phyloseq object.
#' @export
create_phyloseq_object <- function(otu_file, tax_file, meta_file, tree_file) {
  # Import the data
  otu_table <- phyloseq::otu_table(read.delim(otu_file, row.names = 1), taxa_are_rows = TRUE)
  taxonomy_table <- phyloseq::tax_table(as.matrix(read.delim(tax_file, row.names = 1)))
  metadata_table <- phyloseq::sample_data(read.delim(meta_file, row.names = 1))
  phylo_tree <- phyloseq::read_tree(tree_file)

  # Combine into a phyloseq object
  physeq <- phyloseq::phyloseq(otu_table, taxonomy_table, metadata_table, phylo_tree)

  return(physeq)
}


#' @title Download SRA Data
#' @description Downloads FASTQ data from SRA using the SRA Toolkit.
#' @param sra_ids A vector of SRA Run IDs to download.
#' @param output_dir The directory where the data will be downloaded.
#' @return A vector with the paths to the downloaded files.
#' @export
download_sra_data <- function(sra_ids, output_dir = "data/raw/metagenomic") {
  # This function requires the SRA Toolkit to be installed and in the system's PATH.
  # Check if fastq-dump command is available
  if (Sys.which("fastq-dump") == "") {
    stop(
      "SRA Toolkit not found. Please install it and ensure 'fastq-dump' is in your PATH.\n",
      "Installation (Debian/Ubuntu): sudo apt-get install sra-toolkit"
    )
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  message("Downloading SRA data. This may take a significant amount of time...")

  final_file_paths <- c()

  for (sra_id in sra_ids) {
    # Define expected output file. We'll assume single-end reads for now based on metadata.
    # The tool will create _1.fastq.gz for paired-end, but our metadata points to one file.
    expected_file <- file.path(output_dir, paste0(sra_id, "_1.fastq.gz"))
    
    # Also check for the single-end naming convention
    expected_file_single <- file.path(output_dir, paste0(sra_id, ".fastq.gz"))

    # If the file already exists, skip downloading
    if (file.exists(expected_file) || file.exists(expected_file_single)) {
      message("✓ File for ", sra_id, " already exists. Skipping download.")
      if(file.exists(expected_file)) {
        final_file_paths <- c(final_file_paths, expected_file)
      } else {
        final_file_paths <- c(final_file_paths, expected_file_single)
      }
      next
    }

    message("-> Downloading ", sra_id, "...")
    
    # Use system2 to call fastq-dump
    # --split-files: creates _1.fastq, _2.fastq etc. for paired-end reads
    # --gzip: compresses the output
    # -O: specifies the output directory
    tryCatch({
      system2(
        "fastq-dump",
        args = c(
          "--split-files",
          "--gzip",
          "-O", output_dir,
          sra_id
        ),
        stdout = TRUE, # Capture standard output
        stderr = TRUE  # Capture standard error
      )
      
      # Verify download
      if (file.exists(expected_file) || file.exists(expected_file_single)) {
        message("✓ Successfully downloaded ", sra_id)
        if(file.exists(expected_file)) {
          final_file_paths <- c(final_file_paths, expected_file)
        } else {
          final_file_paths <- c(final_file_paths, expected_file_single)
        }
      } else {
         # Sometimes fastq-dump creates files without the _1 suffix for single-end reads
         # Let's rename it to match our convention for consistency
         unsplit_file <- file.path(output_dir, paste0(sra_id, ".fastq.gz"))
         if (file.exists(unsplit_file)) {
            file.rename(unsplit_file, expected_file)
            message("✓ Successfully downloaded and renamed ", sra_id)
            final_file_paths <- c(final_file_paths, expected_file)
         } else {
            warning("✗ Download for ", sra_id, " may have failed. Expected file not found.")
         }
      }

    }, error = function(e) {
      warning("✗ Failed to download ", sra_id, ": ", e$message)
    })
  }
  
  # The file names in the metadata are CORAL_001_R1.fastq.gz, let's rename them
  # to match what the pipeline expects.
  # ERR2043338 -> CORAL_001_R1.fastq.gz
  # ERR2043339 -> CORAL_002_R1.fastq.gz
  
  # This part is tricky because it relies on the metadata.
  # For now, let's assume the downloaded files are what we need.
  # The pipeline might need adjustment if file names are a problem.
  
  return(final_file_paths)
}


#' @title Run FastQC for quality control
#' @description This function runs FastQC on the raw sequencing reads.
#' @param fastq_files A character vector of paths to the FASTQ files.
#' @param qc_output_dir The directory where FastQC reports will be saved.
#' @return A character vector of paths to the generated FastQC zip files.
#' @export
run_fastqc <- function(fastq_files, qc_output_dir = "data/processed/fastqc_reports") {
  if (!requireNamespace("fastqcr", quietly = TRUE)) {
    stop("Package 'fastqcr' is required. Please install it via install.packages('fastqcr').")
  }

  if (length(fastq_files) == 0) {
    stop("No FASTQ files provided to run_fastqc.")
  }

  # The fastqc function works on a directory of fastq files.
  # We'll get the directory from the file paths.
  fastq_dir <- unique(dirname(fastq_files))
  if (length(fastq_dir) > 1) {
    stop("All FASTQ files must be in the same directory to run fastqc.")
  }

  if (!dir.exists(qc_output_dir)) {
    dir.create(qc_output_dir, recursive = TRUE)
  }

  # Run FastQC
  fastqcr::fastqc(
    fq.dir = fastq_dir,
    qc.dir = qc_output_dir,
    threads = 4 # Use 4 threads as an example
  )

  # Return the paths to the generated .zip files for tracking
  list.files(qc_output_dir, pattern = "_fastqc\\.zip$", full.names = TRUE)
}

#' @title Aggregate FastQC reports
#' @description This function aggregates multiple FastQC reports into a single object.
#' @param qc_files A character vector of paths to the FastQC zip files.
#' @param qc_summary_path Path to save the aggregated QC summary report.
#' @return The aggregated QC data frame.
#' @export
aggregate_fastqc_reports <- function(qc_files, qc_summary_path = "data/processed/fastqc_summary.csv") {
    if (!requireNamespace("fastqcr", quietly = TRUE)) {
    stop("Package 'fastqcr' is required. Please install it via install.packages('fastqcr').")
  }
  
  qc_dir <- unique(dirname(qc_files))
   if (length(qc_dir) > 1) {
    stop("All FastQC zip files must be in the same directory.")
  }

  aggregated_qc <- fastqcr::qc_aggregate(qc.dir = qc_dir)
  
  # Save the summary table
  write.csv(aggregated_qc, qc_summary_path, row.names = FALSE)
  
  return(qc_summary_path)
}