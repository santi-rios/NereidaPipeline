# Evolutionary Transcriptomics Analysis using myTAI
# For phylostratigraphic analysis and evolutionary studies in marine organisms

#' Prepare expression data for phylostratigraphic analysis
#' @param expression_data Data frame with gene expression data (genes x samples)
#' @param phylostrata_data Data frame mapping genes to phylostratigraphic ages
#' @param gene_id_col Character. Column name for gene identifiers
#' @return myTAI-compatible PhyloExpressionSet
prepare_marine_expression_set <- function(expression_data,
                                        phylostrata_data,
                                        gene_id_col = "gene_id") {
  
  if (!requireNamespace("myTAI", quietly = TRUE)) {
    stop("Package 'myTAI' is required for evolutionary transcriptomics analysis.")
  }
  
  # Ensure required columns exist
  if (!gene_id_col %in% names(phylostrata_data)) {
    stop("Gene ID column not found in phylostrata data")
  }
  
  if (!"phylostratum" %in% names(phylostrata_data)) {
    stop("'phylostratum' column not found in phylostrata data")
  }
  
  # Match genes between expression and phylostratum data
  common_genes <- intersect(rownames(expression_data), 
                           phylostrata_data[[gene_id_col]])
  
  if (length(common_genes) == 0) {
    stop("No common genes found between expression data and phylostrata")
  }
  
  message("Found ", length(common_genes), " genes with phylostratigraphic information")
  
  # Filter expression data to common genes
  filtered_expression <- expression_data[common_genes, , drop = FALSE]
  
  # Create phylostratum mapping
  phylo_mapping <- phylostrata_data[
    phylostrata_data[[gene_id_col]] %in% common_genes,
    c(gene_id_col, "phylostratum")
  ]
  
  # Ensure proper ordering
  phylo_mapping <- phylo_mapping[match(common_genes, phylo_mapping[[gene_id_col]]), ]
  
  # Create PhyloExpressionSet
  phylo_expression_set <- data.frame(
    Phylostratum = phylo_mapping$phylostratum,
    GeneID = common_genes,
    filtered_expression,
    stringsAsFactors = FALSE
  )
  
  # Validate the set
  if (myTAI::is.ExpressionSet(phylo_expression_set)) {
    message("Successfully created PhyloExpressionSet with ", nrow(phylo_expression_set), " genes")
    return(phylo_expression_set)
  } else {
    stop("Failed to create valid PhyloExpressionSet")
  }
}

#' Create phylostrata mapping for marine organisms
#' @param gene_annotations Data frame with gene annotations
#' @param organism_name Character. Name of the marine organism
#' @param reference_phylogeny Character. Reference phylogenetic framework
#' @return Data frame with gene-to-phylostratum mapping
create_marine_phylostrata <- function(gene_annotations,
                                     organism_name = "marine_organism",
                                     reference_phylogeny = "eukaryota") {
  
  # Define basic phylostratigraphic levels for marine organisms
  # This reference can be used for documentation and validation
  marine_phylostrata <- list(
    "1" = "Cellular_organisms",
    "2" = "Eukaryota",
    "3" = "Opisthokonta",
    "4" = "Metazoa",
    "5" = "Eumetazoa",
    "6" = "Bilateria",
    "7" = "Protostomia_or_Deuterostomia",
    "8" = "Phylum_specific",
    "9" = "Class_specific", 
    "10" = "Order_specific",
    "11" = "Family_specific",
    "12" = "Genus_specific",
    "13" = "Species_specific"
  )
  
  message("Using ", length(marine_phylostrata), " phylostratigraphic levels for ", organism_name)
  
  # Create phylostratum assignment based on gene annotations
  # This is a simplified example - in practice, you would use orthology databases
  
  phylostrata_assignment <- data.frame(
    gene_id = gene_annotations$gene_id,
    phylostratum = NA,
    stringsAsFactors = FALSE
  )
  
  # Assign phylostrata based on gene descriptions/annotations
  for (i in seq_len(nrow(phylostrata_assignment))) {
    gene_desc <- tolower(gene_annotations$description[i])
    
    # Assign based on keywords in gene descriptions
    if (grepl("ribosomal|ribosome|rpl|rps|rrna", gene_desc)) {
      phylostrata_assignment$phylostratum[i] <- 1  # Ancient cellular machinery
    } else if (grepl("histone|chromatin|nucleus", gene_desc)) {
      phylostrata_assignment$phylostratum[i] <- 2  # Eukaryotic
    } else if (grepl("actin|tubulin|cytoskeleton", gene_desc)) {
      phylostrata_assignment$phylostratum[i] <- 2  # Eukaryotic cytoskeleton
    } else if (grepl("collagen|extracellular matrix", gene_desc)) {
      phylostrata_assignment$phylostratum[i] <- 4  # Metazoan innovation
    } else if (grepl("nervous|neural|neuron", gene_desc)) {
      phylostrata_assignment$phylostratum[i] <- 5  # Neural system
    } else if (grepl("immune|defense|antimicrobial", gene_desc)) {
      phylostrata_assignment$phylostratum[i] <- 6  # Immune system
    } else if (grepl("species.specific|unique|novel", gene_desc)) {
      phylostrata_assignment$phylostratum[i] <- 13  # Species-specific
    } else {
      # Default assignment for unclassified genes
      phylostrata_assignment$phylostratum[i] <- 8  # Phylum-specific (conservative)
    }
  }
  
  # Add metadata
  phylostrata_assignment$organism <- organism_name
  phylostrata_assignment$reference_phylogeny <- reference_phylogeny
  phylostrata_assignment$assignment_method <- "annotation_based"
  phylostrata_assignment$assignment_date <- Sys.Date()
  
  return(phylostrata_assignment)
}

#' Calculate transcriptome age index (TAI) for marine development
#' @param phylo_expression_set PhyloExpressionSet from myTAI
#' @param developmental_stages Character vector of developmental stage names
#' @return Data frame with TAI values across developmental stages
calculate_marine_tai <- function(phylo_expression_set, 
                                developmental_stages = NULL) {
  
  if (!requireNamespace("myTAI", quietly = TRUE)) {
    stop("Package 'myTAI' is required for TAI calculation.")
  }
  
  # Validate input
  if (!myTAI::is.ExpressionSet(phylo_expression_set)) {
    stop("Input must be a valid PhyloExpressionSet")
  }
  
  # Set developmental stage names if not provided
  if (is.null(developmental_stages)) {
    n_stages <- ncol(phylo_expression_set) - 2  # Subtract Phylostratum and GeneID columns
    developmental_stages <- paste0("Stage_", seq_len(n_stages))
  }
  
  # Calculate TAI using myTAI
  tai_values <- myTAI::TAI(phylo_expression_set)
  
  # Create result data frame
  tai_result <- data.frame(
    developmental_stage = developmental_stages,
    TAI = tai_values,
    stringsAsFactors = FALSE
  )
  
  # Add summary statistics
  tai_result$TAI_normalized <- (tai_result$TAI - min(tai_result$TAI)) / 
                             (max(tai_result$TAI) - min(tai_result$TAI))
  
  return(tai_result)
}

#' Perform phylostratigraphic analysis for marine conservation
#' @param phylo_expression_set PhyloExpressionSet with marine expression data
#' @param comparison_groups Character vector specifying groups for comparison
#' @param analysis_type Character. Type of analysis ("development", "stress", "comparative")
#' @return List with phylostratigraphic analysis results
analyze_marine_phylostratigraphy <- function(phylo_expression_set,
                                           comparison_groups = NULL,
                                           analysis_type = "development") {
  
  if (!requireNamespace("myTAI", quietly = TRUE)) {
    stop("Package 'myTAI' is required for phylostratigraphic analysis.")
  }
  
  results <- list()
  
  # Basic TAI calculation
  tai_values <- myTAI::TAI(phylo_expression_set)
  results$tai <- data.frame(
    sample = names(tai_values),
    TAI = tai_values,
    stringsAsFactors = FALSE
  )
  
  # Phylostratum contribution analysis
  phylo_contrib <- myTAI::PlotContribution(phylo_expression_set, 
                                          plot.it = FALSE)
  results$phylostratum_contribution <- phylo_contrib
  
  # Age-based expression patterns
  if (analysis_type == "development") {
    # Developmental analysis
    if (ncol(phylo_expression_set) > 4) {  # At least 3 developmental stages
      
      # Test for hourglass pattern
      hourglass_test <- tryCatch({
        myTAI::FlatLineTest(phylo_expression_set)
      }, error = function(e) {
        warning("Hourglass test failed: ", e$message)
        NULL
      })
      
      if (!is.null(hourglass_test)) {
        results$hourglass_test <- hourglass_test
      }
      
      # Early conservation test
      early_conservation <- tryCatch({
        myTAI::EarlyConservationTest(phylo_expression_set)
      }, error = function(e) {
        warning("Early conservation test failed: ", e$message)
        NULL
      })
      
      if (!is.null(early_conservation)) {
        results$early_conservation_test <- early_conservation
      }
    }
    
  } else if (analysis_type == "stress") {
    # Stress response analysis
    if (!is.null(comparison_groups) && length(comparison_groups) == ncol(phylo_expression_set) - 2) {
      
      # Group samples by treatment
      control_samples <- which(comparison_groups == "control")
      stress_samples <- which(comparison_groups == "stress")
      
      if (length(control_samples) > 0 && length(stress_samples) > 0) {
        # Calculate TAI difference between control and stress
        control_tai <- mean(tai_values[control_samples])
        stress_tai <- mean(tai_values[stress_samples])
        
        results$stress_response <- list(
          control_TAI = control_tai,
          stress_TAI = stress_tai,
          TAI_difference = stress_tai - control_tai,
          interpretation = if (stress_tai > control_tai) 
            "Stress activates younger genes" else "Stress activates older genes"
        )
        
        # Phylostratum-specific stress response
        phylo_stress_response <- list()
        
        for (ps in unique(phylo_expression_set$Phylostratum)) {
          ps_genes <- phylo_expression_set$Phylostratum == ps
          ps_expression <- phylo_expression_set[ps_genes, 3:ncol(phylo_expression_set)]
          
          if (nrow(ps_expression) > 0) {
            control_expr <- rowMeans(ps_expression[, control_samples, drop = FALSE])
            stress_expr <- rowMeans(ps_expression[, stress_samples, drop = FALSE])
            
            phylo_stress_response[[paste0("PS", ps)]] <- list(
              n_genes = nrow(ps_expression),
              mean_control = mean(control_expr),
              mean_stress = mean(stress_expr),
              log2_fold_change = log2(mean(stress_expr) / mean(control_expr))
            )
          }
        }
        
        results$phylostratum_stress_response <- phylo_stress_response
      }
    }
  }
  
  # Expression variance by phylostratum
  phylo_variance <- list()
  
  for (ps in unique(phylo_expression_set$Phylostratum)) {
    ps_genes <- phylo_expression_set$Phylostratum == ps
    ps_expression <- phylo_expression_set[ps_genes, 3:ncol(phylo_expression_set)]
    
    if (nrow(ps_expression) > 0) {
      phylo_variance[[paste0("PS", ps)]] <- list(
        n_genes = nrow(ps_expression),
        mean_expression = mean(as.matrix(ps_expression)),
        expression_variance = var(as.matrix(ps_expression)),
        cv = sd(as.matrix(ps_expression)) / mean(as.matrix(ps_expression))
      )
    }
  }
  
  results$phylostratum_variance <- phylo_variance
  
  return(results)
}

#' Compare evolutionary patterns between marine species
#' @param expression_sets List of PhyloExpressionSets for different species
#' @param species_names Character vector of species names
#' @return List with comparative phylostratigraphic analysis
compare_marine_species_evolution <- function(expression_sets,
                                           species_names) {
  
  if (!requireNamespace("myTAI", quietly = TRUE)) {
    stop("Package 'myTAI' is required for comparative analysis.")
  }
  
  if (length(expression_sets) != length(species_names)) {
    stop("Number of expression sets must match number of species names")
  }
  
  comparison_results <- list()
  
  # Calculate TAI for each species
  species_tai <- list()
  
  for (i in seq_along(expression_sets)) {
    species_name <- species_names[i]
    expression_set <- expression_sets[[i]]
    
    if (myTAI::is.ExpressionSet(expression_set)) {
      tai_values <- myTAI::TAI(expression_set)
      species_tai[[species_name]] <- tai_values
    } else {
      warning("Invalid expression set for species: ", species_name)
    }
  }
  
  comparison_results$species_tai <- species_tai
  
  # Compare phylostratum usage between species
  phylostratum_usage <- list()
  
  for (i in seq_along(expression_sets)) {
    species_name <- species_names[i]
    expression_set <- expression_sets[[i]]
    
    if (myTAI::is.ExpressionSet(expression_set)) {
      ps_counts <- table(expression_set$Phylostratum)
      phylostratum_usage[[species_name]] <- ps_counts
    }
  }
  
  comparison_results$phylostratum_usage <- phylostratum_usage
  
  # Calculate evolutionary conservation scores
  if (length(species_tai) >= 2) {
    # Pairwise TAI correlations
    tai_correlations <- matrix(NA, 
                              nrow = length(species_tai), 
                              ncol = length(species_tai),
                              dimnames = list(names(species_tai), names(species_tai)))
    
    for (i in seq_along(species_tai)) {
      for (j in seq_along(species_tai)) {
        if (i != j) {
          tai_i <- species_tai[[i]]
          tai_j <- species_tai[[j]]
          
          # Ensure same length for correlation
          min_length <- min(length(tai_i), length(tai_j))
          if (min_length > 1) {
            correlation <- cor(tai_i[1:min_length], tai_j[1:min_length])
            tai_correlations[i, j] <- correlation
          }
        }
      }
    }
    
    diag(tai_correlations) <- 1
    comparison_results$tai_correlations <- tai_correlations
  }
  
  return(comparison_results)
}

#' Generate conservation-focused evolutionary insights
#' @param phylo_analysis_results Results from phylostratigraphic analysis
#' @param species_name Character. Name of the marine species
#' @param conservation_context Character. Conservation context ("endangered", "invasive", "commercial")
#' @return List with conservation-relevant evolutionary insights
generate_conservation_insights <- function(phylo_analysis_results,
                                         species_name,
                                         conservation_context = "endangered") {
  
  insights <- list()
  insights$species <- species_name
  insights$conservation_context <- conservation_context
  insights$analysis_date <- Sys.Date()
  
  # Extract TAI patterns
  if ("tai" %in% names(phylo_analysis_results)) {
    tai_data <- phylo_analysis_results$tai
    
    insights$tai_summary <- list(
      mean_tai = mean(tai_data$TAI),
      tai_range = range(tai_data$TAI),
      tai_variance = var(tai_data$TAI)
    )
    
    # Conservation interpretations based on TAI patterns
    if (conservation_context == "endangered") {
      if (mean(tai_data$TAI) < 5) {
        insights$conservation_interpretation <- paste(
          "Low TAI suggests reliance on ancient, conserved genes.",
          "This may indicate limited adaptive potential to environmental changes."
        )
      } else {
        insights$conservation_interpretation <- paste(
          "Higher TAI suggests active use of newer genes.",
          "This may indicate better adaptive potential but also higher vulnerability to genetic drift."
        )
      }
      
    } else if (conservation_context == "invasive") {
      if (var(tai_data$TAI) > 1) {
        insights$conservation_interpretation <- paste(
          "High TAI variance suggests flexible gene expression.",
          "This may contribute to invasive success across different environments."
        )
      }
      
    } else if (conservation_context == "commercial") {
      insights$conservation_interpretation <- paste(
        "TAI patterns can inform selective breeding programs.",
        "Ancient genes (low phylostrata) may be targets for maintaining robustness."
      )
    }
  }
  
  # Phylostratum-specific insights
  if ("phylostratum_variance" %in% names(phylo_analysis_results)) {
    ps_variance <- phylo_analysis_results$phylostratum_variance
    
    # Find most variable phylostrata
    variance_values <- sapply(ps_variance, function(x) x$expression_variance)
    most_variable_ps <- names(variance_values)[which.max(variance_values)]
    
    insights$most_variable_phylostratum <- most_variable_ps
    insights$phylostratum_insights <- paste(
      "Phylostratum", most_variable_ps, "shows highest expression variance.",
      "This may represent genes under selection or environmental responsiveness."
    )
  }
  
  # Stress response insights
  if ("stress_response" %in% names(phylo_analysis_results)) {
    stress_data <- phylo_analysis_results$stress_response
    
    insights$stress_insights <- list(
      tai_stress_direction = if (stress_data$TAI_difference > 0) "younger" else "older",
      stress_interpretation = stress_data$interpretation,
      adaptive_capacity = if (abs(stress_data$TAI_difference) > 0.5) "high" else "low"
    )
  }
  
  # Conservation recommendations
  recommendations <- character()
  
  if (conservation_context == "endangered") {
    recommendations <- c(
      "Monitor expression of ancient genes (PS1-3) as indicators of physiological stress",
      "Assess genetic diversity in rapidly evolving gene families (high PS numbers)",
      "Consider captive breeding programs if TAI variance is low (reduced adaptability)"
    )
    
  } else if (conservation_context == "invasive") {
    recommendations <- c(
      "Target management efforts during life stages with high TAI (active gene expression)",
      "Monitor for rapid evolution in species-specific genes (PS13)",
      "Consider evolutionary rescue in native range if TAI patterns are disrupted"
    )
    
  } else if (conservation_context == "commercial") {
    recommendations <- c(
      "Select breeding lines with stable expression of ancient genes for robustness",
      "Monitor TAI changes in response to aquaculture conditions",
      "Develop stress-resistant strains by targeting appropriate phylostrata"
    )
  }
  
  insights$conservation_recommendations <- recommendations
  
  return(insights)
}