# Enhanced Marine Biodiversity Research Environment
# Build script for generating Nix environment with rix
#
# IMPORTANT: This script requires the rix package.
# 
# To run this script, use one of these methods:
#
# Method 1 - Use rix from nixpkgs:
#   nix-shell -p R rPackages.rix --run "R --vanilla -e 'source(\"build_env.R\")'"
#
# Method 2 - Use rix development version:
#   nix-shell --expr "$(curl -sl https://raw.githubusercontent.com/ropensci/rix/master/inst/extdata/default.nix)" --run "R --vanilla -e 'source(\"build_env.R\")'"
#
# Method 3 - Use the helper scripts:
#   ./generate_nix_env.sh       # Uses stable rix from nixpkgs
#   ./generate_nix_env_dev.sh   # Uses latest rix from GitHub

# Check if rix is available
if (!requireNamespace("rix", quietly = TRUE)) {
  stop(
    "❌ The 'rix' package is not available.\n\n",
    "Please run this script using one of these methods:\n\n",
    "1. nix-shell -p R rPackages.rix --run \"R --vanilla -e 'source(\\\"build_env.R\\\")'\" \n\n",
    "2. ./generate_nix_env.sh\n\n",
    "3. ./generate_nix_env_dev.sh\n"
  )
}

library(rix)

# Core packages always needed
core_packages <- c(
  "targets",
  "tarchetypes", 
  "dplyr",
  "quarto",
  "languageserver"  # For VSCode
)

# Packages from _targets.R (update this list when you modify _targets.R)
targets_packages <- c(
  "robis",
  "rgbif",
  "CoordinateCleaner",
  "taxa",
  "terra",
  "nodbi",
  ## bioconductor
  "BiocManager",
  "biomaRt",
  # "BarcodingR",
  "wikitaxa",
  "ggplot2",
  "plotly",
  "jsonlite",
  "DBI",
  "duckdb"
)

# Optional packages (only if available)
optional_packages <- c(
  "biomaRt",
  # "wikitaxa",
  "prism",
  "myTAI",
  "geotargets"
)

# Combine all packages
all_packages <- unique(c(core_packages, targets_packages))

cat("Generating Nix environment with packages:\n")
cat(paste(" -", all_packages, collapse = "\n"), "\n\n")

rix(
  r_ver = "latest-upstream",
  r_pkgs = all_packages,
  system_pkgs = c(
    "git", 
    # "git-filter-repo",
    "python3",
    "quarto"
    # "bfg-repo-cleaner" # alternative to git-filter-branch
  ),
  tex_pkgs = c("amsmath"),
  ide = "none",
  shell_hook = "",
  project_path = ".",
  overwrite = TRUE,
  print = TRUE
)

cat("\n✅ Environment configuration generated!\n")
cat("Next: Run './update_workflow.sh' to build and test\n")