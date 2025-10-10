# Marine Biodiversity Research Environment Setup

## Overview
This document provides comprehensive instructions for setting up the enhanced marine biodiversity research environment using Nix and the `rix` R package.

## Prerequisites

### 1. Install Nix Package Manager
```bash
# Install Nix (multi-user installation recommended for Linux)
curl -L https://nixos.org/nix/install | sh -s -- --daemon

# Source the Nix configuration
source /etc/profile.d/nix.sh

# Verify installation
nix --version
```

### 2. Enable Nix Flakes (Required for rix)
```bash
# Add experimental features to Nix configuration
mkdir -p ~/.config/nix
echo "experimental-features = nix-command flakes" > ~/.config/nix/nix.conf

# Or edit existing configuration
echo "experimental-features = nix-command flakes" >> ~/.config/nix/nix.conf
```

## Environment Setup

### 1. Generate Nix Environment
```bash
# Navigate to project directory
cd /home/santi/Projects/biologia-marina-reproducible

# Generate the Nix environment using the enhanced build_env.R
R --vanilla -e "source('build_env.R')"
```

This will create:
- `default.nix`: Complete environment specification
- Enhanced package dependencies for all pipeline components

### 2. Enter the Nix Environment
```bash
# Enter the Nix development shell
nix-shell

# Alternative: Use direnv for automatic environment loading
echo "use nix" > .envrc
direnv allow
```

### 3. Verify Package Installation
Once inside the Nix shell, verify that all packages are available:

```r
# Start R and check package availability
R

# Core packages
library(targets)
library(tarchetypes)
library(dplyr)

# Biodiversity data sources
library(robis)
library(rgbif)
library(biomaRt)      # May need BiocManager installation
library(wikitaxa)

# Spatial and environmental
library(terra)
library(sf)
library(prism)
library(geotargets)

# Data cleaning
library(CoordinateCleaner)

# Taxonomic management
library(taxa)

# Database integration
library(nodbi)
library(DBI)
library(duckdb)

# Evolutionary analysis
library(myTAI)         # May need BiocManager installation

# Visualization
library(ggplot2)
library(plotly)
library(leaflet)

# Utilities
library(jsonlite)
library(quarto)
```

## Handling Bioconductor Packages

Some packages (like `biomaRt` and `myTAI`) are from Bioconductor and may need special handling:

### Method 1: Install within Nix environment
```r
# Inside the Nix R session
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("biomaRt", "myTAI"))
```

### Method 2: Update build_env.R for Bioconductor
If needed, you can specify Bioconductor packages explicitly:

```r
# Add to build_env.R r_pkgs section:
"BiocManager",
# Then install Bioconductor packages manually
```

## Python Integration Setup

The environment includes Python3 for machine learning components:

```bash
# Verify Python installation
python3 --version

# Install Python packages for the notebook analysis
pip install --user pandas numpy matplotlib seaborn scikit-learn plotly jupyter
```

## Environment Activation Workflow

### Daily Usage
```bash
# 1. Navigate to project
cd /home/santi/Projects/biologia-marina-reproducible

# 2. Enter Nix environment
nix-shell

# 3. Verify environment (first time)
echo $NIX_PATH

# 4. Start analysis
R  # For R analysis
# or
python3  # For Python ML components
# or
quarto render reportes/  # For report generation
```

### Automated Activation with direnv (Recommended)
```bash
# Install direnv
nix-env -iA nixpkgs.direnv

# Setup in project directory
echo "use nix" > .envrc
direnv allow

# Now the environment activates automatically when entering the directory
```

## Troubleshooting

### Common Issues and Solutions

#### 1. Package Not Found
```bash
# Update Nix channels
nix-channel --update

# Regenerate environment
R --vanilla -e "source('build_env.R')"
nix-shell
```

#### 2. Bioconductor Package Issues
```r
# Force Bioconductor installation
BiocManager::install(c("biomaRt", "myTAI"), force = TRUE)

# Check Bioconductor version compatibility
BiocManager::valid()
```

#### 3. System Dependencies Missing
Some packages may need additional system libraries:

```bash
# For spatial packages (already included in system_pkgs)
# gdal, proj, geos should be available

# For database packages
# sqlite, included in system_pkgs

# Verify system packages are available
which gdal-config
which proj
which sqlite3
```

#### 4. Memory Issues with Large Datasets
```bash
# Increase memory limits in R
R --vanilla --max-mem-size=8G

# Or set environment variable
export R_MAX_MEM_SIZE=8G
```

## Performance Optimization

### 1. Parallel Processing
```r
# Enable parallel targets execution
# In _targets.R, add:
tar_option_set(
  packages = c("targets", "tarchetypes", /* all your packages */),
  format = "rds",
  workspace_on_error = TRUE,
  memory = "transient",
  garbage_collection = TRUE
)

# Use parallel workers
library(future)
plan(multisession, workers = 4)  # Adjust based on your CPU cores
```

### 2. Caching Strategy
```bash
# Enable Nix binary cache for faster builds
echo "substituters = https://cache.nixos.org/ https://nix-community.cachix.org" >> ~/.config/nix/nix.conf
echo "trusted-public-keys = cache.nixos.org-1:6NCHdD59X431o0gWypbMrAURkbJ16ZPMQFGspcDShjY= nix-community.cachix.org-1:mB9FSh9qf2dCimDSUo8Zy7bkq5CX+/rkCWyvRCYg3Fs=" >> ~/.config/nix/nix.conf
```

## Integration Testing

### Test the Complete Pipeline
```bash
# Enter environment
nix-shell

# Run pipeline test
R --vanilla -e "
library(targets)
tar_make()  # Execute the complete pipeline
tar_visnetwork()  # Visualize pipeline dependencies
"

# Test Python integration
python3 notebooks/analisis_python_integracion.py  # If converted from .ipynb

# Test Quarto reports
quarto render reportes/analisis_biodiversidad_marina.qmd
```

## Maintenance

### Updating the Environment
```bash
# Update Nix channels
nix-channel --update

# Regenerate environment with latest packages
R --vanilla -e "source('build_env.R')"

# Clean old generations (optional, saves disk space)
nix-collect-garbage -d
```

### Version Pinning (for Reproducibility)
If you need exact reproducibility, consider pinning package versions:

```r
# In build_env.R, specify exact R version
r_ver = "4.3.2"  # Instead of "latest-upstream"

# Or use specific package versions
# This requires more advanced Nix configuration
```

## Environment Variables

The enhanced environment sets useful variables:

```bash
# Check environment variables after entering nix-shell
echo $R_LIBS_SITE
echo $PYTHONPATH
echo $PATH
```

## Support Resources

- **Nix Documentation**: https://nixos.org/manual/nix/stable/
- **rix Package Documentation**: https://b-rodrigues.github.io/rix/
- **Targets Package**: https://books.ropensci.org/targets/
- **Marine Biodiversity Packages**:
  - robis: https://robis.github.io/
  - rgbif: https://docs.ropensci.org/rgbif/
  - biomaRt: https://docs.ropensci.org/biomaRt/

## Quick Reference Card

```bash
# Essential commands:
nix-shell                    # Enter environment
R                           # Start R with all packages
python3                     # Start Python
quarto render              # Generate reports
tar_make()                 # Run targets pipeline (in R)
tar_visnetwork()           # Visualize pipeline (in R)
exit                       # Leave nix-shell
```

This setup provides a completely reproducible environment for your enhanced marine biodiversity research pipeline with all dependencies properly managed through Nix.