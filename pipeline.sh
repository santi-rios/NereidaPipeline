#!/bin/bash
# run pipeline

echo ""
echo "Running pipeline"
nix-shell --run "Rscript -e 'targets::tar_make()'"