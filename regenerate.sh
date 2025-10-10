#!/bin/bash

# Simple script to regenerate Nix environment
# Just runs rix with the configuration from build_env.R

echo "🔧 Generating Nix Environment..."
echo "================================="

# Use rix from nixpkgs to run build_env.R
nix-shell -p R rPackages.rix --run "R --vanilla -e 'source(\"build_env.R\")'"

if [ $? -eq 0 ]; then
    echo ""
    echo "✅ Environment generated successfully!"
    echo ""
    echo "Next steps:"
    echo "  1. Enter environment: nix-shell"
    echo "  2. Test packages: R -e 'library(targets)'"
    echo ""
else
    echo ""
    echo "❌ Environment generation failed"
    echo "Check the error messages above"
    exit 1
fi
