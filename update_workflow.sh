#!/bin/bash
# Streamlined workflow for updating Nix environment and testing

set -e  # Exit on error

echo "🔄 Nix Environment Update Workflow"
echo "=================================="
echo ""

# Step 1: Show current packages (simple version)
echo "Step 1: Verificando paquetes en _targets.R..."
echo ""
echo "Los paquetes se verificarán automáticamente en el paso 5."
echo ""

# Step 2: Ask if build_env.R needs updating
echo "Step 2: Checking build_env.R..."
if [ -f "build_env.R" ]; then
    echo "✅ build_env.R found"
    echo ""
    echo "Make sure these packages are in build_env.R:"
    echo "  $TARGETS_PACKAGES"
    echo ""
    read -p "Have you updated build_env.R? (y/n): " -n 1 -r
    echo ""
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "ℹ️  Please update build_env.R with required packages and run this script again"
        exit 1
    fi
else
    echo "❌ build_env.R not found"
    exit 1
fi

# Step 3: Regenerate Nix environment
echo ""
echo "Step 3: Regenerating Nix environment..."
./regenerate.sh

if [ ! -f "default.nix" ]; then
    echo "❌ Failed to generate default.nix"
    exit 1
fi

# Step 4: Build the environment
echo ""
echo "Step 4: Building Nix environment..."
nix-build

if [ $? -ne 0 ]; then
    echo "❌ Environment build failed"
    exit 1
fi

# Step 5: Test packages
echo ""
echo "Step 5: Testing package availability..."
./test_environment.sh

# Step 6: Optional - Run targets pipeline
echo ""
read -p "Run targets pipeline now? (y/n): " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "🎯 Running targets pipeline..."
    nix-shell --run "Rscript -e 'targets::tar_make()'"
fi

echo ""
echo "✅ Workflow complete!"