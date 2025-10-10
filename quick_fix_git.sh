#!/bin/bash
# Complete repository cleanup script

echo "🧹 Repository Deep Cleanup"
echo "=========================="
echo ""
echo "This script will:"
echo "  1. Create a backup of your repository"
echo "  2. Remove large CSV files from Git history"
echo "  3. Create proper .gitignore entries"
echo "  4. Set up a clean repository ready for GitHub"
echo ""

read -p "Continue? (y/n): " -n 1 -r
echo ""
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Operation cancelled."
    exit 0
fi

# Step 1: Create backup
BACKUP_DIR="../biologia-marina-backup-$(date +%Y%m%d-%H%M%S)"
echo "📦 Creating backup in $BACKUP_DIR..."
cp -r . "$BACKUP_DIR"
echo "✓ Backup created"
echo ""

# Step 2: Remove specific large files using git filter-repo
echo "🔧 Removing large files from history..."

# Ensure git-filter-repo is available
if ! command -v git-filter-repo &> /dev/null; then
    echo "⚠️ git-filter-repo not found. Installing..."
    pip install --user git-filter-repo
    
    if ! command -v git-filter-repo &> /dev/null; then
        echo "❌ Failed to install git-filter-repo. Try adding it manually:"
        echo "    pip install --user git-filter-repo"
        exit 1
    fi
fi

# Add explicit patterns for CSV files
echo "Removing large CSV files..."
git filter-repo --path 'data/raw/obis/orites_astreoides_obis_data.csv' --invert-paths --force
git filter-repo --path 'data/raw/acropora_data.csv' --invert-paths --force

# Clean additional data files
echo "Removing other large data files..."
git filter-repo --path-glob 'data/raw/*.csv' --invert-paths --force
git filter-repo --path-glob '*.rds' --invert-paths --force
git filter-repo --path '_targets/objects' --invert-paths --force
git filter-repo --path-glob '*.RData' --invert-paths --force

# Step 3: Update .gitignore
echo "📄 Creating proper .gitignore..."
cat > .gitignore << 'GITIGNORE'
# R and RStudio files
.Rproj.user/
.Rhistory
.RData
.Ruserdata
*.Rproj

# targets
_targets/
!_targets.R

# Large data files
*.rds
*.RData
data/raw/*.csv
data/raw/obis/*.csv
data/processed/*.csv
data/raw/environmental/*.tif
data/raw/environmental/*.asc

# Nix
result
result-*
*.backup

# System
.DS_Store
*.swp
*~
GITIGNORE
git add .gitignore
git commit -m "Update .gitignore to exclude large data files"

# Step 4: Update build_env.R to include git-filter-repo
echo "📝 Updating build_env.R..."
sed -i 's/"git",/"git", "git-filter-repo",/' build_env.R
git add build_env.R
git commit -m "Add git-filter-repo to system packages"

# Step 5: Clean and optimize repo
echo "🧹 Optimizing repository..."
git gc --aggressive --prune=now

# Step 6: Instructions
echo ""
echo "✅ Repository cleaned successfully!"
echo ""
echo "Repository size: $(du -sh .git | cut -f1)"
echo ""
echo "To push to GitHub:"
echo "  1. Create a new empty repository on GitHub (don't initialize it)"
echo "  2. Run these commands:"
echo ""
echo "     git remote add origin https://github.com/YOUR-USERNAME/biologia-marina-reproducible.git"
echo "     git push -u --force origin main"
echo ""
echo "If issues persist, consider a fresh start with a new repository:"
echo "  1. Create a new empty directory"
echo "  2. Copy all your current files (except .git folder)"
echo "  3. Initialize a new Git repository there"
echo ""