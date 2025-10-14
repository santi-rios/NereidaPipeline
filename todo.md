## Files

- # TODO: consider using git-lfs for database

```bash
# Install Git LFS
sudo apt-get install git-lfs  # For Ubuntu/Debian
# or
brew install git-lfs  # For macOS

# Initialize Git LFS
git lfs install

# Track the DuckDB file with LFS
git lfs track "data/marine_biodiversity.duckdb"

# Create a .gitattributes file (or add to existing)
git add .gitattributes
```


pre-commit hook (optional)


# Create .git/hooks/pre-commit
cat > .git/hooks/pre-commit << 'EOF'
#!/bin/bash

# Find files larger than 50MB
large_files=$(find . -type f -size +50M -not -path "*.git*" | grep -v ".gitignore")

if [ -n "$large_files" ]; then
  echo "Error: Attempting to commit large files:"
  echo "$large_files"
  echo "Please remove these files or add to .gitignore"
  exit 1
fi
EOF

chmod +x .git/hooks/pre-commit