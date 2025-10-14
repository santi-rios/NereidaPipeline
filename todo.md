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