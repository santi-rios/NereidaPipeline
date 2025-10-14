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

```bash
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
```

## quarto dirs

---https://github.com/quarto-dev/quarto-cli/discussions/10535#discussioncomment-10626220

> I'll try to clarify a few things
> 
> > Quarto projects don't seem to be completely self-contained when using quarto::render(execute_dir = ".") though and that's confusing.
> 
> `execute_dir = "."` in the R function is equivalent to calling `quarto render <file> --execute-dir .`. Quarto will normalize the path `.` from where the `quarto render` command is called. So in here,
> 
> ```r
> quarto::quarto_render(
>   input = "manuscript/index.qmd",
>   execute_dir = "."
> )
> ```
> 
> it will be from root of you RStudio project. If you call that from another folder it would be different. Using `.` is not the safest to insure correct path is used.
> 
> For knitr engine, `--execute-dir` will set the value for `knit_root_dir` argument in `rmarkdown::render()`, which is setting `knitr::opts_knit$set(root.dir = )`. So it really controls what R engine is doing.
> 
> > You can change the execute directory via Quarto CLI: [quarto.org/docs/projects/code-execution.html#working-dir](https://quarto.org/docs/projects/code-execution.html#working-dir)
> > (...)
> > Unless I'm missing something, this will set the working directory to the Quarto project, not the R one.
> 
> If no `--execute-dir` is passed at command line, then `execute-dir` option in project YAML will be used to set the value to the project path when `execute-dir: project`. Otherwise, the directory of the source file being rendered is used.
> 
> > After a bit more testing, it seems that the only issue is with .Rprofile. If I copy the .Rprofile into ms_project/, the document will render properly. That means the executing directory is the root of the R project in the end. Although it's executed, the .Rprofile at the root of the R project seems to have no effect.
> 
> I'll check with your exemple but here is how it is supposed to work.
> 
> Quarto will call R from working dir where `quarto render` is called, unless there is a project, and in which case the project dir will be used.
> 
> In your case, the project is inside `ms_project`. So when you call `quarto render ms_project/index.qmd`, it will find this is a project and so call R from there. An `.Rprofile` there will be loaded. Doing `quarto render ms_standalone/index.qmd` is a single file render. In this case, there is no project so the working dir where `quarto render` is called will be used, i.e. your project root.
> 
> > I'll use the pre-render and post-render options to copy the .Rprofile to the Quarto project directory and delete it when the rendering is complete for now.
> 
> This seems like a good workflow for your compendium organization. You could also see if having a .Rprofile in your subfolder can source a script in the root folder.
> 
> Other comment about .Rprofile. The file used can be specified by `R_PROFILE_USER` environment variable. ([from man page](https://stat.ethz.ch/R-manual/R-devel/library/base/html/Startup.html)). So maybe there is way to tell Quarto to always load a specific one no matter where is called. You would have to deal with path resolution to make that safe.
> 
> A related idea for you specific project organization. Quarto project inside your compendium overall structure does not know it is part of such compendium structure. You may set some environment variable or other config to tell Quarto project absolute path to your compendium root maybe this could be useful. Or at least for R code you could leverage **here** package probably.
> 
> Quite long answer, but I hope this will help make your workflow more understandable regarding Quarto's behavior.



---