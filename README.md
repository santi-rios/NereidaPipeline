### Project Repository Structure

```
marine-biodiversity-project/
│
├── .github/
│   └── workflows/
│       └── ci.yml                # GitHub Actions workflow for CI/CD
│
├── data/
│   ├── raw/                       # Raw data files
│   └── processed/                 # Processed data files
│
├── R/
│   ├── analysis.R                 # Main analysis script
│   ├── data_acquisition.R         # Script for data acquisition
│   ├── data_visualization.R        # Script for data visualization
│   └── machine_learning.py         # Python script for ML applications
│
├── reports/
│   ├── technical_report.Rmd        # Technical report in R Markdown
│   └── scientific_uses_es.md       # Scientific uses written in Spanish
│
├── scripts/
│   ├── create_env.R                # Script to create Nix environment
│   └── run_pipeline.R              # Script to run the targets pipeline
│
├── default.nix                     # Nix environment configuration
├── README.md                       # Project overview and instructions
└── targets.R                       # Targets pipeline configuration
```

### Step-by-Step Instructions

1. **Initialize the Repository**:
   - Create a new GitHub repository named `marine-biodiversity-project`.
   - Clone the repository to your local machine.

InstaLLING niX

```bash
curl --proto '=https' --tlsv1.2 -sSf \
    -L https://install.determinate.systems/nix | \
     sh -s -- install
```

Then, install the cachix client and configure our rstats-on-nix cache: this will install binary versions of many R packages which will speed up the building process of environments:

```bash
nix-env -iA cachix -f https://cachix.org/api/v1/install
```

then use the cache:

```bash
cachix use rstats-on-nix
```
Running the following line in a terminal will drop you in an interactive R session that you can use to start generating expressions:

```bash
nix-shell -p R rPackages.rix
```

or if you prefer the development version of rix:

```bash
nix-shell --expr "$(curl -sl https://raw.githubusercontent.com/ropensci/rix/master/inst/extdata/default.nix)"
```

This should immediately start an R session inside your terminal. You can now run something like this:

library(rix)

rix(
  r_ver = "4.4.2",
  r_pkgs = c("dplyr", "ggplot2"),
  system_pkgs = NULL,
  git_pkgs = NULL,
  ide = "none",
  project_path = ".",
  overwrite = TRUE
)

to generate a default.nix, and then use that file to generate an environment with R

nix-build

To start using this environment, open a terminal in the folder containing default.nix and use the following Nix command:

nix-build

nix-build is a Nix command that builds an environment according to the specifications found in a default.nix file. Once the environment is done building, you should find a new file called result next to the default.nix file.

To now use the environment, type in the same terminal as before:

nix-shell

This will activate the environment.

If you wish to run the pipeline whenever you drop into the Nix shell, you could add a Shell-hook to the generated default.nix file:

```r
path_default_nix <- tempdir()

rix(
  r_ver = "4.2.2",
  r_pkgs = c("targets", "tarchetypes", "rmarkdown"),
  system_pkgs = NULL,
  git_pkgs = list(
    package_name = "housing",
    repo_url = "https://github.com/rap4all/housing/",
    commit = "1c860959310b80e67c41f7bbdc3e84cef00df18e"
  ),
  ide = "none",
  shell_hook = "Rscript -e 'targets::tar_make()'",
  project_path = path_default_nix,
  overwrite = TRUE
)
```

Now, each time you drop into the Nix shell for that project using nix-shell, the pipeline gets automatically executed. rix also features a function called tar_nix_ga() that adds a GitHub Actions workflow file to make the pipeline run automatically on GitHub Actions. The GitHub repository linked above has such a file, so each time changes get pushed, the pipeline runs on GitHub Actions and the results are automatically pushed to a branch called targets-runs. 

We are now going to bootstrap this environment. First, run the following line to drop into a temporary shell with R and rix:

```r
# nix-shell --expr "$(curl -sl https://raw.githubusercontent.com/ropensci/rix/main/inst/extdata/default.nix)"
nix-shell --expr "$(curl -sl https://raw.githubusercontent.com/ropensci/rix/master/inst/extdata/default.nix)"

```

Once the build process is done, you can simply start R by typing R and then, within the R session, run source("gen_env.R"). 


```bash
santi@pop-os:~/Projects/NereidaPipeline$ nix-shell --expr "$(curl -sl https://raw.githubusercontent.com/ropensci/rix/master/inst/extdata/default.nix)"
warning: unknown setting 'eval-cores'
warning: unknown setting 'lazy-trees'
unpacking 'https://github.com/rstats-on-nix/nixpkgs/archive/2025-04-29.tar.gz' into the Git cache...
warning: ignoring untrusted substituter 'https://rstats-on-nix.cachix.org', you are not a trusted user.
Run `man nix.conf` for more information on the `substituters` configuration option.
warning: ignoring the client-specified setting 'trusted-public-keys', because it is a restricted setting and you are not a trusted user

[nix-shell:~/Projects/NereidaPipeline]$ nix-shell
warning: unknown setting 'eval-cores'
warning: unknown setting 'lazy-trees'
unpacking 'https://github.com/NixOS/nixpkgs/archive/3de11c4504efc3ac58c7bee62f9391e80f7b2910.tar.gz' into the Git cache...
warning: ignoring untrusted substituter 'https://rstats-on-nix.cachix.org', you are not a trusted user.
Run `man nix.conf` for more information on the `substituters` configuration option.
warning: ignoring the client-specified setting 'trusted-public-keys', because it is a restricted setting and you are not a trusted user
unpacking 'https://flakehub.com/f/DeterminateSystems/nixpkgs-weekly/%2A.tar.gz' into the Git cache...

[nix-shell:~/Projects/NereidaPipeline]$ make help

📚 Comandos Disponibles para Biología Marina
=============================================

  help                 Mostrar esta ayuda
  setup                Configuración inicial del proyecto
  permissions          Establecer permisos de scripts
  regenerate           Regenerar ambiente Nix desde build_env.R
  test                 Probar que todos los paquetes estén disponibles
  update               Flujo completo: regenerar + construir + probar
  git-fix              Arreglar archivos grandes en el historial de Git
  nix-warnings         Arreglar advertencias de configuración de Nix
  clean                Limpiar archivos temporales y objetos de targets


[nix-shell:~/Projects/NereidaPipeline]$ make setup
🔐 Configurando permisos de scripts...
✅ Permisos actualizados
✅ Proyecto configurado y listo para usar
```

--------

```bash
echo "trusted-users = root santi" | sudo tee -a /etc/nix/nix.conf && sudo pkill nix-daemon
cachix use rstats-on-nix
```

2. **Set Up Nix Environment**:
   - Create a `build_env.R` in the `root` directory to define the Nix environment using `r-rix`:



```r
library(rix)

rix(
  r_ver = "latest-upstream",
  r_pkgs = c("dplyr", "ggplot2", "targets", "lubridate", "raster", "robis", "rgbif"),
  system_pkgs = c("git", "python3"),
  ide = "none",
  project_path = ".",
  overwrite = TRUE,
  print = TRUE
)
```

Then use the following command to bootstrap an enivronment with R and rix only (from the same directory):

```bash
nix-shell -p R rPackages.rix
# nix-shell --expr "$(curl -sl https://raw.githubusercontent.com/ropensci/rix/main/inst/extdata/default.nix)" # developement version

nix-build
nix-shell
```

   - Run the script to generate the `default.nix` file.

```r
Rscript your_script.R
```

> Troubleshoot

Extra steps for VScode

- install a piece of software called direnv: direnv will automatically load Nix shells when you open a project that contains a default.nix file in an editor. 
- 


```bash
sudo apt install direnv
nix-env -f '<nixpkgs>' -iA direnv
nix-env -f '<nixpkgs>' -iA nix-direnv
```

then in code install "https://github.com/direnv/direnv-vscode"

Then, in the project root, create a '.envrc' file:

```bash
~/projects/my-nix-r-project/.envrc
```

and put:

```
use nix
mkdir $TMP
```

you will see a pop-up stating direnv: /PATH/TO/PROJECT/.envrc is blocked and a button to allow it.

hen open an R script. You might get another pop-up asking you to restart the extension, so click Restart. Be aware that at this point, direnv will run nix-shell and so will start building the environment. If that particular environment hasn’t been built and cached yet, it might take some time before Code will be able to interact with it. You might get yet another popup, this time from the R Code extension complaining that R can’t be found. In this case, simply restart Code and open the project folder again: now it should work every time. For a new project, simply repeat this process:

    Generate the project’s default.nix file;
    Build it using nix-build;
    Create an .envrc and write the two lines from above in it;
    Open the project’s folder in Code and click allow when prompted;
    Restart the extension and Code if necessary.

if you prefer to use the Code you have installed system-wide, then you need to add the languageserver package and use ide = "none":

```r
rix(
  date = ...,
  r_pkgs = c("languageserver", ...), # languageserver is needed
  ide = "none",
  ...
)
```

> ⚠️ si sale un mensaje en VScode de instalar languageserver, no insttalarlo por este medio. Declinar.
>
> 


1. **Create the Targets Pipeline**:
   - In `targets.R`, define the pipeline for data acquisition, processing, and analysis:
     ```r
     library(targets)

     # Load required packages
     library(robis)
     library(rgbif)
     library(dplyr)

     # Define targets
     list(
       tar_target(raw_data, get_obis_data()),
       tar_target(processed_data, process_data(raw_data)),
       tar_target(visualization, create_visualization(processed_data)),
       tar_target(model, run_machine_learning(processed_data))
     )
     ```

2. **Data Acquisition Script**:
   - In `data_acquisition.R`, write functions to fetch marine biodiversity data using `robis` and `rgbif`:
     ```r
     get_obis_data <- function() {
       library(robis)
       # Fetch data from OBIS
       obis_data <- occurrence(scientificname = "Acropora")
       return(obis_data)
     }
     ```

3. **Data Visualization Script**:
   - In `data_visualization.R`, create functions to visualize the data:
     ```r
     create_visualization <- function(data) {
       library(ggplot2)
       ggplot(data, aes(x = decimalLongitude, y = decimalLatitude)) +
         geom_point() +
         labs(title = "Distribución de Acropora")
     }
     ```

4. **Machine Learning Script**:
   - In `machine_learning.py`, implement Python code for machine learning applications:
     ```python
     import pandas as pd
     from sklearn.ensemble import RandomForestClassifier

     # Load processed data
     data = pd.read_csv('data/processed/enriched_acropora_data.csv')
     # Implement ML model
     model = RandomForestClassifier()
     model.fit(data[['feature1', 'feature2']], data['target'])
     ```

5. **Technical Report**:
   - Create `technical_report.Rmd` in the `reports/` directory to document the workflow, including sections on data acquisition, analysis, and results.

6. **Scientific Uses in Spanish**:
   - Write `scientific_uses_es.md` to explain how to obtain marine data, plot and analyze it with R, and utilize Python for machine learning applications. Include examples and explanations in Spanish.

7. **GitHub Actions Workflow**:
   - Create a CI/CD workflow in `.github/workflows/ci.yml` to automate testing and deployment:
     ```yaml
     name: CI

     on:
       push:
         branches:
           - main

     jobs:
       build:
         runs-on: ubuntu-latest
         steps:
           - name: Checkout code
             uses: actions/checkout@v2

           - name: Set up Nix
             run: |
               curl -L https://nixos.org/nix/install | sh
               . /home/runner/.nix-profile/etc/profile.d/nix.sh

           - name: Build Nix environment
             run: nix-build

           - name: Run targets pipeline
             run: nix-shell --run "Rscript -e 'targets::tar_make()'"
     ```

8.  **Documentation**:
    - Write a `README.md` file to provide an overview of the project, installation instructions, and how to run the pipeline.

### Conclusion

This project structure and setup will allow you to create a reproducible workflow for studying marine biodiversity using R and Python. The use of `r-rix` for environment management, `targets` for pipeline management, and GitHub Actions for CI/CD will ensure that your analyses are reproducible and well-documented.

## License

The code in this repository, including all R scripts, the `targets` pipeline, and Quarto setup, is licensed under the **MIT License**. See the `LICENSE` file for full details.

The generated academic report (e.g., `report.html`) is licensed under a **Creative Commons Attribution 4.0 International (CC BY 4.0)** license.

The input marine occurrence data is sourced from [Name of Provider, e.g., GBIF] and is subject to its own licensing terms, which can be found [here (link to data license)].

 It's one of the most popular licenses in the open-source world, especially for software and code-based projects.

Why the MIT License is a Great Fit for You:

    Maximizes Reuse and Collaboration: It allows anyone—other researchers, students, institutions, or even companies—to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of your code without restriction. This aligns perfectly with the open-science ethos.

    Simple and Understandable: The license is short, clear, and doesn't require a legal expert to understand. This lowers the barrier for others to use your work.

    Compatible with Your Dependencies: The R packages you mentioned (rix, targets, etc.) are almost certainly licensed under permissive licenses like MIT or similar (Apache 2.0, BSD). Using MIT ensures perfect compatibility.

    Ideal for Academic Building Blocks: Your project is a toolkit for analysis. The MIT License allows others to take your code, adapt it for their own marine data (or even other domains), and build upon it, which can lead to citations and collaborations.

What the MIT License Requires:

The only requirement is that you include the original copyright notice and the MIT License text in any substantial portion of the software that is redistributed. In practice, this means anyone who uses your code must keep your LICENSE file intact.

Important Clarifications and Other Components

Your project has multiple components (code, data, report). The license above applies to your code. Let's break it down:
1. Your Code (R scripts, targets pipeline, Quarto files)

This is what the MIT License covers.
2. The Input Data (Public Marine Occurrence Data)

You generally do not license data you do not own.
Since the data is publicly available from other sources, you must respect its original license. Your responsibility is to clearly attribute the data sources in your documentation (e.g., in the README).

    Example: If you use data from the Global Biodiversity Information Facility (GBIF), their data is typically licensed under CC0 1.0 (Public Domain Dedication) or CC BY 4.0, which requires attribution. You would state in your README: "The occurrence data in this repository was sourced from GBIF and is made available under the CC BY 4.0 license."

3. The Output Report (Generated Quarto HTML/PDF)

The academic report you generate is a separate creative work. For this, a Creative Commons license is more appropriate.

    Recommended: CC BY 4.0
    This license allows others to distribute, remix, adapt, and build upon your report, even commercially, as long as they credit you for the original creation. This is the standard for open-access academic work and ensures you get credit.