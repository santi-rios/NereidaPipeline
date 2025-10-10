# check_packages.R

# Leer el archivo _targets.R
targets_file <- readLines("_targets.R")

# Extraer los paquetes de las llamadas a library(), require(), etc.
# Esta expresión regular busca:
# - library(pkg), require(pkg), requireNamespace(pkg)
# - tar_quarto(name, path = "...") -> asume 'quarto'
library_lines <- grep("library\\(|require\\(|requireNamespace\\(|tar_quarto\\(", targets_file, value = TRUE)
packages_from_libs <- gsub(".*(?:library|require|requireNamespace)\\s*\\(([^,)]+)\\).*", "\\1", library_lines)
packages_from_libs <- gsub("[\"']", "", packages_from_libs) # Limpiar comillas

# Detectar 'quarto' de tar_quarto
packages_from_quarto <- if (any(grepl("tar_quarto\\(", library_lines))) "quarto" else character(0)

# Combinar todos los paquetes encontrados
packages <- unique(c(packages_from_libs, packages_from_quarto))
packages <- trimws(packages) # Limpiar espacios en blanco
packages <- packages[packages != "" & !grepl("=", packages)] # Eliminar cadenas con '=' que no son paquetes

# Imprimir paquetes encontrados
cat("\nPaquetes encontrados en _targets.R:\n")
cat(paste("  -", packages, collapse = "\n"), "\n\n")

# Probar si los paquetes están disponibles
cat("Probando disponibilidad de paquetes:\n")
success_count <- 0
failed_packages <- c()

for (pkg in packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  ✓ %s\n", pkg))
    success_count <- success_count + 1
  } else {
    cat(sprintf("  ✗ %s - FAILED\n", pkg))
    failed_packages <- c(failed_packages, pkg)
  }
}

# Imprimir resumen
cat(sprintf("\n%d/%d paquetes disponibles (%.1f%%)\n", 
            success_count, length(packages), 
            (success_count / length(packages)) * 100))

# Si faltan paquetes, mostrar advertencia y salir con error
if (length(failed_packages) > 0) {
  cat("\n⚠️  Librerías no disponibles:\n")
  cat(paste("  -", failed_packages, collapse = "\n"), "\n")
  cat("\nAgrega estos paquetes a 'build_env.R' y ejecuta './regenerate.sh'\n")
  quit(status = 1) # Salir con código de error
}

cat("\n✨ Prueba del ambiente finalizada con éxito.\n")