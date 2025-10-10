#!/bin/bash
# Enhanced test script that syncs with _targets.R

echo "🧪 Probando y verificando el ambiente (software, librerías, dependencias, etc.)"
echo "================================================"
echo ""

echo "Paso 1: Verificando que default.nix exista..."
if [ -f "default.nix" ]; then
    echo "✅ default.nix encontrado"
else
    echo "⚠️  default.nix no disponible. Ejecutando regenerate.sh..."
    ./regenerate.sh
fi

echo ""
echo "Paso 2: Obteniendo paquetes utilizados en la pipeline de _targets.R..."
echo ""

# Extract library() calls from _targets.R and create a test script
nix-shell --run '
R --vanilla --quiet -e "
# Extract packages from _targets.R
targets_file <- readLines(\"_targets.R\")

# Find library() calls
library_lines <- grep(\"library\\(\", targets_file, value = TRUE)
packages <- gsub(\".*library\\(([^)]+)\\).*\", \"\\1\", library_lines)
packages <- unique(packages)
packages <- packages[packages != \"\"]

cat(\"\\nPackages found in _targets.R:\\n\")
cat(paste(\"  -\", packages, collapse = \"\\n\"), \"\\n\\n\")

cat(\"Testing package availability:\\n\")
success_count <- 0
failed_packages <- c()

for (pkg in packages) {
  pkg_clean <- gsub(\"[\\\"\\']\", \"\", pkg)  # Remove quotes
  if (requireNamespace(pkg_clean, quietly = TRUE)) {
    cat(sprintf(\"  ✓ %s\\n\", pkg_clean))
    success_count <- success_count + 1
  } else {
    cat(sprintf(\"  ✗ %s - FAILED\\n\", pkg_clean))
    failed_packages <- c(failed_packages, pkg_clean)
  }
}

cat(sprintf(\"\\n%d/%d packages available (%.1f%%)\\n\", 
            success_count, length(packages), 
            (success_count/length(packages))*100))

if (length(failed_packages) > 0) {
  cat(\"\\n⚠️ Librerías no disponibles Missing packages:\\n\")
  cat(paste(\"  -\", failed_packages, collapse = \"\\n\"), \"\\n\")
  cat(\"\\nAgregar estos paquetes al archivo build_env.R y ejecutar en la consola ./regenerate.sh\\n\")
}

cat(\"\\n✨ Prueba del ambiente finalizada\\n\")
"
'

echo ""
echo "================================================"
echo "Siguientes pasos para comenzar a trabajar:"
echo "  1. nix-shell              # Ejecutar en la consola para acceder al ambiente"
echo "  2. R                      # Escribir R en la consola"  
echo "  3. targets::tar_make()    # Ejecutar la pipeline en la consola de R"
echo ""