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
echo "Paso 2: Probando disponibilidad de paquetes..."
echo ""

# Simplemente ejecutar el script R que ya funciona
nix-shell --run "Rscript R/check_packages.R"

echo ""
echo "================================================"
echo "Siguientes pasos para comenzar a trabajar:"
echo "  1. nix-shell              # Ejecutar en la consola para acceder al ambiente"
echo "  2. R                      # Escribir R en la consola"  
echo "  3. targets::tar_make()    # Ejecutar la pipeline en la consola de R"
echo ""