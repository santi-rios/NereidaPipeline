#!/bin/bash
# filepath: /home/santi/Projects/biologia-marina-reproducible/fix_nix_warnings.sh
# Arreglar advertencias de configuración de Nix
# modifica configuración a nivel de usuario ~/.config/nix/nix.conf
# crea un backup por si hay que revertir

echo "🔧 Arreglando Advertencias de Nix"
echo "=================================="
echo ""

NIX_CONF="$HOME/.config/nix/nix.conf"

# Crear directorio si no existe
mkdir -p "$(dirname "$NIX_CONF")"

# Limpiar configuraciones obsoletas
echo "Limpiando configuraciones obsoletas..."

if [ -f "$NIX_CONF" ]; then
    # Respaldar configuración actual
    cp "$NIX_CONF" "$NIX_CONF.backup"
    echo "✅ Respaldo creado: $NIX_CONF.backup"
    
    # Eliminar configuraciones obsoletas
    sed -i '/eval-cores/d' "$NIX_CONF"
    sed -i '/lazy-trees/d' "$NIX_CONF"
    
    echo "✅ Configuraciones obsoletas eliminadas"
else
    touch "$NIX_CONF"
    echo "✅ Archivo de configuración creado"
fi

echo ""
echo "Configuración actual de Nix:"
echo "----------------------------"
cat "$NIX_CONF"
echo ""

echo "✅ Advertencias corregidas!"
echo ""
echo "Las advertencias sobre 'untrusted substituter' son normales"
echo "y no afectan el funcionamiento del proyecto."
echo ""