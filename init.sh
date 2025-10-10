#!/bin/bash
# filepath: /home/santi/Projects/biologia-marina-reproducible/init.sh
# Script de inicialización única - Ejecutar solo la primera vez

set -e

echo ""
echo "🌊 Inicialización del Proyecto de Biología Marina"
echo "==================================================="
echo ""

# 1. Arreglar advertencias de Nix
echo "Paso 1: Limpiando configuración de Nix..."
if [ -f "fix_nix_warnings.sh" ]; then
    chmod +x fix_nix_warnings.sh
    ./fix_nix_warnings.sh
else
    # Crear script si no existe
    cat > fix_nix_warnings.sh << 'FIXNIX'
#!/bin/bash
NIX_CONF="$HOME/.config/nix/nix.conf"
mkdir -p "$(dirname "$NIX_CONF")"
if [ -f "$NIX_CONF" ]; then
    cp "$NIX_CONF" "$NIX_CONF.backup"
    sed -i '/eval-cores/d' "$NIX_CONF"
    sed -i '/lazy-trees/d' "$NIX_CONF"
fi
echo "✅ Configuración de Nix limpiada"
FIXNIX
    chmod +x fix_nix_warnings.sh
    ./fix_nix_warnings.sh
fi

# 2. Establecer permisos
echo ""
echo "Paso 2: Configurando permisos de scripts..."
chmod +x *.sh 2>/dev/null || true
echo "✅ Permisos configurados"

# 3. Crear .gitignore si no existe
echo ""
echo "Paso 3: Configurando Git..."
if [ ! -f ".gitignore" ]; then
    cat > .gitignore << 'GITIGNORE'
# R y RStudio
.Rproj.user/
.Rhistory
.RData
.Ruserdata
*.Rproj

# targets
_targets/
!_targets.R

# Datos grandes
*.rds
*.RData
data/raw/*
data/processed/*
!data/raw/.gitkeep
!data/processed/.gitkeep

# Nix
result
result-*
*.backup

# Sistema
.DS_Store
*.swp
*~
GITIGNORE
    echo "✅ .gitignore creado"
else
    echo "✅ .gitignore ya existe"
fi

# 4. Crear estructura de directorios
echo ""
echo "Paso 4: Creando estructura de directorios..."
mkdir -p data/raw data/processed R docs
touch data/raw/.gitkeep data/processed/.gitkeep
echo "✅ Estructura creada"

# 5. Generar ambiente Nix
echo ""
echo "Paso 5: Generando ambiente Nix inicial..."
if [ ! -f "default.nix" ]; then
    ./regenerate.sh
else
    echo "✅ default.nix ya existe"
fi

echo ""
echo "✅ ¡Inicialización Completada!"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "📖 Guía Rápida de Uso"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Para trabajar en el proyecto:"
echo ""
echo "  1. Entrar al ambiente:"
echo "     $ nix-shell"
echo ""
echo "  2. Ver comandos disponibles:"
echo "     $ make help"
echo ""
echo "  3. Iniciar R y ejecutar pipeline:"
echo "     $ R"
echo "     > targets::tar_make()"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "🔧 Comandos Útiles con Make"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "  make help        → Ver todos los comandos"
echo "  make regenerate  → Regenerar ambiente Nix"
echo "  make test        → Probar paquetes instalados"
echo "  make update      → Actualizar todo (regenerar + construir + probar)"
echo "  make clean       → Limpiar archivos temporales"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""