#!/bin/bash
# filepath: ~/Projects/biologia-marina-reproducible/fix_git_large_files.sh
# Script para limpiar archivos grandes del historial de Git
# chmod +x fix_git_large_files.sh

set -e

echo "🔧 Limpieza de Archivos Grandes en Git"
echo "======================================="
echo ""
echo "⚠️  ADVERTENCIA: Este script modificará el historial de Git"
echo "Asegúrate de tener un respaldo antes de continuar"
echo ""

read -p "¿Continuar con la limpieza? (s/n): " -n 1 -r
echo ""
if [[ ! $REPLY =~ ^[Ss]$ ]]; then
    echo "Operación cancelada"
    exit 0
fi

# Verificar que git-filter-repo esté disponible
if ! command -v git-filter-repo &> /dev/null; then
    echo "❌ git-filter-repo no está instalado"
    echo ""
    echo "Para instalar:"
    echo "  1. Agrega 'git-filter-repo' a system_pkgs en build_env.R"
    echo "  2. Ejecuta ./update_workflow.sh"
    echo "  3. Entra al ambiente con: nix-shell"
    echo "  4. Ejecuta este script nuevamente"
    exit 1
fi

echo ""
echo "Paso 1: Encontrando archivos grandes en el historial"
echo "-----------------------------------------------------"

# Encontrar los 20 archivos más grandes en el historial
git rev-list --objects --all |
  git cat-file --batch-check='%(objecttype) %(objectname) %(objectsize) %(rest)' |
  sed -n 's/^blob //p' |
  sort --numeric-sort --key=2 --reverse |
  head -n 20 |
  awk '{print $2/1024/1024 " MB\t" $3}' > /tmp/large_files.txt

echo "Archivos grandes encontrados:"
cat /tmp/large_files.txt
echo ""

# Paso 2: Crear lista de archivos a eliminar
echo "Paso 2: Seleccionar archivos para eliminar"
echo "-------------------------------------------"
echo "Escribe los patrones de archivos a eliminar (uno por línea)"
echo "Ejemplos:"
echo "  *.rds"
echo "  data/large_dataset.csv"
echo "  _targets/objects/*"
echo ""
echo "Presiona Ctrl+D cuando termines"
echo ""

# Leer archivos a eliminar
PATHS_FILE="/tmp/paths_to_remove.txt"
cat > "$PATHS_FILE"

if [ ! -s "$PATHS_FILE" ]; then
    echo "No se especificaron archivos. Saliendo..."
    exit 0
fi

echo ""
echo "Archivos/patrones a eliminar:"
cat "$PATHS_FILE"
echo ""

read -p "¿Continuar con la eliminación? (s/n): " -n 1 -r
echo ""
if [[ ! $REPLY =~ ^[Ss]$ ]]; then
    echo "Operación cancelada"
    exit 0
fi

# Paso 3: Hacer respaldo
echo ""
echo "Paso 3: Creando respaldo"
echo "------------------------"
BACKUP_DIR="../biologia-marina-backup-$(date +%Y%m%d-%H%M%S)"
echo "Creando respaldo en: $BACKUP_DIR"
cp -r . "$BACKUP_DIR"
echo "✅ Respaldo creado"

# Paso 4: Limpiar archivos del historial
echo ""
echo "Paso 4: Limpiando historial de Git"
echo "-----------------------------------"
echo "Esto puede tomar varios minutos..."

# Leer cada patrón y ejecutar git-filter-repo
while IFS= read -r path; do
    if [ ! -z "$path" ]; then
        echo "Eliminando: $path"
        git filter-repo --path "$path" --invert-paths --force
    fi
done < "$PATHS_FILE"

echo "✅ Historial limpiado"

# Paso 5: Limpiar y optimizar repositorio
echo ""
echo "Paso 5: Optimizando repositorio"
echo "--------------------------------"
git reflog expire --expire=now --all
git gc --prune=now --aggressive
echo "✅ Repositorio optimizado"

# Paso 6: Mostrar tamaño del repositorio
echo ""
echo "Paso 6: Verificando tamaño"
echo "--------------------------"
du -sh .git
echo ""

echo "✅ ¡Limpieza completada!"
echo ""
echo "Próximos pasos:"
echo "  1. Verifica que el repositorio funcione correctamente"
echo "  2. Agrega el remote nuevamente:"
echo "     git remote add origin https://github.com/santi-rios/biologia-marina-reproducible.git"
echo "  3. Haz push forzado (⚠️ cuidado):"
echo "     git push --force --all"
echo "     git push --force --tags"
echo ""
echo "⚠️  IMPORTANTE: Informa a tu equipo sobre el push forzado"
echo ""