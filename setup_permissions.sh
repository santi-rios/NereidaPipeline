#!/bin/bash
# filepath: ~/Projects/biologia-marina-reproducible/setup_permissions.sh
# Establecer permisos ejecutables para todos los scripts

echo "🔐 Configurando Permisos de Scripts"
echo "===================================="
echo ""

# Lista de scripts que necesitan ser ejecutables
SCRIPTS=(
    "regenerate.sh"
    "update_workflow.sh"
    "test_environment.sh"
    "quick_fix_git.sh"
    "fix_git_large_files.sh"
    "fix_nix_warnings.sh"
    "setup_permissions.sh"
)

echo "Otorgando permisos de ejecución a:"
for script in "${SCRIPTS[@]}"; do
    if [ -f "$script" ]; then
        chmod +x "$script"
        echo "  ✅ $script"
    else
        echo "  ⚠️  $script (no encontrado)"
    fi
done

# Guardar permisos en Git
echo ""
echo "Actualizando permisos en Git..."
git add --chmod=+x "${SCRIPTS[@]}" 2>/dev/null

echo ""
echo "✅ Permisos configurados!"
echo ""
echo "Los permisos se mantendrán en Git para futuros clones"