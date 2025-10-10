# Makefile para gestión del proyecto de Biología Marina
# IMPORTANTE: Este archivo usa TABS, no espacios

.PHONY: help setup permissions regenerate test update clean git-fix nix-warnings

help: ## Mostrar esta ayuda
	@echo ""
	@echo "📚 Comandos Disponibles para Biología Marina"
	@echo "============================================="
	@echo ""
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
        awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}'
	@echo ""

setup: permissions ## Configuración inicial del proyecto
	@echo "✅ Proyecto configurado y listo para usar"

permissions: ## Establecer permisos de scripts
	@echo "🔐 Configurando permisos de scripts..."
	@chmod +x *.sh 2>/dev/null || true
	@echo "✅ Permisos actualizados"

regenerate: permissions ## Regenerar ambiente Nix desde build_env.R
	@echo "🔧 Regenerando ambiente Nix..."
	@./regenerate.sh

test: ## Probar que todos los paquetes estén disponibles
	@./test_environment.sh

update: permissions ## Flujo completo: regenerar + construir + probar
	@./update_workflow.sh

git-fix: ## Arreglar archivos grandes en el historial de Git
	@./quick_fix_git.sh

nix-warnings: ## Arreglar advertencias de configuración de Nix
	@./fix_nix_warnings.sh

clean: ## Limpiar archivos temporales y objetos de targets
	@echo "🧹 Limpiando archivos temporales..."
	@rm -rf result result-* _targets/objects/* 2>/dev/null || true
	@echo "✅ Limpieza completada"

# deploy-docs: 
#     mkdir -p docs
# 	cp reportes/analisis_biodiversidad_marina.html docs/index.html
# 	cp -r reportes/analisis_biodiversidad_marina_files docs/ 2>/dev/null || true
# 	cp reportes/styles.css docs/ 2>/dev/null || true
# 	cp -r figures docs/ 2>/dev/null || true
# 	git add docs/
# 	git commit -m "Update GitHub Pages documentation"

.DEFAULT_GOAL := help