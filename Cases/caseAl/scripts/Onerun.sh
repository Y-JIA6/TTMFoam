#!/bin/bash

# Arrêter le script en cas d'erreur fatale
set -e

# =================================================================
# Lancement du chronomètre global
# =================================================================
START_TOTAL=$(date +%s)

echo "================================================="
echo "   Compilation et Lancement TTM - Ablation       "
echo "================================================="

# 1. NETTOYAGE (Sécurisé et Manuel)
echo "[1/7] Préparation du dossier"
#rm -f Res.txt
rm -f log.*
foamListTimes -rm
rm -rf processor* 2>/dev/null || echo "Dossiers processeurs verrouillés, on continue..."
rm -rf postProcessing/ 2>/dev/null || true

# 2. MAILLAGE
echo "[2/7] Création du maillage (blockMesh)"
blockMesh > log.blockMesh 2>&1

# 3. DÉCOMPOSITION
echo "[3/7] Décomposition du domaine"
decomposePar -force > log.decomposePar 2>&1

# 4. CALCUL EN PARALLÈLE
echo "[4/7] Résolution avec TTm en parallèle"
mpirun -np 4 TTmAl -parallel | tee log.TTmAl | grep --line-buffered -E "Time =|PGC|ExecutionTime"

# 5. RECONSTRUCTION
echo "[5/7] Reconstruction des champs (reconstructPar)"
START_RECON=$(date +%s)
reconstructPar > log.reconstructPar 2>&1 || true
END_RECON=$(date +%s)
TIME_RECON=$((END_RECON - START_RECON))
echo "      Temps de reconstruction : ${TIME_RECON} s"

# 6. POST-PROCESSING (Extraction des sondes)
echo "[6/7] Extraction des données de sondes (postProcess)"
START_POST=$(date +%s)
postProcess -func probes > log.postProcess
END_POST=$(date +%s)
TIME_POST=$((END_POST - START_POST))
echo "      Temps d'extraction : ${TIME_POST} s"

echo "[7/7] Traçage des courbes de pics de température"
gnuplot plots/plot_Temperature.gp 

# =================================================================
# Arrêt du chronomètre global
# =================================================================
END_TOTAL=$(date +%s)
TIME_TOTAL=$((END_TOTAL - START_TOTAL))

echo "================================================="
echo "   Calcul terminé ! Données prêtes dans probes/  "
echo "   TEMPS TOTAL D'EXÉCUTION : ${TIME_TOTAL} s   "
echo "================================================="