#!/bin/bash

# Arret du script en cas d'erreur
set -e

# 1. Creation du dossier resultats
rm -rf resultats_fluences
mkdir -p resultats_fluences
rm -f Res.txt

# 2. Liste des fluences a tester (en J/cm2)
#    A_rel(F) est lue directement depuis Fig A1 dans calculateLaser.H
fluences=(1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0)

for F in "${fluences[@]}"
do
    echo "================================================="
    echo " POUR F = $F J/cm2 "
    echo "================================================="

    # Conversion J/cm2 -> J/m2
    F_m2=$(awk "BEGIN {print $F * 10000}")

    # On modifie SEULEMENT F_laser dans le dictionnaire OpenFOAM
    foamDictionary constant/laserProperties -entry F_laser -set $F_m2

    # Lancement de la simulation
    ./scripts/Onerun.sh

    # Sauvegarde
    cp postProcessing/probes/0/Te resultats_fluences/Te_${F}
    cp postProcessing/probes/0/Tl resultats_fluences/Tl_${F}
done

echo " Tracage des courbes de pics de temperature et d'ablation"
gnuplot plots/plot_All.gp 
gnuplot plots/plot_Ablation.gp

echo "================================================="
echo " FIN! "
echo "================================================="