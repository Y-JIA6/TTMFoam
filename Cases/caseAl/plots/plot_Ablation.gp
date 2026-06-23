# =================================================================
# Script Gnuplot : Comparaison des résultats d'ablation (Fig 4)
# =================================================================

set terminal pngcairo size 1200,600 enhanced font 'Arial,14'
set output 'Ablation.png'

# Utilisation d'un affichage en 2 colonnes (Multiplot)
set multiplot layout 1,2 title "Diameter and Depth Comparison" font ",14"

# =================================================================
# 1. DÉFINITION DES FONCTIONS THÉORIQUES (Tableau 2)
# =================================================================

# Gnuplot utilise log() pour le logarithme népérien (ln).
# L'expression (x > seuil) ? valeur : 1/0 permet de ne pas tracer 
# la courbe quand la fluence est inférieure au seuil d'ablation (1/0 = erreur ignorée).

# a) Profondeur (Eq. 11)
delta_th = 95.66
F_th_depth = 0.49
Depth(x) = (x > F_th_depth) ? delta_th * log(x / F_th_depth) : 1/0

# b) Diamètre au carré (Eq. 12)
w0_th = 9.97
F_th_diam = 0.53
Diam2(x) = (x > F_th_diam) ? 2 * (w0_th**2) * log(x / F_th_diam) : 1/0


# =================================================================
# 2. BLOCS DE DONNÉES INTÉGRÉS
# Colonnes : 1=Fluence | 2=Profondeur(nm) | 3=Diamètre_Carré(µm2)
# =================================================================

$ExpData << EOD
1.0 65 130
2.0 140 250
3.0 180 350
4.0 200 390
5.0 225 425
6.0 250 490
7.0 280 500
8.0 310 510
9.0 340 530
EOD

$SimData << EOD
0.5 0 0
1.2 90 150
2.0 125 240
2.6 150 310
3.2 175 360
3.8 195 410
4.4 215 440
5.2 230 470
5.8 250 490
6.8 265 510
7.2 280 530
8.0 295 550
8.5 305 560
9.2 320 580
10.0 335 590
10.8 345 600
EOD

# =================================================================
# 3. TRACÉ DU GRAPHIQUE GAUCHE : PROFONDEUR
# =================================================================
set title "a) Ablation depths"
set xlabel "Fluence (J/cm^2)"
set ylabel "Depth (nm)"
set xrange [0:11]
set yrange [-10:400]
set grid lc rgb "#E0E0E0"
set key top left

# Ligne de séparation des régimes
set arrow from 5.0, graph 0 to 5.0, graph 1 nohead dt 4 lc rgb "black"

plot \
    $ExpData using 1:2 with points pt 7 ps 1.2 lc rgb "blue" title "Experimental", \
    $SimData using 1:2 with points pt 5 ps 0.8 lc rgb "black" title "Omeñaca simulation", \
    'Res.txt' using 1:3 with points pt 11 ps 2.0 lc rgb "red" title "OpenFOAM"

# =================================================================
# 4. TRACÉ DU GRAPHIQUE DROIT : DIAMÈTRE AU CARRÉ
# =================================================================
set title "b) Squared Diameters"
set xlabel "Fluence (J/cm^2)"
set ylabel "Diameter^2 ({/Symbol m}m^2)"
set xrange [0:11]
set yrange [-20:650]

set arrow from 5.0, graph 0 to 5.0, graph 1 nohead dt 4 lc rgb "black"

# Attention ici : pour MyData, on utilise ($3**2) pour élever 
# votre diamètre mesuré au carré dynamiquement !
plot \
    $ExpData using 1:3 with points pt 7 ps 1.2 lc rgb "blue" title "Experimental", \
    $SimData using 1:3 with points pt 5 ps 0.8 lc rgb "black" title "Omeñaca simulation", \
    'Res.txt' using 1:($5**2) with points pt 11 ps 2.0 lc rgb "red" title "OpenFOAM"

unset multiplot