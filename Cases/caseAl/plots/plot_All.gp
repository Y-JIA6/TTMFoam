# =================================================================
# Script Gnuplot : Validation TTM - Multi-Fluences
# =================================================================

set terminal pngcairo size 1000,700 enhanced font 'Arial,14'
set output 'TTM_All.png'

# Titres et labels
set title "Evolution of Te and Tl for different fluences (280 fs)"
set xlabel "Time (ps)"
set ylabel "Temperature (Kelvin)"
set y2label "Normalized Laser Intensity"

# Configuration des axes
set xrange [-1:21]
set yrange [-10000:200000]
set ytics 25000 nomirror

# Configuration de l'Axe Y2 (Laser)
set y2range [-0.06:1.1]       # Laisse un peu d'espace au-dessus du pic laser
set y2tics 0.2 nomirror

set grid lc rgb "#E0E0E0"

# Legende : decommentez si elle chevauche les courbes
# set key outside right

# =================================================================
# PALETTE : une couleur par fluence, en hexadecimal.
#
# NB : les noms internes de gnuplot s'ecrivent avec un tiret
# ("dark-red", "dark-blue", "dark-orange"...) ; "darkred" ou
# "darkeorange" => "unrecognized color name". Et meme les noms
# corrects ne donnent pas les teintes attendues (le "purple" de
# gnuplot est un mauve clair #C080FF, son "dark-orange" une brique
# #C04000). Les codes hexa donnent exactement la couleur voulue.
#
# Numerotation : 1er chiffre = fluence,
#                dernier chiffre = 1 (Te, pointille) / 2 (Tl, plein)
# =================================================================
set style line 11  lc rgb "#00008B" lw 3 dashtype (3,3)   # Te 1 J/cm2 - bleu fonce
set style line 12  lc rgb "#00008B" lw 4                  # Tl 1 J/cm2
set style line 31  lc rgb "#800080" lw 3 dashtype (3,3)   # Te 3 J/cm2 - violet
set style line 32  lc rgb "#800080" lw 4                  # Tl 3 J/cm2
set style line 51  lc rgb "#006400" lw 3 dashtype (3,3)   # Te 5 J/cm2 - vert fonce
set style line 52  lc rgb "#006400" lw 4                  # Tl 5 J/cm2
set style line 71  lc rgb "#FF8C00" lw 3 dashtype (3,3)   # Te 7 J/cm2 - orange fonce
set style line 72  lc rgb "#FF8C00" lw 4                  # Tl 7 J/cm2

# References
set style line 100 lc rgb "black"   lw 3 dashtype (8,2,1,2)   # Seuil d'ablation
set style line 101 lc rgb "#8B0000" lw 3 dashtype (8,2)       # Laser - rouge fonce

# Definition de l'impulsion laser
tp = 0.28
t0 = 2 * tp
Laser(x) = exp(-2.77 * ((x - t0)/tp)**2)

# Seuil de fusion/ablation de l'aluminium
# (attention : l'ancien commentaire disait 6700 K, la valeur est 6030)
T_fusion = 6030

# ---------------------------------------------------------------------------------
# TRACE DES COURBES MULTIPLES
# ---------------------------------------------------------------------------------
plot \
    T_fusion with lines ls 100 title "Ablation Threshold", \
    Laser(x) axes x1y2 with lines ls 101 title "Laser", \
    'resultats_fluences/Te_1.0' using ($1*1e12):2 with lines ls 11 title "Te 1 J/cm²", \
    'resultats_fluences/Tl_1.0' using ($1*1e12):2 with lines ls 12 title "Tl 1 J/cm²", \
    'resultats_fluences/Te_3.0' using ($1*1e12):2 with lines ls 31 title "Te 3 J/cm²", \
    'resultats_fluences/Tl_3.0' using ($1*1e12):2 with lines ls 32 title "Tl 3 J/cm²", \
    'resultats_fluences/Te_5.0' using ($1*1e12):2 with lines ls 51 title "Te 5 J/cm²", \
    'resultats_fluences/Tl_5.0' using ($1*1e12):2 with lines ls 52 title "Tl 5 J/cm²", \
    'resultats_fluences/Te_9.0' using ($1*1e12):2 with lines ls 71 title "Te 9 J/cm²", \
    'resultats_fluences/Tl_9.0' using ($1*1e12):2 with lines ls 72 title "Tl 9 J/cm²"

# Pour ajouter 9 J/cm2 : definir ls 91 / ls 92 avec une nouvelle
# couleur (ex. "#B22222") et ajouter deux lignes Te_9.0 / Tl_9.0
# sur le meme modele.