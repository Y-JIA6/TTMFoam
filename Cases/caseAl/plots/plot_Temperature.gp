# =================================================================
# Script Gnuplot : Validation TTM - Ablation Laser Aluminium
# =================================================================

# Configuration de l'image de sortie
set terminal pngcairo size 800,600 enhanced font 'Arial,14'
set output 'TTM.png'

# Titres et labels
set title "TTM : Te, Tl, and Laser Pulse (280 fs)"
set xlabel "Time (ps)"
set ylabel "Temperature (K)"
set y2label "Normalized Laser Intensity"

# Configuration de l'Axe Y1 (Températures)
set xrange [-1.0:21.0]    # Réduit pour mieux voir le pic (à ajuster selon besoin)
set yrange [-10000:200000]     # Pas de températures négatives
set ytics 25000
set xtics 2.5
set grid lc rgb "#E0E0E0"

# Configuration de l'Axe Y2 (Laser)
set y2range [-0.06:1.1]       # Laisse un peu d'espace au-dessus du pic laser
set y2tics 0.2                # Active les graduations à droite
set y2tics nomirror           
set ytics nomirror            

# Définition de l'impulsion laser
tp = 0.28
t0 = 2 * tp               
Laser(x) = exp(-2.77 * ((x - t0)/tp)**2)

# Seuil de fusion
T_fusion = 6700

# Tracé des courbes
plot \
    'postProcessing/probes/0/Te' using ($1*1e12):2 axes x1y1 with lines lw 3 lc rgb "blue" title "Te ", \
    'postProcessing/probes/0/Tl' using ($1*1e12):2 axes x1y1 with lines lw 3 dashtype (3,3) lc rgb "blue" title "Tl ", \
    T_fusion axes x1y1 with lines lw 3 dt (8, 2, 1, 2) lc rgb "black" title "Ablation", \
    Laser(x) axes x1y2 with lines dashtype (8, 2) lc rgb "red" lw 3 title "Laser"