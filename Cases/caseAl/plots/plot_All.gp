# =========================================================================
# Script Gnuplot: TTM transient temperatures at representative fluences
# Compares OpenFOAM electron and lattice temperatures over the first 20 ps.
# =========================================================================

set terminal pngcairo size 1400,800 enhanced font 'Arial,14'
set output 'TTM_All.png'
set title "OpenFOAM temperature histories (280 fs pulse): solid T_e, dashed T_l, markers peak T_e"
set xlabel "Time (ps)"
set ylabel "Temperature (K)"
set xrange [-1:21]
set yrange [-10000:250000]
set xtics 5
set ytics 25000
set grid lc rgb "#E0E0E0"
set key outside right vertical font ',11'

# One colour is assigned to each fluence; Te is solid and Tl is dashed.
set style line 1 lc rgb "#1f77b4" lw 2.5
set style line 3 lc rgb "#2ca02c" lw 2.5
set style line 5 lc rgb "#ff7f0e" lw 2.5
set style line 9 lc rgb "#d62728" lw 2.5

stats 'resultats_fluences/Te_1.0' using 2 name 'Te1' nooutput
stats 'resultats_fluences/Te_3.0' using 2 name 'Te3' nooutput
stats 'resultats_fluences/Te_5.0' using 2 name 'Te5' nooutput
stats 'resultats_fluences/Te_9.0' using 2 name 'Te9' nooutput

plot \
    'resultats_fluences/Te_1.0' using (($1*1e12 <= 20) ? ($1*1e12) : 1/0):2 with lines ls 1 title 'OpenFOAM, 1 J cm^{-2}', \
    'resultats_fluences/Tl_1.0' using (($1*1e12 <= 20) ? ($1*1e12) : 1/0):2 with lines ls 1 dashtype (6,3) notitle, \
    'resultats_fluences/Te_1.0' using (($2 == Te1_max) ? ($1*1e12) : 1/0):2 with points pt 7 ps 1.1 lc rgb '#1f77b4' notitle, \
    'resultats_fluences/Te_3.0' using (($1*1e12 <= 20) ? ($1*1e12) : 1/0):2 with lines ls 3 title 'OpenFOAM, 3 J cm^{-2}', \
    'resultats_fluences/Tl_3.0' using (($1*1e12 <= 20) ? ($1*1e12) : 1/0):2 with lines ls 3 dashtype (6,3) notitle, \
    'resultats_fluences/Te_3.0' using (($2 == Te3_max) ? ($1*1e12) : 1/0):2 with points pt 7 ps 1.1 lc rgb '#2ca02c' notitle, \
    'resultats_fluences/Te_5.0' using (($1*1e12 <= 20) ? ($1*1e12) : 1/0):2 with lines ls 5 title 'OpenFOAM, 5 J cm^{-2}', \
    'resultats_fluences/Tl_5.0' using (($1*1e12 <= 20) ? ($1*1e12) : 1/0):2 with lines ls 5 dashtype (6,3) notitle, \
    'resultats_fluences/Te_5.0' using (($2 == Te5_max) ? ($1*1e12) : 1/0):2 with points pt 7 ps 1.1 lc rgb '#ff7f0e' notitle, \
    'resultats_fluences/Te_9.0' using (($1*1e12 <= 20) ? ($1*1e12) : 1/0):2 with lines ls 9 title 'OpenFOAM, 9 J cm^{-2}', \
    'resultats_fluences/Tl_9.0' using (($1*1e12 <= 20) ? ($1*1e12) : 1/0):2 with lines ls 9 dashtype (6,3) notitle, \
    'resultats_fluences/Te_9.0' using (($2 == Te9_max) ? ($1*1e12) : 1/0):2 with points pt 9 ps 1.1 lc rgb '#d62728' notitle
