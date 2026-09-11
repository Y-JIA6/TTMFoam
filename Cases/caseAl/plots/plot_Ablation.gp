# =========================================================================
# Script Gnuplot: Aluminum ablation validation against Omeñaca et al.
# Compares depth and squared diameter with experimental uncertainty bars.
# =========================================================================

set terminal pngcairo size 1200,600 enhanced font 'Arial,14'
set output 'Ablation.png'

# Use a two-panel layout.
set multiplot layout 1,2 title "Ablation Depth and Squared Diameter" font ",14"
set bars 1.5

# Experimental columns: fluence, depth, depth_low, depth_high, diameter2, diameter2_low, diameter2_high.
$ExpData << EOD
1.0  65  40  90  130 122 138
2.0 140 120 158  250 240 260
3.0 180 162 200  350 337 367
4.0 200 170 230  390 360 418
5.0 225 202 245  425 407 445
6.0 250 227 276  490 470 510
7.0 280 260 302  500 480 520
8.0 310 273 350  510 495 525
9.0 340 295 377  530 510 550
EOD

# Omeñaca simulation columns: fluence, depth, squared diameter.
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

# Plot the ablation depth.
set title "a) Ablation depths"
set xlabel "Fluence (J/cm^2)"
set ylabel "Depth (nm)"
set xrange [0:11]
set yrange [-10:400]
set grid lc rgb "#E0E0E0"
set key top left

set arrow 1 from 5.0, graph 0 to 5.0, graph 1 nohead dt 4 lw 1.5 lc rgb "black"
set label 1 "Regime 1" at graph 0.14,0.06 font ",12"
set label 2 "Regime 2" at graph 0.68,0.06 font ",12"

plot \
    $ExpData using 1:2:3:4 with yerrorbars pt 5 ps 1.1 lw 1.5 lc rgb "#2045D8" title "Experimental", \
    $SimData using 1:2 with points pt 7 ps 0.8 lc rgb "black" title "Omeñaca simulation", \
    'Res.txt' using 1:3 with points pt 11 ps 2.0 lc rgb "red" title "OpenFOAM"

# Plot the squared ablation diameter.
set title "b) Squared Diameters"
set xlabel "Fluence (J/cm^2)"
set ylabel "Diameter^2 ({/Symbol m}m^2)"
set xrange [0:11]
set yrange [-20:650]

plot \
    $ExpData using 1:5:6:7 with yerrorbars pt 5 ps 1.1 lw 1.5 lc rgb "#2045D8" title "Experimental", \
    $SimData using 1:3 with points pt 7 ps 0.8 lc rgb "black" title "Omeñaca simulation", \
    'Res.txt' using 1:($5**2) with points pt 11 ps 2.0 lc rgb "red" title "OpenFOAM"

unset multiplot
