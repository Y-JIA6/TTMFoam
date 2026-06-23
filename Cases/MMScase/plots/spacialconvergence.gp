#==============================================================================
# convergence.gp - MMS convergence study for the TTM solver
#
# INPUT FILE: mms_conv.dat  (you create this by hand after the runs)
#   Format - one line per mesh, whitespace separated:
#     Nz    L2_Te        L2_Tl
#     20    1.00e+0      8.0e-1
#     40    2.50e-1      2.0e-1
#     80    6.25e-2      5.0e-2
#     160   1.56e-2      1.25e-2
#
# RUN:  gnuplot convergence.gp     -> produces convergence.png
#       (or: gnuplot -p convergence.gp  to keep the window open on screen)
#==============================================================================

set datafile separator whitespace

# grid spacing h is proportional to 1/Nz
h(Nz) = 1.0/Nz

# ---- fit log(error) = p*log(h) + c  ; slope p = observed order -------------
f_te(x) = p_te*x + c_te
f_tl(x) = p_tl*x + c_tl

fit f_te(x) 'mms_conv.dat' using (log(h($1))):(log($2)) via p_te, c_te
fit f_tl(x) 'mms_conv.dat' using (log(h($1))):(log($3)) via p_tl, c_tl

# ---- print observed orders to the terminal ---------------------------------
print "============================================"
print sprintf(" Observed order  p(Te) = %.3f", p_te)
print sprintf(" Observed order  p(Tl) = %.3f", p_tl)
print "============================================"

# ---- log-log plot ----------------------------------------------------------
set terminal pngcairo size 900,700 enhanced font 'Helvetica,12'
set output 'convergence.png'

set logscale xy
set grid
set key top left
set xlabel "grid spacing  h ~ 1/Nz"
set ylabel "L2 error"
set title sprintf("MMS convergence:  p(Te)=%.2f   p(Tl)=%.2f", p_te, p_tl)

# reference 2nd-order triangle (anchored to the finest Te point)
# slope-2 guide line: y = C * h^2
stats 'mms_conv.dat' using (h($1)):2 nooutput
Cref = STATS_max_y / (STATS_max_x**2)
g2(x) = Cref * x**2

plot 'mms_conv.dat' using (h($1)):2 with linespoints pt 7 ps 1.5 lc rgb 'red'   title 'L2(Te)', \
     'mms_conv.dat' using (h($1)):3 with linespoints pt 9 ps 1.5 lc rgb 'blue'  title 'L2(Tl)', \
     g2(x) with lines dt 2 lc rgb 'black' title 'slope 2 (ref)'

print "Wrote convergence.png"
