# convergence_time.gp - temporal order from mms_time.dat (deltaT L2Te L2Tl)
set datafile separator whitespace
f_te(x)=p_te*x+c_te
f_tl(x)=p_tl*x+c_tl
fit f_te(x) 'mms_time.dat' using (log($1)):(log($2)) via p_te,c_te
fit f_tl(x) 'mms_time.dat' using (log($1)):(log($3)) via p_tl,c_tl
print "============================================"
print sprintf(" Temporal order  p(Te) = %.3f", p_te)
print sprintf(" Temporal order  p(Tl) = %.3f", p_tl)
print "============================================"
set terminal pngcairo size 900,700 enhanced font 'Helvetica,12'
set output 'convergence_time.png'
set logscale xy
set grid
set key top left
set xlabel "time step  deltaT [s]"
set ylabel "L2 error"
set title sprintf("MMS temporal convergence:  p(Te)=%.2f  p(Tl)=%.2f", p_te, p_tl)
stats 'mms_time.dat' using 1:2 nooutput
C1 = STATS_max_y/STATS_max_x
g1(x)=C1*x
plot 'mms_time.dat' u 1:2 w lp pt 7 ps 1.5 lc rgb 'red'  t 'L2(Te)', \
     'mms_time.dat' u 1:3 w lp pt 9 ps 1.5 lc rgb 'blue' t 'L2(Tl)', \
     g1(x) w l dt 2 lc rgb 'black' t 'slope 1 (ref)'
