# =========================================================================
# Script Gnuplot: TTM transient-temperature validation for aluminum
# Generates the manuscript figure of Te, Tl, the laser pulse, and thresholds.
# =========================================================================

# Select the probe directory containing the optical diagnostics from the solver run.
P = system("for d in postProcessing/probes/*/; do [ -f \"$d/R_field\" ] && { printf %s \"$d\"; break; }; done")
if (strlen(P) == 0) { print "ERREUR : sondes introuvables (bloc functions{} actif ?)" ; exit }
print "Lecture des sondes dans : ".P

COL = 2      # Column 2 is probe 0 and column 3 is probe 1.

# Use the same vertical scale as plot_All.gp.
YMAX = 250000
YTIC = 25000      # Vertical tick interval in K.

# Keep these physical parameters synchronized with the solver.
tp      = 0.28        # Pulse duration in ps.
t0      = 2 * tp      # Time of the laser maximum in ps.
T_crit  = 6700.0      # Aluminum critical temperature in K.
T_abl   = 0.9*T_crit  # Phase-explosion threshold used by the solver.
T_melt  = 933.47      # Aluminum melting temperature in K.

# Reproduce the normalized temporal Gaussian used by calculateLaser.H.
Laser(x) = exp(-2.7726 * ((x - t0)/tp)**2)

# Measure and annotate the electron-temperature maximum.
stats P.'Te' using ($1*1e12):COL nooutput
TePeak = STATS_max_y
tPeak  = STATS_pos_max_y
Tmax   = (YMAX > 0) ? YMAX : 1.18*TePeak

# Add margins on every side of the linear-axis figure.
Ymin   = -10000
Xmin   = -0.6
Xmax   = 20.6

print sprintf("Pic Te = %.0f K a t = %.2f ps", TePeak, tPeak)

# Generate the linear-axis overview.
set terminal pngcairo size 900,620 enhanced font 'Arial,14'
set output 'TTM.png'

set title "OpenFOAM TTM: T_e, T_l and laser pulse ({/Symbol t}_p = 280 fs)"
set xlabel "Time (ps)"
set ylabel "Temperature (K)"
set y2label "Normalized laser intensity"

set xrange [Xmin:Xmax]
set yrange [Ymin:Tmax]
set xtics 2.5
set ytics YTIC
set grid lc rgb "#E0E0E0"

# Keep margins below and above the normalized laser pulse.
set y2range [-0.05:1.18]
set y2tics 0.2 nomirror
set ytics nomirror

set key top right box opaque height 0.4 width -2

# Annotate the electron-temperature maximum.
set label 1 sprintf("T_e^{max} = %.0f K", TePeak) \
    at tPeak, TePeak offset 1.2,0.6 font ",11" tc rgb "#B22222"
set arrow 1 from tPeak, TePeak to tPeak, TePeak nohead
set arrow 1 from tPeak, Ymin to tPeak, TePeak nohead dt 3 lw 1 lc rgb "#B22222"

plot \
    P.'Te' using ($1*1e12):COL axes x1y1 with lines lw 3 lc rgb "#B22222" \
        title "OpenFOAM T_e", \
    P.'Tl' using ($1*1e12):COL axes x1y1 with lines lw 3 lc rgb "#1F4E9C" \
        title "OpenFOAM T_l", \
    T_abl  axes x1y1 with lines lw 2 dt (8,3) lc rgb "black" \
        title sprintf("Ablation : 0.9 T_c = %.0f K", T_abl), \
    T_melt axes x1y1 with lines lw 1.5 dt (2,3) lc rgb "#666666" \
        title sprintf("Melting : %.0f K", T_melt), \
    Laser(x) axes x1y2 with lines lw 2.5 dt (8,2) lc rgb "#E08214" \
        title "Laser"

unset label 1
unset arrow 1

print "Figure generated: TTM.png"
