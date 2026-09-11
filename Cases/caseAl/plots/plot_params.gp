#!/usr/bin/gnuplot
# =========================================================================
# Script Gnuplot: Free-electron-gas property curves for aluminum
# Generates Ce, G, chemical potential, conductivity, and relaxation-time plots.
# =========================================================================

set terminal pdfcairo enhanced font "Helvetica,12" size 6,4

DATAFILE = "FEG_param.dat"

# Keep reference values synchronized with constant/ttmProperties.
gamma_S  = 76.7          # Sommerfeld coefficient in J/m3/K2.
EF_eV    = 11.7          # Fermi energy in eV.
n_FEG    = 2.0*gamma_S*(EF_eV*1.602176634e-19)/(pi**2*(1.380649e-23)**2)
Ce_class = 1.5*n_FEG*1.380649e-23    # Classical limit 3/2 n kB in J/m3/K.
T_F      = EF_eV*1.602176634e-19/1.380649e-23

set xlabel "T_e (K)"
set grid
set key top left
set format x "%g"
set format y "%g"

# Compare OpenFOAM heat capacity with the Sommerfeld and classical limits.
set output "Ce_vs_Te.pdf"
set ylabel "C_e (J m^{-3} K^{-1})"
set key top left
set label 1 sprintf("T_F = %.0f K", T_F) at T_F, graph 0.30 left offset 0.5,0 tc rgb "gray40"
set arrow 1 from T_F, graph 0 to T_F, graph 0.28 nohead dt 3 lc rgb "gray40"
plot DATAFILE using 1:2 with lines lw 2.5 lc rgb "blue" \
         title "OpenFOAM C_e (Fermi-Dirac)", \
     gamma_S*x with lines lw 1.5 dt 2 lc rgb "dark-red" \
         title sprintf("Sommerfeld  {/Symbol g}T_e  ({/Symbol g} = %.1f)", gamma_S), \
     Ce_class with lines lw 1.2 dt 4 lc rgb "gray30" \
         title "classical limit  3/2 n k_B"
unset label 1
unset arrow 1

# Use logarithmic axes to emphasize the departure from the Sommerfeld relation.
set output "Ce_vs_Te_log.pdf"
set logscale xy
set key top left
plot DATAFILE using 1:2 with lines lw 2.5 lc rgb "blue" \
         title "OpenFOAM C_e (Fermi-Dirac)", \
     gamma_S*x with lines lw 1.5 dt 2 lc rgb "dark-red" title "Sommerfeld  {/Symbol g}T_e", \
     Ce_class with lines lw 1.2 dt 4 lc rgb "gray30" title "3/2 n k_B"
unset logscale xy

# Plot the electron-phonon coupling factor.
set output "G_vs_Te.pdf"
set ylabel "G (W m^{-3} K^{-1})"
set key top left
plot DATAFILE using 1:3 with lines lw 2.5 lc rgb "red" title "OpenFOAM G(T_e)"

# Plot the chemical potential obtained from particle-number conservation.
set output "mu_vs_Te.pdf"
set ylabel "{/Symbol m} (eV)"
set key top right
set arrow 2 from graph 0, first 0 to graph 1, first 0 nohead dt 3 lc rgb "gray40"
plot DATAFILE using 1:4 with lines lw 2.5 lc rgb "forest-green" title "OpenFOAM {/Symbol m}(T_e)", \
     EF_eV with lines lw 1.2 dt 2 lc rgb "gray30" title sprintf("E_F = %.1f eV", EF_eV)
unset arrow 2

# Plot the electron thermal conductivity.
set output "Ke_vs_Te.pdf"
set ylabel "k_e (W m^{-1} K^{-1})"
set title "Electron thermal conductivity (T_l = 300 K)"
set key top right
plot DATAFILE using 1:5 with lines lw 2.5 lc rgb "dark-violet" title "OpenFOAM k_e(T_e)"
unset title

# Plot the electron relaxation time.
set output "tau_vs_Te.pdf"
set ylabel "{/Symbol t}_e (s)"
set title "Electron relaxation time (T_l = 300 K)"
set logscale y
set key top right
plot DATAFILE using 1:6 with lines lw 2.5 lc rgb "orange-red" title "OpenFOAM {/Symbol t}_e(T_e)"
unset logscale y
unset title

# Generate the combined property panel for the manuscript.
set terminal pdfcairo enhanced font "Helvetica,10" size 7,5
set output "FEG_properties.pdf"
set multiplot layout 2,2

set key top left font ",8"
set ylabel "C_e (J m^{-3} K^{-1})"
plot DATAFILE using 1:2 with lines lw 2 lc rgb "blue" title "OpenFOAM C_e", \
     gamma_S*x with lines lw 1.2 dt 2 lc rgb "dark-red" title "{/Symbol g}T_e"

set ylabel "G (W m^{-3} K^{-1})"
plot DATAFILE using 1:3 with lines lw 2 lc rgb "red" notitle

set ylabel "{/Symbol m} (eV)"
plot DATAFILE using 1:4 with lines lw 2 lc rgb "forest-green" notitle

set ylabel "k_e (W m^{-1} K^{-1})"
plot DATAFILE using 1:5 with lines lw 2 lc rgb "dark-violet" notitle

unset multiplot

print "Figures generees : Ce_vs_Te.pdf, Ce_vs_Te_log.pdf, G_vs_Te.pdf,"
print "                   mu_vs_Te.pdf, Ke_vs_Te.pdf, tau_vs_Te.pdf,"
print "                   FEG_properties.pdf (panneau 2x2 pour l'article)"
