#!/usr/bin/gnuplot
# =========================================================================
# Script Gnuplot: Transient optical and thermophysical response
# Generates surface-probe plots of optical, thermal, and kinetic properties.
# =========================================================================

COL = 3          # Column 2 is the surface-axis probe and column 3 is the first-cell-center probe.

# Select the probe directory that contains the optical diagnostics.
P = system("for d in postProcessing/probes/*/; do [ -f \"$d/R_field\" ] && { printf %s \"$d\"; break; }; done")

if (strlen(P) == 0) {
    print "ERREUR : aucun dossier postProcessing/probes/*/ ne contient R_field."
    print "         Verifier que le bloc functions{} de controlDict est actif"
    print "         et que le solveur a bien ete recompile."
    exit
}
print "Lecture des sondes dans : ".P

set terminal pdfcairo enhanced font "Helvetica,12" size 6,4

# Use logarithmic time in ps to resolve sub-picosecond dynamics.
set logscale x
set xlabel "t (ps)"
set xrange [8e-4:20.6]
set grid
set key top left
set format x "10^{%T}"

# Plot reflectivity and penetration depths.
set output "Reflectivity_vs_t.pdf"
set ylabel "R (-)"
set key bottom left
plot P."R_field" using ($1*1e12):COL with lines lw 2.5 lc rgb "dark-orange" \
        title "OpenFOAM R(t)"

set output "Penetration_vs_t.pdf"
set ylabel "penetration depth (nm)"
set key top left
plot P."deltaOpt_field" using ($1*1e12):(column(COL)*1e9) with lines lw 2.5 \
        lc rgb "dark-green" title "OpenFOAM {/Symbol d}_{opt}", \
     P."deltaEff_field" using ($1*1e12):(column(COL)*1e9) with lines lw 2.5 dt 2 \
        lc rgb "dark-violet" title "OpenFOAM {/Symbol d}_{eff}"

# Plot thermal and kinetic properties.
set output "Properties_vs_t.pdf"
set ylabel "C_e (J m^{-3} K^{-1})   /   k_e (W m^{-1} K^{-1})"
set logscale y
set format y "10^{%T}"
plot P."Ce_var" using ($1*1e12):COL with lines lw 2.5 lc rgb "blue"        title "OpenFOAM C_e", \
     P."ke_var" using ($1*1e12):COL with lines lw 2.5 lc rgb "dark-violet" title "OpenFOAM k_e"
unset logscale y
set format y "%g"

set output "Coupling_vs_t.pdf"
set ylabel "G (W m^{-3} K^{-1})"
plot P."G_var" using ($1*1e12):COL with lines lw 2.5 lc rgb "red" title "OpenFOAM G(T_e)"

set output "tauE_vs_t.pdf"
set ylabel "{/Symbol t}_e (s)"
set logscale y
set format y "10^{%T}"
plot P."tauE_field" using ($1*1e12):COL with lines lw 2.5 lc rgb "orange-red" \
        title "OpenFOAM {/Symbol t}_e(T_e,T_l)"
unset logscale y
set format y "%g"

# Generate the combined optical and thermal panel for the manuscript.
set terminal pdfcairo enhanced font "Helvetica,10" size 7,5
set output "Transient_response.pdf"
set multiplot layout 2,2
set key top left font ",8"

set ylabel "R (-)"
set key bottom left font ",8"
plot P."R_field" using ($1*1e12):COL with lines lw 2 lc rgb "dark-orange" \
        title "OpenFOAM reflectivity"

set ylabel "{/Symbol d} (nm)"
set key top left font ",8"
plot P."deltaOpt_field" using ($1*1e12):(column(COL)*1e9) with lines lw 2 \
        lc rgb "dark-green" title "OpenFOAM {/Symbol d}_{opt}", \
     P."deltaEff_field" using ($1*1e12):(column(COL)*1e9) with lines lw 2 dt 2 \
        lc rgb "dark-violet" title "OpenFOAM {/Symbol d}_{eff}"

set ylabel "C_e , k_e"
set logscale y ; set format y "10^{%T}"
plot P."Ce_var" using ($1*1e12):COL with lines lw 2 lc rgb "blue"        title "OpenFOAM C_e", \
     P."ke_var" using ($1*1e12):COL with lines lw 2 lc rgb "dark-violet" title "OpenFOAM k_e"
unset logscale y ; set format y "%g"

set ylabel "{/Symbol t}_e (s)"
set logscale y ; set format y "10^{%T}"
plot P."tauE_field" using ($1*1e12):COL with lines lw 2 lc rgb "orange-red" \
        title "OpenFOAM {/Symbol t}_e"
unset logscale y ; set format y "%g"

unset multiplot

print "Figures generees : Reflectivity_vs_t.pdf, Penetration_vs_t.pdf,"
print "                   Properties_vs_t.pdf, Coupling_vs_t.pdf, tauE_vs_t.pdf,"
print "                   Transient_response.pdf (panneau 2x2 pour l'article)"
