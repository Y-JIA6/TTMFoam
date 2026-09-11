#!/bin/bash
# ---------------------------------------------------------------------------
# File: Allrun.sh
# Purpose: Run the adaptive aluminium campaign from 1 to 9 J/cm2 and retain
#          the accepted crater result and probe histories for each fluence.
# ---------------------------------------------------------------------------

# Stop immediately if a command fails.
set -e

# Run from the case root regardless of the caller's current directory.
SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
CASE_DIR=$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)
cd "$CASE_DIR"

# Source OpenFOAM when it is not already available in the current shell.
if ! command -v foamDictionary >/dev/null 2>&1; then
    if [ -f /opt/openfoam10/etc/bashrc ]; then
        # shellcheck disable=SC1091
        source /opt/openfoam10/etc/bashrc
    fi
fi

for command_name in foamDictionary blockMesh decomposePar mpirun; do
    if ! command -v "$command_name" >/dev/null 2>&1; then
        echo "ERROR: required command not found: $command_name" >&2
        echo "Source the OpenFOAM environment before running this case." >&2
        exit 1
    fi
done

if ! command -v TTmAl >/dev/null 2>&1; then
    if [ -x "$CASE_DIR/../../platforms/linux64GccDPInt32Opt/bin/TTmAl" ]; then
        export PATH="$CASE_DIR/../../platforms/linux64GccDPInt32Opt/bin:$PATH"
    fi
fi

if ! command -v TTmAl >/dev/null 2>&1; then
    echo "ERROR: TTmAl executable was not found in PATH." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# 1. Prepare campaign outputs
# ---------------------------------------------------------------------------
rm -rf resultats_fluences
mkdir -p resultats_fluences
rm -f Res.txt

# ---------------------------------------------------------------------------
# 2. Define the fluence sequence and convergence settings
# ---------------------------------------------------------------------------
# Absorption is purely collisional: A = 1 - R(Te,Tl), with no empirical,
# fluence-dependent absorption correction.
fluences=(1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0)

# The first fluence begins at 20 ps. Each following fluence starts ten
# picoseconds later than the final endpoint retained for the preceding case.
# A non-converged case is extended in ten-picosecond increments up to 200 ps.
initial_end_ps=20
end_increment_ps=10
quiet_time_ps=10
max_end_ps=200

for F in "${fluences[@]}"
do
    echo "================================================="
    echo " FLUENCE = $F J/cm2 | initial endTime = ${initial_end_ps} ps "
    echo "================================================="

    # Convert J/cm2 to J/m2 for the OpenFOAM dictionary.
    F_m2=$(awk "BEGIN {print $F * 10000}")

    # Apply the fluence; ConvergedRun.sh selects the final end time automatically.
    foamDictionary constant/laserProperties -entry F_laser -set $F_m2

    # Run the current fluence until the crater is quiet, then pass its final
    # endpoint plus ten picoseconds to the following fluence.
    bash "$SCRIPT_DIR/ConvergedRun.sh" "$initial_end_ps" "$end_increment_ps" "$quiet_time_ps" "$max_end_ps"
    if [ ! -s convergence_history/final_endpoint_ps ]; then
        echo "ERROR: final endpoint was not written by ConvergedRun.sh." >&2
        exit 1
    fi
    final_end_ps=$(< convergence_history/final_endpoint_ps)
    initial_end_ps=$(awk -v final="$final_end_ps" -v increment="$end_increment_ps" 'BEGIN { printf "%.12g", final+increment }')
    echo " Final endpoint for ${F} J/cm2: ${final_end_ps} ps; next fluence starts at ${initial_end_ps} ps."

    # -----------------------------------------------------------------------
    # 3. Archive accepted probe histories and the corresponding solver log
    # -----------------------------------------------------------------------
    # Probe folders use the formatted initial time name (for example,
    # 0.000000e+00). Select the folder created by the active probe function.
    PROBEDIR=""
    for d in postProcessing/probes/*/; do
        [ -f "$d/R_field" ] && { PROBEDIR="$d"; break; }
    done
    if [ -z "$PROBEDIR" ]; then
        echo "WARNING: probe files were not found for F = $F (is functions{} active?)."
    else
        for f in Te Tl Qlaser Ce_var ke_var G_var R_field deltaOpt_field deltaEff_field tauE_field
        do
            [ -f "$PROBEDIR/$f" ] && cp "$PROBEDIR/$f" resultats_fluences/${f}_${F}
        done
    fi
    cp log.TTmAl resultats_fluences/log_${F}
done

if command -v gnuplot >/dev/null 2>&1; then
    echo " Plotting temperature and ablation figures"
    gnuplot plots/plot_All.gp
    gnuplot plots/plot_Ablation.gp
else
    echo "WARNING: gnuplot was not found; calculations completed without figures."
fi

echo "================================================="
echo " CAMPAIGN COMPLETED "
echo "================================================="
