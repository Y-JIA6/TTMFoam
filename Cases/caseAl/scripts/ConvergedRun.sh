#!/bin/bash
# ---------------------------------------------------------------------------
# ConvergedRun.sh
# Repeats complete simulations for one fluence until the crater is quiet and
# its geometry changes by less than fixed, user-visible tolerances between two
# successive endpoints. Each trial restarts at t = 0 for reproducibility.
# ---------------------------------------------------------------------------

set -euo pipefail

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
CASE_DIR=$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)
cd "$CASE_DIR"

# ---------------------------------------------------------------------------
# 1. Read and validate convergence settings
# ---------------------------------------------------------------------------

CONTROL_END_TIME_S=$(awk '$1 == "endTime" { value=$2; sub(/;.*/, "", value); print value; exit }' system/controlDict)
if [ -z "$CONTROL_END_TIME_S" ]; then
    echo "ERROR: endTime is missing from system/controlDict." >&2
    exit 2
fi
INITIAL_END_PS=$(awk -v value="$CONTROL_END_TIME_S" 'BEGIN { printf "%.12g", value*1e12 }')
if [ -n "${1:-}" ]; then
    INITIAL_END_PS=$1
fi
INCREMENT_PS=${2:-10}
QUIET_TIME_PS=${3:-1.25}
MAX_END_PS=${4:-200}
DEPTH_TOLERANCE_CELLS=${5:-1.5}
DIAMETER_TOLERANCE_CELLS=${6:-1.0}

is_positive_number()
{
    awk -v value="$1" 'BEGIN { exit !(value + 0 > 0) }'
}

for value in "$INITIAL_END_PS" "$INCREMENT_PS" "$QUIET_TIME_PS" "$MAX_END_PS" "$DEPTH_TOLERANCE_CELLS" "$DIAMETER_TOLERANCE_CELLS"; do
    if ! is_positive_number "$value"; then
        echo "ERROR: all arguments must be positive numbers." >&2
        exit 2
    fi
done

if ! awk -v initial="$INITIAL_END_PS" -v maximum="$MAX_END_PS" 'BEGIN { exit !(initial <= maximum) }'; then
    echo "ERROR: the initial end time must not exceed the safety limit." >&2
    exit 2
fi

# ---------------------------------------------------------------------------
# 2. Prepare convergence history and protect previous accepted results
# ---------------------------------------------------------------------------
mkdir -p convergence_history
HISTORY_FILE="convergence_history/crater_geometry_convergence.tsv"
STATE_FILE="convergence_history/final_endpoint_ps"
if [ ! -f "$HISTORY_FILE" ]; then
    printf '# fluence_J_cm2\tendTime_ps\tdepth_nm\tdiameter_um\tdepthMeshResolution_nm\tdiameterMeshResolution_um\tdeltaDepth_nm\tdeltaDiameter_um\tdepthTolerance_nm\tdiameterTolerance_um\tlastDepthGrowth_ps\tlastRadiusGrowth_ps\tquietTime_ps\tstatus\n' > "$HISTORY_FILE"
fi

# TTmAl appends a line to Res.txt at every trial. Preserve previous fluences
# and append only the final accepted result of this convergence sequence.
RES_BACKUP=$(mktemp "${TMPDIR:-/tmp}/ttmfoam-res-backup.XXXXXX")
FINAL_RESULT=$(mktemp "${TMPDIR:-/tmp}/ttmfoam-final-result.XXXXXX")
cleanup()
{
    rm -f "$RES_BACKUP" "$FINAL_RESULT"
}
trap cleanup EXIT

if [ -f Res.txt ]; then
    cp Res.txt "$RES_BACKUP"
else
    : > "$RES_BACKUP"
fi

restore_results()
{
    cp "$RES_BACKUP" Res.txt
}

# ---------------------------------------------------------------------------
# 3. Run fresh trials until the crater remains quiet for the required window
# ---------------------------------------------------------------------------
end_ps=$INITIAL_END_PS
previous_depth_nm=""
previous_diameter_um=""

while :; do
    end_time_s=$(awk -v value="$end_ps" 'BEGIN { printf "%.12g", value*1e-12 }')

    echo "================================================="
    echo " CRATER-CONVERGENCE RUN: endTime = ${end_ps} ps "
    echo "================================================="

    foamDictionary system/controlDict -entry startFrom -set startTime
    foamDictionary system/controlDict -entry startTime -set 0
    foamDictionary system/controlDict -entry endTime -set "$end_time_s"

    bash "$SCRIPT_DIR/Onerun.sh"

    if [ ! -s Res.txt ]; then
        echo "ERROR: TTmAl did not write Res.txt." >&2
        exit 3
    fi
    tail -n 1 Res.txt > "$FINAL_RESULT"
    restore_results

    fluence=$(awk '/Fluence laser simulee/ {value=$(NF - 1)} END {print value}' log.TTmAl)
    depth_nm=$(awk '/Profondeur NUMERIQUE/ {value=$(NF - 1)} END {print value}' log.TTmAl)
    diameter_um=$(awk '/Diametre NUMERIQUE/ {value=$(NF - 1)} END {print value}' log.TTmAl)
    last_depth_ps=$(awk '/Derniere croissance en profondeur/ {value=$6} END {print value}' log.TTmAl)
    last_radius_ps=$(awk '/Derniere croissance en rayon/ {value=$6} END {print value}' log.TTmAl)
    depth_resolution_nm=$(awk '/Depth mesh resolution/ {value=$(NF - 1)} END {print value}' log.TTmAl)
    diameter_resolution_um=$(awk '/Diameter mesh resolution/ {value=$(NF - 1)} END {print value}' log.TTmAl)

    if [ -z "$fluence" ] || [ -z "$depth_nm" ] || [ -z "$diameter_um" ] || [ -z "$last_depth_ps" ] || [ -z "$last_radius_ps" ] || [ -z "$depth_resolution_nm" ] || [ -z "$diameter_resolution_um" ]; then
        echo "ERROR: final crater diagnostics could not be read from log.TTmAl." >&2
        exit 3
    fi

    last_growth_ps=$(awk -v depth="$last_depth_ps" -v radius="$last_radius_ps" 'BEGIN { print (depth > radius ? depth : radius) }')
    quiet_time_ps=$(awk -v end="$end_ps" -v growth="$last_growth_ps" 'BEGIN { printf "%.12g", end-growth }')

    # Compare the current geometry with the preceding endpoint. 
    endpoint_status="baseline"
    delta_depth_nm="NA"
    delta_diameter_um="NA"
    depth_tolerance_nm="NA"
    diameter_tolerance_um="NA"
    endpoint_stable=0
    if [ -n "$previous_depth_nm" ]; then
        delta_depth_nm=$(awk -v current="$depth_nm" -v previous="$previous_depth_nm" 'BEGIN { value=current-previous; print (value < 0 ? -value : value) }')
        delta_diameter_um=$(awk -v current="$diameter_um" -v previous="$previous_diameter_um" 'BEGIN { value=current-previous; print (value < 0 ? -value : value) }')
        depth_tolerance_nm=$(awk -v cells="$DEPTH_TOLERANCE_CELLS" -v resolution="$depth_resolution_nm" 'BEGIN { printf "%.6g", cells*resolution }')
        diameter_tolerance_um=$(awk -v cells="$DIAMETER_TOLERANCE_CELLS" -v resolution="$diameter_resolution_um" 'BEGIN { printf "%.6g", cells*resolution }')

        if awk -v depth="$delta_depth_nm" -v diameter="$delta_diameter_um" -v depthTol="$depth_tolerance_nm" -v diameterTol="$diameter_tolerance_um" 'BEGIN { exit !(depth <= depthTol && diameter <= diameterTol) }'; then
            endpoint_status="geometry-stable"
            endpoint_stable=1
        else
            endpoint_status="geometry-changing"
        fi
    fi

    if [ "$endpoint_stable" -eq 1 ] && awk -v quiet="$quiet_time_ps" -v required="$QUIET_TIME_PS" 'BEGIN { exit !(quiet >= required) }'; then
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tconverged\n' \
            "$fluence" "$end_ps" "$depth_nm" "$diameter_um" "$depth_resolution_nm" "$diameter_resolution_um" "$delta_depth_nm" "$delta_diameter_um" "$depth_tolerance_nm" "$diameter_tolerance_um" "$last_depth_ps" "$last_radius_ps" "$quiet_time_ps" >> "$HISTORY_FILE"
        sed -n '$p' "$FINAL_RESULT" >> Res.txt
        cp log.TTmAl "convergence_history/log_${fluence}J_${end_ps}ps_converged"
        printf '%s\n' "$end_ps" > "$STATE_FILE"
        echo "CRATER CONVERGED: final endpoint is geometry-stable within ${depth_tolerance_nm} nm (${DEPTH_TOLERANCE_CELLS} depth cells) and ${diameter_tolerance_um} um (${DIAMETER_TOLERANCE_CELLS} diameter cells), and the crater was quiet for ${quiet_time_ps} ps."
        echo "Final endTime retained in system/controlDict: ${end_ps} ps."
        exit 0
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\textended-%s\n' \
        "$fluence" "$end_ps" "$depth_nm" "$diameter_um" "$depth_resolution_nm" "$diameter_resolution_um" "$delta_depth_nm" "$delta_diameter_um" "$depth_tolerance_nm" "$diameter_tolerance_um" "$last_depth_ps" "$last_radius_ps" "$quiet_time_ps" "$endpoint_status" >> "$HISTORY_FILE"
    cp log.TTmAl "convergence_history/log_${fluence}J_${end_ps}ps_extended"
    echo "Crater status: ${endpoint_status}; quiet time = ${quiet_time_ps} ps. Extending by ${INCREMENT_PS} ps."

    previous_depth_nm="$depth_nm"
    previous_diameter_um="$diameter_um"

    next_end_ps=$(awk -v end="$end_ps" -v increment="$INCREMENT_PS" 'BEGIN { printf "%.12g", end+increment }')

    if ! awk -v candidate="$next_end_ps" -v maximum="$MAX_END_PS" 'BEGIN { exit !(candidate <= maximum) }'; then
        echo "ERROR: crater endpoints did not converge before the ${MAX_END_PS} ps safety limit." >&2
        exit 4
    fi
    end_ps=$next_end_ps
done
