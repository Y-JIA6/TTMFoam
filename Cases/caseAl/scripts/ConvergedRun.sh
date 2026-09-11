#!/bin/bash
# ---------------------------------------------------------------------------
# ConvergedRun.sh
# Repeats complete simulations for one fluence until neither crater dimension
# has changed significantly during the final quiet-time window. Each trial
# restarts at t = 0, so its final result is directly reproducible.
# ---------------------------------------------------------------------------

set -euo pipefail

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
CASE_DIR=$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)
cd "$CASE_DIR"

# ---------------------------------------------------------------------------
# 1. Read and validate convergence settings
# ---------------------------------------------------------------------------
# Arguments in ps: initial end time, extension, quiet-time requirement, limit.
INITIAL_END_PS=${1:-20}
INCREMENT_PS=${2:-10}
QUIET_TIME_PS=${3:-10}
MAX_END_PS=${4:-200}

is_positive_number()
{
    awk -v value="$1" 'BEGIN { exit !(value + 0 > 0) }'
}

for value in "$INITIAL_END_PS" "$INCREMENT_PS" "$QUIET_TIME_PS" "$MAX_END_PS"; do
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
HISTORY_FILE="convergence_history/crater_quiet_convergence.tsv"
STATE_FILE="convergence_history/final_endpoint_ps"
if [ ! -f "$HISTORY_FILE" ]; then
    printf '# fluence_J_cm2\tendTime_ps\tdepth_nm\tdiameter_um\tlastDepthGrowth_ps\tlastRadiusGrowth_ps\tquietTime_ps\tstatus\n' > "$HISTORY_FILE"
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

    if [ -z "$fluence" ] || [ -z "$depth_nm" ] || [ -z "$diameter_um" ] || [ -z "$last_depth_ps" ] || [ -z "$last_radius_ps" ]; then
        echo "ERROR: final crater diagnostics could not be read from log.TTmAl." >&2
        exit 3
    fi

    last_growth_ps=$(awk -v depth="$last_depth_ps" -v radius="$last_radius_ps" 'BEGIN { print (depth > radius ? depth : radius) }')
    quiet_time_ps=$(awk -v end="$end_ps" -v growth="$last_growth_ps" 'BEGIN { printf "%.12g", end-growth }')

    if awk -v quiet="$quiet_time_ps" -v required="$QUIET_TIME_PS" 'BEGIN { exit !(quiet >= required) }'; then
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\tconverged\n' \
            "$fluence" "$end_ps" "$depth_nm" "$diameter_um" "$last_depth_ps" "$last_radius_ps" "$quiet_time_ps" >> "$HISTORY_FILE"
        sed -n '$p' "$FINAL_RESULT" >> Res.txt
        cp log.TTmAl "convergence_history/log_${fluence}J_${end_ps}ps_converged"
        printf '%s\n' "$end_ps" > "$STATE_FILE"
        echo "CRATER CONVERGED: no significant depth or radius growth during the final ${quiet_time_ps} ps."
        echo "Final endTime retained in system/controlDict: ${end_ps} ps."
        exit 0
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\textended\n' \
        "$fluence" "$end_ps" "$depth_nm" "$diameter_um" "$last_depth_ps" "$last_radius_ps" "$quiet_time_ps" >> "$HISTORY_FILE"
    cp log.TTmAl "convergence_history/log_${fluence}J_${end_ps}ps_extended"
    echo "Crater remains active after ${quiet_time_ps} ps of quiet time; extending by ${INCREMENT_PS} ps."

    next_end_ps=$(awk -v end="$end_ps" -v increment="$INCREMENT_PS" 'BEGIN { printf "%.12g", end+increment }')

    if ! awk -v candidate="$next_end_ps" -v maximum="$MAX_END_PS" 'BEGIN { exit !(candidate <= maximum) }'; then
        echo "ERROR: crater endpoints did not converge before the ${MAX_END_PS} ps safety limit." >&2
        exit 4
    fi
    end_ps=$next_end_ps
done
