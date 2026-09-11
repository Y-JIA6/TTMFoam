#!/bin/bash
# ---------------------------------------------------------------------------
# File: Serierun.sh
# Purpose: Clean, mesh, and solve one aluminium TTM case in serial. This is
#          intended for debugging and comparison with the MPI execution.
# ---------------------------------------------------------------------------

# Stop immediately if a command fails.
set -e

START_TOTAL=$(date +%s)

echo "================================================="
echo "   Running TTM ablation case in serial           "
echo "================================================="

# ---------------------------------------------------------------------------
# 1. Clean the previous run and build the mesh
# ---------------------------------------------------------------------------
echo "[1/4] Preparing the case directory"
rm -f log.*
foamListTimes -rm
rm -rf processor* 2>/dev/null || echo "Locked processor directories were kept."
rm -rf postProcessing/ 2>/dev/null || true

echo "[2/4] Creating the mesh (blockMesh)"
blockMesh > log.blockMesh 2>&1

# ---------------------------------------------------------------------------
# 2. Solve the case on one core
# ---------------------------------------------------------------------------
echo "[3/4] Solving with TTmAl on one core"
TTmAl | tee log.TTmAl | grep --line-buffered -E "Time =|PGC|ExecutionTime"

# ---------------------------------------------------------------------------
# 3. Extract probe data and report the elapsed time
# ---------------------------------------------------------------------------
echo "[4/4] Extracting probe data (postProcess)"
START_POST=$(date +%s)
postProcess -func probes > log.postProcess
END_POST=$(date +%s)
TIME_POST=$((END_POST - START_POST))
echo "      Probe extraction time: ${TIME_POST} s"

#echo "[optional] Plotting peak temperature histories"
#gnuplot plots/plot_Temperature.gp

END_TOTAL=$(date +%s)
TIME_TOTAL=$((END_TOTAL - START_TOTAL))

echo "================================================="
echo "   Serial calculation completed. Data are ready. "
echo "   TOTAL EXECUTION TIME: ${TIME_TOTAL} s         "
echo "================================================="
