#!/bin/bash
# ---------------------------------------------------------------------------
# File: Onerun.sh
# Purpose: Clean, mesh, decompose, and solve one aluminium TTM case on four
#          MPI ranks. The solver must be compiled before this script is run.
# ---------------------------------------------------------------------------

# Stop immediately if a command fails.
set -e

# ---------------------------------------------------------------------------
# 1. Start wall-clock timing
# ---------------------------------------------------------------------------
START_TOTAL=$(date +%s)

echo "================================================="
echo "   TTM ALUMINIUM ABLATION RUN                    "
echo "================================================="

# ---------------------------------------------------------------------------
# 2. Clean the previous run
# ---------------------------------------------------------------------------
echo "[1/4] Preparing the case directory"
#rm -f Res.txt
rm -f log.*
foamListTimes -rm
rm -rf processor* 2>/dev/null || echo "Locked processor directories were kept."
rm -rf postProcessing/ 2>/dev/null || true

# ---------------------------------------------------------------------------
# 3. Build and decompose the mesh
# ---------------------------------------------------------------------------
echo "[2/4] Creating the mesh (blockMesh)"
blockMesh > log.blockMesh 2>&1

# 3. DÉCOMPOSITION
echo "[3/4] Decomposing the domain"
decomposePar -force > log.decomposePar 2>&1

# ---------------------------------------------------------------------------
# 4. Solve in parallel
# ---------------------------------------------------------------------------
echo "[4/4] Solving with TTmAl in parallel"
mpirun -np 4 TTmAl -parallel | tee log.TTmAl | grep --line-buffered -E "Time =|PGC|ExecutionTime"

# Reconstruction is intentionally disabled. The probe function object writes
# the histories used by the supplied plotting scripts during the solver run.
#echo "[optional] Reconstructing fields (reconstructPar)"
#START_RECON=$(date +%s)
#reconstructPar > log.reconstructPar 2>&1 || true
#END_RECON=$(date +%s)
#TIME_RECON=$((END_RECON - START_RECON))
#echo "      Temps de reconstruction : ${TIME_RECON} s"

#echo "[optional] Extracting probe data (postProcess)"
#START_POST=$(date +%s)
#postProcess -func probes > log.postProcess
#END_POST=$(date +%s)
#TIME_POST=$((END_POST - START_POST))
#echo "      Temps d'extraction : ${TIME_POST} s"

#echo "[optional] Plotting temperature histories"
#gnuplot plots/plot_Temperature.gp

# ---------------------------------------------------------------------------
# 5. Report elapsed wall-clock time
# ---------------------------------------------------------------------------
END_TOTAL=$(date +%s)
TIME_TOTAL=$((END_TOTAL - START_TOTAL))

echo "================================================="
echo "   Calculation completed. Probe data are available."
echo "   TOTAL EXECUTION TIME: ${TIME_TOTAL} s        "
echo "================================================="
