#!/bin/bash
#------------------------------------------------------------------------------
#  runMMS_time.sh — TEMPORAL convergence study (fixed mesh, refine deltaT).
#  Backward-Euler (ddt default in OpenFOAM) is 1st order  -> expect p ~ 1.
#  If your fvSchemes uses CrankNicolson -> expect p ~ 2.
#-----------------------------------------------------------------------------
set -e
SOLVER=TTmAl

# fixed, well-resolved mesh
sed -i -E "s/hex \(0 1 2 0 3 4 5 3\) \([0-9]+ 1 [0-9]+\)/hex (0 1 2 0 3 4 5 3) (200 1 300)/" system/blockMeshDict
sed -i -E "s/simpleGrading \([0-9.]+ [0-9.]+ [0-9.]+\)/simpleGrading (1 1 1)/" system/blockMeshDict
blockMesh > log.blockMesh 2>&1

> mms_time.dat

# endTime is fixed; we shrink deltaT. List = successive halvings.
for DT in 2e-13 1e-13 5e-14 2.5e-14 1.25e-14
do
    echo "--- deltaT = $DT ---"
    # set deltaT and make sure adjustTimeStep is off
    sed -i -E "s/^deltaT[[:space:]]+.*/deltaT          $DT;/" system/controlDict
    sed -i -E "s/^adjustTimeStep[[:space:]]+.*/adjustTimeStep  no;/" system/controlDict

    foamListTimes -rm > /dev/null 2>&1 || true
    rm -f mms_error.dat
    $SOLVER > log.$SOLVER.dt$DT 2>&1

    LAST=$(tail -1 mms_error.dat)
    L2TE=$(echo "$LAST" | awk '{print $2}')
    L2TL=$(echo "$LAST" | awk '{print $3}')
    echo "$DT   $L2TE   $L2TL" >> mms_time.dat
    echo "  L2(Te)=$L2TE  L2(Tl)=$L2TL"
done

echo ""
echo "--- mms_time.dat (deltaT  L2Te  L2Tl) ---"
cat mms_time.dat
echo ""
echo "Observed temporal order p (L2 Te) :"
awk 'NR>1{ p=log(prevE/$2)/log(prevDt/$1);
           printf "   dt %s -> %s : p = %.3f\n", prevDt,$1,p }
     {prevDt=$1; prevE=$2}' mms_time.dat