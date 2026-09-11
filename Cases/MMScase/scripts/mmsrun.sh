#!/bin/bash
# ---------------------------------------------------------------------------
# File: mmsrun.sh
# Purpose: Perform the MMS spatial-convergence study by refining the axial
#          mesh and recording the L2 errors for Te and Tl.
# ---------------------------------------------------------------------------

# Stop immediately if a command fails.
set -e
START_TOTAL=$(date +%s)

SOLVER=TTmAl
NZ_LIST="16 32 64 128"

echo "-------------------------------------------------"
echo "   MMS VERIFICATION - SPATIAL L2 CONVERGENCE     "
echo "-------------------------------------------------"

# ---------------------------------------------------------------------------
# 1. Confirm that MMS mode is active
# ---------------------------------------------------------------------------
if ! grep -q "mmsVerification true" system/controlDict ; then
    echo "ERROR: 'mmsVerification true;' is missing from system/controlDict."
    echo "       Enable it before launching the verification."
    exit 1
fi
echo "[OK] mmsVerification is active in controlDict"

> mms_conv.dat

for NZ in $NZ_LIST
do
    echo ""
    echo "-------------------------------------------------"
    echo "   Nz = $NZ"
    echo "-------------------------------------------------"

    # -----------------------------------------------------------------------
    # 2. Clean the previous mesh and MMS output
    # -----------------------------------------------------------------------
    rm -f log.*
    foamListTimes -rm > /dev/null 2>&1 || true
    rm -rf processor* postProcessing/ 2>/dev/null || true
    rm -f mms_error.dat

    # -----------------------------------------------------------------------
    # 3. Update the axial mesh resolution and rebuild the mesh
    # -----------------------------------------------------------------------
    sed -i -E "s/hex \(0 1 2 0 3 4 5 3\) \([0-9]+ 1 [0-9]+\)/hex (0 1 2 0 3 4 5 3) (200 1 $NZ)/" system/blockMeshDict
    sed -i -E "s/simpleGrading \([0-9.]+ [0-9.]+ [0-9.]+\)/simpleGrading (1 1 1)/" system/blockMeshDict

    echo "   [mesh] blockMesh"
    blockMesh > log.blockMesh 2>&1

    # -----------------------------------------------------------------------
    # 4. Solve the MMS case in serial
    # -----------------------------------------------------------------------
    echo "   [solve]  $SOLVER (serial)"
    $SOLVER > log.$SOLVER 2>&1

    # Confirm that the solver reported the MMS mode in its log.
    if ! grep -q "MODE MMS IS ACTIVE" log.$SOLVER ; then
        echo "   Warning : MMS banner not found in solver log."
        echo "              Verify mmsVerification in controlDict."
    fi

    # -----------------------------------------------------------------------
    # 5. Extract the final L2 errors
    # -----------------------------------------------------------------------
    if [ -f mms_error.dat ]; then
        LAST=$(tail -1 mms_error.dat)
        L2TE=$(echo "$LAST" | awk '{print $2}')
        L2TL=$(echo "$LAST" | awk '{print $3}')
        echo "$NZ   $L2TE   $L2TL" >> mms_conv.dat
        echo "   [results] L2(Te)=$L2TE   L2(Tl)=$L2TL"
    else
        echo "   ERROR : mms_error.dat not found for Nz=$NZ"
    fi
done

echo ""
echo "-------------------------------------------------"
echo "   CONVERGENCE TABLE (mms_conv.dat)"
echo "-------------------------------------------------"
cat mms_conv.dat

echo ""
echo "   Local convergence order p (L2 Te) :"
awk 'NR>1{
        p = log(prevE/$2)/log($1/prevN);
        printf "     Nz %4d -> %4d : p = %.3f\n", prevN, $1, p
     }
     {prevN=$1; prevE=$2}' mms_conv.dat

# ---------------------------------------------------------------------------
# 6. Generate the spatial-convergence plot when gnuplot is available
# ---------------------------------------------------------------------------
if command -v gnuplot >/dev/null 2>&1 && [ -f spacialconvergence.gp ]; then
    echo ""
    echo "   [plot] gnuplot spacialconvergence.gp -> spacialconvergence.png"
    gnuplot plots/spacialconvergence.gp 2>/dev/null || echo "   (Warning) gnuplot failed to generate the plot."
fi

END_TOTAL=$(date +%s)
echo ""
echo "-------------------------------------------------"
echo "   MMS VERIFICATION COMPLETED"
echo "   TOTAL TIME: $((END_TOTAL - START_TOTAL)) s"
echo "   RESULTS: mms_conv.dat (+ spacialconvergence.png)"
echo "-------------------------------------------------"
