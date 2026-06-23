#!/bin/bash

set -e
START_TOTAL=$(date +%s)

SOLVER=TTmAl
NZ_LIST="16 32 64 128"

echo "-------------------------------------------------"
echo "   VERIFICATION MMS - Etude de convergence L2    "
echo "-------------------------------------------------"

# --- MMS switch activated --------------------------------
if ! grep -q "mmsVerification true" system/controlDict ; then
    echo "ERREUR: 'mmsVerification true;' absent from system/controlDict."
    echo "        Add it before launching the verification."
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

    # [1] Cleanning 
    rm -f log.*
    foamListTimes -rm > /dev/null 2>&1 || true
    rm -rf processor* postProcessing/ 2>/dev/null || true
    rm -f mms_error.dat

    # [2] Here we modify dynamically the mesh size in blockMeshDict (NZ) and rebuild the mesh
    sed -i -E "s/hex \(0 1 2 0 3 4 5 3\) \([0-9]+ 1 [0-9]+\)/hex (0 1 2 0 3 4 5 3) (200 1 $NZ)/" system/blockMeshDict
    sed -i -E "s/simpleGrading \([0-9.]+ [0-9.]+ [0-9.]+\)/simpleGrading (1 1 1)/" system/blockMeshDict

    # [3] Mesh series
    echo "   [maillage] blockMesh"
    blockMesh > log.blockMesh 2>&1

    # [4] Resolution 
    echo "   [calcul]   $SOLVER (serie)"
    $SOLVER > log.$SOLVER 2>&1

    # Verification of the presence of the MMS banner in the solver log
    if ! grep -q "MODE MMS IS ACTIVE" log.$SOLVER ; then
        echo "   Warning : MMS banner not found in solver log."
        echo "              Verify mmsVerification in controlDict."
    fi

    # [5] Check if mms_error.dat exists and extract the last line for L2 errors
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
echo "   CONVERGENCE Table (mms_conv.dat)"
echo "-------------------------------------------------"
cat mms_conv.dat

echo ""
echo "   Local convergence order p (L2 Te) :"
awk 'NR>1{
        p = log(prevE/$2)/log($1/prevN);
        printf "     Nz %4d -> %4d : p = %.3f\n", prevN, $1, p
     }
     {prevN=$1; prevE=$2}' mms_conv.dat

# --- Plotting ---------------------------------------------------
if command -v gnuplot >/dev/null 2>&1 && [ -f spacialconvergence.gp ]; then
    echo ""
    echo "   [plot] gnuplot spacialconvergence.gp -> spacialconvergence.png"
    gnuplot plots/spacialconvergence.gp 2>/dev/null || echo "   (Warning) gnuplot failed to generate the plot."
fi

END_TOTAL=$(date +%s)
echo ""
echo "-------------------------------------------------"
echo "   VERIFICATION MMS - FINISHED"
echo "   TOTAL TIME : $((END_TOTAL - START_TOTAL)) s"
echo "   Results   : mms_conv.dat  (+ spacialconvergence.png)"
echo "-------------------------------------------------"