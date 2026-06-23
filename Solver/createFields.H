// ============================================================================
//  createFields.H
//  Declares all fields used by the solver:
//    - primary unknowns Te, Tl (read from 0/, written each output)
//    - the laser source Qlaser
//    - per-cell property fields (*_var) updated each step
//    - geometric helper dz_field (precomputed once)
//    - MMS verification sources Qe_mms, Ql_mms
//  Property fields are seeded from the dictionary values read in
//  readProperties.H, then overwritten per-cell during the Picard loop.
// ============================================================================

Info<< "Initialisation des champs (Te, Tl, Qlaser)\n" << endl;

// --- 1. Primary temperature fields (read from 0/, auto-written) --------------
volScalarField Te
(
    IOobject("Te", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

volScalarField Tl
(
    IOobject("Tl", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
    mesh
);

// --- 2. Laser source term [W/m^3] --------------------------------------------
volScalarField Qlaser
(
    IOobject("Qlaser", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("Qlaser", dimPower/dimVolume, 0.0)   // W/m^3
);

// --- 3. Per-cell dynamic property fields -------------------------------------
// Seeded from dictionary constants; overwritten each Picard iteration.
volScalarField Ce_var
(
    IOobject("Ce_var", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("Ce_init", dimEnergy/dimVolume/dimTemperature, 0.0) // overwritten by table
);

volScalarField Cl_var
(
    IOobject("Cl_var", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh, Cl
);

volScalarField ke_var
(
    IOobject("ke_var", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh, ke
);

volScalarField kl_var
(
    IOobject("kl_var", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh, kl
);

volScalarField G_var
(
    IOobject("G_var", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh, G
);

// --- 4. Geometric helper: per-cell z-thickness, computed ONCE ----------------
// Used by the Beer-Lambert traversal in calculateLaser.H. Cost is paid only at
// startup, so the brute-force point scan below is acceptable.
volScalarField dz_field
(
    IOobject("dz_field", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("dz", dimLength, 1e-12)
);

const pointField& pts = mesh.points();
forAll(mesh.cells(), cellI)
{
    const cell& currentCell = mesh.cells()[cellI];
    scalar minZ_i =  GREAT;
    scalar maxZ_i = -GREAT;

    // span of the cell in z = max(point z) - min(point z)
    forAll(currentCell, fI)
    {
        const face& f = mesh.faces()[currentCell[fI]];
        forAll(f, pI)
        {
            scalar pointZ = pts[f[pI]].z();
            if (pointZ < minZ_i) minZ_i = pointZ;
            if (pointZ > maxZ_i) maxZ_i = pointZ;
        }
    }
    dz_field[cellI] = max(maxZ_i - minZ_i, 1.0e-12);  // floor avoids divide-by-zero
}

Info << ">>Champ dz calculé et stocké dans dz_field\n" << endl;

// --- 5. MMS verification sources [W/m^3] -------------------------------------
// (dimEnergy/dimVolume/dimTime == dimPower/dimVolume == W/m^3; same dimension,
//  just written differently. Only nonzero when mmsVerification is true.)
volScalarField Qe_mms
(
    IOobject("Qe_mms", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("z", dimEnergy/dimVolume/dimTime, 0.0)
);
volScalarField Ql_mms("Ql_mms", Qe_mms);   // same dimensions, value-copied

scalar initialisationTime = runTime.elapsedCpuTime();
Info << ">>Initialisation des champs terminée en " << initialisationTime << " s\n" << endl;