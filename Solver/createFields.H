// ===========================================================================
//  createFields.H
//  Creates the primary temperature fields, laser source, material-property
//  fields, surface optical diagnostics, and MMS verification sources.
//
//  Te and Tl are read from the selected start time. Thermophysical fields are
//  seeded from ttmProperties and updated in every Picard iteration. dz_field
//  is computed once and used by calculateLaser.H for Beer--Lambert integration.
// ===========================================================================

Info<< "Initialisation des champs (Te, Tl, Qlaser)\n" << endl;

// ---------------------------------------------------------------------------
// 1. Primary temperature fields
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// 2. Volumetric laser source
// ---------------------------------------------------------------------------
volScalarField Qlaser
(
    IOobject("Qlaser", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("Qlaser", dimPower/dimVolume, 0.0)   // Volumetric power density in W/m3.
);

// ---------------------------------------------------------------------------
// 3. Temperature-dependent material-property fields
// ---------------------------------------------------------------------------
volScalarField Ce_var
(
    IOobject("Ce_var", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh, Ce
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
    IOobject("G_var", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh, G
);

// ---------------------------------------------------------------------------
// 4. Surface optical diagnostics sampled by the probe function object
// ---------------------------------------------------------------------------
volScalarField R_field
(
    IOobject("R_field", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("R", dimless, 0.0)
);

volScalarField deltaOpt_field
(
    IOobject("deltaOpt_field", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("deltaOpt", dimLength, 0.0)
);

volScalarField deltaEff_field
(
    IOobject("deltaEff_field", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("deltaEff", dimLength, 0.0)
);

volScalarField tauE_field
(
    IOobject("tauE_field", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("tauE", dimTime, 0.0)
);

// ---------------------------------------------------------------------------
// 5. Per-cell vertical thickness for Beer--Lambert integration
// ---------------------------------------------------------------------------
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

    // Compute the vertical span from the cell vertices.
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
    dz_field[cellI] = max(maxZ_i - minZ_i, 1.0e-12);  // Enforce a positive minimum thickness.
}

Info << ">>Champ dz calculé et stocké dans dz_field\n" << endl;

// ---------------------------------------------------------------------------
// 6. Method-of-Manufactured-Solutions source fields
// ---------------------------------------------------------------------------
volScalarField Qe_mms
(
    IOobject("Qe_mms", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("z", dimEnergy/dimVolume/dimTime, 0.0)
);
volScalarField Ql_mms("Ql_mms", Qe_mms);   // Copy dimensions and initial values from Qe_mms.

scalar initialisationTime = runTime.elapsedCpuTime();
Info << ">>Initialisation des champs terminée en " << initialisationTime << " s\n" << endl;
