// ------------------------------------------------------------------------------
//  calculateLaser.H
//  Builds the laser source term Qlaser [W/m^3] each time step.
//
//    Qlaser(r,z,t) = S(r,z) * T(t)
//
//  T(t) : normalized Gaussian temporal envelope (FWHM = tau_p)
//  S(r,z): Gaussian beam profile x depth-resolved Beer-Lambert absorption
//
//  Absorption = Drude (1 - R) plus an empirical fluence-dependent term
//  A_rel(F) that lumps the non-Drude (collisionless) channels, following the
//  phenomenological approach of Omenaca et al. (their Fig. A1).
//
//  Depth deposition uses an exact cell-averaged Beer-Lambert factor
//  (sinh(x)/x correction), walked column-by-column from the surface down via
//  an O(Nz) topological face traversal.
// ----------------------------------------------------------------------------

scalar t   = runTime.value();
scalar tp  = tau.value();
scalar ln2 = Foam::log(2.0);

// Laser off beyond t = 4*tp (Gaussian tail negligible): zero the source.
if (t > 4.0 * tp)
{
    Qlaser = 0.0 * Qlaser;
}
else
{
    // --- Temporal Gaussian: peak at t=2tp, FWHM=tp, integral over t = 1 ---
    scalar T_t = Foam::sqrt(4.0 * ln2 / (constant::mathematical::pi * Foam::sqr(tp)))
               * Foam::exp(-4.0 * ln2 * Foam::sqr((t - 2.0*tp)/tp));

    const volVectorField& centres = mesh.C();
    scalar surfaceZ    = mesh.bounds().max().z();

    scalar w0_val      = w0.value();
    scalar lambda_val  = lambda.value();
    scalar F_laser_val = F_laser.value();          // [J/m^2]
    scalar zR          = constant::mathematical::pi * Foam::sqr(w0_val) / lambda_val; // Rayleigh length

    Qlaser = 0.0 * Qlaser;

    // --- Surface patch (laser-facing boundary) ---
    label topPatchID = mesh.boundaryMesh().findPatchID("top");
    const polyPatch& topPatch        = mesh.boundaryMesh()[topPatchID];
    const fvPatchScalarField& Te_top = Te.boundaryField()[topPatchID];
    const fvPatchScalarField& Tl_top = Tl.boundaryField()[topPatchID];

    forAll(topPatch, faceI)
    {
        // --- Surface optical properties (Drude, at surface Te) ---
        scalar Te_face   = Te_top[faceI];
        scalar Tl_face   = Tl_top[faceI];
        scalar tau_e_loc = calculateTauE(Te_face, Tl_face, Ae_val, Bl_val);
        OpticalProps opt_loc = calculateOpticalProperties(tau_e_loc, lambda_val, vF_val, Te_face);

        scalar A_drude = 1.0 - opt_loc.R;   // Drude absorbed fraction

        label currentCell = topPatch.faceCells()[faceI];
        scalar integ_d = 0.0;   // cumulative optical depth tau = int dz/delta

        while (currentCell != -1)   // descend the column to the bottom boundary
        {
            scalar r_i = Foam::sqrt(Foam::sqr(centres[currentCell].x())
                                  + Foam::sqr(centres[currentCell].y()));
            scalar z_i = centres[currentCell].z();
            scalar z_depth = max(0.0, surfaceZ - z_i);

            // Gaussian beam: divergence w(z) and energy-conserving factor (w0/w)^2
            scalar wz              = w0_val * Foam::sqrt(1.0 + Foam::sqr(z_depth/zR));
            scalar spatial_profile = Foam::exp(-2.0 * Foam::sqr(r_i/wz));
            scalar geom_factor     = Foam::sqr(w0_val)/Foam::sqr(wz);

            scalar dz_cell        = dz_field[currentCell];
            scalar delta_eff_cell = opt_loc.delta_eff;

            // --- Beer-Lambert with exact cell-average (sinh(x)/x) correction ---
            scalar d_tau_d = dz_cell / delta_eff_cell ;       // cell effective depth
            scalar tau_c_d = integ_d + 0.5*d_tau_d ;          // depth at cell center
            scalar arg_d   = 0.5*d_tau_d ;
            scalar cor_d   = (arg_d > 1.0e-4) ? Foam::sinh(arg_d)/arg_d : 1.0;
            scalar S_d     = A_drude * (F_laser_val / delta_eff_cell)
                           * Foam::exp(-tau_c_d) * cor_d;

            // --- Source [W/m^3] = depth profile x beam profile x temporal ---
            scalar S_r_z = S_d * geom_factor * spatial_profile;
            Qlaser[currentCell] = S_r_z * T_t;

            integ_d += d_tau_d;

            // --- Topological descent: pick the most downward internal face ---
            const cell& c = mesh.cells()[currentCell];
            label nextCell = -1;
            scalar minDirZ = GREAT;
            label bestNextCell = -1;
            forAll(c, i)
            {
                label internalFaceID = c[i];
                if (mesh.isInternalFace(internalFaceID))
                {
                    vector sf = mesh.Sf()[internalFaceID];
                    scalar dirZ = 0.0;
                    label potentialNext = -1;
                    if (mesh.faceOwner()[internalFaceID] == currentCell)
                    {
                        dirZ = sf.z();
                        potentialNext = mesh.faceNeighbour()[internalFaceID];
                    }
                    else
                    {
                        dirZ = -sf.z();
                        potentialNext = mesh.faceOwner()[internalFaceID];
                    }
                    scalar area = Foam::mag(sf);
                    if (area > VSMALL)
                    {
                        scalar unitZ = dirZ / area;     // normalized z-component
                        if (unitZ < minDirZ)
                        {
                            minDirZ = unitZ;
                            bestNextCell = potentialNext;
                        }
                    }
                }
            }
            if (minDirZ < -0.5) nextCell = bestNextCell;   // accept only downward
            currentCell = nextCell;                         // -1 ends the column
        }
    }
}