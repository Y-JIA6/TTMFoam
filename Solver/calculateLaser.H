// ---------------------------------------------------------------------------
//  calculateLaser.H
//  Builds the volumetric laser source Qlaser [W/m^3] at every time step.
//
//    Qlaser(r,z,t) = S(r,z) T(t)
//
//  T(t) is a normalized Gaussian temporal envelope with FWHM tau_p. S(r,z)
//  combines the Gaussian beam profile with depth-resolved Beer--Lambert
//  absorption. The absorbed surface fraction is Drude A = 1 - R, without a
//  fluence-dependent empirical correction. Reflection is evaluated on the
//  irradiated surface, while delta_eff = delta_opt + vF*tau_e is evaluated
//  locally in every cell for the attenuation path.
//
//  Cell-averaged Beer--Lambert absorption uses a sinh(x)/x correction. Each
//  mesh column is traversed from the top surface through internal faces.
// ---------------------------------------------------------------------------

scalar t   = runTime.value();
scalar tp  = tau.value();
scalar ln2 = Foam::log(2.0);

const label topPatchID = mesh.boundaryMesh().findPatchID("top");
const polyPatch& topPatch = mesh.boundaryMesh()[topPatchID];
const fvPatchScalarField& Te_top = Te.boundaryField()[topPatchID];
const fvPatchScalarField& Tl_top = Tl.boundaryField()[topPatchID];
List<scalar> tauE_surface(topPatch.size());
List<OpticalProps> optical_surface(topPatch.size());

// ---------------------------------------------------------------------------
// 1. Surface optical diagnostics
// ---------------------------------------------------------------------------
{
    const scalar lambda_srf = lambda.value();

    forAll(topPatch, faceI)
    {
        tauE_surface[faceI] = calculateTauE(Te_top[faceI], Tl_top[faceI], Ae_val, Bl_val);
        optical_surface[faceI] = calculateOpticalProperties
        (
            tauE_surface[faceI], lambda_srf, tauERef
        );

        const label c = topPatch.faceCells()[faceI];
        R_field[c]        = optical_surface[faceI].R;
        deltaOpt_field[c] = optical_surface[faceI].delta_opt;
        deltaEff_field[c] = optical_surface[faceI].delta_opt + vF_val*tauE_surface[faceI];
        tauE_field[c]     = tauE_surface[faceI];
    }
}

// ---------------------------------------------------------------------------
// 2. Normalized Gaussian temporal envelope
// ---------------------------------------------------------------------------
if (t > 4.0 * tp)
{
    Qlaser = 0.0 * Qlaser;
}
else
{
    // Center the normalized Gaussian at two pulse widths.
    scalar T_t = Foam::sqrt(4.0 * ln2 / (constant::mathematical::pi * Foam::sqr(tp)))
               * Foam::exp(-4.0 * ln2 * Foam::sqr((t - 2.0*tp)/tp));

    const volVectorField& centres = mesh.C();
    scalar surfaceZ    = mesh.bounds().max().z();

    scalar w0_val      = w0.value();
    scalar lambda_val  = lambda.value();
    scalar F_laser_val = F_laser.value();          // Peak fluence in J/m2.
    scalar zR          = constant::mathematical::pi * Foam::sqr(w0_val) / lambda_val; // Rayleigh length.

    Qlaser = 0.0 * Qlaser;

    forAll(topPatch, faceI)
    {
        // -------------------------------------------------------------------
        // 3. Surface absorption and cell-local penetration depth
        // -------------------------------------------------------------------
        // Reuse the surface response to set the reflected and absorbed power.
        const OpticalProps& opt_surface = optical_surface[faceI];
        scalar A_drude = 1.0 - opt_surface.R;

        label currentCell = topPatch.faceCells()[faceI];

        scalar integ_d = 0.0;   // Accumulate the optical depth along the column.

        while (currentCell != -1)   // Descend from the surface to the bottom boundary.
        {
            scalar r_i = Foam::sqrt(Foam::sqr(centres[currentCell].x())
                                  + Foam::sqr(centres[currentCell].y()));
            scalar z_i = centres[currentCell].z();
            scalar z_depth = max(0.0, surfaceZ - z_i);

            // Apply Gaussian-beam divergence and its energy-conserving amplitude.
            scalar wz              = w0_val * Foam::sqrt(1.0 + Foam::sqr(z_depth/zR));
            scalar spatial_profile = Foam::exp(-2.0 * Foam::sqr(r_i/wz));
            scalar geom_factor     = Foam::sqr(w0_val)/Foam::sqr(wz);

            scalar dz_cell = dz_field[currentCell];

            const scalar tau_e_cell = calculateTauE(Te[currentCell], Tl[currentCell], Ae_val, Bl_val);
            const OpticalProps opt_cell = calculateOpticalProperties
            (
                tau_e_cell, lambda_val, tauERef
            );

            const scalar delta_opt_cell = opt_cell.delta_opt;

            const scalar delta_ball_cell = vF_val * tau_e_cell;
            const scalar delta_eff_cell = delta_opt_cell + delta_ball_cell;

            // -------------------------------------------------------------------
            // 4. Exact cell-averaged Beer--Lambert attenuation
            // -------------------------------------------------------------------
            scalar d_tau_d = dz_cell / delta_eff_cell ;       // Optical-depth increment across the cell.
            scalar tau_c_d = integ_d + 0.5*d_tau_d ;          // Optical depth at the cell center.
            scalar arg_d   = 0.5*d_tau_d ;
            scalar cor_d   = (arg_d > 1.0e-4) ? Foam::sinh(arg_d)/arg_d : 1.0;
            scalar S_d     = A_drude * (F_laser_val / delta_eff_cell)
                           * Foam::exp(-tau_c_d) * cor_d;

            // Combine depth, radial, and temporal profiles into the volumetric source.
            scalar S_r_z = S_d * geom_factor * spatial_profile;
            Qlaser[currentCell] = S_r_z * T_t;

            integ_d += d_tau_d;

            // -------------------------------------------------------------------
            // 5. Continue through the most downward internal face
            // -------------------------------------------------------------------
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
                        scalar unitZ = dirZ / area;     // Normalize the face-normal z component.
                        if (unitZ < minDirZ)
                        {
                            minDirZ = unitZ;
                            bestNextCell = potentialNext;
                        }
                    }
                }
            }
            if (minDirZ < -0.5) nextCell = bestNextCell;   // Accept only a downward connection.
            currentCell = nextCell;                         // A value of -1 ends the column.
        }
    }
}
