// ===========================================================================
//  TTmAl.C
//  Solves the two-temperature model for ultrashort-pulse laser heating of
//  aluminum on a 2D-axisymmetric OpenFOAM mesh.
//
//  The physical mode couples Fermi--Dirac electron properties, Shomate lattice
//  thermodynamics, local Drude--Beer--Lambert absorption, and phase-explosion
//  crater diagnostics. The MMS mode verifies the discretization independently.
// ===========================================================================

#include "fvCFD.H"
#include <fstream>
#include "updateThermo.H"
#include "updateOptical.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info << "Lecture des dictionnaires de propriétés\n" << endl;

    #include "readProperties.H"
    #include "createFields.H"

    // Cache dictionary values used repeatedly in cell loops.
    scalar gamma_val = gamma.value();
    scalar Ae_val = Ae.value();
    scalar Bl_val = Bl.value();
    scalar vF_val = vF.value();

    // Compute the ambient relaxation time used by the optical model.
    const scalar tauERef = calculateTauE(T0.value(), T0.value(), Ae_val, Bl_val);
    Info<< ">> tau_e(T0) de reference pour le modele optique = "
        << tauERef << " s" << endl;

    scalar shA_sol_val = ShA_sol.value();
    scalar shB_sol_val = ShB_sol.value();
    scalar shC_sol_val = ShC_sol.value();
    scalar shD_sol_val = ShD_sol.value();
    scalar shE_sol_val = ShE_sol.value();

    scalar shA_liq_val = ShA_liq.value();
    scalar shB_liq_val = ShB_liq.value();
    scalar shC_liq_val = ShC_liq.value();
    scalar shD_liq_val = ShD_liq.value();
    scalar shE_liq_val = ShE_liq.value();

    scalar cl_mult_val = Cl_mult.value();

    #include "createTables.H"   // Build the Ce, G, and chemical-potential tables once.

    // Read the Picard iteration controls from controlDict.
    const label maxIter = runTime.controlDict().lookupOrDefault<label>("picardMaxIter", 15);
    const scalar picardAbsTol = runTime.controlDict().lookupOrDefault<scalar>("picardAbsTol", 1.0e-3);  // K
    const scalar picardRelTol = runTime.controlDict().lookupOrDefault<scalar>("picardRelTol", 1.0e-6);  // -

    Info<< ">> Picard controls : maxIter = " << maxIter
        << ", absTol = " << picardAbsTol << " K"
        << ", relTol = " << picardRelTol << endl;

    // Initialize the phase-explosion ablation diagnostics.
    scalar T_ablation = 6700.0; // critical temperature Tc [K]
    scalar profondeur = 0.0;    // crater depth  [m]
    // Record the latest depth and radius growth times.
    scalar profPrevGlobal    = 0.0;
    scalar rayPrevGlobal     = 0.0;
    scalar tLastDepthGrowth  = 0.0;
    scalar tLastRadiusGrowth = 0.0;
    scalar depthMeshResolution = 0.0; // Axial cell size at the latest advancing crater front.

    // Record every newly resolved advance of the crater front.
    const scalar craterGrowthTol = SMALL;
    const scalar craterQuietTime = 1.25e-12; // Require 1.25 ps without a new resolved front cell.
    scalar rayon = 0.0;         // crater radius [m]
    scalar surfaceZ_mesh = mesh.bounds().max().z();
    reduce(surfaceZ_mesh, maxOp<scalar>()); // Use the physical surface height on every MPI rank.

    // Derive the ablation sampling bands from the surface cell size.
    scalar dz_surf = 0.0;
    scalar dr_axis = 0.0;
    {
        const volVectorField& centres = mesh.C();

        const vector pointImpact(0.0, 0.0, surfaceZ_mesh - 1e-18);
        const label impactCell = mesh.findCell(pointImpact);

        if (impactCell != -1)
        {
            const vector& cc = centres[impactCell];
            const scalar r_c = Foam::sqrt(Foam::sqr(cc.x()) + Foam::sqr(cc.y()));
            dz_surf = 1.5 * 2.0 * (surfaceZ_mesh - cc.z());   // Apply a 1.5 safety factor.
            dr_axis = 1.5 * 2.0 * r_c;
        }

        reduce(dz_surf, maxOp<scalar>());
        reduce(dr_axis, maxOp<scalar>());

        if (dz_surf < SMALL || dr_axis < SMALL)
        {
            dz_surf = 1.5 * 5.0e-9;     // Use fallback bands when no rank finds the impact cell.
            dr_axis = 1.5 * 25.0e-9;
            Info<< ">> ATTENTION : maille d'impact introuvable, tolerances par"
                   " defaut utilisees." << endl;
        }
    }

    // A one-index radial advance changes the reported diameter by twice the
    // radial cell width. This is the meaningful diameter resolution.
    const scalar diameterMeshResolution = 2.0*dr_axis/1.5;
    Info<< ">> Diagnostic ablation : bande de surface dz = " << dz_surf*1e9
        << " nm, bande d'axe dr = " << dr_axis*1e9 << " nm" << endl;

    // Read the MMS verification switch from controlDict.
    Switch mmsOn
    (
        runTime.controlDict().lookupOrDefault<Switch>("mmsVerification", false)
    );

    // Declare the MMS constants shared by the source and error calculations.
    #include "mmsParams.H"

    // Initialize both temperatures with the manufactured profile in MMS mode.
    if (mmsOn)
    {
        forAll(mesh.C(), i)
        {
            scalar z   = mesh.C()[i].z() - z0_mms;
            scalar shp = Foam::cos(kz_mms*z);          // Use the wave number defined in mmsParams.H.
            Te[i] = T0_mms + Ae_mms*shp;               // The exponential decay equals one at t = 0.
            Tl[i] = T0_mms + Al_mms*shp;
        }
        Te.correctBoundaryConditions();
        Tl.correctBoundaryConditions();
        Info << ">> MMS ACTIVE " << endl;
    }

    scalar peakTe = 300, peakTl = 300, peakKe = 0;

    Info << "\n>>Début de la simulation..." << endl;
    while (runTime.loop())
    {
        // Use adaptive time steps for physics and the controlDict value for MMS.
        if (!mmsOn)
        {
            scalar t_actuel = runTime.value();
            if (t_actuel <= 2.0e-12) runTime.setDeltaT(1.0e-15);  // Resolve the laser pulse.
            else                     runTime.setDeltaT(10.0e-15); // Advance faster after the pulse.
        }

        if (runTime.timeIndex() % 500 == 0)
            Info << "Time = " << runTime.timeName() << " s" << endl;

        // Perform the nonlinear Picard iterations.
        for (label iter = 0; iter < maxIter; iter++)
        {
            volScalarField Te_prev = Te;
            volScalarField Tl_prev = Tl;

            if (mmsOn)
            {
                // Freeze all material properties at the MMS test values.
                Ce_var = dimensionedScalar("Ce", Ce_var.dimensions(), 3.0e4);
                Cl_var = dimensionedScalar("Cl", Cl_var.dimensions(), 2.42e6);
                ke_var = dimensionedScalar("ke", ke_var.dimensions(), 235.0);
                kl_var = dimensionedScalar("kl", kl_var.dimensions(), 1.0e8);
                G_var  = dimensionedScalar("G",  G_var.dimensions(),  0.0);

                #include "mmS.H"   // Build the manufactured sources at the current time.
            }
            else
            {
                // Update temperature-dependent physical properties in every cell.
                forAll(Te, celli)
                {
                    scalar te_loc = Te[celli];
                    scalar tl_loc = Tl[celli];

                    // Interpolate the property tables with a clamped index.
                    scalar te_safe = max(300.0, min(te_loc, T_max - dT - 1.0));
                    int idx = floor(te_safe / dT);
                    scalar reste = (te_safe - idx * dT) / dT;

                    Ce_var[celli] = Ce_table[idx]*(1.0 - reste) + Ce_table[idx + 1]*reste;
                    G_var[celli]  = G_table[idx] *(1.0 - reste) + G_table[idx + 1] *reste;

                    // Evaluate the solid or liquid Shomate lattice heat capacity.
                    Cl_var[celli] = calculateShomateCl(
                        tl_loc,
                        shA_sol_val, shB_sol_val, shC_sol_val, shD_sol_val, shE_sol_val,
                        shA_liq_val, shB_liq_val, shC_liq_val, shD_liq_val, shE_liq_val,
                        cl_mult_val
                    );

                    // Evaluate the relaxation time and thermal conductivities.
                    scalar tau_e_loc = calculateTauE(te_loc, tl_loc, Ae_val, Bl_val);
                    ke_var[celli] = calculateKe(tau_e_loc, Ce_var[celli], vF_val);
                    kl_var[celli] = calculateKl(ke_var[celli]);

                }
            }

            // Update the boundary values of all property fields.
            Ce_var.correctBoundaryConditions();
            Cl_var.correctBoundaryConditions();
            ke_var.correctBoundaryConditions();
            kl_var.correctBoundaryConditions();
            G_var.correctBoundaryConditions();

            // Build the laser source only in physical mode.
            if (!mmsOn)
            {
                #include "calculateLaser.H"
            }

            // Select manufactured sources in MMS mode and the laser source otherwise.
            volScalarField Se(mmsOn ? Qe_mms : Qlaser);
            volScalarField Sl("Sl", Ql_mms);
            if (!mmsOn) { Sl = 0.0 * Sl; }   // The laser does not heat the lattice directly.

            // Solve the electron energy equation.
            fvScalarMatrix TeEqn
            (
                  Ce_var * fvm::ddt(Te)
                - fvm::laplacian(ke_var, Te)
                + fvm::Sp(G_var, Te)
                ==
                  Se + G_var * Tl
            );
            SolverPerformance<scalar> pTe = TeEqn.solve();

            // Solve the lattice energy equation.
            fvScalarMatrix TlEqn
            (
                  Cl_var * fvm::ddt(Tl)
                - fvm::laplacian(kl_var, Tl)
                + fvm::Sp(G_var, Tl)
                ==
                  Sl + G_var * Te
            );
            SolverPerformance<scalar> pTl = TlEqn.solve();

            // Compute absolute and relative Picard increments from the previous iterate.
            const scalar absErrTe = gMax(mag(Te.primitiveField() - Te_prev.primitiveField()));
            const scalar absErrTl = gMax(mag(Tl.primitiveField() - Tl_prev.primitiveField()));

            const scalar normTe = max(gMax(mag(Te.primitiveField())), SMALL);
            const scalar normTl = max(gMax(mag(Tl.primitiveField())), SMALL);

            const scalar relErrTe = absErrTe/normTe;
            const scalar relErrTl = absErrTl/normTl;

            // Retain linear-solver residuals for diagnostics only.
            const scalar resTe = pTe.initialResidual();
            const scalar resTl = pTl.initialResidual();

            if (runTime.timeIndex() % 500 == 0)
                Info<< "  Picard " << iter
                    << " | absTe=" << absErrTe << " K  absTl=" << absErrTl << " K"
                    << " | relTe=" << relErrTe << "  relTl=" << relErrTl
                    << " | linResTe=" << resTe << " linResTl=" << resTl << endl;

            // Advance only when both fields satisfy an absolute or relative tolerance.
            const bool convTe = (absErrTe < picardAbsTol) || (relErrTe < picardRelTol);
            const bool convTl = (absErrTl < picardAbsTol) || (relErrTl < picardRelTol);

            if (convTe && convTl)
            {
                break;                                  // Picard iteration has converged.
            }
            else if (iter == maxIter - 1)
            {
                Info<< "  WARNING: Picard not converged @ t=" << runTime.value()
                    << " s | absTe=" << absErrTe << " K  absTl=" << absErrTl << " K"
                    << " | relTe=" << relErrTe << "  relTl=" << relErrTl << endl;
            }

        } // End the Picard loop.
        peakTe = max(peakTe, gMax(Te.primitiveField()));
        peakTl = max(peakTl, gMax(Tl.primitiveField()));
        peakKe = max(peakKe, gMax(ke_var.primitiveField()));

        // Evaluate the phase-explosion crater where Tl is at least 0.9 Tc.
        if (!mmsOn)
        {
            const volVectorField& centres = mesh.C();

            // Reuse the globally reduced surface and axis bands computed at startup.

            forAll(mesh.cells(), cellI)
            {
                if (Tl[cellI] >= 0.9 * T_ablation)
                {
                    scalar z_i = centres[cellI].z();
                    scalar r_i = Foam::sqrt(Foam::sqr(centres[cellI].x()) + Foam::sqr(centres[cellI].y()));
                    scalar current_depth = surfaceZ_mesh - z_i;

                    // Measure depth along the symmetry axis.
                    if (r_i < dr_axis && current_depth > profondeur) profondeur = current_depth;
                    // Measure radius inside the surface band.
                    if (current_depth < dz_surf && r_i > rayon)      rayon = r_i;
                }
            }

            // Track only the last significant crater growth after global reduction.
            {
                scalar profGlobal = profondeur;
                scalar rayGlobal  = rayon;
                reduce(profGlobal, maxOp<scalar>());
                reduce(rayGlobal,  maxOp<scalar>());

                if (profGlobal > profPrevGlobal + craterGrowthTol)
                {
                    profPrevGlobal = profGlobal;
                    tLastDepthGrowth = runTime.value();

                    // Identify the cell size at the newly reached axial front.
                    // The global maximum can belong to any MPI rank, so each
                    // rank contributes only when it owns the matching cell.
                    scalar localFrontResolution = 0.0;
                    forAll(mesh.cells(), cellI)
                    {
                        if (Tl[cellI] >= 0.9*T_ablation)
                        {
                            const scalar r_i = Foam::sqrt
                            (
                                Foam::sqr(centres[cellI].x())
                              + Foam::sqr(centres[cellI].y())
                            );
                            const scalar depth_i = surfaceZ_mesh - centres[cellI].z();
                            if
                            (
                                r_i < dr_axis
                             && Foam::mag(depth_i - profGlobal) <= 0.5*dz_field[cellI]
                            )
                            {
                                localFrontResolution = max
                                (
                                    localFrontResolution,
                                    dz_field[cellI]
                                );
                            }
                        }
                    }
                    reduce(localFrontResolution, maxOp<scalar>());
                    if (localFrontResolution > SMALL)
                    {
                        depthMeshResolution = localFrontResolution;
                    }
                }
                if (rayGlobal > rayPrevGlobal + craterGrowthTol)
                {
                    rayPrevGlobal = rayGlobal;
                    tLastRadiusGrowth = runTime.value();
                }
            }
        }

        // Evaluate the MMS L2 error when verification is enabled.
        if (mmsOn)
        {
            #include "mmSL2.H"
        }

        runTime.write();

    } // End the time loop.

    // Report peak values and the diffusion length based on the reference coupling.
    {
        scalar Ldiff = (G_ref_val > VSMALL) ? Foam::sqrt(peakKe / G_ref_val)*1e9 : 0.0;   // [nm]
        Info<< "peak Te=" << peakTe << " K  peak Tl=" << peakTl
            << " K  peak ke=" << peakKe << " W/m/K"
            << "  L=sqrt(ke/G)=" << Ldiff << " nm" << endl;
    }

    Info<< "Fin du calcul.\n" << endl;
    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s\n"
        << "ClockTime = " << runTime.elapsedClockTime() << " s" << endl;

    // Report the final ablation dimensions in physical mode.
    if (!mmsOn)
    {
        scalar F_Jcm2 = F_laser.value() / 10000.0;   // Convert J/m2 to J/cm2.

        // Evaluate the theoretical logarithmic depth law.
        scalar delta_fit = 98.25;
        scalar F_th_depth = 0.49;
        scalar depth_th = (F_Jcm2 > F_th_depth)
                        ? delta_fit * Foam::log(F_Jcm2 / F_th_depth) : 0.0;

        // Evaluate the theoretical Gaussian-beam diameter law.
        scalar w0_fit = 9.69;
        scalar F_th_diam = 0.53;
        scalar diametre_th = (F_Jcm2 > F_th_diam)
                           ? Foam::sqrt(2.0 * Foam::sqr(w0_fit) * Foam::log(F_Jcm2 / F_th_diam)) : 0.0;

        // Reduce the crater dimensions across all MPI ranks.
        reduce(profondeur, maxOp<scalar>());
        reduce(rayon, maxOp<scalar>());
        scalar diametre = 2.0 * rayon;

        // Check whether the crater remained unchanged for the required quiet time.
        const scalar tEnd = runTime.value();
        const scalar tLastGrowth = max(tLastDepthGrowth, tLastRadiusGrowth);
        const scalar quietTime = max(scalar(0.0), tEnd - tLastGrowth);
        const scalar reportedDepthResolution = max(depthMeshResolution, dz_surf/1.5);
        Info << nl << " CONVERGENCE DU CRATERE" << endl;
        Info << "Crater-growth event                : new resolved front cell" << endl;
        Info << "Duree minimale sans croissance      : " << craterQuietTime*1e12
             << " ps" << endl;
        Info << "Derniere croissance en profondeur : " << tLastDepthGrowth*1e12
             << " ps   (fin du calcul : " << tEnd*1e12 << " ps)" << endl;
        Info << "Derniere croissance en rayon      : " << tLastRadiusGrowth*1e12
             << " ps" << endl;
        Info << "Duree observee sans croissance    : " << quietTime*1e12
             << " ps" << endl;
        Info << "Depth mesh resolution             : " << reportedDepthResolution*1e9
             << " nm" << endl;
        Info << "Diameter mesh resolution          : " << diameterMeshResolution*1e6
             << " um" << endl;
        if (quietTime < craterQuietTime)
        {
            Info << "CRATER STATUS: NOT CONVERGED." << nl
                 << "    No new resolved front cell was observed for only " << quietTime*1e12
                 << " ps; the required quiet time is " << craterQuietTime*1e12 << " ps." << nl
                 << "    Increase endTime or run ConvergedRun.sh." << endl;
        }
        else
        {
            Info << "CRATER STATUS: CONVERGED." << nl
                 << "    No new resolved front cell was observed for " << quietTime*1e12
                 << " ps (required: " << craterQuietTime*1e12 << " ps)." << endl;
        }

        Info << nl << " BILAN ABLATION" << endl;
        Info << "Fluence laser simulee         : " << F_Jcm2 << " J/cm2" << endl;
        Info << "------------------------------------------------------" << endl;
        Info << "Profondeur THEORIQUE : " << depth_th << " nm" << endl;
        Info << "Profondeur NUMERIQUE : " << profondeur * 1e9 << " nm" << endl;
        Info << "------------------------------------------------------" << endl;
        Info << "Diametre THEORIQUE  : " << diametre_th << " um" << endl;
        Info << "Diametre NUMERIQUE  : " << diametre * 1e6 << " um" << endl;
        Info << "------------------------------------------------------" << nl << endl;

        if (Pstream::master())
        {
            std::ofstream outFile("Res.txt", std::ios_base::app);
            if (outFile.is_open())
            {
                // Write fluence, theoretical and numerical depth, and theoretical and numerical diameter.
                outFile << F_Jcm2 << "   "
                        << depth_th << "   "
                        << (profondeur * 1e9) << "   "
                        << diametre_th << "   "
                        << (diametre * 1e6) << "\n";
                outFile.close();
            }
            else
            {
                Info << "Attention : Impossible d'ouvrir le fichier Res.txt" << endl;
            }
        }
    }

    Info<< "End\n" << endl;
    return 0;
}
