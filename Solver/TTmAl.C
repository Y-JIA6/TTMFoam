// ----------------------------------------------------------------------------
//  TTmAl.C
//  Two-temperature model (TTM) solver for ultrashort-pulse laser ablation of
//  aluminum, on a 2D-axisymmetric OpenFOAM mesh.
//
//  Solves the coupled electron/lattice energy equations with
//  temperature-dependent properties (Fermi-Dirac Ce, G; Drude optics; Shomate
//  Cl). A Picard loop linearizes the nonlinear coupling each time step.
//
//  Two run modes, selected by 'mmsVerification' in controlDict:
//    false : physics (laser source, ablation diagnostics)
//    true  : Method of Manufactured Solutions verification (no physics)
// ----------------------------------------------------------------------------

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

    // Cache dictionary scalars once (avoids repeated .value() in the cell loop)
    scalar gamma_val = gamma.value();
    scalar Ae_val = Ae.value();
    scalar Bl_val = Bl.value();
    scalar vF_val = vF.value();

    scalar ke_base = ke.value();
    scalar kl_base = kl.value();
    scalar Cl_base = Cl.value();

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

    #include "createTables.H"   // builds Ce_table, G_table, mu_table (once)

    const int maxIter = 15;     // Picard iteration cap per time step

    // Ablation diagnostic parameters
    scalar T_ablation = 6700.0; // critical temperature Tc [K]
    scalar profondeur = 0.0;    // crater depth  [m]
    scalar rayon = 0.0;         // crater radius [m]
    scalar surfaceZ_mesh = mesh.bounds().max().z();

    // --- MMS verification switch (read from controlDict) ---
    Switch mmsOn
    (
        runTime.controlDict().lookupOrDefault<Switch>("mmsVerification", false)
    );

    // MMS constants, shared by mmS.H and mmSL2.H (declared once here)
    #include "mmsParams.H"

    // In MMS mode, set Te/Tl to the manufactured profile at t=0
    if (mmsOn)
    {
        forAll(mesh.C(), i)
        {
            scalar z   = mesh.C()[i].z() - z0_mms;
            scalar shp = Foam::cos(kz_mms*z);          // kz_mms from mmsParams.H
            Te[i] = T0_mms + Ae_mms*shp;               // t=0 -> decay = 1
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
        // Adaptive time step (physics only). MMS uses the fixed controlDict dt
        if (!mmsOn)
        {
            scalar t_actuel = runTime.value();
            if (t_actuel <= 2.0e-12) runTime.setDeltaT(1.0e-15);  // during pulse
            else                     runTime.setDeltaT(10.0e-15); // after pulse
        }

        if (runTime.timeIndex() % 500 == 0)
            Info << "Time = " << runTime.timeName() << " s" << endl;

        // ------- Picard loop -------------------
        for (int iter = 0; iter < maxIter; iter++)
        {
            volScalarField Te_prev = Te;
            volScalarField Tl_prev = Tl;

            if (mmsOn)
            {
                // freeze properties to the TEST constants.
                Ce_var = dimensionedScalar("Ce", Ce_var.dimensions(), 3.0e4);
                Cl_var = dimensionedScalar("Cl", Cl_var.dimensions(), 2.42e6);
                ke_var = dimensionedScalar("ke", ke_var.dimensions(), 235.0);
                kl_var = dimensionedScalar("kl", kl_var.dimensions(), 1.0e8);
                G_var  = dimensionedScalar("G",  G_var.dimensions(),  0.0);

                #include "mmS.H"   // builds Qe_mms, Ql_mms at current time
            }
            else
            {
                // PHYSICS: temperature-dependent properties, per cell.
                forAll(Te, celli)
                {
                    scalar te_loc = Te[celli];
                    scalar tl_loc = Tl[celli];

                    // Table lookup with linear interpolation (clamped index)
                    scalar te_safe = max(300.0, min(te_loc, T_max - dT - 1.0));
                    int idx = floor(te_safe / dT);
                    scalar reste = (te_safe - idx * dT) / dT;

                    Ce_var[celli] = Ce_table[idx]*(1.0 - reste) + Ce_table[idx + 1]*reste;
                    G_var[celli]  = G_table[idx] *(1.0 - reste) + G_table[idx + 1] *reste;

                    // Lattice heat capacity (Shomate, solid/liquid)
                    Cl_var[celli] = calculateShomateCl(
                        tl_loc,
                        shA_sol_val, shB_sol_val, shC_sol_val, shD_sol_val, shE_sol_val,
                        shA_liq_val, shB_liq_val, shC_liq_val, shD_liq_val, shE_liq_val,
                        cl_mult_val, Cl_base
                    );

                    // Relaxation time, conductivities
                    scalar tau_e_loc = calculateTauE(te_loc, tl_loc, Ae_val, Bl_val);
                    ke_var[celli] = calculateKe(tau_e_loc, Ce_var[celli], vF_val, ke_base);
                    kl_var[celli] = calculateKl(ke_var[celli], kl_base);

                }
            }

            // Boundary conditions on all property fields
            Ce_var.correctBoundaryConditions();
            Cl_var.correctBoundaryConditions();
            ke_var.correctBoundaryConditions();
            kl_var.correctBoundaryConditions();
            G_var.correctBoundaryConditions();

            // Laser source (physics only)
            if (!mmsOn)
            {
                #include "calculateLaser.H"
            }

            // Source selection: MMS -> Qe_mms/Ql_mms ; physics -> Qlaser/0
            volScalarField Se(mmsOn ? Qe_mms : Qlaser);
            volScalarField Sl("Sl", Ql_mms);
            if (!mmsOn) { Sl = 0.0 * Sl; }   // lattice has no laser source

            // --- Electron energy equation ---
            fvScalarMatrix TeEqn
            (
                  Ce_var * fvm::ddt(Te)
                - fvm::laplacian(ke_var, Te)
                + fvm::Sp(G_var, Te)
                ==
                  Se + G_var * Tl
            );
            SolverPerformance<scalar> pTe = TeEqn.solve();

            // --- Lattice energy equation ---
            fvScalarMatrix TlEqn
            (
                  Cl_var * fvm::ddt(Tl)
                - fvm::laplacian(kl_var, Tl)
                + fvm::Sp(G_var, Tl)
                ==
                  Sl + G_var * Te
            );
            SolverPerformance<scalar> pTl = TlEqn.solve();

            // Convergence monitor (field change + linear residual)
            scalar deltaTe = gMax(mag(Te.primitiveField() - Te_prev.primitiveField()));
            scalar deltaTl = gMax(mag(Tl.primitiveField() - Tl_prev.primitiveField()));
            scalar resTe = pTe.initialResidual();
            scalar resTl = pTl.initialResidual();

            if (runTime.timeIndex() % 500 == 0)
                Info<< "  Picard " << iter
                    << " | dTe=" << deltaTe << " dTl=" << deltaTl
                    << " | resTe=" << resTe << " resTl=" << resTl << endl;

            const scalar outerResTol = 1e-4;
            if (resTe < outerResTol && resTl < outerResTol)
            {
                break;                                  // Picard converged
            }
            else if (iter == maxIter - 1)
            {
                Info<< "  ATTENTION: Picard non-converge @ t=" << runTime.value()
                    << " | dTe=" << deltaTe << " dTl=" << deltaTl
                    << " | resTe=" << resTe << " resTl=" << resTl << endl;
            }

        } // end Picard

        peakTe = max(peakTe, gMax(Te.primitiveField()));
        peakTl = max(peakTl, gMax(Tl.primitiveField()));
        peakKe = max(peakKe, gMax(ke_var.primitiveField()));

        // ------- Ablation diagnostics (physics only) -------------------
        // Crater = region where Tl >= 0.9*Tc (phase-explosion criterion).
        if (!mmsOn)
        {
            const volVectorField& centres = mesh.C();

            // Tolerances sized from the surface cell at the beam center
            scalar dz_surf = 0.0;
            scalar dr_axis = 0.0;

            vector pointImpact(0.0, 0.0, surfaceZ_mesh - 1e-12);
            label impactCell = mesh.findCell(pointImpact);

            if (impactCell != -1)
            {
                scalar z_c = centres[impactCell].z();
                scalar r_c = Foam::sqrt(Foam::sqr(centres[impactCell].x()) + Foam::sqr(centres[impactCell].y()));
                scalar dz_cell = 2.0 * (surfaceZ_mesh - z_c);
                scalar dr_cell = 2.0 * r_c;
                dz_surf = 1.5 * dz_cell;     // 1.5x safety band
                dr_axis = 1.5 * dr_cell;
            }
            else
            {
                dz_surf = 1.5 * 5.0e-9;      // fallback if center cell not found
                dr_axis = 1.5 * 25.0e-9;
            }

            forAll(mesh.cells(), cellI)
            {
                if (Tl[cellI] >= 0.9 * T_ablation)
                {
                    scalar z_i = centres[cellI].z();
                    scalar r_i = Foam::sqrt(Foam::sqr(centres[cellI].x()) + Foam::sqr(centres[cellI].y()));
                    scalar current_depth = surfaceZ_mesh - z_i;

                    // depth: max along the axis (r < dr_axis)
                    if (r_i < dr_axis && current_depth > profondeur) profondeur = current_depth;
                    // radius: max at the surface (depth < dz_surf)
                    if (current_depth < dz_surf && r_i > rayon)      rayon = r_i;
                }
            }
        }

        // ------- L2 error monitor (MMS only) ---------------------------
        if (mmsOn)
        {
            #include "mmSL2.H"
        }

        runTime.write();

    } // end time loop

    // Final peak summary. Diffusion length uses the dictionary G (guarded).
    {
        scalar Ldiff = (G_ref_val > VSMALL) ? Foam::sqrt(peakKe / G_ref_val)*1e9 : 0.0;   // [nm]
        Info<< "peak Te=" << peakTe << " K  peak Tl=" << peakTl
            << " K  peak ke=" << peakKe << " W/m/K"
            << "  L=sqrt(ke/G)=" << Ldiff << " nm" << endl;
    }

    Info<< "Fin du calcul.\n" << endl;
    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s\n"
        << "ClockTime = " << runTime.elapsedClockTime() << " s" << endl;

    // --------- Final ablation report (physics only) -----------------------
    if (!mmsOn)
    {
        scalar F_Jcm2 = F_laser.value() / 10000.0;   // J/m^2 -> J/cm^2

        // Theoretical depth law: h = delta * ln(F/Fth)
        scalar delta_fit = 98.25;
        scalar F_th_depth = 0.49;
        scalar depth_th = (F_Jcm2 > F_th_depth)
                        ? delta_fit * Foam::log(F_Jcm2 / F_th_depth) : 0.0;

        // Theoretical diameter law: D = sqrt(2 w0^2 ln(F/Fth))
        scalar w0_fit = 9.69;
        scalar F_th_diam = 0.53;
        scalar diametre_th = (F_Jcm2 > F_th_diam)
                           ? Foam::sqrt(2.0 * Foam::sqr(w0_fit) * Foam::log(F_Jcm2 / F_th_diam)) : 0.0;

        // Parallel reduction of the crater extents
        reduce(profondeur, maxOp<scalar>());
        reduce(rayon, maxOp<scalar>());
        scalar diametre = 2.0 * rayon;

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
                // Fluence | depth_theory | depth_num(nm) | diam_theory | diam_num(um)
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