// ===========================================================================
//  createTables.H
//  Builds the FEG lookup tables ONCE at startup: mu(Te), Ce(Te), G(Te).
//  Included AFTER createFields.H and updateThermo.H.
//
//  Inputs:
//     gamma_val : scalar, Sommerfeld coefficient [J/m^3/K^2]
//     G, Ce     : dimensionedScalar, 300 K reference values
//  Uses the same functions/argument order as updateThermo.H:
//     calculateMu(Te,EF,gamma); calculateCe(Te,EF,gamma,mu);
//     calculateIntegralG(Te,EF,gamma,mu) -> shape integral I(Te)
//
//  The solver then interpolates these tables per cell each step, which is far
//  cheaper than re-integrating Fermi-Dirac in every cell.
// ===========================================================================

Info<< "\n>> Generation des tables FEG : mu(Te), Ce(Te), G(Te)" << endl;

// --- Temperature discretization --------------------------------------------
const scalar T_max = 2.0e5;                 // K (covers Te up to plasma range)
const scalar dT    = 10.0;                  // K resolution
const label  num_entries = label(T_max/dT) + 10;

// --- FEG parameters (must match updateThermo.H and the paper) --------------
const scalar EF_val = 11.7*1.602176634e-19; // J  

// --- 300 K references / floors (from ttmProperties) ------------------------
const scalar G_ref_val = G.value();         // W/m^3/K
const scalar Ce_min    = Ce.value();        // J/m^3/K 

// --- Tables -----------------------------------------------------------------
List<scalar> mu_table(num_entries, 0.0);
List<scalar> Ce_table(num_entries, 0.0);
List<scalar> G_table (num_entries, 0.0);

// --- Reference coupling integral at 300 K (denominator of the G ratio) ------
const scalar mu_ref  = calculateMu(300.0, EF_val, gamma_val);
const scalar I_G_ref = calculateIntegralG(300.0, EF_val, gamma_val, mu_ref);

// --- Fill loop --------------------------------------------------------------
for (label i = 0; i < num_entries; ++i)
{
    const scalar Te = max(scalar(1.0), i*dT);

    const scalar mu = calculateMu(Te, EF_val, gamma_val);
    mu_table[i] = mu;
    Ce_table[i] = max(Ce_min, calculateCe(Te, EF_val, gamma_val, mu));
    G_table[i]  = G_ref_val*calculateIntegralG(Te, EF_val, gamma_val, mu)/I_G_ref;

    if (i % 2000 == 0)
    {
        Info<< "  -- Te = " << Te << " K :"
            << "  Ce = " << Ce_table[i]
            << "  G = "  << G_table[i] << endl;
    }
}

Info<< ">> Tables FEG generees (" << num_entries
    << " entrees, pas " << dT << " K)." << endl;

// --- Verification export -------
if (Pstream::master())
{
    OFstream propFile("FEG_param.dat");
    propFile << "# Te(K)  Ce(J/m3.K)  G(W/m3.K)  mu(eV)" << nl;
    for (label i = 0; i < num_entries; i += 10)
    {
        const scalar Te = max(scalar(1.0), i*dT);
        propFile << Te << "  " << Ce_table[i] << "  " << G_table[i]
                 << "  " << mu_table[i]/1.602176634e-19 << nl;
    }
}

Info<< ">> Generation terminee en " << runTime.elapsedCpuTime() << " s\n" << endl;