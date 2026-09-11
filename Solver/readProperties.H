// ===========================================================================
//  readProperties.H
//  Reads laser and material constants from the case dictionaries at run time.
//  Parameters can therefore be changed without recompiling the solver.
//
//    - laserProperties : constant/laserProperties
//    - ttmProperties   : constant/ttmProperties
//
//  Ambient Ce, Cl, ke, and kl initialize the physical fields. Temperature-
//  dependent values are subsequently updated by the thermophysical, optical,
//  and laser-source modules.
// ===========================================================================

Info<< "Lecture des propriétés laser et matériaux\n" << endl;

// ---------------------------------------------------------------------------
// 1. Laser dictionary
// ---------------------------------------------------------------------------
IOdictionary laserProperties
(
    IOobject("laserProperties", runTime.constant(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE)
);

dimensionedScalar F_laser("F_laser", dimEnergy/dimArea,  laserProperties); // Peak fluence in J/m2.
dimensionedScalar w0     ("w0",      dimLength,          laserProperties); // Beam waist at 1/e2 intensity in m.
dimensionedScalar tau    ("tau",     dimTime,            laserProperties); // Pulse FWHM in s.
dimensionedScalar lambda ("lambda",  dimLength,          laserProperties); // Wavelength in m.

// ---------------------------------------------------------------------------
// 2. Two-temperature-model material dictionary
// ---------------------------------------------------------------------------
IOdictionary ttmProperties
(
    IOobject("ttmProperties", runTime.constant(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE)
);

// ---------------------------------------------------------------------------
// 3. Ambient property references and electron thermodynamics
// ---------------------------------------------------------------------------
dimensionedScalar T0("T0", dimTemperature, ttmProperties);                     // Temperature in K.

dimensionedScalar Ce("Ce", dimEnergy/dimVolume/dimTemperature, ttmProperties); // Electron heat capacity in J/m3/K.
dimensionedScalar Cl("Cl", dimEnergy/dimVolume/dimTemperature, ttmProperties); // Lattice heat capacity in J/m3/K.
dimensionedScalar ke("ke", dimPower/dimLength/dimTemperature,  ttmProperties); // Electron conductivity in W/m/K.
dimensionedScalar kl("kl", dimPower/dimLength/dimTemperature,  ttmProperties); // Lattice conductivity in W/m/K.
dimensionedScalar G ("G",  dimPower/dimVolume/dimTemperature,  ttmProperties); // Coupling factor in W/m3/K.

dimensionedScalar gamma("gamma", dimEnergy/dimVolume/dimTemperature/dimTemperature, ttmProperties);
dimensionedScalar mu   ("mu",    dimEnergy, ttmProperties); // Fermi energy in J.

// ---------------------------------------------------------------------------
// 4. Transport, Drude, and electron--phonon parameters
// ---------------------------------------------------------------------------
dimensionedScalar vF("vF", dimLength/dimTime, ttmProperties);                 // Fermi velocity in m/s.
dimensionedScalar Ae("Ae", dimless/dimTime/dimTemperature/dimTemperature, ttmProperties); // Electron-electron coefficient in 1/(K2 s).
dimensionedScalar Bl("Bl", dimless/dimTime/dimTemperature, ttmProperties);    // Electron-phonon coefficient in 1/(K s).

// ---------------------------------------------------------------------------
// 5. Solid and liquid Shomate heat-capacity coefficients
// ---------------------------------------------------------------------------
dimensionedScalar ShA_sol("ShA_sol", dimless, ttmProperties);
dimensionedScalar ShB_sol("ShB_sol", dimless, ttmProperties);
dimensionedScalar ShC_sol("ShC_sol", dimless, ttmProperties);
dimensionedScalar ShD_sol("ShD_sol", dimless, ttmProperties);
dimensionedScalar ShE_sol("ShE_sol", dimless, ttmProperties);

dimensionedScalar ShA_liq("ShA_liq", dimless, ttmProperties);
dimensionedScalar ShB_liq("ShB_liq", dimless, ttmProperties);
dimensionedScalar ShC_liq("ShC_liq", dimless, ttmProperties);
dimensionedScalar ShD_liq("ShD_liq", dimless, ttmProperties);
dimensionedScalar ShE_liq("ShE_liq", dimless, ttmProperties);

dimensionedScalar Cl_mult("Cl_mult", dimEnergy/dimVolume/dimTemperature, ttmProperties); // Convert Shomate values to volumetric heat capacity.
