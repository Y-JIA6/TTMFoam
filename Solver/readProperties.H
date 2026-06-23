// ============================================================================
//  readProperties.H
//  Reads all laser and material constants from the case dictionaries at run
//  time. Because everything is read here, parameters can be changed WITHOUT
//  recompiling the solver (only the case files change).
//    - laserProperties : constant/laserProperties
//    - ttmProperties   : constant/ttmProperties
//  Room-temperature optical values (R, delta_eff, delta_opt) are only SEED
//  values; they are overwritten per-cell each step by updateOptical.H.
// ============================================================================

Info<< "Lecture des propriétés laser et matériaux\n" << endl;

// --- 1. Laser dictionary -----------------------------------------------------
IOdictionary laserProperties
(
    IOobject("laserProperties", runTime.constant(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE)
);

dimensionedScalar F_laser("F_laser", dimEnergy/dimArea,  laserProperties); // peak fluence [J/m^2]
dimensionedScalar w0     ("w0",      dimLength,          laserProperties); // 1/e^2 beam waist [m]
dimensionedScalar tau    ("tau",     dimTime,            laserProperties); // pulse FWHM [s]
dimensionedScalar lambda ("lambda",  dimLength,          laserProperties); // wavelength [m]

// --- 2. Material (TTM) dictionary --------------------------------------------
IOdictionary ttmProperties
(
    IOobject("ttmProperties", runTime.constant(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE)
);

// Base properties and coupling (room-temperature reference values)
dimensionedScalar Ce("Ce", dimEnergy/dimVolume/dimTemperature, ttmProperties); // [J/m^3/K]
dimensionedScalar Cl("Cl", dimEnergy/dimVolume/dimTemperature, ttmProperties); // [J/m^3/K]
dimensionedScalar ke("ke", dimPower/dimLength/dimTemperature,  ttmProperties); // [W/m/K]
dimensionedScalar kl("kl", dimPower/dimLength/dimTemperature,  ttmProperties); // [W/m/K]
dimensionedScalar G ("G",  dimPower/dimVolume/dimTemperature,  ttmProperties); // G_0 [W/m^3/K]

// Electron thermodynamics (used by the table builder / Sommerfeld reference)
dimensionedScalar gamma("gamma", dimEnergy/dimVolume/dimTemperature/dimTemperature, ttmProperties);
dimensionedScalar mu   ("mu",    dimEnergy, ttmProperties); // Fermi energy [J]

// Drude kinetics
dimensionedScalar vF("vF", dimLength/dimTime, ttmProperties);                 // Fermi velocity [m/s]
dimensionedScalar Ae("Ae", dimless/dimTime/dimTemperature/dimTemperature, ttmProperties); // e-e coeff [1/(K^2 s)]
dimensionedScalar Bl("Bl", dimless/dimTime/dimTemperature, ttmProperties);    // e-ph coeff [1/(K s)]

// Shomate coefficients: two sets (solid + liquid) to capture the melt transition
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

dimensionedScalar Cl_mult("Cl_mult", dimEnergy/dimVolume/dimTemperature, ttmProperties); // Shomate->volumetric scale