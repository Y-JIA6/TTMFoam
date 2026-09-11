// ===========================================================================
//  updateOptical.H
//  Temperature-dependent Drude optical response of aluminum.
//  Returns the reflectivity R and optical penetration depth delta_opt.
//
//  Calibration: the raw Drude reflectivity is scaled at the reference
//  relaxation time tau_e_ref to match R = 0.91 at 1032 nm. The same
//  fluence-independent factor is then applied at every time step, preserving
//  the Drude trend while anchoring the ambient value to experiment.
//
//  Ballistic transport is not an optical property. calculateLaser.H forms
//  delta_eff = delta_opt + vF*tau_e locally in every cell.
//
//  SI units: tau_e, tau_e_ref [s]; lambda, delta_opt [m].
// ===========================================================================
#ifndef UPDATEOPTICAL_H
#define UPDATEOPTICAL_H

#include "fvCFD.H"

struct OpticalProps
{
    scalar R;           // Surface reflectivity.
    scalar delta_opt;   // Optical penetration depth in m.
};

// ---------------------------------------------------------------------------
// Evaluate the temperature-dependent Drude optical properties
// ---------------------------------------------------------------------------
inline OpticalProps calculateOpticalProperties
(
    scalar tau_e,
    scalar lambda,
    scalar tau_e_ref
)
{
    OpticalProps props;

    // -----------------------------------------------------------------------
    // 1. Drude constants and ambient-reflectivity calibration
    // -----------------------------------------------------------------------
    const scalar c        = 299792458.0;        // Speed of light in m/s.
    const scalar eps_inf  = 8.5116;             // Aluminum high-frequency permittivity.
    const scalar e_charge = 1.602176634e-19;    // Elementary charge in C.
    const scalar m_e      = 9.1093837e-31;      // Electron mass in kg.
    const scalar eps_0    = 8.8541878e-12;      // Vacuum permittivity in F/m.
    const scalar n_e      = 1.808e29;           // Aluminum valence-electron density in 1/m3.

    // Set the measured ambient aluminum reflectivity at 1032 nm.
    const scalar R_dict = 0.91;

    // Compute the plasma and laser angular frequencies.
    scalar omega_p = Foam::sqrt((n_e * Foam::pow(e_charge, 2)) / (m_e * eps_0));
    scalar omega   = (2.0 * constant::mathematical::pi * c) / lambda;

    // Evaluate the Drude reflectivity at the reference relaxation time.
    const scalar gamma_init = 1.0 / tau_e_ref;
    const scalar denom_init = Foam::pow(omega, 2) + Foam::pow(gamma_init, 2);
    const scalar eps1_init = eps_inf - Foam::pow(omega_p, 2)/denom_init;
    const scalar eps2_init = Foam::pow(omega_p, 2)*gamma_init/(omega*denom_init);
    const scalar mod_eps_init = Foam::sqrt(Foam::sqr(eps1_init) + Foam::sqr(eps2_init));
    const scalar n_init = Foam::sqrt((eps1_init + mod_eps_init)/2.0);
    const scalar k_init = Foam::sqrt((-eps1_init + mod_eps_init)/2.0);
    const scalar R_drude_init = (Foam::sqr(n_init - 1.0) + Foam::sqr(k_init))
                              / (Foam::sqr(n_init + 1.0) + Foam::sqr(k_init));
    const scalar scaleFactor = R_dict / R_drude_init;

    // -----------------------------------------------------------------------
    // 2. Temperature-dependent Drude dielectric response
    // -----------------------------------------------------------------------
    scalar gamma_D = 1.0 / tau_e;
    scalar denom = Foam::pow(omega, 2) + Foam::pow(gamma_D, 2);

    scalar eps1 = eps_inf - (Foam::pow(omega_p, 2) / denom);
    scalar eps2 = (Foam::pow(omega_p, 2) * gamma_D) / (omega * denom);
    scalar mod_eps = Foam::sqrt(Foam::pow(eps1, 2) + Foam::pow(eps2, 2));

    scalar n = Foam::sqrt((eps1 + mod_eps) / 2.0);
    scalar k = Foam::sqrt((-eps1 + mod_eps) / 2.0);

    const scalar R_drude = (Foam::sqr(n - 1.0) + Foam::sqr(k))
                         / (Foam::sqr(n + 1.0) + Foam::sqr(k));
    props.R = min(R_drude * scaleFactor, R_dict);

    // -----------------------------------------------------------------------
    // 3. Optical penetration depth from the extinction coefficient
    // -----------------------------------------------------------------------
    scalar alpha     = (2.0 * omega * k) / c;
    scalar delta_opt = 1.0 / alpha;

    props.delta_opt = delta_opt;

    return props;
}

#endif
