#ifndef UPDATEOPTICAL_H
#define UPDATEOPTICAL_H

#include "fvCFD.H"

// ===========================================================================
//  updateOptical.H
//  Temperature-dependent optical response of aluminum via the Drude model.
//  Returns: reflectivity R, optical skin depth delta_opt, and effective
//  penetration depth delta_eff = delta_opt + delta_ball.
//
//  Calibration: the raw Drude reflectivity is scaled so that, at the reference
//  electron relaxation time (tau_e_init, room-temperature value), the model
//  matches the experimental R_target = 0.91. The same scaleFactor is then
//  applied at all temperatures, preserving the Drude trend while anchoring the
//  absolute level to experiment.
//
//  SI units: tau_e [s], lambda [m], vF [m/s], Te_local [K].
// ===========================================================================

struct OpticalProps
{
    scalar R;           // reflectivity
    scalar delta_eff;   // effective penetration depth (optical + ballistic) [m]
    scalar delta_opt;   // pure optical penetration depth [m]
};

inline OpticalProps calculateOpticalProperties(scalar tau_e, scalar lambda, scalar vF, scalar Te_local)
{
    OpticalProps props;

    // --- Physical constants ---
    const scalar c        = 299792458.0;        // speed of light [m/s]
    const scalar eps_inf  = 8.5116;             // Al high-frequency permittivity [-]
    const scalar e_charge = 1.602176634e-19;    // elementary charge [C]
    const scalar m_e      = 9.1093837e-31;      // electron mass [kg]
    const scalar eps_0    = 8.8541878e-12;      // vacuum permittivity [F/m]
    const scalar n_e      = 1.808e29;           // Al valence electron density [m^-3]

    // Experimental reflectivity for room temperature
    const scalar R_dict = 0.91;   // 0.91 for Al at 1032 nm

    // Plasma and laser angular frequencies
    scalar omega_p = Foam::sqrt((n_e * Foam::pow(e_charge, 2)) / (m_e * eps_0));
    scalar omega   = (2.0 * constant::mathematical::pi * c) / lambda;

    // -----------------------------------------------------------------------
    // Step 1: Drude reflectivity at the reference (room-temperature) state.
    //         Used only to compute the factor scaleFactor.
    // -----------------------------------------------------------------------
    const scalar tau_e_init = 5.626667e-15;     // reference relaxation time [s]
    scalar gamma_init = 1.0 / tau_e_init;       // reference collision frequency [1/s]
    scalar denom_init = Foam::pow(omega, 2) + Foam::pow(gamma_init, 2);

    scalar eps1_init = eps_inf - (Foam::pow(omega_p, 2) / denom_init);
    scalar eps2_init = (Foam::pow(omega_p, 2) * gamma_init) / (omega * denom_init);
    scalar mod_eps_init = Foam::sqrt(Foam::pow(eps1_init, 2) + Foam::pow(eps2_init, 2));

    scalar n_init = Foam::sqrt((eps1_init + mod_eps_init) / 2.0);
    scalar k_init = Foam::sqrt((-eps1_init + mod_eps_init) / 2.0);

    scalar R_drude_init = (Foam::pow(n_init - 1.0, 2) + Foam::pow(k_init, 2))
                        / (Foam::pow(n_init + 1.0, 2) + Foam::pow(k_init, 2));

    // Multiplicative calibration factor (anchors level, keeps Drude trend)
    const scalar scaleFactor = R_dict / R_drude_init;

    // -----------------------------------------------------------------------
    // Step 2: Temperature-dependent Drude reflectivity (uses current tau_e).
    // -----------------------------------------------------------------------
    scalar gamma_D = 1.0 / tau_e;
    scalar denom = Foam::pow(omega, 2) + Foam::pow(gamma_D, 2);

    scalar eps1 = eps_inf - (Foam::pow(omega_p, 2) / denom);
    scalar eps2 = (Foam::pow(omega_p, 2) * gamma_D) / (omega * denom);
    scalar mod_eps = Foam::sqrt(Foam::pow(eps1, 2) + Foam::pow(eps2, 2));

    scalar n = Foam::sqrt((eps1 + mod_eps) / 2.0);
    scalar k = Foam::sqrt((-eps1 + mod_eps) / 2.0);

    scalar R_drude = (Foam::pow(n - 1.0, 2) + Foam::pow(k, 2))
                   / (Foam::pow(n + 1.0, 2) + Foam::pow(k, 2));

    // -----------------------------------------------------------------------
    // Step 3: Apply calibration factor to the Drude reflectivity.
    // -----------------------------------------------------------------------
    props.R = R_drude * scaleFactor;
    if (props.R > R_dict) props.R = R_dict;

    // -----------------------------------------------------------------------
    // Step 4: Penetration depths.
    // -----------------------------------------------------------------------
    scalar alpha     = (2.0 * omega * k) / c;
    scalar delta_opt = 1.0 / alpha;

    scalar delta_ball = vF * tau_e;

    props.delta_opt = delta_opt;
    props.delta_eff = delta_opt + delta_ball;

    return props;
}

#endif