#ifndef UPDATE_THERMO_H
#define UPDATE_THERMO_H

#include "fvCFD.H"

// ===========================================================================
//  updateThermo.H
//  Free-electron-gas (FEG) thermophysical models for the aluminum TTM.
//
//  Analytic DOS:  g(E) = C*sqrt(E),  C = (3/2)*n/EF^{3/2}
//                                      = 3*gamma/(pi^2 kB^2 sqrt(EF))
//
//  SI units:
//     Te, Tl : K        EF, mu : J  (EF_Al ~ 11.7 eV ~ 1.87e-18 J)
//     gamma  : J/m^3/K^2 (Ce(300)=2.3e4 -> gamma=76.7)
//     Ce, Cl : J/m^3/K   ke, kl : W/m/K   tau_e : s   vF : m/s (~2.03e6)
//
//  CONSISTENCY: for a FEG, gamma, EF and n are NOT independent:
//     n = 2*gamma*EF/(pi^2 kB^2).  Choose gamma and EF consistently.
// ===========================================================================

namespace ttm
{
    const scalar kB = 1.380649e-23;     // J/K
    const scalar eV = 1.602176634e-19;  // J
}

// ---------------------------------------------------------------------------
// 1. Chemical potential mu(Te) at conserved electron number.
// ---------------------------------------------------------------------------
inline scalar calculateMu(scalar Te, scalar EF, scalar gamma)
{
    const scalar kB = ttm::kB;
    const scalar pi = constant::mathematical::pi;
    const scalar kT = kB*Te;

    const scalar r   = kT/EF;
    const scalar muS = EF*(1.0 - (pi*pi/12.0)*r*r);   // Sommerfeld estimate
    if (Te <= 13000.0) return muS;

    const scalar n = 2.0*gamma*EF/(pi*pi*kB*kB);       // electron density (FEG)
    const scalar C = 1.5*n*Foam::pow(EF, -1.5);        // DOS prefactor

    // bisection: find mu such that integral of g*f equals n
    scalar lo = -15.0*kT, hi = EF, mu = EF;
    const int M = 1000;
    for (int it = 0; it < 60; ++it)
    {
        mu = 0.5*(lo + hi);
        const scalar Emax = max(mu, scalar(0.0)) + 50.0*kT;
        const scalar dE = Emax/M;
        scalar nCalc = 0.0;
        for (int i = 0; i <= M; ++i)
        {
            const scalar E = i*dE;
            const scalar X = (E - mu)/kT;
            scalar f = 0.0;
            if (X < -50.0)      f = 1.0;                   // avoid exp overflow
            else if (X < 50.0)  f = 1.0/(Foam::exp(X) + 1.0);
            const scalar w = (i == 0 || i == M) ? 0.5 : 1.0; // trapezoid ends
            nCalc += w*C*Foam::sqrt(E)*f*dE;
        }
        if (nCalc > n) hi = mu; else lo = mu;
    }

    if (Te >= 14000.0) return mu;

    // Hann blend between Sommerfeld (muS) and exact (mu) over 13000-14000 K
    const scalar wgt = 0.5*(1.0 - Foam::cos(pi*(Te - 13000.0)/1000.0));
    return (1.0 - wgt)*muS + wgt*mu;
}

// ---------------------------------------------------------------------------
// 2. Electron heat capacity
// ---------------------------------------------------------------------------
inline scalar calculateCe(scalar Te, scalar EF, scalar gamma, scalar mu)
{
    const scalar kB = ttm::kB;
    const scalar pi = constant::mathematical::pi;
    const scalar kT = kB*Te;
    if (Te < 1.0) return gamma*Te;       // low-T limit

    const scalar n = 2.0*gamma*EF/(pi*pi*kB*kB);
    const scalar C = 1.5*n*Foam::pow(EF, -1.5);

    const scalar Elo = max(1e-12, mu - 60.0*kT);   // window +/-60 kT around mu
    const scalar Ehi = mu + 60.0*kT;
    const int N = 4000;
    const scalar dE = (Ehi - Elo)/N;

    scalar A0 = 0.0, A1 = 0.0, A2 = 0.0;   // 0th,1st,2nd moments of g*f(1-f)
    for (int i = 0; i <= N; ++i)
    {
        const scalar E = Elo + i*dE;
        if (E <= 0.0) continue;
        const scalar X = (E - mu)/kT;
        scalar f1mf = 0.0;                  // f(1-f), peaked near mu
        if (X > -50.0 && X < 50.0)
        {
            const scalar ex = Foam::exp(X);
            f1mf = ex/Foam::sqr(ex + 1.0);
        }
        const scalar g      = C*Foam::sqrt(E);
        const scalar wquad  = (i == 0 || i == N) ? 0.5 : 1.0;
        const scalar weight = wquad * g * f1mf * dE;

        A0 += weight;
        A1 += weight * (E - mu);
        A2 += weight * Foam::sqr(E - mu);
    }

    if (A0 < 1e-300) return gamma*Te;    // numerical guard only
    return (A2 - A1*A1/A0) / (kB*Te*Te); // variance form / (kB Te^2)
}

// ---------------------------------------------------------------------------
// 3. Electron-phonon coupling SHAPE integral (FEG)
// ---------------------------------------------------------------------------
inline scalar calculateIntegralG(scalar Te, scalar EF, scalar gamma, scalar mu)
{
    const scalar kB = ttm::kB;
    const scalar pi = constant::mathematical::pi;
    const scalar kT = kB*Te;

    const scalar n = 2.0*gamma*EF/(pi*pi*kB*kB);
    const scalar C = 1.5*n*Foam::pow(EF, -1.5);

    const scalar Elo = max(1e-6*ttm::eV, mu - 60.0*kT);
    const scalar Ehi = mu + 60.0*kT;
    const int N = 4000;
    const scalar dE = (Ehi - Elo)/N;

    scalar I = 0.0;
    for (int i = 0; i <= N; ++i)
    {
        const scalar E = Elo + i*dE;
        if (E <= 0.0) continue;
        const scalar X = (E - mu)/kT;
        scalar mdf = 0.0;                          // -df/dE
        if (X > -50.0 && X < 50.0)
        {
            const scalar ex = Foam::exp(X);
            mdf = ex/(kT*Foam::sqr(ex + 1.0));
        }
        const scalar g = C*Foam::sqrt(E);
        const scalar w = (i == 0 || i == N) ? 0.5 : 1.0;
        I += w*Foam::sqr(g)*mdf*dE;
    }
    return I;
}

// ---------------------------------------------------------------------------
// 4. Lattice heat capacity Cl (Shomate, solid/liquid branches)
// ---------------------------------------------------------------------------
inline scalar calculateShomateCl
(
    scalar Tl,
    scalar A_sol, scalar B_sol, scalar C_sol, scalar D_sol, scalar E_sol,
    scalar A_liq, scalar B_liq, scalar C_liq, scalar D_liq, scalar E_liq,
    scalar cl_mult, scalar cl_base
)
{
    scalar Tl_safe = min(Tl, 2743.0);     // saturate liquid polynomial
    const scalar t = Tl_safe / 1000.0;    // Shomate uses T/1000

    if (Tl_safe < 933.47)   // solid (Al melting point 933.47 K)
    {
        const scalar cp = A_sol + B_sol*t + C_sol*Foam::sqr(t)
                        + D_sol*Foam::pow(t,3) + E_sol/Foam::sqr(t);
        return max(cl_base, cp*cl_mult);
    }
    else                    // liquid
    {
        const scalar cp = A_liq + B_liq*t + C_liq*Foam::sqr(t)
                        + D_liq*Foam::pow(t,3) + E_liq/Foam::sqr(t);
        return max(cl_base, cp*cl_mult);
    }
}

// ---------------------------------------------------------------------------
// 5. Electron relaxation time  1/tau_e = Ae*Te^2 + Bl*Tl.
// ---------------------------------------------------------------------------
inline scalar calculateTauE(scalar Te, scalar Tl, scalar Ae, scalar Bl)
{
    return 1.0/(Ae*Foam::sqr(Te) + Bl*(Tl));
}

// ---------------------------------------------------------------------------
// 6. Electron thermal conductivity  ke = (1/3) vF^2 tau_e Ce  (kinetic theory).
// ---------------------------------------------------------------------------
inline scalar calculateKe(scalar tau_e, scalar Ce, scalar vF, scalar ke_base)
{
    const scalar ke = (1.0/3.0)*Foam::sqr(vF)*tau_e*Ce;
    return max(ke_base, ke);
}

// ---------------------------------------------------------------------------
// 7. Lattice thermal conductivity  kl ~ ke/99 
// ---------------------------------------------------------------------------
inline scalar calculateKl(scalar ke, scalar kl_base)
{
    return ke/99.0;
}

#endif