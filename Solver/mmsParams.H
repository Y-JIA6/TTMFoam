//==============================================================================
//  mmsParams.H — MMS constants (declared ONCE before the time loop).
//
//  Manufactured solution (mode m in z):
//    Te*(z,t) = T0 + Ae*cos(m*pi*z/Lz)*exp(-t/tau)
//    Tl*(z,t) = T0 + Al*cos(m*pi*z/Lz)*exp(-t/tau)
//
//  m=1 : single half-wave  -> spatial error tiny (time error dominates)
//  m=3 : three half-waves  -> spatial error large -> reveals 2nd-order in space
//
//  cos(m*pi*z/Lz) has zero z-gradient at z=0 and z=Lz for ANY integer m,
//  so the zeroGradient BCs at top/bottom remain consistent.
//==============================================================================
const scalar T0_mms  = 300.0;     // K
const scalar Ae_mms  = 1000.0;    // K
const scalar Al_mms  = 500.0;     // K
const scalar tau_mms = 1.0e-12;   // s
const scalar mMode   = 3.0;       // spatial mode (set 1 for time study, 3 for space study)
const scalar Lz_mms  = mesh.bounds().max().z() - mesh.bounds().min().z();
const scalar z0_mms  = mesh.bounds().min().z();
const scalar pi_m    = constant::mathematical::pi;
const scalar kz_mms  = mMode*pi_m/Lz_mms;     // wavenumber

// Frozen property constants (MUST match the mmsOn freeze in TTmAl.C)
const scalar Ce_c = 3.0e4;
const scalar Cl_c = 2.42e6;
const scalar ke_c = 235.0;
const scalar kl_c = 1.0e8;;//10.0;
const scalar G_c  = 0.0; //2.5e17;

// Source bracket coefficients (Laplacian term now uses kz^2, not (pi/L)^2)
const scalar coefQe = -Ae_mms*Ce_c/tau_mms
                    +  kz_mms*kz_mms*Ae_mms*ke_c
                    +  G_c*(Ae_mms - Al_mms);

const scalar coefQl = -Al_mms*Cl_c/tau_mms
                    +  kz_mms*kz_mms*Al_mms*kl_c
                    -  G_c*(Ae_mms - Al_mms);