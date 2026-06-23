# TTMFoam

**An OpenFOAM solver for two-temperature modeling of ultrashort-pulse laser ablation of metals with temperature-dependent quantum thermophysical properties.**

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![OpenFOAM 10](https://img.shields.io/badge/OpenFOAM-10-brightgreen.svg)](https://openfoam.org/version/10/)

TTMFoam is a finite-volume solver, built on OpenFOAM 10, for the
**Two-Temperature Model (TTM)** of femtosecond laser heating and ablation of metals. It evaluates the electron heat capacity `Ce(Te)`, the electron–phonon coupling factor `G(Te)`, and the dynamic optical response from direct numerical integration of the Fermi–Dirac distribution**, so it remains valid from room temperature up to dense-plasma electron temperatures — beyond the range where the Sommerfeld linearization holds.

---

## Table of contents

- [Features](#features)
- [Physics overview](#physics-overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Repository structure](#repository-structure)
- [Running a case](#running-a-case)
- [Configuring a simulation](#configuring-a-simulation)
- [Outputs](#outputs)
- [Code verification (MMS)](#code-verification-mms)
- [Validation](#validation)
- [Citing TTMFoam](#citing-ttmfoam)
- [License](#license)
- [Contact](#contact)

---
## Features

- **Quantum thermophysical properties.** `Ce(Te)` and `G(Te)` computed from the Fermi-Dirac distribution, valid to plasma temperatures.
- **Dynamic optical model.** Temperature-dependent Drude reflectivity `R(Te)` and effective penetration depth `δ_eff = δ_opt + δ_ball`.
- **Depth-resolved laser source.** Beer–Lambert deposition with exact cell-averaged absorption (`sinh` correction) and a Gaussian beam, traversed column-by-column with  an `O(Nz)` face-walking algorithm.
- **2D-axisymmetric** (wedge) geometry on a stretched, structured mesh.
- **Ablation diagnostics.** Crater depth and radius extracted on the fly using the phase-explosion criterion `Tl >= 0.9 Tc`.
- **Built-in verification mode.** A run-time switch activates a Method of Manufactured Solutions (MMS) source and an L2-error monitor.
- **Modular headers.** Each physical model lives in its own `.H` file and can be replaced without touching the solver loop.

---

## Physics overview

The solver advances the coupled energy equations :

```
Ce(Te) dTe/dt = div(ke grad Te) - G(Te)(Te - Tl) + S(r,z,t)
Cl(Tl) dTl/dt = div(kl grad Tl) + G(Te)(Te - Tl)
```

- `Te`, `Tl` — electron and lattice temperatures
- `Ce`, `Cl` — volumetric heat capacities
- `ke`, `kl` — thermal conductivities
- `G`        — electron–phonon coupling factor
- `S`        — laser source term

Boundaries are adiabatic (zero-flux), valid over the sub-20 ps timescale where radiative and convective losses are negligible.

---

## Requirements

| Dependency | Version | Notes |
|---|---|---|
| OpenFOAM | 10 | [openfoam.org](https://openfoam.org/version/10/) |
| C++ compiler | C++14 | bundled with OpenFOAM toolchain |
| MPI | any OpenFOAM-supported | for parallel runs |
| gnuplot | optional | convergence / result plots |

**Operating systems:** Linux, macOS, or Windows via WSL. OpenFOAM 10 must be installed and its environment sourced (`source /opt/openfoam10/etc/bashrc`) before building.

---

## Installation

```bash
# 1. Source the OpenFOAM environment
source /opt/openfoam10/etc/bashrc

# 2. Clone the repository
cd TTMFoam

# 3. Build the solver
cd /TTmAl
wclean
wmake
```

A successful build places the `TTmAl` executable in `$FOAM_USER_APPBIN`. Verify with:

```bash
which TTmAl
```

---

## Repository structure

```
TTMFoam/
├── TTmAl/
│   ├── TTmAl.C              # main solver: time loop + Picard loop
│   ├── createFields.H       # field declarations
│   ├── readProperties.H     # reads dictionaries at run time
│   ├── createTables.H       # precomputes Ce(Te), G(Te) tables
│   ├── updateThermo.H       # Cl, tau_e, ke, kl
│   ├── updateOptical.H      # Drude R(Te), delta_eff
│   ├── calculateLaser.H     # Beer–Lambert Gaussian source
│   ├── mmsParams.H          # MMS constants (verification mode)
│   ├── mmS.H                # MMS manufactured source
│   ├── mmSL2.H              # MMS L2-error monitor
│   └── Make/                # wmake build files
│       ├── files            
│       └── options
├── cases/
│   ├── caseAl/                  # validation case (Omeñaca 2024)
│   │   ├── 0/                   # initial/boundary fields (Te, Tl)
│   │   ├── constant/            # material & laser dictionaries
│   │   ├── plots/
│   │   │   ├── plot_Ablation.gp           # plotting ablation
│   │   │   ├── plot_All.gp                # Plot multiples Fluences peaks of temperatures Tl and Te 
│   │   │   └── plot_Temperatures.gp       # Plot one Fluence peak of temperatures Tl and Te
│   │   ├── scripts/
│   │   │   ├── Onerun.sh                  # run a physics case end to end
│   │   │   └── Allrun.sh                  # run a multiples phisics cases 
│   │   └── system/                        # blockMeshDict, controlDict, fvSchemes ...
│   ├── MMScase/                 # verification case (MMS)
│   │   ├── 0/                   # initial/boundary fields (Te, Tl)
│   │   ├── constant/            # material & laser dictionaries
│   │   ├── plots/
│   │   │   ├── spacialconvergence.gp                                   # plotting spacial convergence
│   │   │   └── plot_Temperatures.gp  timeconvergence.gp                # Plotting time convergence
│   │   ├── scripts/
│   │   │   ├── mmsrun.sh                  # spatial convergence sweep
│   │   │   └── timerun.sh                 # temporal convergence sweep 
│   └── └── system/                        # blockMeshDict, controlDict, fvSchemes ...
├── docs/                        # figures and additional documentation
├── README.md
└── License.txt
```

---

## Running a case

### Physics (ablation) case

```bash
cd cases/caseAl
blockMesh                        # generate the mesh
TTmAl                            # run in serial
# or in parallel:
decomposePar -force
mpirun -np 4 TTmAl -parallel
reconstructPar
```

At the end, the solver prints an ablation report and appends a line to
`Res.txt`:

```
Fluence | Depth_theory | Depth_numeric (nm) | Diameter_theory | Diameter_numeric (um)
```

### Convenience script

```bash
cd cases/ablationAl
../../scripts/Onerun.sh
```

---

## Configuring a simulation

All physical parameters are read from dictionaries — **no recompilation is needed** to change a case. Key entries (in `constant/` and `system/`):

| Quantity | Symbol | Example | Location |
|---|---|---|---|
| Wavelength | λ | 1032 nm | laser dict |
| Pulse duration | τ_p | 280 fs | laser dict |
| Beam waist | w0 | 10 µm | laser dict |
| Peak fluence | F | 1–9 J/cm² | laser dict |
| Ablation temperature | Tc | 6700 K | solver/dict |
| Domain radius / depth | Rmax / Zmax | 15 / 2 µm | blockMeshDict |
| Mesh | Nr × Nz | 200 × 300 | blockMeshDict |

The mesh uses a geometric grading toward the surface (≈1 nm cells at the surface, ≈25 nm in the bulk) to resolve the optical skin depth.

---

## Outputs

- **Time directories** (`0/`, `1e-13/`, …) — `Te`, `Tl`, and property fields, viewable in ParaView (`paraFoam` or `touch case.foam`).
- **`Res.txt`** : ablation depth and diameter vs. fluence.
- **Console report** : peak temperatures and crater geometry per run.

---

## Code verification (MMS)

TTMFoam includes a built-in Method of Manufactured Solutions mode for verifying the discretization independently of the physics.

Enable it by adding to `system/controlDict` :

```c
mmsVerification  true ;
```

Then run a convergence study :

```bash
cd cases/MMScase

# temporal order (expect ~ 1, backward Euler)
../../scripts/runMMS_time.sh
gnuplot ../../plots/timeconvergence.gp

# spatial order (expect ~ 2, Gauss linear)
../../scripts/runMMS.sh
gnuplot ../../plots/spacialconvergence.gp
```

**Reference results** (this solver): temporal order `p ≈ 0.99`; spatial order `p ≈ 1.97` (electron eq.) and `p ≈ 2.01` (lattice eq.), matching the theoretical orders of the backward-Euler and Gauss-linear schemes.

The manufactured solution and parameters (mode, amplitudes, frozen properties) are set in `mmsParams.H`.

---

## Validation

The `cases/caseAl` case reproduces the femtosecond aluminum ablation of Omeñaca et al. (Opt. Laser Technol. 170, 2024, 110283): λ = 1032 nm, τ_p = 280 fs, w0 = 10 µm. Simulated crater depth and squared diameter follow the expected logarithmic law

```
h(F)  = δ   · ln(F / Fth_h)
D²(F) = 2w0² · ln(F / Fth_D)
```

and agree with the reference experimental and numerical data over F = 1–9 J/cm². A simple run for one fluencce a few minutes and run completes in a few hours on four cores.

---

## Citing TTMFoam

If you use TTMFoam in your research, please cite the SoftwareX article :

```bibtex
@article{kpelly2026ttmfoam,
  title   = {TTMFoam: An OpenFOAM solver for two-temperature modeling of ultrashort-pulse laser ablation of metals with temperature-dependent quantum thermophysical properties},
  author  = {Kpelly, Koffi and Jia, Yabo},
  journal = {SoftwareX},
  year    = {2026},
  note    = {Submitted}
}
```

---

## License

Distributed under the **GNU General Public License v3 (GPLv3)** : see [`License.txt`](License.txt). This is compatible with OpenFOAM, which is itself GPLv3.

---

## Contact

**Koffi Kpelly** — koffi.kpelly@uph.fr
**Yabo Jia** - yabo.jia@uphf.fr
Issues and feature requests : please use the GitHub [issue tracker](https://github.com/<user>/TTMFoam/issues).