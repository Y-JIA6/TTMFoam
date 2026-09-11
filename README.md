# TTMFoam

**An OpenFOAM solver for two-temperature modelling of ultrashort-pulse laser heating and ablation diagnostics in aluminium.**

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)
[![OpenFOAM 10](https://img.shields.io/badge/OpenFOAM-10-brightgreen.svg)](https://openfoam.org/version/10/)

TTMFoam is a finite-volume solver for the two-temperature model (TTM) of ultrashort-pulse laser--metal interaction. It is implemented in OpenFOAM 10 and is supplied with an axisymmetric aluminium validation case and a Method of Manufactured Solutions (MMS) verification case. The code evaluates temperature-dependent electron properties from free-electron Fermi--Dirac integrals rather than assuming a Sommerfeld-linear heat capacity.

## Main features

- Fermi--Dirac lookup tables for the chemical potential `mu(Te)`, electron heat capacity `Ce(Te)`, and electron--phonon coupling factor `G(Te)`.
- Temperature-dependent electron and lattice properties, including a Shomate lattice heat capacity and relaxation-time-based thermal conductivities.
- A calibrated Drude optical response at 1032 nm. The ambient reflectivity is anchored to `R = 0.91`; the temperature dependence follows the Drude relaxation time without a fluence-dependent empirical absorption correction.
- A Gaussian, depth-resolved Beer--Lambert laser source with exact cell averaging.
- A 2D axisymmetric wedge mesh and MPI parallel execution.
- On-the-fly crater diagnostics based on the threshold `Tl >= 0.9 x 6700 K` (6030 K).
- An automated multi-fluence campaign that extends each run until the crater diagnostic is quiet.
- MMS spatial and temporal convergence studies.

## Governing model

TTMFoam solves

```text
Ce(Te) dTe/dt = div(ke grad Te) - G(Te) (Te - Tl) + Qlaser
Cl(Tl) dTl/dt = div(kl grad Tl) + G(Te) (Te - Tl)
```

where `Te` and `Tl` are the electron and lattice temperatures, `Ce` and `Cl` are volumetric heat capacities, `ke` and `kl` are thermal conductivities, `G` is the electron--phonon coupling factor, and `Qlaser` is the volumetric laser source.

For the supplied aluminium case, `Ce`, `G`, and `mu` are tabulated once at start-up on a 10 K grid from 300 K to 300 000 K and then linearly interpolated in each cell. The simulation log reports whether the peak electron temperature remains inside this table range.

The laser source uses the surface Drude absorptivity, `A = 1 - R`, and a local attenuation length

```text
delta_eff = delta_opt + vF tau_e.
```

Here, `R` is evaluated at the irradiated surface and `delta_eff` is evaluated cell by cell along the Beer--Lambert path. The ballistic contribution is therefore treated as a transport-length correction to the deposition source, not as an independent optical property.

## Scope and limitations

The supplied model is intended for ultrashort-pulse aluminium simulations within the assumptions below.

- The TTM presumes that the electron population can be represented by an electron temperature. Non-thermal electron kinetics and hyperbolic TTM physics are not implemented.
- The Fermi--Dirac properties use a free-electron density of states. Material-specific *ab initio* densities of states, especially for transition metals with d bands, are not included.
- The optical response is a calibrated Drude model. Interband transitions, Drude--Lorentz terms, and explicit collisionless absorption are outside the present implementation.
- The code diagnoses a crater from a lattice-temperature threshold; it does not remove material or solve hydrodynamic expansion.
- The distributed material parameters and validation case are for aluminium. Applying the solver to another metal requires appropriate material parameters and validation.

## Requirements

| Dependency | Requirement | Purpose |
|---|---:|---|
| OpenFOAM | 10 | finite-volume framework and build system |
| C++ compiler | C++14-compatible | compilation through `wmake` |
| MPI | OpenFOAM-compatible | parallel validation runs |
| Bash and awk | standard POSIX tools | automation scripts |
| gnuplot | optional | figures from the supplied plotting scripts |

OpenFOAM must be installed and its environment sourced before compilation, for example:

```bash
source /opt/openfoam10/etc/bashrc
```

## Installation

Clone or download the repository, then compile the solver from the repository root:

```bash
cd TTMFoam
source /opt/openfoam10/etc/bashrc
cd Solver
wclean
wmake
which TTmAl
```

The executable is installed in `$FOAM_USER_APPBIN`. The final command should print the path to `TTmAl`.

## Repository structure

```text
TTMFoam/
├── Solver/
│   ├── TTmAl.C                 # main time loop, Picard iteration, diagnostics
│   ├── readProperties.H         # case dictionaries
│   ├── createFields.H           # OpenFOAM fields
│   ├── createTables.H           # FEG property tables
│   ├── updateThermo.H           # thermophysical properties and relaxation time
│   ├── updateOptical.H          # calibrated Drude R and delta_opt
│   ├── calculateLaser.H         # Beer--Lambert laser source
│   ├── mmsParams.H              # manufactured solution and MMS parameters
│   ├── mmS.H                    # manufactured source terms
│   ├── mmSL2.H                  # MMS L2-error evaluation
│   └── Make/                    # wmake configuration
├── Cases/
│   ├── caseAl/                  # aluminium validation case
│   │   ├── 0/, constant/, system/
│   │   ├── plots/               # gnuplot scripts
│   │   └── scripts/
│   │       ├── Onerun.sh        # one parallel run on four MPI ranks
│   │       ├── Serierun.sh      # one serial run
│   │       ├── ConvergedRun.sh  # adaptive end-time convergence for one fluence
│   │       └── Allrun.sh        # adaptive 1--9 J/cm2 campaign
│   └── MMScase/                 # Method of Manufactured Solutions case
│       ├── 0/, constant/, system/, plots/
│       └── scripts/mmsrun.sh, scripts/timerun.sh
├── manuscript/                  # manuscript and response-letter sources
└── README.md
```

## Running the aluminium case

All commands below assume that `TTmAl` has already been compiled and OpenFOAM has been sourced.

### Parallel run (four ranks)

```bash
cd TTMFoam/Cases/caseAl
bash scripts/Onerun.sh
```

`Onerun.sh` clears previous time directories, builds the mesh, decomposes the domain, and runs `TTmAl` on four MPI ranks. It does not recompile the solver.

### Serial run

```bash
cd TTMFoam/Cases/caseAl
bash scripts/Serierun.sh
```

This is useful for debugging, MPI comparisons, or machines without an MPI launch configuration.

### Adaptive multi-fluence campaign

```bash
cd TTMFoam/Cases/caseAl
bash scripts/Allrun.sh
```

The campaign evaluates fluences from 1 to 9 J/cm2. The 1 J/cm2 case begins at 20 ps. If crater growth remains significant, the same case is repeated from `t = 0` with an end time increased by 10 ps. Once accepted, the next fluence begins at the previous accepted end time plus 10 ps. The default safety limit is 200 ps.

The crater diagnostic uses a 0.5 nm event-detection threshold and requires 10 ps without significant depth or radius growth. It also compares two successive endpoints: the depth and diameter changes must be smaller than the corresponding local mesh resolutions. Only the accepted result for each fluence is retained in `Res.txt`; all attempted endpoints are recorded in `convergence_history/crater_geometry_convergence.tsv`.

To run the same procedure for one already configured fluence:

```bash
bash scripts/ConvergedRun.sh 20 10 10 200
```

The four arguments are, respectively, the initial end time, the extension increment, the required quiet-time window, and the maximum end time, all in ps.

## Configuring a simulation

Most case inputs are read at run time and can be changed without recompiling:

| Input | Location |
|---|---|
| Fluence, wavelength, pulse duration, beam waist | `Cases/caseAl/constant/laserProperties` |
| Aluminium thermophysical and relaxation parameters | `Cases/caseAl/constant/ttmProperties` |
| End time, write interval, Picard tolerances | `Cases/caseAl/system/controlDict` |
| Mesh dimensions and grading | `Cases/caseAl/system/blockMeshDict` |
| Linear solver controls | `Cases/caseAl/system/fvSolution` |
| Surface and axial probes | `Cases/caseAl/system/probes` |

Changing a solver-level model, the FEG table range, or the crater threshold requires modifying the corresponding source file in `Solver/` and rebuilding with `wmake`.

## Outputs and figures

The principal outputs of the aluminium case are:

- `postProcessing/probes/` for histories of `Te`, `Tl`, `Qlaser`, `Ce_var`, `ke_var`, `G_var`, `R_field`, `deltaOpt_field`, `deltaEff_field`, and `tauE_field`.
- `FEG_param.dat` for the generated FEG property tables.
- `Res.txt` for the accepted fluence, crater depth, and crater diameter values.
- `convergence_history/` for individual endpoint logs and the convergence history.
- `resultats_fluences/` for probe histories retained by `Allrun.sh`.

After a completed run, selected figures can be generated with:

```bash
gnuplot plots/plot_Temperature.gp
gnuplot plots/plot_All.gp
gnuplot plots/plot_Ablation.gp
gnuplot plots/plot_params.gp
gnuplot plots/plot_Optical.gp
```

## Code verification with MMS

The MMS case verifies the transient and diffusion discretization independently of the laser-ablation physics. Its `system/controlDict` enables `mmsVerification true;`.

```bash
cd TTMFoam/Cases/MMScase
bash scripts/timerun.sh
gnuplot plots/timeconvergence.gp

bash scripts/mmsrun.sh
gnuplot plots/spacialconvergence.gp
```

The temporal study refines `deltaT` on a fixed mesh; the spatial study refines the axial mesh with fixed time settings. Both scripts modify their case dictionaries during the sweep. Restore the tracked case files afterwards if you need the original configuration:

```bash
git restore system/controlDict system/blockMeshDict
```

The expected orders are approximately one for backward Euler in time and two for the Gauss-linear diffusion discretization in space.

## Aluminium validation case

The supplied case uses the experimental conditions reported by Omeñaca *et al.* for femtosecond laser ablation of aluminium: wavelength 1032 nm, pulse duration 280 fs, and beam waist 10 um. It is intended as a reproducible validation and comparison case; the interpretation of crater depth remains subject to the scope limitations stated above.

## Citation

If TTMFoam contributes to your work, please cite the associated SoftwareX manuscript. Replace the placeholder with the final bibliographic information and DOI once available.

```bibtex
@article{kpelly2026ttmfoam,
  title   = {TTMFoam: Two-temperature modeling of ultrashort-pulse laser ablation of metals with temperature-dependent quantum thermophysical properties},
  author  = {Kpelly, K, Ivan and Jia, Yabo},
  journal = {SoftwareX},
  year    = {2026},
  note    = {Manuscript under revision}
}
```

## License

TTMFoam is distributed under the GNU General Public License v3.0 (GPL-3.0-or-later).

## Contact

- Koffi Kpelly — ivan.kpelly@icloud.com
- Yabo Jia — yabo.jia@uphf.fr

For questions, bug reports, or feature requests, please open an issue in the repository.
