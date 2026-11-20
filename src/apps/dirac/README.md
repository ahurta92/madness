# Dirac Hartree-Fock Solver (`DFdriver`)

`DFdriver` is a 4-component molecular Dirac Hartree-Fock (DHF) solver implemented in the MADNESS (Multiresolution Adaptive Numerical Environment for Scientific Simulation) framework. It computes relativistic electronic structures using multiwavelet adaptive representations.

---

## Overview

`DFdriver` solves the 4-component Dirac-Fock equations:
$$\begin{pmatrix} V + J - K & c\, \boldsymbol{\sigma} \cdot \mathbf{p} \\ c\, \boldsymbol{\sigma} \cdot \mathbf{p} & V + J - K - 2mc^2 \end{pmatrix} \begin{pmatrix} \psi^L \\ \psi^S \end{pmatrix} = \varepsilon \begin{pmatrix} \psi^L \\ \psi^S \end{pmatrix}$$

Key features:

- **Initial Guess from `moldft`**: Initializes 4-component spinors from a converged non-relativistic ground-state calculation in `moldft` via the kinetic balance condition:
  $$\psi^L = \phi_{\text{NR}}, \quad \psi^S \approx \frac{-i}{2c} (\boldsymbol{\sigma} \cdot \boldsymbol{\nabla}) \phi_{\text{NR}}$$
- **BSH Integral Operator Formulation**: Employs the bound-state Helmholtz (BSH) integral representation to solve for orbital updates without finite-difference grid artifacts.
- **KAIN Acceleration**: Uses the Krylov Accelerated Inexact Newton (KAIN) method for robust convergence.
- **Flexible Nuclear Models**: Supports Gaussian, Fermi, and point-charge nuclear distributions.
- **Restart Capability**: Can restart from `moldft` restart files or resume previous `DFdriver` runs.

---

## Compilation

From your build directory:

```bash
cmake ../madness -GNinja -DENABLE_MPI=OFF \
  -DLAPACK_LIBRARIES="-L/home/rjh/Devel/mTxm/lib -lmTxmq -L/opt/arm/armpl_23.10_gcc-12.2/lib/ -larmpl_lp64 -lamath -lastring"

ninja src/apps/dirac/DFdriver
```

---

## Basic Workflow

Running a Dirac-Fock calculation is a two-step process:

### Step 1: Run `moldft` for the Non-Relativistic Ground State

Run `moldft` to generate the molecular orbital restart archive (default prefix: `mad.restartdata.00000`):

```bash
# Example: Helium atom
MAD_NUM_THREADS=20 ./src/apps/moldft/moldft --geometry=he
```

Or for a custom geometry using an input deck (e.g. Beryllium):
```
dft
  xc hf
end

geometry
  be 0.0 0.0 0.0
end
```
```bash
MAD_NUM_THREADS=20 ./src/apps/moldft/moldft
```

### Step 2: Run `DFdriver`

Create an `input` file in the working directory:

```
DiracFock
  archive mad.restartdata
  small 1e-5
  thresh 1e-6
  k 8
  kain
  maxsub 5
  no_save
end
```

Run `DFdriver`:

```bash
MAD_NUM_THREADS=20 ./src/apps/dirac/DFdriver
```

By default, `DFdriver` reads `input` in the current working directory. You can specify a different input file via:
```bash
MAD_NUM_THREADS=20 ./src/apps/dirac/DFdriver --input=my_input
```

---

## Input Options

All parameters are specified within the `DiracFock ... end` block. Lines starting with `#` are treated as comments.

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `archive <file>` | `string` | *(required)* | Path or name of the input restart archive. Accepts `mad.restartdata`, `mad.restartdata.00000`, or `mad`. |
| `job <n>` | `int` | `0` | Calculation type. `0`: DF on occupied orbitals only. |
| `thresh <val>` | `double` | `1e-6` | Wavelet truncation/refinement threshold. Supplies the default for `dconv` and `thresh_mul`. |
| `dconv <val>` | `double` | `thresh` | Density convergence threshold ($\lVert \Delta \rho \rVert$). Omitted, it inherits the final value of `thresh`; an explicit value wins regardless of keyword order. |
| `thresh_mul <val>` | `double` | `thresh` | Screening threshold for sparse multiplication (`mul_sparse`). Omitted, it inherits the final value of `thresh`; an explicit value wins regardless of keyword order. |
| `k <order>` | `int` | `8` | Multiwavelet polynomial order. If different from archive, orbitals are automatically projected. |
| `small <val>` | `double` | `1e-5` | Smallest length scale to be resolved. |
| `max_iter <n>` | `int` | `20` | Maximum number of SCF iterations. Reaching it stops the run *without* convergence; see [Termination](#termination). |
| `min_iter <n>` | `int` | `2` | Minimum number of SCF iterations. Must be less than or equal to `max_iter`. |
| `convergence_criteria <val>` | `string` | `bsh_residual` | Which quantities must converge to stop iterating. `bsh_residual`: max BSH residual &le; `thresh`. `energy_density_residual`: relative total-energy change &le; `thresh`, *and* $\lVert \Delta \rho \rVert$ &le; `dconv` &times; N<sub>elec</sub>, *and* max BSH residual &le; 100 &times; `thresh`. |
| `kain` | flag | `false` | Enable KAIN nonlinear accelerator. |
| `maxsub <n>` | `int` | `10` | Maximum subspace size for KAIN. |
| `maxrotn <val>` | `double` | `0.25` | Maximum orbital rotation step allowed by KAIN. |
| `nucleus <n>` | `int` | `0` | Nuclear charge distribution: `0` = Gaussian, `1` = Fermi, `2` = Point. |
| `restart` | flag | `false` | Resume from a previous `DFdriver` run rather than from `moldft`. |
| `savefile <file>`| `string` | `DFrestartdata`| Name of the archive to write restart data to. |
| `no_save` | flag | `false` | Disable saving restart data at the end of each iteration. |
| `lb_iter <n>` | `int` | `20` | Number of initial unconverged SCF iterations during which load balancing runs (iteration_number <= lb_iter). |
| `speed_of_light`| `double` | `137.035999...`| Speed of light in atomic units (CODATA 2022). |
| `bohr_rad <val>` | `double` | `52917.72...`  | Bohr radius in fm (CODATA 2022). |
| `Krestricted` | flag | `false` | Enforce time-reversal symmetry (Kramers pairs) for closed-shell systems. Recommended for closed shells. |
| `lineplot` | flag | `false` | Generate 1D lineplots of large and small spinor components along the x-axis. |
| `no_compute` | flag | `false` | Skip SCF iterations and exit after setup. |

### Termination

The SCF loop ends in one of two ways, and the final line of the iteration log says
which:

- **Converged.** The quantities selected by `convergence_criteria` all met their
  tolerances. The run prints `Converged due to residuals` or
  `Converged due to energy, density, and residuals`.
- **Iteration cap reached.** `max_iter` iterations ran without meeting the
  criteria. The run prints
  `WARNING: maximum iterations reached without convergence`, followed by the
  number of iterations performed.

An iteration-limited calculation is **not** converged, and its energies and
properties should not be reported as such. The latest orbitals and the restart
archive are still written (unless `no_save` is set), so the calculation can be
continued with `restart` and a larger `max_iter` rather than started over.

---

## Closed-Shell Calculations & Kramers Pairs

In relativistic electronic structure, time-reversal symmetry ($\mathcal{K}$) plays the role of spin-degeneracy in non-relativistic closed-shell systems. Each spatial orbital splits into a Kramers doublet $(\psi, \bar{\psi})$, where the time-reversed partner is:
$$\bar{\psi} = \mathcal{K} \psi = \begin{pmatrix} -\psi_1^* \\ \psi_0^* \\ -\psi_3^* \\ \psi_2^* \end{pmatrix}$$

The performance observations below are historical benchmarks; they depend on the system, numerical settings, and execution environment and were not remeasured during this documentation bootstrap. By specifying `Krestricted` in the input deck, `DFdriver`:

- **Reduces computational cost significantly**: Solves and iterates over only $N_{\text{pairs}}$ independent 4-spinors (half the total electrons). This yields ~2× fewer 3D Poisson solves in exchange ($N(N+1)$ vs $2N^2+N$), an exact 2× reduction in BSH operator applications, and 50% lower memory usage compared to an unrestricted calculation.
- **Generates time-reversed partners algebraically**: $\bar{\psi}$ is constructed on the fly without numerical integration, and the quaternionic Fock matrix is formed and diagonalized directly.
- **Preserves Kramers degeneracies**: Exact partner degeneracy ($\varepsilon_i = \varepsilon_{\bar{i}}$) is maintained to $< 10^{-6}$ Ha.
- **Preserves spherical symmetry**: For atoms with closed subshells (e.g. $p^6$ in neon and magnesium), it preserves spherical isotropy and correctly resolves the physical spin-orbit splitting ($p_{1/2}$ doublet vs $p_{3/2}$ quartet) without unphysical $z$-polarization.

Without `Krestricted`, a closed-shell system is treated as Kramers-unrestricted: both members of each
Kramers doublet are stored and iterated explicitly, so a guess from a spin-restricted `moldft` run gives
$2N$ singly-occupied 4-spinors rather than $N$ doubly-occupied ones. This costs about twice as much as
`Krestricted` but does not impose time-reversal symmetry on the solution.

---

## Restarting Calculations

### From `moldft`
Provide the base name or chunk filename of the `moldft` calculation:
```
DiracFock
  archive mad.restartdata
  Krestricted
  ...
end
```
Both `archive mad.restartdata` and `archive mad.restartdata.00000` are accepted.

### From a Previous `DFdriver` Calculation
To continue an interrupted or unconverged `DFdriver` calculation:
```
DiracFock
  archive DFrestartdata
  restart
  Krestricted
  kain
  max_iter 20
end
```
Ensure the previous run did not specify `no_save`.

DF restart archives begin with an explicit format version field. Archives written before this
field was introduced carry no version and are rejected with an actionable error; they are not
migrated. Regenerate them by restarting from a `moldft` archive. The incompatibility is detected
in both directions: an older `DFdriver` pointed at a new archive reads the version field where it
expects the total energy, so rank 0 reports an archive type mismatch instead of misparsing the
file. Multi-rank runs may hang at that point, because the other ranks wait for a broadcast from
rank 0 that never comes.

---

## Example Inputs

### 1. Standard Closed-Shell Calculation (`DFinput_sample`)
```
DiracFock
  archive mad.restartdata
  Krestricted
  small 1e-5
  thresh 1e-6
  k 8
  kain
  maxsub 5
  no_save
end
```

### 2. High-Precision Calculation with Fermi Nuclear Model
```
DiracFock
  archive /path/to/calc.restartdata
  thresh 1e-7
  dconv 1e-7
  k 10
  nucleus 1
  kain
  maxsub 8
  savefile my_df_checkpoint
end
```

<!-- EGM architecture contract: begin -->

## Architecture overview

Active four-component Dirac–Fock application with versioned restart/convergence helpers and dormant Douglas–Kroll research code.

## Interface

- DFdriver builds DFdriver.cc, DF.cc, fcwf.cc and links MADchem; DF::solve dispatches job=0 occupied solve. [CMakeLists.txt:add_mad_executable](./CMakeLists.txt) [DF.cc:DF::solve](./DF.cc)
- DFParameters reads DiracFock block; InitParameters reads moldft or DF archive; Fcwf supplies four-component arithmetic, density, KramersPair, transform and KAIN allocator. [DFParameters.h:DFParameters](./DFParameters.h) [InitParameters.h:InitParameters](./InitParameters.h) [fcwf.h:Fcwf](./fcwf.h)
- DKops.h,relops.cc,rk.cc are dormant experimental relativistic operator/fitting/driver sources, excluded from DFdriver. [CMakeLists.txt:add_mad_executable](./CMakeLists.txt) [DKops.h:Ebar](./DKops.h) [relops.cc:bshrel_fit](./relops.cc) [rk.cc:main](./rk.cc)

## Invariants

- DF construction is collective; min_iter<=max_iter always checked; dconv and thresh_mul inherit final thresh unless explicitly set. [DF.cc:DF::DF](./DF.cc) [DFParameters.h:DFParameters::read](./DFParameters.h)
- Fcwf stores small components multiplied by c, so norms/inner products/densities divide their contribution by c²; assignment is shallow and copy(Fcwf) is explicit deep copy. [fcwf.cc:Fcwf::norm2](./fcwf.cc) [fcwf.cc:Fcwf::operator=](./fcwf.cc) [fcwf.cc:copy](./fcwf.cc)
- DF restart version=1 required; rank zero probes then broadcasts version so all ranks reject together; stored Kramers restriction must match request; moldft archives representing nemo/znemo are rejected. [DFRestart.h:DF_RESTART_VERSION](./DFRestart.h) [InitParameters.h:InitParameters::read](./InitParameters.h)
- bsh_residual requires max residual<=thresh; combined criterion additionally requires relative energy change<=thresh and density change<=dconv*nelec, with residual<=100*thresh; max_iter exit is explicitly unconverged. [DFConvergence.h:df_iteration_converged](./DFConvergence.h) [DF.cc:DF::solve_occupied](./DF.cc)

## Side Effects

- Reads input and orbital archives, sets FunctionDefaults and reprojects to requested k, balances orbital/potential trees, saves DFrestartdata (or savefile) each iteration unless no_save. [DF.cc:DF::DF](./DF.cc) [DF.cc:DF::saveDF](./DF.cc) [DF.cc:DF::solve_occupied](./DF.cc)
- Prints orbital residuals, energies/components, radial expectation, coefficient counts/depth and timing; lineplot writes large_component_lineplots and small_component_lineplots. [DF.cc:DF::iterate](./DF.cc) [DF.cc:DF::solve](./DF.cc)
- Dormant DKops factories read absolute /gpfs/home/jscanderson/DKproject CSV paths and exit on missing data; relops includes external ../../../../mpfrc++-3/mpreal.h; rk writes logarithmic plots. [DKops.h:Pbar](./DKops.h) [relops.cc:bshrel_omega](./relops.cc) [rk.cc:logplot](./rk.cc)

## Source files

| File | Role |
| --- | --- |
| [DF.h](./DF.h) | DF solver API. |
| [DFConvergence.h](./DFConvergence.h) | Pure criterion/termination reporting. |
| [DFParameters.h](./DFParameters.h) | Input parsing and normalization. |
| [DFRestart.h](./DFRestart.h) | Restart format validation. |
| [InitParameters.h](./InitParameters.h) | Archive import/kinetic-balance initial spinors; readnw explicitly marked obsolete. |
| [DFdriver.cc](./DFdriver.cc) | Active CLI. |
| [DF.cc](./DF.cc) | SCF algorithm, energies, restart and plot effects. |
| [fcwf.h](./fcwf.h) | Spinor and allocator declarations. |
| [fcwf.cc](./fcwf.cc) | Scaled-small-component spinor operations. |
| [DKops.h](./DKops.h) | Dormant hard-coded CSV operator factories. |
| [relops.cc](./relops.cc) | Dormant Gaussian expansion fitting using MPFR. |
| [rk.cc](./rk.cc) | Dormant one-electron relativistic kinetic/DK prototype. |

## Verification

Commands are run from the repository root against an existing configured build. Replace `<build-dir>` with that build directory.

```sh
cmake --build <build-dir> --target DFdriver test_dirac
```

```sh
ctest --test-dir <build-dir> -R "madness/test/dirac/test_dirac/run"
```

<!-- EGM architecture contract: end -->
