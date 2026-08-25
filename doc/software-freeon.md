# FreeON — An Exhaustive Review

## 1. Overview

**FreeON** is an experimental, open-source (GPL) suite of programs for linear-scaling (*O(N)*) quantum chemistry, targeted at large molecular and periodic systems. It was formerly known as **MondoSCF**, and was rebranded and released as FreeON in 2010–2011. It grew out of computational chemistry work at Los Alamos National Laboratory (LANL), primarily under the leadership of **Matt Challacombe**, with major contributions from **Eric Schwegler**, **C.J. Tymczak**, **Anders M.N. Niklasson**, **Valéry Weber**, **Karoly Nemeth**, **Chee Kwan Gan**, and others.

| Attribute | Detail |
|---|---|
| Formerly known as | MondoSCF |
| First public FreeON release | ~2010–2011 (rebrand of MondoSCF) |
| Latest stable release | 1.0.8 (November 8, 2013) |
| License | GNU GPL (v2/v3 depending on source) |
| Languages | Fortran 95 (majority), C (HDF5/interfacing), some Python tooling |
| Build system | Autotools (`autogen.sh`, `configure`, `make`) and a CMake path |
| Parallelism | MPI (distributed memory) |
| I/O | HDF5 (platform-independent binary storage) |
| Primary developer / origin | Los Alamos National Laboratory (LANL); Matt Challacombe et al. |
| Repository | github.com/FreeON/freeon (mirrored from the original Savannah/Launchpad project) |
| Homepage (legacy) | nongnu.org/freeon |
| Project status | Largely dormant/archival since ~2014; no active maintenance |

## 2. Scientific Scope and Methods

FreeON performs ground-state electronic structure calculations using Gaussian-basis LCAO methods, with the specific design goal of achieving **linear, or near-linear, scaling with system size** for non-metallic (insulating/gapped) systems. Supported levels of theory and features include:

- **Hartree–Fock (HF)** theory
- **Pure density functional theory (DFT)** (local and gradient-corrected functionals)
- **Hybrid HF/DFT** functionals, e.g. B3LYP
- **Cartesian-Gaussian LCAO basis sets**
- **Periodic boundary conditions** in 1, 2, and 3 dimensions, implemented via a Lorentz-field (Γ-point-only) treatment of the periodic Coulomb and exchange problems
- **Geometry optimization**, including a full atom + unit-cell relaxation using an internal/curvilinear coordinate optimizer with analytic energy gradients with respect to both atomic positions and cell parameters
- **Born–Oppenheimer and time-reversible ab initio molecular dynamics**
- **Time-dependent SCF / perturbation-theory response properties** (e.g., polarizabilities, excited states) via density-matrix-based (wavefunction-free) response theory rather than conventional TD-DFT formalisms
- **Density matrix purification / density matrix minimization** as the central O(N) replacement for canonical diagonalization of the Fock/Kohn–Sham matrix

### 2.1 Core algorithmic strategy

FreeON's linear scaling is achieved through several coupled techniques rather than a single trick:

1. **QCTC (Quantum Chemical Tree Code)** — a hierarchical, octree-like fast multipole/tree method for assembling the Coulomb matrix, analogous in spirit to N-body tree codes from astrophysics, adapted to Gaussian charge distributions.
2. **ONX (O(N) Exchange)** — rigorous distance/multipole-based bounds and screening to construct the Hartree–Fock exchange matrix with linear (or N log N) cost, including permutational symmetry exploitation and incremental Fock builds.
3. **HiCu (Hierarchical Cubature)** — a hierarchical numerical integration scheme for the exchange–correlation matrix in DFT, built for efficient, adaptive quadrature over the Gaussian charge/orbital density.
4. **Density matrix purification** (McWeeny-type and trace-resetting variants), avoiding full diagonalization of the Fock matrix — replacing it with sparse matrix–matrix multiplications, which is the step that ultimately fixes the overall O(N) or O(N log N) scaling for gapped systems.
5. **Sparse-blocked parallel matrix algebra** — a general parallel sparse-blocked matrix multiplication library underlying the purification and Fock-build steps, enabling distributed-memory scaling.
6. **Density matrix perturbation theory (DMPT)** — a purification-based analogue of coupled-perturbed SCF theory, allowing linear-scaling computation of response properties (polarizabilities, NMR-type response, TD perturbation theory) without explicit virtual orbitals.
7. **Curvilinear/internal coordinate geometry optimization**, including a "quasi-independent curvilinear coordinate approximation" designed to make optimization scale well for large and periodic systems (crystals).
8. **Time-reversible Born–Oppenheimer MD**, an integration scheme that preserves long-term energy conservation in ab initio MD without full SCF convergence at every step.

### 2.2 Internal module structure (as reflected in the source tree)

The codebase is organized into functionally distinct Fortran/C modules, each corresponding to a stage of the SCF pipeline:

- `FreeON` — main driver/front end
- `OneE` — one-electron integrals
- `TwoE` — two-electron integral infrastructure
- `QCTC` — Coulomb matrix tree-code assembly
- `ONX` — exchange matrix construction
- `HiCu` — hierarchical cubature for XC integration
- `SCFeqs` — SCF equations / purification / density matrix machinery
- `Modules` — shared utility and data-structure modules
- `BasisSets` — bundled Gaussian basis set library
- `Validate` — regression/validation test suite
- `MPI_examples`, `tests`, `examples` — parallel usage examples and test inputs

## 3. Distinguishing Characteristics

- **Genuinely linear-scaling, not just "large-system capable."** Unlike many quantum chemistry packages that bolt O(N) methods onto a conventional diagonalization-based SCF engine, FreeON was designed from the ground up around O(N)/O(N log N) algorithms for every major SCF step (Coulomb, exchange, XC integration, density matrix construction).
- **Purification-first design.** Density matrix purification is the default route to the ground-state density, rather than an optional accelerator, which is what allows the code to bypass the O(N³) diagonalization bottleneck.
- **Analytic full-cell relaxation.** Its Γ-point periodic implementation includes analytic gradients with respect to *both* atomic coordinates and unit-cell parameters, enabling full crystal structure relaxation — a nontrivial feature for a Gaussian-basis periodic code of this era.
- **Response theory without virtual orbitals.** Its time-dependent/perturbation-theory formalism is built as an extension of density-matrix perturbation theory, making it naturally compatible with O(N) purification rather than requiring canonical MO-based coupled-perturbed equations.
- **Applications record.** Beyond methodology papers, the code was used for real large-scale/periodic applications, including density-functional studies of hydrostatic compression of the energetic material PETN (pentaerythritol tetranitrate), and even an application to helium conductivity in cool white dwarf atmospheres (an unusual excursion into astrophysical materials modeling).
- **Platform-independent I/O.** Uses HDF5 for binary storage/checkpointing rather than ad hoc formats, aiding portability and large-dataset handling.

## 4. Practical / Engineering Aspects

- **Build system:** GNU Autotools workflow (`./reconfigure.sh` → `./configure` → `make` → `make install`), with `--prefix` for install location; a partial CMake-based build path also exists in the tree (`CMakeLists.txt`, `cmake-tests`).
- **Cleaning targets:** graduated cleanup levels (`quickclean`, `clean`, `distclean`, `maintainer-clean`) reflecting a fairly heavyweight, legacy-style Autotools build.
- **Execution model:** invoked as `FreeON something.inp`, optionally specifying output/log/geometry filenames explicitly; runs under MPI for parallel execution.
- **Supported platforms:** Linux and Unix-like systems (including FreeBSD); OS X support was historically claimed as well.
- **Dependencies:** MPI implementation, HDF5 libraries, a Fortran 95 compiler, standard linear algebra (BLAS/LAPACK-class) libraries for local dense-block operations within the sparse-blocked matrix framework.
- **Distribution packaging:** historically packaged for Debian-family distributions (a `debian/` packaging directory exists in the repository).

## 5. Project History and Current Status

- **Origins:** development traces to the mid-1990s, growing out of Matt Challacombe and coworkers' work at the University of Minnesota and subsequently Los Alamos National Laboratory, initially released and known as **MondoSCF**.
- **Rebranding:** the project was renamed **FreeON** and hosted on Savannah (nongnu.org) and later mirrored to Launchpad and GitHub.
- **Peak activity:** roughly mid-1990s through the mid-2000s for core methodology, with continued development (time-reversible MD, response theory, DMPT) through the late 2000s.
- **Last stable release:** version 1.0.8, dated November 8, 2013.
- **Present state:** the GitHub mirror (`FreeON/freeon`) shows a low level of ongoing activity (on the order of single-digit stars/forks and open issues), consistent with an archival/legacy project rather than an actively maintained one. It is generally regarded in the electronic-structure community as one of the earlier open-source entrants in the O(N) quantum chemistry space, alongside codes such as CONQUEST, ONETEP, SIESTA, OpenMX, and — closer in spirit (Gaussian-basis, all-electron, purification-based) — **ErgoSCF**, which is often viewed as FreeON's closest methodological relative and, to some extent, spiritual successor in terms of continued maintenance.

## 6. Comparison to Related Linear-Scaling Codes

| Code | Basis type | Primary O(N) mechanism | Status (relative) |
|---|---|---|---|
| **FreeON** (MondoSCF) | Gaussian (all-electron) | Density matrix purification + tree-code Coulomb + O(N) exchange | Archival/legacy |
| **ErgoSCF** | Gaussian (all-electron) | Trace-correcting purification + fast multipole + sparse matrix algebra | Actively maintained |
| **CONQUEST** | Numerical/pseudo-atomic orbitals | Density matrix / localized orbital methods | Actively maintained |
| **ONETEP** | Psinc/plane-wave-like local basis | Localized orbitals, linear-scaling DFT | Actively maintained |
| **SIESTA** | Numerical atomic orbitals | Efficient diagonalization + O(N) options | Actively maintained |
| **OpenMX** | Pseudo-atomic orbitals | Divide-and-conquer / recursion methods | Actively maintained |

FreeON and ErgoSCF are frequently grouped together in the literature as the two Gaussian-basis, quantum-chemistry-lineage entries in the O(N) electronic structure landscape, in contrast to the numerical-atomic-orbital or plane-wave-derived codes more common in the condensed-matter/materials community.

## 7. Practical Caveats for Prospective Users

- The code is **experimental** by its own description, and the last tagged stable release is over a decade old (1.0.8, 2013) — expect friction building against modern compilers, MPI stacks, and HDF5 versions without patching.
- Documentation is comparatively sparse relative to mainstream packages (Gaussian, ORCA, Psi4, NWChem, Q-Chem); most authoritative technical detail exists in the primary literature (Section 8) rather than in user manuals.
- Periodic DFT/HF is restricted to the **Γ-point** — there is no full Brillouin-zone k-point sampling, which limits its applicability for band-structure-sensitive periodic/solid-state problems relative to plane-wave codes.
- Best suited to **large, non-metallic (gapped) molecular or periodic systems**, since the linear-scaling machinery (density matrix purification, integral screening) relies on the exponential decay of the density matrix that holds for insulators/molecules with a gap, and degrades for metallic or near-metallic systems.
- Given its dormant maintenance status, it is more likely to be encountered today as a **historical/methodological reference implementation** for O(N) SCF techniques than as a first choice for new production calculations.

## 8. Publications Related to FreeON's Theory

The list below compiles the primary methodology papers underlying FreeON (as MondoSCF/FreeON), covering Coulomb/exchange matrix construction, density matrix purification and perturbation theory, periodic boundary conditions, geometry optimization, molecular dynamics, and applications, drawn from the project's own publication list and cross-referenced against the literature.

### 8.1 Foundational Coulomb/Exchange matrix methods

- Challacombe, M.; Schwegler, E.; Almlöf, J. **"Fast assembly of the Coulomb matrix: A quantum chemical tree code."** *Journal of Chemical Physics*, 104, 4685 (1996).
- Challacombe, M.; Schwegler, E. **"Linear scaling computation of the Fock matrix."** *Journal of Chemical Physics*, 106, 5526 (1997).
- Schwegler, E.; Challacombe, M.; Head-Gordon, M. **"Linear scaling computation of the Fock matrix. II. Rigorous bounds on exchange integrals and incremental Fock build."** *Journal of Chemical Physics*, 106, 9708 (1997).
- Schwegler, E.; Challacombe, M. **"Linear scaling computation of the Hartree–Fock exchange matrix."** *Journal of Chemical Physics*, 105, 2726 (1996).
- Schwegler, E.; Challacombe, M.; Head-Gordon, M. **"A multipole acceptability criterion for electronic structure theory."** *Journal of Chemical Physics*, 109, 8764 (1998).
- Schwegler, E.; Challacombe, M. **"Linear scaling computation of the Fock matrix. IV. Multipole accelerated formation of the exchange matrix."** *Journal of Chemical Physics*, 111, 6223 (1999).
- Schwegler, E.; Challacombe, M. **"Linear scaling computation of the Fock matrix. III. Formation of the exchange matrix with permutational symmetry."** *Theoretical Chemistry Accounts*, 104, 344 (2000).

### 8.2 Density matrix purification, minimization, and perturbation theory

- Challacombe, M. **"A simplified density matrix minimization for linear scaling self-consistent field theory."** *Journal of Chemical Physics*, 110, 2332 (1999).
- Niklasson, A. M. N.; Tymczak, C. J.; Challacombe, M. **"Trace resetting density matrix purification in O(N) self-consistent-field theory."** *Journal of Chemical Physics*, 118, 8611 (2003).
- Niklasson, A. M. N.; Challacombe, M. **"Density matrix perturbation theory."** *Physical Review Letters*, 92, 193001 (2004).
- Weber, V.; Niklasson, A. M. N.; Challacombe, M. **"Ab initio linear scaling response theory: Electric polarizability by perturbed projection."** *Physical Review Letters*, 92, 193002 (2004).
- Weber, V.; Niklasson, A. M. N.; Challacombe, M. **"Higher-order response in O(N) by perturbed projection."** *Journal of Chemical Physics*, 123, 044106 (2005).
- Niklasson, A. M. N.; Weber, V.; Challacombe, M. **"Nonorthogonal density-matrix perturbation theory."** *Journal of Chemical Physics*, 123, 044107 (2005).

### 8.3 Numerical integration and parallel matrix algebra

- Challacombe, M. **"Linear scaling computation of the Fock matrix. V. Hierarchical Cubature for numerical integration of the exchange-correlation matrix."** *Journal of Chemical Physics*, 113, 10037 (2000).
- Challacombe, M. **"General parallel sparse-blocked matrix multiply for linear scaling SCF theory."** *Computer Physics Communications*, 128, 93 (2000).
- Gan, C. K.; Challacombe, M. **"Linear scaling computation of the Fock matrix. VI. Data parallel computation of the exchange-correlation matrix."** *Journal of Chemical Physics*, 118, 9128 (2003).
- Gan, C. K.; Tymczak, C. J.; Challacombe, M. **"Linear scaling computation of the Fock matrix. VII. Parallel computation of the Coulomb matrix."** *Journal of Chemical Physics*, 121, 6608 (2004).
- Weber, V.; Challacombe, M. **"Parallel algorithm for the computation of the Hartree-Fock exchange matrix: gas phase and periodic parallel ONX."** *Journal of Chemical Physics*, 125, 104110 (2006).

### 8.4 Periodic boundary conditions

- Challacombe, M.; White, C.; Head-Gordon, M. **"Periodic boundary conditions and the fast multipole method."** *Journal of Chemical Physics*, 107, 10131 (1997).
- Tymczak, C. J.; Challacombe, M. **"Linear scaling computation of the Fock matrix. VII. Periodic density functional theory at the Gamma point."** *Journal of Chemical Physics*, 122, 134102 (2005). *(Note: this and the item above share sequence number VII in the original numbering as published.)*
- Tymczak, C. J.; Weber, V.; Schwegler, E.; Challacombe, M. **"Linear scaling computation of the Fock matrix. VIII. Periodic boundaries for exact exchange at the Gamma point."** *Journal of Chemical Physics*, 122, 124105 (2005).
- Weber, V.; Daul, C.; Challacombe, M. **"Exchange energy gradients with respect to atomic positions and cell parameters within the Hartree-Fock Γ-point approximation."** *Journal of Chemical Physics*, 124, 214105 (2006).
- Weber, V.; Tymczak, C. J.; Challacombe, M. **"Energy gradients with respect to atomic positions and cell parameters for the Kohn-Sham density-functional theory at the Gamma point."** *Journal of Chemical Physics*, 124, 224107 (2006).

### 8.5 Geometry optimization

- Nemeth, K.; Challacombe, M. **"The quasi-independent curvilinear coordinate approximation for geometry optimization."** *Journal of Chemical Physics*, 121, 2877 (2004).
- Nemeth, K.; Challacombe, M. **"Geometry optimization of crystals by the quasi-independent curvilinear coordinate approximation."** *Journal of Chemical Physics*, 123, 194112 (2005).

### 8.6 Ab initio molecular dynamics

- Niklasson, A. M. N.; Tymczak, C. J.; Challacombe, M. **"Time-reversible Born-Oppenheimer molecular dynamics."** *Physical Review Letters*, 97, 123001 (2006).
- Niklasson, A. M. N.; Tymczak, C. J.; Challacombe, M. **"Time-reversible ab initio molecular dynamics."** *Journal of Chemical Physics*, 126, 144103 (2007).

### 8.7 Time-dependent / excited-state and response theory

- Challacombe, M. **"Linear Scaling Solution of the Time-Dependent Self-Consistent-Field Equations."** arXiv:1001.2586 (2010).
- Lucero, M.; Niklasson, A. M. N.; Tretiak, S.; Challacombe, M. **"Molecular orbital free algorithm for excited states in time-dependent perturbation theory."** *Journal of Chemical Physics*, 129, 064114 (2008).
- Tretiak, S.; Isborn, C.; Niklasson, A. M. N.; Challacombe, M. **"Representation independent algorithms for molecular response calculations in time-dependent self-consistent field theories."** *Journal of Chemical Physics*, 130, 054111 (2009).

### 8.8 Applications

- Gan, C. K.; Sewell, T.; Challacombe, M. **"All-electron density-functional studies of hydrostatic compression of pentaerythritol tetranitrate C(CH₂ONO₂)₄."** *Physical Review B*, 69, 035116 (2004).
- Mazevet, S.; Challacombe, M.; Kowalski, P. M.; Saumon, D. **"He conductivity in cool white dwarf atmospheres."** *Astrophysics and Space Science*, 307, 273 (2007).

---

*Sources consulted: the FreeON/MondoSCF project homepage (nongnu.org/freeon), the project's own Publications page, the FreeON GitHub repository (github.com/FreeON/freeon), Wikipedia's FreeON article, and review literature on O(N) electronic structure methods (e.g., "O(N) methods in electronic structure calculations," arXiv:1108.5976, and "Challenges in Large Scale Quantum Mechanical Calculations," arXiv:1609.00252).*

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the FreeON 	Open-source, linear-scaling quantum chemistry/DFT code designed for very large molecular and periodic systems. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
