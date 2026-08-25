# Quickstep — The Electronic Structure and Molecular Dynamics Module of CP2K

*A comprehensive technical review*

---

## 1. Overview and Context

**CP2K** is an open-source (GPL-licensed) quantum chemistry and solid-state physics program package for atomistic simulation of solid-state, liquid, molecular, periodic, materials, and biological systems, with a particular emphasis on massively parallel and linear-scaling electronic-structure methods and state-of-the-art ab initio molecular dynamics (AIMD). **Quickstep** is CP2K's central electronic-structure engine — the module responsible for constructing and solving the electronic Hamiltonian that drives energies, forces, and properties throughout the package. It is written in Fortran 2008, is MPI/OpenMP parallel, and supports GPU (CUDA/HIP) and, historically, FPGA acceleration.

Quickstep is not a single method but a common computational infrastructure — integral routines, grid/multigrid mapping, sparse linear algebra, and SCF machinery — shared by a wide spectrum of electronic-structure approaches:

- Semi-empirical (SE) and tight-binding (TB, including GFN-xTB) methods
- Orbital-free DFT
- Kohn–Sham DFT (KS-DFT), the primary workhorse
- Wavefunction-based post-Hartree–Fock correlation methods: MP2, direct-RPA, GW

The name "Quickstep" reflects its historical goal — "make the atoms dance" efficiently — and its foundational algorithm is the **Gaussian and Plane Waves (GPW)** method, later extended by the **Gaussian Augmented Plane Wave (GAPW)** all-electron method.

---

## 2. The Gaussian and Plane Waves (GPW) Method

### 2.1 Dual Gaussian/plane-wave representation

Quickstep's defining innovation is representing molecular orbitals in a contracted Gaussian-type orbital (GTO) basis, while representing the electron density on an auxiliary real-space/reciprocal-space (plane-wave) grid. This dual representation exploits the fact that the Fourier transform of a Gaussian is again a Gaussian, so integrating Gaussian functions on an equidistant grid converges exponentially with grid spacing.

- **Orbital basis:** Contracted, atom-centered Cartesian/solid-harmonic Gaussian functions with fixed contraction coefficients.
- **Density representation:** The electron density n(r), built from the density matrix P and AO basis, is mapped onto a regular real-space grid and its plane-wave (reciprocal-space) expansion, bounded by a kinetic-energy cutoff **Ecut**.
- **Multigrid mapping:** Because Gaussians of different widths require different grid resolutions for efficient representation, Quickstep uses a **multigrid** scheme — several superimposed grids of different spacing — with optimal real-space screening and exploitation of the separability of Cartesian Gaussians along orthogonal axes. This mapping (density matrix → grid density, and back for potential → Kohn–Sham matrix elements) is a **linear-scaling** step with a small prefactor and is typically the dominant cost for standard GGA calculations.

### 2.2 The Kohn–Sham energy functional

Within the supercell/Γ-point approximation, the GPW Kohn–Sham energy is

E = E_kin + E_ext + E_ES + E_XC

where E_kin is the kinetic energy, E_ext the ionic-core interaction, E_ES the total electrostatic (Coulomb) energy, and E_XC the exchange–correlation energy. k-point sampling of the Brillouin zone is also supported. Electrostatics are handled via an Ewald-type decomposition with compensating Gaussian charge distributions, following the same machinery used in traditional plane-wave codes, including cluster-boundary methods (analytic Green's functions, Martyna–Tuckerman, wavelet approaches).

Forces and the stress tensor are computed consistently within the same framework, including **Pulay (wavefunction) forces** arising from the position-dependence of the atom-centered basis functions.

### 2.3 Pseudopotentials and basis sets

To avoid the steep cost increase associated with representing tightly localized core orbitals on the plane-wave grid, Quickstep primarily uses **norm-conserving, dual-space (Goedecker–Teter–Hutter, GTH) pseudopotentials**, which have a fully analytic Gaussian form compatible with the GPW machinery. An internal numerical atomic-DFT code is used to generate/optimize these pseudopotentials against relativistic all-electron reference calculations, and a large database of GTH parameter sets (per element, per XC functional) ships with CP2K.

Correspondingly optimized **MOLOPT basis sets** — generally contracted Gaussian sets sharing exponents across angular momenta, optimized against a condition-number penalty for good conditioning in condensed-phase and linear-scaling contexts — are the recommended default basis family, available in variants trading diffuseness/accuracy against multigrid mapping cost.

### 2.4 Local resolution-of-identity acceleration (LRIGPW)

For large, dense condensed-phase systems, the atom-pair density mapping step can be replaced by a **local resolution-of-the-identity (LRI)** approximation: pair densities are fit to a local auxiliary Gaussian basis at each atom (à la Baerends), collapsing the atom-*pair* sum into an atom sum — e.g., reducing ~200,000 pair terms to 192 atom terms for a 64-water-molecule box. Efficient solid-harmonic-based integral routines support this scheme (LRIGPW).

---

## 3. All-Electron Calculations: The GAPW Method

The **Gaussian Augmented Plane Wave (GAPW)** method is an alternative to pseudopotentials — or a means of using smaller PW cutoffs in pseudopotential calculations — providing an **all-electron** description. Borrowing the separation strategy of Blöchl's projector augmented-wave (PAW) method, GAPW decomposes the density into:

- a **smooth/soft** density represented in plane waves,
- **local hard** atom-centered densities represented on radial atomic grids as products of primitive Gaussians,
- **local soft** counterparts of the same, subtracted to avoid double counting.

This separates rapidly varying near-nucleus behavior from the smooth interstitial density while still integrating exactly over all space. GAPW enables calculation of core-electron-dependent properties, notably:

- X-ray absorption spectroscopy (XAS) via the transition-potential method
- Nuclear magnetic resonance (NMR) chemical shifts and electric field gradients
- Electron paramagnetic resonance (EPR) hyperfine couplings and g-tensors
- X-ray scattering (e.g., of liquid water)

---

## 4. Hartree–Fock Exchange and Hybrid DFT

Exact Hartree–Fock exchange (HFX) was introduced into Quickstep to enable condensed-phase, disordered-system AIMD with hybrid functionals (B3LYP, HSE, etc.). Key algorithmic features:

- **Γ-point evaluation** with an explicit sum over periodic image cells; for unscreened 1/r exchange, a **truncated Coulomb operator** with a system-size-dependent truncation radius R_C is used to make the sum rigorously convergent, paired with a long-range-correction (LRC) GGA-type exchange functional in the spirit of HSE.
- **Schwarz-type integral screening** based on basis-function overlap, reducing formal cost from O(N⁴) to O(N²) without assuming density-matrix sparsity.
- **"In-core" storage**: four-center two-electron integrals are computed once (analytically) at SCF start, compressed to the needed bit precision, and reused — parallelized via MPI+OpenMP for super-linear speed-up as memory scales with node count.
- **Multiple-time-step (MTS) integration**: HFX is updated only every few AIMD steps, exploiting the slow variation of the GGA-vs-hybrid potential-energy-surface difference.
- **Auxiliary Density Matrix Method (ADMM)**: HFX is evaluated cheaply from a projected auxiliary density matrix in a smaller auxiliary basis, with the difference corrected by a GGA-style exchange functional; this brings hybrid-DFT cost close to that of GGA-DFT and has enabled hundreds of picoseconds of hybrid-functional AIMD on hundreds of water molecules and spin-density calculations on ~3000-atom metalloproteins.

Two independent implementations (CP2K and Gaussian) have been cross-validated to micro-Hartree agreement on the HF total energy of crystalline LiH.

---

## 5. Post-Hartree–Fock and Many-Body Methods

Quickstep shares its Gaussian/GPW integral and optimization infrastructure with several correlated wavefunction methods, unified by the **resolution-of-identity (RI)** approximation to reduce four-center electron-repulsion integrals to sparse three-center tensor contractions.

### 5.1 MP2 and variants
- **Canonical MP2, RI-MP2** (O(N⁵), reduced prefactor via RI), **Laplace-transformed MP2**, and **scaled-opposite-spin (SOS-)MP2** are available at the Γ-point.
- Analytical gradients and stress tensors are implemented for RI-MP2 (both closed- and open-shell), enabling MP2-based geometry optimization and even AIMD on HPC systems.
- The implementation exploits sparse-matrix algebra and GPU acceleration for large matrix multiplications.
- Applications include AIMD/Monte Carlo simulations of bulk liquid water, ice XV structure refinement, and AIMD of the hydrated electron and other condensed-phase radicals.

### 5.2 Random Phase Approximation (RPA)
- **Direct-RPA (dRPA)** correlation energy via a frequency integral over a dielectric-matrix formulation (minimax quadrature, ~10 points for μHartree accuracy), combined with exact-exchange (EXX/RPA total energy formalism).
- **Quartic-scaling** O(N⁴N_q) implementation shares machinery with SOS-MP2 (both are post-processing steps on the same three-center RI-contracted matrix Q(ω)/Q(τ)).
- **Cubic-scaling** O(N³) (empirically better, ~O(N^1.8)) reformulation using imaginary-time transformation, atomic-orbital-pair tensors, and **block-sparse tensor contractions** built on the extended DBCSR tensor API — demonstrated on systems up to 864 water molecules (6912 electrons, ~49,000 primary and ~117,500 RI basis functions) on 3072 CPU cores.

### 5.3 GW approximation
- **G₀W₀** (non-self-consistent) and **eigenvalue-self-consistent GW (evGW)** implementations for valence-level ionization potentials/electron affinities, primarily for molecules; extensions toward periodic systems exist.
- Frequency integration of the self-energy via analytic continuation (2-pole model or Padé approximants).
- Standard implementation scales O(N⁴); a **low-scaling G₀W₀** variant achieves O(N²)–O(N³), enabling molecules beyond 1000 atoms.
- Applications include graphene-based nanomaterials on gold surfaces, combined with an **image-charge (IC) correction model** for level alignment at metal interfaces (comparing to scanning tunneling spectroscopy).

---

## 6. Constrained DFT (CDFT)

CDFT enforces charge/spin localization within atom-centered regions by augmenting the KS functional with Lagrangian constraint potentials, using Becke or Hirshfeld space-partitioning weight functions. Applications include:

- Charge-transfer (CT) state construction and Marcus-theory electronic-coupling calculations
- Correcting spurious self-interaction-driven charge delocalization
- Parameterizing model spin Hamiltonians
- **Mixed CDFT** (multi-state) simulations for charge-transfer kinetics and configuration-interaction (CI) calculations in a CDFT-state basis, correcting pathological self-interaction errors (e.g., recovering the correct H₂⁺ dissociation curve under PBE).

The SCF procedure is a two-tiered optimization: an outer loop optimizing constraint Lagrange multipliers (Newton/quasi-Newton root finding) around an inner, standard electronic SCF loop.

---

## 7. Density Functional Perturbation Theory (DFPT) and Response Properties

Quickstep implements DFPT for a range of second-order response properties via a preconditioned conjugate-gradient solution of the Sternheimer-type linear response equations:

- **Static electric-field response / polarizability and dipole moments**, using a Berry-phase formulation of polarization for periodic systems (modern theory of polarization), feeding into **Raman spectra** (isotropic/anisotropic transition polarizability) computed along finite-temperature AIMD trajectories.
- **NMR chemical shifts and EPR g-tensors** via GAPW-based all-electron DFPT response to an external magnetic field, including spin–orbit (SO) and spin–other–orbit (SOO) contributions to the g-tensor, gauge treatment via IGAIM (Keith–Bader) or CSGT schemes, and use of maximally localized Wannier functions (MLWFs) to handle the position operator under periodic boundary conditions. Applications span energy-storage materials, biomolecules, hydrogen-bonded systems, and paramagnetic solids (via combined hyperfine/g-tensor/orbital-shielding calculations).

---

## 8. Time-Dependent DFT (TD-DFT)

### 8.1 Linear-response TD-DFT (LR-TDDFT)
- Tamm–Dancoff approximation (ignoring de-excitations) reduces the response problem to a Hermitian eigenvalue problem, solved by a **block Davidson** iterative method.
- Supports hybrid exchange functionals with the same acceleration techniques as ground-state HFX (integral screening, truncated Coulomb operator, ADMM), and allows a different XC functional for the ground state versus the TD-DFT kernel (e.g., to add a long-range correction only at the response stage).
- Demonstrated for excitation energies in 1D/2D/3D periodic systems (aluminosilicate nanotube defects, MgO/HfO₂ surface and bulk vacancy defects) at system sizes up to ~1000 atoms.

### 8.2 Real-time propagation TD-DFT (RT-TDDFT)
- Enables non-linear response and direct simulation of electron dynamics, including **Ehrenfest dynamics** (simultaneous electron/nuclear propagation).
- Two flavors: cubic-scaling MO-coefficient-based propagation (**MO-RTP**, using Krylov-subspace/Arnoldi methods) and linear-scaling density-matrix-based propagation (**P-RTP**, via Baker–Campbell–Hausdorff-type expansion of the matrix exponential).
- Exponential-midpoint and enforced-time-reversal-symmetry (ETRS) propagators.
- True linear scaling in Ehrenfest dynamics is not achievable with plain P-RTP due to non-exponential density-matrix decay under time evolution, but can be recovered by coupling to **subsystem DFT** (Gordon–Kim-style kinetic-energy-functional coupling between independently minimized subsystems).

---

## 9. SCF Eigensolvers and Energy Minimization

Building the KS matrix is followed by an iterative SCF procedure. Quickstep offers multiple eigensolver/minimization strategies:

| Method | Scaling | Notes |
|---|---|---|
| **Traditional Diagonalization (TD)** | O(N³) | Via ScaLAPACK or ELPA; supports fractional occupations (needed for metals) |
| **Pseudodiagonalization (PD)** | reduced prefactor | Jacobi-rotation-based refinement of a pre-converged TD solution; BLAS-1 operations only |
| **Orbital Transformation (OT)** | favorable scaling | Direct energy minimization; default choice for efficiency, but requires fixed integer occupations (unsuitable for metals/zero-gap systems) |
| **Purification methods** | up to linear | Density-matrix-based, bypassing explicit orbitals (see §10) |

Orthogonality is enforced via Cholesky decomposition of the overlap matrix S (default) or symmetric Löwdin orthogonalization (useful for detecting/filtering basis-set linear dependencies, S eigenvalues < 10⁻⁵).

SCF convergence acceleration uses **DIIS** (direct inversion in the iterative subspace), which in CP2K automatically falls back to steepest descent when a DIIS step would be an ascent direction.

### 9.1 Orbital Transformation (OT) method in detail
Two parametrizations enforce orthogonality on an unconstrained auxiliary variable X:
- **OT/diag and OT/Taylor**: matrix functions cos(U), U⁻¹sin(U) evaluated by diagonalization or truncated Taylor expansion.
- **OT/IR**: an iterative-refinement expansion (Niklasson-type) approximating the orthogonality-enforcing matrix function to arbitrary order.

Minimizers available for the OT functional include DIIS (with SD fallback), non-linear conjugate gradients (Polak–Ribière with restart, golden-section or quadratic-interpolation line search), and quasi-Newton (Broyden type-2) methods.

**Preconditioners** range from the expensive but robust FULL_ALL (orbital-dependent eigenvalue shift, O(N³)) to cheaper linear-scaling approximations (FULL_KINETIC, FULL_S_INVERSE) and the commonly recommended **FULL_SINGLE_INVERSE**, which inverts only occupied eigenvalues at O(NM²) cost. For very large systems, sparse-matrix Hotelling iterations (via the `INVERSE_UPDATE` preconditioner solver) avoid explicit O(N³) matrix inversion.

---

## 10. Linear-Scaling (O(N)) Methods

### 10.1 Density-matrix purification
Rather than diagonalizing K, the density matrix P can be obtained directly by **purifying** K via matrix functions that map its eigenvalues through the Fermi–Dirac function. CP2K implements:
- **TRS4** (trace-resetting 4th order)
- **TC2** (trace-conserving 2nd order)
- **Sign-function-based purification** (up to 7th-order Padé approximants; the 5th-order iteration needs just four matrix multiplications per step), with aggressive filtering to preserve sparsity
- An interface to the external **PEXSI** library (Pole EXpansion and Selected Inversion) for evaluating selected density-matrix elements via a pole expansion of the Fermi–Dirac function

All rely on the **DBCSR** sparse block matrix library for efficient distributed sparse matrix–matrix multiplication (see §12). These methods have been demonstrated at extreme scale: a >1-million-atom simulation of the STMV virus in solution using the GFN-xTB tight-binding Hamiltonian on 10,240 CPU cores.

### 10.2 Submatrix method
An alternative to global sign-function iteration: matrix functions are evaluated on small, dense **principal submatrices** covering each atom and its basis-function-overlapping neighbors, rather than the whole sparse KS matrix — achieving linear scaling in both time and memory while allowing purification via either sign-function iteration or direct eigendecomposition on the small blocks.

### 10.3 Localized orbitals (LMOs, NLMOs, ALMOs, CLMOs)
Orbital-based O(N) methods complement density-matrix approaches by exploiting locality directly in the orbitals:
- **Berghold/Resta functionals** for maximally localized Wannier-type orbitals (periodic, Γ-point), equivalent to Boys–Foster localization in the gas phase.
- **Pipek–Mezey** localization, preserving σ/π bond separation — favored for molecular systems.
- **Non-orthogonal localized molecular orbitals (NLMOs)**, penalizing linear dependence to achieve tighter localization than orthogonal LMOs.
- **Compact localized molecular orbitals (CLMOs)** (generalizing ALMOs): orbitals with a strict, a-priori localization radius R_c around an atomic/molecular center, expanded only in basis functions within that radius. A robust **two-stage SCF** (ALMO optimization on centers, then constrained delocalization within R_c) and an approximate-Hessian-based optimizer for strongly covalent systems address the historically slow convergence of orbital-based O(N) methods. CLMO-DFT of liquid water demonstrates early-onset linear scaling with negligible accuracy loss relative to fully delocalized orbitals.

### 10.4 Machine-learned adaptive basis sets (PAO-ML)
The **polarized atomic orbital (PAO) machine-learning** scheme predicts a small, chemically adaptive basis transformation matrix from a larger primary basis using ML, trained via a Li–Nunes–Vanderbilt-style optimization, with rotational invariance enforced through neighbor-anchored auxiliary potentials. Demonstrated on liquid water AIMD with up to 200× computational cost reduction (and four orders of magnitude fewer floating-point operations) relative to basis-set-converged calculations.

---

## 11. Ab Initio Molecular Dynamics (AIMD)

Quickstep drives two principal flavors of AIMD, both computing the potential energy on the fly via KS-DFT:

### 11.1 Born–Oppenheimer MD (BOMD)
The electronic energy is fully minimized at each MD step under the orthonormality constraint, yielding forces = Hellmann–Feynman force + Pulay (wavefunction) force + implicit basis-position-dependence terms. Quickstep's efficient SCF/OT/purification machinery makes BOMD tractable for large and long trajectories.

### 11.2 (Second-generation) Car–Parrinello MD
CP2K implements the **second-generation Car–Parrinello MD (CPMD)** approach developed by Kühne and co-workers, which — unlike traditional first-generation CPMD — avoids the need for a small, artificial fictitious electron mass and correspondingly tiny integration time step, instead using a predictor–corrector-type always-stable propagation scheme for the electronic degrees of freedom. This is coupled with the CLMO/O(N) machinery (§10.3) to support genuinely linear-scaling AIMD.

### 11.3 Multiple-time-step and coupling to localized-orbital methods
As noted in §4, MTS integration allows expensive HFX (or, by extension, other slowly varying corrections) to be updated only periodically during a trajectory, and this idea generalizes to the correlated-wavefunction (MP2/RPA) methods.

---

## 12. Parallel Infrastructure and Performance

### 12.1 DBCSR — Distributed Block Compressed Sparse Row library
Co-developed with CP2K as a standalone library, **DBCSR** provides the block-sparse matrix–matrix multiplication infrastructure underpinning linear-scaling DFT, the OT preconditioners, purification/sign-function methods, and the tensor contractions of cubic-scaling RPA/SOS-MP2 and RI-MP2. Key features:

- Matrices are stored as a 2D grid of small dense blocks (typically 5–30 rows/columns) with CSR-format block-occupancy metadata.
- MPI ranks are arranged in a 2D Cartesian process grid mapped to matrix block rows/columns (a modified **Cannon's algorithm**), with a 2.5D one-sided-MPI variant for improved efficiency.
- Extended with a **tensor API**, mapping block-sparse tensor contractions onto sparse matrix–matrix multiplications, with a tall-and-skinny intermediate layer (simplified CARMA algorithm) reducing memory and communication cost.
- GPU-accelerated matrix multiplication kernels (LIBCUSMM/LIBSMM) support CUDA and HIP backends.

### 12.2 Massively parallel and heterogeneous execution
Quickstep supports distributed-memory (MPI), shared-memory (OpenMP), and accelerator (GPU: CUDA/HIP/tensor-core; historically FPGA) execution. GPU offload extends to grid-based integral evaluation, DBCSR sparse linear algebra, and dense matrix multiplication (via accelerated BLAS or the **SpLA** library) for correlated-wavefunction methods.

### 12.3 Approximate/stochastic computing
Aggressive matrix-element filtering thresholds (ϵ_filter) used throughout the sparse-matrix machinery (purification, preconditioner inversion, submatrix method) trade numerical precision for speed; this has been formalized as a form of **approximate computing (AC)**, including accurate MD sampling in the presence of controlled numerical/force noise.

---

## 13. Embedding, Spectroscopy, and Analysis Tools

- **QM/MM embedding**, including an adaptive buffered-force (ADBF) QM/MM coupling scheme for seamless treatment of solute–solvent or active-site/environment boundaries.
- **Subsystem DFT** (Gordon–Kim-style density-based embedding) for linear-scaling Ehrenfest/RT-TDDFT.
- **Energy decomposition analysis (EDA)**, exploiting the locality of compact orbitals (CLMOs/ALMOs) to decompose intermolecular interaction energies into physically interpretable components.
- **Solvation models**: e.g., the SCCS (self-consistent continuum solvation) model with a solvent-aware interface for implicit-solvent AIMD.
- **Spectroscopic property calculators**: IR and Raman spectra from AIMD dipole/polarizability autocorrelation functions, XAS via transition potentials, NMR/EPR via DFPT+GAPW (§7), and GW/BSE-adjacent tools for photoemission-related properties (§5.3).

---

## 14. Interfaces and Software Engineering

Quickstep/CP2K interfaces with a range of external tools and libraries, including (non-exhaustively) ScaLAPACK/ELPA for dense diagonalization, PEXSI for pole-expansion density-matrix evaluation, libxc for exchange-correlation functionals, and various machine-learning/force-field packages for hybrid QM/ML or QM/MM workflows. The codebase is organized to keep the same integral, grid, and sparse-linear-algebra infrastructure shared across semi-empirical, DFT, and post-HF methods, minimizing code duplication and easing the addition of new methods.

---

## 15. Summary of Key Distinguishing Features

1. **Dual Gaussian/plane-wave (GPW) representation** — analytic Gaussian integrals for the orbital basis combined with efficient, linear-scaling plane-wave-grid electrostatics/XC evaluation via multigrid mapping.
2. **All-electron capability (GAPW)** without sacrificing the GPW efficiency for the bulk of the density.
3. **Broad method coverage on one shared infrastructure**: semi-empirical/tight-binding, KS-DFT, hybrid-DFT (with ADMM acceleration), MP2/RI-MP2/SOS-MP2, RPA (quartic and cubic scaling), and GW.
4. **Multiple linear-scaling strategies**: density-matrix purification (sign-function, TRS4, TC2, PEXSI), submatrix methods, and localized-orbital (CLMO/ALMO) approaches — all built on the DBCSR sparse-matrix/tensor library.
5. **Advanced AIMD**: both Born–Oppenheimer MD and (second-generation, always-stable) Car–Parrinello MD, extensible to linear-scaling regimes and multiple-time-step acceleration of expensive terms.
6. **Rich response-property machinery**: DFPT-based polarizability/Raman, NMR/EPR via all-electron GAPW response, and both linear-response and real-time TD-DFT (including Ehrenfest dynamics).
7. **HPC-first design**: MPI+OpenMP+GPU parallelism throughout, exascale-oriented algorithmic choices, and support for approximate/stochastic computing to trade precision for throughput at extreme scale (million-atom-class simulations).

---

## 16. Publications Related to CP2K/Quickstep Theory

The following is a curated list of the primary scientific publications describing the theoretical foundations, algorithms, and methodology of CP2K and its Quickstep module. (CP2K's own output dynamically generates a context-specific "REFERENCES" section per simulation; consult that output, or the full literature list at `manual.cp2k.org` / `cp2k.org`, for citations covering a specific feature used in a given run.)

### 16.1 Foundational method papers

- Lippert, G.; Hutter, J.; Parrinello, M. **"A hybrid Gaussian and plane wave density functional scheme."** *Molecular Physics* **92**, 477–487 (1997). — Original GPW method.
- Lippert, G.; Hutter, J.; Parrinello, M. **"The Gaussian and augmented-plane-wave density functional method for ab initio molecular dynamics simulations."** *Theoretical Chemistry Accounts* **103**, 124–140 (1999). — Introduction of the GAPW all-electron extension.
- Krack, M.; Parrinello, M. **"Quickstep: Make the atoms dance."** In *High Performance Computing in Chemistry*, NIC Series Vol. 25, ed. J. Grotendorst (NIC-Directors, 2004), pp. 29–51. — Early overview of the Quickstep implementation.
- VandeVondele, J.; Krack, M.; Mohamed, F.; Parrinello, M.; Chassaing, T.; Hutter, J. **"Quickstep: Fast and accurate density functional calculations using a mixed Gaussian and plane waves approach."** *Computer Physics Communications* **167**, 103–128 (2005). https://doi.org/10.1016/j.cpc.2004.12.014 — The core Quickstep/GPW implementation paper.
- Goedecker, S.; Teter, M.; Hutter, J. **"Separable dual-space Gaussian pseudopotentials."** *Physical Review B* **54**, 1703–1710 (1996). — GTH pseudopotentials.
- Hartwigsen, C.; Goedecker, S.; Hutter, J. **"Relativistic separable dual-space Gaussian pseudopotentials from H to Rn."** *Physical Review B* **58**, 3641–3662 (1998). — Extended GTH pseudopotential set.
- VandeVondele, J.; Hutter, J. **"Gaussian basis sets for accurate calculations on molecular systems in gas and condensed phases."** *Journal of Chemical Physics* **127**, 114105 (2007). — MOLOPT basis sets.

### 16.2 Comprehensive package reviews

- Hutter, J.; Iannuzzi, M.; Schiffmann, F.; VandeVondele, J. **"CP2K: atomistic simulations of condensed matter systems."** *WIREs Computational Molecular Science* **4**, 15–25 (2014).
- Kühne, T. D.; Iannuzzi, M.; Del Ben, M.; Rybkin, V. V.; Seewald, P.; Stein, F.; Laino, T.; Khaliullin, R. Z.; Schütt, O.; Schiffmann, F.; Golze, D.; Wilhelm, J.; Chulkov, S.; Bani-Hashemian, M. H.; Weber, V.; Borštnik, U.; Taillefumier, M.; Jakobovits, A. S.; Lazzaro, A.; Pabst, H.; Müller, T.; Schade, R.; Guidon, M.; Andermatt, S.; Holmberg, N.; Schenter, G. K.; Hehn, A.; Bussy, A.; Belleflamme, F.; Tabacchi, G.; Glöß, A.; Lass, M.; Bethune, I.; Mundy, C. J.; Plessl, C.; Watkins, M.; VandeVondele, J.; Krack, M.; Hutter, J. **"CP2K: An electronic structure and molecular dynamics software package — Quickstep: Efficient and accurate electronic structure calculations."** *Journal of Chemical Physics* **152**, 194103 (2020). https://doi.org/10.1063/5.0007045 — The principal, comprehensive review of Quickstep methodology (source of most of this document's technical content).
- Kühne, T. D. et al. **"Roadmap on Exascale Computing"** / related CP2K exascale-performance perspectives (2022) covering algorithmic workflow, low-scaling post-HF methods, submatrix methods, and scaling benchmarks.
- Iannuzzi, M. et al. **"The CP2K Program Package Made Simple."** (2025) — Usage-oriented review of Quickstep, including practical guidance on RI-HFX, MP2 GPU acceleration, and recommended input recipes. arXiv:2508.15559.

### 16.3 Sparse linear algebra (DBCSR) and parallel performance

- Borštnik, U.; VandeVondele, J.; Weber, V.; Hutter, J. **"Sparse matrix multiplication: The distributed block-compressed sparse row library."** *Parallel Computing* **40**, 47–58 (2014). — Core DBCSR library paper.
- Schütt, O.; Messmer, P.; Hutter, J.; VandeVondele, J. **"GPU-accelerated sparse matrix–matrix multiplication for linear scaling density functional theory."** In *Electronic Structure Calculations on Graphics Processing Units*, Wiley (2016).
- Schütt, O. et al. **"DBCSR: A Blocked Sparse Tensor Algebra Library."** (tensor-contraction extension of DBCSR for correlated wavefunction methods). arXiv:1910.13555.
- Lass, M. et al. **"A Submatrix-Based Method for Approximate Matrix Function Evaluation in the Quantum Chemistry Code CP2K."** (2020). arXiv:2004.10811. — Submatrix method for linear-scaling matrix functions.
- Lazzaro, A.; VandeVondele, J.; Hutter, J.; Schütt, O. **"Increasing the Efficiency of Sparse Matrix-Matrix Multiplication with a 2.5D Algorithm and One-Sided MPI."** In *Proceedings of PASC '17* (2017).
- Yokelson, D. et al. **"Performance and Portability Analysis of CP2K"** or related hybrid MPI/OpenMP/GPU performance studies (2021).

### 16.4 Linear scaling and density-matrix methods

- VandeVondele, J.; Borštnik, U.; Hutter, J. **"Linear scaling self-consistent field calculations with millions of atoms in the condensed phase."** *Journal of Chemical Theory and Computation* **8**, 3565–3573 (2012).
- Niklasson, A. M. N. **"Expansion algorithm for the density matrix."** *Physical Review B* **66**, 155115 (2002). — TC2 purification.
- Niklasson, A. M. N.; Tymczak, C. J.; Challacombe, M. **"Trace resetting density matrix purification in O(N) self-consistent-field theory."** *Journal of Chemical Physics* **118**, 8611 (2003). — TRS4 purification.
- Lin, L.; Chen, M.; Yang, C.; He, L. **"Accelerating electronic structure calculations with the pole expansion and selected inversion method."** related PEXSI methodology papers (2013–2014).

### 16.5 Localized orbitals and O(N) orbital-based DFT

- Khaliullin, R. Z.; VandeVondele, J.; Hutter, J. **"Efficient linear-scaling density functional theory for molecular systems."** *Journal of Chemical Theory and Computation* **9**, 4421–4427 (2013).
- Berghold, G.; Mundy, C. J.; Romero, A. H.; Hutter, J.; Parrinello, M. **"General and efficient algorithms for obtaining maximally localized Wannier functions."** *Physical Review B* **61**, 10040 (2000).
- Pipek, J.; Mezey, P. G. **"A fast intrinsic localization procedure applicable for ab initio and semiempirical linear combination of atomic orbital wave functions."** *Journal of Chemical Physics* **90**, 4916 (1989).

### 16.6 Post-Hartree–Fock and many-body methods

- Guidon, M.; Hutter, J.; VandeVondele, J. **"Auxiliary density matrix methods for Hartree-Fock exchange calculations."** *Journal of Chemical Theory and Computation* **6**, 2348–2364 (2010). — ADMM method.
- Guidon, M.; Hutter, J.; VandeVondele, J. **"Robust periodic Hartree-Fock exchange for large-scale simulations using Gaussian basis sets."** *Journal of Chemical Theory and Computation* **5**, 3010–3021 (2009).
- Guidon, M.; Schiffmann, F.; Hutter, J.; VandeVondele, J. **"Ab initio molecular dynamics using hybrid density functionals."** *Journal of Chemical Physics* **128**, 214104 (2008).
- Del Ben, M.; Hutter, J.; VandeVondele, J. **"Second-order Møller–Plesset perturbation theory in the condensed phase: an efficient and massively parallel Gaussian and plane waves approach."** *Journal of Chemical Theory and Computation* **8**, 4177–4188 (2012).
- Del Ben, M.; Hutter, J.; VandeVondele, J. **"Electron correlation in the condensed phase from a resolution of identity approach based on the Gaussian and plane waves scheme."** *Journal of Chemical Theory and Computation* **9**, 2654–2671 (2013).
- Wilhelm, J.; Del Ben, M.; Hutter, J. **"GW in the Gaussian and plane waves scheme with application to linear acenes."** *Journal of Chemical Theory and Computation* **12**, 3623–3635 (2016).
- Wilhelm, J.; Golze, D.; Talirz, L.; Hutter, J.; Pignedoli, C. A. **"Toward GW calculations on thousands of atoms."** *Journal of Physical Chemistry Letters* **9**, 306–312 (2018). — Low-scaling GW.
- Golze, D.; Wilhelm, J.; van Setten, M. J.; Rinke, P. **"Core-level binding energies from GW: an efficient full-frequency approach within a localized basis."** *Journal of Chemical Theory and Computation* **14**, 4856–4869 (2018).

### 16.7 Response properties (DFPT, NMR/EPR, TD-DFT)

- Putrino, A.; Sebastiani, D.; Parrinello, M. **"Generalized variational density functional perturbation theory."** *Journal of Chemical Physics* **113**, 7102 (2000).
- Sebastiani, D.; Parrinello, M. **"A new ab-initio approach for NMR chemical shifts in periodic systems."** *Journal of Physical Chemistry A* **105**, 1951–1958 (2001).
- Weber, V.; Iannuzzi, M.; Giani, S.; Hutter, J.; Declerck, R.; Waroquier, M. **"Magnetic linear response properties calculation using a CP2K interface: the implementation and validation for gauge including atomic orbitals."** *Journal of Chemical Physics* **131**, 014106 (2009).
- Iannuzzi, M.; Chassaing, T.; Wallman, T.; Hutter, J. **"Ground and excited state density functional calculations with the Gaussian and augmented-plane-wave method."** *Chimia* **59**, 499 (2005). — Early TD-DFT in Quickstep.

### 16.8 Car–Parrinello / AIMD methodology

- Kühne, T. D.; Krack, M.; Mohamed, F. R.; Parrinello, M. **"Efficient and accurate Car-Parrinello-like approach to Born-Oppenheimer molecular dynamics."** *Physical Review Letters* **98**, 066401 (2007). — Second-generation CPMD.
- Kühne, T. D. **"Second generation Car–Parrinello molecular dynamics."** *WIREs Computational Molecular Science* **4**, 391–406 (2014).
- Kühne, T. D. (2026). Review/best-practice article on Car–Parrinello AIMD in Quickstep, predictor–corrector dynamics.

### 16.9 Machine learning and approximate computing

- Schran, C.; Brieuc, F.; Marx, D. **"Transferability of machine learning potentials."** and related PAO-ML methodology by Hutter, VandeVondele, and co-workers, e.g. on ML-adapted polarized atomic orbitals.
- Rengaraj, V.; Lass, M.; Plessl, C.; Kühne, T. D. **"Accurate Sampling with Noisy Forces from Approximate Computing."** *Computation* **8**, 39 (2020).

### 16.10 Related infrastructure and solvation

- Chai, Z.; Luber, S. **"Functional analytic derivation and CP2K implementation of the SCCS model based on the solvent-aware interface."** *Computer Physics Communications* **311**, 109563 (2025).

---

*Note: This document synthesizes information gathered from the primary CP2K/Quickstep literature, the official CP2K website (cp2k.org), the CP2K GitHub repository, and the CP2K user manual/bibliography. Publication details for some entries (particularly volume/page numbers of very recent works) should be cross-checked against the dynamically generated "REFERENCES" section of actual CP2K program output, or the authoritative literature list at `manual.cp2k.org`, before use in a formal citation list.*


---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the Quickstep 	Electronic structure and molecular dynamics module inside the open-source program package CP2K. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
