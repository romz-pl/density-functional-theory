# PARSEC: An Exhaustive Review

**Pseudopotential Algorithm for Real-Space Electronic Structure Calculations**

---

## 1. Overview

PARSEC (Pseudopotential Algorithm for Real-Space Electronic Calculations) is a real-space, grid-based electronic-structure code that solves the Kohn–Sham (KS) equations of density functional theory (DFT) by expanding spatial derivatives with **high-order finite differences (FD)** on a uniform Cartesian grid, rather than expanding wavefunctions in a plane-wave (Fourier) basis or in localized atomic-orbital basis sets. Electron wavefunctions are represented directly by their numerical values at discrete grid points, so there is no explicit basis set at all. PARSEC uses **norm-conserving pseudopotentials** (Troullier–Martins and related forms) to replace core electrons with a smooth effective potential for the valence electrons.

Created in 1994, PARSEC was the first practical code to solve the Kohn–Sham problem on a real-space grid using high-order finite differences, making it one of the founding codes of the real-space DFT family (a family that today also includes codes such as SPARC, DFT‑FE, RESCU, Octopus, and GPAW's finite-difference mode). It is written primarily in Fortran 95, is released under the GNU GPL (current public release: version 1.4.x), and is developed and maintained chiefly by James R. Chelikowsky's group (University of Texas at Austin, Center for Computational Materials / Oden Institute), with long-standing numerical-methods collaboration from Yousef Saad's group (University of Minnesota).

Project home: **real-space.org** · Source code: **github.com/PARSEC-real-space-code/PARSEC**

---

## 2. Historical Development

| Period | Milestone |
|---|---|
| Early 1990s | Motivating ideas developed by James R. Chelikowsky, Norm Troullier, and Yousef Saad at the University of Minnesota, seeking an alternative to plane-wave pseudopotential DFT that would scale better on massively parallel machines and handle non-periodic systems without artificial supercells. |
| 1994 | First formal presentation of the finite-difference-pseudopotential method (Chelikowsky, Troullier & Saad, *Phys. Rev. Lett.* **72**, 1240) and its higher-order extension (Chelikowsky, Troullier, Wu & Saad, *Phys. Rev. B* **50**, 11355), applied to diatomic molecules (Si₂, C₂, O₂, CO). |
| Mid-to-late 1990s | Development of efficient sparse iterative eigensolvers tailored to the method (Saad, Stathopoulos, Chelikowsky, Wu, Öğüt and others), enabling application to atomic clusters and early nanostructures. |
| 2000s | Extension to periodic and partially periodic (1D/2D/3D) boundary conditions; development of the **Chebyshev-filtered subspace iteration (CheFSI)** method (Zhou, Saad, Tiago & Chelikowsky) to bypass expensive explicit diagonalization at most self-consistency steps; addition of time-dependent DFT (TDDFT), GW/Bethe–Salpeter excited-state post-processing, spin-orbit coupling, and non-collinear magnetism. |
| 2006 | Comprehensive review/reference paper published in *Physica Status Solidi (b)* **243**, 1063, describing the formalism, numerics, and applications of the mature code — still one of the standard citations for PARSEC. |
| 2010s–2020s | Continued algorithmic modernization: polynomial/spectrum-slicing filtering, space-filling-curve data layouts for sparse matrix-vector products, finite-difference interpolation (FDI) to reduce "egg-box" grid errors, orbital-based force schemes, and demonstrations of electronic-structure calculations exceeding 100,000 atoms on leadership-class supercomputers (e.g., Frontera at TACC). |
| Present | PARSEC remains an actively maintained open-source (GPL-3.0) real-space DFT code, hosted on GitHub, and is featured as one of 14 codes surveyed in the 2022 community "Roadmap on Electronic Structure Codes in the Exascale Era." |

---

## 3. Core Methodology

### 3.1 Real-space finite-difference discretization

Rather than expanding the Kohn–Sham orbitals in plane waves or atomic orbitals, PARSEC samples all physical quantities — wavefunctions, electron density, potentials — on a uniform Cartesian grid covering the physical domain. The kinetic-energy (Laplacian) operator acting on the wavefunctions is approximated by a **high-order finite-difference stencil**, typically a centered difference formula of order 4–12, which controls the accuracy/cost trade-off of the discretization independent of the grid spacing itself. This is in contrast to standard second-order finite-difference schemes, and was one of the key innovations distinguishing PARSEC's approach from earlier, lower-accuracy real-space attempts.

For non-periodic (cluster/molecule/nanostructure) calculations, the domain is a finite region — most commonly a sphere — enclosing the physical system, with wavefunctions constrained to vanish outside a specified boundary radius. This "natural" confined-boundary-condition treatment avoids the need for artificial periodic supercells and vacuum padding that plane-wave codes require for finite systems, and it makes charged and defect systems especially natural and efficient to treat. PARSEC also supports full 3D periodic boundary conditions, as well as partially periodic systems (1D wires/nanotubes, 2D slabs/surfaces).

### 3.2 Pseudopotentials

PARSEC uses **norm-conserving pseudopotentials**, most commonly of the Troullier–Martins form, to remove chemically inert core electrons from the explicit calculation while preserving scattering properties of the valence electrons. Nonlocal (Kleinman–Bylander) separable projector forms are used for the angular-momentum-dependent parts of the pseudopotential. Because the pseudopotentials are evaluated on the same real-space grid as the wavefunctions, aliasing/positional ("egg-box") errors can arise when atoms move relative to the fixed grid; PARSEC addresses this with Fourier-filtering of the pseudopotentials (mask-function methods) and, more recently, with finite-difference interpolation (FDI) schemes that exploit the pseudopotential's high intrinsic resolution to reduce this error at low added cost.

### 3.3 Solving the Kohn–Sham eigenproblem

Discretization on the grid turns the Kohn–Sham equation into a large, **sparse** matrix eigenvalue problem (the finite-difference Hamiltonian is sparse because each grid point couples only to a limited stencil of neighbors). PARSEC exploits this sparsity in several complementary ways:

- **Direct sparse diagonalization** using iterative Krylov-subspace eigensolvers (e.g., Lanczos/Davidson-type methods), historically the default approach and still used to "seed" other schemes.
- **Chebyshev-filtered subspace iteration (CheFSI)**: rather than diagonalizing the Hamiltonian explicitly at every self-consistent-field (SCF) iteration, a low-degree Chebyshev polynomial filter is applied to the subspace from the previous iteration to damp components outside the occupied eigenspace. Explicit diagonalization is then needed only at the first SCF step, dramatically reducing cost for later iterations. This nonlinear subspace-iteration approach (Zhou, Saad, Tiago, Chelikowsky, 2006) is one of PARSEC's signature numerical contributions and has since been adopted or reproduced in other electronic-structure codes (e.g., DFT-FE).
- **Polynomial/spectrum-slicing filtering** (Liou, Yang & Chelikowsky, 2020): a hybrid scheme combining CheFSI with spectrum slicing, in which the eigenvalue spectrum is partitioned into slices that are filtered (band-pass) and solved somewhat independently, improving parallel scalability for very large problems.
- **Sparse-matrix data layout optimizations**, including blockwise Hilbert space-filling curves to improve locality and speed of sparse matrix–vector multiplication, which is the dominant kernel in the filtering-based solvers.

These developments have allowed PARSEC to converge self-consistent electronic structure for systems well beyond 100,000 atoms (e.g., a Si₁₀₇,₆₄₁H₉,₀₈₄ silicon nanocluster ~16 nm in diameter) on leadership machines such as TACC's Frontera, using thousands of MPI ranks.

### 3.4 Parallelization

PARSEC is parallelized with the MPI message-passing model. Parallelization is achieved primarily by domain decomposition of the real-space grid across processes, along with parallel dense/sparse linear algebra (e.g., ScaLAPACK) for the subspace-projected problems arising in CheFSI-type solvers. Because the finite-difference Hamiltonian is sparse and local, real-space methods like PARSEC avoid the global data communication associated with the fast Fourier transforms central to plane-wave codes, which is one of the method's principal scalability advantages on massively parallel and (more recently) GPU-unified-memory architectures.

---

## 4. Capabilities and Features

- **Ground-state DFT** within LDA and GGA exchange-correlation functionals.
- **Boundary conditions**: confined/non-periodic (finite droplet of grid points, natural for clusters, molecules, nanocrystals, and defects/charged systems) or periodic in 1, 2, or 3 dimensions.
- **Structural relaxation** (geometry optimization) and **simulated annealing**.
- **Born–Oppenheimer / Langevin molecular dynamics**.
- **Spin polarization**, **spin-orbit coupling**, and **non-collinear magnetism**.
- **Static polarizability** calculations (confined/non-periodic boundary conditions).
- **Excited-state workflows**: PARSEC's ground-state Kohn–Sham orbitals serve as input to post-processing excited-state methods used by the Chelikowsky group and collaborators, including:
  - **Time-dependent DFT (TDDFT / TDLDA)** for optical absorption spectra.
  - **GW approximation** for quasiparticle (ionization potential/electron affinity) energies.
  - **Bethe–Salpeter equation (BSE)** solutions for excitonic optical spectra, including efficient Lanczos-based full-frequency GW implementations designed to exploit the same real-space framework.
- **Massively parallel, sparse eigensolvers** as described above (CheFSI, polynomial/spectrum-slicing filtering).
- Demonstrated scalability to systems of **~100,000+ atoms / ~1,000,000 electrons** on modern HPC platforms.

PARSEC was specifically designed with **nanostructures** in mind — quantum dots, nanocrystals, clusters, surfaces, and point defects — where the natural confined-domain treatment, direct access to real-space wavefunctions/densities, and efficient handling of charged or symmetry-broken systems are particularly advantageous relative to supercell plane-wave approaches.

---

## 5. Comparison with Plane-Wave DFT

| Aspect | PARSEC (real-space FD) | Typical plane-wave DFT |
|---|---|---|
| Basis / representation | None; orbitals sampled on a real grid | Fourier/plane-wave expansion |
| Non-periodic systems | Natural: finite domain, no periodic images | Requires artificial supercell + vacuum |
| Charged / defect systems | Handled directly and efficiently | Needs charge-compensation schemes |
| Parallel communication | Local (sparse stencil operations); avoids global FFTs | Requires global FFTs (all-to-all communication) |
| Boundary conditions | Confined or periodic (1D/2D/3D), flexibly mixed | Inherently periodic |
| Convergence parameter | Grid spacing (and FD order) | Plane-wave cutoff energy |
| Known artifact | "Egg-box" grid-positional error (mitigated via filtering/FDI) | Pulay stress in variable-cell relaxation |

---

## 6. Key Publications

### 6.1 Foundational / Official Citation Papers
(These are the citations the PARSEC developers formally request users acknowledge.)

1. J. R. Chelikowsky, N. Troullier, and Y. Saad, "Finite-difference-pseudopotential method: electronic structure calculations without a basis," *Physical Review Letters* **72**, 1240 (1994). https://doi.org/10.1103/PhysRevLett.72.1240
2. J. R. Chelikowsky, N. Troullier, K. Wu, and Y. Saad, "Higher-order finite-difference pseudopotential method: An application to diatomic molecules," *Physical Review B* **50**, 11355 (1994). https://doi.org/10.1103/PhysRevB.50.11355
3. J. R. Chelikowsky, "The pseudopotential-density functional method applied to nanostructures," *Journal of Physics D: Applied Physics* **33**, R33 (2000). https://doi.org/10.1088/0022-3727/33/8/201
4. L. Kronik, A. Makmal, M. L. Tiago, M. M. G. Alemany, M. Jain, X. Huang, Y. Saad, and J. R. Chelikowsky, "PARSEC – the pseudopotential algorithm for real-space electronic structure calculations: recent advances and novel applications to nano-structures," *Physica Status Solidi (b)* **243**, 1063 (2006). https://doi.org/10.1002/pssb.200541463
5. Y. Saad, J. R. Chelikowsky, and S. M. Shontz, "Numerical methods for electronic structure calculations of materials," *SIAM Review* **52**, 3 (2010). https://doi.org/10.1137/060651653

### 6.2 Eigensolver / Numerical Linear Algebra Theory

6. Y. Saad, A. Stathopoulos, J. Chelikowsky, K. Wu, and S. Öğüt, "Solution of large eigenvalue problems in electronic structure calculations," *BIT Numerical Mathematics* **36**, 563 (1996). https://doi.org/10.1007/BF01731934
7. Y. Zhou, Y. Saad, M. L. Tiago, and J. R. Chelikowsky, "Self-consistent-field calculations using Chebyshev-filtered subspace iteration," *Journal of Computational Physics* **219**, 172 (2006). https://doi.org/10.1016/j.jcp.2006.03.017
8. Y. Zhou, Y. Saad, M. L. Tiago, and J. R. Chelikowsky, "Parallel self-consistent-field calculations via Chebyshev-filtered subspace acceleration," *Physical Review E* **74**, 066704 (2006). https://doi.org/10.1103/PhysRevE.74.066704
9. Y. Saad, Y. Zhou, C. Bekas, M. Tiago, and J. R. Chelikowsky, "Diagonalization methods in PARSEC," *Physica Status Solidi (b)* **243**, 2188 (2006). https://doi.org/10.1002/pssb.200666816
10. Y. Zhou, J. R. Chelikowsky, and Y. Saad, "Chebyshev-filtered subspace iteration method free of sparse diagonalization for solving the Kohn–Sham equation," *Journal of Computational Physics* **274**, 770 (2014). https://doi.org/10.1016/j.jcp.2014.06.056
11. G. Schofield, J. R. Chelikowsky, and Y. Saad, "A spectrum slicing method for the Kohn–Sham problem," *Computer Physics Communications* **183**, 497 (2012).
12. K.-H. Liou, A. Biller, L. Kronik, and J. R. Chelikowsky, "Space-filling curves for real-space electronic structure calculations," *Journal of Chemical Theory and Computation* **17**, 4039 (2021). https://doi.org/10.1021/acs.jctc.1c00237
13. K.-H. Liou, C. Yang, and J. R. Chelikowsky, "Scalable implementation of polynomial filtering for density functional theory calculation in PARSEC," *Computer Physics Communications* **254**, 107330 (2020). https://doi.org/10.1016/j.cpc.2020.107330

### 6.3 Boundary Conditions and Method Extensions

14. M. M. G. Alemany, M. Jain, L. Kronik, and J. R. Chelikowsky, "Real-space pseudopotential method for computing the electronic properties of periodic systems," *Physical Review B* **69**, 075101 (2004). https://doi.org/10.1103/PhysRevB.69.075101
15. A. Natan, A. Benjamini, D. Naveh, L. Kronik, M. L. Tiago, S. P. Beckman, and J. R. Chelikowsky, "Real-space pseudopotential method for first principles calculations of general periodic and partially periodic systems," *Physical Review B* **78**, 075109 (2008). https://doi.org/10.1103/PhysRevB.78.075109
16. J. Han, M. L. Tiago, T.-L. Chan, and J. R. Chelikowsky, "Real space method for the electronic structure of one-dimensional periodic systems," *Journal of Chemical Physics* **129**, 144109 (2008). https://doi.org/10.1063/1.2988316

### 6.4 Force Precision and Grid ("Egg-Box") Error Reduction

17. D. Roller, A. M. Rappe, L. Kronik, and O. Hellman, "Finite Difference Interpolation for Reduction of Grid-Related Errors in Real-Space Pseudopotential Density Functional Theory," *Journal of Chemical Theory and Computation* **19**, 3889 (2023). https://doi.org/10.1021/acs.jctc.3c00217
18. D. Roller et al., "Improving the precision of forces in real-space pseudopotential density functional theory," *Journal of Chemical Physics* **161**, 074113 (2024). https://doi.org/10.1063/5.0219847

### 6.5 Large-Scale / Exascale-Era Performance

19. M. Dogan, K.-H. Liou, and J. R. Chelikowsky, "Solving the electronic structure problem for over 100 000 atoms in real space," *Physical Review Materials* **7**, L063001 (2023). https://doi.org/10.1103/PhysRevMaterials.7.L063001
20. M. Dogan, K.-H. Liou, and J. R. Chelikowsky, "Real-space solution to the electronic structure problem for nearly a million electrons" (companion large-scale-performance study).
21. K.-H. Liou, M. Dogan, and J. R. Chelikowsky, "PARSEC: Real-Space Pseudopotential Density Functional Theory Code," Section 9 in *"Roadmap on Electronic Structure Codes in the Exascale Era,"* arXiv:2209.12747 (2022).

### 6.6 Excited-State (TDDFT / GW / BSE) Applications Built on PARSEC

22. M. L. Tiago and J. R. Chelikowsky, "Optical excitations in organic molecules, clusters, and defects studied by first-principles Green's function methods," *Physical Review B* **73**, 205334 (2006). https://doi.org/10.1103/PhysRevB.73.205334
23. M. L. Tiago, J. C. Idrobo, S. Öğüt, J. Jellinek, and J. R. Chelikowsky, "Electronic and optical excitations in clusters: Comparison of density-functional and many-body theories," *Physical Review B* **79**, 155419 (2009). https://doi.org/10.1103/PhysRevB.79.155419
24. M. Lopez del Puerto, M. L. Tiago, and J. R. Chelikowsky, "Excitonic effects and optical properties of passivated CdSe clusters," *Physical Review Letters* **97**, 096401 (2006).
25. L. Hung, F. H. da Jornada, J. Souto-Casares, J. R. Chelikowsky, S. G. Louie, and S. Öğüt, "Excitation spectra of aromatic molecules within a real-space GW-BSE formalism: Role of self-consistency and vertex corrections," *Physical Review B* **94**, 085125 (2016).

### 6.7 Early Cluster/Application Papers Establishing the Method

26. J. R. Chelikowsky, N. Troullier, X. Jing, D. Dean, N. Binggeli, K. Wu, and Y. Saad, "Algorithms for the structural properties of clusters," *Computer Physics Communications* **85**, 325 (1995).
27. X. Jing, N. Troullier, J. R. Chelikowsky, K. Wu, and Y. Saad, "Vibrational modes of silicon nanostructures," *Solid State Communications* **96**, 231 (1995).
28. A. Stathopoulos, S. Öğüt, Y. Saad, J. Chelikowsky, and H. Kim, "Parallel methods and tools for predicting material properties," *Computing in Science & Engineering* **2**, 19 (2000). https://doi.org/10.1109/5992.852388

---

## 7. Summary Assessment

PARSEC occupies a distinctive and historically important niche in the electronic-structure software ecosystem: it was the pioneering demonstration that high-order finite-difference discretization on a uniform real-space grid could match plane-wave accuracy while offering more natural treatment of non-periodic, charged, and defect-containing systems, and while avoiding the global-communication bottleneck (FFTs) that limits plane-wave codes' scalability on massively parallel machines. Three decades of subsequent development — particularly the Chebyshev-filtered subspace iteration family of eigensolvers, space-filling-curve data layouts, and finite-difference interpolation schemes for reducing grid-positional ("egg-box") force and energy errors — have kept PARSEC numerically competitive and allowed it to be applied to systems exceeding 100,000 atoms on leadership supercomputers. Its tight coupling to real-space GW/BSE and TDDFT post-processing methods has also made it a workhorse for excited-state and optical-property studies of nanostructures (quantum dots, passivated clusters, nanocrystals) where quantum confinement and finite-system boundary conditions are central to the physics. The code remains open source (GPL) and under active development, and is recognized in the community's 2022 exascale-era roadmap as one of the established real-space DFT codes charting a path toward exascale computing.

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the PARSEC 	Real-space pseudopotential DFT code using finite differences on a grid rather than plane waves. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
