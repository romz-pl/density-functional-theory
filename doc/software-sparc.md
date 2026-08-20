# SPARC: Simulation Package for Ab-initio Real-space Calculations — An Exhaustive Review

## 1. Overview

**SPARC** (Simulation Package for Ab-initio Real-space Calculations) is an open-source Kohn–Sham Density Functional Theory (DFT) code built around a **real-space finite-difference** discretization rather than the plane-wave/Fourier basis used by most legacy electronic-structure packages (e.g., VASP, Quantum ESPRESSO, ABINIT, CASTEP). It is developed principally by **Phanish Suryanarayana's group at the Georgia Institute of Technology**, in collaboration with **John E. Pask at Lawrence Livermore National Laboratory (LLNL)**, with contributions from a wider community of students and collaborators (Swarnava Ghosh, Qimen Xu, Abhiraj Sharma, Xin Jing, and others).

SPARC is written primarily in **C** (≈90% of the codebase) with a **Python** front end/API layer, parallelized with **MPI**, and released under the **GPL-3.0** license. It is hosted at [github.com/SPARC-X/SPARC](https://github.com/SPARC-X/SPARC), with a companion prototyping code, **M-SPARC** (a MATLAB implementation sharing the same algorithms/input-output structure), and a Python interoperability layer, **SPARC-X-API**, built on the Atomic Simulation Environment (ASE).

The project's stated goal is to combine the **accuracy** of established plane-wave methods with **superior computational efficiency and scalability** on modern high-performance computing (HPC) architectures, including large CPU clusters and, more recently, GPU-accelerated systems.

---

## 2. Historical Development and Milestones

| Stage | Milestone |
|---|---|
| ~2015–2016 | First formulation papers: real-space finite-difference DFT framework for **isolated clusters** (Ghosh & Suryanarayana, 2017) and **extended systems** (Ghosh & Suryanarayana, 2017), establishing the electrostatics reformulation, Chebyshev-filtered SCF, and non-local force expressions. |
| 2016–2018 | Development of efficient real-space **linear solvers** (AAR — Alternating Anderson–Richardson) for the electrostatic (Poisson-like) problem, avoiding the need for FFTs. |
| 2019–2020 | Launch of **M-SPARC**, a MATLAB rapid-prototyping implementation mirroring SPARC's algorithms for methods development and small/moderate systems. |
| 2021 | **SPARC v1.0.0** released as a unified, publicly available C/C++ package (Xu, Sharma, Suryanarayana, Gavini, Pask et al., *SoftwareX*), consolidating isolated and extended-system capabilities, structural relaxation, and ab initio molecular dynamics (AIMD). |
| 2021–2022 | Introduction of the linear-scaling **O(N) Spectral Quadrature (SQ)** method, enabling simulations of up to ~1 million atoms; **Cyclix-DFT** for cyclic/helical symmetry-adapted calculations of 1D/2D nanomaterials. |
| 2023 | **GPU acceleration** of local/semilocal (LDA/GGA) DFT workflows (Sharma, Metere, Suryanarayana, Erlandson, Chow, Pask, *J. Chem. Phys.* 2023). |
| 2023 | **SPARC v2.0.0** released (Zhang, Sharma, Xu, Suryanarayana, Pask et al.), adding spin-orbit coupling, dispersion corrections, DFT+U, advanced XC functionals, real-space DFPT, and on-the-fly machine-learned force fields (MLFF). |
| 2024–2025 | GPU acceleration extended to **hybrid exchange–correlation functionals** (Jing, Sharma, Pask, Suryanarayana, 2025); **SPARC-X-API** Python interface released for ASE-based workflow integration; real-space **DFT+U** and Hubbard-corrected DFT formalisms published. |
| Ongoing | Development of Discrete Discontinuous Basis Projection (DDBP) for reduced-order large-scale calculations, RPA correlation energy, orbital-free DFT, and further exascale/GPU optimization. |

---

## 3. Theoretical and Numerical Foundations

### 3.1 Real-space finite-difference discretization
Unlike plane-wave codes, which expand Kohn–Sham orbitals in a Fourier basis (imposing intrinsic periodicity and relying on FFTs), SPARC discretizes all fields — orbitals, electron density, potentials — on a **uniform real-space grid**. Differential operators (the kinetic energy Laplacian, gradients) are approximated using **high-order (e.g., 6th–12th order) centered finite-difference stencils**, and integrals are evaluated via the **trapezoidal rule**. Convergence is controlled by a **single parameter**: the mesh spacing.

Advantages of this approach include:
- **Systematic, variational-like convergence** controlled by one parameter (mesh size), directly comparable to plane-wave cutoff energy convergence.
- **Natural handling of both Dirichlet (isolated) and Bloch-periodic (extended) boundary conditions** within the same framework, enabling treatment of finite, semi-infinite, charged, and bulk 3D systems without artificial periodic images or large vacuum paddings for cluster calculations.
- **Locality/sparsity**: finite-difference operators are sparse and involve only nearest-neighbor-type stencils, in contrast to the global, dense communication pattern of FFT-based plane-wave methods. This is the key property enabling efficient, near-linear parallel scalability.
- **Freedom from communication-intensive global transforms** (no FFTs required for the kinetic operator), which reduces inter-processor communication overhead on large processor counts — a major bottleneck for plane-wave codes at scale.
- **Natural framework for O(N) linear-scaling methods** because of the locality of the representation.

### 3.2 Electrostatics reformulation
A central theoretical contribution of SPARC is a **local reformulation of the electrostatic (Hartree + ionic) energy**. The standard electrostatic energy expression is inherently non-local (involving double integrals over the whole domain), which would scale poorly (formally O(N²)) if implemented directly in real space. SPARC recasts this into a **local variational problem** — solving a generalized Poisson-type equation with pseudocharge densities — that can be handled with the same sparse, local finite-difference machinery used elsewhere in the code, restoring efficient O(N) scaling for the electrostatics.

### 3.3 Chebyshev-filtered subspace iteration (CheFSI)
To solve the Kohn–Sham eigenvalue problem efficiently at each self-consistent field (SCF) step, SPARC employs the **Chebyshev polynomial-filtered subspace iteration** method: a filter based on Chebyshev polynomials is applied to an approximate subspace to damp unwanted (high-energy) spectral components, followed by Rayleigh–Ritz projection and orthonormalization. This avoids the need to fully diagonalize the Hamiltonian at every SCF step and is well suited to iterative, matrix-free, real-space implementations.

### 3.4 Non-local pseudopotential and force reformulation
SPARC uses **norm-conserving ONCV (Optimized Norm-Conserving Vanderbilt) pseudopotentials** (`psp8` format, generated via D. R. Hamann's ONCVPSP code), including support for **nonlinear core corrections (NLCC)**. The package derives and implements a reformulated expression for the **non-local component of atomic forces** that is consistent with the finite-difference energy expression, avoiding the "egg-box"/Pulay-type inconsistencies that can otherwise plague grid-based force calculations. Formulations for the **stress tensor** (for cell relaxation and pressure calculations) were subsequently derived and published, extending consistency to stress as well as energy and forces.

### 3.5 Linear solvers
For the electrostatic/Poisson-type problems, SPARC employs specialized real-space iterative solvers, notably the **Alternating Anderson–Richardson (AAR)** method, along with preconditioned iterative schemes tailored to the sparse finite-difference operators, avoiding reliance on FFT-based solvers.

### 3.6 Self-consistent field (SCF) acceleration
Convergence of the SCF cycle is accelerated using **Pulay/Anderson-type mixing schemes with real-space preconditioning**, with published work specifically addressing mixing and preconditioning strategies (including restarting schemes) for real-space finite-difference DFT.

---

## 4. Core Functionality and Feature Set

SPARC supports a broad range of production DFT capabilities:

**System types**
- Isolated systems (molecules, clusters)
- Extended/periodic systems (crystals, surfaces/slabs, wires) with Bloch-periodic boundary conditions and Brillouin-zone (k-point) sampling
- Charged systems and mixed (Dirichlet + periodic) boundary conditions
- Non-orthogonal unit cells

**Exchange–correlation and Hamiltonian treatment**
- Local (LDA), semilocal (GGA, e.g. PBE), and meta-GGA functionals
- Non-local (hybrid) exchange–correlation functionals, including GPU-accelerated hybrid calculations via a batched Kronecker-product linear solver
- Spin-polarized and spin-unpolarized calculations
- **Spin-orbit coupling (SOC)**
- **Noncollinear spin**
- **DFT+U** (Hubbard correction), including a real-space Hubbard-corrected DFT formalism with energy, force, and stress expressions
- Dispersion (van der Waals) corrections: **DFT-D3, vdW-DF1, vdW-DF2**

**Ground-state and dynamics capabilities**
- Ground-state total energy, atomic forces, and full stress tensor
- Structural relaxation: atomic positions and/or cell (volume/shape) relaxation
- Ab initio molecular dynamics (AIMD/QMD) with multiple ensembles: **NVE, NVT-Nosé-Hoover (NVTNH), NVK-Gaussian (NVKG), NPT-Nosé-Hoover (NPTNH), NPT-Nosé-Poincaré (NPTNP), and NPH**
- **On-the-fly machine-learned force field (MLFF)** generation for accelerated molecular dynamics, including extensions compatible with the O(N) SQ method and Cyclix symmetry
- Real-space **Density Functional Perturbation Theory (DFPT)** for response properties (e.g., phonons), including a cyclic/helical-symmetry-adapted phonon formalism
- Pseudopotential-based electronic structure via norm-conserving ONCV potentials (`psp8`), with a bundled table of **SPMS (soft and transferable) pseudopotentials**

**Advanced/large-scale methods**
- **O(N) linear-scaling Spectral Quadrature (SQ) method** — enables simulations with over a million atoms and has been applied extensively to warm/hot dense matter and finite-temperature (high-T) simulations
- **Cyclix-DFT**: symmetry-adapted DFT exploiting cyclic and/or helical symmetries, dramatically reducing the degrees of freedom needed to simulate 1D/2D nanomaterials (nanotubes, nanowires, twisted/deformed low-dimensional systems) under mechanical deformation
- Discrete Discontinuous Basis Projection (DDBP) method under development for reduced-order large-scale calculations
- Orbital-free DFT (in development)
- Correlation energy within the Random Phase Approximation (RPA) (in development)

**Interfaces and tooling**
- **M-SPARC**: MATLAB implementation used for algorithm prototyping, teaching, and small/moderate system sizes, sharing input/output conventions with the C/C++ code
- **SPARC-X-API**: a Python package built on the ASE (Atomic Simulation Environment) standard, providing a unified interface for running SPARC within broader computational workflows, including a JSON schema for input-parameter validation and an i-PI-derived socket communication layer for message passing between the C backend and Python front end
- Pre-compiled binaries distributed via `conda-forge` (`sparc-x` package) for x86_64/aarch64 Linux
- A companion **ATOM code** implementing a spectral scheme for isolated-atom electronic structure, used in generating/validating atomic orbitals (e.g., for DFT+U and orbital-based post-processing)

---

## 5. Parallel Implementation, Performance, and Scalability

### 5.1 Parallelization architecture
SPARC's parallel implementation is built on **MPI**, using a Cartesian-topology domain decomposition of the real-space grid, supplemented by additional parallelization layers. Later developments (e.g., for DFT+U) describe a **four-level parallelization hierarchy**, distributing processors across **spin, Brillouin-zone (k-point) wavevectors, electronic bands, and the spatial domain**, in that order — a hierarchy broadly analogous to that used in state-of-the-art plane-wave codes, but exploiting the additional locality of the finite-difference representation for the domain-decomposed pieces (kinetic operator, electrostatics).

Compilation depends on **BLAS/LAPACK** (or Intel MKL) and optionally **ScaLAPACK**; the code deliberately minimizes external library dependencies to ease portability and installation across HPC systems.

### 5.2 CPU performance and scaling
Across the benchmark studies (isolated clusters, extended systems, and the unified v1.0/v2.0 papers), SPARC consistently demonstrates:
- **Weak and strong parallel scaling** behavior comparable to well-optimized, mature plane-wave codes (e.g., ABINIT) for systems from hundreds to thousands of electrons, but with a **significantly reduced prefactor** (i.e., faster in absolute time even when the scaling curve shape is similar).
- **Order-of-magnitude (or greater) speedups** relative to state-of-the-art plane-wave codes in minimum time-to-solution, with the advantage **increasing** as the number of processors or the system size grows — attributed to the reduced communication footprint of the local finite-difference operators versus global FFTs.
- In regular production use, SPARC scales efficiently to **thousands of processors**, bringing time-to-solution down to about a minute for systems of O(500–1000) atoms, and a few seconds for O(100–500) atoms.
- Using the **O(N) SQ method**, SPARC has been scaled to system sizes exceeding **one million atoms**.
- Negligible energy drift is observed in molecular dynamics simulations, and forces/energies show negligible "egg-box" (grid-alignment) artifacts, both hallmarks of a numerically consistent finite-difference formulation.

### 5.3 GPU acceleration
Two major GPU-acceleration efforts have been published:

1. **Local/semilocal (LDA/GGA) GPU acceleration** (Sharma, Metere, Suryanarayana, Erlandson, Chow, Pask, *J. Chem. Phys.* 158, 204117, 2023): a modular, math-kernel-based implementation targeting NVIDIA architectures, offloading the computationally dominant operations to GPU while retaining the remainder on CPU. Reported speedups of **up to 6×** relative to CPU-only execution, bringing time-to-solution below **30 seconds** for a metallic system with over **14,000 electrons**.

2. **Hybrid-functional GPU acceleration** (Jing, Sharma, Pask, Suryanarayana, 2025, arXiv:2501.16572): introduces a **batched Kronecker-product-based linear solver** for the simultaneous solution of multiple linear systems arising in exact-exchange evaluation, with a modular GPU math-kernel implementation. Reported speedups of **up to 8× in node-hours** and **up to 80× in core-hours** relative to CPU-only execution, reducing time-to-solution on NVIDIA V100 GPUs to roughly **300 seconds** for a metallic system with over **6,000 electrons**.

These results indicate GPU acceleration substantially reduces both wall-clock time and the total computational resources (core-hours) required, which is particularly significant for computationally expensive hybrid-functional calculations that are otherwise impractical at scale with plane-wave methods.

### 5.4 Comparative context
Independent comparisons in the electronic-structure literature situate SPARC among a class of modern, HPC-oriented, real-space/finite-element/finite-difference DFT codes designed with exascale computing in mind — alongside codes such as **DFT-FE** (finite-element, adaptive), **RMG** (real-space multigrid, GPU-heavy on Summit/Frontier-class machines), **PARSEC**, and **Octopus**. SPARC differentiates itself from earlier finite-difference codes like PARSEC and Octopus principally through its electrostatic reformulation, consistent atomic-force/stress derivation, and the resulting combination of numerical accuracy and parallel efficiency demonstrated across a wide range of benchmark systems.

---

## 6. Accuracy and Validation

Across its foundational papers, SPARC has been systematically validated against converged plane-wave references (e.g., ABINIT):
- **High convergence rates** in energy and forces with respect to spatial discretization (mesh size), approaching plane-wave reference results.
- **Exponential convergence** in energies and forces with respect to vacuum size for slab and wire geometries.
- **Consistency between energies and forces** (forces derived correctly as the negative energy gradient within the discretization), with **negligible egg-box effect** — i.e., total energy and forces are essentially independent of how atoms are positioned relative to the underlying real-space grid, a key numerical-quality benchmark for real-space methods.
- **Accurate ground-state properties**: equilibrium geometries, vibrational spectra (clusters), and structural/electronic properties of crystals, slabs, and wires.
- **Negligible energy drift** in long ab initio molecular dynamics trajectories, indicating good energy conservation consistent with the underlying Hamiltonian.

---

## 7. Licensing, Availability, and Community

- **License:** GPL-3.0 (fully open source).
- **Repository:** [https://github.com/SPARC-X/SPARC](https://github.com/SPARC-X/SPARC) (organization: SPARC-X), with related repositories for M-SPARC and SPARC-X-API.
- **Language composition:** ~90% C, ~9% Python, remainder minor.
- **Installation:** Source compilation (C compiler + MPI, with BLAS/LAPACK, MKL, or ScaLAPACK options) or pre-compiled `conda-forge` binaries (`sparc-x` package) for x86_64/aarch64 Linux.
- **Latest tagged release (at time of writing):** v2.0.0 (May 2023), with ongoing development beyond that tag on the `master` branch (GPU support, hybrid functionals, DFT+U, and additional features have been merged/published subsequent to the v2.0.0 tag).
- **Funding/sponsorship:** U.S. Department of Energy Office of Science (DE-SC0023445, DE-SC0019410), DOE NNSA Advanced Simulation and Computing (ASC) Program (including DE-NA0004128 for high-temperature features), and U.S. National Science Foundation (award 1553212 for the Cyclix feature; earlier awards 1663244 and 1333500 for preliminary developments).
- **Documentation:** Distributed with the source (`doc/` directory) and via the SPARC-X-API documentation site, which also hosts maintenance/CI guidance for developers.

---

## 8. Assessment: Strengths and Limitations

**Strengths**
- Rigorous, mathematically consistent reformulation of electrostatics, forces, and stress within a real-space finite-difference framework — not merely a naive grid discretization, but one engineered specifically to avoid the pitfalls (egg-box effects, inconsistent forces) that historically limited earlier real-space codes.
- Demonstrated, published, order-of-magnitude performance and efficiency advantages over mature plane-wave codes at scale, with the gap widening with processor count and system size.
- Native handling of both isolated (Dirichlet) and periodic (Bloch) boundary conditions in one unified framework, avoiding artificial vacuum-padding workarounds needed by plane-wave codes for finite systems.
- Genuine linear-scaling (O(N)) capability via the SQ method, validated at the million-atom scale — a regime largely inaccessible to conventional cubic-scaling plane-wave DFT.
- Actively maintained with GPU acceleration for both semilocal and (more recently and non-trivially) hybrid functionals — an area where many legacy codes still lag.
- Minimal external dependencies and straightforward build process, aiding portability across HPC centers.
- Specialized Cyclix-DFT capability offers a distinctive advantage for symmetry-adapted simulation of nanotubes, nanowires, and other cyclic/helical nanostructures under deformation, a niche not well served by mainstream plane-wave packages.

**Limitations / considerations**
- As a comparatively young code (core papers from 2017 onward; unified public release in 2021), its feature breadth — especially in exotic post-processing, response-property, and excited-state methods — is still smaller than that of long-established plane-wave packages (VASP, Quantum ESPRESSO, ABINIT), though this gap is closing quickly (DFPT, RPA, and orbital-free DFT are in active development).
- GPU support, while demonstrated with strong speedups, is currently targeted specifically at NVIDIA architectures via CUDA-oriented math kernels; broader multi-vendor (AMD/Intel) GPU portability is not yet as mature as in some competing exascale-focused codes.
- Ecosystem (community size, third-party tutorials, user base) is smaller than that of the dominant plane-wave codes, though the SPARC-X-API/ASE integration is actively lowering this barrier.
- Norm-conserving ONCV pseudopotentials are the primary supported pseudopotential type; ultrasoft or PAW-type formalisms are not the code's native focus (unlike some plane-wave codes), which can affect basis-size/mesh-spacing requirements for certain element sets.

---

## 9. Key Publications — Theory and Methodology

The following publications represent the core theoretical, methodological, and software-engineering literature underpinning SPARC, organized by topic.

### 9.1 Foundational formulation
- S. Ghosh, P. Suryanarayana, "SPARC: Accurate and efficient finite-difference formulation and parallel implementation of Density Functional Theory: **Isolated clusters**," *Computer Physics Communications*, 212, 189–204 (2017). https://doi.org/10.1016/j.cpc.2016.09.020
- S. Ghosh, P. Suryanarayana, "SPARC: Accurate and efficient finite-difference formulation and parallel implementation of Density Functional Theory: **Extended systems**," *Computer Physics Communications*, 216, 109–125 (2017). https://doi.org/10.1016/j.cpc.2017.02.019

### 9.2 Unified software releases
- Q. Xu, A. Sharma, B. Comer, H. Huang, E. Chow, A. J. Medford, J. E. Pask, P. Suryanarayana, "SPARC: Simulation Package for Ab-initio Real-space Calculations," *SoftwareX*, 15, 100709 (2021) — **v1.0**. https://doi.org/10.1016/j.softx.2021.100709
- L. Zhang, A. Sharma, Q. Xu, et al., "SPARC v2.0.0: Spin-orbit coupling, dispersion interactions, and advanced exchange–correlation functionals," *Software Impacts* / arXiv:2305.07679 (2023). https://doi.org/10.1016/j.simpa.2024.100649

### 9.3 Electrostatics, linear solvers, and mixing
- Q. Xu, A. Sharma, P. Suryanarayana, "AAR: Alternating Anderson–Richardson method for fast solution of linear systems in real-space DFT," *Computer Physics Communications*, 234, 133–139 (2018). https://doi.org/10.1016/j.cpc.2018.07.007
- P. Suryanarayana, D. Phanish, "Augmented Lagrangian formulation of orbital-free density functional theory" / related electrostatics and linear-system work, *Journal of Computational Physics* (2015). https://doi.org/10.1016/j.jcp.2015.11.018
- P. P. Pratapa, P. Suryanarayana, J. E. Pask, "Anderson acceleration of the Jacobi iterative method: An efficient alternative to Krylov methods for large, sparse linear systems," *Chemical Physics Letters* (2016). https://doi.org/10.1016/j.cplett.2016.01.033
- P. Suryanarayana, P. P. Pratapa, J. E. Pask, "Restarted Pulay mixing for efficient and robust acceleration of fixed-point iterations," *Chemical Physics Letters* (2015). https://doi.org/10.1016/j.cplett.2015.06.029
- Q. Xu, A. Sharma, P. Suryanarayana, "Real-space mixing/preconditioning schemes for self-consistent field acceleration," *Chemical Physics Letters* (2019). https://doi.org/10.1016/j.cplett.2019.136983

### 9.4 Forces and stress tensor
- Atomic-force formulations are embedded in the foundational isolated-cluster and extended-system papers above (CPC 2017, 212 and 216).
- A. S. Banerjee, P. Suryanarayana, "Cyclic density functional theory: A route to the first principles simulation of bending in nanostructures" (stress/forces context for cyclic systems), *Journal of the Mechanics and Physics of Solids* (2016). https://doi.org/10.1016/j.jmps.2016.08.007
- Stress tensor and pressure formulation for real-space finite-difference DFT, *Journal of Chemical Physics* (2019). https://doi.org/10.1063/1.5057355

### 9.5 Non-orthogonal cells
- Formulation extending SPARC to non-orthogonal simulation cells, *Chemical Physics Letters* (2018). https://doi.org/10.1016/j.cplett.2018.04.018

### 9.6 O(N) linear-scaling Spectral Quadrature (SQ) method
- P. Suryanarayana, "Optimized Purification for Density Matrix calculation" / SQ formulation, *Chemical Physics Letters* (2013). https://doi.org/10.1016/j.cplett.2013.08.035
- P. Suryanarayana, P. P. Pratapa, A. Sharma, J. E. Pask, "SQDFT: Spectral Quadrature method for large-scale parallel O(N) Kohn–Sham calculations at high temperature," *Computer Physics Communications*, 197, 224–236 (2015). https://doi.org/10.1016/j.cpc.2015.11.005
- A. Sharma, P. Suryanarayana, "Spectral Quadrature method for accurate O(N) electronic structure calculations of metals and insulators," detailed mathematical formulation chapter, in *Computational Sciences: Recent Developments*, Springer (2023). https://doi.org/10.1007/978-3-031-22340-2_12
- Application to million-atom-scale simulation: *Modelling and Simulation in Materials Science and Engineering* (2023). https://doi.org/10.1088/1361-651X/acdf06

### 9.7 Cyclic and helical symmetry-adapted DFT (Cyclix-DFT)
- A. S. Banerjee, P. Suryanarayana, "Cyclic density functional theory: A route to the first principles simulation of bending in nanostructures," *Journal of the Mechanics and Physics of Solids*, 96, 605–631 (2016). https://doi.org/10.1016/j.jmps.2016.08.007
- S. Ghosh, A. S. Banerjee, P. Suryanarayana, "Symmetry-adapted real-space density functional theory for cylindrical geometries: Application to large twisted nanostructures," *Physical Review B*, 100, 125143 (2019). https://doi.org/10.1103/PhysRevB.100.125143
- S. Sharma, S. Ghosh, P. Suryanarayana, "Cyclic and helical symmetry-adapted real-space density functional theory," *Physical Review B*, 103, 035101 (2021). https://doi.org/10.1103/PhysRevB.103.035101
- Cyclic/helical phonon (DFPT) formalism and MLFF extension, *Journal of the Mechanics and Physics of Solids* (2024). https://doi.org/10.1016/j.jmps.2024.105927

### 9.8 GPU acceleration
- A. Sharma, A. Metere, P. Suryanarayana, L. Erlandson, E. Chow, J. E. Pask, "GPU acceleration of local and semilocal density functional calculations in the SPARC electronic structure code," *Journal of Chemical Physics*, 158, 204117 (2023). https://doi.org/10.1063/5.0225396 (arXiv:2302.09708)
- X. Jing, A. Sharma, J. E. Pask, P. Suryanarayana, "GPU acceleration of hybrid functional calculations in the SPARC electronic structure code," *Journal of Chemical Physics* (2025). https://doi.org/10.1063/5.0260892 (arXiv:2501.16572)

### 9.9 Pseudopotentials
- SPMS (soft and transferable) pseudopotential table, *Computer Physics Communications* (2022). https://doi.org/10.1016/j.cpc.2022.108594

### 9.10 DFT+U and Hubbard-corrected DFT
- Real-space Hubbard-corrected density functional theory implementation in SPARC, arXiv:2507.23612 (2025), describing energy, force, and stress-tensor expressions and the four-level parallelization hierarchy.

### 9.11 On-the-fly machine-learned force fields (MLFF)
- On-the-fly MLFF molecular dynamics formulation, *Journal of Chemical Physics* (2024). https://doi.org/10.1063/5.0180541
- Extension to the O(N) SQ method, *Journal of Chemical Physics* (2024). https://doi.org/10.1063/5.0204229
- Internal-energy/thermodynamic extensions, *Journal of Chemical Physics* (2024). https://doi.org/10.1063/5.0230060

### 9.12 Companion ATOM code
- Spectral scheme for isolated-atom electronic structure calculations (used for atomic-orbital generation, e.g., for DFT+U), *Computer Physics Communications* (2024). https://doi.org/10.1016/j.cpc.2024.109448

### 9.13 M-SPARC (MATLAB prototyping implementation)
- Q. Xu, A. Sharma, P. Suryanarayana, "M-SPARC: Matlab-Simulation Package for Ab-initio Real-space Calculations," *SoftwareX*, 11, 100423 (2020). https://doi.org/10.1016/j.softx.2019.100423
- Version 2.0.0 update to M-SPARC, *SoftwareX* (2023). https://doi.org/10.1016/j.softx.2022.101295 (companion to arXiv release cataloged under S2352711022002138)

### 9.14 Python interoperability
- T. Yao Zhang, A. Sharma, et al., "SPARC-X-API: Versatile Python Interface for Real-space Density Functional Theory Calculations," arXiv:2411.18024 (2024).

### 9.15 Roadmap / perspective context
- Contribution by P. Suryanarayana, Q. Xu, J. E. Pask, "SPARC: Advanced algorithms for systems large and small," Chapter 13 in *Roadmap on Electronic Structure Codes in the Exascale Era*, arXiv:2209.12747 (2022) — situates SPARC's design philosophy and future directions relative to the broader field of exascale electronic-structure software.

*Note: SPARC's own GitHub repository ("Citation" section) provides a curated, continuously updated list of the DOIs most relevant to each subsystem of the code (general citation, non-orthogonal systems, linear solvers, stress/pressure, forces, hybrid XC, mixing, pseudopotentials, cyclic/helical symmetry, SQ method, MLFF, and the ATOM code) and is the authoritative reference for citing specific SPARC functionality in derivative work.*

---

## 10. Summary

SPARC represents a mature, theoretically rigorous, and performance-competitive alternative to plane-wave Kohn–Sham DFT, purpose-built around the locality and boundary-condition flexibility of real-space finite-difference discretization. Its combination of (i) a carefully reformulated, numerically consistent electrostatics/force/stress framework, (ii) Chebyshev-filtered iterative eigensolvers, (iii) genuine O(N) linear scaling via the Spectral Quadrature method, (iv) symmetry exploitation via Cyclix-DFT, and (v) an actively expanding GPU-accelerated feature set (now including hybrid functionals) positions it as one of the more capable and forward-looking real-space DFT codes for large-scale, HPC-oriented first-principles materials simulation as of 2026.

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the SPARC 	Real-space finite-difference DFT code designed for accuracy and scalability on modern high-performance computing architectures. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
