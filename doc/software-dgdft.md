# DGDFT: A Massively Parallel Discontinuous Galerkin Method for Large-Scale Electronic Structure Calculations — An Exhaustive Review

## 1. Overview

DGDFT (**D**iscontinuous **G**alerkin **D**ensity **F**unctional **T**heory) is a massively parallel electronic structure software package developed primarily by Wei Hu, Lin Lin, Chao Yang (Lawrence Berkeley National Laboratory / UC Berkeley) and John E. Pask (Lawrence Livermore National Laboratory), with contributions from Amartya S. Banerjee, Gaigong Zhang, Mathias Jacquelin, Eugene Vecharynski, and collaborators including Jinlong Yang's group at USTC. It is written in C++ and parallelized with MPI (with later hybrid MPI/OpenMP and GPU-oriented extensions), and is designed specifically to push Kohn–Sham density functional theory (KS-DFT) calculations to system sizes and processor counts far beyond the reach of conventional planewave codes — from thousands to millions of atoms, and from tens of thousands to tens of millions of processor cores.

The defining idea of DGDFT is to combine two normally competing strengths of standard basis sets:

- **Planewave basis sets**: systematically improvable (convergence controlled by a single cutoff parameter), but delocalized and therefore expensive to scale to large systems.
- **Localized basis sets** (Gaussians, numerical atomic orbitals): compact, but not systematically improvable and often require empirical tuning for new chemical environments.

DGDFT resolves this tension by constructing an **adaptive local basis (ALB)** set that is generated *on-the-fly*, self-consistently, during the SCF iteration, from local Kohn–Sham problems, and then coupling these local bases together using the mathematical machinery of the **discontinuous Galerkin (DG) finite element method**.

---

## 2. Theoretical and Methodological Foundations

### 2.1 Domain decomposition and the adaptive local basis (ALB)

DGDFT partitions the global simulation domain Ω into a set of non-overlapping rectangular subdomains ("elements") 𝒯 = {E₁, …, E_M}. For each element E_k, an **extended element** Q_k is defined that includes E_k plus a buffer region (typically its nearest-neighbor elements). A small, local Kohn–Sham eigenvalue problem is solved on Q_k (with periodic boundary conditions) using an auxiliary basis — typically planewaves, though the local problem can equally be solved with finite differences or other discretizations. The resulting local eigenfunctions are then **restricted and truncated back onto the central element E_k**, producing a set of basis functions that are:

- Strictly localized within each element (compact support),
- Discontinuous across element boundaries,
- Orthonormal within each element,
- Adaptive — they change from one SCF iteration to the next as the density/potential evolves, automatically encoding local atomic and environmental (chemical bonding) information into the basis.

Because the local problems capture the relevant local physics, chemical accuracy in total energies and forces can typically be reached with only a **few tens of ALB functions per atom** (and, with further optimization, as few as ~4–40 depending on the element/species), versus the hundreds of planewaves per atom that a global planewave calculation would require at the same real-space resolution.

### 2.2 Discontinuous Galerkin discretization of the Kohn–Sham Hamiltonian

Once the ALB set is generated, DGDFT discretizes the global Kohn–Sham Hamiltonian using the **interior penalty discontinuous Galerkin (IPDG) method**. This weak formulation adds boundary/penalty terms across element interfaces that enforce approximate continuity of the wavefunctions and their derivatives in a variational sense, without requiring the basis functions themselves to be continuous. As the number of ALBs per element is increased, the DG solution converges systematically to the exact (infinite-basis, e.g. planewave-limit) Kohn–Sham solution, and convergence can be monitored rigorously via **a posteriori error estimators** derived from finite-element analysis.

A key structural payoff of this construction is that the resulting DG Hamiltonian matrix is:

- **Sparse and block-structured** (block-tridiagonal-like, reflecting nearest-neighbor element coupling),
- Regular in structure, resembling a finite-difference-type discretization,
- This structure holds even for **metallic systems**, unlike many linear-scaling methods that rely on the nearsightedness/locality of the density matrix (which fails for metals) — DGDFT's sparsity comes from the basis construction itself, not from decay of the density matrix, so it is applicable to both insulators and metals.

This sparse block structure is what enables efficient parallel Hamiltonian construction, matrix-vector operations, and scalable eigensolvers/density-matrix solvers.

### 2.3 Electron density, energy, and force evaluation: diagonalization alternatives

Given the sparse DG Hamiltonian, DGDFT supports multiple strategies for computing the density matrix (and hence electron density, total energy, atomic forces):

- **Direct diagonalization (DIAG)**: standard dense/parallel eigensolvers (e.g., ScaLAPACK-based); simplest but scales cubically (O(N³)) with system size and becomes the bottleneck for large systems.
- **Chebyshev-polynomial Filtered Subspace Iteration (CheFSI)**: uses Chebyshev polynomial filters to amplify the occupied-subspace eigencomponents, avoiding full diagonalization at every SCF step; substantially faster than DIAG (reported >2 orders of magnitude speed-up) while retaining the accuracy benefits of the DG Hamiltonian's sparsity, though the cubic-scaling subspace (Rayleigh–Ritz) step still eventually dominates at very large scale.
- **Pole Expansion and Selected Inversion (PEXSI)**: represents the Fermi–Dirac operator via a rational (pole) expansion and computes only the selected entries of the resulting matrix inverses needed to construct the density matrix, energy, and forces — entirely avoiding diagonalization. Combined with the DG Hamiltonian's sparsity, this yields sub-cubic complexity: at most **quadratic (O(N²))** scaling in general and as low as **O(N^1.5)** for quasi-2D systems, and is highly parallel-scalable (demonstrated to >100,000 processors independent of DGDFT).
- **Two-level Chebyshev filter / complementary subspace method**: a refinement that further pushes scalability by separating the eigenspectrum into "top" and "bottom" complementary subspaces treated with different numerical strategies, extending applicability to larger systems and higher core counts than single-level CheFSI.
- **ELSI interoperability**: DGDFT has been integrated with the **ELSI** (ELectronic Structure Infrastructure) solver interface, giving it access to a broader family of interchangeable eigensolvers/density-matrix solvers (ELPA, PEXSI, NTPoly, etc.) through a common API shared with other codes such as FHI-aims, DFTB+, and SIESTA.

### 2.4 Two-level parallelization strategy

DGDFT's massive scalability rests on a **two-level (nested) MPI parallelization** scheme:

1. **Outer level — element-based parallelization**: the global domain's elements are distributed across groups of MPI processes (a 2D process grid), with each element/group responsible for constructing its local Hamiltonian block and communicating boundary information only with neighboring elements ("fixed neighbor communication"), keeping most communication local and minimizing global collectives.
2. **Inner level — within-element parallelization**: within each element/extended-element group, further parallelization is applied to the local Kohn–Sham solve used to generate the ALB (band parallelization with column-cyclic FFT distribution, and grid parallelization with row-block GEMM distribution for the tall-and-skinny wavefunction matrices).

This nested strategy, combined with the block-sparse DG Hamiltonian, is what has allowed DGDFT to demonstrate strong/weak scaling to extraordinarily large processor counts (see benchmarks in §4).

### 2.5 Beyond ground-state energetics: forces, dynamics, and excited states

- **Analytic atomic forces**: DGDFT computes Hellmann–Feynman-type atomic forces consistently within the DG/ALB framework (including Pulay-type correction terms arising from the basis being atom/environment-dependent), enabling geometry optimization.
- **Ab initio molecular dynamics (AIMD)**: forces computed at each MD step allow Born–Oppenheimer AIMD simulations of large systems (thousands to tens of thousands of atoms) at finite temperature.
- **Real-time time-dependent DFT (rt-TDDFT)** via the **parallel transport (PT) gauge**: DGDFT implements a gauge-optimal formulation of the time-dependent Kohn–Sham equations that removes fast oscillatory phases associated with the usual gauge choice, permitting much larger time steps in explicit propagation; this has been used to reach electron-dynamics simulations of unprecedented scale (up to million-atom regimes in later "DG-TDDFT" extensions built on the DGDFT framework).
- **Hybrid exchange–correlation functionals**: because exact-exchange (Fock exchange) evaluation is normally the dominant cost for hybrid functionals (e.g., HSE06), DGDFT/PWDFT implement acceleration techniques including the **Adaptively Compressed Exchange (ACE)** operator formulation and **Interpolative Separable Density Fitting (ISDF)**, including an ISDF variant based on **Centroidal Voronoi Tessellation (CVT)** selection of interpolation points, and a **Projected Commutator DIIS (PC-DIIS)** method for accelerating SCF convergence specifically under the nonlinear eigenvalue problem hybrid functionals induce.

### 2.6 PWDFT — the embedded planewave module

DGDFT ships with a **self-contained planewave DFT module, PWDFT**, used mainly for validation, benchmarking, and as the local-problem solver, but also functional as a standalone standard planewave DFT code (Γ-point, periodic boundary conditions, norm-conserving Hartwigsen–Goedecker–Hutter (HGH) pseudopotentials, LibXC-based exchange-correlation, real-space pseudo-charge formulation). Its accuracy has been cross-validated against established planewave codes such as Quantum ESPRESSO. Most of the hybrid-functional acceleration algorithms above (ACE, ISDF, PC-DIIS) were developed and demonstrated first within PWDFT before/alongside their use in the full DG framework.

---

## 3. Software Architecture and Implementation

| Aspect | Details |
|---|---|
| **Language** | C++ |
| **Parallelization** | MPI (two-level/nested); later work adds hybrid MPI + heterogeneous (many-core / GPU-style) parallelism |
| **Core numerical dependencies** | FFTW (planewave/FFT operations), BLAS/LAPACK/ScaLAPACK (dense linear algebra), PEXSI library (sparse pole-expansion/selected-inversion solver), LibXC (exchange-correlation functionals) |
| **Solver interoperability** | ELSI infrastructure (ELPA, PEXSI, NTPoly, etc.) |
| **Pseudopotentials** | Norm-conserving pseudopotentials (e.g., HGH-type), real-space pseudo-charge formulation |
| **Modules** | Core DG/ALB engine; **PWDFT** (embedded planewave module); DG-TDDFT extension for real-time dynamics |
| **Deployment/scale demonstrated** | From modest clusters up to leadership-class systems: NERSC Edison (up to 128,000 cores), Sunway TaihuLight (millions of cores, SW26010 heterogeneous processor), and newer Sunway systems (tens of millions of cores) |
| **License / distribution** | Developed and distributed as a research HPC code out of LBNL (Computational Research Division) under the U.S. DOE SciDAC program, with source historically hosted via LBNL's code repository infrastructure; typically obtained by request/collaboration from the developing groups |

### Computational workflow (per SCF iteration)

The main SCF loop in DGDFT consists of four computationally significant stages:
1. **ALB generation** (inner SCF): solve local KS problems on each extended element to produce the adaptive basis.
2. **DG Hamiltonian construction**: assemble the global sparse Hamiltonian from element/interface contributions using the IPDG formulation.
3. **Diagonalization / density-matrix solve** (outer SCF): via DIAG, CheFSI, or PEXSI depending on system size/type and available resources.
4. **Density, energy, and force evaluation**, feeding back into the next SCF cycle.

---

## 4. Performance, Scalability, and Landmark Demonstrations

DGDFT's central design goal — and its most distinctive achievement relative to other DFT codes — is extreme parallel scalability paired with retained planewave-level accuracy. Key reported benchmarks include:

- **Basis compactness**: chemical accuracy (errors on the order of 10⁻⁴ Hartree/atom in energy) achieved with roughly **15–40 ALB functions per atom**, versus hundreds needed for equivalent-accuracy planewave calculations.
- **NERSC Edison benchmark**: **80% parallel efficiency at 128,000 CPU cores** studying 2D phosphorene systems of 3,500–14,000 atoms, using the two-level parallelization scheme.
- **Sunway TaihuLight benchmark**: parallel efficiency of **32.3%** (speedup ~42,383×) using **8,519,680 processor cores** (131,072 core groups) for 2D metallic graphene systems (2,880 and 11,520 carbon atoms), leveraging the SW26010 heterogeneous master–slave many-core architecture.
- **Extreme-scale metallic heterostructure simulation**: on a newer Sunway supercomputer, DGDFT + PEXSI computed the electronic structure of a system with **2.5 million atoms (17.2 million electrons)** using **35.9 million cores**, with PEXSI's sparse direct solver reaching a peak of **~64 PFLOPS** (~5% of theoretical peak) — described as unprecedented for sparse direct solvers at that scale, opening a path toward mesoscale quantum-mechanical device simulation.
- **AIMD throughput**: for an 8,000-atom bulk silicon system at finite temperature, an average SCF step wall time of ~51 seconds on 34,560 processors was reported, enabling ~1.0 ps of ab initio MD in ~28 hours of wall time (via the parallel-transport-gauge-enabled dynamics implementation).
- **Million-atom real-time TDDFT**: the DG-TDDFT extension (built on the DGDFT framework, adding hybrid-functional and excited-state support) has been demonstrated for electron-dynamics simulations at the million-atom scale, exploiting the same fixed-neighbor sparse communication pattern.

Complexity scaling achieved: at most **O(N²)** overall with DIAG/CheFSI-class solvers replaced or supplemented by **O(N^1.5)–O(N)**-scaling PEXSI-based density-matrix evaluation for quasi-2D and sufficiently sparse systems, applicable uniformly to **both insulating and metallic** systems (a notable distinction from many linear-scaling DFT methods that rely on density-matrix locality/nearsightedness and hence struggle with metals).

---

## 5. Representative Application Domains

DGDFT has been applied to a range of large-scale materials and molecular systems that stress-test both accuracy and scalability:

- **2D materials**: graphene nanoflakes (electronic structure and aromaticity of large hexagonal flakes), 2D phosphorene monolayers and nanoribbons (including armchair-edge reconstruction physics and edge-controlled phosphorene nanoflake heterojunctions).
- **Bulk semiconductors and alloys**: bulk silicon systems of varying size, disordered silicon–aluminum alloys, used extensively as benchmark/validation systems for hybrid-functional acceleration methods (ACE, ISDF, PC-DIIS).
- **Point defects**: defect electronic structure in silicon (e.g., vacancy defects) under hybrid functionals.
- **Metallic heterostructures**: large-scale (multi-million-atom) metallic heterostructure electronic structure for next-generation electronic device design at mesoscopic scale.
- **Finite-temperature ab initio molecular dynamics**: large bulk silicon systems (thousands of atoms) as a demonstration of the parallel-transport-gauge dynamics scheme.
- **Real-time electron dynamics**: ultrafast/excited-state electron dynamics simulations at unprecedented (million-atom) system sizes via DG-TDDFT.

---

## 6. Comparison to Related Codes and Methodological Context

DGDFT belongs to a family of "systematically improvable, reduced-basis, real-space-adjacent" DFT methods aimed at large-scale HPC, alongside codes such as:

- **BigDFT** and **ONETEP** — also generate localized basis functions on-the-fly, but via wavelets (BigDFT) or continuous, optimization-refined localized orbitals (ONETEP), as opposed to DGDFT's strictly local, orthonormal, discontinuous ALB functions defined through the DG framework.
- **DFT-FE** — a finite-element-discretization DFT code with comparable large-scale ambitions and hybrid CPU-GPU acceleration; conceptually related to DGDFT's DG/finite-element lineage but uses a continuous (not discontinuous) finite-element basis and has emphasized GPU acceleration more explicitly and earlier in its public releases.
- **Linear-scaling (O(N)) codes** such as CONQUEST, SIESTA (linear-scaling mode), LS3DF, RSDFT — these typically exploit the nearsightedness/decay of the density matrix, which is effective for insulators/large-gap systems but breaks down for metals; DGDFT's sparsity instead derives from the ALB/DG construction itself, so it remains applicable — and has been explicitly benchmarked — on metallic systems.
- **PEXSI and ELSI** — DGDFT was one of the original and principal use cases motivating the PEXSI sparse solver library, and both projects share core developers (notably Lin Lin and Chao Yang); DGDFT is also one of the codes integrated with the broader ELSI solver-interoperability infrastructure alongside FHI-aims, DFTB+, and SIESTA.
- **A more recent, related line of work** ("Fast adaptive discontinuous basis sets for electronic structure," 2025) explicitly builds on and extends the DGDFT approach by allowing hybrid combinations of atom-centered Gaussian-type functions and polynomials within the same DG framework, aiming to reduce the overhead of diagonalizing local Hamiltonians that the original ALB construction requires.

---

## 7. Summary Assessment

**Strengths**
- Combines systematic improvability (planewave-like convergence control) with compact basis size (localized-orbital-like efficiency).
- Rigorous mathematical foundation (finite-element/DG theory) enabling a posteriori error control.
- Sparse, regular Hamiltonian structure that holds for both metals and insulators, unlike many linear-scaling alternatives.
- Demonstrated scalability that is, at the time of its landmark publications, among the most extreme reported in DFT software (tens of millions of cores, millions of atoms).
- Multiple interchangeable density-matrix/eigensolver back-ends (DIAG, CheFSI, PEXSI, ELSI-mediated solvers) letting users trade accuracy/robustness against speed depending on system size and metallicity.
- Extensions covering forces, AIMD, hybrid functionals, and real-time TDDFT — not merely a ground-state energy engine.

**Limitations / considerations**
- The local-problem solves needed to generate ALBs (diagonalizing local Hamiltonians on each extended element every SCF iteration) add non-trivial overhead compared to fixed, precomputed basis sets — a design trade-off explicitly noted even by follow-on work seeking to reduce it.
- As a specialized HPC research code (rather than a general-purpose, broadly distributed community package like Quantum ESPRESSO or VASP), DGDFT's user base, documentation, and general accessibility are narrower; it is primarily used by its developing groups and close collaborators, and demonstrated benchmarks emphasize periodic, largely bulk/2D crystalline and metallic systems rather than the full breadth of molecular chemistry use cases.
- Best performance advantages materialize specifically at very large scale (thousands+ atoms, large core counts); for small/medium systems the ALB-generation overhead and two-level parallelization machinery are less likely to outweigh simpler planewave or Gaussian-basis approaches.

---

## 8. Publications Related to DGDFT's Theory and Methodology

### Foundational method and code papers
1. L. Lin, J. Lu, L. Ying, W. E, "Adaptive local basis set for Kohn–Sham density functional theory in a discontinuous Galerkin framework I: Total energy calculation," *J. Comput. Phys.* **231**, 2140–2154 (2012).
2. G. Zhang, L. Lin, W. Hu, C. Yang, J. E. Pask, "Adaptive local basis set for Kohn–Sham density functional theory in a discontinuous Galerkin framework II: Force, vibration, and molecular dynamics calculations," *J. Comput. Phys.* **335**, 426–443 (2017).
3. W. Hu, L. Lin, C. Yang, "DGDFT: A massively parallel method for large scale density functional theory calculations," *J. Chem. Phys.* **143**, 124110 (2015).

### Eigensolver / density-matrix methodology
4. A. S. Banerjee, L. Lin, W. Hu, C. Yang, J. E. Pask, "Chebyshev polynomial filtered subspace iteration in the discontinuous Galerkin method for large-scale electronic structure calculations," *J. Chem. Phys.* **145**, 154101 (2016).
5. A. S. Banerjee, L. Lin, P. Suryanarayana, C. Yang, J. E. Pask, "Two-level Chebyshev filter based complementary subspace method: pushing the envelope of large-scale electronic structure calculations" (arXiv:1712.04439; *J. Chem. Theory Comput.*).
6. L. Lin, J. Lu, L. Ying, R. Car, W. E, "Fast algorithm for extracting the diagonal of the inverse matrix with application to the electronic structure analysis of metallic systems," *Commun. Math. Sci.* **7**, 755–777 (2009). *(PEXSI diagonal-extraction precursor.)*
7. L. Lin, C. Yang, J. C. Meza, J. Lu, L. Ying, W. E, "SelInv—An algorithm for selected inversion of a sparse symmetric matrix," *ACM Trans. Math. Softw.* **37**, 40 (2011).
8. L. Lin, M. Chen, C. Yang, L. He, "Accelerating atomic orbital-based electronic structure calculation via pole expansion and selected inversion," *J. Phys.: Condens. Matter* **25**, 295501 (2013). *(PEXSI method core reference.)*

### Hybrid functional acceleration
9. L. Lin, "Adaptively compressed exchange operator," *J. Chem. Theory Comput.* **12**, 2242–2249 (2016) (arXiv:1601.07159).
10. K. Dong, W. Hu, L. Lin, "Interpolative separable density fitting through centroidal Voronoi tessellation with applications to hybrid functional electronic structure calculations," *J. Chem. Theory Comput.* **14**, 1311–1320 (2018) (arXiv:1711.01531).
11. W. Hu, L. Lin, C. Yang, "Interpolative separable density fitting decomposition for accelerating hybrid density functional calculations with applications to defects in silicon," *J. Chem. Theory Comput.* **13**, 5420–5431 (2017) (arXiv:1707.09141).
12. W. Hu, L. Lin, C. Yang, "Projected commutator DIIS method for accelerating hybrid functional electronic structure calculations," *J. Chem. Theory Comput.* **13**, 5458–5467 (2017) (arXiv:1708.06485).

### Dynamics (AIMD and real-time TDDFT)
13. J. Zhang, X. Wang, L. Lin, C. Yang, "Quantum dynamics with the parallel transport gauge," (parallel transport gauge for TDDFT — foundational formulation).
14. G. Zhang, L. Lin, "Fast real-time time-dependent density functional theory calculations with the parallel transport gauge," *J. Chem. Theory Comput.* (2018).
15. W. Hu et al., "Million-Atom Ab Initio Electron Dynamics: Discontinuous Galerkin Real-Time Time-Dependent Density Functional Theory," *Proceedings of SC'25 (International Conference for High Performance Computing, Networking, Storage and Analysis)* (2025).

### Applications demonstrating/validating the methodology
16. W. Hu, L. Lin, C. Yang, J. Yang, "Electronic structure and aromaticity of large-scale hexagonal graphene nanoflakes," *J. Chem. Phys.* **141**, 214704 (2014).
17. W. Hu, L. Lin, C. Yang, "Edge reconstruction in armchair phosphorene nanoribbons revealed by discontinuous Galerkin density functional theory," *Phys. Chem. Chem. Phys.* **17**, 31397–31404 (2015).
18. W. Hu, L. Lin, C. Yang, J. Yang, "Edge-controlled phosphorene nanoflake heterojunctions," (companion phosphorene nanoflake study).

### Extreme-scale HPC implementation papers
19. Z. Wu et al. (incl. W. Hu), "Extreme-scale density functional theory high performance computing of DGDFT for tens of thousands of atoms using millions of cores on Sunway TaihuLight," (arXiv:2003.00407; SC'20 proceedings).
20. W. Hu et al., "High performance computing of DGDFT for tens of thousands of atoms using millions of cores on Sunway TaihuLight," *Science Bulletin* (2021).
21. Q. Liu et al. (incl. W. Hu), "2.5 Million-Atom Ab Initio Electronic-Structure Simulation of Complex Metallic Heterostructures Using DGDFT" — Gordon Bell-class paper, *Proceedings of SC'22*.

### Solver interoperability infrastructure
22. V. W. Yu et al., "ELSI: A unified software interface for Kohn–Sham electronic structure solvers," *Comput. Phys. Commun.* **222**, 267–285 (2018), and subsequent ELSI infrastructure papers describing DGDFT/PEXSI/ELPA/NTPoly interoperability.

### Related/extending recent work
23. "Fast adaptive discontinuous basis sets for electronic structure" (arXiv:2510.21213, 2025/2026) — extends the DGDFT-style discontinuous Galerkin adaptive basis idea to hybrid Gaussian/polynomial primitive functions.

*Note: Item numbering above groups papers thematically rather than strictly chronologically; several titles (items 13, 18) are cited as they consistently appear as companion/in-preparation references across the primary DGDFT papers and the developers' own project pages, and readers wishing to pin exact volume/page/DOI details should cross-check against the developers' publication list or each journal's official record, as some entries (especially very recent conference papers) were still finalizing bibliographic details at the time of writing.*

---

## 9. Suggested Entry Points for Further Reading

- Start with the two **Adaptive Local Basis** papers (Lin et al. 2012; Zhang et al. 2017) for the mathematical formulation.
- Read the **DGDFT: A massively parallel method** paper (Hu, Lin, Yang, *J. Chem. Phys.* 2015) for the software/parallelization architecture.
- For solver internals, see the **CheFSI** and **PEXSI**-related papers.
- For the extreme-scale HPC engineering, see the **Sunway TaihuLight** and **Gordon Bell (SC'22)** papers.
- For hybrid-functional and dynamics capabilities, see the **ACE/ISDF/PC-DIIS** and **parallel transport gauge / DG-TDDFT** papers respectively.

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the DGDFT 	Discontinuous Galerkin DFT code designed for massively parallel, large-scale electronic structure calculations. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
