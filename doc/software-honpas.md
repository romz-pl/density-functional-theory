# HONPAS: An Exhaustive Technical Review

## 1. Overview

**HONPAS** (**H**efei **O**rder-**N** **P**ackages for **A**b initio **S**imulations) is an open-source ab initio electronic structure code for linear-scaling (*O(N)*) first-principles density functional theory (DFT) calculations of large-scale periodic systems. It is built directly within the **SIESTA methodology**: standard norm-conserving pseudopotentials, strictly localized numerical atomic orbital (NAO) basis sets, and a real-space grid for the Hartree and exchange-correlation terms, under periodic boundary conditions.

HONPAS is developed and maintained by the group of **Jinlong Yang** (with **Xinming Qin** as lead/main developer and **Honghui Shang** as lead developer of the hybrid-functional module) at the **Hefei National Laboratory for Physical Sciences at the Microscale, University of Science and Technology of China (USTC)**. Its defining scientific contribution is bringing computationally tractable, linear-scaling **hybrid exchange-correlation functionals** (exact/Hartree-Fock exchange mixed with semi-local DFT) to large periodic systems expressed in an NAO basis — a regime where such calculations are traditionally prohibitively expensive.

- **Homepage:** http://honpas.ustc.edu.cn/
- **Source code:** public repository on GitHub (`xmqin/HONPAS`) and Bitbucket
- **License:** GPL-3.0 (inherited from/aligned with SIESTA's licensing)
- **Base version:** currently built on top of **SIESTA 4.1.5**
- **Language/platform:** Fortran (SIESTA codebase) with an interfaced C/C++ integral library; MPI + OpenMP hybrid parallelism

## 2. Relationship to SIESTA

HONPAS is not a from-scratch DFT code — it is a **derivative/extension of SIESTA** that inherits essentially the entire SIESTA feature set (self-consistent-field DFT, geometry relaxation, molecular dynamics, band structures, DOS, LDA/GGA functionals, k-point sampling, pseudopotentials, etc.) and adds a hybrid-functional engine on top. Practically:

- All standard SIESTA input parameters and calculation modules work unmodified in HONPAS, so it can be used as a conventional (semi-local) DFT code as well.
- Enabling a hybrid functional (e.g., HSE06) in a HONPAS run typically requires changing only a single input flag, `xc.authors HSE06`, in an otherwise standard SIESTA `.fdf` input file — all other hybrid-specific parameters default sensibly for most systems.
- Because it tracks SIESTA releases, HONPAS periodically resynchronizes with upstream SIESTA to inherit bug fixes and new features, while carrying its own hybrid-functional and linear-scaling additions forward.

This lineage is significant to note: as of 2026, the core NAO2GTO methodology pioneered in HONPAS has also been folded directly into mainline **SIESTA** itself (via an interface to the `libint` library), meaning HONPAS's central algorithmic contribution has influenced — and is converging with — the broader SIESTA ecosystem.

## 3. Core Methodology

### 3.1 Basis sets and pseudopotentials
- Strictly localized, finite-range **numerical atomic orbitals (NAOs)**, following SIESTA's pseudo-atomic-orbital (PAO) construction (single/multiple-ζ, polarized bases such as DZP).
- **Norm-conserving pseudopotentials** (Troullier–Martins type) for core-valence separation, generated with the standard SIESTA/ATOM toolchain.
- NAO strict locality gives Hamiltonian/overlap matrices natural sparsity, which is the structural property HONPAS exploits for linear scaling.

### 3.2 The NAO2GTO scheme (central innovation)
The central bottleneck in bringing hybrid functionals to an NAO-based code is the evaluation of **electron repulsion integrals (ERIs)**, which are needed for the nonlocal Hartree–Fock exchange (HFX) term but have no closed analytical form over numerical (tabulated) orbitals. HONPAS solves this with the **NAO2GTO** scheme:

1. Each numerical atomic orbital is fitted by a short, compact expansion of **Gaussian-type orbitals (GTOs)** (e.g., a "543"/"544" contraction pattern — 5, 4, 3 or 5, 4, 4 Gaussians per s-, p-, d-type shell).
2. The four-center ERIs and their analytic derivatives are then evaluated over the *fitted* GTOs using the mature, highly optimized **LIBINT** integral library (originally LIBINT 1.1.5, bundled with HONPAS).
3. Because the original NAOs are localized and short-ranged, the resulting shell-pair/shell-quartet ERI list is sparse and can be aggressively screened (Schwarz-type and distance-based screening), giving near-linear-scaling cost in system size rather than the conventional O(N⁴) formal scaling of exact exchange.
4. A later refinement (2023) replaces the *original* NAOs entirely with the *fitted* GTOs as the working numerical basis for the HFX/force terms, removing the small basis-mismatch error introduced by NAO2GTO fitting while preserving full analytic-derivative machinery.

### 3.3 Hybrid exchange-correlation functionals supported
Using the NAO2GTO engine, HONPAS implements total-energy **and** atomic-force (later, also stress/analytic-gradient) evaluation for:

- **PBE0** (global hybrid, 25% exact exchange)
- **B3LYP** (Becke three-parameter hybrid)
- **HSE06** (Heyd–Scuseria–Ernzerhof screened/range-separated hybrid) — the most extensively benchmarked and most widely used functional in the package, since screening the long-range part of exchange is what enables true linear scaling for extended/periodic systems.

### 3.4 Linear-scaling density matrix solvers
Rather than diagonalizing the Kohn–Sham Hamiltonian (formally O(N³)), HONPAS implements **density-matrix purification** methods that operate directly on sparse matrices, exploiting Kohn's "nearsightedness" principle:

- **Trace-preserving canonical purification** (Palser–Manolopoulos, PM)
- **Trace-correcting purification (TC / TC2)**, based on Niklasson's algorithm
- **Trace-resetting purification (TRS)**
- Spin-polarized extensions: **PSUTC2** and **SUTC2**, for spin-unrestricted systems with and without predetermined spin multiplicity, respectively
- Sparse matrix–matrix multiplication in compressed sparse row (CSR) format, parallelized via `MPI_Allgather`-based schemes, shown to scale to hundreds of cores for systems of tens of thousands of atoms (e.g., boron-nitride nanotubes).

### 3.5 Parallelization strategy
HONPAS uses a **two-level MPI + OpenMP** hybrid parallelization:
- **Static distribution algorithms**: ERIs distributed either over shell pairs or over shell quartets, with each ERI computed independently to minimize inter-process communication.
- **Dynamic master–worker distribution**: a two-level master–worker scheme dynamically load-balances ERI batches across cores, storing density and HFX matrices in sparse format and communicating only ERI indices and the final sparse HFX matrix — minimizing network traffic.
- Demonstrated scalability on systems including the **Tianhe-2** supercomputer, with benchmarks on molecular and periodic (solid-state, multi-thousand-atom) systems, and separately on iterative eigensolvers (Davidson, LOBPCG, PPCG, CheFSI) scaling to tens of thousands of processing cores for silicon systems with >10,000 atoms.

### 3.6 Accelerated exchange algorithms
- **Interpolative Separable Density Fitting (ISDF)** decomposition: constructs a low-rank approximation of the HFX matrix, avoiding explicit ERI evaluation and cutting the cost of PBE0-type HFX matrix construction by roughly two orders of magnitude in benchmark tests (e.g., benzene, polycyclic aromatic hydrocarbons), at controlled accuracy.
- **Machine-learning acceleration via DeepH**: a 2025 development interfaces HONPAS with **DeepH**, a deep-learning Hamiltonian-prediction framework, to bypass explicit SCF iteration. This lets HSE06-level hybrid calculations be extended to systems of **>10,000 atoms**, demonstrated on twisted bilayer graphene and twisted bilayer MoS₂ (moiré/incommensurate van der Waals systems), substantially cutting the time-to-solution for otherwise prohibitively expensive hybrid-functional calculations on large or aperiodic-in-approximation systems.

### 3.7 Post-Hartree–Fock capability (stated roadmap)
The project has indicated (per its group homepage) an intention to extend beyond DFT hybrids toward post-HF methods — **second-order Møller–Plesset perturbation theory (MP2)** and **coupled-cluster theory** — though these were described as "under development" rather than production features, and no completed, peer-reviewed implementation of full MP2/CC in HONPAS was identified as of this review.

## 4. Distinguishing Features vs. Other NAO/Hybrid Codes

| Aspect | HONPAS approach | Contrast |
|---|---|---|
| ERI evaluation | Fit NAOs → GTOs (NAO2GTO), then use LIBINT analytic GTO integrals | FHI-aims/ABACUS use resolution-of-identity (RI) with auxiliary NAO-like basis; CP2K uses Gaussian-and-plane-wave (GPW/GAPW) with its own RI/ADMM schemes; SIESTA (native, post-2025) now also interfaces libint directly |
| Scaling target | Explicit linear-scaling (O(N)) density-matrix purification + sparse ERI screening | Many hybrid-DFT LCAO codes remain cubic-scaling (direct diagonalization) unless specifically extended |
| Base infrastructure | Built as a SIESTA fork/extension, so inherits SIESTA's full conventional-DFT feature set for free | Purpose-built hybrid-DFT codes (e.g., CRYSTAL, Gaussian, FHI-aims) are not SIESTA derivatives |
| Force/stress support | Fully analytic HFX forces (2020) and, later, higher-accuracy analytic forces on the fitted-GTO basis (2023); more recent work also targets stress tensors for cell relaxation | Analytic hybrid-functional forces/stresses under periodic boundary conditions remain a nontrivial, actively developed capability across the NAO-DFT field generally |

## 5. Typical Application Domains

Based on the published benchmark and application studies, HONPAS (and its NAO2GTO methodology) has been used or benchmarked for:

- Bulk semiconductors and insulators — lattice constants, bulk moduli, and **band gaps** compared against experiment, addressing the well-known gap underestimation of semi-local (LDA/GGA) DFT via HSE06/PBE0.
- Two-dimensional materials.
- Point-defect and polaron physics — e.g., geometry optimization and **small-polaron formation from an excess electron in rutile TiO₂**, studied with HSE06 to validate the analytic-gradient implementation.
- Large-scale nanostructures — boron-nitride nanotubes with tens of thousands of atoms (linear-scaling density-matrix purification benchmark).
- Graphene quantum dots — room-temperature magnetism and tunable energy gaps in edge-passivated zigzag graphene quantum dots.
- Twisted van der Waals heterostructures / moiré systems — twisted bilayer graphene and twisted bilayer MoS₂, enabled at >10,000-atom scale via the DeepH–HONPAS interface.
- Molecular test systems (benzene, polycyclic aromatic hydrocarbons) for validating the ISDF-accelerated HFX approach.

## 6. Installation and Usage Notes

- **Dependency:** HONPAS requires the external **LIBINT** library (version 1.1.5 is bundled with the source distribution) for analytic two-electron integral evaluation over the fitted GTOs.
- **Build process:** LIBINT is compiled first (`./configure` / `make` / `make install`), then HONPAS itself is built from `Obj/` using an `arch.make` file essentially identical to a standard SIESTA build, but with an added `LIBINT_LIBS` path pointing to `libderiv.a` and `libint.a`.
- **Conventional DFT use:** because all SIESTA 4.1.5 input parameters and modules are supported unmodified, HONPAS can run ordinary LDA/GGA calculations exactly as SIESTA would.
- **Hybrid-functional use:** switching on a hybrid functional (e.g., HSE06) is done via a single `.fdf` keyword (`xc.authors HSE06`); other NAO2GTO-related parameters have working defaults for most systems.
- **Examples provided:** the repository ships a `HONPAS_Examples` directory with HSE06 band-structure input/output examples for several semiconductors, plus general `Examples/`, `Tests/`, and a `Tutorials/` directory, and a bundled PDF user manual (`HONPAS_manual.pdf`).
- **Repository activity (GitHub, `xmqin/HONPAS`):** a modest, low-traffic public repository (single-digit stars/forks range at time of writing), consistent with a specialized academic research code rather than a large community-maintained project; development activity and releases are better tracked through the associated journal publications than through GitHub release cadence.

## 7. Strengths and Limitations

**Strengths**
- Genuinely linear-scaling hybrid-functional total energies **and** analytic forces (and now progressing toward stresses), which is uncommon among NAO-based codes offering hybrid DFT.
- Reuses LIBINT, a mature, well-optimized GTO integral engine, rather than reinventing ERI evaluation from scratch.
- Full inheritance of SIESTA's mature conventional-DFT infrastructure (pseudopotentials, k-sampling, MD, structural relaxation, etc.).
- Demonstrated scalability on leadership-class HPC systems (Tianhe-2) to systems with thousands to tens of thousands of atoms.
- Active, ongoing methodological development into the mid-2020s (ISDF acceleration, fitted-GTO force refinement, DeepH machine-learning acceleration), and its methodology has been influential enough to be absorbed into mainline SIESTA.

**Limitations / caveats**
- The NAO2GTO fitting step introduces an approximation relative to using the "true" numerical orbitals directly; this has been an active source of algorithmic refinement (the 2023 fitted-NAO analytic-gradient work exists specifically to mitigate it).
- Community size and maintenance cadence appear modest relative to mainstream general-purpose DFT codes (VASP, Quantum ESPRESSO, CP2K); documentation is comparatively sparse outside the published papers and bundled manual/examples.
- Post-HF capabilities (MP2, coupled cluster) are described by the developers as under development rather than complete, production-ready features.
- Requires an external, somewhat dated dependency (LIBINT 1.1.5) with compiler-specific build instructions (Intel `icc`/`icpc`/`ifort` shown in the official build recipe), which may complicate builds on modern toolchains.
- As a derivative of a specific SIESTA version (4.1.5 at time of writing), HONPAS's feature set lags behind the very latest upstream SIESTA releases until a resync is performed.

## 8. Summary

HONPAS is a specialized, research-grade extension of SIESTA whose primary contribution to the field is the **NAO2GTO scheme**: an efficient, linear-scaling method for evaluating Hartree–Fock-type exact exchange (and its analytic derivatives) within a strictly localized numerical-atomic-orbital basis, by fitting NAOs to compact Gaussian expansions and leveraging the LIBINT integral library. This allows hybrid functionals (PBE0, B3LYP, and especially the screened HSE06 functional) — normally very expensive for large periodic systems — to be applied to systems ranging from small molecules to multi-thousand-atom (and, with recent machine-learning acceleration, >10,000-atom) periodic solids and van der Waals heterostructures, while retaining SIESTA's full conventional-DFT feature set. It remains an actively developed, narrowly-scoped academic code centered at USTC's Hefei National Laboratory, with its core methodology increasingly cross-pollinating into mainline SIESTA itself.

---

## 9. Publications Related to HONPAS Theory and Methodology

The following publications constitute the primary theoretical and methodological literature underpinning HONPAS, listed roughly in chronological/developmental order:

1. **Qin, X.; Shang, H.; Xu, L.; Hu, W.; Yang, J.; Xiang, H.** "Implementation of Exact Exchange with Numerical Atomic Orbitals." — earliest work implementing Hartree-Fock-type exact exchange in the SIESTA-based NAO framework using real-space Poisson-equation/interpolating-scaling-function evaluation of ERIs (precursor to NAO2GTO); validated against Gaussian03 and Crystal06.

2. **Shang, H.; Li, Z.; Yang, J.** "Implementation of Screened Hybrid Density Functional for Periodic Systems with Numerical Atomic Orbitals: Basis Function Fitting and Integral Screening." *J. Chem. Phys.* (2011) — foundational NAO2GTO fitting/screening paper referenced as the origin of the NAO2GTO scheme adopted by HONPAS.

3. **Qin, X.; Shang, H.; Xiang, H.; Li, Z.; Yang, J.** "HONPAS: A linear scaling open-source solution for large system simulations." *International Journal of Quantum Chemistry*, 115(10), 647–655 (2015). DOI: 10.1002/qua.24837 — the principal HONPAS package-description paper: introduces the NAO2GTO implementation of the HSE screened hybrid functional, ERI screening for linear scaling, and the PM/TC/TRS/(P)SUTC2 density-matrix purification algorithms.

4. **Shang, H.; Xu, L.; Wu, B.; Qin, X.; Zhang, Y.; Yang, J.** "The dynamic parallel distribution algorithm for hybrid density-functional calculations in HONPAS package." *Computer Physics Communications* (2020). arXiv: 2009.03555. — introduces the two-level master–worker dynamic parallel scheme for ERI/HFX-matrix distribution, benchmarked on Tianhe-2.

5. **Qin, X.; Shang, H.; Xu, L.; Hu, W.; Yang, J.; Li, S.; Zhang, Y.** "The static parallel distribution algorithms for hybrid density-functional calculations in HONPAS package." arXiv: 2009.03559 (2020). — presents the static shell-pair and shell-quartet ERI distribution algorithms as an alternative/complement to the dynamic scheme.

6. **Luo, Z.; Qin, X.; Wan, L.; Hu, W.; Yang, J.** "Parallel Implementation of Large-Scale Linear Scaling Density Functional Theory Calculations With Numerical Atomic Orbitals in HONPAS." *Frontiers in Chemistry*, 8:589910 (2020). DOI: 10.3389/fchem.2020.589910. — parallel sparse-matrix (CSR, `MPI_Allgather`) implementation of the trace-correcting (TC) density-matrix purification algorithm, demonstrated on boron-nitride nanotubes with tens of thousands of atoms.

7. **Zhang, Y.; Qin, X.; Hu, W.; Yang, J.; et al.** "Interpolative separable density fitting decomposition for accelerating Hartree-Fock exchange calculations within numerical atomic orbitals." arXiv: 2003.01654. — introduces the ISDF-based low-rank acceleration of the PBE0 HFX matrix construction in HONPAS.

8. **Qin, X.; Shang, H.; Yang, J.** "Efficient implementation of analytical gradients for periodic hybrid functional calculations within fitted numerical atomic orbitals from NAO2GTO." *Frontiers in Chemistry*, 11:1232425 (2023). DOI: 10.3389/fchem.2023.1232425. — replaces original NAOs with fitted GTOs as the working basis for HFX forces; analytic HSE06 force implementation validated on semiconductor lattice constants, bulk moduli, band gaps, and TiO₂ polaron formation.

9. **Ke, Y.; Qin, X.; Hu, W.; Yang, J.** "Combining DeepH with HONPAS for accurate and efficient hybrid functional electronic structure calculations with ten thousand atoms." *Digital Discovery*, 4, 2627–2638 (2025). DOI: 10.1039/D5DD00128E. — interfaces the DeepH machine-learning Hamiltonian framework with HONPAS to bypass SCF iterations, extending HSE06-level hybrid calculations to >10,000-atom twisted bilayer graphene/MoS₂ systems.

### Closely related methodological literature (broader NAO/hybrid-DFT context)

10. **Soler, J. M.; Artacho, E.; Gale, J. D.; García, A.; Junquera, J.; Ordejón, P.; Sánchez-Portal, D.** "The SIESTA method for ab initio order-N materials simulation." *J. Phys.: Condens. Matter*, 14(11), 2745–2779 (2002). — the foundational SIESTA methodology paper that HONPAS is built upon.

11. **García, A.; Papior, N.; Akhtar, A.; et al.** "Siesta: Recent developments and applications." *J. Chem. Phys.*, 152, 204108 (2020). DOI: 10.1063/5.0005077. — comprehensive update on modern SIESTA capabilities, forming the current baseline (SIESTA 4.1.x) that HONPAS extends.

12. **Pouillon, Y.; Oyomo, B. C.; Sifuna, J.; Camarasa-Gómez, M.; Qin, X.; Beltrán, C.; Gómez-Ortiz, F.; Shang, H.; Junquera, J.** "Implementation of the hybrid exchange-correlation functionals in the siesta code." *Computer Physics Communications*, 323, 110086 (2026). DOI: 10.1016/j.cpc.2026.110086. — describes the integration of a libint-based hybrid-functional (HFX) implementation directly into mainline SIESTA, explicitly positioned alongside CP2K's and HONPAS's NAO2GTO approach as prior art; involves a HONPAS co-developer (Xinming Qin) and signals convergence of the two codebases' hybrid-DFT infrastructure.

13. **Cao, Y.; Zhang, M.-Y.; Lin, P.; Chen, M.; Ren, X.** "Efficient Hybrid-Functional-Based Force and Stress Calculations for Periodic Systems with Thousands of Atoms." *J. Chem. Theory Comput.* (2025). DOI: 10.1021/acs.jctc.4c01635. — a parallel, independent NAO-based (ABACUS-family) implementation of linear-scaling analytic HFX forces/stresses using localized resolution-of-identity (LRI), useful as a direct comparison point to HONPAS's NAO2GTO force methodology.

14. **Lin, P.; Ren, X.; He, L.** "Efficient hybrid density functional calculations for large periodic systems using numerical atomic orbitals." *J. Chem. Theory Comput.*, 17(1), 222 (2021). — RI-based (rather than NAO2GTO-based) alternative approach to large-scale periodic hybrid-functional NAO calculations, relevant as comparative methodology.

---

*Note: Several of the arXiv preprints above (e.g., 2009.03555, 2009.03559, 2003.01654) correspond to work subsequently published in peer-reviewed venues (*Computer Physics Communications*, RSC/ACS journals); where the peer-reviewed DOI was confirmed it is given, otherwise the arXiv identifier is provided as the traceable reference.*

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the HONPAS 	Numerical atomic orbital DFT code based on the SIESTA methodology with hybrid functional support. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
