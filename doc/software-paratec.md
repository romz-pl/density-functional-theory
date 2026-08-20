# PARATEC (PARAllel Total Energy Code) — Exhaustive Technical Review

## 1. Overview

PARATEC (**PARA**llel **T**otal **E**nergy **C**ode) is an *ab initio* density functional theory (DFT) software package that performs quantum-mechanical total-energy calculations using **norm-conserving pseudopotentials** and a **plane-wave basis set**. It was developed principally at **Lawrence Berkeley National Laboratory (LBNL)**, with contributing collaborators from the Université Pierre et Marie Curie, the University of Montreal, and the University of Cambridge, and additional contributions from F. Mauri, M. Côté, Y. Yoon, C. Pickard, and P. Haynes.

PARATEC is designed primarily for **massively parallel computing platforms** (clusters and supercomputers) but can also run on serial machines. It is one of a family of plane-wave pseudopotential DFT codes (alongside VASP, Quantum ESPRESSO, ABINIT, Qbox, CASTEP, CPMD, and others) developed in the 1990s–2000s HPC materials-science community.

- **Primary authors / maintainers:** Bernd Pfrommer, David Raczkowski, Andrew Canning, Steven G. Louie (LBNL)
- **Institution:** Lawrence Berkeley National Laboratory / NERSC
- **Language:** Fortran 90, with MPI for parallelism
- **Domain:** Condensed matter physics, materials science, computational chemistry
- **Official project page (historical):** nersc.gov/projects/paratec (NERSC-hosted)

## 2. Scientific Method

### 2.1 Theoretical Foundation
PARATEC solves the **Kohn–Sham equations of Density Functional Theory (DFT)**, the standard quantum-mechanical framework for computing the ground-state electronic structure of atoms, molecules, and periodic solids. DFT recasts the many-electron Schrödinger problem into a tractable set of single-particle equations governed by the electron density.

### 2.2 Basis Set and Pseudopotentials
- Electron wavefunctions are expanded in a **plane-wave basis set**, well suited to periodic (crystalline) systems.
- Core electrons are removed from the explicit calculation via **norm-conserving pseudopotentials**, which replace the strong Coulombic core potential with a smoother effective potential, reducing the number of plane waves needed for convergence.

### 2.3 Total Energy Minimization
PARATEC's defining algorithmic contribution is its **all-band, unconstrained conjugate-gradient (CG) minimization** of the Kohn–Sham total energy functional, described in Pfrommer, Demmel & Simon (1999), *"Unconstrained Energy Functionals for Electronic Structure Calculations,"* *J. Comput. Phys.* **150**(1), 287–298. Two solution strategies are supported:

1. **Self-consistent field (SCF) method** — the conventional iterative approach of alternating density and potential updates until convergence.
2. **Direct energy minimization** — an all-band CG approach minimizing the total energy functional directly; at the time of its introduction this was implemented only for systems with a band gap (insulators/semiconductors), avoiding some of the convergence difficulties SCF mixing can encounter.

### 2.4 Derived Properties
- **Forces** on atoms can be computed directly from the converged wavefunctions, enabling **structural relaxation** (geometry optimization) toward equilibrium atomic positions.
- **XANES (X-ray Absorption Near-Edge Structure)** calculations have been implemented within PARATEC's full-potential framework, used for simulating core-level X-ray absorption spectra (e.g., K-edge spectra in diamond and α-quartz, per Taillefumier et al. 2002, *Phys. Rev. B* **66**, 195107).

## 3. Parallel Computing Architecture

PARATEC was explicitly engineered as an **HPC/scientific-computing benchmark-class application** and is frequently cited in supercomputer evaluation literature.

### 3.1 Computational Profile
Reported profiling shows PARATEC spends the bulk of its runtime in:
- **Vendor-supplied BLAS3 routines** (~30% of runtime) — dense linear algebra (matrix–matrix operations)
- **FFT routines** (~30% of runtime) — transformations between real and reciprocal space

This BLAS3/FFT-dominated profile makes PARATEC well suited to a wide range of architectures, since both operation classes are typically heavily vendor-optimized.

### 3.2 Parallel FFT
A major engineering contribution associated with PARATEC is a **custom scalable parallel 3D FFT** designed specifically for electronic-structure codes, minimizing communication overhead in the global transpositions required by distributed 3D transforms. This FFT implementation has been benchmarked at scale (thousands of processors) on IBM SP, Cray XT, and NEC SX platforms. A related optimization — skipping FFT operations over the zero-valued rectangular region surrounding the non-zero sphere in reciprocal space — has since been adopted or compared against by other codes (e.g., CPMD, BigDFT, ONETEP).

### 3.3 Historical Performance Studies
PARATEC has served as a representative workload in several vector/parallel-architecture evaluation studies, notably:
- **"Scientific Computations on Modern Parallel Vector Systems"** (SC'04) — benchmarked PARATEC on the Earth Simulator (ES), Cray X1, and IBM Power3, using Si bulk systems of 432 and 686 atoms. Results showed the Earth Simulator outperforming the Cray X1 by more than 3.5× at 256 processors on the larger system, attributed to the X1's weaker interconnect bisection bandwidth limiting FFT transpose performance at higher processor counts.
- Cited as a reference workload in Knights Landing (Intel Xeon Phi) optimization studies for related codes such as BerkeleyGW.

## 4. Software Capabilities Summary

| Capability | Description |
|---|---|
| Ground-state DFT total energy | Plane-wave pseudopotential Kohn–Sham DFT |
| Pseudopotentials | Norm-conserving |
| SCF solver | Traditional self-consistent field mixing |
| Direct minimization | All-band conjugate-gradient, gap systems |
| Forces / relaxation | Atomic force calculation and structural relaxation |
| Spectroscopy | XANES (core-level X-ray absorption) |
| Parallelism | MPI-based; custom scalable 3D FFT |
| Target hardware | Massively parallel HPC systems; also runs serial |

## 5. Integration with BerkeleyGW

PARATEC has particularly close, purpose-built integration with **BerkeleyGW**, the many-body perturbation theory (GW/Bethe-Salpeter) package for excited-state properties:

- PARATEC is described in the BerkeleyGW manual as *"a simple DFT code optimized for small- and medium-sized systems"* with tight BerkeleyGW integration.
- PARATEC output relevant to BerkeleyGW is controlled by dedicated flags: `gw_shift`, `gwc`, `gwr`, `gwscreening`, `gwcscreening`, and `vxc_matrix_elements` (combinable, e.g. `gwr_gwscreening`; note `gwr`/`gwc` and `gwscreening`/`gwcscreening` are mutually incompatible pairs).
- Utilities **`bgw2para`** and **`rho2cd`** convert BerkeleyGW `WFN`/`RHO` files into PARATEC-format inputs (`WFN$n.$s`, `BAND`, `CD`), enabling workflows where plane waves generated elsewhere are diagonalized in PARATEC.
- Example DFT input files for PARATEC ship within the BerkeleyGW `examples/DFT` directory tree.
- As of the cited BerkeleyGW 2.2 (alpha) manual, the referenced PARATEC release is **version 5.1.12**.

## 6. Position in the DFT Software Landscape

PARATEC belongs to the plane-wave pseudopotential family of periodic DFT codes, alongside:

| Code | Basis/Method | Notes |
|---|---|---|
| **VASP** | Plane wave, PAW/pseudopotential | Widely used commercial code |
| **Quantum ESPRESSO (PWscf)** | Plane wave, pseudopotential | Open-source, GPL |
| **ABINIT** | Plane wave / wavelet, pseudopotential | Open-source |
| **Qbox** | Plane wave, pseudopotential | Designed for very large parallel scaling |
| **CASTEP** | Plane wave, pseudopotential | Commercial, UK academic origin |
| **CPMD** | Plane wave, pseudopotential | Car-Parrinello MD focus |
| **PARSEC** | Real space, pseudopotential | GPL |
| **BigDFT** | Wavelet basis | GPL, massively parallel |

Distinguishing PARATEC among these is its early emphasis on (a) an **unconstrained, all-band CG total-energy minimization scheme** as an alternative to standard diagonalization/SCF mixing, and (b) a **hand-tuned scalable parallel FFT** developed in tandem with vector- and massively-parallel supercomputer evaluation work at NERSC/LBNL — making it historically significant as much for its HPC engineering as for its DFT methodology.

## 7. Licensing and Availability

Available source material describes PARATEC principally through NERSC/LBNL project pages and its role as a companion code to BerkeleyGW, rather than through a standalone modern open public code repository (e.g., no actively maintained mainstream GitHub organization repository was found analogous to those of VASP, Quantum ESPRESSO, or ABINIT). Users seeking current access or licensing terms should consult NERSC/LBNL directly or the BerkeleyGW project, since PARATEC is distributed there as the reference DFT front end for BerkeleyGW workflows.

## 8. Key Publications

1. Pfrommer, B. G.; Demmel, J.; Simon, H. (1999). "Unconstrained Energy Functionals for Electronic Structure Calculations." *Journal of Computational Physics*, 150(1), 287–298. doi:10.1006/jcph.1998.6181
2. Taillefumier, M.; Cabaret, D.; Flank, A.-M.; Mauri, F. (2002). "X-ray absorption near-edge structure calculations with the pseudopotentials: Application to the K edge in diamond and α-quartz." *Physical Review B*, 66(19), 195107. doi:10.1103/PhysRevB.66.195107
3. Canning, A.; Wang, L. W.; Williamson, A.; Zunger, A. (SC'04 and related). "Scientific Computations on Modern Parallel Vector Systems," *Proceedings of the 2004 ACM/IEEE Conference on Supercomputing*.
4. Pfrommer, B.; Raczkowski, D.; Canning, A.; Louie, S. G. *PARATEC (PARAllel Total Energy Code)*. Lawrence Berkeley National Laboratory (with contributions from F. Mauri, M. Côté, Y. Yoon, C. Pickard, P. Haynes).

## 9. Summary Assessment

**Strengths**
- Mature, well-validated plane-wave pseudopotential DFT methodology with a documented and cited energy-minimization algorithm.
- Strong historical HPC pedigree: purpose-built parallel FFT, extensively benchmarked on leading-edge supercomputers of its era (Earth Simulator, Cray X1, IBM Power3/SP).
- Tight, well-documented interoperability with BerkeleyGW for GW/BSE excited-state calculations — a significant practical advantage for researchers running combined ground-state + many-body-perturbation-theory workflows.
- Efficient use of vendor BLAS3/FFT libraries gives portable performance across diverse HPC architectures.

**Limitations / Considerations**
- Documented primarily as optimized for **small- and medium-sized systems**, in contrast to codes explicitly engineered for very large-scale (thousands-of-atoms) simulation (e.g., PWDFT, ONETEP, CONQUEST).
- Direct energy-minimization mode historically limited to systems with a band gap (insulators/semiconductors), not general metals.
- Public documentation and independent tooling (packaging, active issue trackers, recent releases) are comparatively sparse relative to more actively promoted contemporary open-source codes (Quantum ESPRESSO, ABINIT, VASP), and it appears today mainly as the DFT front end bundled alongside BerkeleyGW rather than as a broadly, independently maintained standalone project.
- Norm-conserving-only pseudopotential support (no native PAW/ultrasoft formalism referenced in available material) may require larger plane-wave cutoffs for some elements compared to PAW-based codes.

**Overall**: PARATEC is a scientifically well-grounded, historically influential plane-wave DFT code notable for its energy-functional minimization algorithm and its role as a high-performance, tightly-coupled front end to the BerkeleyGW excited-state package. It remains most relevant today within GW/BSE workflows and within the historical literature on HPC performance engineering for electronic-structure codes.

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the PARATEC 	Parallel Total Energy Code is an ab initio density functional theory (DFT) software package. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
