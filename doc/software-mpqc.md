# MPQC — Massively Parallel Quantum Chemistry: An Exhaustive Review

## 1. Overview

MPQC (Massively Parallel Quantum Chemistry) is an open-source *ab initio* computational chemistry software package designed from its inception for parallel execution. It computes properties of atoms, molecules, and (in its current version) periodic solids from first principles by solving the time-independent Schrödinger equation. Development is led by the **Valeev Group** at **Virginia Tech**, with historical roots at **Sandia National Laboratories**, where the project originated.

MPQC is distinguished from legacy quantum chemistry packages (e.g., Gaussian, GAMESS) by three defining design choices:

- **Object-oriented architecture** implemented in C++, intended to make the codebase modular and extensible for method developers.
- **Parallelism as a first-class design goal**, rather than an afterthought bolted onto a serial code — MPQC was "written from the beginning to run in parallel," scaling from individual workstations to distributed-memory massively parallel supercomputers.
- **Fully open-source availability**, released under the **GNU General Public License (GPL)**, in contrast to most commercial or restricted-access quantum chemistry codes.

MPQC has gone through several major design generations, summarized below.

---

## 2. Version History

| Version | Era | Core Technology | Notes |
|---|---|---|---|
| MPQC 1–2 | 1990s–2000s | Scientific Computing Toolkit (SC): custom C++ core libraries for memory management, object state I/O, and parallel communication | Original Sandia-developed platform; established the object-oriented, parallel-native design philosophy |
| MPQC 3 | ~2000s–2010s | Continuation of SC-based `mpqc.Core` / `mpqc.Chemistry` library split | Split into a general-purpose C++ toolkit and QC-specific libraries; packaged for Linux distributions (Debian, Ubuntu) |
| MPQC 4 ("MPQC4") | 2018–present | **TiledArray** (massively-parallel block-sparse tensor framework), **MADWorld/MADNESS** (distributed task-based programming runtime), **Libint** (Gaussian integral evaluation library) | A ground-up reenvisioning of the original design on modern distributed task-based infrastructure; the release referenced by the current GitHub repository and most recent literature |

The final tagged release of the "classic" line is commonly cited as **MPQC 2.3 / v4, September 24, 2018**, with the codebase written in a mix of **C++, C, and FORTRAN 77**, targeting **Linux (Debian, Ubuntu) on x86-64**. The GitHub repository (`ValeevGroup/mpqc`) represents the actively developed MPQC4 codebase (81 stars, 27 forks, ~9,850 commits as of recent record).

---

## 3. Architectural Design

### 3.1 Two-Layer Component Structure

MPQC has historically comprised two components:

1. **A set of C++ libraries** for molecular electronic structure (quantum chemistry), of interest to software developers who want to implement new methods.
2. **An end-user program** for performing production QC computations, of interest to computational chemists applying standard methods.

The developer-facing libraries are organized as:

- **`mpqc.Core`** (formerly known as the **Scientific Computing Toolkit, "SC"**) — foundational C++ classes for memory management, object state serialization/restoration, input-file parsing, and parallel communication.
- **`mpqc.Chemistry`** — QC-specific libraries built atop `mpqc.Core`, used to implement both existing and new electronic-structure methods.

### 3.2 MPQC4 Technology Stack

The current (4th) generation rebuilds this architecture around three external, purpose-built parallel-computing frameworks, all also developed in close collaboration with the Valeev Group:

- **TiledArray** — a massively-parallel, block-sparse tensor framework that provides the fundamental distributed tensor data structures and contraction operations underlying the many-body methods in MPQC4.
- **MADWorld (MADNESS runtime)** — a distributed, task-based programming model and runtime system that manages asynchronous task scheduling and data movement across compute nodes.
- **Libint** — a library for the evaluation of molecular integrals of many-body operators over Gaussian basis functions, used for computing one- and two-electron (and higher-order) integrals.

This combination allows MPQC4 to express complex many-body methods (e.g., coupled-cluster) as compositions of tensor operations that are automatically distributed and scheduled across a cluster, rather than requiring hand-coded MPI parallelization for every new method — substantially lowering the barrier to implementing new high-accuracy methods that still scale to large machines.

### 3.3 Licensing

MPQC is distributed under the terms of the **GPL v3+**. Alternative licensing terms are available by contacting the Valeev Group directly. Earlier community packaging (e.g., Gentoo) lists **GPL-2** for older releases.

---

## 4. Computational Methods and Capabilities

### 4.1 Wavefunction / Mean-Field Methods

- **Hartree–Fock (HF)** theory:
  - Closed-shell (restricted, RHF)
  - Unrestricted (UHF)
  - General restricted open-shell (ROHF)
  - Energies and analytic gradients for all variants

- **Reduced-scaling / large-system Hartree–Fock**: a reduced-scaling LCAO formalism for periodic and molecular systems, using multipole-accelerated real-space summation and density fitting for the Coulomb potential, and concentric atomic density fitting for the exchange potential — enabling HF calculations on very large and periodic systems.

### 4.2 Density Functional Theory (DFT)

- **Kohn–Sham DFT**, implemented in parallel with the HF machinery:
  - Closed-shell, unrestricted, and general restricted open-shell Kohn–Sham DFT energies and analytic gradients.
- DFT is treated within the same modular integral/tensor infrastructure as HF, allowing exchange-correlation functional evaluation to be combined with the same distributed Coulomb/exchange machinery used for large-scale HF.

### 4.3 Correlated / Many-Body Methods

- **Møller–Plesset perturbation theory (MP2)**:
  - Second-order closed-shell MP2 energies and gradients.
  - **Explicitly correlated MP2-R12** methods (linear R12 correlation factor with an auxiliary basis set, ABS-MP2-R12), improving basis-set convergence for correlation energies.
  - Open-shell variants: **OPT2[2]** (second-order open-shell perturbation theory) and **ZAPT2** (Z-averaged perturbation theory).
- **Coupled-Cluster (CC) methods**:
  - Coupled-Cluster Singles and Doubles (**CCSD**) and related variants, implemented as a massively parallel tensor contraction workflow atop TiledArray.
  - **Explicitly correlated coupled-cluster (CCSD-F12/R12)** methods for accelerated basis-set convergence, including formulations suitable for **open-shell molecules with hundreds of atoms**.
- **Pair-natural-orbital (PNO) and reduced-scaling correlated methods**, including real-space/grid-based approaches for determining optimal pair-natural orbitals (e.g., for MP2-type wavefunctions).

### 4.4 Periodic / Solid-State Extensions

- Electronic structure methods extended to **periodic solids**, including reduced-scaling Hartree–Fock with Monkhorst–Pack k-point sampling.
- **Robust Pipek–Mezey (PM) orbital localization** for periodic systems, enabling localized Wannier-function-like orbitals in solids.
- **Concentric atomic density fitting** for efficient evaluation of exact exchange in periodic systems.

### 4.5 Geometry Optimization and Other Utilities

- An internal-coordinate **geometry optimizer**.
- Tensor-network approximation tools (e.g., robust approximation of tensor networks for grid-free factorization of the Coulomb interaction), used to accelerate integral and correlation-energy evaluation.

### 4.6 Summary Table of Methods

| Category | Methods |
|---|---|
| Mean-field | RHF, UHF, ROHF (energies + gradients); reduced-scaling/periodic HF |
| DFT | Closed-shell, unrestricted, and ROHF Kohn–Sham DFT (energies + gradients) |
| Perturbation theory | MP2, MP2-R12 (ABS), OPT2[2], ZAPT2 |
| Coupled cluster | CCSD, explicitly correlated CCSD-F12/R12 (incl. open-shell, large systems) |
| Localization | Robust Pipek–Mezey localization (molecular and periodic) |
| Solid state | Periodic Hartree–Fock, concentric atomic density fitting for exact exchange |
| Structural | Internal-coordinate geometry optimization |

---

## 5. Parallel Performance

MPQC (and MPQC4 specifically) targets architectures ranging from **individual workstations and symmetric multiprocessors to massively parallel distributed-memory supercomputers**. Demonstrated strong-scaling performance for the CCSD wavefunction solver (uracil trimer, 6-31G* basis, "BlueRidge" cluster at Virginia Tech) shows favorable parallel speed-up from a single 16-core node baseline (wall time ≈1290 s) across increasing node counts, illustrating the benefit of the TiledArray/MADWorld task-based design for many-body methods.

The MPQC4 platform has also been used as a testbed for large-scale multi-GPU tensor contraction research (block-sparse tensor contraction across distributed multi-GPU systems), extending its parallel model beyond CPU-only clusters.

---

## 6. Development History and Funding

- MPQC is described in its own literature as a **"30-year-old project"** as of 2020, tracing back to early-1990s development, historically associated with **Sandia National Laboratories** (operated for the U.S. DOE/NNSA).
- Present development is centered in the **Valeev Group**, Department of Chemistry, **Virginia Tech**.
- Development has been supported by:
  - **U.S. National Science Foundation** awards (CHE-0847295, CHE-0741927, OCI-1047696, CHE-1362655, ACI-1450262, ACI-1550456)
  - The **Alfred P. Sloan Foundation**
  - The **Camille and Henry Dreyfus Foundation**
  - The **U.S. Department of Energy Exascale Computing Project** (as part of the **NWChemEx** subproject)
  - The **DOE INCITE Program**

---

## 7. Availability and Distribution

- **Source code**: hosted at `https://github.com/ValeevGroup/mpqc` (MPQC4) and historically mirrored/archived on SourceForge (`mpqc.sourceforge.net`) and other community forges (e.g., `qsnake/mpqc` for MPQC 2.3).
- **Official site**: `www.mpqc.org`.
- **Linux packaging**: available as prebuilt/source packages for **Debian**, **Ubuntu**, and **Gentoo** (`sci-chemistry/mpqc`).
- **Interoperability**: supported as a backend by third-party GUI/analysis tools such as **Gabedit**, a graphical front end that also interfaces with GAMESS, Molpro, NWChem, Orca, Psi4, and other packages.
- MPQC is also referenced in academic surveys of open-source molecular modeling software as a representative parallel-native, object-oriented ab initio QC package.

---

## 8. Distinguishing Characteristics — Summary

1. **Parallel-first design**: unlike codes that retrofit MPI or GPU support onto a serial core, MPQC's architecture (both the classic SC-based design and the MPQC4 TiledArray/MADWorld design) was conceived around distributed execution from the start.
2. **Open-source and GPL-licensed**, unusual among historically dominant QC packages.
3. **Object-oriented C++ design**, intended to make new-method development (rather than just new-method *use*) tractable — a key differentiator for a "research platform" rather than solely a production tool.
4. **Modern tensor/task-based reengineering (MPQC4)**: reflects a broader trend in HPC quantum chemistry toward expressing many-body methods as distributed tensor algebra (shared conceptually with efforts like NWChemEx, of which MPQC4 development is formally a subproject).
5. **Breadth from mean-field to explicitly correlated coupled-cluster and periodic solids**, spanning both routine DFT/HF production calculations and cutting-edge high-accuracy correlated methods research.

---

## 9. Key Publications — MPQC Theory and Methods

The following publications document the theoretical and software-engineering foundations of MPQC and its core numerical methods, as cited in the official repository, developer publication lists, and related literature.

### 9.1 MPQC Platform / Software Architecture

- C. Peng, C. A. Lewis, X. Wang, M. C. Clement, K. Pierce, V. Rishi, F. Pavošević, S. Slattery, J. Zhang, N. Teke, A. Kumar, C. Masteran, A. Asadchev, J. A. Calvin, E. F. Valeev, "Massively Parallel Quantum Chemistry: A High-Performance Research Platform for Electronic Structure," *J. Chem. Phys.* **153**, 044120 (2020). https://doi.org/10.1063/5.0005889
- J. A. Calvin, C. Peng, V. Rishi, A. Kumar, E. F. Valeev, "Many-Body Quantum Chemistry on Massively Parallel Computers," *Chem. Rev.* (2021). https://doi.org/10.1021/acs.chemrev.0c00006

### 9.2 Tensor Infrastructure and Integral Evaluation

- C. A. Lewis, J. A. Calvin, E. F. Valeev, "Clustered Low-Rank Tensor Format: Introduction and Application to Fast Construction of Hartree–Fock Exchange," *J. Chem. Theory Comput.* **12**, 5868–5880 (2016). https://doi.org/10.1021/acs.jctc.6b00884
- T. Herault, Y. Robert, G. Bosilca, R. J. Harrison, C. A. Lewis, E. F. Valeev, J. J. Dongarra, "Distributed-Memory Multi-GPU Block-Sparse Tensor Contraction for Electronic Structure," *35th IEEE International Parallel & Distributed Processing Symposium (IPDPS)*, 537–546 (2021). https://doi.org/10.1109/IPDPS49936.2021.00062
- K. Pierce, V. Rishi, E. F. Valeev, "Robust Approximation of Tensor Networks: Application to Grid-Free Tensor Factorization of the Coulomb Interaction," *J. Chem. Theory Comput.* **17**, 2217–2230 (2021). https://doi.org/10.1021/acs.jctc.0c01310

### 9.3 Coupled-Cluster and Explicitly Correlated Methods

- C. Peng, J. A. Calvin, F. Pavošević, J. Zhang, E. F. Valeev, "Massively Parallel Implementation of Explicitly Correlated Coupled-Cluster Singles and Doubles Using the TiledArray Framework," *J. Phys. Chem. A* (2017). https://doi.org/10.1021/acs.jpca.6b10150
- A. Kumar, F. Neese, E. F. Valeev, "Explicitly Correlated Coupled Cluster Method for Accurate Treatment of Open-Shell Molecules with Hundreds of Atoms," *J. Chem. Phys.* **153**, 094105 (2020). https://doi.org/10.1063/5.0012753

### 9.4 Periodic / Solid-State Methods

- M. C. Clement, X. Wang, E. F. Valeev, "Robust Pipek–Mezey Orbital Localization in Periodic Solids," *J. Chem. Theory Comput.* **17**, 7406–7415 (2021). https://doi.org/10.1021/acs.jctc.1c00238
- X. Wang, C. A. Lewis, E. F. Valeev, "Efficient Evaluation of Exact Exchange for Periodic Systems via Concentric Atomic Density Fitting," *J. Chem. Phys.* **153**, 124116 (2020). https://doi.org/10.1063/5.0016856

### 9.5 Historical / Contextual References

- C. L. Janssen, I. M. B. Nielsen (Eds.), *Parallel Computing in Quantum Chemistry*, CRC Press, Boca Raton, FL (2008). ISBN 978-1-4200-5164-3. (Discusses MPQC's parallel, object-oriented design in the context of the broader field.)
- S. Pirhadi, J. Sunseri, D. R. Koes, "Open Source Molecular Modeling," *J. Mol. Graph. Model.* **69**, 127–143 (2016). https://doi.org/10.1016/j.jmgm.2016.07.008 (Surveys open-source computational chemistry software, including MPQC.)

*Note: earlier MPQC 1–3 development (1990s–2000s, Sandia National Laboratories era) is documented primarily in conference proceedings and technical reports (e.g., Sandia fact sheets, ClusterWorld articles by J. P. Kenny and C. L. Janssen) rather than in indexed journal articles; DOIs for these are not consistently available.*

---

## 10. Suggested Citation

Per the MPQC repository's `CITATION` file, users of the software are asked to cite the 2020 *Journal of Chemical Physics* platform paper (Peng et al., §9.1) as the primary reference, in addition to the specific method papers relevant to the calculation performed (e.g., the coupled-cluster or R12 papers in §9.3 if using those methods).

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the MPQC (Massively Parallel Quantum Chemistry) 	Open-source quantum chemistry package designed for parallel computation, including DFT methods. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
