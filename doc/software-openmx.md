# OpenMX: A Comprehensive Review

**Open source package for Material eXplorer — DFT with pseudo-atomic localized basis functions and norm-conserving pseudopotentials**

---

## 1. Overview

OpenMX ("Open source package for Material eXplorer") is a software package for nanoscale material simulations based on density functional theory (DFT), norm-conserving pseudopotentials, and variationally optimized pseudo-atomic localized basis functions (PAOs). Development was initiated by Taisuke Ozaki in 2000 (originally at JRCAT/AIST, later ISSP, University of Tokyo), and the code has since been extended by a distributed group of contributors. It is distributed under the GNU General Public License version 3 (GPLv3) and is freely downloadable, with source code, precompiled pseudopotential/basis databases, and full documentation hosted at openmx-square.org.

The defining architectural choice of OpenMX is its use of a **strictly localized numerical atomic-orbital (LCPAO) basis** combined with **norm-conserving pseudopotentials**, evaluated on a **real-space grid** (via FFT-based Poisson solves). This combination — shared conceptually with codes such as SIESTA, CONQUEST, FHI-aims, and ONETEP — makes the Hamiltonian and overlap matrices sparse, which in turn enables genuinely large-scale (and, with O(N) solvers, linear-scaling) electronic structure calculations on parallel machines, while retaining an all-electron-like or plane-wave-comparable accuracy once the local basis is properly converged (as demonstrated in the Δ-gauge benchmark).

---

## 2. Theoretical and Methodological Foundations

### 2.1 Basis set: variationally optimized pseudo-atomic orbitals (PAOs)
- Basis functions are numerical, strictly localized (finite-support) pseudo-atomic orbitals, generated as eigenstates of a confined atomic Kohn–Sham problem and then variationally optimized (via a force-theorem-based procedure) to minimize the total energy for representative reference systems.
- Basis quality is systematically improved by increasing the number of radial functions per angular momentum channel (e.g., "s2p2d1"), analogous to multiple-zeta plus polarization Gaussian/Slater basis sets.
- A comprehensive, systematic basis library spans elements from H to Kr and beyond, with cutoff radii and confinement schemes tuned per element and pre-optimized/benchmarked sets ("Standard", "Precise") distributed with the code.

### 2.2 Pseudopotentials
- Norm-conserving pseudopotentials in the Bachelet–Hamann–Schlüter / Troullier–Martins / Kleinman–Bylander separable forms, plus fully and scalar relativistic treatments (Morrison–Bylander–Kleinman-type forms for relativistic/SOC-including cases).
- Efficient projector expansion techniques reduce the cost of evaluating non-local pseudopotential projectors in the LCAO representation.

### 2.3 Exchange-correlation and extensions
- LDA/LSDA and GGA (PBE) exchange-correlation functionals.
- DFT+U for correlated d/f-electron systems, formulated on the non-orthogonal PAO basis.
- Non-collinear DFT (including generalized Bloch theorem for spin-spiral states) and constrained non-collinear spin/orbital DFT for complex magnetic structures.
- DFT-D2/D3 dispersion (van der Waals) corrections.

### 2.4 Numerical methods
- **Real-space grid** techniques for numerical integration and for solving the Poisson equation via FFT (energy-cutoff-controlled grid density).
- **Eigensolvers:** conventional dense diagonalization (parallel ELPA-based solver and ScaLAPACK, scaling to thousands of cores) for cluster/band-structure calculations of up to roughly 1000-atom systems.
- **O(N) / linear- and low-order-scaling methods** for much larger systems (10,000+ atoms with sufficient CPU resources):
  - **O(N) divide-and-conquer (DC) method**
  - **O(N) divide-and-conquer with localized natural orbitals (DC-LNO)** — a more recent, more efficient variant
  - **O(N) Krylov subspace method**
  - **Numerically exact low-order-scaling method** (contour-integration based)
- **Charge (density) mixing schemes:** simple mixing, Kerker preconditioning, Pulay/GR-Pulay, RMM-DIIS, and RMM-DIIS variants applied directly to the Hamiltonian matrix.
- **Geometry/cell optimization:** structural relaxation and variable-cell optimization (with and without constraints), and the nudged elastic band (NEB) method for reaction pathways/transition states.
- **Molecular dynamics:** NVE and NVT (velocity scaling, Nosé–Hoover) ensembles.
- **GPU acceleration:** more recent work offloads matrix multiplication and eigenvalue-solve steps of collinear and non-collinear DFT to GPUs (cuBLAS/cuSOLVER, OpenACC), giving substantial speedups on modern GPU-equipped clusters.

### 2.5 Property calculations and post-processing
- Electric/nonequilibrium transport via the non-equilibrium Green's function (NEGF) method for two-probe device geometries under finite bias.
- Maximally localized Wannier function construction, plus interfaces to Wannier90.
- Berry-phase macroscopic polarization, Z₂ invariant/Chern number/Berry curvature analysis, spin-texture analysis.
- Magnetic exchange coupling parameters (Liechtenstein-type mapping), magnetic anisotropy energy (second-variational SOC treatment).
- Band unfolding for supercell-to-primitive-cell band structure recovery.
- Absolute core-level binding energies (comparable to XPS), using fully relativistic pseudopotentials and a variational penalty functional with exact Coulomb cutoff.
- Density of states (DOS)/projected DOS, optical conductivity and dielectric function, effective screening medium (ESM) method for surfaces/electrochemical interfaces, scanning tunneling microscopy (STM) simulation via Tersoff–Hamann theory.
- Population/charge analysis: Mulliken, Voronoi, ESP fitting, and natural population/natural atomic orbital analysis.
- Interfaces to BoltzTrap (transport coefficients), Wannier90, XCrySDen, and its own OpenMX Viewer visualization tool.

---

## 3. Parallelism and Scalability

- Hybrid **MPI/OpenMP** parallelism is used throughout, essential to control the fairly high per-rank memory footprint of the code.
- A **three-dimensional domain decomposition** scheme for large-scale calculations improves scalability of the linear-scaling Krylov subspace method, demonstrated on machines up to the K computer (>100,000 cores, systems of ~262,000 atoms).
- Conventional (diagonalization-based) calculations scale efficiently to roughly 1000-atom systems on machines with hundreds to a few thousand cores; O(N) methods extend applicability toward 10,000+ atoms given sufficient core counts.
- Independent benchmarking (e.g., by Peter Larsson, NSC) confirms that OpenMX gives accurate delta-gauge results comparable to plane-wave codes, but notes that memory replication across MPI ranks makes the MPI/OpenMP hybrid mode important for practical large jobs.

---

## 4. Typical Application Domains

OpenMX is widely used across:
- Bulk crystalline solids, surfaces, interfaces, and low-dimensional/2D materials
- Molecular and cluster systems, including biomolecules (proteins, DNA, polysaccharides) via pre-optimized basis sets tuned for biological elements
- Magnetic materials — collinear and non-collinear magnetism, spin-orbit coupling, magnetic anisotropy, exchange coupling constants, spin-spiral structures
- Electronic transport / two-probe nanodevices via the NEGF-DFT method
- Strongly correlated materials via DFT+U
- Liquids and disordered/amorphous systems, often in molecular dynamics workflows
- High-throughput and machine-learning training-data generation pipelines (e.g., via ASE's OpenMX calculator interface), owing to its efficient, well-parametrized basis/pseudopotential databases

---

## 5. Software Ecosystem

- **License:** GPLv3, free and open source.
- **Language:** written in C (with Fortran components), MPI/OpenMP parallel.
- **Data:** curated VPS (pseudopotential) and PAO (basis orbital) databases for most of the periodic table, generated with the companion `ADPACK` tool and benchmarked using the Δ-gauge protocol against reference all-electron/plane-wave results.
- **Interfaces:** Atomic Simulation Environment (ASE) calculator, Wannier90, BoltzTraP, XCrySDen, VESTA-compatible output formats, and its own OpenMX Viewer.
- **Availability:** packaged on numerous HPC clusters (e.g., Tetralith/NSC) and in scientific Linux distributions/Copr repositories, in addition to source-only distribution from the official site.

---

## 6. Strengths and Limitations

**Strengths**
- Strictly localized basis and pseudopotentials give sparse matrices, enabling both efficient cluster/band diagonalization and true O(N) scaling for very large systems.
- Extensive, well-validated pseudopotential/PAO database lowers the barrier to starting new calculations.
- Broad functionality: magnetism (collinear/non-collinear/SOC), transport (NEGF), strong correlation (DFT+U), and Berry-phase/topological property calculations are all natively supported.
- Fully open source (GPLv3), in contrast to several other localized-basis O(N) DFT codes.
- Actively maintained, with continuing additions (DC-LNO, spin-spiral generalized Bloch theorem, GPU acceleration, XPS binding energies) into recent versions (3.9.x).

**Limitations**
- Finite, localized basis sets require basis-convergence testing per system/property, unlike plane-wave (fully systematic) basis sets.
- Documented high per-MPI-rank memory replication can make naive flat-MPI runs on large jobs memory-limited; hybrid MPI/OpenMP configuration is effectively required for large systems.
- As with other LCAO/O(N) codes, achieving both linear scaling and high numerical accuracy simultaneously requires careful choice of cutoff radii, Krylov subspace dimensions, or DC/DC-LNO buffer regions — these are non-trivial convergence parameters.

---

## 7. Publications Related to OpenMX Theory and Methodology

The following list follows the official "Related Papers" citation guide maintained by the OpenMX developers, organized by functionality, supplemented with foundational methodology references.

### General / Core Method
- T. Ozaki, *Variationally optimized atomic orbitals for large-scale electronic structures*, **Phys. Rev. B 67**, 155108 (2003).
- T. Ozaki and H. Kino, *Numerical atomic basis orbitals from H to Kr*, **Phys. Rev. B 69**, 195113 (2004).
- T. Ozaki and H. Kino, *Efficient projector expansion for the ab initio LCAO method*, **Phys. Rev. B 72**, 045121 (2005).
- T. Ozaki and H. Kino, *Variationally optimized basis orbitals for biological molecules*, **J. Chem. Phys. 121**, 10879 (2004).
- K. Lejaeghere et al., *Reproducibility in density functional theory calculations of solids*, **Science 351**, aad3000 (2016). *(Δ-gauge benchmark validating OpenMX among other codes)*

### Large-Scale Parallel Calculation
- T. V. T. Duy and T. Ozaki, *A three-dimensional domain decomposition method for large-scale DFT electronic structure calculations*, **Comput. Phys. Commun. 185**, 777 (2014).
- T. V. T. Duy and T. Ozaki, *A decomposition method with minimum communication amount for parallelization of multi-dimensional FFTs*, **Comput. Phys. Commun. 185**, 153 (2014).

### O(N) Divide-and-Conquer with Localized Natural Orbitals (DC-LNO)
- T. Ozaki, M. Fukuda, and G. Jiang, *Efficient O(N) method with localized natural orbitals*, **Phys. Rev. B 98**, 245137 (2018).

### O(N) Krylov Subspace Method
- T. Ozaki, *O(N) Krylov-subspace method for large-scale ab initio electronic structure calculations*, **Phys. Rev. B 74**, 245101 (2006).

### Numerically Exact Low-Order Scaling Method
- T. Ozaki, *Continued fraction representation of the Fermi–Dirac function for large-scale electronic structure calculations*, **Phys. Rev. B 82**, 075131 (2010).

### DFT+U Method
- M. J. Han, T. Ozaki, and J. Yu, *O(N) LDA+U electronic structure calculation method based on the nonorthogonal pseudoatomic orbital basis*, **Phys. Rev. B 73**, 045110 (2006).
- S. Ryee and M. J. Han, *The effect of double counting, spin density, and Hund interaction in the different DFT+U functionals*, **J. Phys.: Condens. Matter 30**, 275802 (2018).
- S. Ryee and M. J. Han, *SolvedU: SOLar Values of the effective on-site Coulomb interaction using the constrained DFT method*, **Scientific Reports 8**, 9559 (2018).

### Exchange Coupling Parameter
- M. J. Han, T. Ozaki, and J. Yu, *Exchange coupling and magnetic anisotropy in Fe nanowires*, **Phys. Rev. B 70**, 184421 (2004).
- A. Terasawa, M. Matsumoto, T. Ozaki, and Y. Gohda, *Efficient algorithm based on Liechtenstein method for computing exchange coupling constants*, **J. Phys. Soc. Jpn. 88**, 114706 (2019).

### NEGF (Electron Transport) Method
- T. Ozaki, K. Nishio, and H. Kino, *Efficient implementation of the nonequilibrium Green function method for electronic transport calculations*, **Phys. Rev. B 81**, 035116 (2010).
- T. Ozaki, *Efficient recursion method for inverting an overlap matrix*, **Phys. Rev. B 75**, 035123 (2007).

### Effective Screening Medium (ESM) Method
- T. Ohwaki, M. Otani, T. Ikeshoji, and T. Ozaki, *Selective formation of hydrogen and hydronium ions*, **J. Chem. Phys. 136**, 134101 (2012).

### Generation of Wannier Functions
- H. Weng, T. Ozaki, and K. Terakura, *Revisiting magnetic coupling in transition-metal-benzene complexes with maximally localized Wannier functions*, **Phys. Rev. B 79**, 235118 (2009).

### Generation of Natural Atomic Orbitals
- T. Ohwaki, M. Otani, and T. Ozaki, *A method to construct isolated orbitals for large systems using locally constructed Wannier functions*, **J. Chem. Phys. 140**, 244105 (2014).

### Band Unfolding Method
- C.-C. Lee, Y. Yamada-Takamura, and T. Ozaki, *Unfolding method for first-principles LCAO electronic structure calculations*, **J. Phys.: Condens. Matter 25**, 345501 (2013).

### XPS Core-Level Binding Energies
- T. Ozaki and C.-C. Lee, *Absolute binding energies of core levels in solids from first principles*, **Phys. Rev. Lett. 118**, 026401 (2017).

### BoltzTraP Interface (Transport Coefficients)
- M. Miyata, T. Ozaki, T. Takeuchi, S. Nishino, M. Inukai, and M. Koyano, *High-throughput screening of Sb-based thermoelectric materials using calculations of transport properties*, **J. Electronic Materials 47**, 3254 (2017).

### FermiSurfer
- M. Kawamura, *FermiSurfer: Fermi-surface viewer providing multiple representation schemes*, **Comput. Phys. Commun. 239**, 197 (2019).

### Spin Texture Analysis in k-Space
- H. Kotaka, F. Ishii, and M. Saito, *Rashba effect on thin bismuth film from first-principles*, **Jpn. J. Appl. Phys. 52**, 035204 (2013).
- N. Yamaguchi and F. Ishii, *Spin-orbit splitting in the vacuum region of surface systems*, **Appl. Phys. Express 10**, 123003 (2017).

### Spin Spiral Calculations
- T. B. Prayitno and F. Ishii, *First-principles calculation of the generalized Bloch theorem in the presence of spin-orbit coupling*, **J. Phys. Soc. Jpn. 87**, 114709 (2018).
- T. B. Prayitno and F. Ishii, *Effect of the Hubbard U on the generalized Bloch theorem*, **J. Phys. Soc. Jpn. 88**, 054701 (2019).

### Z₂ Invariant, Chern Number, and Berry Curvature
- H. Sawahata, N. Yamaguchi, H. Kotaka, and F. Ishii, *First-principles study of the Z₂ topological invariants*, **Jpn. J. Appl. Phys. 57**, 030309 (2018).

### OpenMX Viewer
- Y.-T. Lee and T. Ozaki, *OpenMX Viewer: An open-source project for scientific data visualization*, **J. Molecular Graphics and Modelling 89**, 192 (2019).

### Additional Foundational and Methodological References (widely cited alongside the above)
- G. B. Bachelet, D. R. Hamann, and M. Schlüter, *Pseudopotentials that work: From H to Pu*, **Phys. Rev. B 26**, 4199 (1982). *(norm-conserving pseudopotential formalism)*
- N. Troullier and J. L. Martins, *Efficient pseudopotentials for plane-wave calculations*, **Phys. Rev. B 43**, 1993 (1991).
- L. Kleinman and D. M. Bylander, *Efficacious form for model pseudopotentials*, **Phys. Rev. Lett. 48**, 1425 (1982).
- P. E. Blöchl, *Generalized separable potentials for electronic-structure calculations*, **Phys. Rev. B 41**, 5414 (1990).
- A. Morrison, D. M. Bylander, and L. Kleinman, *Nonlocal Hermitian norm-conserving Vanderbilt pseudopotential* (MBK-type relativistic pseudopotential formalism).
- J. P. Perdew, K. Burke, and M. Ernzerhof, *Generalized gradient approximation made simple*, **Phys. Rev. Lett. 77**, 3865 (1996).
- S. Grimme, *Semiempirical GGA-type density functional constructed with a long-range dispersion correction*, **J. Comput. Chem. 27**, 1787 (2006). *(DFT-D2)*
- S. Grimme, J. Antony, S. Ehrlich, and H. Krieg, *A consistent and accurate ab initio parametrization of density functional dispersion correction (DFT-D) for the 94 elements H–Pu*, **J. Chem. Phys. 132**, 154104 (2010). *(DFT-D3)*

---

## 8. Key References for Further Reading

- Official OpenMX website and documentation: <https://www.openmx-square.org/>
- User's Manual (versioned, e.g. Ver. 3.9): openmx-square.org/openmx_man3.9/
- Pseudopotential and PAO database with Δ-gauge validation: openmx-square.org (VPS/PAO database pages)
- ASE OpenMX calculator documentation: <https://wiki.fysik.dtu.dk/ase/ase/calculators/openmx.html>

---

*This review is a synthesis of publicly available OpenMX documentation, the official OpenMX "Related Papers" citation list, and peer-reviewed literature describing and applying the code. It is intended as a technical orientation document, not an exhaustive line-by-line reproduction of the OpenMX manual.*

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the OpenMX 	Open-source DFT code using pseudo-atomic localized basis functions and norm-conserving pseudopotentials, efficient for large systems. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
