# Density Functional Theory (DFT) Software Packages

An exhaustive reference list of software packages that implement Density Functional Theory, organized by primary basis-set/methodology family. Packages are open-source or free-for-academic-use unless noted as commercial.

---

## Plane-Wave / Pseudopotential Codes

| Package | Description |
|---|---|
| **VASP** (Vienna Ab initio Simulation Package) | Widely used commercial plane-wave PAW/pseudopotential code for periodic solids, surfaces, and molecules; a mainstay in materials science and catalysis. |
| **Quantum ESPRESSO** | Open-source integrated suite (PWscf, CP, PHonon, etc.) for plane-wave pseudopotential DFT, phonons, and molecular dynamics. |
| **ABINIT** | Open-source plane-wave/PAW code with strong capabilities in DFPT (density-functional perturbation theory), phonons, and many-body perturbation theory (GW, BSE). |
| **CASTEP** | Plane-wave pseudopotential code for solids, surfaces, and interfaces, distributed with/without BIOVIA Materials Studio; commercial with academic licensing. |
| **CPMD** (Car-Parrinello Molecular Dynamics) | Plane-wave code specialized in Car-Parrinello and Born-Oppenheimer ab initio molecular dynamics. |
| **CP2K** | Open-source code using a mixed Gaussian and plane-wave (GPW) approach, excelling at large-scale condensed-phase and QM/MM simulations. |
| **PWmat** | GPU-accelerated plane-wave pseudopotential DFT code designed for high-throughput and machine-learning-assisted materials simulation. |
| **JDFTx** | Open-source plane-wave DFT code emphasizing joint density-functional theory for solvation and electrochemical interfaces. |
| **Qbox** | Scalable plane-wave/pseudopotential code designed for large-scale parallel first-principles molecular dynamics on supercomputers. |
| **PARSEC** | Real-space pseudopotential DFT code using finite differences on a grid rather than plane waves. |
| **Octopus** | Real-space, real-time DFT/TDDFT code for finite and periodic systems, focused on excited-state and spectroscopic properties. |
| **BigDFT** | DFT code using Daubechies wavelets as a basis, enabling adaptive resolution and linear-scaling algorithms. |
| **ONETEP** | Linear-scaling plane-wave-based DFT code using local orbitals ("psinc" functions), designed for very large systems. |
| **CONQUEST** | Linear-scaling DFT code for simulations of very large systems (up to millions of atoms) using numerical atomic orbitals or blip basis functions. |
| **SIESTA** | Open-source, order-N (linear-scaling capable) DFT code using numerical atomic orbital basis sets and norm-conserving pseudopotentials. |
| **OpenMX** | Open-source DFT code using pseudo-atomic localized basis functions and norm-conserving pseudopotentials, efficient for large systems. |
| **FHI-aims** | All-electron, numeric atom-centered orbital DFT code offering high accuracy across molecules, clusters, and periodic solids. |
| **DFTB+** | Fast approximate DFT code implementing density-functional tight-binding (DFTB) and extended tight-binding (xTB) methods for large systems. |
| **DMol3** | Numerical atomic orbital DFT code for molecules, surfaces, and solids, distributed within BIOVIA Materials Studio (commercial). |
| **GPAW** | Real-space/plane-wave/LCAO DFT code implementing the projector augmented-wave (PAW) method, built on the Python ASE framework. |
| **PROFESS** | Orbital-free DFT code for large-scale simulations using kinetic-energy density functionals instead of Kohn-Sham orbitals. |
| **DFTK.jl** (Density-Functional Toolkit) | Julia-based, plane-wave pseudopotential DFT code emphasizing algorithmic research, verification, and mathematical rigor. |
| **RESCU** | GPU-accelerated linear combination of atomic orbitals (LCAO) and plane-wave DFT code for large-scale materials simulations. |
| **SPARC** | Real-space finite-difference DFT code designed for accuracy and scalability on modern high-performance computing architectures. |
| **PWDFT / PWDFT.jl** | Plane-wave DFT codes (C++/Julia implementations) developed for large-scale and GPU-accelerated electronic structure calculations. |
| **Abacus** | Open-source DFT package supporting both plane-wave and numerical atomic orbital bases, developed for large-scale and high-throughput calculations. |
| **DGDFT** | Discontinuous Galerkin DFT code designed for massively parallel, large-scale electronic structure calculations. |
| **HONPAS** | Numerical atomic orbital DFT code based on the SIESTA methodology with hybrid functional support. |

---

## All-Electron / Augmented Plane-Wave & LMTO Codes

| Package | Description |
|---|---|
| **WIEN2k** | Widely used all-electron full-potential linearized augmented plane-wave plus local orbitals (FP-LAPW+lo) code for solids; commercial/academic license. |
| **Elk** | Open-source, full-potential all-electron LAPW code, notable for extensive support of advanced DFT and beyond-DFT methods. |
| **FLEUR** | Open-source full-potential linearized augmented plane-wave (FLAPW) code, strong in magnetism and spin-orbit coupling studies. |
| **Exciting** | Full-potential all-electron LAPW+lo code with a strong focus on excited-state properties (GW, BSE, TDDFT). |
| **FHI96 / FLAIR** | All-electron LAPW-based codes with historical roots in surface and thin-film electronic structure studies. |
| **RSPt** (Relativistic Spin Polarized toolkit) | Full-potential linear muffin-tin orbital (FP-LMTO) all-electron DFT code for solids. |
| **Questaal (LMTO/ASA)** | Suite of full-potential and atomic-sphere-approximation LMTO codes for electronic structure, including DFT+DMFT and GW. |
| **AkaiKKR (MACHIKANEYAMA)** | All-electron Korringa-Kohn-Rostoker (KKR) Green's-function DFT code, well suited to disordered alloys via the coherent potential approximation. |
| **SPR-KKR** | Fully relativistic spin-polarized KKR Green's-function DFT code for solids, surfaces, and disordered systems. |

---

## Quantum Chemistry Codes (Molecular, Gaussian-Type Orbitals)

| Package | Description |
|---|---|
| **Gaussian** | Widely used commercial quantum chemistry package supporting a vast range of DFT functionals plus HF and post-HF methods for molecules. |
| **ORCA** | Free-for-academic-use quantum chemistry package popular for DFT, TD-DFT, and correlated wavefunction calculations, known for its efficiency and usability. |
| **TURBOMOLE** | Commercial/academic quantum chemistry package known for fast, efficient DFT and post-HF calculations via RI/density-fitting techniques; also supports periodic DFT (Riper module). |
| **NWChem** | Open-source, highly scalable quantum chemistry suite supporting DFT, TDDFT, and periodic (plane-wave and Gaussian) calculations on HPC systems. |
| **Psi4** | Open-source quantum chemistry package with extensive DFT and post-HF methods, popular for method development and Python scripting. |
| **Q-Chem** | Commercial quantum chemistry package offering a broad range of DFT functionals plus advanced excited-state and solvation methods. |
| **GAMESS (US)** | Free general quantum chemistry package supporting DFT, HF, and many post-HF methods, including QM/MM via the QuanPol module. |
| **GAMESS (UK)** | Related but independently developed quantum chemistry package with DFT and ab initio capabilities, historically UK-centered. |
| **MOLPRO** | Commercial quantum chemistry package strong in high-accuracy correlated methods, also offering DFT functionality. |
| **Molcas / OpenMolcas** | Quantum chemistry package specializing in multiconfigurational methods, with DFT capabilities included; OpenMolcas is the open-source branch. |
| **DALTON** | Open-source quantum chemistry package with a strong focus on molecular properties and response theory, including DFT/TDDFT. |
| **deMon2k** | DFT-focused quantum chemistry code using Gaussian-type auxiliary density functional methods, efficient for large molecular systems. |
| **DeMonNano** | DFT-based tight-binding code (auxiliary density functional tight-binding) for large biomolecular and nanoscale systems. |
| **Firefly (PC GAMESS)** | Quantum chemistry package derived from GAMESS(US), optimized for x86 architectures, supporting DFT calculations. |
| **Jaguar** | Commercial quantum chemistry package (Schrödinger suite) offering fast DFT methods for pharma and materials applications. |
| **deal.II-based codes / DFTFE** | Real-space finite-element DFT code (DFT-FE) designed for large-scale, high-accuracy all-electron and pseudopotential calculations on HPC/GPU systems. |
| **MPQC** (Massively Parallel Quantum Chemistry) | Open-source quantum chemistry package designed for parallel computation, including DFT methods. |
| **PySCF** | Open-source, Python-based quantum chemistry package supporting DFT, post-HF, and periodic calculations, popular for method prototyping. |
| **CFOUR** | Quantum chemistry package specializing in coupled-cluster methods, also including DFT functionality. |
| **ADF** (Amsterdam Density Functional) | Commercial DFT-centric package using Slater-type orbitals, part of the Amsterdam Modeling Suite (AMS), strong in relativistic and spectroscopic calculations. |
| **BAND** | Periodic-system DFT code (Slater-type orbitals) within the Amsterdam Modeling Suite, for solids, surfaces, and semiconductors. |
| **deMon2k / StoBe** | Gaussian-basis DFT codes historically used for X-ray absorption spectroscopy simulations. |
| **Spartan** | Commercial molecular modeling package offering DFT alongside HF/semi-empirical methods, aimed at broad usability. |
| **HyperChem** | Commercial molecular modeling software including basic DFT functionality alongside molecular mechanics and semi-empirical methods. |
| **Turbomole Riper** | Periodic-DFT extension module within TURBOMOLE using Gaussian-type orbitals for solids and surfaces. |
| **MOPAC** | Primarily semi-empirical quantum chemistry package; re-released as open source, sometimes used alongside DFT workflows for large systems. |
| **xtb (Grimme group)** | Open-source semi-empirical extended tight-binding package built to approximate DFT-level accuracy at very low computational cost. |
| **ACES / ACES III** | Quantum chemistry package specializing in coupled-cluster and many-body methods, with parallel implementations, also supporting DFT. |
| **Molden / GaussSum** (auxiliary) | Not DFT engines themselves but common companion tools for visualizing/post-processing DFT output (included for completeness of the ecosystem). |

---

## Linear-Scaling / Order-N and Fragment-Based Codes

| Package | Description |
|---|---|
| **ONETEP** | *(see above)* Linear-scaling plane-wave-quality DFT for very large systems. |
| **CONQUEST** | *(see above)* Order-N DFT for million-atom-scale simulations. |
| **FreeON** | Open-source, linear-scaling quantum chemistry/DFT code designed for very large molecular and periodic systems. |
| **ERGOSCF (Ergo)** | Open-source linear-scaling HF/DFT code focused on very large molecular systems using density matrix purification techniques. |

---

## Real-Time / Excited-State & TDDFT-Focused Codes

| Package | Description |
|---|---|
| **Octopus** | *(see above)* Real-time TDDFT and ground-state DFT for spectroscopy and non-equilibrium dynamics. |
| **NWChem TDDFT module** | Time-dependent DFT capability within NWChem for excited-state properties. |
| **SALMON** | Open-source real-time TDDFT code for light-matter interaction simulations in molecules and solids. |
| **Yambo** | Many-body perturbation theory code (GW, BSE) built to work alongside ground-state DFT codes (e.g., Quantum ESPRESSO, Abinit) for excited states. |
| **BerkeleyGW** | Many-body perturbation theory package (GW/BSE) that uses mean-field DFT wavefunctions from codes like Quantum ESPRESSO, PARATEC, or Abinit as input. |

---

## Orbital-Free / Tight-Binding / Semi-Empirical DFT Approximations

| Package | Description |
|---|---|
| **DFTB+** | *(see above)* Density-functional tight-binding package. |
| **xtb** | *(see above)* Semi-empirical extended tight-binding methods. |
| **PROFESS** | *(see above)* Orbital-free DFT. |
| **DFTB (original code)** | Predecessor implementation of the density-functional tight-binding method, foundational to DFTB+. |
| **hotbit** | Open-source density-functional tight-binding code implemented in Python, integrated with ASE. |

---

## Materials-Informatics / High-Throughput & Workflow Frameworks (DFT-Adjacent)

*(Not DFT engines themselves, but widely used to drive DFT codes for high-throughput or automated workflows — included for completeness.)*

| Package | Description |
|---|---|
| **AiiDA** | Open-source workflow and provenance-tracking framework for automating and managing high-throughput DFT calculations across many codes. |
| **ASE** (Atomic Simulation Environment) | Python library providing a common interface to set up, run, and analyze calculations from dozens of DFT codes. |
| **Pymatgen / Atomate** | Python materials analysis library and associated high-throughput DFT workflow framework, closely tied to the Materials Project. |
| **DFTTK** (DFT ToolKit) | Python-based high-throughput toolkit built on Atomate for finite-temperature thermodynamic properties from DFT phonon calculations. |
| **VASPilot** | Multi-agent, LLM-assisted automation framework for running and managing VASP DFT simulations. |

---

## Notes

- **Commercial vs. open-source** varies by package and sometimes by module; check current vendor/project pages for licensing terms.
- Several codes overlap categories (e.g., CP2K and TURBOMOLE support both molecular and periodic calculations); placement above reflects their most characteristic use.
- For an actively maintained, continuously updated master list, see Wikipedia's *"List of quantum chemistry and solid-state physics software"* and *"Category: Density functional theory software."*

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive list Density Functional Theory (DFT) Software Packages. For each package provide the short description of the package. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
