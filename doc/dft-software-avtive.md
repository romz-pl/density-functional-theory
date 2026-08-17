# Actively Developed DFT Software Packages

An exhaustive reference list of density functional theory (DFT) software packages currently under active development, organized by primary basis-set/methodology category. "Actively developed" here means the package has had a release, commit, or documented update within the last ~1–2 years.

---

## Plane-Wave / Pseudopotential Codes (Periodic Systems)

| Package | Description |
|---|---|
| **VASP** (Vienna Ab initio Simulation Package) | Commercial plane-wave DFT code using PAW/ultrasoft pseudopotentials; the most widely used code in condensed matter/materials science for periodic solids and surfaces. |
| **Quantum ESPRESSO** | Open-source (GPL) integrated suite (PWscf, CP, PHonon, etc.) for plane-wave pseudopotential and PAW DFT, DFPT, and ab initio molecular dynamics. |
| **ABINIT** | Open-source plane-wave DFT/DFPT package with strong support for response-function properties (phonons, dielectric tensors, elastic constants) and PAW. |
| **CASTEP** | Academic/commercial plane-wave pseudopotential code widely used for solid-state materials, spectroscopy, and NMR properties. |
| **CPMD** | Plane-wave/pseudopotential code specialized in Car–Parrinello and Born–Oppenheimer ab initio molecular dynamics. |
| **Qbox** | Scalable plane-wave DFT code designed for large-scale parallel ab initio molecular dynamics on supercomputers. |
| **JDFTx** | Plane-wave DFT code (successor to DFT++) with strong support for solvation models and electrochemistry/joint density-functional theory. |
| **Quantum Espresso GPU / QE-GPU** | GPU-accelerated fork/extension of Quantum ESPRESSO. |
| **PWmat** | GPU-accelerated commercial plane-wave DFT code developed in China, used for large-scale materials simulation. |
| **SPARC** | Real-space and plane-wave finite-difference DFT code emphasizing high-performance, large-scale, GPU-accelerated calculations. |
| **INQ** | Modern, GPU-native (time-dependent) DFT framework built for exascale/GPU-first architectures. |
| **DFTpy** | Python-based orbital-free DFT (OF-DFT) code built on a plane-wave grid. |

## All-Electron / Full-Potential Codes

| Package | Description |
|---|---|
| **FHI-aims** | All-electron, full-potential code using numeric atom-centered orbitals (NAOs); supports molecules and periodic solids with high accuracy, plus GW/RPA post-DFT methods. |
| **WIEN2k** | Full-potential linearized augmented plane-wave (FLAPW) code, a long-standing standard for high-accuracy all-electron solid-state DFT. |
| **exciting** | Full-potential all-electron FLAPW/LAPW+lo code with strong excited-state (GW, BSE) capabilities. |
| **FLEUR** | Full-potential linearized augmented plane-wave (FLAPW) code from the Jülich group, strong in magnetism and spin-orbit coupling. |
| **Elk** | Open-source full-potential LAPW code with extensive support for advanced DFT extensions including DFT+U, hybrid functionals, and TDDFT. |
| **RSPt** | Full-potential linear muffin-tin orbital (LMTO) code for electronic structure of solids. |
| **FPLO** | Full-potential local-orbital minimum-basis code for electronic structure of solids. |

## Gaussian-Type / Atomic-Orbital Codes

| Package | Description |
|---|---|
| **Gaussian** | Widely used commercial quantum chemistry package with comprehensive DFT functionality for molecules; a long-standing industry/academic standard. |
| **ORCA** | Free-for-academic-use quantum chemistry package with extensive DFT (including DLPNO coupled-cluster and other post-DFT correlation methods), popular for molecular and moderate periodic systems. |
| **Q-Chem** | Commercial quantum chemistry package with broad DFT functional coverage and advanced excited-state/response-theory methods. |
| **TURBOMOLE** | Commercial/academic package known for high-performance molecular DFT (RI approximations) with growing periodic (Riper module) capability. |
| **NWChem** | Open-source, highly scalable quantum chemistry suite for both molecular and periodic (plane-wave) DFT on large HPC systems. |
| **GAMESS (US)** | Free general ab initio/DFT quantum chemistry package for molecular systems, developed at Iowa State. |
| **Psi4** | Open-source quantum chemistry package with modern DFT and post-HF methods, popular in education and Python-integrated workflows. |
| **PySCF** | Python-based, modular open-source package for molecular and periodic (PBC) DFT/post-HF methods, widely used as a research and development platform. |
| **CRYSTAL** (CRYSTAL23) | Gaussian-basis periodic ab initio code specializing in crystalline solids, surfaces, and polymers with hybrid functional support. |
| **deMon2k** | Gaussian-basis DFT code (auxiliary density functional theory, ADFT) for molecular systems. |
| **MOLGW** | Open-source many-body Green's function code (GW/BSE) built on a Gaussian-basis DFT foundation. |
| **ERKALE** | Open-source Gaussian-basis DFT/HF code, notable for X-ray spectroscopy and finite-element extensions. |
| **deMon-Nano** | DFT-based tight-binding/DFT package extension for large systems and nanostructures. |
| **BAGEL** | Open-source package for accurate multireference and relativistic electronic structure, including relativistic DFT. |
| **OpenMolcas** | Open-source multireference quantum chemistry package with DFT and multiconfigurational methods. |
| **CFOUR** | High-accuracy coupled-cluster/DFT quantum chemistry package for molecular properties. |
| **MRCC** | High-accuracy quantum chemistry code with DFT and systematically improvable correlation methods. |
| **Dalton/LSDalton** | Open-source molecular electronic structure program with strong response-theory (spectroscopic property) DFT capabilities. |
| **Firefly (PC GAMESS)** | Free quantum chemistry package (GAMESS-US derivative) optimized for x86 architectures. |
| **HORTON** | Python-based electronic structure development platform, used for DFT method development and density-based analysis. |
| **DIRAC** | Relativistic (four-component) quantum chemistry code including relativistic DFT. |

## Numerical Atomic Orbital / Localized-Basis Codes

| Package | Description |
|---|---|
| **SIESTA** | Open-source DFT code using numerical atomic orbitals and norm-conserving pseudopotentials, optimized for large-scale (linear-scaling) simulations. |
| **OpenMX** | Open-source DFT code using pseudo-atomic localized basis functions, efficient for large-scale materials and nanostructure simulations. |
| **CP2K** | Open-source package combining Gaussian and plane-wave (GPW/GAPW) methods for DFT, ab initio MD, and QM/MM, scalable to large periodic systems. |
| **ONETEP** | Linear-scaling DFT code using optimized local orbitals, designed for very large systems (thousands of atoms) on HPC systems. |
| **CONQUEST** | Linear-scaling DFT code for very large-scale (million-atom) simulations using local orbitals. |
| **BigDFT** | Open-source DFT code using Daubechies wavelets as basis functions, supporting linear scaling and free/periodic boundary conditions. |
| **RMG** (Real-space Multigrid) | Open-source, GPU-accelerated real-space multigrid DFT code for large-scale materials and molecular simulation. |
| **DGDFT** | Discontinuous Galerkin density functional theory code for large-scale, massively parallel electronic structure calculations. |

## Real-Space Grid / Finite-Element Codes

| Package | Description |
|---|---|
| **Octopus** | Open-source real-space grid TDDFT/DFT code specializing in excited-state and time-dependent phenomena in finite and periodic systems. |
| **GPAW** | Open-source Python-based DFT code using the projector-augmented wave (PAW) method, with real-space grid, plane-wave, and LCAO modes. |
| **HelFEM** | Finite-element method code for high-accuracy DFT/HF calculations on atoms and diatomic molecules. |
| **DFTK** (Density-Functional ToolKit) | Modern, open-source Julia-based plane-wave DFT code focused on flexibility, verification, and method-development research. |

## Tight-Binding / Semi-Empirical DFT-Based Codes

| Package | Description |
|---|---|
| **DFTB+** | Open-source density functional tight-binding (DFTB) code enabling fast approximate DFT-quality simulations of large systems. |
| **xtb** | Open-source semi-empirical extended tight-binding package (GFN-xTB family) from the Grimme group, widely used for fast conformer/geometry screening. |

## Orbital-Free and Embedding / Multiscale DFT

| Package | Description |
|---|---|
| **PROFESS** | Orbital-free DFT code for large-scale simulations using kinetic energy density functionals in place of explicit orbitals. |
| **eQE** (embedded Quantum ESPRESSO) | Subsystem/embedding DFT extension of Quantum ESPRESSO for multiscale simulations. |

## Machine-Learning-Accelerated / Hybrid DFT Frameworks

| Package | Description |
|---|---|
| **DeePKS / DeePMD-kit (DFT-coupled workflows)** | ML-based frameworks that learn DFT-quality potential energy surfaces/corrections, actively developed for accelerating and extending DFT-level accuracy to large systems. |
| **ABACUS** | Open-source Chinese-developed plane-wave and numerical-atomic-orbital DFT code with growing GPU and machine-learning integration. |

## Commercial/Industrial Materials Modeling Suites (DFT engines included)

| Package | Description |
|---|---|
| **Materials Studio / DMol3 / CASTEP** (BIOVIA) | Commercial materials modeling suite bundling multiple DFT engines (DMol3 numerical-orbital code and CASTEP plane-wave code) for industrial and academic use. |
| **QuantumATK** (Synopsys) | Commercial atomistic simulation platform combining DFT, tight-binding, and NEGF transport methods for nanoelectronics and materials design. |
| **Amsterdam Modeling Suite / ADF** (SCM) | Commercial DFT suite (Amsterdam Density Functional, BAND for periodic systems) using Slater-type orbitals, strong in molecular property prediction. |
| **Jaguar** (Schrödinger) | Commercial quantum chemistry package with DFT, widely used in pharmaceutical and materials industry workflows. |
| **Spartan** | Commercial molecular modeling software with integrated DFT for educational and research use. |
| **TeraChem** | Commercial GPU-accelerated quantum chemistry/DFT package optimized for large molecular systems and ab initio MD. |

## High-Throughput / Workflow Frameworks (DFT orchestration, not DFT engines themselves)

| Package | Description |
|---|---|
| **ASE** (Atomic Simulation Environment) | Open-source Python library for scripting, driving, and analyzing DFT calculations across many backend codes. |
| **AiiDA** | Open-source workflow and data-provenance engine for automating and managing high-throughput DFT calculations. |
| **pymatgen / Atomate2** | Open-source Python materials-analysis library and workflow suite for automated high-throughput DFT (built around VASP and other codes). |
| **DFTTK** (DFT ToolKit) | Python package for high-throughput DFT-based lattice-dynamics and thermodynamic-property calculations. |

---

### Notes
- This list emphasizes packages with demonstrable recent activity (recent releases, active repositories, or recent papers/preprints as of 2025–2026); truly abandoned or historical-only codes (e.g., original DACAPO, CADPAC, AMPAC-era tools) were excluded.
- "DFT software" spans several distinct methodological families (plane-wave, all-electron, Gaussian-basis, real-space, tight-binding, orbital-free); the categorization above reflects the primary basis/representation each package uses, though many support multiple modes.
- Licensing varies widely — from fully open-source (GPL/MIT/Apache) to free-for-academic to fully commercial; license status changes over time and should be verified against each project's current site before use.


---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive list of all actively developed density functional theory (DFT) software packages. For each package provide the short description of the package. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
