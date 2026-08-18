# CPMD (Car–Parrinello Molecular Dynamics): An Exhaustive Review

## 1. Overview

CPMD is a plane-wave/pseudopotential implementation of Kohn–Sham Density Functional Theory (DFT), purpose-built for *ab initio* molecular dynamics (AIMD). Its defining feature is efficient support for the **Car–Parrinello (CP) extended-Lagrangian scheme**, in which the electronic wavefunction coefficients are propagated as classical dynamical variables alongside the nuclei, avoiding an explicit self-consistent-field (SCF) diagonalization at every MD step. CPMD also supports conventional **Born–Oppenheimer molecular dynamics (BOMD)**, in which the electronic ground state is fully converged at each nuclear configuration before the forces are computed.

The code originated in 1990 at IBM Research – Zurich, and between 1994 and 2001 was developed jointly by IBM and the Max Planck Institute for Solid State Research, Stuttgart. In 2001 a CPMD Consortium was formed to coordinate worldwide development and distribution, founded by Prof. Michele Parrinello (USI) and Prof. Wanda Andreoni (EPFL/IBM Zurich), and coordinated by Alessandro Curioni (IBM). Core contributors across its history include Jürg Hutter, Mauro Boero, Alessandro Curioni, Mark Tuckerman, Ivano Tavernelli, Axel Kohlmeyer, Nisanth Nair, Ari Paavo Seitsonen, Wolfram Quester, Dominik Marx, Valery Weber, Teodoro Laino, Takashi Ikeda, Ali Alavi, Michiel Sprik, Marcella Iannuzzi, and Jochen Blumberger, among others.

In **2022, IBM released the CPMD source code freely on GitHub under the MIT License** (organization: `CPMD-code`), ending its previous registration/license-agreement distribution model. The code is copyrighted 1990–2023 by IBM Corp. and 1994–2001 by the Max Planck Institute, Stuttgart. Following this release, a community initiative called **OpenCPMD** formed to continue collaborative development, documentation, and tutorials outside IBM.

CPMD is written primarily in **Fortran**, parallelized with **MPI** (and OpenMP/hybrid schemes in modern versions), and historically has been noted as one of the most HPC-efficient codes for quantum molecular dynamics via the Car–Parrinello scheme, though simulations are typically restricted to systems of up to a few hundred atoms (extendable via QM/MM coupling).

---

## 2. Theoretical and Methodological Foundations

### 2.1 Car–Parrinello Molecular Dynamics
The core idea (Car & Parrinello, 1985) treats the Kohn–Sham orbital coefficients as fictitious dynamical degrees of freedom with an associated fictitious electron mass, coupled to a classical Lagrangian for the nuclei:

$$
\mathcal{L}_{CP} = \sum_I \frac{1}{2} M_I \dot{R}_I^2 + \sum_i \mu \int |\dot{\psi}_i(\mathbf{r})|^2\, d\mathbf{r} - E_{KS}[\{\psi_i\}, \{R_I\}] + \text{constraints}
$$

This avoids re-diagonalizing the Kohn–Sham Hamiltonian at every step: the electrons stay close to the Born–Oppenheimer surface via fictitious dynamics rather than explicit minimization, provided the fictitious mass $\mu$ and timestep are chosen to maintain adiabatic separation between ionic and electronic motion.

### 2.2 Born–Oppenheimer Molecular Dynamics
CPMD also implements direct BOMD, where the electronic ground state is optimized (via direct minimization or diagonalization schemes) at every nuclear step, which is more robust for systems where the CP adiabaticity condition is difficult to satisfy (e.g., small HOMO–LUMO gaps, metals) but computationally more expensive per step.

### 2.3 Electronic Structure Machinery
- **Basis set:** plane waves with periodic boundary conditions (natural for condensed-phase/crystalline systems); isolated/cluster systems handled via appropriate boundary condition treatments (e.g., Poisson solvers/cutoff Coulomb techniques).
- **Pseudopotentials:** norm-conserving pseudopotentials (including Goedecker–Teter–Hutter, GTH, separable dual-space pseudopotentials) and ultrasoft pseudopotentials.
- **Exchange-correlation:** LDA, LSD, and the most widely used GGA functionals, plus hybrid functionals in later versions.
- **k-point sampling** for periodic/crystalline systems, and exploitation of molecular/crystal symmetry.
- **Wavefunction optimization:** direct minimization and iterative diagonalization routes to obtain Kohn–Sham energies/orbitals.

### 2.4 Structural and Dynamical Capabilities
- Geometry optimization (local optimization, simulated annealing).
- Molecular dynamics ensembles: constant energy (NVE), constant temperature (Nosé–Hoover thermostats), constant pressure (Parrinello–Rahman-type variable-cell Lagrangian).
- Path integral molecular dynamics (nuclear quantum effects), including efficient two-level parallelization schemes.
- Free-energy functional methods for metallic/finite-temperature electronic systems (Alavi/Deutsch-type free energy functional with k-points).
- Metadynamics (including coarse-grained, non-Markovian variants) for free-energy surface exploration and rare-event sampling.
- Response/perturbation-theory functionality (e.g., vibrational spectra, polarizabilities) and calculation of a range of electronic/spectroscopic properties.
- Time-dependent DFT (TDDFT) for excited states, and excited-state (e.g., Born–Oppenheimer surface-hopping-type) molecular dynamics.
- External time-dependent/static field application, including specialized periodic-boundary-condition-compatible field implementations (e.g., sinusoidal external potentials) contributed by the user community.

### 2.5 Multiscale / QM-MM Extensions
- A hybrid **QM/MM interface** was historically built around routines from the GROMOS96 MD code, enabling extension to larger biomolecular systems beyond the few-hundred-atom regime feasible for pure QM treatment.
- **MiMiC** (Multiscale Modeling in Computational Chemistry) is the modern, actively maintained QM/MM interface coupling CPMD (QM region) with **GROMACS** (MM region), implementing electrostatic embedding with a scheme developed by the group of Prof. Ursula Röthlisberger. MiMiC treats bonded/non-bonded cross terms at the QM/MM boundary and requires CPMD version 4.1 or later built with MiMiC support; setup involves a Python preprocessing script (`prepare-qmmm.py`) that maps a GROMACS index/topology to the CPMD QM/MM input.
- BioExcel (a Centre of Excellence for Computational Biomolecular Research) supported development of this modernized QM/MM coupling to improve on the performance limitations of the original GROMOS96-based interface.

---

## 3. Software Architecture and Distribution

| Aspect | Detail |
|---|---|
| Primary language | Fortran (with MPI parallelization; hybrid MPI/OpenMP in newer builds) |
| License (current) | MIT License (since the 2022 GitHub release) |
| License (historical) | Free registration-based license for non-profit/non-commercial users; separate commercial licensing route via the CPMD Consortium |
| Copyright holders | IBM Corp. (1990–2023); Max Planck Institute, Stuttgart (1994–2001) |
| Repository host | GitHub organization `CPMD-code` (main code, `Tests`, `Regtests`, `Additional_Packages`, `CPQA` quality-assurance tooling) |
| Community fork/initiative | **OpenCPMD** — community-driven continuation for documentation, tutorials, and collaborative development post-2022 |
| Typical scale | Systems up to a few hundred atoms at the pure-QM (plane-wave DFT) level; larger systems addressed via QM/MM (GROMOS96 interface or MiMiC/GROMACS) |
| Pseudopotential ecosystem | GTH (Goedecker–Teter–Hutter) dual-space pseudopotentials widely used and distributed (e.g., via the OpenCPMD `GTH-pseudopotentials` repository) |
| Typical execution | `cpmd.x file.in [PP_path] > file.out`, where `file.in` is the CPMD input file and pseudopotential files are supplied via a specified path |
| Supported platforms | Historically ported across many HPC platforms; particular historical optimization for Cray systems (shmem interface) and other major HPC architectures |

### 3.1 Brief Version History Highlights
- **Early 1990s:** Original implementation at IBM Zurich; first major update (1993) added keyword-driven input, atomic pseudo-wavefunction initial guesses, geometry optimization, additional MD types, Nosé–Hoover thermostats, and a diagonalization routine (code size ≈17,000 lines).
- **1994:** Improved communication layer and an MPI library interface; version 2.5 reached a stable milestone; a separate ab initio path-integral branch developed by Dominik Marx.
- **Version 3.0:** Improved portability (including a Cray shmem interface), constant-pressure MD via a Parrinello–Rahman-type Lagrangian, symmetry constraints, and Stefan Goedecker's dual-space pseudopotentials; code grew to ≈55,000 lines.
- **Version 3.1:** Minor update, but the starting point for the free-energy functional (with k-points) development (Alavi and Deutsch) and an efficient two-level-parallel path-integral implementation.
- **Version 4.x series:** Continued expansion (response functions, excited states, TDDFT, metadynamics, improved QM/MM); version 4.1 is the commonly cited baseline for modern QM/MM (MiMiC) support.
- **2022:** Public release on GitHub under the MIT License, opening the codebase to unrestricted community access and contribution.

---

## 4. Typical Application Domains

- Liquids and solvation phenomena (notably liquid water and aqueous solutions), where CPMD has been extensively used historically.
- Surfaces and heterogeneous catalysis.
- Crystalline solids and semiconductors.
- Biomolecular systems via QM/MM (enzymatic reaction mechanisms, active-site chemistry, metalloproteins) where universal classical force fields are insufficient because bond-breaking/forming or electronic effects matter.
- Reaction mechanism studies and free-energy sampling of chemical reactions (via metadynamics and constrained/thermodynamic-integration approaches).
- Response-property and spectroscopic predictions (vibrational and electronic spectra).
- Interfacial and electrochemical phenomena, including studies incorporating external static/time-dependent electric fields.

---

## 5. Related and Complementary Software

- **CP2K** — an alternative open-source (GPL) plane-wave/Gaussian mixed-basis DFT and Car–Parrinello MD code, often compared to CPMD; historically CPMD's Quality Assurance (CPQA) tooling has been shared/adapted with the CP2K distribution.
- **GROMACS** — coupled to CPMD via the MiMiC interface for QM/MM simulations.
- **GROMOS96** — source of the original QM/MM coupling routines (separately licensed).
- **pwtools** — a third-party Python pre-/post-processing toolkit supporting CPMD (alongside Quantum ESPRESSO, CP2K, and LAMMPS) output/input handling.

---

## 6. Key Publications on CPMD Theory and Methodology

### 6.1 Foundational Method Papers
- R. Car and M. Parrinello, *"Unified Approach for Molecular Dynamics and Density-Functional Theory,"* **Physical Review Letters**, 55, 2471 (1985). — The original Car–Parrinello paper introducing the extended-Lagrangian *ab initio* MD scheme underlying CPMD.
- P. Hohenberg and W. Kohn, *"Inhomogeneous Electron Gas,"* **Physical Review**, 136, B864 (1964). — Foundational DFT theorem paper.
- W. Kohn and L. J. Sham, *"Self-Consistent Equations Including Exchange and Correlation Effects,"* **Physical Review**, 140, A1133 (1965). — Kohn–Sham DFT formalism underlying all CPMD electronic-structure calculations.

### 6.2 The CPMD Code Reference Papers
- J. Hutter and M. Iannuzzi, *"CPMD: Car–Parrinello molecular dynamics,"* **Zeitschrift für Kristallographie – Crystalline Materials**, 220(5–6), 549–551 (2005). — The standard citable reference paper describing the CPMD code, its plane-wave/pseudopotential Kohn–Sham DFT foundation, feature set, and applicability to liquids, surfaces, crystals, and biomolecules.
- CPMD, http://www.cpmd.org/ and https://github.com/CPMD-code — Copyright IBM Corp. 1990–2023, Copyright MPI für Festkörperforschung Stuttgart 1997–2001 (the code's own standard self-citation used across the literature).

### 6.3 Pseudopotentials
- S. Goedecker, M. Teter, and J. Hutter, *"Separable dual-space Gaussian pseudopotentials,"* **Physical Review B**, 54, 1703 (1996). — The GTH pseudopotential scheme central to CPMD's standard pseudopotential library.
- C. Hartwigsen, S. Goedecker, and J. Hutter, *"Relativistic separable dual-space Gaussian pseudopotentials from H to Rn,"* **Physical Review B**, 58, 3641 (1998). — Extension of the GTH pseudopotentials across the periodic table.

### 6.4 Second-Generation Car–Parrinello / Methodological Extensions
- J. VandeVondele, M. Krack, F. Mohamed, M. Parrinello, T. Chassaing, and J. Hutter, *"Quickstep: Fast and accurate density functional calculations using a mixed Gaussian and plane waves approach,"* **Computer Physics Communications**, 167, 103 (2005). — Related mixed-basis approach informing the wider CP-based code ecosystem (implemented in CP2K, methodologically connected to CPMD's lineage).
- P. Tangney, *"On the theory underlying the Car–Parrinello method and the role of the fictitious mass parameter,"* **Journal of Chemical Physics**, 124, 044111 (2006). — Analysis of the fictitious-mass/adiabaticity foundations of the CP method as implemented in codes like CPMD.
- A. M. N. Niklasson, C. J. Tymczak, and M. Challacombe, *"Time-Reversible Born–Oppenheimer Molecular Dynamics,"* **Physical Review Letters**, 97, 123001 (2006). — Development relevant to the BOMD side of CPMD-type codes, improving energy conservation in direct BOMD schemes.

### 6.5 Free Energy Functional / Metallic Systems
- S. de Gironcoli, *related finite-temperature DFT works*; and A. Alavi, J. Kohanoff, M. Parrinello, D. Frenkel, *"Ab Initio Molecular Dynamics with Excited Electrons,"* **Physical Review Letters**, 73, 2599 (1994). — Free-energy functional approach for finite-electronic-temperature/metallic ab initio MD as later incorporated into CPMD.

### 6.6 Nosé–Hoover Thermostats and Extended-System Methods
- S. Nosé, *"A unified formulation of the constant temperature molecular dynamics methods,"* **Journal of Chemical Physics**, 81, 511 (1984).
- W. G. Hoover, *"Canonical dynamics: Equilibrium phase-space distributions,"* **Physical Review A**, 31, 1695 (1985). — Together, the Nosé–Hoover thermostat formalism used for constant-temperature CPMD simulations.
- G. J. Martyna, M. L. Klein, and M. Tuckerman, *"Nosé–Hoover chains: The canonical ensemble via continuous dynamics,"* **Journal of Chemical Physics**, 97, 2635 (1992). — Nosé–Hoover chain thermostatting as implemented for robust temperature control in CPMD.

### 6.7 Constant-Pressure / Variable-Cell Dynamics
- M. Parrinello and A. Rahman, *"Polymorphic transitions in single crystals: A new molecular dynamics method,"* **Journal of Applied Physics**, 52, 7182 (1981). — The Parrinello–Rahman variable-cell Lagrangian underlying CPMD's constant-pressure MD.
- M. Parrinello and A. Rahman, *"Crystal Structure and Pair Potentials: A Molecular-Dynamics Study,"* **Physical Review Letters**, 45, 1196 (1980).

### 6.8 Path Integral Ab Initio Molecular Dynamics
- D. Marx and M. Parrinello, *"Ab initio path integral molecular dynamics: Basic ideas,"* **Journal of Chemical Physics**, 104, 4077 (1996). — Path-integral CPMD methodology for nuclear quantum effects.
- D. Marx, M. E. Tuckerman, and G. J. Martyna, *"Quantum dynamics via adiabatic ab initio centroid molecular dynamics,"* **Computer Physics Communications**, 118, 166 (1999).

### 6.9 QM/MM and Multiscale Modeling
- A. Laio, J. VandeVondele, and U. Röthlisberger, *"A Hamiltonian electrostatic coupling scheme for hybrid Car–Parrinello molecular dynamics simulations,"* **Journal of Chemical Physics**, 116, 6941 (2002). — The Roethlisberger-group electrostatic embedding scheme underlying CPMD's (and later MiMiC's) QM/MM coupling.
- V. Kabadshow, ... (general MiMiC-related literature) — describing the MiMiC coupling of CPMD with GROMACS via a shared multiple-program-multiple-data (MPMD) MPI framework and electrostatic embedding for cross terms at the QM/MM boundary (see GROMACS/MiMiC reference-manual documentation for coupling formalism citation numbering).

### 6.10 Metadynamics
- A. Laio and M. Parrinello, *"Escaping free-energy minima,"* **Proceedings of the National Academy of Sciences**, 99, 12562 (2002). — Foundational metadynamics method, later incorporated (including coarse-grained, non-Markovian variants) into CPMD for free-energy surface exploration.

### 6.11 Excited States / TDDFT in CPMD
- I. Tavernelli, U. F. Röthlisberger, and collaborators — CPMD-implemented linear-response and real-time TDDFT methodology papers underlying excited-state and non-adiabatic/excited-state molecular dynamics functionality in CPMD (see CPMD manual bibliography for the full, version-specific citation set).

### 6.12 Historical/Contextual Review
- D. Marx and J. Hutter, *Ab Initio Molecular Dynamics: Basic Theory and Advanced Methods*, Cambridge University Press (2009). — The standard textbook-level treatment of the Car–Parrinello and Born–Oppenheimer ab initio MD methodology as embodied in CPMD, written by two of the method's/code's principal developers.

> **Note on completeness:** CPMD's official manual (distributed with the source, e.g., the CPMD 4.3 manual) contains an extensive, version-specific bibliography (typically 150+ entries) covering every implemented method (pseudopotentials, functionals, response theory, TDDFT, path integrals, QM/MM, metadynamics, etc.) with exact equation-level attributions. The list above collects the primary, most load-bearing theoretical references; for exhaustive, section-by-section citations, the CPMD manual's bibliography (packaged with each release on GitHub) is the authoritative source.

---

## 7. Access and Further Resources

- Source code: GitHub organization **CPMD-code** (`CPMD`, `Tests`, `Regtests`, `Additional_Packages`, `CPQA`).
- Community initiative: **OpenCPMD** (GitHub organization `OpenCPMD`) — tutorials, GTH pseudopotential collections, regression tests.
- QM/MM: **MiMiC** documentation is maintained as part of the GROMACS reference manual (special topics section).
- Historical/manual documentation: CPMD manuals (e.g., version 3.17.1 and 4.3) contain full keyword references, theoretical background chapters, and bibliographies.

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the CPMD (Car-Parrinello Molecular Dynamics) 	Plane-wave code specialized in Car-Parrinello and Born-Oppenheimer ab initio molecular dynamics. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
