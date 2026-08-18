# Octopus: A Real-Space, Real-Time DFT/TDDFT Code for Finite and Periodic Systems

## 1. Overview

Octopus is a free, open-source (GPL-licensed) electronic-structure package built around **real-space grid discretization** of the Kohn–Sham (KS) equations, with **real-time propagation** of the time-dependent Kohn–Sham (TDKS) equations as its defining methodological signature. It was conceived specifically to compute excited-state and spectroscopic properties of matter — optical absorption, nonlinear optical response, and light-driven electron/ion dynamics — for both finite systems (atoms, molecules, clusters, nanostructures) and periodic systems (crystals, surfaces, 1D/2D/3D solids).

The project originated around 2000 in Ángel Rubio's group (then at the University of Valladolid, Spain), with the first Octopus-based publication appearing in 2001. Its original purpose was to implement real-time TDDFT, at the time a newly proposed route to excited-state properties of molecules, as an alternative to the linear-response (Casida) formulation dominant in quantum chemistry. Over more than two decades the code has grown into a general-purpose electronic-structure and quantum-dynamics platform (on the order of hundreds of thousands of lines of Fortran/C, contributed to by dozens of developers across many countries) while retaining its founding emphasis on excited states driven far from equilibrium.

## 2. Core Methodology

### 2.1 Real-space grid representation
- Kohn–Sham orbitals, densities, and potentials are represented by their numerical values on a set of points in real space, rather than expanded in a basis of analytic functions (plane waves, Gaussians, atomic orbitals).
- Spatial derivatives (kinetic energy operator, Laplacians) are evaluated with high-order **finite-difference stencils**.
- Grids can be uniform or, in more recent developments, adaptive/non-uniform, with support for arbitrary user-defined simulation-box shapes (spheres, unions of spheres around atoms, cylinders, minimum-box shapes, parallelepipeds).
- The real-space approach gives systematic, controllable convergence with a single parameter (grid spacing), natural treatment of open (finite) boundary conditions, and straightforward extension to periodic boundary conditions in 1, 2, or 3 dimensions — including mixed periodicity (e.g., slabs, wires).
- Core electrons are eliminated via **norm-conserving pseudopotentials** (historically Troullier–Martins-type; later extended to Hartwigsen–Goedecker–Hutter and other separable forms), and the projected augmented-wave (PAW) formalism is also supported for improved accuracy near nuclei.

### 2.2 Ground-state DFT
Octopus solves the static Kohn–Sham equations self-consistently as a foundation for excited-state calculations, supporting:
- LDA, GGA, meta-GGA, and hybrid functionals (via LibXC, the exchange-correlation functional library co-developed by the Octopus team and now used across many electronic-structure codes).
- Exact-exchange / optimized-effective-potential (OEP) and Kohn–Sham exact-exchange (KLI) methods.
- Spin-polarized (collinear and non-collinear) and spin–orbit coupling calculations.
- Structural relaxation and geometry optimization.

### 2.3 Excited-state and spectroscopic methods
Octopus offers **three principal routes to excited-state properties within TDDFT**, plus additional post-processing/embedding schemes:

1. **Real-time TDDFT (RT-TDDFT).** The signature capability: KS orbitals are propagated explicitly in time under the time-dependent KS Hamiltonian, typically after an impulsive (delta-function) electric-field "kick" (for linear absorption spectra) or under an explicit laser pulse (for nonlinear/strong-field response). The optical absorption spectrum is obtained from the Fourier transform of the induced time-dependent dipole moment. This approach:
   - Scales favorably with system size relative to matrix-diagonalization methods, since it avoids computing/storing the Casida response matrix.
   - Naturally captures nonlinear and non-perturbative phenomena (high-harmonic generation, above-threshold phenomena, strong-field ionization) inaccessible to linear response.
   - Extends directly to periodic systems using time-dependent Bloch orbitals and gauge choices (length gauge with a time-dependent vector potential, or "kick" plus crystal-momentum evolution) suited to extended systems, enabling real-time optical spectra, high-harmonic generation, and other nonlinear responses in solids.
   - Supports coupling to classical light fields, and, through Octopus's quantum-electrodynamical DFT (QEDFT) extension, quantized photon modes for polaritonic/strong-coupling chemistry and cavity-modified spectroscopy.

2. **Casida linear-response TDDFT** (frequency-domain). Excitation energies and oscillator strengths are obtained by diagonalizing the Casida matrix eigenvalue equation built from occupied–unoccupied Kohn–Sham orbital pairs and the exchange-correlation (and Hartree) kernel. This is the standard route for a small number of well-resolved discrete excitations in finite systems (as in quantum-chemistry codes), and Octopus supports:
   - Singlet and triplet excitations for closed-shell references.
   - Excited-state (analytic) forces/gradients (experimental), enabling excited-state geometry optimization and vibrational analysis on excited potential-energy surfaces.
   - Selection of active occupied/unoccupied KS-state windows to truncate the particle–hole basis.

3. **Sternheimer linear-response TDDFT** (density-functional perturbation theory, DFPT). Rather than diagonalizing the full Casida matrix or summing over unoccupied states, the frequency-dependent Sternheimer equation is solved directly (a linear-response Kohn–Sham equation at fixed frequency), avoiding explicit reference to unoccupied-state sums. This scales more favorably for larger systems and dense/continuous spectra, and gives access to:
   - Static and dynamical polarizabilities and hyperpolarizabilities.
   - Van der Waals (dispersion) coefficients.
   - Vibrational (phonon) properties and infrared intensities.
   - Resonant and non-resonant response beyond the lowest few discrete excitations, complementing Casida in the regime where dense excitation spectra or continua make matrix diagonalization impractical.

4. **Bethe–Salpeter-equation (BSE) interoperability.** Octopus can export ground-state and screening data required to perform GW quasiparticle and Bethe–Salpeter calculations with external many-body-perturbation-theory codes (e.g., the BerkeleyGW code), extending Octopus-based excited-state workflows beyond (TD)DFT into explicit many-body treatments of excitons in periodic systems.

### 2.4 Dynamics beyond the electronic problem
- **Ehrenfest and Born–Oppenheimer molecular dynamics**, allowing coupled electron–ion (nonadiabatic) real-time dynamics — important for excited-state relaxation, photochemistry, and vibronic coupling.
- A modified/multiple-time-step Ehrenfest scheme for adiabatic dynamics of larger systems with favorable scaling.
- **Quantum optimal control theory (QOCT)**, used to design laser pulses that steer the system toward a target excited-state population or observable — directly tied to spectroscopic and photochemical objectives.
- **Real-time and linear-response treatments of core electrons/X-ray spectroscopies**, including simulation of core-level excitations relevant to X-ray absorption.
- Time-dependent current-density functional theory (TDCDFT)-adjacent capabilities and current-density observables.

### 2.5 Coupled light–matter / polaritonic and photoemission features
- Quantum-electrodynamical DFT (QEDFT) coupling electrons to quantized photon modes, used for cavity-QED and polaritonic-chemistry simulations, including modified exchange-correlation functionals derived for strongly coupled light–matter systems.
- Simulation of photoemission spectra and photoelectron distributions via time-dependent wavefunction propagation into the continuum (mask/absorbing-boundary techniques).
- Plasmonic-system simulations (metallic nanoparticle optical response).
- Electronic circular dichroism (ECD) spectra from real-time propagation.
- Electronic-stopping-power calculations (energy loss of fast charged projectiles in electronic systems), a distinctive Octopus strength connecting excited-state electron dynamics to radiation-damage/ion-beam physics.
- High-harmonic generation (HHG) in solids and molecules, including spin-polarized and topological materials.

### 2.6 Observables, output, and interoperability
Octopus can compute and export densities, orbitals, current densities, and the time-dependent electron localization function (TD-ELF), among other quantities, for visualization and post-processing, and interfaces with external packages (e.g., for GW/BSE follow-up calculations) to broaden the excited-state toolchain beyond native (TD)DFT.

## 3. Numerical and Computational Infrastructure

- **Discretization:** finite-difference real-space grids (uniform or adaptive), arbitrary simulation-cell shapes, and support for arbitrary spatial dimensionality (used e.g. for reduced-dimensionality model systems as well as standard 1D/2D/3D physical systems).
- **Boundary conditions:** open (zero) boundary conditions for finite systems and periodic boundary conditions (1D/2D/3D) for extended systems, plus absorbing boundaries/complex absorbing potentials for photoemission and ionization simulations.
- **Time propagators:** a range of unitary/norm-conserving propagators for the TDKS equations (e.g., approximate enforced time-reversal symmetry (AETRS), Magnus expansion-based, exponential midpoint, and Crank–Nicolson-type schemes), chosen to balance accuracy, stability, and cost, particularly under strong time-dependent fields.
- **Parallelization:** large-scale hybrid MPI/OpenMP (and GPU-accelerated) parallelization across real-space domains, k-points, electronic states, and independent time-dependent "runs" (e.g., multiple laser-kick directions/frequencies), enabling use on massively parallel HPC architectures — a focus of dedicated methodological papers (see §5).
- **Exchange-correlation functionals:** delivered through LibXC, giving access to the large standard library of LDA/GGA/meta-GGA/hybrid functionals used across the broader electronic-structure community.

## 4. Scope of Systems and Typical Applications

Octopus targets systems ranging from small molecules to large organic compounds, metallic and semiconductor clusters of hundreds of atoms, low-dimensional materials (2D materials, nanowires, nanotubes), surfaces, and bulk periodic solids. Representative application areas driven by its excited-state/spectroscopic machinery include:
- Linear optical absorption spectra of molecules, clusters, and solids.
- Nonlinear and strong-field phenomena: high-harmonic generation, above-threshold ionization, multiphoton processes.
- Photoemission and time-resolved photoelectron spectroscopy.
- Electronic circular dichroism.
- Plasmonics and nanophotonics.
- Polaritonic chemistry / strong light–matter coupling (cavity QED).
- Nonadiabatic (Ehrenfest) excited-state molecular dynamics and photochemistry.
- Electronic stopping power for ion-beam and radiation-damage studies.
- Spin dynamics and magnon physics, spin-polarized HHG.
- Topological and Floquet-engineered electronic states under periodic driving.
- Quantum optimal control of excited-state populations and observables.

## 5. Development History and Milestones

| Period | Milestone |
|---|---|
| ~2000 | Project initiated in Á. Rubio's group (Univ. of Valladolid, Spain), aimed at real-time TDDFT for excited states. |
| 2001 | First scientific publication using Octopus. |
| 2003 | First dedicated code paper, "Octopus: a first-principles tool for excited electron-ion dynamics" (Comput. Phys. Commun.), establishing RT-TDDFT + Ehrenfest dynamics for finite systems. |
| 2006 | Major methodology update paper (Phys. Status Solidi B), introducing adaptive coordinates, extension of the real-space technique to periodic systems, and large-scale parallelization. |
| 2012 | "Time-dependent density-functional theory in massively parallel computer architectures: the Octopus project" (J. Phys.: Condens. Matter), documenting HPC-scale parallelization strategies. |
| 2015 | "Real-space grids and the Octopus code as tools for the development of new simulation approaches for electronic systems" (Phys. Chem. Chem. Phys.), broadening scope to plasmonics, quantum optimal control, and new real-space methodology. |
| 2020 | "Octopus, a computational framework for exploring light-driven phenomena and quantum dynamics in extended and finite systems" (J. Chem. Phys., JCP Special Topic on Electronic Structure Software) — the current comprehensive reference paper, covering QEDFT, periodic RT-TDDFT, magnons, and much-expanded functionality. |
| Ongoing | Continued development (GPU acceleration, PAW, machine-learning-assisted functionals, extended polaritonic/cavity-QED features, TDHF via adaptively compressed exchange, core-level/X-ray spectroscopy). |

## 6. Licensing and Availability

Octopus is distributed free of charge under the **GNU General Public License (GPL)**. Source code, documentation, tutorials, and a database of publications using the code are maintained at the project website (octopus-code.org; historically also hosted at tddft.org/programs/octopus).

## 7. Summary Assessment

Octopus occupies a distinctive niche among first-principles excited-state codes: rather than treating linear-response (Casida-type) TDDFT as the primary tool, it is architected from the ground up around **real-time propagation on real-space grids**, giving it particular strength in (i) nonlinear and non-perturbative light-driven phenomena, (ii) systems and processes where a dense or continuous excitation spectrum makes matrix-diagonalization approaches costly, (iii) coupled electron-ion nonadiabatic dynamics, and (iv) periodic as well as finite systems within a single unified framework. Its three complementary excited-state engines (real-time, Casida, Sternheimer/DFPT) let users trade off between these regimes within one code, and its GW/BSE export pathway and QEDFT extensions connect it to explicit many-body and cavity-QED treatments of excitations. This combination — real-space flexibility, native periodicity, strong-field/nonlinear capability, and open-source HPC scalability — explains its wide adoption for spectroscopic simulation across molecular, nano, and solid-state physics/chemistry communities.

## 8. Key Publications on Octopus Theory and Methodology

1. Marques, M. A. L.; Castro, A.; Bertsch, G. F.; Rubio, A. *Octopus: a first-principles tool for excited electron-ion dynamics.* Comput. Phys. Commun. **151**, 60–78 (2003).
2. Castro, A.; Appel, H.; Oliveira, M.; Rozzi, C. A.; Andrade, X.; Lorenzen, F.; Marques, M. A. L.; Gross, E. K. U.; Rubio, A. *Octopus: a tool for the application of time-dependent density functional theory.* Phys. Status Solidi B **243**(11), 2465–2488 (2006).
3. Andrade, X.; Alberdi-Rodriguez, J.; Strubbe, D. A.; Oliveira, M. J. T.; Nogueira, F.; Castro, A.; Muguerza, J.; Arruabarrena, A.; Louie, S. G.; Aspuru-Guzik, A.; Rubio, A.; Marques, M. A. L. *Time-dependent density-functional theory in massively parallel computer architectures: the Octopus project.* J. Phys.: Condens. Matter **24**, 233202 (2012).
4. Andrade, X.; Strubbe, D. A.; De Giovannini, U.; Larsen, A. H.; Oliveira, M. J. T.; Alberdi-Rodriguez, J.; Varas, A.; Theophilou, I.; Helbig, N.; Verstraete, M. J.; Stella, L.; Nogueira, F.; Aspuru-Guzik, A.; Castro, A.; Marques, M. A. L.; Rubio, A. *Real-space grids and the Octopus code as tools for the development of new simulation approaches for electronic systems.* Phys. Chem. Chem. Phys. **17**, 31371–31396 (2015).
5. Tancogne-Dejean, N.; Oliveira, M. J. T.; Andrade, X.; Appel, H.; Borca, C. H.; Le Breton, G.; Buchholz, F.; Castro, A.; Corni, S.; Correa, A. A.; De Giovannini, U.; Delgado, A.; Eich, F. G.; Flick, J.; Gil, G.; Gomez, A.; Helbig, N.; Hübener, H.; Jestädt, R.; Jornet-Somoza, J.; Larsen, A. H.; Lebedeva, I. V.; Lüders, M.; Marques, M. A. L.; Ohlmann, S. T.; Pipolo, S.; Rampp, M.; Rozzi, C. A.; Strubbe, D. A.; Sato, S. A.; Schäfer, C.; Theophilou, I.; Welden, A.; Rubio, A. *Octopus, a computational framework for exploring light-driven phenomena and quantum dynamics in extended and finite systems.* J. Chem. Phys. **152**, 124119 (2020).

These five papers form the core methodological lineage of the code (2003 → 2006 → 2012 → 2015 → 2020), each documenting a successive expansion of theoretical scope (finite-system RT-TDDFT and Ehrenfest dynamics → periodic systems and adaptive coordinates → HPC parallelization → new real-space simulation approaches, plasmonics, and optimal control → QEDFT, magnons, and the fully generalized light-driven-phenomena framework).

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the Octopus 	Real-space, real-time DFT/TDDFT code for finite and periodic systems, focused on excited-state and spectroscopic properties. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
