# ABINIT — Exhaustive Technical Review

**Category:** Open-source plane-wave / PAW electronic-structure package
**Core strengths:** Density-Functional Perturbation Theory (DFPT), phonons and lattice dynamics, many-body perturbation theory (GW, Bethe–Salpeter Equation)
**License:** GNU General Public License (GPL)
**Website / repository:** [abinit.org](https://www.abinit.org) · [github.com/abinit/abinit](https://github.com/abinit/abinit) · [docs.abinit.org](https://docs.abinit.org)

---

## 1. Overview

ABINIT is one of the earliest electronic-structure packages to have been released under an open-source license (its first public distribution dates to 1998, with the project's academic roots going back to the early 1990s at UCL/Louvain-la-Neuve and Corning Inc.). It implements Kohn–Sham Density Functional Theory (DFT), Density-Functional Perturbation Theory (DFPT), Many-Body Perturbation Theory (MBPT: GW approximation and the Bethe–Salpeter equation), and a growing set of more specialized formalisms such as Dynamical Mean-Field Theory (DMFT) and the temperature-dependent effective potential (TDEP) approach to anharmonic lattice dynamics.

The code relies on a **plane-wave basis set** for wavefunctions, densities, and other space-dependent quantities, combined with **norm-conserving pseudopotentials** or the **Projector-Augmented Wave (PAW)** method to treat the core–valence separation. This makes ABINIT naturally suited to periodic crystalline solids, while molecules, surfaces, and nanostructures can also be treated via the supercell technique.

ABINIT is developed by a large, distributed international collaboration (originally centered at the Université catholique de Louvain, with major contributing groups at CEA, Liège, ETH Zürich, MIT, and many other institutions worldwide) and is maintained under a rolling development model with a public Git/GitLab repository, an extensive automated test suite, and structured documentation (input-variable database, "topics" pages, and tutorials).

---

## 2. Theoretical & Methodological Capabilities

### 2.1 Ground-state DFT
- Kohn–Sham DFT with LDA, GGA (PBE, PBEsol, etc.), meta-GGA, hybrid functionals (PBE0, HSE), and DFT+U / DFT+DMFT for strongly correlated systems.
- van der Waals corrections (DFT-D/Grimme schemes, vdW-DF functionals).
- Constrained DFT (CDFT) via Lagrange-multiplier methods for pinning charges/magnetic moments and probing excited/local states.
- Non-collinear magnetism and spin-orbit coupling.
- Finite electric fields and Berry-phase polarization (modern theory of polarization).
- Structural relaxation, molecular dynamics (including path-integral MD, extended first-principles MD for warm dense matter/plasma regimes).
- Real-time TDDFT (RT-TDDFT) for non-adiabatic dynamics and laser-pulse response (actively expanding).

### 2.2 Density-Functional Perturbation Theory (DFPT) — a flagship capability
ABINIT's DFPT implementation is one of the most complete of any open-source plane-wave/PAW code, allowing linear- and higher-order response to be computed self-consistently without the need for large supercells or finite differences. Response properties available include:

- **Phonons** at arbitrary **q**-points, dynamical matrices, interatomic force constants, and Fourier interpolation of phonon band structures/DOS (via the companion `anaddb` post-processor).
- **Dielectric tensors**, **Born effective charges**, and **LO–TO splitting**.
- **Elastic and piezoelectric tensors**, and **internal-strain (force-response) tensors** — response to homogeneous strain.
- **Raman scattering efficiencies** and the **electro-optic effect** (static non-linear response).
- **Non-linear optical susceptibilities** (second-harmonic generation and related third-order responses).
- **Flexoelectricity** and **dynamical (Born) quadrupoles** — spatial-dispersion (long-wave) DFPT.
- **Electron–phonon coupling**, including mobilities, temperature dependence of the electronic structure (zero-point renormalization of band gaps), and superconducting properties (Eliashberg-type quantities).
- **Effective masses** computed directly via DFPT (rather than finite differences of the band structure).
- Response to a **homogeneous magnetic field** and general spin-perturbation DFPT.
- DFPT combined with **PAW** (not just norm-conserving pseudopotentials), including cross terms in oscillator strengths (`pawcross`) and careful handling of the acoustic sum rule.
- DFPT is compatible with **van der Waals interactions** and **non-collinear magnetism**, which is unusually general among plane-wave DFPT implementations.
- **Anharmonic lattice dynamics** via the `a-TDEP` post-processor (Temperature-Dependent Effective Potential), giving temperature-renormalized phonon spectra, free energies, and thermodynamic properties beyond the harmonic approximation.
- Lattice and spin models fitted from DFPT/DFT data are handled by the **Multibinit** second-principles engine (effective Hamiltonians for ferroelectric/multiferroic lattice dynamics, spin dynamics, and lattice Wannier functions).

### 2.3 Many-Body Perturbation Theory: GW
- **G₀W₀** and iterative/self-consistent GW variants (including QSGW-like self-consistency on top of the self-energy) for quasiparticle band structures and band gaps.
- Compatible with **plane-wave/pseudopotential** and **PAW** formalisms (PAW-GW requires special care for the frozen-core contribution to the self-energy).
- Conventional GW implementation formally scaling as O(N⁴_atoms); the plasmon-pole and contour-deformation/imaginary-frequency approaches are supported for the dynamically-screened Coulomb interaction W.
- **Completeness-relation** technique to accelerate convergence with respect to the number of empty bands (`gwcomp`, `gwencomp`).
- A newer **GWR** (real-space/imaginary-time, low-scaling) implementation, and an ongoing **GW Perturbation Theory (GWPT)** development for electron-phonon-renormalized quasiparticle spectra.
- Coupling to the **GreenX** library for efficient minimax time/frequency grids used in low-scaling RPA and GW calculations.
- **cRPA** (constrained RPA) for computing effective (screened) Hubbard U/J parameters used in DFT+U and DFT+DMFT.
- Massive parallelism (band/FFT/k-point/frequency distributed parallelism) for GW workloads on HPC systems.

### 2.4 Many-Body Perturbation Theory: Bethe–Salpeter Equation (BSE)
- Solves the BSE on top of a GW (or model-dielectric-function) quasiparticle spectrum to obtain **optical absorption spectra including excitonic effects** — essential for accurately describing insulators and semiconductors where electron–hole binding is significant.
- Supports both **direct diagonalization** and **iterative (Haydock-type) solvers** for the BSE Hamiltonian, allowing access to large systems where full diagonalization is prohibitive.
- Provides the macroscopic dielectric function along multiple **q**-directions (reciprocal-lattice vectors and Cartesian axes) for optical-limit response, plus finite-momentum-transfer response relevant to EELS/inelastic X-ray scattering.
- PAW-compatible, enabling BSE spectra for systems requiring the augmentation-region description (transition metals, first-row elements).
- Has served as the electronic-structure engine behind external core-level BSE codes (e.g., the OCEAN package for XAS/XES/NRIXS combining ABINIT wavefunctions with the NIST BSE solver), illustrating the reach of ABINIT's MBPT outputs beyond the native code.

### 2.5 Other advanced/specific capabilities
- **Dynamical Mean-Field Theory (DMFT)**, including a self-consistent DFT+DMFT scheme within PAW, with support for external CT-HYB solvers via TRIQS.
- **DFT+U**, LDA−1/2-type schemes.
- **Positron annihilation** (two-component DFT for electron-positron states).
- **Nuclear-site properties**: electric field gradients (EFG), Mössbauer isomer shifts, orbital magnetization, NMR-relevant quantities.
- **Wannier functions** via interface to Wannier90, and native lattice Wannier functions in Multibinit.
- **Topological/geometric-phase properties**: Berry phases, orbital magnetization, interface to Z2Pack for topological invariants.
- **Electron-phonon transport** (mobilities, resistivities) through the `eph` driver and the `gstore` electron-phonon matrix-element storage framework.
- **High-pressure / warm dense matter** physics: extended first-principles molecular dynamics for very high temperature/pressure regimes relevant to planetary and stellar interiors.
- **GPU acceleration**: ground-state SCF loops, FFTs, and dense linear algebra offloaded via OpenMP and vendor libraries (cuFFT/cuBLAS on NVIDIA, rocFFT/rocBLAS on AMD), giving roughly an order-of-magnitude speed-up over CPU-only runs on systems of ~1000 atoms.
- **High-throughput infrastructure**: community-curated databases of phonon band structures, second-harmonic-generation coefficients, and GW band-gap benchmarks, plus the **AbiPy** Python post-processing/workflow ecosystem and integration with AiiDA-style workflow managers.
- **LIBPAW**: a standalone, reusable PAW library shared with other electronic-structure codes.

---

## 3. Software Architecture

- **Language/build:** Fortran (with C components), Autotools/CMake-style build system, extensive `configure` options for HPC environments; EasyBuild recipes available.
- **Parallelism:** MPI-based, with multiple simultaneous levels of parallel distribution — over k-points, spin/spinor components, plane-wave/FFT grids, electronic bands, perturbations (for DFPT), and (for GW/BSE) frequencies and transitions. GPU offloading (OpenMP target + vendor math libraries) supplements MPI parallelism for the costliest kernels.
- **Auxiliary programs bundled with the package**: `anaddb` (analysis of Derivative Data Bases — phonons, dielectric/elastic/piezoelectric properties, Fourier interpolation), `mrgddb` (merging DDB files), `cut3d` (post-processing densities/wavefunctions), `aim` (Bader/AIM analysis), `optic` (linear/non-linear optical spectra from band-structure data), `fold2Bloch` (band unfolding), `a-TDEP` (anharmonic lattice dynamics), `Multibinit` (second-principles lattice/spin dynamics), `abitk` (utility toolkit), `spacegroup` (symmetry utilities).
- **Pseudopotentials/PAW datasets:** ABINIT supports multiple pseudopotential formats (including Hartwigsen–Goedecker–Hutter, and importantly **ONCVPSP**-generated norm-conserving pseudopotentials, widely distributed via the **PseudoDojo** project) and PAW datasets (generated with tools such as **AtomPAW**), with community-curated tables optimized for accuracy/efficiency trade-offs.
- **Testing/verification:** a large regression-test suite (YAML-based test infrastructure) underpins each release, and the ABINIT team participates in cross-code DFT-precision verification efforts (the "Δ-factor"/reproducibility benchmarking exercises).
- **Companion Python ecosystem:** **AbiPy** for automating and post-processing ABINIT workflows (band structures, DFPT/phonon analysis, GW/BSE post-processing, high-throughput screening).

---

## 4. Typical Application Domains

- Vibrational/thermodynamic properties of solids (phonon dispersions, thermal expansion, heat capacity, anharmonic effects).
- Dielectric, piezoelectric, ferroelectric, and multiferroic materials (via DFPT strain/electric-field response and Multibinit second-principles modeling).
- Non-linear optical materials (Raman, SHG, electro-optic coefficients).
- Quasiparticle band structures and optical absorption spectra of semiconductors, insulators, and nanostructures (GW/BSE).
- Strongly correlated materials (actinides, lanthanides, transition-metal oxides) via DFT+U and DFT+DMFT.
- Electron-phonon-driven phenomena: temperature-dependent band gaps, carrier mobilities, superconductivity.
- Extreme conditions: high-pressure phase transitions, warm dense matter, planetary/stellar interior modeling.
- Defect physics (point defects, luminescence lineshapes via electron-phonon coupling, e.g., phosphor materials).

---

## 5. Key Publications on the Underlying Theory and Software

The list below follows the structure suggested by the **official ABINIT "Acknowledgments" documentation page** (docs.abinit.org/theory/acknowledgments), which distinguishes general project-overview articles ("Gen.") from papers describing specific theoretical/methodological capabilities ("Spe."). This is the authoritative citation guidance maintained by the ABINIT development team itself.

### 5.1 General project / software overview articles

1. **Verstraete, M. J. et al.** (2025). *Abinit 2025: New Capabilities for the Predictive Modeling of Solids and Nanomaterials.* — the most recent general overview article; the current default citation for any ABINIT-based work. (arXiv: 2507.08578)
2. **Gonze, X., Amadon, B., Antonius, G., et al.** (2020). *The Abinit project: Impact, environment and recent developments.* *Computer Physics Communications* **248**, 107042.
3. **Romero, A. H., Allan, D. C., Amadon, B., et al.** (2020). *ABINIT: Overview and focus on selected capabilities.* *Journal of Chemical Physics* **152**, 124102.
4. **Gonze, X., Jollet, F., Abreu Araujo, F., et al.** (2016). *Recent developments in the ABINIT software package.* *Computer Physics Communications* **205**, 106–131.
5. **Gonze, X., Amadon, B., Anglade, P.-M., et al.** (2009). *ABINIT: First-principles approach to material and nanosystem properties.* *Computer Physics Communications* **180**, 2582–2615.
6. **Gonze, X. et al.** (2005). *A brief introduction to the ABINIT software package.* *Zeitschrift für Kristallographie* **220**, 558–562.
7. **Gonze, X., Beuken, J.-M., Caracas, R., et al.** (2002). *First-principles computation of material properties: the ABINIT software project.* *Computational Materials Science* **25**, 478–492.

### 5.2 Methodology-specific references (PAW, DFPT, GW, BSE, and related capabilities)

**Projector-Augmented Wave (PAW) implementation**
- **Torrent, M., Jollet, F., Bottin, F., Zérah, G., Gonze, X.** (2008). *Implementation of the projector augmented-wave method in the ABINIT code: Application to the study of iron under pressure.* *Computational Materials Science* **42**, 337–351.
- **Blöchl, P. E.** (1994). *Projector augmented-wave method.* *Physical Review B* **50**, 17953–17979. (Foundational PAW formalism.)
- **Audouze, C., Jollet, F., Torrent, M., Gonze, X.** (2006). *Projector augmented-wave approach to density-functional perturbation theory.* *Physical Review B* **73**, 235101.
- **Audouze, C., Jollet, F., Torrent, M., Gonze, X.** (2008). *Comparison between projector augmented-wave and ultrasoft pseudopotential formalisms at the density-functional perturbation theory level.* *Physical Review B* **78**, 035105.

**DFPT: phonons, dielectric response, effective charges**
- **Gonze, X., Lee, C.** (1997). *Dynamical matrices, Born effective charges, dielectric permittivity tensors, and interatomic force constants from density-functional perturbation theory.* *Physical Review B* **55**, 10355–10368.
- **Gonze, X.** (1997). *First-principles responses of solids to atomic displacements and homogeneous electric fields: Implementation of a conjugate-gradient algorithm.* *Physical Review B* **55**, 10337–10354.
- **Baroni, S., de Gironcoli, S., Dal Corso, A., Giannozzi, P.** (2001). *Phonons and related crystal properties from density-functional perturbation theory.* *Reviews of Modern Physics* **73**, 515–562. (General DFPT review, widely cited alongside ABINIT's own DFPT papers.)
- **Baroni, S., Giannozzi, P., Testa, A.** (1987). *Green's-function approach to linear response in solids.* *Physical Review Letters* **58**, 1861–1864. (Foundational DFPT paper.)

**Strain / elastic and piezoelectric response**
- **Hamann, D. R., Wu, X., Rabe, K. M., Vanderbilt, D.** (2005). *Metric tensor formulation of strain in density-functional perturbation theory.* *Physical Review B* **71**, 035117.

**Non-linear optical response (Raman, electro-optic effect)**
- **Veithen, M., Gonze, X., Ghosez, P.** (2005). *Nonlinear optical susceptibilities, Raman efficiencies, and electro-optic tensors from first-principles density functional perturbation theory.* *Physical Review B* **71**, 125107.

**Thermal/phonon integration (thermodynamic properties)**
- **Lee, C., Gonze, X.** (1995). *Ab initio calculation of the thermodynamic properties and atomic temperature factors of SiO₂ α-quartz and stishovite.* *Physical Review B* **51**, 8610–8613.

**Flexoelectricity and dynamical (Born) quadrupoles**
- **Royo, M., Stengel, M.** (2019). *First-principles theory of spatial dispersion: Dynamical quadrupoles and flexoelectricity.* *Physical Review X* **9**, 021050.

**Effective masses via DFPT**
- **Laflamme Janssen, J., Rousseau, B., Côté, M.** (2016). *Efficient effective-mass calculations using density-functional perturbation theory.* *Physical Review B* **93**, 205147.

**GW approximation and self-consistency**
- **Bruneval, F., Vast, N., Reining, L.** (2006). *Effect of self-consistency on quasiparticles in solids.* *Physical Review B* **74**, 045102. (Cited when self-consistent GW/COHSEX/HF options are used.)
- **Bruneval, F., Gonze, X.** (2008). *Accurate GW self-energies in a plane-wave basis using only a few empty states: Towards large systems.* *Physical Review B* **78**, 085125. (Completeness-relation acceleration of GW convergence.)
- **Hybertsen, M. S., Louie, S. G.** (1986). *Electron correlation in semiconductors and insulators: Band gaps and quasiparticle energies.* *Physical Review B* **34**, 5390–5413. (Foundational GW methodology.)
- **Hedin, L.** (1965). *New method for calculating the one-particle Green's function with application to the electron-gas problem.* *Physical Review* **139**, A796–A823. (Original Hedin's-equations/GW formalism.)
- **Azizi, M., Wilhelm, J., Golze, D., et al.** (2024). *Validation of the GreenX library time-frequency component for efficient GW and RPA calculations.* *Physical Review B* **109**, 245101. (Low-scaling GW/RPA via minimax grids.)

**Bethe–Salpeter Equation and optical/excitonic spectra**
- **Albrecht, S., Reining, L., Del Sole, R., Onida, G.** (1998). *Ab initio calculation of excitonic effects in the optical spectra of semiconductors.* *Physical Review Letters* **80**, 4510–4513.
- **Benedict, L. X., Shirley, E. L., Bohn, R. B.** (1998). *Optical absorption of insulators and the electron-hole interaction: An ab initio calculation.* *Physical Review Letters* **80**, 4514–4517. (Foundational BSE-for-optics papers frequently cited alongside ABINIT's BSE module.)
- **Onida, G., Reining, L., Rubio, A.** (2002). *Electronic excitations: density-functional versus many-body Green's-function approaches.* *Reviews of Modern Physics* **74**, 601–659. (Standard GW/BSE review.)

**Massive parallelism**
- **Bottin, F., Leroux, S., Knyazev, A., Zérah, G.** (2008). *Large-scale ab initio calculations based on three levels of parallelization.* *Computational Materials Science* **42**, 329–336. (arXiv: 0707.3405)

**DFT+U**
- **Amadon, B., Lechermann, F., Georges, A., Jollet, F., Wehling, T. O., Lichtenstein, A. I.** (2008). *Plane-wave based electronic structure calculations for correlated materials using dynamical mean-field theory and projected local orbitals.* *Physical Review B* **77**, 205112.

**Norm-conserving (ONCVPSP) pseudopotentials**
- **Hamann, D. R.** (2013). *Optimized norm-conserving Vanderbilt pseudopotentials.* *Physical Review B* **88**, 085117.

**Van der Waals DFT-D (Grimme-type) functionals**
- **Van Troeye, B., Torrent, M., Gonze, X.** (2016). *Interatomic force constants including the DFT-D dispersion contribution.* *Physical Review B* **93**, 144304.

**Temperature dependence of electronic structure / zero-point renormalization**
- **Ponce, S., Antonius, G., Boulanger, P., Cannuccia, E., Gonze, X.** (2014). *Temperature dependence of electronic eigenenergies in the adiabatic harmonic approximation.* *Physical Review B* **90**, 214304.
- **Ponce, S., Antonius, G., Gillet, Y., Boulanger, P., Laflamme Janssen, J., Marini, A., Côté, M., Gonze, X.** (2015). *Temperature dependence of the electronic structure of semiconductors and insulators.* *Journal of Chemical Physics* **143**, 102813.

**Electron-phonon coupling driver (`optdriver = 7`)**
- **Brunin, G., Miglio, A., Ponce, S., et al.** (2020). *Electron-phonon beyond Fröhlich: Dynamical quadrupoles in polar and covalent solids.* *Physical Review Letters* **125**, 136601.
- **Brunin, G., Miglio, A., Ponce, S., et al.** (2020). *Phonon-limited electron mobility in Si, GaAs, and GaP with exact treatment of dynamical quadrupoles.* *Physical Review B* **102**, 094308.

**Linear-response U/J determination**
- **MacEnulty, L., et al.** (2024). *Reference paper for the Linear Response UJ (lrUJ) utility and renovated UJ Determination (UJdet) implementation.*
- **MacEnulty, L., et al.** (2023). *Related earlier paper on UJ determination methodology.*

*(Full, continuously updated BibTeX entries for all of the above — and for the several hundred additional method- and application-specific references cited throughout the ABINIT documentation — are maintained at [docs.abinit.org/theory/bibliography](https://docs.abinit.org/theory/bibliography/), with the master `.bib` file linked from that page.)*

---

## 6. Access to Tutorials and Documentation (for reference)

- Input-variable database: [docs.abinit.org/variables](https://docs.abinit.org/variables/)
- Topic-organized documentation (DFPT, GW, GWR, BSE, PAW, EPH, etc.): [docs.abinit.org/topics](https://docs.abinit.org/topics/)
- Tutorials, including dedicated DFPT (RF1/RF2), GW1/GW2, GWR, and Bethe-Salpeter tutorials: [docs.abinit.org/tutorial](https://docs.abinit.org/tutorial/)
- Theory pages (PAW, MBPT, BSE formalism write-ups): [docs.abinit.org/theory](https://docs.abinit.org/theory/mbt/)

---

*Document compiled from ABINIT's official documentation (docs.abinit.org), the ABINIT project's general and specific citation guidance, and peer-reviewed literature on the DFPT/GW/BSE methodology as implemented in the code.*

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the ABINIT 	Open-source plane-wave/PAW code with strong capabilities in DFPT (density-functional perturbation theory), phonons, and many-body perturbation theory (GW, BSE). Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
