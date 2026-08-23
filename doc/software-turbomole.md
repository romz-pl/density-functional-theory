# TURBOMOLE: An Exhaustive Review

## 1. Overview

TURBOMOLE is a modular, high-performance quantum chemistry program package for *ab initio* electronic structure calculations on molecules, clusters, periodic systems (solids, surfaces, polymers), and solutions. It was originally developed in the late 1980s by Reinhart Ahlrichs and coworkers at the University of Karlsruhe, based on one of the first widely used direct-SCF implementations. Development continued at the University of Karlsruhe and Forschungszentrum Karlsruhe from 1989–2007, after which commercial stewardship and continued development passed to TURBOMOLE GmbH, a spin-off founded by the original academic developers. The code remains a multi-national academic–industrial collaborative project, with a large group of university and research-institute contributors (KIT, TU Berlin, UC Irvine, DTU, University of Augsburg, FSU Jena, and others) feeding new methods into the GmbH-maintained codebase.

TURBOMOLE is distributed under a dual licensing model: free-of-charge (or low-cost) licenses for academic/non-commercial research use, and commercial licenses for industrial users, the latter often bundled through BIOVIA (Dassault Systèmes), which also supplies the TmoleX graphical front end.

## 2. Design Philosophy and Target Hardware

- **Efficiency on modest hardware.** From its inception, TURBOMOLE was designed to run efficient, stable calculations on affordable, widely available hardware — multi-core workstations and small clusters — rather than requiring large-scale HPC allocations.
- **Gaussian-type orbital (GTO) basis.** All methods use segmented-contracted Gaussian basis functions (the Ahlrichs def2/def-family basis sets were developed largely in tandem with the code).
- **Accuracy-to-cost ratio.** The stated design goal is not the highest possible accuracy in isolation, but the best trade-off of accuracy versus computational cost, achieved through aggressive but well-controlled approximations (RI/density fitting, multipole acceleration, Laplace-transform and pair-natural-orbital techniques, integral-direct algorithms).
- **Symmetry exploitation.** TURBOMOLE supports use of all finite point groups (not just Abelian subgroups), which is comparatively rare among mainstream QC codes and can give large speedups (up to roughly two orders of magnitude for high-symmetry systems such as icosahedral clusters).

## 3. Core Electronic-Structure Methods

| Category | Methods available |
|---|---|
| Mean-field | Restricted/unrestricted Hartree–Fock (HF); Kohn–Sham DFT |
| DFT functional classes | LDA, GGA, global and range-separated hybrids, meta-GGA, double hybrids, local hybrid functionals (LHs) with position-dependent exact-exchange admixture |
| Dispersion corrections | DFT-D3, DFT-D4 (also usable with HF) |
| Post-HF / correlation | MP2 (RI-MP2, Laplace-transform SOS/scaled-opposite-spin variants), explicitly correlated (F12) coupled-cluster methods, CCSD, CCSD(T), and reduced-scaling variants (e.g., pair natural orbital and DLPNO-type approaches) |
| Many-body / Green's function | Random Phase Approximation (RPA) total energies, GW and GW–Bethe–Salpeter Equation (BSE) methods for quasiparticle and excited-state/optical properties |
| Excited states | Linear-response TDDFT/TDA, ADC-type and CC-based excited-state methods, RT-TDDFT (real-time propagation, including high-harmonic-generation simulations) |
| Relativistic effects | Effective core potentials (ECPs), scalar-relativistic and two-component (X2C) treatments, spin–orbit coupling for HF/DFT and periodic systems |
| Solvation | COSMO (conductor-like screening model), interfacing to COSMO-RS/COSMOtherm for free energies in solution |
| Embedding | Frozen-density embedding and projection-based embedding schemes coupling DFT with correlated wavefunction methods, molecule-in-molecule, molecule-in-periodic, and periodic-in-periodic embedding |
| Nuclear/other | Nuclear electronic orbital (NEO) methods (under active development), Periodic Electrostatic Embedded Cluster Method (PEECM) |

Typical calculated properties include ground- and excited-state geometries and transition states, vibrational frequencies and IR/Raman/VCD spectra, NMR shifts and coupling constants, UV-Vis/CD/two-photon absorption spectra, polarizabilities and hyperpolarizabilities, and thermodynamic functions.

## 4. Resolution-of-the-Identity (RI) / Density-Fitting Machinery

RI (also called density fitting, DF) is central to TURBOMOLE's speed advantage:

- **RI-J**: fits the electron density in an auxiliary Gaussian basis to accelerate evaluation of the Coulomb (J) term, turning formally $O(N^4)$ integral evaluation into an effectively much cheaper procedure while introducing only small, well-controlled errors (on the order of tens of microhartree per atom with standard auxiliary basis sets).
- **Multipole-Accelerated RI (MARI-J)**: combines RI with a multipole expansion of the long-range part of the Coulomb interaction, extending efficient treatment to much larger systems.
- **RI-K / seminumerical exchange**: analogous fitting/acceleration techniques applied to the Hartree–Fock exchange term and to local hybrid exchange–correlation contributions.
- **RI-MP2 / RI-CC2**: RI is likewise used to accelerate post-HF correlation methods (MP2, CC2, ADC(2)), which is one of the historical hallmarks that made TURBOMOLE distinctive among QC codes in the 1990s–2000s.
- **Continuous Fast Multipole Method (CFMM)**: used in conjunction with density fitting (DF-CFMM) to push Kohn–Sham matrix construction toward near-linear, $O(N)$, scaling for both molecular and periodic systems.
- **Low-memory iterative density fitting**: a variant that trades a modest increase in computation time for up to a ~15-fold reduction in memory footprint, enabling DFT calculations on molecular systems with thousands of atoms on a single workstation.

## 5. Periodic DFT: The Riper Module

For nearly 30 years TURBOMOLE was restricted to finite (molecular/cluster) systems. Periodic boundary conditions were introduced with version 7.0 (2015) through a newly written module, **riper**, rather than by retrofitting the molecular code.

Key characteristics of riper:

- **Unified treatment of dimensionality.** riper treats molecules, chains, polymers, slabs/surfaces, and 3-D periodic solids on an equal footing within the same Gaussian-type-orbital framework, reusing only the most efficient low-level integral routines from the molecular code.
- **Density fitting for periodic Coulomb sums.** Because the Coulomb interaction is long-ranged, periodic systems require special density-fitting techniques distinct from the molecular case; riper implements a projection-based extension of the molecular RI/DF scheme, partitioned into near-field (real-space, fitted) and far-field (multipole-expanded) contributions, again via a DF-CFMM-type approach operating entirely in real space (not reciprocal/plane-wave space).
- **Scaling.** The combination of RI, CFMM, and a hierarchical numerical integration scheme for the exchange–correlation term targets $O(N)$ scaling for Kohn–Sham matrix formation and nuclear gradients; unit cells with several hundred to ~640 atoms and up to ~19,000 GTO basis functions have been demonstrated.
- **Low-memory RI variant with iterative solvers.** A preconditioned conjugate-gradient-based, low-memory modification of RI is available in riper, again substantially reducing memory demand relative to the standard RI implementation.
- **Extended functionality.** riper has since been extended to periodic Hartree–Fock and hybrid-functional exchange (including two-component/spin–orbit variants), relativistic effective core potentials, real-time TDDFT for periodic systems, and DFT-based embedding (molecule-in-periodic and periodic-in-periodic).
- **Tooling.** A community-built web application, *RIPER-Tools* (developed by M. Sharma with contributions from Y.-F. Chen, under M. Sierka's group), assists with generating riper input files and retrieving crystal structures for periodic DFT setup.

## 6. Parallelization and Performance

TURBOMOLE supports several parallelization models depending on the module:

- **Fork-SMP**: simple shared-memory parallelism via forked processes for some modules.
- **OpenMP**: thread-based shared-memory parallelism.
- **MPI**: distributed-memory parallelism for multi-node clusters.
- **Hybrid OpenMP/MPI**: combined shared- and distributed-memory parallelism for selected modules.
- **GPU acceleration**: available for a subset of modules/methods in recent versions.

Not every module supports every parallelization mode; the manual documents which combinations apply to which executable (e.g., `dscf`, `ridft`, `riper`, `ricc2`, etc.).

## 7. Software Architecture

TURBOMOLE is organized as a set of independent, single-purpose command-line executables sharing a common set of input/output files (notably the `control` file), rather than one monolithic binary — a modular UNIX-philosophy design dating to its 1989 origins. Representative modules include:

- `define` — interactive setup of geometry, basis set, and calculation type, producing the `control` file
- `dscf` / `ridft` — conventional and RI-accelerated HF/DFT SCF energies
- `grad` / `rdgrad` — analytical gradients
- `aoforce` / `NumForce` — analytical/numerical Hessians and vibrational frequencies
- `escf` / `egrad` — excited-state (TDDFT/TDA) energies and gradients
- `ricc2` — RI-based MP2, CC2, ADC(2), CCSD, CCSD(T), and related correlation methods
- `riper` — periodic DFT/HF
- `mpshift` — NMR chemical shifts
- Auxiliary utilities for geometry optimization (`jobex`), COSMO-RS interfacing, and visualization export

A graphical front end, **TmoleX** (developed and distributed by BIOVIA/Dassault Systèmes), provides input generation, job submission/monitoring, and results visualization on top of the native TURBOMOLE executables, including a native Windows-compatible build of the core package.

## 8. Licensing, Distribution, and Community

- **Licensing model**: Dual academic/commercial licensing administered by TURBOMOLE GmbH; commercial access and the TmoleX GUI are distributed via BIOVIA (Dassault Systèmes).
- **Source availability**: Not open source; distributed as licensed binaries/source to license holders.
- **Documentation**: An extensive official manual (several hundred pages, versioned per release, e.g., "Turbomole Manual 7.9") documents theory, keywords, and a bibliography of underlying method papers; the TURBOMOLE website also hosts white papers, HOWTOs (e.g., on the Periodic Electrostatic Embedded Cluster Method), and links to instructional material on TDDFT and post-HF methods for large systems.
- **Release cadence**: Versioned major releases (e.g., V7.0 in 2015 introducing periodic DFT; V7.2 introducing local hybrids; V7.7 covered in the 2023 "Today and Tomorrow" perspective; V8.0 released in 2026 with expanded two-component spin–orbit and current-dependent meta-GGA/local-hybrid response theory, among other changes).
- **User community**: Widely used in academic and industrial computational chemistry/materials science — catalysis (homogeneous and heterogeneous), inorganic and organic chemistry, spectroscopy, biochemistry, and (via riper) surface science and solid-state materials modeling.

## 9. Typical Application Domains

- Catalyst structure/reactivity modeling in homogeneous and heterogeneous catalysis (including cluster models of surfaces and full periodic slab models via riper)
- Reaction path and transition-state searches, potential energy surface scans
- Prediction of vibrational, UV-Vis, CD, VCD, and NMR spectra for structure elucidation
- Solvation free energies via COSMO/COSMO-RS for reaction and partitioning studies
- Organic electronics and photophysics (e.g., OLED emitter design, excited-state and non-radiative decay studies)
- Nonlinear optics and strong-field phenomena via RT-TDDFT (e.g., high-harmonic generation)
- Low-dimensional nanostructures, surfaces, and bulk crystalline materials via periodic DFT (Riper)

## 10. Strengths and Limitations

**Strengths**
- Mature, extensively validated RI/density-fitting infrastructure giving strong accuracy-per-CPU-hour on standard hardware
- Full finite-point-group symmetry exploitation
- Unified Gaussian-basis treatment spanning molecules through 3-D periodic solids in a single package (riper)
- Broad coverage of correlated wavefunction methods (RI-MP2, CC2, CCSD(T), explicitly correlated F12, RPA/GW-BSE) alongside DFT
- Long-standing, well-documented theoretical basis with a large, traceable primary-literature bibliography

**Limitations (as stated by the developers)**
- Not intended to cover highly parameterized semi-empirical methods
- Does not target the highest-end correlation treatments for spectroscopic-accuracy benchmarking
- Limited support for most multi-reference methods and model Hamiltonians
- Commercial/dual license (not open source), which can be a barrier relative to freely redistributable codes
- Plane-wave-based periodic codes (e.g., VASP) remain more established for some classes of bulk/metallic solid-state problems, though riper's GTO-based approach is competitive for many molecular-crystal, surface, and hybrid-functional periodic applications

## 11. Version History Highlights

| Year | Milestone |
|---|---|
| 1989 | Original direct-SCF TURBOMOLE code published (Ahlrichs, Bär, Häser, Horn, Kölmel) |
| 1990s | Introduction of RI-J (Coulomb) resolution-of-the-identity approximation |
| 2003 | Multipole-accelerated RI (MARI-J) |
| 2007 | Transition from University of Karlsruhe/Forschungszentrum Karlsruhe stewardship to TURBOMOLE GmbH |
| ~2013 | Explicitly correlated (F12) coupled-cluster methods added |
| 2015 (V7.0) | Periodic boundary conditions introduced via the new riper module |
| V7.2 | Local hybrid functionals (LHs) added |
| 2020 | "Modular program suite" reference publication (J. Chem. Phys. 152, 184107) |
| 2023 (~V7.7) | "Today and Tomorrow" perspective; NEO methods, HF-based adiabatic connection models, simplified TDDFT under development |
| 2024 | Two-component spin–orbit periodic exchange in riper; Scalmani–Frisch current-dependent meta-GGA/local-hybrid response formalism |
| 2026 (V8.0) | Major release: updated COSMO cavity construction defaults, expanded two-component/spin–orbit formalism across functional classes, current-dependent response theory in `mpshift` and RT-TDDFT |

---

## 12. Key Theory Publications

### 12.1 Foundational and Reference Papers

1. Ahlrichs, R.; Bär, M.; Häser, M.; Horn, H.; Kölmel, C. **"Electronic structure calculations on workstation computers: The program system TURBOMOLE."** *Chem. Phys. Lett.* **1989**, *162*, 165–169.
2. Furche, F.; Ahlrichs, R.; Hättig, C.; Klopper, W.; Sierka, M.; Weigend, F. **"Turbomole."** *WIREs Comput. Mol. Sci.* **2014**, *4*, 91–100. https://doi.org/10.1002/wcms.1162
3. Balasubramani, S. G.; Chen, G. P.; Coriani, S.; Diedenhofen, M.; Frank, M. S.; Franzke, Y. J.; Furche, F.; Grotjahn, R.; Harding, M. E.; Hättig, C.; Hellweg, A.; Helmich-Paris, B.; Holzer, C.; Huniar, U.; Kaupp, M.; Marefat Khah, A.; Karbalaei Khani, S.; Müller, T.; Mack, F.; Nguyen, B. D.; Parker, S. M.; Perlt, E.; Rappoport, D.; Reiter, K.; Roy, S.; Rückert, M.; Schmitz, G.; Sierka, M.; Tapavicza, E.; Tew, D. P.; van Wüllen, C.; Voora, V. K.; Weigend, F.; Wodyński, A.; Yu, J. M. **"TURBOMOLE: Modular program suite for *ab initio* quantum-chemical and condensed-matter simulations."** *J. Chem. Phys.* **2020**, *152*, 184107. https://doi.org/10.1063/5.0004635
4. Franzke, Y. J.; Holzer, C.; Andersen, J. H.; Begušić, T.; Bruder, F.; Coriani, S.; Della Sala, F.; Fabiano, E.; Fedotov, D. A.; Fürst, S.; Gillhuber, S.; Grotjahn, R.; Kaupp, M.; Kehry, M.; Krstić, M.; Mack, F.; Majumdar, S.; Nguyen, B. D.; Parker, S. M.; Pauly, F.; Pausch, A.; Perlt, E.; Phun, G. S.; Rajabi, A.; Rappoport, D.; Samal, B.; Schrader, T.; Sharma, M.; Tapavicza, E.; Treß, R. S.; Voora, V.; Wodyński, A.; Yu, J. M.; Zerulla, B.; Furche, F.; Hättig, C.; Sierka, M.; Tew, D. P.; Weigend, F. **"TURBOMOLE: Today and Tomorrow."** *J. Chem. Theory Comput.* **2023**, *19*, 6859–6890. https://doi.org/10.1021/acs.jctc.3c00347

### 12.2 Periodic DFT / Riper Module

5. Burow, A. M.; Sierka, M.; Mohamed, F. **"Resolution of identity approximation for the Coulomb term in molecular and periodic systems."** *J. Chem. Phys.* **2009**, *131*, 214101. https://doi.org/10.1063/1.3267858
6. Burow, A. M.; Sierka, M. **"Linear scaling hierarchical integration scheme for the exchange–correlation term in molecular and periodic systems."** *J. Chem. Theory Comput.* **2011**, *7*, 3097–3104. https://doi.org/10.1021/ct200412r
7. Łazarski, R.; Burow, A. M.; Sierka, M. **"Density Functional Theory for Molecular and Periodic Systems Using Density Fitting and Continuous Fast Multipole Methods."** *J. Chem. Theory Comput.* **2015**, *11*, 3029–3041. https://doi.org/10.1021/acs.jctc.5b00252
8. Grajciar, L. **"Low-memory iterative density fitting."** *J. Comput. Chem.* **2015**, *36*, 1521–1535. https://doi.org/10.1002/jcc.23961
9. Sharma, M.; Franzke, Y. J.; Holzer, C.; Pauly, F.; Sierka, M. **"Density Functional Theory for Molecular and Periodic Systems in TURBOMOLE: Theory, Implementation, and Applications."** *J. Phys. Chem. A* **2025**, *129*, 9062–9083. https://doi.org/10.1021/acs.jpca.5c02937

### 12.3 Resolution-of-the-Identity (RI) / Density Fitting — Molecular

10. Sierka, M.; Hogekamp, A.; Ahlrichs, R. **"Fast evaluation of the Coulomb potential for electron densities using multipole accelerated resolution of identity approximation."** *J. Chem. Phys.* **2003**, *118*, 9136–9148.
11. Weigend, F. **"Accurate Coulomb-fitting basis sets for H to Rn."** *Phys. Chem. Chem. Phys.* **2006**, *8*, 1057–1065. (RI-J auxiliary basis sets underpinning TURBOMOLE's RI machinery.)

### 12.4 Relativistic Effects and Spin–Orbit Coupling

12. Armbruster, M. K.; Weigend, F.; van Wüllen, C.; Klopper, W. **"Self-consistent treatment of spin–orbit interactions with efficient Hartree–Fock and density functional methods."** *Phys. Chem. Chem. Phys.* **2008**, *10*, 1748–1756.
13. (2024) Two-component spin–orbit Fock exchange for HF, global hybrids, and range-separated hybrids under periodic boundary conditions in riper. *Phys. Rev. B* **2024**, *109*, 165144. https://doi.org/10.1103/PhysRevB.109.165144
14. (2024) Scalmani–Frisch formalism with current-dependent meta-GGAs and local hybrids in `ridft`, `rdgrad`, and `riper`. *J. Chem. Phys.* **2024**. https://doi.org/10.1063/5.0246433

### 12.5 Correlated Wavefunction Methods

15. Hättig, C.; Weigend, F. **"CC2 excitation energy calculations on large molecules using the resolution of the identity approximation."** *J. Chem. Phys.* **2000**, *113*, 5154–5161.
16. Hättig, C. **"Geometry optimizations with the coupled-cluster model CC2 using the resolution-of-the-identity approximation."** *J. Chem. Phys.* **2003**, *118*, 7751–7761.
17. Schmitz, G.; Hättig, C.; Tew, D. P. **"Explicitly correlated PNO-MP2 and PNO-CCSD and their application to the S66 set and large molecular systems."** *Phys. Chem. Chem. Phys.* **2014**, *16*, 22167–22178.

### 12.6 Excited States / TDDFT / GW-BSE

18. Furche, F.; Ahlrichs, R. **"Adiabatic time-dependent density functional methods for excited state properties."** *J. Chem. Phys.* **2002**, *117*, 7433–7447.
19. Kehry, M.; Franzke, Y. J.; Holzer, C.; Klopper, W. **"Quasirelativistic two-component core excitations and polarizabilities from a damped-response formulation of the Bethe–Salpeter equation."** *Mol. Phys.* **2020**, *118*, e1755064.

### 12.7 Local Hybrid Functionals

20. Kaupp, M.; Bahmann, H.; Arbuznikov, A. V. **"Local hybrid functionals: An assessment for thermochemical kinetics."** *J. Chem. Phys.* **2007**, *127*, 194102. (Foundational local-hybrid methodology later implemented efficiently in TURBOMOLE from V7.2 onward via seminumerical integration.)

### 12.8 Solvation (COSMO / COSMO-RS)

21. Klamt, A.; Schüürmann, G. **"COSMO: A new approach to dielectric screening in solvents with explicit expressions for the screening energy and its gradient."** *J. Chem. Soc., Perkin Trans. 2* **1993**, 799–805.
22. Klamt, A. **"Conductor-like Screening Model for Real Solvents: A New Approach to the Quantitative Calculation of Solvation Phenomena."** *J. Phys. Chem.* **1995**, *99*, 2224–2235. (Basis of the COSMO-RS methodology later interfaced with TURBOMOLE via COSMOtherm.)

### 12.9 Companion/Reference Bibliography

23. TURBOMOLE GmbH. **"Turbomole Program Package for *ab initio* Electronic Structure Calculations"** — the official TURBOMOLE user manual, which contains a numbered, comprehensive bibliography (currently >150 entries) of every method paper underlying implemented features; consult the current manual PDF distributed with each release for the complete, version-matched list.

*Note: entries 10, 11, 15–20, and 22 are well-established, canonical papers underlying TURBOMOLE's implementations that are consistently cited in the official manual bibliography and the review articles above; consult the official TURBOMOLE manual's numbered reference list for the complete and version-exact set of theory citations, as it is updated with every release.*

---

**Primary sources consulted:** turbomole.org (official site, feature list, release notes, documentation page), the 2020 *J. Chem. Phys.* reference paper, the 2023 "Today and Tomorrow" *JCTC* perspective, the 2025 *J. Phys. Chem. A* Riper/periodic-DFT review, the official TURBOMOLE manual (v7.9), and the BIOVIA/Dassault Systèmes TURBOMOLE product pages.

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the TURBOMOLE 	Commercial/academic quantum chemistry package known for fast, efficient DFT and post-HF calculations via RI/density-fitting techniques; also supports periodic DFT (Riper module). Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
