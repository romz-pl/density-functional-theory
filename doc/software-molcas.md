# Molcas / OpenMolcas: An Exhaustive Review

## 1. Overview

MOLCAS/OpenMolcas is an *ab initio* electronic structure software suite whose defining specialty is **multiconfigurational (multireference) quantum chemistry** — methods designed for molecular systems where a single Slater determinant or single Kohn–Sham reference is inadequate (bond breaking, transition metals, lanthanides/actinides, excited states, conical intersections, diradicals, etc.). Alongside this multireference core, the package also implements Hartree–Fock and Kohn–Sham Density Functional Theory (DFT), and — distinctively — hybrid wave-function/DFT approaches such as Multiconfiguration Pair-Density Functional Theory (MC-PDFT).

- **Lineage**: Development traces to the late 1980s under Björn O. Roos at Lund University. The name combines *Mol*ecule (integral code by Jan Almlöf) and *CAS* (Complete Active Space code by Roos).
- **Two branches today**:
  - **MOLCAS** — the historical, licensed (academic-fee) package, still maintained.
  - **OpenMolcas** — since 2019, the free, open-source fork released under the **LGPL v2.1**. It is not a reimplementation but shares the bulk of the MOLCAS codebase; a small number of components remain proprietary and are excluded because their original authors could not release them or chose not to.
- **Community model**: OpenMolcas is developed openly on GitLab/GitHub by a large international consortium (over 100 contributing authors across the flagship papers), with community-driven feature contributions, in contrast to MOLCAS's more centralized, license-fee-funded model.
- **Platforms**: Linux, Unix, macOS (MOLCAS additionally lists Windows support historically); built with CMake; can link to external linear algebra libraries (MKL, OpenBLAS) or use a bundled LAPACK submodule.

## 2. Architecture and Design Philosophy

- **Modular program suite**: [Open]Molcas is not a monolithic binary but a collection of many individual, largely independent programs/modules (e.g., SEWARD for integrals, RASSCF/CASSCF, CASPT2, RASSI, SLAPAF for geometry optimization, etc.), orchestrated by a driver script (`pymolcas`).
- **Inter-module communication** historically occurs through shared files (binary/text), such as `RunFile`, `InpOrb`, and `JobIph`, which store orbitals, wave function, and job-state information between steps of a calculation.
- **Cholesky decomposition (CD) / RI infrastructure**: A defining technical innovation is that essentially every method in the package is reformulated in terms of **Cholesky vectors** rather than conventional two-electron integrals. This drastically reduces storage requirements and lowers the computational scaling of most tensor contractions used in both mean-field and correlated-method energy/property evaluations, and underlies the package's ability to treat comparatively large multiconfigurational active spaces and basis sets efficiently.
- **Developer platform**: The codebase is explicitly designed to be a platform for method developers as well as an end-user tool, which is reflected in the very large number of named contributors (100+) to the flagship overview articles.

## 3. Core Electronic-Structure Methods

### 3.1 Mean-field and DFT
- Hartree–Fock (HF)
- Kohn–Sham DFT (KS-DFT), including a range of exchange-correlation functionals

### 3.2 Multiconfigurational self-consistent field (MCSCF) family — the specialty
- **CASSCF** (Complete Active Space SCF) — the historic core method, using the Super-CI approach with a density-matrix formulation for orbital optimization.
- **RASSCF** (Restricted Active Space SCF) — partitions the active space into RAS1/RAS2/RAS3 subspaces with excitation-level restrictions, extending accessible active-space sizes beyond CASSCF's exponential wall.
- **GASSCF** (Generalized Active Space SCF) — a further generalization allowing arbitrary numbers of active subspaces with flexible occupation constraints; includes the related **SplitGAS** approach.
- **State-Average (SA)** variants of the above for simultaneous treatment of multiple electronic states, including analytic gradients (with density fitting) for SA-CASSCF.

### 3.3 Circumventing exponential active-space scaling
Because exact CASSCF-type diagonalization scales combinatorially with active-space size (practically limiting CASSCF to roughly CAS(18,18) without special techniques), OpenMolcas has pioneered/integrated several advanced strategies:
- **DMRG-SCF** (Density Matrix Renormalization Group as an active-space solver), including relativistic and 4-component DMRG.
- **Stochastic-CASSCF**, based on **Full Configuration Interaction Quantum Monte Carlo (FCIQMC)**, via interface to the NECI code — one- and two-electron reduced density matrices (RDMs) are stochastically sampled and fed back into the Super-CI orbital-relaxation step.
- **Selected CI (SCI)** methods.
- Interfaces to external, specialized codes for large active spaces: **Block**, **CheMPS2**, **QCMaquis**, **Dice**, **GronOR**, and **NECI**.
- A parallel interface for computing 4-particle RDM (4-RDM) elements on clusters, extending feasible active-space sizes for methods that require them (with N ≲ 22 orbitals at publication of the 2020 review, targeted to expand toward N ≲ 30).

### 3.4 Dynamic-correlation / perturbation-theory methods
- **CASPT2** (Complete Active Space second-order Perturbation Theory) — historically one of Molcas's signature contributions to the field, including its **multistate (MS-CASPT2)** and **extended multistate (XMS-CASPT2)** variants.
- **RASPT2** and **GASPT2** — PT2 built on RASSCF/GASSCF references.
- Newer **quasi-degenerate CASPT2 variants** combining the strengths of MS- and XMS-CASPT2 in a single formulation, useful both for accurate relative energies and for near-degenerate regions of potential energy surfaces.
- **σp-CASPT2**: an exponential-regularization scheme for the first-order perturbative amplitudes that suppresses the intruder-state problem with minimal sensitivity to the regularization parameter, replacing older empirical level-shift techniques.
- **NEVPT2** (n-Electron Valence State Perturbation Theory) — including via DMRG references (DMRG-NEVPT2).
- **Frozen Natural Orbital (FNO)** truncation schemes: FNO-CASPT2 and, more recently, **FNO-RASPT2** / **FNO-GASPT2**, which reduce virtual-space cost via orbital deletion criteria robust to intruder-state-like singularities in restricted/generalized active spaces.
- **Large-scale parallel Multireference Configuration Interaction (MRCI)**, including analytic gradients, via interface to the **Columbus** package.
- **Transcorrelated (TC) Hamiltonians** for treating dynamic (cusp) correlation, accessible via imaginary-time-propagated time-dependent DMRG (TD-DMRG) or via TC-FCIQMC.

### 3.5 Multiconfiguration Pair-Density Functional Theory (MC-PDFT) — the wave-function/DFT hybrid
A hallmark "third way" method combining a multiconfigurational wave function (CASSCF/RASSCF/GASSCF or DMRG) with a density functional evaluated on the pair density rather than solely the one-electron density:
- **CAS-PDFT / RAS-PDFT / GAS-PDFT / DMRG-PDFT** variants.
- **Hybrid MC-PDFT (HMC-PDFT)**, mixing translated functionals with wave-function-derived energy components.
- **State-Interaction PDFT (SI-PDFT)** for strongly interacting electronic states.
- Analytic **gradients** for state-specific CAS-PDFT (constructed via a fully variational Lagrangian), enabling geometry optimization at MC-PDFT cost that is markedly lower than CASPT2 while giving comparable equilibrium geometries.
- Translated KS-DFT exchange-correlation functionals adapted for on-top pair-density use, with tunable/scaled exchange-correlation contributions in more recent releases.

### 3.6 State-interaction and property methods
- **RASSI** (RAS State Interaction) module: computes matrix elements (including transition properties, spin–orbit coupling, and Dyson orbitals for ionization processes) between CI/CASSCF/RASSCF wave functions of possibly different symmetry, spin, or active space.
- **MPSSI**: a non-orthogonal state-interaction approach extended to matrix product state (DMRG) wave functions.
- Spin–orbit coupling (SOC) treatments for heavy-element and lanthanide/actinide chemistry.
- Magnetic property calculations: g-tensors, zero-field splitting, exchange couplings (J), magnetic circular dichroism, natural/spin natural orbitals for magnetic-property analysis, exploiting optimized DMRG wave functions among others.

## 4. Relativistic Treatments

- **Douglas–Kroll–Hess (DKH)** scalar-relativistic Hamiltonian, including local and linear-scaling implementations for large systems.
- Spin–orbit coupling via RASSI-based approaches (state-interaction with SOC operators), applicable to lanthanide, actinide, and heavy transition-metal systems.
- Relativistic/four-component DMRG implementations for strongly relativistic multireference problems.

## 5. Potential Energy Surfaces, Gradients, and Dynamics

- **Analytic gradients**: available for SA-CASSCF (with density fitting), state-specific and multistate CASPT2/RASPT2 variants, and state-specific CAS-PDFT; more recently extended to multiscale QM/MM contexts (e.g., MCSCF coupled to polarizable fluctuating-charge force fields).
- **Conical intersection optimization**: dedicated algorithms for locating and characterizing conical intersections and other PES crossing seams, essential for photochemical mechanism studies.
- **Ab initio molecular dynamics**: both adiabatic and **nonadiabatic** dynamics support, including surface-hopping approaches (interface to **SHARC** — Surface Hopping including ARbitrary Couplings).
- **Constrained fragment optimization** and diabatization methods for characterizing electron/energy-transfer pathways.
- Interfaces to external nuclear-dynamics tools for semiclassical and full quantum mechanical nuclear motion simulations.

## 6. Spectroscopy and Light–Matter Interaction

- Simulation of UV/Vis absorption and emission spectra (vertical and vibronic) from multiconfigurational wave functions.
- X-ray spectroscopy: core-ionization and core-excitation methods (e.g., via RASSCF/RASPT2 on localized core holes), simulation of X-ray absorption (NEXAFS) and related processes.
- **Semiclassical light–matter interaction** via Gauss–Hermite quadrature, offering an alternative to the standard multipole expansion for simulating spectroscopic transitions.
- Magnetic circular dichroism (MCD) property calculations.
- Complex Absorbing Potential (CAP) combined with multistate CASPT2 (CAP/MS-CASPT2) for metastable/resonance electronic states (via external interface, e.g., PyOpenCAP), used to extract resonance energies and lifetimes.
- Vibronic spectral simulation workflows, including via the wave function/DMRG-based approaches noted above.

## 7. Multiscale and Embedding Methods

- **QM/MM** interfaces (e.g., coupling to Tinker for the MM part).
- **Frozen-Density Embedding Theory (FDET)** for multiscale simulations embedding a quantum-mechanically treated subsystem in a frozen density environment.
- Coupling of MCSCF wave functions to **polarizable force fields** (e.g., Fluctuating Charges, FQ) with analytic nuclear gradients, used for solvated-system vibronic spectroscopy.

## 8. Basis Sets and Integrals

- Wide range of standard Gaussian basis sets, plus specialized **electronic and muonic basis sets**.
- Density fitting / resolution-of-identity and Cholesky-based auxiliary-basis approaches for two-electron integral approximation, unified under the CD infrastructure described in Section 2.
- Efficient integral-derivative evaluation schemes (e.g., Rys quadrature-based methods) inherited from the historical Molcas codebase and refined over successive releases.

## 9. Analysis, Visualization, and Ecosystem Interfaces

- Built-in and add-on tools for post-calculation wave function analysis (e.g., natural orbitals, binatural orbitals) and visualization.
- Wave function analysis via external libraries (e.g., libwfa) for excited-state character analysis.
- Automated/semi-automated **active-space selection** tools, addressing one of the historically most user-skill-dependent steps of multiconfigurational calculations (from small-basis orbital selection heuristics to schemes such as the approximate pair coefficient, APC, method for larger systems).
- Interfaces to numerous external packages beyond those already mentioned: **Columbus** (large-scale parallel MRCI), **SHARC** (nonadiabatic dynamics), **GronOR** (massively parallel non-orthogonal CI on systems such as molecular crystals, demonstrated at scale on leadership-class HPC systems), and quantum Monte Carlo codes.

## 10. Licensing and Availability

| Aspect | MOLCAS | OpenMolcas |
|---|---|---|
| License | Academic/commercial license fee | LGPL v2.1 (free, open source) |
| Source availability | Restricted | Public (GitLab: `Molcas/OpenMolcas`, mirrored on GitHub) |
| Support model | Vendor/institution-based | Community-driven, web-based |
| Codebase relationship | Superset (includes some components not released) | Large subset of MOLCAS codebase; not a fork/reimplementation |

- **Primary repository**: `https://gitlab.com/Molcas/OpenMolcas` (GitHub is a mirror).
- **Build system**: CMake-based; a `pymolcas` script drives job execution.
- **Documentation**: HTML/PDF manual (much of it predates OpenMolcas's creation and may be partly outdated in details, per the project's own README).

## 11. Typical Application Domains

Based on the applications reviewed across the flagship papers and associated literature:
- Transition-metal and lanthanide/actinide chemistry (spin states, magnetic coupling, single-molecule magnets).
- Photochemistry and photophysics: conical intersections, excited-state proton transfer, isomerization (e.g., retinal chromophores, azobenzene photoswitching).
- X-ray and core-level spectroscopy simulation.
- Strongly correlated molecular clusters (e.g., cuprate-type exchange-coupled systems relevant to high-*T*c superconductivity discussions, magneto-electric coupling).
- Molecular qubit candidates (spin-orbit and spin-decoherence properties of transition-metal complexes).
- Benchmark studies of excited-state methods (e.g., large-scale comparisons of CASSCF/NEVPT2/MC-PDFT/HMC-PDFT excitation energies against reference databases such as QUESTDB).

## 12. Strengths and Limitations (Synthesis)

**Strengths**
- Unmatched breadth of multireference/multiconfigurational methodology within a single open-source package, spanning wave-function theory, DMRG/stochastic approaches, and MC-PDFT hybrids.
- Strong analytic-gradient and PES-exploration toolset for photochemistry-oriented research (conical intersections, nonadiabatic dynamics).
- Cholesky/RI infrastructure gives favorable scaling and memory footprint relative to conventional integral-based implementations.
- Very active, broad developer community (>100 contributing authors on core papers), enabling rapid uptake of new methods (e.g., σp-CASPT2, FNO-RASPT2/GASPT2).
- Rich set of interfaces to specialized external codes (DMRG solvers, FCIQMC, MRCI, nonadiabatic dynamics, non-orthogonal CI at HPC scale) rather than reinventing every capability internally.

**Limitations**
- Exact CASSCF/CASCI-type diagonalization remains combinatorially limited without invoking DMRG/FCIQMC/selected-CI extensions; active-space selection still requires user expertise, notwithstanding automated-selection tools.
- Documentation is acknowledged by the developers themselves to be partly outdated in places, since much of it predates the OpenMolcas fork.
- Not a "black-box" DFT-only package: KS-DFT capabilities exist but are not the primary focus or differentiator relative to dedicated mainstream DFT codes.
- Some legacy MOLCAS components remain proprietary and are unavailable in the OpenMolcas branch.

---

# Key Publications on the Package's Theory and Methodology

### Flagship overview / release papers (chronological)
1. Karlström, G. *et al.* **"MOLCAS: A program package for computational chemistry."** *Comput. Mater. Sci.* **28** (2003) 222. doi:10.1016/s0927-0256(03)00109-5
2. Aquilante, F. *et al.* **"MOLCAS 7: The next generation."** *J. Comput. Chem.* **31** (2010) 224. doi:10.1002/jcc.21318
3. Aquilante, F.; Autschbach, J.; Carlson, R. K.; Chibotaru, L. F.; Delcey, M. G.; De Vico, L.; Fdez. Galván, I.; Ferré, N.; Frutos, L. M.; Gagliardi, L. *et al.* **"MOLCAS 8: New Capabilities for Multiconfigurational Quantum Chemical Calculations across the Periodic Table."** *J. Comput. Chem.* **37** (2016) 506–541. doi:10.1002/jcc.24221
4. Fdez. Galván, I.; Vacher, M.; Alavi, A.; Angeli, C.; Aquilante, F.; Autschbach, J. *et al.* **"OpenMolcas: From Source Code to Insight."** *J. Chem. Theory Comput.* **15** (2019) 5925–5964. doi:10.1021/acs.jctc.9b00532
5. Aquilante, F.; Autschbach, J.; Baiardi, A.; Battaglia, S.; Borin, V. A.; Chibotaru, L. F. *et al.* **"Modern quantum chemistry with [Open]Molcas."** *J. Chem. Phys.* **152** (2020) 214117. doi:10.1063/5.0004835
6. Li Manni, G.; Fdez. Galván, I.; Alavi, A.; Aleotti, F.; Aquilante, F.; Autschbach, J. *et al.* **"The OpenMolcas Web: A Community-Driven Approach to Advancing Computational Chemistry."** *J. Chem. Theory Comput.* **19** (2023) 6933–6991. doi:10.1021/acs.jctc.3c00182

### Foundational method papers cited across the above reviews
7. Roos, B. O.; Taylor, P. R.; Siegbahn, P. E. M. **"A Complete Active Space SCF Method (CASSCF) Using a Density Matrix Formulated Super-CI Approach."** *Chem. Phys.* **48** (1980) 157–173.
8. Siegbahn, P. E. M.; Almlöf, J.; Heiberg, A.; Roos, B. O. **"The Complete Active Space SCF (CASSCF) Method in a Newton–Raphson Formulation with Application to the HNO Molecule."** *J. Chem. Phys.* **74** (1981) 2384–2396.
9. Finley, J.; Malmqvist, P. Å.; Roos, B. O.; Serrano-Andrés, L. **"The multi-state CASPT2 method."** *Chem. Phys. Lett.* **288** (1998) 299–306.
10. Malmqvist, P. Å.; Rendell, A.; Roos, B. O. **"The restricted active space self-consistent-field method, implemented with a split graph unitary group approach."** *J. Phys. Chem.* **94** (1990) 5477–5482.
11. Angeli, C.; Cimiraglia, R.; Malrieu, J.-P. **"Introduction of n-electron valence states for multireference perturbation theory."** *J. Chem. Phys.* **114** (2001) 10252.
12. Manni, G. L.; Carlson, R. K.; Luo, S.; Ma, D.; Olsen, J.; Truhlar, D. G.; Gagliardi, L. **"Multiconfiguration Pair-Density Functional Theory."** *J. Chem. Theory Comput.* **10** (2014) 3669–3680.
13. Knecht, S.; Legeza, Ö.; Reiher, M. **"Communication: Four-component density matrix renormalization group."** *J. Chem. Phys.* **140** (2014) 041101.
14. Kurashige, Y.; Yanai, T. **"Second-order perturbation theory with a density matrix renormalization group self-consistent field reference function: Theory and application to the study of chromium dimer."** *J. Chem. Phys.* **135** (2011) 094104.
15. Freitag, L.; Knecht, S.; Angeli, C.; Reiher, M. **"Multireference Perturbation Theory with Cholesky Decomposition for the Density Matrix Renormalization Group."** *J. Chem. Theory Comput.* **13** (2017) 451–459.
16. Battaglia, S.; Keller, S.; Knecht, S. **"Efficient Relativistic Density-Matrix Renormalization Group Implementation in a Matrix-Product Formulation."** *J. Chem. Theory Comput.* **14** (2018) 2353–2369.
17. Sharma, P.; Bao, J. J.; Truhlar, D. G.; Gagliardi, L. **"Density Matrix Renormalization Group Pair-Density Functional Theory (DMRG-PDFT): Singlet–Triplet Gaps in Polyacenes and Polyacetylenes."** *Chem. Sci.* **10** (2019) 1716–1723.
18. Bao, J. J.; Dong, S. S.; Gagliardi, L.; Truhlar, D. G. **"Automatic Selection of an Active Space for Calculating Electronic Excitation Spectra by MS-CASPT2 or MC-PDFT."** *J. Chem. Theory Comput.* **14** (2018).
19. Dong, S. S.; Gagliardi, L.; Truhlar, D. G. **"Excitation spectra of retinal by multiconfiguration pair-density functional theory."** *Phys. Chem. Chem. Phys.* **20** (2018) 7265–7276.
20. Vancoillie, S.; Malmqvist, P. Å.; Veryazov, V. **"Parallelization of a multiconfigurational perturbation theory."** *J. Comput. Chem.* **34** (2013) 1937.
21. Aquilante, F.; Pedersen, T. B.; Lindh, R. **"Low-cost evaluation of the exchange Fock matrix from Cholesky and density fitting representations of the electron repulsion integrals."** (Cholesky infrastructure foundation.)
22. Taube, A. G.; Bartlett, R. J. **"Frozen natural orbital coupled-cluster theory: Forces and application to decomposition of nitroethane."** *J. Chem. Phys.* (basis for FNO deletion criteria later adapted to FNO-CASPT2/RASPT2/GASPT2.)

*Note:* Items 7–22 represent a representative, non-exhaustive selection of the underlying method-development literature most frequently cross-referenced by the three flagship OpenMolcas release papers (items 4–6); the flagship papers themselves each carry several hundred references and constitute the most complete bibliographic entry points into the full body of Molcas/OpenMolcas-related theoretical literature.

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the Molcas / OpenMolcas 	Quantum chemistry package specializing in multiconfigurational methods, with DFT capabilities included; OpenMolcas is the open-source branch. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
