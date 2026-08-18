# CASTEP: An Exhaustive Technical and Practical Review

**CASTEP** (Cambridge Serial Total Energy Package) is a first-principles plane-wave pseudopotential code for calculating the properties of crystalline solids, surfaces, interfaces, molecules, liquids, and amorphous materials from quantum mechanics, within the framework of density functional theory (DFT). It is jointly developed by the CASTEP Developers' Group (CDG) — a consortium spanning the Universities of Cambridge, Durham, York, and other UK institutions, historically linked to the Theory of Condensed Matter (TCM) group at Cambridge — and is commercially distributed by Dassault Systèmes BIOVIA as part of the *Materials Studio* modelling suite.

---

## 1. Identity, History, and Governance

| Attribute | Detail |
|---|---|
| Full name | Cambridge Serial Total Energy Package |
| Category | Plane-wave pseudopotential DFT electronic-structure code |
| Origin | Lineage traces to the CASTEP code developed from the late 1980s/1990s at Cambridge (M.C. Payne and co-workers), building on the total-energy pseudopotential methods reviewed in Payne *et al.* (1992) |
| Modern rewrite | CASTEP was substantially rewritten (~2000s) in modern, modular Fortran95/2003 by the "new CASTEP" developer consortium (Clark, Segall, Pickard, Hasnip, Probert, Refson, Payne, and others), distinct from the earlier Fortran77 code base |
| Maintained by | CASTEP Developers' Group (CDG), in partnership with Cambridge Enterprise (University of Cambridge) and STFC (Science and Technology Facilities Council), which administers academic licensing via CoSeC |
| Commercial distributor | Dassault Systèmes BIOVIA, as a module of BIOVIA Materials Studio |
| Website | castep.org (community/academic) and 3ds.com/biovia (commercial) |
| Primary reference citation | Clark, S.J., Segall, M.D., Pickard, C.J., Hasnip, P.J., Probert, M.I.J., Refson, K., Payne, M.C. *First principles methods using CASTEP*, Z. Kristallogr. **220**, 567–570 (2005) |

CASTEP exists in two parallel forms that share the same underlying source and physics engine:

1. **Academic/standalone CASTEP** — a command-line Fortran code, source-licensed free of charge to academic research groups worldwide via STFC/CoSeC (superseding the earlier UK Academic Licence and European Source-Code Licence). This licence provides only the command-line interface; it does **not** include the BIOVIA Materials Studio graphical front end, and redistribution/public deployment beyond the licensed research group is restricted.
2. **BIOVIA Materials Studio CASTEP** — the same computational engine wrapped in Materials Studio's graphical "Materials Visualizer" environment, sold commercially by BIOVIA (with module-based licensing, e.g., a separate purchasable NMR CASTEP module for solid-state NMR calculations).

---

## 2. Underlying Physics and Methodology

### 2.1 Core electronic-structure approach
CASTEP solves the Kohn–Sham equations of density functional theory using:

- **Plane-wave basis sets**, exploiting periodic boundary conditions and Bloch's theorem, so the natural domain of application is 3D-periodic crystalline systems (solids, surfaces via slab models, and interfaces); finite/isolated systems are treated with supercell approaches.
- **Pseudopotentials** to represent the electron–ion interaction, avoiding explicit treatment of core electrons:
  - Norm-conserving pseudopotentials
  - Ultrasoft (Vanderbilt-type) pseudopotentials
  - An in-built **"on-the-fly generation" (OTFG)** scheme that constructs pseudopotentials at run time from simple specification strings, avoiding the need for external pseudopotential files; standard OTFG libraries (e.g., the "C19" ultrasoft library) ship as defaults
  - Non-self-consistent projector-augmented-wave (PAW)-style reconstruction for property calculations
- **Self-consistent-field solvers**: block Davidson diagonalisation combined with density mixing; conjugate-gradient direct energy minimisation; and finite-temperature/"ensemble DFT" methods for metals and systems with partial orbital occupancies (important for robust convergence in metallic and small-gap systems).

### 2.2 Exchange–correlation functionals and Hamiltonians
- **LDA**: standard local-density approximation
- **GGA family**: PW91, PBE, RPBE, PBEsol, Wu–Cohen (WC)
- **Meta-GGA**: rSCAN (regularised SCAN)
- **Hybrid functionals**: PBE0 (Hartree–Fock exchange + PBE), B3LYP, screened-exchange sX-LDA, and the HSE family (including user-defined screening/mixing parameters)
- **DFT+U**: Hubbard-corrected functionals for strongly correlated (e.g., transition-metal oxide, f-electron) systems
- **Dispersion corrections**: semi-empirical DFT+D schemes (Tkatchenko–Scheffler, many-body dispersion/MBD, Grimme D2/D3/D4, XDM)
- **LibXC interface**, giving access to the broader library of exchange-correlation functionals maintained by the wider DFT community

### 2.3 Structural, dynamical and transition-state methods
- Full variable-cell geometry optimisation (BFGS, LBFGS, TPSD algorithms), including optimisation in internal (delocalised) coordinates and via damped molecular dynamics
- Transition-state search using the LST/QST (linear/quadratic synchronous transit) method and the nudged elastic band (NEB) method
- Molecular dynamics: fixed- and variable-cell MD; NVE, NVT, NPH, and NPT ensembles; path-integral MD for treating quantum nuclear motion (e.g., proton tunnelling/zero-point effects)

### 2.4 Vibrational, elastic and thermodynamic properties
- Phonon dispersion and density of states via density-functional perturbation theory (DFPT/linear response) and via finite-displacement supercell methods
- Infrared and Raman intensities
- Elastic constants and derived polycrystalline moduli (bulk/shear modulus, Poisson ratio, etc.) via the finite-strain (stress–strain) method, together with piezoelectric and internal-strain tensors

### 2.5 Magnetic resonance (solid-state NMR)
CASTEP includes the **GIPAW** (gauge-including projector-augmented-wave) formalism for computing, under periodic boundary conditions:
- NMR chemical shielding/shift tensors
- Electric field gradient (EFG) tensors (quadrupolar parameters)
- J-couplings
- Hyperfine and *g*-tensors (EPR-relevant parameters)

In Materials Studio, this functionality is packaged as the separately licensed **NMR CASTEP** module.

### 2.6 Spectroscopy: core-level, optical, and electronic structure
- Electron energy-loss/electron energy-loss near-edge structure (EELS/ELNES) and X-ray absorption near-edge structure (XANES) spectra, including core-hole (excited-state) treatments for X-ray photoemission and absorption spectroscopy simulation
- Optical matrix elements and derived optical spectra (complex dielectric function, absorption, reflectivity, refractive index) from interband transition matrix elements
- Electronic band structure and density of states
- Population/bonding analysis: Mulliken and Hirshfeld population analysis
- Electron localisation function (ELF) for real-space bonding visualisation
- Time-dependent DFT (TDDFT) excitation-energy and state-projection analysis

### 2.7 Implicit solvation and surfaces/interfaces
CASTEP supports implicit (continuum) solvation models, allowing automated preparation of vacuum vs. solvated calculations and extraction of solvation free energies — directly relevant to modelling solid–liquid interfaces, in addition to its native strength in modelling free/periodic surfaces via slab geometries and heterostructure interfaces via appropriately constructed supercells.

---

## 3. Software Architecture and Performance

- Written in modular **Fortran 2003**, following the "new CASTEP" rewrite that replaced the older Fortran77 code base, improving maintainability and enabling systematic parallel scaling
- **Parallelism**: hybrid MPI + OpenMP; data and workload are distributed over plane-wave coefficients, k-points, electronic bands, and (for high-throughput property calculations) over independent "property farms" of parallel tasks
- Scales from a single desktop/laptop workstation up to large HPC clusters/supercomputers, with an "intelligent default" parallel decomposition chosen automatically based on the number of processors and k-points
- Companion/auxiliary tools in the ecosystem include **OptaDOS** (advanced density-of-states, ELNES/core-loss and optical-spectra post-processing), **MagresView/Soprano** (NMR tensor visualisation and analysis, reading CASTEP's `.magres` output), and **c2x** (a third-party visualisation/format-conversion tool supporting CASTEP input/output alongside other electronic-structure codes)
- Key native output file types include `.bands` (band structure), `.phonon`/`.phonon_dos` (vibrational data), `.elastic` (elastic/piezoelectric tensors), `.magres` (NMR tensors), `.tddft` (excitation data), `.elnes_bin` (core-loss matrix elements), and `.efield` (IR oscillator strengths/permittivity)

---

## 4. Licensing Model

CASTEP's licensing is bifurcated between free academic source-code access and paid commercial distribution:

- **Academic use**: A worldwide, cost-free "research-group" source-code licence is available via the STFC/CoSeC licence portal, jointly administered with Cambridge Enterprise and the CASTEP Developers' Group. Applicants must demonstrate bona fide academic status and institutional backing (e.g., a university-affiliated signatory), and requests are subject to security/legal vetting; the licence covers CASTEP and NMR-CASTEP source code but is **command-line only** — it does not bundle the BIOVIA Materials Studio GUI. This academic licence superseded older, more restrictive schemes (the former UK Academic Licence and a separate European source-code licence formerly administered by BIOVIA).
- **Commercial use**: Organisations wishing to use CASTEP for commercial/non-academic purposes must license it through **BIOVIA Materials Studio**, purchased from Dassault Systèmes BIOVIA. In Materials Studio, CASTEP functionality is gated by license "features" (e.g., `MS_castep_ui`), with some sub-capabilities (notably NMR calculations) requiring an additional, separately purchased module. Materials Studio also offers auto-trial licensing (typically ~30 days) for limited evaluation.
- Academic computing centres can obtain source access for HPC deployment but may not redistribute it publicly beyond the terms of the licence agreement; direct enquiries are handled via castep@stfc.ac.uk.

---

## 5. Typical Application Domains

Reflecting its plane-wave/periodic-boundary-condition design, CASTEP is most heavily used for:

- **Bulk crystalline solids** — ceramics, semiconductors, metals, minerals, superconductors, and functional oxides: structural, elastic, electronic, vibrational, thermodynamic, and optical property prediction
- **Surfaces and interfaces** — adsorption, surface energetics, heterostructures, and solid–liquid interfaces (via implicit solvation)
- **High-pressure/geophysical mineralogy** — phase stability and equations of state under extreme conditions
- **Spectroscopic interpretation** — connecting computed NMR, EELS/XANES, IR/Raman, and optical spectra to experimental measurements
- **Reaction pathways and defects** — transition-state search (NEB/LST-QST), defect and doping studies in periodic supercells

It is comparatively less suited to strongly non-periodic, large-biomolecular, or linear-scaling problems, where localised-basis or O(N) codes (e.g., ONETEP) are typically preferred.

---

## 6. Key Publications on CASTEP Theory and Methods

The following list gathers the foundational and methodological papers most commonly cited when using CASTEP, grouped by theme.

### 6.1 Primary code/method references
- Clark, S.J., Segall, M.D., Pickard, C.J., Hasnip, P.J., Probert, M.I.J., Refson, K., Payne, M.C. (2005). *First principles methods using CASTEP.* **Zeitschrift für Kristallographie – Crystalline Materials**, 220(5–6), 567–570. https://doi.org/10.1524/zkri.220.5.567.65075
- Segall, M.D., Lindan, P.J.D., Probert, M.J., Pickard, C.J., Hasnip, P.J., Clark, S.J., Payne, M.C. (2002). *First-principles simulation: ideas, illustrations and the CASTEP code.* **Journal of Physics: Condensed Matter**, 14(11), 2717–2744. https://doi.org/10.1088/0953-8984/14/11/301
- Payne, M.C., Teter, M.P., Allan, D.C., Arias, T.A., Joannopoulos, J.D. (1992). *Iterative minimization techniques for ab initio total-energy calculations: molecular dynamics and conjugate gradients.* **Reviews of Modern Physics**, 64(4), 1045–1097. https://doi.org/10.1103/RevModPhys.64.1045

### 6.2 Pseudopotential theory
- Vanderbilt, D. (1990). *Soft self-consistent pseudopotentials in a generalized eigenvalue formalism.* **Physical Review B**, 41, 7892.
- Troullier, N., Martins, J.L. (1991). *Efficient pseudopotentials for plane-wave calculations.* **Physical Review B**, 43, 1993.

### 6.3 Exchange–correlation functionals
- Perdew, J.P., Zunger, A. (1981). *Self-interaction correction to density-functional approximations for many-electron systems.* **Physical Review B**, 23, 5048. (LDA)
- Perdew, J.P., Burke, K., Ernzerhof, M. (1996). *Generalized gradient approximation made simple.* **Physical Review Letters**, 77, 3865–3868. (PBE)
- Perdew, J.P. *et al.* (2008). *Restoring the density-gradient expansion for exchange in solids and surfaces.* **Physical Review Letters**, 100, 136406 / related **Physical Review B**, 79, 155107. (PBEsol)
- Hammer, B., Hansen, L.B., Nørskov, J.K. (1999). *Improved adsorption energetics within density-functional theory using revised Perdew-Burke-Ernzerhof functionals.* **Physical Review B**, 59, 7413. (RPBE)
- Wu, Z., Cohen, R.E. (2006). *More accurate generalized gradient approximation for solids.* **Physical Review B**, 73, 235116. (WC)
- Bartók, A.P., Yates, J.R. (2019). *Regularized SCAN functional.* **Journal of Chemical Physics**, 150, 161101. (rSCAN)
- Becke, A.D. (1993). *A new mixing of Hartree–Fock and local density-functional theories.* **Journal of Chemical Physics**, 98, 1372 / related B3LYP references. (Hybrid functionals)
- Heyd, J., Scuseria, G.E., Ernzerhof, M. (2003, erratum 2006). *Hybrid functionals based on a screened Coulomb potential.* **Journal of Chemical Physics**, 118, 8207; 124, 219906. (HSE family)
- Clark, S.J., Robertson, J. (2010). *Screened exchange density functional applied to solids.* **Physical Review B**, 82, 085208. (sX-LDA)

### 6.4 DFT+U and dispersion corrections
- Anisimov, V.I., Zaanen, J., Andersen, O.K. and related DFT+U formalism papers; O'Regan, D.D. *et al.* on linear-response Hubbard-U in the CASTEP implementation.
- Tkatchenko, A., Scheffler, M. (2009). *Accurate molecular van der Waals interactions from ground-state electron density and free-atom reference data.* **Physical Review Letters**, 102, 073005.
- Grimme, S. (2006, and D3/D4 successors). *Semiempirical GGA-type density functional constructed with a long-range dispersion correction.* **Journal of Computational Chemistry**, 27, 1787.

### 6.5 Structural optimisation and transition-state methods
- Pfrommer, B.G., Côté, M., Louie, S.G., Cohen, M.L. (1997). *Relaxation of crystals with the quasi-Newton method.* **Journal of Computational Physics**, 131, 233–240. (BFGS geometry optimisation)
- Halgren, T.A., Lipscomb, W.N. (1977). *The synchronous-transit method for determining reaction pathways and locating molecular transition states.* **Chemical Physics Letters**, 49, 225. (LST/QST transition-state search)
- Henkelman, G., Jónsson, H. and related papers on the nudged elastic band (NEB) method.

### 6.6 Ensemble/finite-temperature DFT
- Marzari, N., Vanderbilt, D., Payne, M.C. (1997). *Ensemble density-functional theory for ab initio molecular dynamics of metals and finite-temperature insulators.* **Physical Review Letters**, 79, 1337.

### 6.7 Vibrational/phonon theory (DFPT)
- Refson, K., Tulip, P.R., Clark, S.J. (2006). *Variational density-functional perturbation theory for dielectrics and lattice dynamics.* **Physical Review B**, 73, 155114.

### 6.8 Solid-state NMR (GIPAW)
- Pickard, C.J., Mauri, F. (2001). *All-electron magnetic response with pseudopotentials: NMR chemical shifts.* **Physical Review B**, 63, 245101.
- Yates, J.R., Pickard, C.J., Mauri, F. (2007). *Calculation of NMR chemical shifts for extended systems using ultrasoft pseudopotentials.* **Physical Review B**, 76, 024401.
- Profeta, M., Mauri, F., Pickard, C.J. (2003). *Accurate first principles prediction of ¹⁷O NMR parameters in SiO₂: assignment of the zeolite ferrierite spectrum.* **Journal of the American Chemical Society**, 125, 541.

### 6.9 Core-level and optical spectroscopy
- Gao, S.-P., Pickard, C.J., Payne, M.C., Zhu, J., Yuan, J. (2008). *Theory of core-hole effects in 1s core-level spectroscopy of the first-row elements.* **Physical Review B**, 77, 115122. (Core-hole/XANES/ELNES methodology)
- Refson, K. and co-workers, and Yates, J.R. *et al.* on optical matrix-element and dielectric-function implementations underlying CASTEP's optical properties module.

### 6.10 Crystal structure prediction (a major CASTEP-adjacent research theme)
- Pickard, C.J., Needs, R.J. (2011). *Ab initio random structure searching.* **Journal of Physics: Condensed Matter**, 23, 053201. (The "AIRSS" methodology, closely associated with the CASTEP developer community and widely used with CASTEP for structure prediction under pressure.)

---

## 7. Summary Assessment

CASTEP is a mature, widely validated, and continuously developed plane-wave DFT code with a comprehensive property set — spanning structural relaxation, dynamics, elasticity, vibrational spectroscopy, solid-state NMR, and core-level/optical spectroscopy — built around a shared modern Fortran engine. Its dual-track distribution model (free academic source licence vs. commercial BIOVIA Materials Studio bundling) is somewhat distinctive among major plane-wave codes (contrasted with e.g. VASP's purely commercial model or Quantum ESPRESSO's fully open-source model), giving academic users full command-line/HPC access while commercial users obtain a GUI-integrated, professionally supported package with modular licensing for advanced features such as NMR. Its close association with method-defining publications (GIPAW NMR, DFPT lattice dynamics, AIRSS structure prediction) reflects its role as much as a research vehicle for new first-principles methodology as a production simulation tool.

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the CASTEP 	Plane-wave pseudopotential code for solids, surfaces, and interfaces, distributed with/without BIOVIA Materials Studio; commercial with academic licensing. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
