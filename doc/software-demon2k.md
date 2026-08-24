# deMon2k: An Exhaustive Review

## 1. Overview

**deMon2k** ("densité de Montréal") is a specialized quantum chemistry program for **Density Functional Theory (DFT)** simulations, built around the **Linear Combination of Gaussian-Type Orbitals (LCGTO)** approach for solving the self-consistent Kohn–Sham equations. Its defining methodological feature is **Auxiliary Density Functional Theory (ADFT)**: instead of computing the full four-center electron repulsion integrals (ERIs) needed in conventional Kohn–Sham DFT, deMon2k variationally fits the Kohn–Sham (and, for hybrid functionals, Hartree–Fock exchange) density using an auxiliary Gaussian function basis. This reduces the Coulomb-term scaling from formally $\mathcal{O}(N^4)$ to $\mathcal{O}(N^3)$ (asymptotically closer to linear with further approximations), giving substantial speed gains with limited loss of accuracy — hence its particular suitability for **large molecular and nanoscale systems** (hundreds to thousands of atoms).

The code descends from the original **deMon** program developed by Alain St-Amant and Dennis R. Salahub at the Université de Montréal (1990), later merged with the **AllChem** code (Hannover) to form deMon2k in the mid-2000s. Development is coordinated primarily by groups at **Cinvestav (Mexico City)**, the **University of Calgary**, and collaborating laboratories in France, Germany, and elsewhere, under the umbrella of the **International deMon Developers Community (deMonCity)**.

---

## 2. Core Theoretical Framework

### 2.1 Auxiliary Density Functional Theory (ADFT)

- The electron density is expanded in an **auxiliary Gaussian ("Hermite") basis set**, distinct from the atomic-orbital basis, following the original Sambe–Felton density-fitting idea later formalized by Dunlap, Connolly, and Sabin.
- The **Coulomb (Fock) potential** is obtained through **variational fitting**, avoiding costly four-center ERIs and instead requiring only two- and three-center integrals.
- Crucially, deMon2k uses the **fitted (auxiliary) density directly** — not just for the Coulomb term but also for evaluating the **exchange–correlation potential and energy** on the numerical integration grid — which is the central ADFT innovation distinguishing it from earlier semi-numerical fitting schemes (which were prone to instabilities in energy derivatives).
- Special **primitive Hermite Gaussian auxiliary function sets** are used, with automatic, adaptive generation of auxiliary basis sets (GEN-A2*, etc.) tuned to the orbital basis.

### 2.2 Auxiliary Density Perturbation Theory (ADPT)

An extension of ADFT to analytic derivative theory, allowing efficient, stable computation of response properties (polarizabilities, NMR shieldings, vibrational frequencies) using perturbation-dependent auxiliary functions rather than finite-difference approaches.

### 2.3 Electronic Structure Methods Available

- Kohn–Sham DFT (LCGTO-KS) and ADFT
- Hartree–Fock theory
- Restricted, unrestricted, restricted open-shell, and constrained-unrestricted SCF formalisms
- Large functional library: LDA, GGA, meta-GGA, and hybrid exchange-correlation functionals, plus an interface to the **Libxc** functional library
- Empirical dispersion corrections (Grimme-type, DFT-D) for essentially the full periodic table
- Effective core potentials (ECP) and model core potentials (MCP), evaluated via **half-numerical integral recurrence relations without angular-momentum (*l*) limitation**
- Analytical three-center ERI recurrence relations with **double asymptotic expansion** for long-range/large-system efficiency

---

## 3. Key Computational Capabilities

| Category | Capabilities |
|---|---|
| **Energies & gradients** | QM, MM, and hybrid QM/MM energy models within one package; analytic energy gradients |
| **Structure search** | Local geometry optimization and transition-state search via (quasi-)Newton restricted-step algorithms; double-ended saddle interpolation for non-intuitive TS searches; intrinsic reaction coordinate (IRC) calculations |
| **Dynamics** | Born–Oppenheimer molecular dynamics (BOMD, NVE/NVT); QM/MM molecular dynamics; Ehrenfest non-adiabatic dynamics; Real-Time TD-ADFT for attosecond electron dynamics (private-version feature) |
| **Excited states / spectroscopy** | Time-dependent ADFT (TD-ADFT) for UV/Vis absorption spectra; X-ray absorption and emission spectra |
| **Response properties** | ADPT- and NIA-CPKS–based static/dynamic polarizabilities and hyperpolarizabilities (α, β); Fukui functions |
| **Magnetic properties** | NMR chemical shielding and magnetizability via ADFT-GIAO; rotational g-tensors and nuclear spin–rotation constants |
| **Thermochemistry** | Thermodynamic functions from the polyatomic ideal-gas approximation and from BOMD simulations |
| **Population/topological analysis** | Mulliken, Löwdin, NBO, Bader (QTAIM), Voronoi, Becke, and Hirshfeld population schemes; topological analysis of electron density and electrostatic potential; molecular field visualization |
| **QM/MM** | Additive QM/MM with non-polarizable and polarizable force fields (e.g., Drude model); link-atom schemes; implicit continuum treatment of remote environment; interfaces to CHARMM and to the Cuby wrapper; metadynamics (private-version) |
| **Special models** | Constrained (unrestricted) DFT and constrained ADFT for charge/spin localization and charge-transfer studies; multicomponent ADFT for nuclear quantum effects; cyclic cluster model for periodic/extended systems; σ–π energy separation analysis |
| **Basis sets** | DFT-optimized orbital basis sets (e.g., DZVP-GGA) and automatically generated auxiliary function sets |
| **Parallelization** | Fully MPI-parallelized, portable across computing platforms and operating systems |
| **Visualization interfaces** | Molden, Molekel, Vu, XAIM |

Several advanced features (e.g., constrained ADFT, multicomponent ADFT, range-separated functionals, real-time TD-ADFT, revised local hardness kernels, dispersion-corrected variants, polarizable QM/MM) exist primarily in **private/development branches** of the code, accessible to registered members of the deMon developer community.

---

## 4. Performance Characteristics for Large Systems

- The variational density-fitting scheme replaces the dominant $\mathcal{O}(N^4)$ four-center ERI bottleneck with cheaper two- and three-center integral evaluation, and the direct use of the fitted density for exchange-correlation evaluation keeps the numerical-integration cost scaling favorably as system size grows.
- The **LDF-EXX** (linear-scaling Hartree–Fock/hybrid exchange) implementation allows hybrid-functional QM(/MM) calculations on systems of several hundred atoms — e.g., reported SCF iteration times of only a few minutes for a ~250-atom polyalanine α-helix using paired auxiliary function sets (a smaller set for SCF iterations, a larger one for the final non-self-consistent energy evaluation).
- Fully optimized structures without symmetry constraints have been demonstrated for systems as large as **C₅₄₀ fullerenes** and Al-zeolite models using all-electron basis sets, within practical wall-clock times on parallel hardware.
- MPI parallelization supports scaling across many cores/nodes, and adaptive atom-centered integration grids balance accuracy against cost for exchange-correlation evaluation.

---

## 5. Program Ecosystem, Branches, and Availability

- **License/Access**: deMon2k is distributed under specific licensing terms via the official website (registration and agreement to distribution restrictions are required); precompiled binaries and source access are available depending on version.
- **Related/derived codes**:
  - **StoBe** (Stockholm–Berlin version, focused on surface science/core-level spectroscopy)
  - **deMon@Grenoble** (local development branch)
  - **deMon2k-KSCED** (Geneva; orbital-free embedding / subsystem DFT extensions following the Wesolowski–Warshel formalism)
  - **deMonNano** (tight-binding/semi-empirical DFT companion code for larger systems)
- **Interfaces**: CHARMM (QM/MM), Cuby (QM/MM wrapper), PUPIL, Libxc (functional library), and multiple visualization packages.
- **Documentation**: An online user's guide/manual and tutorials are maintained on the official site, along with a hands-on DFT workbook designed for teaching at the Master's/advanced-undergraduate level.
- **Community**: Development and support are coordinated through the deMon2k mailing list and periodic "deMon Developers" meetings.

---

## 6. Required Citation

Per the developers, any publication using deMon2k results must cite the version-appropriate program reference, e.g., for the current public/development releases:

> A.M. Köster, G. Geudtner, A. Alvarez-Ibarra, P. Calaminici, M.E. Casida, J. Carmona-Espindola, V.D. Dominguez, R. Flores-Moreno, G.U. Gamboa, A. Goursot, T. Heine, A. Ipatov, A. de la Lande, F. Janetzko, J.M. del Campo, D. Mejia-Rodriguez, J.U. Reveles, J. Vasquez-Perez, A. Vela, B. Zuniga-Gutierrez, and D.R. Salahub, *deMon2k*, Version 5/6, The deMon developers, Cinvestav, Mexico City (2018).

Specific functionals, basis sets, auxiliary function sets, and program features used must also be cited individually (see reference lists below and on the official site).

---

## 7. Publications Related to the Package's Theory

### 7.1 Foundational DFT Theory (Underlying deMon2k)

1. P. Hohenberg, W. Kohn, "Inhomogeneous Electron Gas," *Phys. Rev.* **136**, B864 (1964).
2. W. Kohn, L.J. Sham, "Self-Consistent Equations Including Exchange and Correlation Effects," *Phys. Rev.* **140**, A1133 (1965).
3. R.M. Dreizler, E.K.U. Gross, *Density Functional Theory* (Springer-Verlag, Berlin, 1990).
4. S.H. Vosko, L. Wilk, M. Nusair, "Accurate spin-dependent electron liquid correlation energies," *Can. J. Phys.* **58**, 1200 (1980).
5. A.D. Becke, "Density-functional exchange-energy approximation with correct asymptotic behavior," *Phys. Rev. A* **38**, 3098 (1988).
6. J.P. Perdew, "Density-functional approximation for the correlation energy of the inhomogeneous electron gas," *Phys. Rev. B* **33**, 8822 (1986).

### 7.2 Density-Fitting / Auxiliary Basis Origins

7. H. Sambe, R.H. Felton, "A new computational approach to Slater's SCF-Xα equation," *J. Chem. Phys.* **62**, 1122 (1975).
8. B.I. Dunlap, J.W.D. Connolly, J.R. Sabin, "On some approximations in applications of Xα theory," *J. Chem. Phys.* **71**, 3396 (1979).
9. B.I. Dunlap, J.W.D. Connolly, J.R. Sabin, "On first-row diatomic molecules and local density models," *J. Chem. Phys.* **71**, 4993 (1979).
10. J.W. Mintmire, B.I. Dunlap, "Fitting the Coulomb potential variationally in linear-combination-of-atomic-orbitals density-functional calculations," *Phys. Rev. A* **25**, 88 (1982).
11. J.W. Mintmire, J.R. Sabin, S.B. Trickey, "Local-density-functional methods in two-dimensionally periodic systems," *Phys. Rev. B* **26**, 1743 (1982).

### 7.3 Core ADFT / ADPT Methodology (deMon2k-Specific)

12. A.M. Köster, J.U. Reveles, J.M. del Campo, "Calculation of exchange-correlation potentials with auxiliary function densities," *J. Chem. Phys.* **121**, 3417–3424 (2004).
13. R. Flores-Moreno, A.M. Köster, "Auxiliary density perturbation theory," *J. Chem. Phys.* **128**, 134105 (2008).
14. A.M. Köster, J.M. del Campo, F. Janetzko, B. Zuñiga-Gutierrez, "A MinMax self-consistent-field approach for auxiliary density functional theory," *J. Chem. Phys.* **130**, 114106 (2009).
15. D. Mejia-Rodriguez, A.M. Köster, "Robust and efficient variational fitting of Fock exchange," *J. Chem. Phys.* **141**, 124114 (2014).
16. B. Zuñiga-Gutierrez, A.M. Köster, "Analytical GGA exchange-correlation kernel calculation in auxiliary density functional theory," *Mol. Phys.* **114**, 1026–1035 (2016).
17. L.I. Hernández-Segura, F.A. Olvera-Rubalcava, R. Flores-Moreno, P. Calaminici, A.M. Köster, "Exchange-correlation kernel for perturbation dependent auxiliary functions in auxiliary density perturbation theory," *J. Mol. Model.* (2024). DOI: 10.1007/s00894-024-06091-z

### 7.4 Program Overview / Review Articles

18. G. Geudtner, P. Calaminici, J. Carmona-Espindola, J.M. del Campo, V.D. Dominguez-Soria, R. Flores-Moreno, G.U. Gamboa, A. Goursot, A.M. Köster, J.U. Reveles, T. Mineva, J.M. Vasquez-Perez, A. Vela, B. Zuñiga-Gutierrez, D.R. Salahub, "deMon2k," *WIREs Comput. Mol. Sci.* **2**, 548–555 (2012). (Primary program overview article.)
19. D.R. Salahub, A. Goursot, J. Weber, A.M. Köster, A. Vela, "Applied Density Functional Theory and the deMon Codes: 1964–2004," in *Theory and Applications of Computational Chemistry: The First 40 Years*, eds. C.E. Dykstra, G. Frenking, K.S. Kim, G.E. Scuseria (Elsevier, 2005).
20. D.R. Salahub, S.Yu. Noskov, B. Lev, R. Zhang, V. Ngo, A. Goursot, P. Calaminici, A.M. Köster, A. Alvarez-Ibarra, D. Mejía-Rodríguez, J. Řezáč, F. Cailliez, A. de la Lande, "QM/MM Calculations with deMon2k," *Molecules* **20**, 4780–4812 (2015).
21. M.E. Casida, A. Alvarez-Ibarra, P. Calaminici, D. Mejía-Rodríguez, L. López-Sosa, G. Geudtner, I. Navizet, C. Garcia-Iriepa, D.R. Salahub, "Molecular Simulations with in-deMon2k QM/MM, a Tutorial-Review," *Molecules* **24**, 1653 (2019).
22. J.D. Samaniego-Rojas et al., "QM/MM with Auxiliary DFT in deMon2k," Chapter 1 in *Multiscale Dynamics Simulations: Nano and Nano-bio Systems in Complex Environments*, Royal Society of Chemistry (2021).
23. A. Ostojić et al. (deMon2k collaborators), "Current status of deMon2k for the investigation of the early stages of matter irradiation by time-dependent DFT approaches," *Eur. Phys. J. Spec. Top.* **232**, 2167 (2023).

### 7.5 Historically Significant Method Developments Originating in deMon/deMon-KS

24. V.G. Malkin, O.L. Malkina, M.E. Casida, D.R. Salahub, "Absolute NMR shielding constants and chemical shifts calculated using the deMon program," *J. Am. Chem. Soc.* **116**, 5898 (1994). (Foundational GIAO-DFT NMR implementation.)
25. C. Jamorski, M.E. Casida, D.R. Salahub, "Dynamic polarizabilities and excitation spectra from a molecular implementation of time-dependent density-functional response theory," *J. Chem. Phys.* **104**, 5134 (1996). (Foundational TD-DFT implementation in deMon.)

### 7.6 Selected PhD Theses on deMon2k Theory and Implementation

| Thesis | Author | Year |
|---|---|---|
| Entwicklung einer LCGTO-Dichtefunktionalmethode mit Hilfsfunktionen | A.M. Köster | 1998 |
| Geometry Optimization in LCGTO-DFT Methods with Auxiliary Functions | J.U. Reveles Ramírez | 2004 |
| Separación de las Energías σ y π en la Teoría de los Funcionales de la Densidad | Z. Gómez Sandoval | 2005 |
| Analytic Derivatives in LCGTO-DFT Pseudo-Potential Methods with Auxiliary Functions | R. Flores Moreno | 2006 |
| Exploring Chemical Reactivity with Auxiliary Density Functional Theory | J. Martín del Campo Ramírez | 2008 |
| Cálculos de Sistemas Grandes Basados en la Teoría de Funcionales de la Densidad Auxiliar | V.D. Domínguez Soria | 2009 |
| Magnetic Shielding, Spin–Spin Coupling Constant and Magnetizability Tensors from ADPT | B.A. Zúñiga Gutiérrez | 2011 |
| Born–Oppenheimer Molecular Dynamics with Auxiliary Density Functional Theory | G.U. Gamboa Martínez | 2011 |
| Propiedades Termodinámicas de Sistemas Finitos a partir de Dinámica Molecular de Born–Oppenheimer | J.M. Vásquez Pérez | 2011 |
| Time-Dependent Auxiliary Density Perturbation Theory: Method, Implementation and Applications | J. Carmona Espíndola | 2012 |
| Asymptotic Expansion of Molecular Integrals in Self-Consistent Auxiliary Density Functional Methods | A. Álvarez Ibarra | 2013 |
| Low-Order Scaling Methods for Auxiliary Density Functional Theory | D. Mejía Rodríguez | 2015 |

*(Full topic-organized publication lists — covering Basis Sets, Constrained DFT, Cyclic Cluster methods, Functionals, Geometry Optimization, Molecular Graphs, Molecular Integrals, Molecular Properties, Numerical Integration, Parallelization, Population Analysis, QM/MM, Scalar Fields, σ/π Energy Separation, and Transition State Search — are maintained on the official site under the "Publications → Development" section.)*

---

## 8. Summary Assessment

**Strengths:**
- Genuine methodological distinctiveness (ADFT/ADPT) rather than a reimplementation of standard Kohn–Sham DFT, giving real computational advantages for medium-to-large systems (tens to hundreds of atoms, and specialized applications up to thousands).
- Broad, mature property toolkit (NMR, TD-DFT, IR/Raman, X-ray spectra, response properties) built natively around the fitted-density formalism rather than bolted on.
- Strong QM/MM infrastructure with multiple interface options (native, CHARMM, Cuby) suited to biomolecular and condensed-phase modeling.
- Long, well-documented lineage with consistent peer-reviewed methodological validation.

**Limitations / considerations:**
- Some advanced capabilities are restricted to private developer branches rather than the public release, requiring community engagement/registration for full access.
- As a specialist academic code (versus large commercial suites), the user base, packaged tutorials, and third-party interoperability are narrower.
- Documentation is comprehensive but distributed across many topic-specific pages/manuals rather than a single unified reference.

---

*Sources consulted: official deMon2k website (demon-software.com), WIREs Comput. Mol. Sci. review, Molecules tutorial-reviews (2015, 2019), RSC book chapter (2021), Eur. Phys. J. Spec. Top. (2023), and primary ADFT/ADPT methodology papers cited above.*


---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the deMon2k 	DFT-focused quantum chemistry code using Gaussian-type auxiliary density functional methods, efficient for large molecular systems. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
