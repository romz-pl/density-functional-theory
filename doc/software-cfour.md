# CFOUR: An Exhaustive Review

**CFOUR** (**C**oupled-Cluster techniques **For** **O**ur **U**nderstanding of **R**eality) is a program system for high-level *ab initio* quantum-chemical calculations, developed principally by the research groups of John F. Stanton (University of Florida) and Jürgen Gauss (Universität Mainz), with numerous contributing authors worldwide. Its defining strength — and the reason for its name — is an exceptionally deep and often unique implementation of coupled-cluster (CC) theory and its associated analytic derivative and property machinery. CFOUR traces its lineage back to 1989 and has been continuously developed since.

---

## 1. Identity and Distribution

- **Full name / acronym:** CFOUR — Coupled-Cluster techniques for Computational Chemistry.
- **Lead developers:** J.F. Stanton, J. Gauss, L. Cheng, M.E. Harding, F. Lipparini, D.A. Matthews, S. Stopkowicz, P.G. Szalay, with dozens of contributing authors (A.A. Auer, R.J. Bartlett, Y.J. Bomble, O. Christiansen, T.-C. Jagau, C. Puzzarini, K. Ruud, J.D. Watts, and many others).
- **Origins:** Development began in 1989; the first public production release appeared in June 2010, with a major update (v2.1) in July 2019.
- **License model:** Not commercial software. It is distributed free of charge for non-commercial/academic use under a signed license agreement; it is *not* open source in the sense of an unrestricted public repository, but source access is granted to licensees.
- **Underlying integral packages:** CFOUR incorporates and builds upon several external integral/derivative packages — MOLECULE (J. Almlöf and P.R. Taylor), PROPS (P.R. Taylor), ABACUS (T. Helgaker, H.J.Aa. Jensen, P. Jørgensen, J. Olsen) — plus effective-core-potential (ECP) routines by A.V. Mitin and C. van Wüllen.
- **External method interfaces:** CFOUR interfaces with the MRCC program (M. Kállay) for very high-order CI/CC methods, and with GIMIC, NEWTON-X, DIRAC, and QCUMBRE for specialized property/dynamics calculations.

---

## 2. Core Philosophy: A Wavefunction-Based, Coupled-Cluster-Centric Code

CFOUR is best understood as sitting at the opposite pole from general-purpose DFT-heavy packages (e.g., Gaussian, ORCA, NWChem). Its entire architecture is organized around systematically improvable, single-reference wavefunction methods — Hartree–Fock, Møller–Plesset perturbation theory, and, above all, coupled-cluster theory — together with a very large catalogue of analytic derivatives for each. The design goal is chemical accuracy and rigorous, well-defined error control rather than breadth of empirical functionals or applicability to very large systems.

This specialization is also why CFOUR is frequently paired with other programs in practice: several published studies use a DFT package (e.g., Gaussian) for geometry pre-optimization or broad functional screening, and then switch to CFOUR specifically to obtain high-accuracy CCSD(T)/EOM-CC energies, gradients, or spectroscopic properties on top of those structures.

---

## 3. Electronic-Structure Methods

### 3.1 Mean-field / reference methods
- Restricted (RHF) and unrestricted (UHF) Hartree–Fock SCF.
- Restricted open-shell HF (ROHF).
- Two-configurational SCF (TCSCF).
- Complete-active-space SCF (CASSCF), including a quadratically/second-order convergent implementation and a Cholesky-decomposition-based variant.
- DIIS-based SCF convergence acceleration.

### 3.2 Perturbation theory (Møller–Plesset / many-body perturbation theory)
- MP2, MP3, MP4 (including the partial SDQ-MP4 variant).
- ROHF-based MBPT/MP formulations for open-shell references.

### 3.3 Configuration-interaction and quadratic CI (present but not the focus)
- CID, CISD (explicitly noted by the developers as "not recommended" relative to CC alternatives).
- QCISD, QCISD(T).

### 3.4 Coupled-Cluster methods — the heart of the package
**a) Truncated cluster operator:**
- CCD, CCSD, CCSDT, and (non-public) CCSDTQ.

**b) Perturbative treatments of higher excitations:**
- CCSD(T), CCSD+T(CCSD), and the non-public CCSDT[Q] / CCSDT(Q).

**c) Approximate/iterative CC hierarchies:**
- CCSDT-1, CCSDT-1b, CCSDT-2, CCSDT-3, CCSDT-4.
- The CCn family: CC2 and CC3.

**d) Linearized coupled-cluster methods:**
- LCCD, LCCSD.

**e) Unitary coupled-cluster:**
- UCC(4), plus a modern quadratic UCCSD implementation.

**f) Brueckner coupled-cluster:**
- B-CCD.

**g) Higher excitations via the MRCC interface (M. Kállay):**
- CI(n), FCI, CCSDTQ, CCSDTQP, CCSDTQPH, general CC(n), CCSDT(Q) variants A/B, CC4, CCSDTQ-1/1a/1b/3, and related iterative/perturbative hierarchies up to arbitrary excitation level.

**h) Open-shell / spin-adaptation variants:**
- UHF- and ROHF-based spin-orbital CC.
- Partially spin-adapted and fully spin-restricted open-shell CC (including CCSD and CCSD with triples).
- Unitary-group-based, fully spin-adapted CC using a combinatoric open-shell ansatz (implemented via an automated non-antisymmetric tensor contraction engine).
- Quasi-restricted HF (QRHF) orbital generation for difficult open-shell references.

**i) Multireference coupled-cluster:**
- State-specific Mukherjee-style multireference CC (Mk-MRCC) at the CCSD (and, with MRCC-interface extensions, CCSDT) level — not part of the standard public release.

**j) Relativistic / spin-orbit CC:**
- Two-component CCSD and CCSD(T) with spin-orbit coupling.
- Spin-free exact two-component (SFX2C-1e) relativistic treatment.

### 3.5 Excited, ionized, and electron-attached states
- Equation-of-motion coupled-cluster (EOM-CC): EOM-EE (singlet/triplet excited states), EOM-IP (ionization), and EOM-EA (electron attachment), closely related to Fock-space multireference CC.
- CC linear-response (CC-LR) formulations equivalent to EOM-CC for many properties.
- Analytic gradients for excited/ionized/electron-attached states, enabling geometry optimization on these surfaces.
- Nonadiabatic coupling within the EOM framework.
- Core-level/X-ray spectroscopy methods, including transition-potential coupled-cluster (TP-CC) approaches for accurate X-ray absorption spectra at cost comparable to EOM-CCSD.

### 3.6 Higher-order and cutting-edge methods
- Full CCSDT analytic second derivatives (Hessians) — a capability essentially unique to CFOUR among general-purpose packages.
- CCSDT(Q) analytic gradient theory (recent development), extending the well-established CCSD(T) analog to a full-CCSDT-based perturbative quadruples correction.
- Composite/automated schemes and basis-set extrapolation protocols (e.g., HEAT-type approaches) built around this CC hierarchy.

---

## 4. DFT Functionality — A Notable Limitation

A careful, accurate accounting must state this clearly: **CFOUR is fundamentally a wavefunction-based *ab initio* program and does not offer a native, general-purpose Kohn–Sham density-functional-theory (DFT) module** in the way that packages such as Gaussian, ORCA, Q-Chem, NWChem, or Turbomole do.

Consistent with this:
- The program's own documentation lists all quantum-chemical methods available via the `CALC_LEVEL`/`CALCULATION` keyword, and this list is composed exclusively of Hartree–Fock, MP-*n*, CI, QCI, and the extensive coupled-cluster hierarchy described above. No exchange-correlation functionals (B3LYP, PBE0, TPSS, M06, etc.) or Kohn–Sham SCF option appear anywhere in this specification.
- The Features and Manual pages of the CFOUR wiki contain no DFT, TD-DFT, or exchange-correlation-functional entries.
- In the broader literature, CFOUR is consistently characterized (including by its own overview paper) as a program for "high-level coupled-cluster theory," in contrast to packages explicitly advertised for DFT breadth (e.g., NWChem, ORCA, Turbomole).
- A common practical pattern in published research is to combine CFOUR with a separate DFT-capable program: geometries and/or exploratory functional benchmarking are performed in a DFT package (frequently Gaussian), while CFOUR is invoked specifically for high-accuracy CCSD(T)/QCISD(T)/EOM-CC energies, properties, or spectroscopic constants on top of those structures — precisely because CFOUR itself does not provide the DFT step.

**Practical implication:** users wanting a genuinely integrated wavefunction-theory + DFT + double-hybrid workflow within a single code should look to programs such as ORCA, Q-Chem, Turbomole, Psi4, or NWChem. CFOUR's role in a mixed workflow is exclusively the high-accuracy correlated wavefunction side (typically CCSD(T) and beyond, or EOM-CC for excited states), not the DFT side. If a task genuinely requires CFOUR to perform DFT, this is not supported by the public program and should not be assumed possible without independent confirmation from the current CFOUR manual or developers.

---

## 5. Molecular Properties and Analytic Derivatives

CFOUR's second major strength is the breadth and rigor of its analytic-derivative infrastructure, extended to essentially the full CC hierarchy:

- **Geometries and energetics:** single-point energies, full geometry optimization (including for excited/ionized/electron-attached states), intrinsic reaction coordinate calculations (non-public).
- **First-order properties:** dipole moments and other first-order response properties; Mulliken population analysis.
- **Vibrational properties:** analytic harmonic vibrational frequencies and IR intensities at essentially all supported levels; Raman intensities (non-public); anharmonic force fields via numerical differentiation of analytic derivatives; vibrationally averaged properties; vibrational perturbation theory (VPT2) and effective Hamiltonians; thermodynamic properties at finite temperature.
- **NMR and magnetic properties:** NMR chemical shifts at MP and CC levels; nuclear spin–rotation constants; magnetizabilities; rotational g-tensors; indirect (J) spin–spin coupling constants; electronic g-tensors and electron spin-rotation tensors (both non-public).
- **Second- and third-order/response properties:** static and frequency-dependent polarizabilities; frequency-dependent hyperpolarizabilities and Verdet constants (largely non-public); linear- and quadratic-response properties (partially non-public).
- **Relativistic and beyond-Born–Oppenheimer corrections:** scalar-relativistic corrections; spin-orbit coupling; the spin-free exact two-component (SFX2C-1e) method; analytic Diagonal Born–Oppenheimer Correction (DBOC) calculations.
- **Symmetry and efficiency infrastructure:** exploitation of point-group symmetry via a direct-product decomposition approach in both energy and gradient evaluation; frozen-core treatments (including efficient frozen-core gradients); partial-AO algorithms; Cholesky decomposition of two-electron integrals for energies, gradients, and NMR shielding calculations, reducing storage and computational cost.
- **Basis-set technology:** support for basis functions up to i-type angular momentum; effective core potentials; automated composite/additivity schemes and basis-set extrapolation; basis-set superposition error (BSSE) corrections; use of localized orbitals.
- **External-field and finite-difference tools:** finite-field and finite-difference calculation modes; calculations with external electric fields or external one-electron perturbations; support for non-standard isotopologues.
- **Structure refinement:** least-squares fitting of equilibrium (r_e) structures to experimental rotational constants.

---

## 6. Parallelization and Computational Infrastructure

- A dedicated parallel implementation exists for CCSD and CCSD(T) analytic first and second derivatives (Harding, Metzroth, Gauss, Auer, 2008), as well as parallel CCSDT and Mk-MRCCSDT energy evaluation (Prochnow, Harding, Gauss, 2010).
- Parallel finite-difference schemes for harmonic-frequency and VPT2 (anharmonic) calculations, including supporting scripts.
- Grid-based/distributed calculation support, including parallel grid execution for select workflows.
- A "MINT" integral package option and a partial-AO algorithm framework to control memory/disk footprint for large calculations.

---

## 7. Interfaces to External Software

- **MRCC** (M. Kállay): extends CFOUR's accessible excitation levels to full CI, CCSDTQ, CCSDTQP, CCSDTQPH, and general CC(n) methods well beyond CCSDT(Q).
- **DIRAC:** interfacing for relativistic four-component quantum chemistry contexts.
- **GIMIC:** for magnetically induced current density analysis.
- **NEWTON-X:** for nonadiabatic/excited-state dynamics simulations.
- **QCUMBRE:** specialized interface (not part of the public release).
- **Visualization:** graphical interfaces to JMOL, MOLDEN, and MOLEKEL for orbital/structure visualization.
- CFOUR is also commonly used as a "backend" quantum-chemistry engine by multiscale/workflow frameworks (e.g., ASH) that combine it with other QM, MM, and machine-learning methods for hybrid QM/MM or QM/ML simulations.

---

## 8. Representative Application Domains

Consistent with its analytic-derivative and spectroscopic-property strengths, CFOUR is widely used for:
- High-accuracy thermochemistry and reaction energetics (e.g., composite/HEAT-type protocols built on CCSD(T)/CCSDT(Q)).
- Rovibrational spectroscopy: harmonic and anharmonic force fields, vibrationally averaged structures, and IR/Raman intensity prediction.
- NMR and other magnetic-property prediction at correlated wavefunction levels.
- Electronic spectroscopy of excited, ionized, and electron-attached states via EOM-CC, including vibronic Hamiltonians.
- Core-level (X-ray) spectroscopy via transition-potential CC methods.
- Studies of open-shell radicals and doublet/triplet species using ROHF- and spin-adapted CC methods.
- Benchmark-quality reference data generation for astrochemistry, atmospheric chemistry, and method-development studies (frequently in combination with DFT pre-screening performed in other packages).

---

## 9. Summary Assessment

| Aspect | CFOUR's Position |
|---|---|
| Coupled-cluster breadth/depth | Best-in-class; among the very few codes offering full CCSDT and CCSDT(Q) analytic derivatives |
| MP-*n* / HF | Solid, complete, well-optimized reference implementations |
| DFT / Kohn–Sham | **Not natively available** — no exchange-correlation functionals or KS-SCF module in the public program |
| Excited/ionized states | Strong via EOM-CC and CC-LR, with analytic gradients |
| Molecular properties | Extremely broad (NMR, magnetic, vibrational, relativistic) at correlated levels |
| Multireference CC | Present (Mk-MRCC), though not part of the public release |
| Scalability to very large systems | Modest relative to DFT-centric or local-correlation-focused codes; strength lies in accuracy, not system size |
| Licensing | Free for academic/non-commercial use under signed agreement; not open-source |

CFOUR should be regarded as a **specialist high-accuracy coupled-cluster and molecular-property engine**, not a general-purpose DFT platform. Users requiring DFT/TD-DFT functionality alongside CC-level benchmarking typically pair CFOUR with a DFT-capable package.

---

## 10. Key Publications on CFOUR's Theory and Implementation

### 10.1 Primary program-overview papers
1. D.A. Matthews, L. Cheng, M.E. Harding, F. Lipparini, S. Stopkowicz, T.-C. Jagau, P.G. Szalay, J. Gauss, J.F. Stanton, "Coupled-cluster techniques for computational chemistry: The CFOUR program package," *J. Chem. Phys.* **152**, 214108 (2020). https://doi.org/10.1063/5.0004837
2. M.E. Harding, T. Metzroth, J. Gauss, A.A. Auer, "Parallel Calculation of CCSD and CCSD(T) Analytic First and Second Derivatives," *J. Chem. Theory Comput.* **4**, 64–74 (2008). https://doi.org/10.1021/ct700152c

### 10.2 Coupled-cluster theory — foundational and CCSDT/CCSDT(Q) methodology
3. J. Čížek, "On the correlation problem in atomic and molecular systems...," *J. Chem. Phys.* **45**, 4256–4266 (1966). https://doi.org/10.1063/1.1727484
4. G.D. Purvis III, R.J. Bartlett, "A full coupled-cluster singles and doubles model: The inclusion of disconnected triples," *J. Chem. Phys.* **76**, 1910–1918 (1982). https://doi.org/10.1063/1.443164
5. K. Raghavachari, G.W. Trucks, J.A. Pople, M. Head-Gordon, "A fifth-order perturbation comparison of electron correlation theories," *Chem. Phys. Lett.* **157**, 479–483 (1989). https://doi.org/10.1016/S0009-2614(89)87395-6
6. J. Noga, R.J. Bartlett, "The full CCSDT model for molecular electronic structure," *J. Chem. Phys.* **86**, 7041–7050 (1987); Erratum **89**, 3401 (1988). https://doi.org/10.1063/1.452353
7. J.F. Stanton, J. Gauss, J.D. Watts, R.J. Bartlett, "A direct product decomposition approach for symmetry exploitation in many-body methods. I. Energy calculations," *J. Chem. Phys.* **94**, 4334–4345 (1991). https://doi.org/10.1063/1.460620
8. Y.J. Bomble, J.F. Stanton, M. Kállay, J. Gauss, "Coupled-cluster methods including noniterative corrections for quadruple excitations," *J. Chem. Phys.* **123**, 054101 (2005). https://doi.org/10.1063/1.1950567
9. M. Kállay, J. Gauss, "Approximate treatment of higher excitations in coupled-cluster theory," *J. Chem. Phys.* **123**, 214105 (2005). https://doi.org/10.1063/1.2121589
10. E. Prochnow, M.E. Harding, J. Gauss, "Parallel Calculation of CCSDT and Mk-MRCCSDT Energies," *J. Chem. Theory Comput.* **6**, 2339–2347 (2010). https://doi.org/10.1021/ct1002016
11. D.A. Matthews, J.F. Stanton, "Non-orthogonal spin-adaptation of coupled cluster methods: A new implementation of methods including quadruple excitations," *J. Chem. Phys.* **142**, 064108 (2015). https://doi.org/10.1063/1.4907278
12. J.J. Eriksen, D.A. Matthews, P. Jørgensen, J. Gauss, "The performance of non-iterative coupled cluster quadruples models," *J. Chem. Phys.* **143**, 041101 (2015). https://doi.org/10.1063/1.4927247

### 10.3 Open-shell and spin-adapted coupled-cluster methods
13. P.G. Szalay, J. Gauss, "Spin-restricted open-shell coupled-cluster theory," *J. Chem. Phys.* **107**, 9028–9038 (1997). https://doi.org/10.1063/1.475220
14. M. Heckert, O. Heun, J. Gauss, P.G. Szalay, "Towards a spin-adapted coupled-cluster theory for high-spin open-shell states," *J. Chem. Phys.* **124**, 124105 (2006). https://doi.org/10.1063/1.2179070
15. D. Datta, J. Gauss, "A non-antisymmetric tensor contraction engine for the automated implementation of spin-adapted coupled cluster approaches," *J. Chem. Theory Comput.* **9**, 2639–2653 (2013). https://doi.org/10.1021/ct400216h

### 10.4 Multireference and unitary coupled-cluster theory
16. U.S. Mahapatra, B. Datta, D. Mukherjee, "A size-consistent state-specific multireference coupled cluster theory: Formal developments and molecular applications," *J. Chem. Phys.* **101**, 6171 (1999). (Journal/page as listed by CFOUR bibliography.)
17. F.A. Evangelista, W.D. Allen, H.F. Schaefer III, "Coupling term derivation and general implementation of state-specific multireference coupled cluster theories," *J. Chem. Phys.* **127**, 024102 (2007). https://doi.org/10.1063/1.2743014
18. J. Liu, D.A. Matthews, L. Cheng, "Quadratic Unitary Coupled-Cluster Singles and Doubles Scheme: Efficient Implementation, Benchmark Study, and Formulation of an Extended Version," *J. Chem. Theory Comput.* **18**, 2281–2291 (2022). https://doi.org/10.1021/acs.jctc.1c01210

### 10.5 CASSCF and reference-function methodology
19. F. Lipparini, J. Gauss, "Cost-effective treatment of scalar relativistic effects for multireference systems: A CASSCF implementation based on the spin-free Dirac–Coulomb Hamiltonian," *J. Chem. Theory Comput.* **12**, 4284–4295 (2016). https://doi.org/10.1021/acs.jctc.6b00609
20. T. Nottoli, F. Lipparini, J. Gauss, "Second-Order CASSCF Algorithm with the Cholesky Decomposition of the Two-Electron Integrals," *J. Chem. Theory Comput.* **17**, 6819–6831 (2021). https://doi.org/10.1021/acs.jctc.1c00327

### 10.6 Relativistic and spin-orbit coupled-cluster theory
21. F. Wang, J. Gauss, C. van Wüllen, "Closed-shell coupled-cluster theory with spin-orbit coupling," *J. Chem. Phys.* **129**, 064113 (2008). https://doi.org/10.1063/1.2968136

### 10.7 Review articles on coupled-cluster theory (background/context, cited by CFOUR's own bibliography)
22. R.J. Bartlett, "Coupled-cluster approach to molecular structure and spectra: a step toward predictive quantum chemistry," *J. Phys. Chem.* **93**, 1697–1708 (1989). https://doi.org/10.1021/j100342a008
23. T.D. Crawford, H.F. Schaefer III, "An Introduction to Coupled Cluster Theory for Computational Chemists," *Rev. Comput. Chem.* **14**, 33 (2000). https://doi.org/10.1002/9780470125915.ch2
24. R.J. Bartlett, M. Musiał, "Coupled-cluster theory in quantum chemistry," *Rev. Mod. Phys.* **79**, 291 (2007). https://doi.org/10.1103/RevModPhys.79.291
25. I. Shavitt, R.J. Bartlett, *Many-Body Methods in Chemistry and Physics* (Cambridge University Press, 2009).

---

*This review is based on CFOUR's official website and manual (cfour.uni-mainz.de), the 2020 CFOUR overview paper in the Journal of Chemical Physics, and cross-checking against independent literature discussing CFOUR's scope and typical usage patterns, as of August 2026.*


---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the CFOUR 	Quantum chemistry package specializing in coupled-cluster methods, also including DFT functionality. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
