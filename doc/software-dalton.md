# DALTON: An Exhaustive Review of the Open-Source Quantum Chemistry Package

## 1. Overview

**DALTON** is a general-purpose *ab initio* electronic-structure program suite developed collaboratively by an international consortium of academic groups (primarily in Scandinavia, but with contributors worldwide). It is one of the longest-running actively maintained quantum chemistry codes, with development tracing back to the early-to-mid 1980s, and it has historically distinguished itself less as a "black-box energy/geometry" code and more as **the reference implementation for molecular response theory** — i.e., the systematic, wavefunction- and DFT-based computation of frequency-dependent electric, magnetic, and mixed electric-magnetic properties.

The Dalton **suite** actually comprises two separate executables sharing common infrastructure and philosophy:

- **Dalton** (often called "the Dalton program" or historically "DALTON2011/2013/2016/2018/2020") — the flagship code, strongest at molecular properties and response theory.
- **LSDalton** — a linear-scaling Hartree–Fock and DFT code designed for very large molecular systems, sharing integral and infrastructure code with Dalton but optimized for massively parallel, near-linear-scaling execution.

Since 2020, both codes can also be driven through the **Dalton Project**, a Python scripting/analysis layer that exposes intermediate quantities (integrals, response vectors, etc.) as NumPy arrays and standardizes common workflows (spectra generation, environment embedding, etc.).

- **License:** GNU Lesser General Public License (LGPL) since August 2017 (previously a free-of-charge but restrictive academic license).
- **Language:** Fortran (with C/C++ and Python components in newer infrastructure).
- **Platform:** Unix/Linux (MPI + OpenMP hybrid parallelism).
- **Distribution:** Source code hosted at `gitlab.com/dalton/dalton`; documentation and releases at `daltonprogram.org`.
- **User base:** Historically cited as having over 2,000 issued licenses, reflecting a large worldwide academic user community, particularly in the NMR/EPR spectroscopy, nonlinear optics, and electronic-circular-dichroism communities.

---

## 2. Historical Development

| Period | Milestone |
|---|---|
| Early 1980s | Origins in the SIRIUS MCSCF program (Jensen & Ågren, Uppsala) and the ABACUS response-property code, later merged into a unified framework. |
| 1997 | First public "DALTON" release combining SIRIUS (MCSCF), ABACUS (response properties), and the HERMIT integral code. |
| ~2005 (v2.0) | Introduction of Kohn–Sham DFT and time-dependent DFT (TDDFT) functionality. |
| 2013 (v2013) | Major expansion of coupled-cluster response and relativistic four-component capability sharing infrastructure with **DIRAC** (a related, separately maintained relativistic code co-founded by several Dalton developers). |
| 2013–2014 | Introduction and consolidation of **LSDalton** as the linear-scaling module for large systems. |
| August 2017 | Dalton becomes fully open source under the LGPL license (previously a gratis-but-closed academic license model). |
| 2020 (v2020) | Release of the **Dalton Project** Python platform; continued extension of multiconfigurational short-range DFT (MC-srDFT) and environment (PCM/PE/QM-MM) models. |
| Ongoing | Active development on GitLab; incremental releases continue to extend property modules, environment models, and relativistic/spin–orbit treatments. |

The two principal review-style "master" papers documenting the whole package are the 2014 *WIREs Comput. Mol. Sci.* paper (Aidas *et al.*) and the 2020 *J. Chem. Phys.* Dalton Project paper (Olsen *et al.*) — both listed in the bibliography below — supplemented by numerous method-specific primary publications.

---

## 3. Electronic-Structure Methods

Dalton supports a broad hierarchy of wavefunction- and density-based methods, essentially all of which are also available (to varying extent) for response/property evaluation, not merely for ground-state energies:

### 3.1 Mean-field and correlated wavefunction methods
- **Hartree–Fock (HF)** — restricted (RHF), unrestricted (UHF), and restricted open-shell (ROHF).
- **Kohn–Sham DFT (KS-DFT)** — a wide range of LDA, GGA, meta-GGA, and hybrid/global-hybrid exchange-correlation functionals, plus range-separated (Coulomb-attenuated, CAM-style) functionals.
- **Møller–Plesset perturbation theory (MP2)**, including RI/Cholesky-approximated variants for large systems.
- **Multiconfigurational SCF (MCSCF)**, including CASSCF and RASSCF, via the SIRIUS module.
- **Multi-reference / multiconfigurational short-range DFT (MC-srDFT)** — a hybrid wavefunction-DFT approach combining long-range MCSCF/CASSCF treatment of static correlation with a short-range DFT description of dynamical correlation, a signature Dalton-specific method.
- **Configuration Interaction (CI)**, including large-scale CI via the LUCITA module.
- **Coupled Cluster (CC)** hierarchy: CCS, CC2, CCSD, CCSD(T), CC3, and the general "CC response" framework — historically one of Dalton's flagship strengths, particularly for high-accuracy frequency-dependent properties and excited states.

### 3.2 Relativistic and heavy-element treatments
- Scalar-relativistic approaches, including one- and two-component treatments.
- Spin–orbit coupling in select modules.
- Shared ancestry/infrastructure with **DIRAC**, the four-component fully relativistic sister program, for users needing rigorous relativistic treatments beyond Dalton's native scope.

### 3.3 Linear-scaling and large-system methods (LSDalton)
- Density-fitting / Cholesky-decomposition-based integral evaluation.
- Linear-scaling Hartree–Fock and DFT via localized/sparse-matrix techniques.
- Massively parallel (MPI/OpenMP hybrid) execution, enabling calculations on molecules with hundreds to over a thousand atoms.

### 3.4 Environment / embedding models
- **Polarizable Continuum Model (PCM)** and related dielectric-continuum solvation models.
- **Polarizable Embedding (PE)** model — QM/MM-type embedding with polarizable classical environments, developed extensively within the Dalton/PyFraME ecosystem for biomolecular and solvated systems.
- Traditional discrete/QM-MM coupling schemes.
- **PyFraME** — an automated fragmentation and parameterization package (part of the Dalton Project) for building embedding potentials for complex environments (e.g., proteins, solvent shells) with minimal manual setup.

---

## 4. Molecular Properties and Response Theory — Dalton's Core Strength

Dalton's defining feature, emphasized across every major review of the package, is its systematic and unusually complete implementation of **analytic response theory** at essentially all supported levels of electronic-structure theory (HF, MCSCF, DFT/TDDFT, and CC). This is the origin of its strong reputation in the spectroscopy and nonlinear-optics communities.

### 4.1 Response theory hierarchy
- **Linear response** — first-order frequency-dependent properties (e.g., polarizabilities, NMR shielding tensors, excitation energies as poles of the linear response function).
- **Quadratic response** — second-order nonlinear properties (e.g., first hyperpolarizabilities, two-photon absorption cross-sections).
- **Cubic response** — third-order nonlinear properties (e.g., second hyperpolarizabilities, three-photon processes).
- Response formulated consistently at the **HF, MCSCF/MC-srDFT, KS-DFT/TDDFT, and CC** levels, with the CC response hierarchy (CC2, CCSD, CC3) providing a systematically improvable route to high accuracy for excited states and frequency-dependent properties.

### 4.2 Excited states and TDDFT
- **Time-Dependent DFT (TDDFT)**, formulated within the linear-response framework, for vertical excitation energies, oscillator strengths, and one-/two-photon absorption properties.
- Restricted and unrestricted (open-shell) TDDFT, including the "restricted-unrestricted" generalized approach for open-shell reference states.
- CC-based excited-state models (CC2, CCSD, CC3) as alternatives/benchmarks to TDDFT.
- MC-srDFT-based excited-state treatment for systems with strong static correlation where single-reference TDDFT is unreliable.

### 4.3 Perturbation operators and spectroscopic properties
Dalton implements an unusually large catalogue of **singlet and triplet perturbation operators**, enabling calculation of:

- **Electric properties:** dipole moments, (hyper)polarizabilities, two- and three-photon absorption/emission cross-sections.
- **Magnetic resonance properties:** NMR shielding constants and spin–spin coupling constants, EPR g-tensors and hyperfine coupling constants — computed in a fully **gauge-origin-invariant** manner (via London atomic orbitals / GIAOs), which historically has been a major Dalton strength and adoption driver in the NMR community.
- **Optical activity:** electronic and vibrational circular dichroism (ECD/VCD), optical rotation, Raman optical activity (ROA) — again gauge-origin-invariant.
- **Vibrational properties:** molecular gradients and Hessians for geometry optimization, harmonic vibrational frequencies, and IR/Raman intensities.
- **Nonlinear/multiphoton spectroscopies:** one-, two-, and three-photon absorption processes via the singlet/triplet operator framework combined with quadratic/cubic response.
- **X-ray and core-level spectroscopies** supported in more recent versions via specialized response approaches for core excitations.

### 4.4 Analytic derivatives
- Analytic molecular gradients and (in many cases) Hessians across supported wavefunction models, used for geometry optimization, transition-state search, vibrational analysis, and molecular dynamics.

---

## 5. Software Architecture

The Dalton codebase is modular, with distinct components historically developed somewhat independently and later integrated:

| Module | Function |
|---|---|
| **HERMIT** | Core one- and two-electron integral evaluation, including all integrals of the Breit–Pauli Hamiltonian needed for magnetic/relativistic properties. |
| **ERI** | Vectorized, distribution-oriented integral generator used especially in integral-direct coupled-cluster calculations. |
| **SIRIUS** | (MC)SCF energy optimization engine (direct, second-order/norm-extended optimization methods). |
| **CC module** | Coupled-cluster energy, gradient, and response-property calculations. |
| **LUCITA** | Large-scale CI engine for general CI expansions; also supports MCSCF with parallel configuration-vector and density-matrix evaluation. |
| **DFT module** | Numerical exchange–correlation quadrature, feeding Kohn–Sham matrix elements to SIRIUS. |
| **LSDalton** | Separate executable for linear-scaling HF/DFT on large systems, sharing integral infrastructure. |
| **Dalton Project (Python layer)** | Post-2020 scripting/analysis interface; exposes integrals and results as NumPy arrays; interfaces with PyFraME for embedding workflows. |

Parallelization is via a hybrid MPI + OpenMP model, historically using a master/worker task-distribution design (a design whose scalability at very high process counts has itself been the subject of dedicated performance-engineering studies — see bibliography).

---

## 6. Typical Application Domains

Reflecting its property/response focus, Dalton is used heavily in:

- **NMR and EPR spectroscopy** — prediction of shielding constants, spin–spin couplings, g-tensors, and hyperfine couplings, especially where gauge-origin invariance is essential.
- **Nonlinear optics (NLO)** — computational screening of hyperpolarizabilities for candidate NLO chromophores (e.g., push–pull Schiff-base and metal-complex systems).
- **Chiroptical spectroscopy** — ECD, VCD, ROA, and optical rotation for stereochemical characterization.
- **Multiphoton absorption** — two- and three-photon cross-section prediction for photonics and bioimaging chromophore design.
- **Solvation and biomolecular environments** — via PCM and Polarizable Embedding, e.g., property calculations for chromophores embedded in proteins or solvent.
- **Benchmark and method-development studies** — as a testbed for new response-theory and CC-response methodology given its transparent, well-documented response infrastructure.
- **Large-system electronic structure** (via LSDalton) — HF/DFT energetics on systems of hundreds to 1000+ atoms.

---

## 7. Distinguishing Features vs. Other Packages

- **Breadth and depth of response theory** across essentially the full hierarchy of electronic-structure models (HF → MCSCF/MC-srDFT → DFT/TDDFT → CC), rather than response capability being bolted onto only one or two reference methods, is widely regarded (including by the developers themselves and independent reviewers) as unmatched or nearly unmatched among general-purpose packages.
- **Gauge-origin invariance** built systematically into magnetic-property and optical-activity modules via London orbitals.
- **MC-srDFT**, a hybrid multiconfigurational/DFT method largely unique to the Dalton ecosystem, addressing static correlation with a DFT-level treatment of dynamic correlation.
- **Historical and ongoing kinship with DIRAC**, giving a conceptual and code-sharing bridge to fully relativistic four-component calculations for heavy-element systems.
- **Open architecture for methodology development** — because response equations are implemented generically, Dalton has long served as a proving ground for new response-theory methods before broader adoption elsewhere.
- Comparatively **less emphasis on "black-box" workflow convenience** (geometry optimization ergonomics, broad functional/basis-set catalogs comparable to Gaussian/ORCA/Q-Chem) relative to its property/response depth — this is a trade-off noted implicitly by its niche adoption pattern (strong in spectroscopy-focused groups, less commonly a first choice for routine organic/ground-state screening).

---

## 8. Access and Practical Notes

- Source code: `https://gitlab.com/dalton/dalton`
- Documentation/manuals/tutorials: `https://daltonprogram.org`
- License: LGPL (open source since August 2017); citation of the appropriate reference paper(s) is requested/required for published work using the program.
- Deployed as a standard module on numerous national/university HPC systems (e.g., HPC2N, UiT, bwHPC), typically built as an MPI+OpenMP hybrid.

---

## 9. Bibliography — Key Publications on DALTON's Theory and Implementation

### 9.1 Primary "master" program-description papers
1. Aidas, K.; Angeli, C.; Bak, K. L.; Bakken, V.; Bast, R.; Boman, L.; Christiansen, O.; Cimiraglia, R.; Coriani, S.; Dahle, P.; Dalskov, E. K.; Ekström, U.; Enevoldsen, T.; Eriksen, J. J.; Ettenhuber, P.; Fernández, B.; Ferrighi, L.; Fliegl, H.; Frediani, L.; Hald, K.; Halkier, A.; Hättig, C.; Heiberg, H.; Helgaker, T.; Hennum, A. C.; Hettema, H.; Hjertenæs, E.; Høst, S.; Høyvik, I.-M.; Iozzi, M. F.; Jansík, B.; Jensen, H. J. Aa.; Jonsson, D.; Jørgensen, P.; Kauczor, J.; Kirpekar, S.; Kjærgaard, T.; Klopper, W.; Knecht, S.; Kobayashi, R.; Koch, H.; Kongsted, J.; Krapp, A.; Kristensen, K.; Ligabue, A.; Lutnæs, O. B.; Melo, J. I.; Mikkelsen, K. V.; Myhre, R. H.; Neiss, C.; Nielsen, C. B.; Norman, P.; Olsen, J.; Olsen, J. M. H.; Osted, A.; Packer, M. J.; Pawlowski, F.; Pedersen, T. B.; Provasi, P. F.; Reine, S.; Rinkevicius, Z.; Ruden, T. A.; Ruud, K.; Rybkin, V. V.; Sałek, P.; Samson, C. C. M.; Sánchez de Merás, A.; Saue, T.; Sauer, S. P. A.; Schimmelpfennig, B.; Sneskov, K.; Steindal, A. H.; Sylvester-Hvid, K. O.; Taylor, P. R.; Teale, A. M.; Tellgren, E. I.; Tew, D. P.; Thorvaldsen, A. J.; Thøgersen, L.; Vahtras, O.; Watson, M. A.; Wilson, D. J. D.; Ziolkowski, M.; Ågren, H. **"The Dalton quantum chemistry program system."** *WIREs Comput. Mol. Sci.* **2014**, *4* (3), 269–284. DOI: 10.1002/wcms.1172.
2. Olsen, J. M. H.; Reine, S.; Vahtras, O.; Kjellgren, E.; Reinholdt, P.; Dundas, K. O. H.; Li, X.; Cukras, J.; Ringholm, M.; Hedegård, E. D.; Di Remigio, R.; List, N. H.; Faber, R.; Tenorio, B. N. C.; Bast, R.; Pedersen, T. B.; Rinkevicius, Z.; Sauer, S. P. A.; Mikkelsen, K. V.; Kongsted, J.; Coriani, S.; Ruud, K.; Helgaker, T.; Jensen, H. J. Aa.; Norman, P. **"Dalton Project: A Python platform for molecular- and electronic-structure simulations of complex systems."** *J. Chem. Phys.* **2020**, *152* (21), 214115. DOI: 10.1063/1.5144298.

### 9.2 Foundational response-theory and MCSCF methodology
3. Jensen, H. J. Aa.; Ågren, H. **"Documentation of SIRIUS: A General Purpose Direct Second Order MCSCF Program."** Technical Notes TN783, TN784, TN785, Dept. of Quantum Chemistry, University of Uppsala, 1986.
4. Jensen, H. J. Aa.; Ågren, H. **"MCSCF optimization using the direct, restricted step, second-order norm-extended optimization method."** *(Foundational Dalton MCSCF optimization method.)*

### 9.3 DFT / TDDFT methodology specific to Dalton
5. Teale, A. M.; De Proft, F.; Tozer, D. J. **"Orbital energies and negative electron affinities from density functional theory: insight from the integer discontinuity."** *J. Chem. Phys.* **2008**, *129*, 044110.
6. Rinkevicius, Z.; Tunell, I.; Sałek, P.; Vahtras, O.; Ågren, H. **"Restricted density functional theory of linear time-dependent properties in open-shell molecules."** *J. Chem. Phys.* **2003**, *119*, 34–46.
7. Oprea, C. I.; Telyatnyk, L.; Rinkevicius, Z.; Vahtras, O.; Ågren, H. **"Time-dependent density functional theory with the generalized restricted-unrestricted approach."** *J. Chem. Phys.* **2006**, *124*, 174103.

### 9.4 Performance / parallelization studies
8. **"Scaling Dalton, a molecular electronic structure program."** Performance and parallel-scalability study addressing the master/worker design bottleneck of Dalton at high process counts (DiVA-portal technical report).
9. Performance characterization and optimization study of the Dalton program system for improved scalability and parallel efficiency (referenced in connection with LDRD "Quantum Vibrations in Molecules" final report).

### 9.5 Related/companion program (shared lineage)
10. Saue, T. *et al.* — **DIRAC**, the four-component relativistic quantum chemistry program sharing historical and methodological lineage with Dalton (relevant for users requiring fully relativistic treatments beyond Dalton's native scope).

> **Citation note:** Users publishing results obtained with the Dalton or LSDalton programs are expected to cite the appropriate master reference(s) above (items 1 and/or 2, depending on version/module used) together with any method-specific papers relevant to the particular property or model employed, per the citation policy maintained at daltonprogram.org.

---

## 10. Summary

DALTON occupies a distinctive niche among general-purpose quantum chemistry packages: rather than competing primarily on breadth of routine ground-state functionality, it has built a 40-year track record as the most systematically complete platform for **analytic response theory** — linear, quadratic, and cubic — across the HF, MCSCF/MC-srDFT, DFT/TDDFT, and coupled-cluster hierarchies. Its gauge-origin-invariant treatment of magnetic and chiroptical properties, its unique MC-srDFT method, its shared lineage with the relativistic DIRAC code, and (via LSDalton) its linear-scaling capability for large systems collectively make it a primary tool of choice for the NMR/EPR, nonlinear-optics, and chiroptical-spectroscopy communities, now fully open source under the LGPL and further extended by the Python-based Dalton Project platform.


---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the DALTON 	Open-source quantum chemistry package with a strong focus on molecular properties and response theory, including DFT/TDDFT. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
