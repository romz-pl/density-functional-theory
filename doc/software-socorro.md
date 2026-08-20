# Socorro: An Object-Oriented Electronic Structure Code for Density-Functional Theory

## 1. Overview

**Socorro** is an open-source, object-oriented software package for performing self-consistent, plane-wave, pseudopotential-based **density-functional theory (DFT)** calculations of the electronic structure of periodically repeated (and, with appropriate cell construction, isolated/molecular) systems. It was developed primarily at **Sandia National Laboratories** (Albuquerque, NM), as a collaborative effort also involving **Vanderbilt University** and (per some records) **Wake Forest University**. It is described as a collaborative effort with Wake Forest and Vanderbilt Universities, providing object-oriented software for performing electronic-structure calculations based on density-functional theory, utilizing libraries such as MPI, BLAS, and LAPACK and developed under the GNU General Public License.

Socorro's core purpose is stated directly in its official software description: Socorro can be used to compute the ground-state electron density for a periodically repeated simulation cell in which the external potential is obtained from norm-conserving pseudopotentials or projector-augmented-wave functions. The ground-state electron density is computed by summation over one-electron orbitals which are obtained using the Kohn-Sham formulation of density-functional theory and which are expanded in a plane-wave basis. Various quantities can be computed from the ground-state solution, including atomic forces which can be used to optimize atom positions in the simulation cell and to perform molecular dynamics simulations.

---

## 2. History and Development Context

- Socorro's development traces back to Sandia's broader, long-running program of electronic-structure code development, sitting alongside sibling Sandia codes such as SeqQuest (a localized-orbital LCAO code) and the earlier QUEST massively-parallel electronic-structure effort. Sandia scientists developed QUEST (QUantum Electronic STructure), an electronic structure program designed for massively parallel supercomputers, useful for understanding systems at a quantum level and important for applications such as catalysis, drug design, microelectronics, new chemicals, and new alloys. Socorro emerged from this lineage as Sandia's modern, object-oriented plane-wave DFT platform.
- The code has long been indexed by Sandia's own DFT/CCR pages as one of its principal in-house electronic-structure tools, alongside SeqQuest, for materials modeling work. This resource area is intended to provide information useful for people using and/or developing density-functional-theory-based tools for electronic structure calculations, with a focus on usage and development of DFT methods within Sandia.
- A dedicated FAQ and "About/Intro" page existed at Sandia's `dft.sandia.gov/Socorro` domain (currently returning access errors, suggesting the original hosting has been retired or restricted).
- The Socorro developer/author list registered with the U.S. Department of Energy's Office of Scientific and Technical Information (OSTI) includes: Foiles, Stephen M.; Leung, Kevin; Lippert, Ross A.; Wright, Alan F.; Modine, Normand A.; Plimpton, Steven J.; Wills, Ann Elisabet (all Sandia); Hatcher, Ryan; Tackett, Alan R. (Vanderbilt-affiliated); and Henkelman, Graeme. This reflects the multi-institution, multi-author character of the project consistent with the Wikipedia summary of a Sandia/Vanderbilt/Wake Forest collaboration.
- The DOE/HydroGEN capability profile lists **Normand Modine** and **Alan Wright** (both of Sandia's Center for Computing Research) as the current capability experts/maintainers for the code, reflecting continued Sandia stewardship of the project into the 2010s–2020s.
- Socorro has been used to underpin methodological work as recently as 2023, appearing in a Physical Review E journal-article listing on OSTI associated with the software record, indicating ongoing scientific output tied to the code.

---

## 3. Licensing and Distribution

- Socorro is an open-source code available under the terms of the GNU Public License (GPL). It is available for download from Vanderbilt University.
- This matches the general Sandia open-source software catalog description: it is developed under the GNU General Public License, consistent with Sandia's broader practice of releasing computational science tools (ParaView, Pyomo, etc.) as open source.
- The primary Sandia-hosted documentation portal (`dft.sandia.gov/Socorro/`) is referenced in the scientific literature as the authoritative source: SOCORRO is developed at Sandia National Laboratories and available from http://dft.sandia.gov/Socorro/. As of this review, that domain's specific Socorro subpages return access-rejection errors, so current users are advised to obtain the code via the Vanderbilt distribution point or by contacting Sandia's CCR/DFT program directly.
- Socorro does **not** currently appear to maintain an actively-updated public GitHub organization repository the way many newer Sandia codes do (e.g., Albany, SEACAS); it predates Sandia's large-scale migration of software to GitHub and appears to remain distributed primarily via direct download/tarball rather than a continuously mirrored public git history.

---

## 4. Physical and Numerical Methodology

### 4.1 Electronic-structure framework
Socorro solves the **Kohn–Sham equations** of DFT self-consistently for the ground-state electron density of a periodic simulation cell. The ground-state electron density is computed by summation over one-electron orbitals obtained using the Kohn-Sham formulation of density-functional theory, expanded in a plane-wave basis.

### 4.2 Treatment of core electrons / ion cores
Socorro supports two complementary approaches for representing the ionic potential and core electrons:
- **Norm-conserving pseudopotentials (NCPs)**
- **Projector-Augmented-Wave (PAW) potentials**

Socorro can be used to perform spin-polarized and non-spin-polarized DFT calculations in a plane-wave basis using norm-conserving pseudopotentials (NCPs) or projector-augmented-wave (PAW) potentials to model the ions and core electrons, and density-dependent functionals from the LibXC library to describe exchange and correlation effects among valence electrons.

### 4.3 Exchange-correlation functionals
Socorro interfaces with the widely used **LibXC** functional library for semi-local (LDA/GGA-type) exchange-correlation treatments, and additionally supports orbital-dependent **hybrid functionals**: In addition, NCP-based calculations can be performed using LibXC hybrid functionals, which combine density- and orbital-dependent (exact) exchange with density-dependent correlation.

A distinguishing methodological contribution of the Socorro project is its scalable treatment of **exact (Fock) exchange** in a plane-wave/pseudopotential hybrid-functional context — historically one of the hardest operations to parallelize efficiently in plane-wave codes: Novel algorithms are used to construct and apply the exact-exchange operator, allowing the use of significantly more processing cores than are currently allowed with commercial and academic electronic structure codes.

### 4.4 Dynamics and structural capabilities
Beyond ground-state total-energy and force evaluation, Socorro implements:
1. **Structural relaxation** — ionic-coordinate relaxation plus cell-shape/size optimization.
2. **Molecular dynamics** — both:
   - **Born–Oppenheimer MD** (standard adiabatic ionic dynamics on the DFT ground-state surface), and
   - **Ehrenfest dynamics** — i.e., time-dependent DFT (TDDFT) coupled to ion dynamics, allowing non-adiabatic electron-ion coupled trajectories.
3. **Transition-state / saddle-point searches** for migrating point defects.

Socorro has state-of-the-art capabilities to: (1) relax the ionic coordinates in a simulation cell and optimize the cell size and shape, (2) run Born-Oppenheimer (molecular-dynamics) or Ehrenfest (time-dependent DFT plus ion dynamics) trajectories, and (3) identify transition states of migrating point defects.

Consistent with this, an independent educational/mapping project at Murray State/Kentucky State University characterized Socorro's TDDFT capability directly: Socorro is a free-ware scientific code that implements Time-Dependent Density Functional Theory (TDDFT), a quantum-mechanical method used to determine the time-dependent properties of materials at the level of atoms and electrons.

### 4.5 Usability design
The code was deliberately designed for a low barrier to entry for users already familiar with commercial/academic plane-wave codes: The input formats for lattice vectors and ion positions are similar to those used in commercial and academic codes and the keywords used to invoke capabilities are easily recognized, enabling short learning curves for experienced users of these codes.

---

## 5. Software Architecture

### 5.1 Object-oriented design
Socorro is explicitly characterized in Sandia's own software catalog and on Wikipedia's Sandia software list as **object-oriented** software for DFT-based electronic-structure calculation — a somewhat distinctive design choice among plane-wave DFT codes of its era, many of which (VASP, ABINIT, Quantum ESPRESSO's PWscf) originated as more traditionally structured Fortran 77/90 codes.

### 5.2 Implementation language
Despite the object-oriented design philosophy, the underlying implementation is in Fortran: Socorro is programmed in FORTRAN 90/95, a language widely used in scientific computing. The Socorro code is complex, and makes use of many advanced programming techniques. The complexity of the theory implemented and the FORTRAN implementation in Socorro make modifying and extending the code extremely difficult. This reflects a common pattern in scientific HPC software: object-oriented design patterns implemented within Fortran 90/95's derived-type and module system, rather than a code written in C++/Python.

Because of this complexity, an NSF EPSCoR-funded undergraduate research project (Kentucky State University / University of Kentucky) built a wiki-based "map" or encyclopedia of the Socorro source tree to aid comprehension: A wiki is a flexible and extensible website that allows easy creation of interlinked web pages, and the map is documented in an online wiki. This independently corroborates both the code's object-oriented complexity and its status as a serious, actively studied research code base.

### 5.3 Parallelism and external libraries
- **MPI** — distributed-memory parallelism across compute nodes, essential for the code's scalability claims.
- **BLAS / LAPACK** — dense linear algebra kernels. It utilizes libraries such as MPI, BLAS, and LAPACK.
- **OpenMP / Intel threaded libraries** — for node-level (shared-memory) threading: Socorro is compatible with existing and near-term capability and capacity computing resources at Sandia, including an emerging capability to utilize OpenMP hardware threads and Intel threaded libraries.
- **LibXC** — external exchange-correlation functional library (see §4.3).
- The code's scaling performance is presented as a headline feature: The code was designed to scale efficiently with the number of cores and exceeds the performance of popular commercial codes in capability class calculations, especially for hybrid-functional calculations.

### 5.4 Position among peer plane-wave DFT codes
Socorro is recognized in the broader HPC/DFT community as one of a canonical set of plane-wave pseudopotential DFT/AIMD codes, appearing alongside VASP, CASTEP, CPMD, ABINIT, PWSCF (Quantum ESPRESSO), DACAPO, DFT++, PARATEC, CP2K, SPHInX, and Qbox in comparative surveys of scalable plane-wave DFT-MD software: VASP, CASTEP, CPMD, ABINIT, PWSCF, DACAPO, SOCORRO, DFT++, PARATEC, DOD-PW, CP2K, SPHINX, QBOX, PEtot.

---

## 6. Applications and Scientific Use Cases

Socorro has been applied to a range of materials-science and defect-physics problems, particularly at Sandia's Center for Integrated Nanotechnologies (CINT) and in point-defect / semiconductor physics:

- **Semiconductor point-defect physics** — e.g., silicon self-interstitial migration:  Thermal and carrier-induced migration processes of a neutral silicon interstitial in bulk silicon have been obtained using Socorro, involving thermal transitions between a ground state and a metastable state through a transition state, and carrier-induced transitions caused by alternating capture of holes and electrons.
- **Exchange-correlation functional development** — Socorro served as the implementation platform for a Sandia-developed surface-specific XC functional: a density-functional-theory (DFT) exchange-correlation functional designed to enable an accurate treatment of systems with electronic surfaces, validated for two metals (Al, Pt) and one semiconductor (Si) via bulk lattice constants, bulk moduli, and vacancy formation energies.
- **Point-defect bounds analysis in semiconductors/insulators**, and III-V compound semiconductor defect surveys (see publication list, §7).
- **Electrolysis / water-splitting materials research**, as part of DOE's HydroGEN Advanced Water Splitting Materials Consortium capability catalog.
- **Cross-validation / benchmarking role**, e.g., serving as a target output format for pseudopotential-generation and XC-functional testing workflows in supplementary materials of related Sandia-authored papers. The software was modified to use functionals presented in a related paper and to obtain the XC potential from the calculation, using Hamann-type and Troullier-Martins-type pseudopotentials for the elements studied.

---

## 7. Related Publications (Theory, Methods, and Applications)

The following publications are directly associated with Socorro's theoretical foundations, methodology, or its use as the computational engine for defect/materials-physics studies. (Citation details reflect what is available from search of the scientific literature; users should verify full bibliographic details against the original journals.)

### Exchange-correlation functional theory
- **A.E. Mattsson, R. Armiento, P.A. Schultz, and T.R. Mattsson**, *"A functional designed to include surface effects in self-consistent density-functional theory,"* *Physical Review B* **72**, 085108 (2005). — Presents a subsystem-functional (surface + interior interpolation) approach implemented and validated in Socorro. SOCORRO is developed at Sandia National Laboratories and available from http://dft.sandia.gov/Socorro/... A subsystem functional approach is used, in which an interpolation index combines a surface functional with a functional for interior regions, validated for Al, Pt, and Si via bulk properties and vacancy formation energies.

### Point-defect and bounds-analysis theory (Socorro as computational engine)
- **N.A. Modine, A.F. Wright, and S.R. Lee**, *"Bounds on the range of density-functional-theory point defect levels in semiconductors and insulators,"* *Computational Materials Science* **91**, 431 (2014).
- **A.F. Wright and N.A. Modine**, *"Application of the bounds-analysis approach to arsenic and gallium antisite defects in gallium arsenide,"* *Physical Review B* **91**, 014110 (2015).
- **S.R. Lee, A.F. Wright, N.A. Modine, C.C. Battaile, S.M. Foiles, J.C. Thomas, and A. Van der Ven**, *"First-principles survey of the structure, formation energies, and transition levels of As-interstitial defects in InGaAs,"* *Physical Review B* **92**, 045205 (2015).

(These three papers are cited together in Sandia's DOE/HydroGEN capability description of Socorro: Bounds on the range of density-functional-theory point defect levels in semiconductors and insulators, N.A. Modine, A.F. Wright, and S.R. Lee, Computational Materials Science 91, 431 (2014). Application of the bounds-analysis approach to arsenic and gallium antisite defects in gallium arsenide, A.F. Wright and N.A. Modine, Physical Review B 91, 014110 (2015). First-principles survey of the structure, formation energies, and transition levels of As-interstitial defects in InGaAs, S.R. Lee, A.F. Wright, N.A. Modine, C.C. Battaile, S.M. Foiles, J.C. Thomas, and A. Van der Ven, Physical Review B 92, 045205 (2015).)

### Pseudopotential optimization methodology (Socorro-adjacent tooling)
- **C.N. Brock, B.C. Paikoff, M.I. Md Sallih, A.R. Tackett, and D.G. Walker**, *"Force-based optimization of pseudopotentials for non-equilibrium configurations,"* *Computer Physics Communications* (year per journal record) — listed alongside the defect-physics papers as a Socorro-associated Computer Physics Communications publication on force-based pseudopotential optimization. (Note: A.R. Tackett is one of the OSTI-registered Socorro co-developers.)

### Software/methodology record
- **OSTI Software Record**: *"Socorro Electronic Structure Software"*, listing authors S.M. Foiles, K. Leung, R.A. Lippert, A.F. Wright, N.A. Modine, S.J. Plimpton, A.E. Wills, R. Hatcher, A.R. Tackett, and G. Henkelman, with an associated 2023 journal article in *Physical Review E* (OSTI ID: 2283284), indicating continued theoretical/methodological output tied to the code. Journal Article, 2023, Physical Review E, OSTI ID:2283284.

### Contextual/comparative literature (Socorro cited as a peer code)
- Surveys of scalable plane-wave pseudopotential DFT-MD codes that list Socorro among comparable production codes (VASP, ABINIT, Quantum ESPRESSO/PWSCF, CPMD, CASTEP, CP2K, DACAPO, PARATEC, Qbox, DFT++). VASP, CASTEP, CPMD, ABINIT, PWSCF, DACAPO, SOCORRO, DFT++, PARATEC, DOD-PW, CP2K, SPHINX, QBOX, PEtot.

*Note on completeness:* Because Socorro's primary documentation domain (`dft.sandia.gov`) is currently inaccessible and the code lacks a continuously indexed public GitHub bibliography/CITATION file, this publication list reflects what is discoverable through indirect literature and DOE/OSTI records rather than an authoritative, code-maintained reference list. Users seeking a definitive, exhaustive bibliography should contact Sandia's Center for Computing Research (CCR) DFT program directly.

---

## 8. Comparison to Sibling Sandia Codes

| Feature | **Socorro** | **SeqQuest** |
|---|---|---|
| Basis set | Plane waves | Localized Gaussian (LCAO) |
| Core treatment | Norm-conserving PP / PAW | Norm-conserving PP |
| Boundary conditions | Periodic cells | Periodic slabs/solids or finite molecules |
| Scaling character | Plane-wave FFT-based, GPL, MPI/OpenMP | Near-linear scaling LCAO, efficient on modest hardware |
| XC treatment | LibXC (LDA/GGA + hybrids w/ exact exchange) | LDA/GGA, with/without spin polarization |
| License | GPL, publicly downloadable | Distributed by request/Sandia Quest program |

SeqQuest is a general-purpose electronic structure code to compute energies and forces for periodic surfaces (slabs) or solids, or finite molecules, using norm-conserving pseudopotentials and high-quality contracted-Gaussian basis sets in an LCAO approach, actively under development within Sandia's Multiscale Science Department. This positions Socorro and SeqQuest as complementary Sandia DFT platforms — Socorro for plane-wave/periodic-cell problems where PAW accuracy or hybrid-functional treatments are desired, and SeqQuest for very large-scale, LCAO-based simulations on modest hardware.

---

## 9. Summary

Socorro is a mature, GPL-licensed, Fortran-based but object-oriented plane-wave DFT code developed principally at Sandia National Laboratories (with Vanderbilt and Wake Forest collaboration), supporting norm-conserving pseudopotentials and PAW, LibXC semi-local and hybrid exchange-correlation functionals with a scalable exact-exchange implementation, structural/cell relaxation, Born–Oppenheimer and Ehrenfest (TDDFT) molecular dynamics, and transition-state searches for point defects. It has been used extensively in Sandia's own defect-physics and materials-science research (particularly semiconductor point defects and surface-sensitive functional development) and is recognized in the broader computational-materials community as a peer of major production plane-wave DFT codes such as VASP, ABINIT, and Quantum ESPRESSO, while remaining a comparatively narrowly-distributed, lab-centric code relative to those larger open community projects.

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the Socorro 	Open-source, object-oriented electronic structure computer code developed primarily at Sandia National Laboratories for performing self-consistent density-functional theory (DFT) calculations. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
