# GAMESS-UK: An Exhaustive Review

## 1. Overview

GAMESS-UK (General Atomic and Molecular Electronic Structure System – UK version) is an ab initio and density functional theory (DFT) quantum chemistry package for molecular electronic structure calculations. It is one of three independently evolved descendants of the original GAMESS code, alongside GAMESS (US) — maintained at Iowa State University — and Firefly (PC GAMESS). Although all three trace back to a common ancestor, GAMESS-UK has followed its own development trajectory for over four decades and today differs substantially from its US cousin in algorithms, code base, and feature set.

GAMESS-UK is historically associated with the UK's Daresbury Laboratory and the academic and industrial computational chemistry community that grew up around the Collaborative Computational Project in Electronic Structure of Molecules (CCP1).

## 2. History and Origins

- **Common ancestor (pre-1981):** The original GAMESS code was developed at the National Resource for Computational Chemistry (NRCC) in the United States by Michel Dupuis, Dave Spangler, and John Wendoloski, and distributed via the NRCC Software Catalog (1980, Program No. QG01).
- **1981 fork:** Martyn F. Guest obtained the NRCC GAMESS code from Michel Dupuis in 1981 while at Daresbury Laboratory. From this point, the UK and US versions diverged into what became two functionally distinct programs sharing only a name and distant common origin.
- **ATMOL heritage:** Much of the UK code's early distinctive character came not from the imported GAMESS base but from **ATMOL**, an earlier British ab initio program developed by Vic R. Saunders and Martyn Guest at Daresbury/Rutherford in the 1970s. ATMOL contributed data-handling and integral technology but, unlike GAMESS, lacked analytic energy gradients — a capability judged essential for practical geometry optimisation, which motivated migration of ATMOL-derived methodology into the GAMESS-derived UK codebase.
- **CCP1 stewardship:** GAMESS-UK became the flagship code of **CCP1**, the UK Collaborative Computational Project for the electronic structure of molecules, originally funded under the Science and Engineering Research Council (SERC) and subsequently the Engineering and Physical Sciences Research Council (EPSRC). Coordinated development was based at CCLRC/STFC Daresbury Laboratory.
- **International consortium model:** Long-term development was coordinated by an international team including Martyn F. Guest (Daresbury Laboratory), Joop H. van Lenthe (Utrecht University, Netherlands), and John Kendrick (University of Bradford / QMolecular Ltd, UK), with substantial contributions from researchers including Ian J. Bush, Huub J. J. van Dam, Paul Sherwood, Jens M. H. Thomas, Remco W. A. Havenith, and many others across UK and European institutions.
- **Continuous development span:** By the time of its major 2005 review publication, the code had been under active development for "nearly a quarter of a century," and it continued to be updated for some years thereafter (stable release 7.0, dated January 2010, is the last widely cited version number).
- **Licensing model:** GAMESS-UK was distributed as proprietary software rather than open source — free of charge for UK academic use, with demonstration, serial, parallel, site-wide, and joint-software-development licence tiers available for other users, administered originally via Daresbury Laboratory and later via the spin-out company Computing for Science Ltd.

## 3. Technical Architecture

- **Languages:** Primarily Fortran, with C components.
- **Platforms:** Historically supported an unusually wide range of architectures for a quantum chemistry code of its era, including x86, x86-64, PowerPC, MIPS, SPARC, and Alpha, and operating systems including Linux, Mac OS X, AIX, Tru64 UNIX, and Windows — reflecting Daresbury's role in benchmarking chemistry codes across successive generations of national supercomputing hardware (e.g. FPS array processors, Cray systems, IBM SP, HPCx).
- **Parallelism:** Supported both message-passing (MPI) and Global Arrays (GA)-based parallel models, including a parallel virtual-memory distributed-data SCF/DFT implementation and parallel post-Hartree–Fock modules; later work explored GPU acceleration of performance-critical integral and SCF kernels.
- **Integral evaluation:** Electron repulsion integrals over s, p, d, f, and g Cartesian Gaussian orbitals were evaluated using Rys polynomial (Rys quadrature) and rotated-axis integral techniques inherited and extended from the ATMOL/GAMESS lineage.

## 4. Core Scientific Capabilities

### 4.1 Wavefunction methods
- Restricted and unrestricted Hartree–Fock (RHF/UHF), and restricted open-shell HF (ROHF)
- Multi-configuration SCF (MCSCF)
- Møller–Plesset perturbation theory (MP2, MP3)
- Coupled cluster methods (CCSD, CCSD(T))
- Configuration interaction (CI), including large-scale Direct-CI for ground and low-lying excited states, using table-driven selection algorithms
- Valence bond theory via the **TURTLE** module, developed by Joop H. van Lenthe, enabling modern non-orthogonal/valence-bond wavefunctions
- Analytic energy gradients throughout most wavefunction levels, enabling efficient geometry optimisation and transition-state location — the capability whose absence in ATMOL had originally motivated adoption of the GAMESS base
- Analytic and numerical second derivatives (Hessians) for vibrational frequency analysis in a number of cases

### 4.2 Density functional theory (DFT)
- A DFT module that was, by the mid-2000s, described as the most extensively developed and actively enhanced part of the package
- Support for a wide range of exchange-correlation functionals (LDA, GGA, hybrid functionals such as B3LYP-type forms), numerical integration grids (e.g. Lebedev-type angular quadrature, Becke-partitioning schemes), and analytic DFT gradients
- Second derivatives of the DFT energy, supporting DFT vibrational analysis
- Access to the Density Functional Repository maintained in collaboration with Daresbury/CCLRC for exchange-correlation functional definitions

### 4.3 Relativistic and heavy-element treatments
- Relativistic effective core potentials (ECPs)
- Implementation of the two-component scaled **zeroth-order regular approximation (ZORA)**, allowing all-electron-like relativistic treatment across essentially the entire periodic table, including the lanthanides

### 4.4 Excited states and spectroscopy
- Time-dependent DFT and correlated excited-state methods supporting prediction of electronic excitation spectra
- Tools for computation of molecular properties relevant to spectroscopy (e.g. NMR-related and other response properties, depending on version and interfaced modules)

### 4.5 Solvation and environment effects
- Continuum solvation models for treating solvent effects on structure and energetics
- Hybrid quantum mechanics/molecular mechanics (QM/MM) capability, primarily realised through tight integration with the **ChemShell** environment (see Section 5), used extensively for enzyme catalysis and zeolite/heterogeneous catalysis problems where a quantum region is embedded in a large classical environment

### 4.6 Applications emphasised in the literature
The package's own review literature repeatedly highlights three flagship application domains:
1. **Enzyme catalysis** — QM/MM modelling of reaction mechanisms in biomolecular active sites
2. **Zeolite and heterogeneous catalysis** — cluster and embedded-cluster modelling of framework materials
3. **Spectroscopy** — prediction of excited-state and vibrational spectra for closed-shell and open-shell species

## 5. Ecosystem and Associated Software

- **ATMOL** — the predecessor UK code (Saunders & Guest) from which early GAMESS-UK integral and SCF technology was adapted.
- **TURTLE** — the valence-bond module integrated into GAMESS-UK, authored by J. H. van Lenthe, allowing modern valence bond (nonorthogonal) wavefunction calculations within the same framework.
- **ChemShell** — a modular, Tcl/Python-scriptable QM/MM environment originally developed at Daresbury by Paul Sherwood and collaborators, for which GAMESS-UK was one of the principal supported QM engines (alongside codes such as GAMESS(US)/interfaced via CHARMM, Gaussian, Turbomole, MNDO99, and others). ChemShell provided the scripting, data-handling, embedding-potential, and optimisation (via the DL-FIND library) infrastructure that let GAMESS-UK be used for large QM/MM studies of enzymes, zeolites, ionic solids, and surfaces. ChemShell was later substantially redeveloped as an open-source, Python-based package ("Py-ChemShell"), which has continued to support GAMESS-UK as one of its ionic-embedding-capable QM backends.
- **CCP1GUI** and **Molden** — visualisation front-ends commonly used alongside GAMESS-UK for input preparation and output visualisation.
- **DL_POLY / GULP** — classical molecular mechanics and lattice-dynamics codes frequently paired with GAMESS-UK within ChemShell-based QM/MM workflows, particularly for ionic solids and catalytic surfaces.

## 6. Relationship to GAMESS (US) and Firefly

| Aspect | GAMESS-UK | GAMESS (US) |
|---|---|---|
| Origin | Forked 1981 from NRCC GAMESS, heavily shaped by ATMOL | Continued NRCC lineage |
| Institutional home | Daresbury Laboratory / CCP1, UK | Iowa State University, Gordon Research Group |
| Valence bond | TURTLE module (van Lenthe) | Not native (no TURTLE-equivalent) |
| Licensing | Proprietary; free for UK academics, tiered licences elsewhere | Source-available freeware, not OSI open-source |
| Distinctive strengths often cited | DFT module depth, ZORA relativistic treatment, QM/MM via ChemShell, valence bond via TURTLE | Broad ab initio method coverage, effective fragment potential (EFP), long-standing very large US user community |

Firefly (formerly PC GAMESS) is a separate, distantly related fork based on the US GAMESS source rather than the UK branch, optimised historically for x86 desktop performance; it is not otherwise connected to the UK codebase.

## 7. Assessment

**Strengths**
- Long, well-documented development history with a stable international consortium model, giving methodological continuity across ~30 years of active work
- Comparatively early and thorough investment in DFT, including relativistic (ZORA) treatments spanning the full periodic table
- Distinctive access to modern valence bond theory (TURTLE) not commonly found integrated into competitor packages
- Deep integration with ChemShell made it a workhorse for QM/MM studies of enzymes and heterogeneous/zeolite catalysis, a niche in which it was particularly influential
- Demonstrated portability across an unusually broad span of historical hardware architectures, reflecting its roots in a national supercomputing centre

**Limitations / considerations**
- Proprietary licensing (rather than open source) limited community-driven extension relative to some contemporaries and successors
- Primary development activity and major version releases appear concentrated up to around 2010, after which the software's visibility in the primary literature and active maintenance signals decline relative to actively maintained modern alternatives (e.g. NWChem, Q-Chem, ORCA, Psi4)
- Much of GAMESS-UK's continuing relevance today comes indirectly, through its role as a legacy QM backend within the actively maintained, open-source Py-ChemShell QM/MM environment, rather than through standalone development
- As with any long-lived scientific code with dozens of contributors, documentation of the very latest (post-2010) internal algorithmic state is comparatively sparse in the open literature relative to earlier, well-reviewed versions

## 8. Summary

GAMESS-UK represents a historically significant, UK-centred branch of the GAMESS family of ab initio/DFT quantum chemistry codes, distinguished by its ATMOL heritage, a strong and early DFT and relativistic (ZORA) capability, native valence bond theory via TURTLE, and a central role in enabling large-scale QM/MM simulation of enzymatic and heterogeneous catalytic systems through its tight coupling with ChemShell. Though its era of most active independent development appears to have concluded around 2010, its methodology and code remain embedded in ongoing computational chemistry workflows via the modern, open-source ChemShell ecosystem.

---

# Publications Related to GAMESS-UK Theory and Development

The following are the principal peer-reviewed and technical publications describing the theoretical methods, algorithms, and applications of GAMESS-UK and its closely coupled software (ATMOL, ChemShell, TURTLE-related QM/MM work). Full bibliographic detail is given where available from search results; some early technical reports are cited as they appear in the secondary literature.

### Core GAMESS-UK descriptive papers
1. Guest, M. F.; Bush, I. J.; van Dam, H. J. J.; Sherwood, P.; Thomas, J. M. H.; van Lenthe, J. H.; Havenith, R. W. A.; Kendrick, J. **"The GAMESS-UK electronic structure package: algorithms, developments and applications."** *Molecular Physics*, 2005, 103 (6–8), 719–747. DOI: 10.1080/00268970512331340592. (Published in a special issue of *Molecular Physics* in honour of Professor Nicholas C. Handy; also available as arXiv:1506.05421.)
2. Guest, M. F.; van Dam, H. J. J. **"Algorithms, developments and applications in molecular modelling: the GAMESS-UK ab initio code."** *AIP Conference Proceedings*, 1999, 479 (1), 9–18. DOI: 10.1063/1.59480.

### Foundational ATMOL / early algorithmic papers underlying GAMESS-UK
3. Saunders, V. R.; Guest, M. F. **"ATMOL3 Part 9."** RL-76-106, Rutherford Laboratory, 1976.
4. Guest, M. F.; Saunders, V. R. **"On methods for converging open-shell Hartree-Fock wave-functions."** *Molecular Physics*, 1974, 28, 819.
5. Guest, M. F.; Wilson, S. Daresbury Laboratory Preprint DL/SCI/P290T; also in *Supercomputers in Chemistry*, eds. P. Lykos and I. Shavitt, ACS Symposium Series 173, 1981, p. 1.
6. Saunders, V. R.; Guest, M. F. *Computer Physics Communications*, 1982, 26, 389.
7. Guest, M. F. **"Performance of Various Computers in Computational Chemistry"** (survey chapter), in *Supercomputer Simulations in Chemistry*, ed. M. Dupuis, Lecture Notes in Chemistry Vol. 44, Springer Verlag, 1986, p. 98.
8. Guest, M. F.; Harrison, R. J.; van Lenthe, J. H.; van Corler, L. C. H. **"Computational Chemistry on the FPS-X64 Scientific Computers: Experience on single- and multi-processor systems."** *Theoretica Chimica Acta*, 1987, 71, 117.

### Original GAMESS lineage (pre-1981 common ancestor)
9. Dupuis, M.; Spangler, D.; Wendoloski, J. **GAMESS**, NRCC Software Catalog, Vol. 1, Program No. QG01, 1980.

### DFT-specific methodology
10. van Dam, H. J. J. **"2nd Derivatives of the Electronic Energy in Density Functional Theory."** CLRC Daresbury Laboratory Technical Report DL-TR-01-002, 2001.
11. Johnson, B. G.; Fisch, M. J. *Journal of Chemical Physics*, 1994, 100, 7429. (DFT gradient/functional methodology cited in GAMESS-UK's DFT implementation.)
12. Becke, A. D. *Journal of Chemical Physics*, 1988, 88, 2547. (Becke exchange functional, used within GAMESS-UK's DFT module.)
13. Murray, C. W.; Handy, N. C.; Laming, G. J. *Molecular Physics*, 1993, 78, 997. (Numerical integration grid methodology.)
14. Mura, M. E.; Knowles, P. J. *Journal of Chemical Physics*, 1996, 104, 9848. (Quadrature grid methodology used in DFT integration.)
15. Lebedev, V. I.; Laikov, D. N. *Doklady Mathematics*, 1999, 59, 477. (Lebedev angular quadrature, used in DFT numerical integration.)
16. Stratmann, R. E.; Scuseria, G. E.; Frisch, M. J. *Chemical Physics Letters*, 1996, 257, 213. (Grid/quadrature scheme referenced in DFT implementation.)
17. Curtiss, L. A.; Raghavachari, K.; Trucks, G. W.; Pople, J. A. *Journal of Chemical Physics*, 1991, 94, 7221. (Gaussian-n composite method context cited alongside GAMESS-UK's thermochemical capability.)

### Relativistic (ZORA) methodology implemented in GAMESS-UK
18. van Lenthe, E.; Baerends, E. J.; Snijders, J. G. (ZORA formalism — foundational relativistic method later implemented within GAMESS-UK by van Lenthe and coworkers.)
19. Faas, S.; van Lenthe, J. H.; Hennum, A. C.; Snijders, J. G. **"An efficient two-component relativistic method to describe conformational stability of a hexacoordinated ruthenium complex"** and related implementation papers describing the two-component scaled ZORA implementation within GAMESS-UK.

### Valence bond theory (TURTLE module)
20. Dijkstra, F.; van Lenthe, J. H.; Havenith, R. W. A.; et al. Papers describing the TURTLE non-orthogonal valence bond program and its integration with GAMESS-UK (multiple methodological papers by van Lenthe and coworkers on Direct-CI and valence bond wavefunction construction).

### QM/MM methodology (ChemShell, tightly coupled to GAMESS-UK)
21. Sherwood, P.; de Vries, A. H.; Guest, M. F.; Schreckenbach, G.; Catlow, C. R. A.; French, S. A.; Sokol, A. A.; Bromley, S. T.; Thiel, W.; Turner, A. J.; Billeter, S.; Terstegen, F.; Thiel, S.; Kendrick, J.; Rogers, S. C.; Casci, J.; Watson, M.; King, F.; Karlsen, E.; Sjøvoll, M.; Fahmi, A.; Schäfer, A.; Lennartz, C. **"QUASI: A general purpose implementation of the QM/MM approach and its application to problems in catalysis."** *Journal of Molecular Structure (THEOCHEM)*, 2003, 632, 1–28. DOI: 10.1016/S0166-1280(03)00285-9.
22. Sherwood, P.; de Vries, A. H.; Collins, S. J.; Greatbanks, S. P.; Burton, N. A.; Vincent, M. A.; Hillier, I. H. *Faraday Discussions*, 1997, 106, 79. (Early QM/MM embedding methodology involving GAMESS-UK.)
23. Sokol, A. A.; Bromley, S. T.; French, S. A.; Catlow, C. R. A.; Sherwood, P. **"Hybrid QM/MM embedding approach for the treatment of localized surface states in ionic materials."** *International Journal of Quantum Chemistry*, 2004, 99, 695–712. DOI: 10.1002/qua.20032.
24. French, S. A.; Sokol, A. A.; Bromley, S. T.; Catlow, C. R. A.; Sherwood, P. *Topics in Catalysis*, 2003, 24, 161.
25. Metz, S.; Kästner, J.; Sokol, A. A.; Keal, T. W.; Sherwood, P. **"ChemShell—a modular software package for QM/MM simulations."** *WIREs Computational Molecular Science*, 2014, 4, 101–110. DOI: 10.1002/wcms.1163.
26. Kästner, J.; Carr, J. M.; Keal, T. W.; Thiel, W.; Wander, A.; Sherwood, P. **"DL-FIND: an open-source geometry optimizer for atomistic simulations."** *Journal of Physical Chemistry A*, 2009, 113, 11856–11865.
27. Lu, Y.; Farrow, M. R.; Fayon, P.; Logsdail, A. J.; Sokol, A. A.; Catlow, C. R. A.; Sherwood, P.; Keal, T. W. **"Open-Source, Python-Based Redevelopment of the ChemShell Multiscale QM/MM Environment."** *Journal of Chemical Theory and Computation*, 2019, 15 (2), 1317–1328. DOI: 10.1021/acs.jctc.8b01036.
28. Lu, Y.; Sen, K.; Yong, C.; Gunn, D. S. D.; Purton, J. A.; Guan, J.; Desmoutier, A.; Nasir, J. A.; Zhang, X.; Zhu, L.; Hou, Q.; Jackson-Masters, J.; Watts, S.; Hanson, R.; Thomas, H. N.; Jayawardena, O.; Logsdail, A. J.; Woodley, S. M.; Senn, H. M.; Sherwood, P.; Catlow, C. R. A.; Sokol, A. A.; Keal, T. W. **"Multiscale QM/MM modelling of catalytic systems with ChemShell."** *Physical Chemistry Chemical Physics*, 2023, 25 (33), 21816–21835. DOI: 10.1039/D3CP00648D.

### GPU/performance-oriented methodology papers
29. **"Acceleration of the GAMESS-UK Electronic Structure Package on Graphical Processing Units"** — study of GPU-based acceleration of Hartree–Fock/SCF integral evaluation and Fock matrix construction within GAMESS-UK, addressing the irregular computational structure of the electron repulsion integral and SCF algorithms.

### Contextual / comparative reviews
30. Wikipedia contributors. **"GAMESS (UK)."** *Wikipedia, The Free Encyclopedia* (overview article summarising licensing, capability list, and core citation).
31. Wikipedia contributors. **"GAMESS."** *Wikipedia, The Free Encyclopedia* (disambiguation/overview of the GAMESS family: UK, US, and Firefly forks).

---

*Note: Items 18–20 and 26 (TURTLE/ZORA implementation specifics) are referenced as they appear cited in the secondary literature (principally the 2005 Molecular Physics review and related QM/MM literature); readers requiring exact original citations for the ZORA and TURTLE implementation papers are advised to consult the reference list of Guest et al. (2005) directly, as the underlying search sources did not yield complete independent bibliographic records for every such item.*


---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the GAMESS (UK) 	Related but independently developed quantum chemistry package with DFT and ab initio capabilities, historically UK-centered. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
