# ACES: An Exhaustive Technical Review

## 1. Overview

**ACES** (Advanced Concepts in Electronic Structure) is a family of ab initio quantum chemistry program systems developed primarily by the Rodney J. Bartlett research group at the Quantum Theory Project (QTP), University of Florida, Gainesville. The suite specializes in **many-body methods** — many-body perturbation theory (MBPT) and, above all, **coupled-cluster (CC) theory** — for the accurate treatment of electron correlation, and has more recently incorporated **Kohn–Sham density functional theory (KS-DFT)**, including the QTP family of "consistent" exchange-correlation functionals developed in-house.

The ACES lineage comprises three principal generations:

| Generation | Character | Released | Status |
|---|---|---|---|
| **ACES I** | Early diagrammatic many-body code | Late 1970s–1980s, Battelle Memorial Institute | Historical / superseded |
| **ACES II** | Serial, vector-oriented ab initio package | Development began 1990 | Legacy; still used; forked into CFOUR |
| **ACES III** | Completely rewritten, massively parallel implementation | Fall 2008 | Actively distributed |
| **Aces4** | Next-generation successor built on the same parallel infrastructure | Development ongoing (GitHub) | In development |

A parallel fork of ACES II, maintained jointly by groups in Mainz, Austin, and Budapest ("ACES II-MAB"), evolved independently into what is now the separate and widely used **CFOUR** package. ACES III is *not* CFOUR; it is the Bartlett-group parallel rewrite that shares ACES II's ancestry but not its codebase.

All three modern programs (ACES II, ACES III, Aces4) are distributed free of charge for academic/research use.

---

## 2. Historical Development

- **Late 1970s:** Diagrammatic many-body correlation methods first implemented at Battelle Memorial Institute, the intellectual root of the ACES lineage (shared ancestry with COLUMBUS).
- **Late 1989:** John F. Stanton, working in the Bartlett group, wrote interfaces to the SCF and integral packages then in use, seeding what became the ACES program system.
- **1990:** Jürgen Gauss joined the Bartlett group as a postdoc; together with Stanton (and John D. Watts), the core MBPT and CC codes were written during 1990–1991, forming the backbone of **ACES II**.
- The original integral infrastructure was not written from scratch — ACES II incorporated the **MOLECULE** integral package (J. Almlöf), the **VPROP** property package (P. R. Taylor), and the **ABACUS** integral-derivative package (T. Helgaker, P. Jørgensen, J. Olsen, H. J. Aa. Jensen), the latter substantially modified for ACES II.
- **1990s:** ACES II bifurcated into two maintained branches — the Florida (Bartlett-group) version, and the "ACES II-MAB" (Mainz–Austin–Budapest) version, which by 2008 had diverged enough in features and architecture to be renamed **CFOUR**.
- **2000s:** Recognizing that ACES II's architecture could not exploit distributed-memory parallel computers, the Bartlett group (principally Erik Deumens, Victor Lotrich, Ajith Perera, Mark Ponton, and Beverly Sanders, with DoD HPCMP funding) undertook a ground-up rewrite targeting petascale machines.
- **Fall 2008:** **ACES III** released — a complete reimplementation of ACES II's computationally demanding kernels (SCF, MBPT(2), CCSD(T) energies/gradients/Hessians) using a new parallel software architecture, the **Super Instruction Architecture (SIA)**.
- **2009–2014:** Successive ACES III releases (3.0.0 through 3.0.7, culminating in the July 2014 release) added new theoretical capabilities — explicitly correlated CC methods (F12), EOM-CC variants, spin-orbit and ESR-property modules, and KS-DFT support — with incremental performance improvements.
- **Aces4:** A further re-engineered successor built on an evolved SIA runtime, addressing scalability limitations identified in ACES III's original CCSD(T) implementation; under continued development and hosted on GitHub.

---

## 3. Architecture: The Super Instruction Architecture (SIA)

ACES III's defining technical contribution is not a new quantum-chemistry method per se but a **new software architecture for writing massively parallel many-body codes**. The rationale: coupled-cluster algorithms involve extremely complex tensor-contraction expressions over distributed, often out-of-core, multidimensional arrays, and hand-coding MPI-level data movement, load balancing, and synchronization for each new method is prohibitively costly in developer time. SIA addresses this by splitting the problem into two layers:

### 3.1 Super Instruction Assembly Language (SIAL)
- A **domain-specific programming language** ("SIAL", pronounced "sail") in which quantum-chemistry algorithms (e.g., the CCSD amplitude equations) are expressed at a high, tensor/block level, rather than in explicit loop nests over MPI ranks.
- Domain scientists write and maintain SIAL code without having to manage low-level parallel data distribution directly — this is the "high-productivity" aspect emphasized throughout the SIA literature.
- The ACES III distribution contains on the order of **580,000 lines of SIAL code** (roughly 200,000 of which are comments).

### 3.2 Super Instruction Processor (SIP)
- The **runtime/virtual machine** that executes compiled SIAL programs. SIP is best understood as an MPMD (multiple-program, multiple-data) parallel virtual machine, functioning conceptually "like a hardware processor" that operates not on scalars but on **blocks of tensor data** ("super instructions").
- SIP manages: distributed/served arrays, out-of-core storage and checkpointing of large arrays, asynchronous communication, block scheduling, and load balancing across potentially tens of thousands of cores.
- The C/C++/Fortran runtime beneath SIAL comprises roughly **230,000 lines** of code (about 62,000 of which are comments).

### 3.3 Supporting Infrastructure
- **OED/ERD atomic integral package** — the parallel integral-evaluation engine feeding ACES III's SIAL-coded methods.
- Demonstrated scalability to **tens of thousands of cores** on petascale DoD/HPCMP systems (showcased at Supercomputing 2010).
- The architecture was later reused/extended in **Aces4**, which reimplemented CCSD and CCSD(T) in an evolved SIA runtime; benchmarking showed Aces4 scales similarly to ACES III but runs its (T) correction roughly 20–30% faster due to runtime-level improvements (elimination of block wait times via asynchronous communication, reduced per-block metadata overhead).
- SIA/SIAL was also explored as a general-purpose parallel-programming platform beyond quantum chemistry (e.g., applications in environmental/computational-chemistry adjacent domains), reflecting its design as a reusable framework rather than an ACES-specific hack.

---

## 4. Theoretical / Methodological Capabilities

ACES III implements the computationally intensive core of ACES II's method set in parallel, with several capabilities original to ACES III. In broad strokes:

### 4.1 Reference wavefunctions
- Restricted Hartree–Fock (RHF), unrestricted Hartree–Fock (UHF), restricted open-shell Hartree–Fock (ROHF) self-consistent field.
- Quasi-restricted Hartree–Fock (QRHF) and Kohn–Sham (KS) references for open-shell/EOM work.

### 4.2 Many-body perturbation theory
- Second-order MBPT [MBPT(2)] energies, analytic gradients, and analytic Hessians — parallelized from the outset in ACES III.

### 4.3 Coupled-cluster theory (the package's core strength)
- **CCSD** — coupled-cluster singles and doubles, with RHF/UHF/ROHF references, energies, gradients, and (via later modules) Hessians.
- **CCSD(T)** — the "gold standard" perturbative-triples correction, with parallel energy and gradient implementations; a principal performance benchmark for the package (petascale scaling demonstrations).
- Higher-rank iterative CC (CCSDT, CCSDTQ, etc.) and combined CC/MBPT hybrid schemes, with automated-derivation and parallel-implementation approaches explored in the broader ACES/SIA literature.
- **Explicitly correlated CC — CCSD(F12)** and its B-approximation variant, using Slater-type correlation factors, Rys-quadrature evaluation of Slater/Yukawa integrals, and Becke fuzzy-cell numerical quadrature for many-electron integrals — implemented specifically in ACES III to accelerate basis-set convergence.
- **Equation-of-motion coupled-cluster (EOM-CC)** family for excited, ionized, and electron-attached states:
  - EE-EOM-CCSD (excitation energies)
  - IP-EOM-CCSD (ionization potentials) and its explicitly correlated **IP-EOM-CCSD(F12)**
  - EA-EOM-CCSD (electron affinities)
  - **DIP-EOM-CCSD(F12)** and **DEA-EOM-CCSD(F12)** — double-ionization-potential and double-electron-attachment EOM-CC, including explicitly correlated F12 variants, implemented in ACES III.
- Massively parallel **linear-response CC** module for first- and second-order static molecular properties (e.g., static polarizabilities) of both closed- and open-shell species over RHF/UHF/ROHF references — demonstrated on systems such as C20 isomers and biphospholylidene dioxide/disulfide oligomers.
- Massively parallel CC modules for **electron spin resonance (ESR) tensors** — isotropic hyperfine coupling constants (A-tensors), with a planned/realized series extending to g-tensors and D-tensors, targeting large open-shell radicals previously inaccessible to CC-level treatment.
- Multi-reference coupled-cluster extensions (e.g., Mk-MRCC-type methods) were an area of active/ongoing development noted in the SIA literature.

### 4.4 Density functional theory
- **Kohn–Sham DFT (KS-DFT)** support, including local/GGA and hybrid/range-separated functionals.
- The **QTP family of "consistent" exchange-correlation functionals** (e.g., CAM-QTP00, CAM-QTP01, CAM-QTP02, QTP17), developed by the Bartlett group specifically so that KS orbital eigenvalues satisfy **Bartlett's IP (ionization-potential) eigenvalue theorem** — i.e., that occupied KS eigenvalues (including core orbitals) closely approximate the corresponding vertical ionization energies, unlike conventional functionals such as B3LYP or PBE0.
- QTP functionals are parameterized as re-optimizations of established range-separated/global hybrid forms (CAM-B3LYP-type, B3LYP-type) against IP/excitation-energy/atomization-energy training sets, and are benchmarked extensively against EOM-CCSD (IP-EOM-CCSD, EA-EOM-CCSD) as a correlated reference.
- Applications include vertical ionization potential and electron affinity benchmarking (including photovoltaic-candidate molecules), core-excitation and X-ray absorption (NEXAFS) spectral simulation via IP-optimized global hybrids, and fundamental-gap/exciton-binding-energy studies of conjugated organic polymers (trans-polyacetylene, polyacenes) combining QTP-DFT with EOM-CC benchmarks computed on ACES III.
- Note: some DFT/QTP-functional benchmarking work has been carried out using ACES III alongside other packages (e.g., NWChem, ACES II) for the DFT side, reflecting that DFT capability is comparatively newer/less central to ACES III than its CC/MBPT core.

### 4.5 Derivative properties and structure determination
- Analytic gradients (RHF/UHF/MBPT(2)/CCSD(T)) for geometry optimization.
- Analytic Hessians (RHF/UHF/MBPT(2)) for vibrational-frequency calculations; CCSD-level Hessian capability was extended over successive releases.
- Static and dynamic molecular properties (polarizabilities, ESR tensors) via linear-response CC.

---

## 5. Software Characteristics

- **Languages:** SIAL (domain-specific tensor algorithm language) for the science layer; C/C++ and Fortran for the SIP runtime and supporting infrastructure.
- **Parallel model:** MPMD execution via the Super Instruction Processor, targeting distributed-memory clusters and demonstrated at petascale (tens of thousands of cores).
- **Licensing/availability:** Free for academic use; installation is nontrivial (the project's own documentation cautions that "the build process has not been smoothed out" relative to ACES II), and the developers historically offered a paid consulting service for installation/validation on user systems.
- **Funding:** ACES III development was principally supported by the U.S. Department of Defense High Performance Computing Modernization Program (HPCMP), specifically the Common High Performance Software Initiative (CHSSI, project CBD-03) and the User Productivity Enhancement and Technology Transfer (PET) program, with continuing support from the U.S. Army Research Office (ARO). Aces4 development received support from an Air Force Research Laboratory HASI grant and internal ENSCO R&D funding.
- **Relationship to CFOUR:** ACES III is architecturally and organizationally distinct from CFOUR (the Mainz-Austin-Budapest descendant of ACES II), though both trace to the same 1990s ACES II codebase and both remain active, separately maintained coupled-cluster packages today.
- **Successor:** **Aces4**, hosted on GitHub, reimplements core CC kernels (CCSD, CCSD(T)) on an improved SIA runtime with better scaling behavior and reduced per-block overhead relative to ACES III, and is presented as the ongoing evolution of the platform.

---

## 6. Representative Application Domains

- Benchmark-quality thermochemistry and CCSD(T)/CCSDT(Q) energetics on systems ranging from small molecules to oligomeric/polymeric chains.
- Static polarizabilities of closed- and open-shell molecules and oligomers via parallel linear-response CC.
- Electron spin resonance parameter prediction (hyperfine couplings) for large organic and inorganic radicals, previously impractical at the CC level due to cost.
- Ionization potentials, electron affinities, and fundamental band gaps of conjugated polymers (e.g., trans-polyacetylene, polyacenes) via combined EOM-CC/QTP-DFT studies, including thermodynamic-limit extrapolations and exciton-binding-energy estimates.
- Core-level (X-ray) ionization and absorption spectroscopy via IP-optimized DFT functionals benchmarked against CC.
- Reaction-mechanism studies involving heavy/actinide-containing systems (e.g., uranium hexafluoride hydrolysis) drawing on the package's correlated-method and DFT capabilities.
- Photovoltaic-candidate molecule screening via IP/EA benchmarking combining DFT/QTP functionals, EOM-CCSD, and G0W0 corrections.

---

## 7. Summary Assessment

ACES III's principal legacy is less a single new quantum-chemical method than a **reusable, high-productivity parallel-programming framework (SIA/SIAL/SIP)** that let the Bartlett group's decades of many-body/coupled-cluster method development be ported to distributed-memory, petascale hardware without a full method-by-method hand-parallelization effort. This enabled genuinely new science — CC-level ESR tensors, CC-level static polarizabilities, and CCSD(T)-quality thermochemistry — on systems well beyond what serial ACES II (or contemporaneous competitor codes) could handle. Its DFT capability, centered on the in-house **QTP functional family**, is a comparatively newer and more specialized contribution, aimed less at general-purpose DFT usage and more at producing KS-DFT methods whose orbital eigenvalues are rigorously interpretable as ionization potentials — a niche but influential idea within the broader DFT literature. The package remains most authoritative, and most cited, for its coupled-cluster and many-body-perturbation-theory core.

---

## 8. Selected Publications (Theory, Software, and Methodology)

### Foundational / Historical
1. H. P. Kelly, *Phys. Rev. Lett.* **23**, 455 (1969). — Early many-body perturbation theory groundwork underlying the ACES lineage's diagrammatic approach.
2. J. F. Stanton, J. Gauss, J. D. Watts, R. J. Bartlett, "A direct product decomposition approach for symmetry exploitation in many-body methods," and related early ACES II papers, *Int. J. Quantum Chem.* **526**, 879 (1992). — Definitive early account of the ACES II program system.

### ACES III Architecture and Software Design
3. V. Lotrich, N. Flocke, M. Ponton, A. Yau, A. Perera, E. Deumens, R. J. Bartlett, "Parallel Implementation of Electronic Structure Energy, Gradient, and Hessian Calculations," *J. Chem. Phys.* **128**, 194104 (2008). https://doi.org/10.1063/1.2920482 — The founding ACES III paper.
4. E. Deumens, V. F. Lotrich, A. Perera, M. J. Ponton, B. A. Sanders, R. J. Bartlett, "Software design of ACES III with the super instruction architecture," *WIREs Comput. Mol. Sci.* **1**, 895–901 (2011). https://doi.org/10.1002/wcms.77 — Definitive description of SIA/SIAL/SIP.
5. "The Super Instruction Architecture: A Framework for High-Productivity Parallel Implementation of Coupled-Cluster Methods on Petascale Computers," in *Petascale Computing: Algorithms and Applications*, Ch. 8 (Elsevier/CRC, 2011). https://doi.org/10.1016/B978-0-444-53835-2.00008-0
6. SIAL Programmer Guide, ACES III Documentation, Quantum Theory Project, University of Florida. https://aces.qtp.ufl.edu
7. A. Perera, R. J. Bartlett, B. A. Sanders, V. F. Lotrich, J. N. Byrd, "Advanced concepts in electronic structure (ACES) software programs," *J. Chem. Phys.* **152**, 184105 (2020). https://doi.org/10.1063/5.0002581 — Comprehensive modern review of ACES II/III/Aces4.
8. S. E. Masters, B. A. Sanders et al., "Leveraging the Super Instruction Architecture to Develop Massively Parallel Computational Chemistry Applications" (Aces4 scaling study). arXiv:2003.01688.

### Coupled-Cluster and EOM-CC Methodology on ACES III
9. P. Verma, A. Perera, R. J. Bartlett, "Massively parallel implementations of coupled-cluster methods for electron spin resonance spectra. I. Isotropic hyperfine coupling tensors in large radicals," *J. Chem. Phys.* **139**, 174103 (2013).
10. A. Perera, J. A. Morales, "Implementation of a Parallel Linear-Response Coupled-Cluster-Theory Module in ACES III: First Application to the Static Polarizabilities of the C20 Isomers and of the Biphospholylidene Dioxide and Disulfide Oligomers," *Adv. Quantum Chem.*, ScienceDirect (2015). https://doi.org/10.1016/bs.aiq.2015.06.001 (approx.)
11. A. Perera, J. A. Morales, "New massively parallel linear-response coupled-cluster module in ACES III: application to static polarisabilities of closed-shell molecules and oligomers and of open-shell radicals," *Mol. Phys.* **114**(3-4), 547–561 (2016). https://doi.org/10.1080/00268976.2015.1126367
12. "Explicitly-correlated double ionization potentials and double electron attachment equation-of-motion coupled cluster methods," implemented and benchmarked in ACES III (CCSD(F12), DIP/DEA-EOM-CCSD(F12)). *Chem. Phys. Lett.* (2017/2018). https://doi.org/10.1016/j.cplett.2017.11.043 (approx.)
13. E. Prochnow, M. E. Harding, J. Gauss, "Parallel calculation of CCSDT and Mk-MRCCSDT energies," *J. Chem. Theory Comput.* **6**, 2339–2347 (2010). — Related parallel higher-order CC context.
14. J. N. Byrd, R. Montgomery Jr., R. J. Bartlett, "Combined Coupled-Cluster and Many-Body Perturbation Theories: Automated Derivation and Parallel Implementation," *J. Chem. Phys.* **121**(24), 12197–12207 (2004).

### QTP Density Functionals and DFT Applications
15. P. Verma, R. J. Bartlett, "Increasing the applicability of density functional theory. III. Do consistent Kohn–Sham density functional methods exist?" and related QTP-functional papers (CAM-QTP00 introduction).
16. "The QTP family of consistent functionals and potentials in Kohn-Sham density functional theory" (foundational QTP functional development), *J. Chem. Phys.* / related QTP series.
17. Y. C. Park, A. Perera, R. J. Bartlett, "Density functionals for core excitations," *J. Chem. Phys.* **157**, 094107 (2022). https://doi.org/10.1063/5.0111095
18. R. J. Bartlett, "Adventures in DFT by a wavefunction theorist," *J. Chem. Phys.* **151**, 160901 (2019). https://doi.org/10.1063/1.5116338 — Perspective covering the IP-eigenvalue theorem and QTP functional rationale.
19. "Vertical valence ionization potential benchmarks from equation-of-motion coupled cluster theory and QTP functionals," *J. Chem. Phys.* **150**, 074108 (2019). https://doi.org/10.1063/1.5084728
20. H. Kim, A. Perera, R. A. Mendes, R. J. Bartlett, "Benchmarking ionization potentials and electron affinities of potential photovoltaic molecules using DFT/QTP functionals and EOM-CC," *J. Chem. Phys.* **163**, 174703 (2025). https://doi.org/10.1063/5.0293131
21. "Examining fundamental and excitation gaps at the thermodynamic limit: A combined (QTP) DFT and coupled cluster study on trans-polyacetylene and polyacene," *J. Chem. Phys.* **156**, 204308 (2022). https://doi.org/10.1063/5.0089829 (approx.)

### Related CFOUR-Lineage Context (for comparison)
22. D. A. Matthews, L. Cheng, M. E. Harding, F. Lipparini, S. Stopkowicz, T.-C. Jagau, P. G. Szalay, J. Gauss, J. F. Stanton, "Coupled-Cluster Techniques for Computational Chemistry: The CFOUR Program Package," *J. Chem. Phys.* **152**, 214108 (2020). — Companion history covering the ACES II → CFOUR branch and the ACES III branch's divergence.

*Note: Some full bibliographic details (exact volume/page/DOI for a few entries) reflect the most complete information available from indexed search sources at the time of writing; readers should verify exact citation details against the publisher record before formal citation.*

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the ACES / ACES III 	Quantum chemistry package specializing in coupled-cluster and many-body methods, with parallel implementations, also supporting DFT. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
