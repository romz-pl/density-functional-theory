# Qbox: A Scalable Plane-Wave/Pseudopotential Code for Large-Scale First-Principles Molecular Dynamics

## 1. Overview

Qbox is a C++/MPI(/OpenMP) implementation of first-principles molecular dynamics (FPMD) based on the plane-wave, pseudopotential formalism, designed for operation on large parallel computers. It solves the Kohn–Sham equations of density functional theory (DFT) self-consistently and propagates ionic positions using forces derived from the resulting electronic structure, enabling ab initio simulation of solids, liquids, nanostructures, and molecules without empirical parameters.

Qbox is developed and maintained by the **Gygi research group** (originally at Lawrence Livermore National Laboratory, now at the **University of California, Davis**), with a growing list of external contributors. It is free, open-source software released under the **GNU General Public License (GPL)**.

| | |
|---|---|
| **Original author** | François Gygi |
| **Initial release** | 2003 |
| **Latest stable release** | 1.78.4 (24 April 2025) |
| **Language** | C++ (with MPI and OpenMP) |
| **License** | GPL |
| **Repository** | github.com/qboxcode |
| **Homepage** | qboxcode.org |
| **Primary developing institution** | University of California, Davis (Gygi Research Group / Electronic Structure Laboratory) |

## 2. Scientific and Numerical Method

- **Basis set:** Plane waves, providing a systematically improvable, unbiased representation of the electronic wavefunctions and controlled convergence via a single energy cutoff parameter.
- **Electron–ion interaction:** Norm-conserving pseudopotentials, including the **SG15** set of Optimized Norm-Conserving Vanderbilt (ONCV) pseudopotentials, which were built to reproduce all-electron results with high accuracy and are distributed for use with Qbox. Pseudopotentials from the PseudoDojo project can also be converted into Qbox's native QSO format via a bundled `upf2qso` utility.
- **Exchange–correlation functionals:** LDA and GGA functionals, together with range-separated and screened hybrid functionals (e.g., HSE, PBE0), BHandHLYP, and a spin-polarized implementation of B3LYP.
- **Ionic/electronic dynamics:**
  - Born–Oppenheimer molecular dynamics in the microcanonical (NVE) or canonical (NVT) ensemble.
  - Car–Parrinello-style fictitious dynamics for the electronic degrees of freedom.
  - Constrained molecular dynamics for thermodynamic integration.
  - Several wavefunction-update/minimization algorithms: steepest descent, preconditioned steepest descent (PSD), PSD with Anderson acceleration (PSDA), and Jacobi–Davidson (JD).
- **Response and analysis capabilities:**
  - Efficient computation of maximally localized Wannier functions (MLWFs), used for polarization, dipole, and localization analysis.
  - Electronic structure calculations in the presence of arbitrary/constant external electric fields, and computation of electronic polarizability.
  - Partial charge/spin population analysis in atom-centered spheres.
  - Constant-pressure (NpT) molecular dynamics.
  - A "client–server" execution mode allowing external drivers to steer Qbox interactively, used, for example, in nudged-elastic-band (NEB) reaction-pathway calculations.
  - Support for executing arbitrary external commands at regular intervals during an MD run, for on-the-fly analysis or coupling to other tools.

## 3. Software Architecture

The architecture of Qbox is described in detail in Gygi's 2008 *IBM Journal of Research and Development* paper, which explains that Qbox is built upon well-optimized parallel numerical libraries such as the Basic Linear Algebra Communication Subprograms (BLACS) and the Scalable Linear Algebra Package (ScaLAPACK), and also features an Extensible Markup Language (XML) interface built on the Apache Xerces-C library. Key architectural characteristics include:

- **Language and design.** Qbox is written in modern C++, applying object-oriented design patterns (as documented in the architecture paper's references to Gamma et al.'s *Design Patterns* and Stroustrup's *The C++ Programming Language*) to separate physics, numerics, and parallel data distribution concerns.
- **Parallel decomposition.** Qbox distributes plane-wave coefficients and electronic states across a two-dimensional (and, in later releases, hierarchical MPI/OpenMP) process grid, using ScaLAPACK/BLACS for dense linear algebra (orthogonalization, subspace diagonalization) and custom parallel FFTs for the transforms between real and reciprocal space that dominate the plane-wave method's cost.
- **I/O and interoperability.** Qbox's native input is a simple line-oriented command interpreter, while output is a well-defined **XML** document conforming to a public XML Schema (hosted at quantum-simulation.org). This format is directly supported by the **NOMAD** materials-data repository, which promotes reuse and long-term archival/interoperability of Qbox results. A companion suite of post-processing and visualization tools is distributed alongside the code.
- **Interactive/client–server operation.** Beyond batch execution from an input script, Qbox can run interactively, reading commands from standard input, and can be driven as a "server" by an external client process — a design that has been reused for workflow tools such as NEB reaction-path searches.
- **Hybrid parallelism.** More recent releases are described as a C++/MPI/OpenMP implementation, reflecting the addition of on-node threading to complement the original pure-MPI design as target machines moved to many-core nodes.

## 4. Scalability and Performance Record

Scalability on very large parallel machines has been Qbox's defining characteristic since its inception.

- **Blue Gene/L era (2005–2006).** Qbox was used to run first-principles simulations of high-Z metallic systems on the full 65,536-node (128k CPU) Blue Gene/L system at Lawrence Livermore National Laboratory, sustaining a peak performance of 207.3 TFlop/s, corresponding to 56.5% of the theoretical peak of the entire machine. This result — "Large-Scale Electronic Structure Calculations of High-Z Metals on the BlueGene/L Platform" — won the **2006 ACM/IEEE Gordon Bell Prize for Peak Performance**.
- **Architecture paper benchmark.** The 2008 IBM J. Res. Dev. architecture paperspecifically discusses the case of the IBM Blue Gene/L platform, on which Qbox was run using up to 65,536 nodes, and also discusses future design challenges for upcoming petascale computers.
- **Subsequent platforms.** Qbox has since been ported to and run on Blue Gene/Q, Cray XT-5, Cray XE-6, Cray XC30, NERSC's Perlmutter, and a variety of Linux/Intel clusters, and is used in ongoing large-scale materials-science and liquid-state simulation projects.
- **2022 Roadmap perspective.** In the community "Roadmap on Electronic Structure Codes in the Exascale Era" (2022), the Qbox contribution states that the codewas designed for scalability on thousands of tasks, involving tens of thousands of processor cores, and is routinely used to perform molecular dynamics simulations of systems including several hundred atoms, with features including constant-temperature (NVT) and constant-pressure (NpT) MD, computation of maximally localized Wannier functions, hybrid-DFT exchange-correlation functionals, and computation of electronic response to arbitrary periodic perturbations.

## 5. Development History and Governance

- **Origins (2003).** Qbox was created as a from-scratch, modern rewrite of an FPMD code, explicitly intended, in the words of its architecture paper's reference list, to seize "the opportunity to create a modern code that profits from current design ideas and programming tools," rather than incrementally extend legacy Fortran-era plane-wave codes.
- **Institutional history.** Development began at LLNL's Center for Applied Scientific Computing, in a context where FPMD was explicitly targeted as a workload for the leadership-class Blue Gene/L system. Development subsequently moved with François Gygi to UC Davis, where the code continues to be developed in the Gygi Research Group / Electronic Structure Laboratory (ESLab).
- **Contributors.** In addition to François Gygi, publicly credited developers include Ivan Duchemin, Jun Wu, Quan Wan, William Dawson, Martin Schlipf, He Ma, and Michael LaCount, among others.
- **Funding.** Development has been supported by the U.S. Department of Energy Office of Basic Energy Sciences (including grant DE-SC0008938) and, more recently, through a collaboration with the DOE **Midwest Integrated Center for Computational Materials (MICCoM)**.
- **Distribution.** Qbox source code is hosted on an ESLAB GitLab server and mirrored/available on GitHub at github.com/qboxcode, under the GPL.
- **Release cadence.** The code has seen continuous, incremental releases over two decades; recent versions include 1.74.4 (Python-compatibility release, 2022), 1.75.0 (band-alignment update, 2022), 1.75.1 ("enable large samples," 2023), 1.76.1 (2023), 1.78.2 (dipoles in arbitrary unit cells; improved MLWF center accuracy in non-cubic cells, 2024), and 1.78.4 (2025).

## 6. Representative Scientific Applications

Qbox has been applied across condensed matter physics, high-pressure/planetary science, liquid-state physics, and materials chemistry, including:

- Simulations of **liquids** (notably water and aqueous interfaces, including ab initio MD of an electrified silicon–water interface).
- **High-pressure/extreme-condition physics**, e.g., melting curves of light elements/hydrides and behavior of matter at planetary-interior conditions.
- **Semiconductor nanostructures and nanocrystals**, e.g., optical properties of silicon nanocrystals.
- **Point defects in wide-gap semiconductors**, e.g., charge-state and entropic effects governing divacancy formation/dynamics in 3C-SiC, using Qbox's client–server NEB capability.
- General **materials science** and solid-state simulations distributed across national laboratory and university groups using the code.

## 7. Position in the Electronic-Structure Software Landscape

Qbox belongs to the family of plane-wave/pseudopotential DFT codes (alongside, e.g., Quantum ESPRESSO, VASP, CP2K/Quickstep, and ABINIT), distinguished within that family primarily by:

- An architecture designed from the outset (2003) for **extreme parallel scalability** on the largest available supercomputers, rather than being retrofitted for parallelism onto an earlier serial/small-cluster design.
- A clean, modern **C++/object-oriented codebase** rather than legacy Fortran, chosen explicitly to leverage contemporary software-engineering practices.
- A native **XML-based output format** with a public schema, enabling structured post-processing and integration with data repositories such as NOMAD.
- A demonstrated, peer-recognized (Gordon Bell Prize) track record of near-full-machine-scale performance on leadership supercomputers, which has historically served as a stress test and proof of concept for scalable FPMD more broadly.

## 8. Publications Related to Qbox's Theory and Implementation

### Core architecture and method papers
1. F. Gygi, **"Architecture of Qbox: A scalable first-principles molecular dynamics code,"** *IBM Journal of Research and Development* **52**(1/2), 137–144 (2008). — The principal reference describing Qbox's design, parallel numerical-library dependencies (BLACS/ScaLAPACK), and XML interface.
2. F. Gygi, E. W. Draeger, M. Schulz, B. R. de Supinski, J. A. Gunnels, V. Austel, J. C. Sexton, F. Franchetti, S. Kral, C. Ueberhuber, and J. Lorenz, **"Large-Scale Electronic Structure Calculations of High-Z Metals on the BlueGene/L Platform,"** *Proceedings of the ACM/IEEE Conference on Supercomputing (SC'06)* (2006). — The Gordon Bell Prize–winning paper reporting 207.3 TFlop/s sustained performance on 65,536 Blue Gene/L nodes.
3. F. Gygi, R. K. Yates, J. Lorenz, E. W. Draeger, F. Franchetti, C. W. Ueberhuber, B. R. de Supinski, S. Kral, J. A. Gunnels, and J. C. Sexton, **"Large-Scale First-Principles Molecular Dynamics Simulations on the BlueGene/L Platform using the Qbox Code,"** *Proceedings of the 2005 ACM/IEEE Conference on Supercomputing (SC'05)*, IEEE Computer Society (2005). LLNL/OSTI technical report also available.
4. F. Gygi, **"Qbox: a large-scale parallel implementation of First-Principles Molecular Dynamics,"** LLNL preprint (2003) — the original code-introduction preprint.
5. F. Gygi, **"Large-Scale First-Principles Molecular Dynamics: Moving from Terascale to Petascale Computing,"** *Journal of Physics: Conference Series* **46**, 268–277 (2006).
6. F. Gygi, **"The Qbox first-principles molecular dynamics code,"** in *Roadmap on Electronic Structure Codes in the Exascale Era*, V. Gavini et al., Section 10, *Electronic Structure* (2023); preprint arXiv:2209.12747. — Current-status and exascale-development-priorities overview.

### Underlying Car–Parrinello / FPMD method
7. R. Car and M. Parrinello, **"Unified Approach for Molecular Dynamics and Density-Functional Theory,"** *Physical Review Letters* **55**, 2471–2474 (1985). — Foundational FPMD method paper underlying Qbox's algorithmic approach.
8. D. G. Anderson, **"Iterative Procedures for Nonlinear Integral Equations,"** *Journal of the ACM* **12**(4), 547–560 (1965). — Basis for the Anderson-acceleration scheme used in Qbox's PSDA wavefunction-update algorithm.

### Pseudopotentials used by Qbox
9. M. Schlipf and F. Gygi, **"Optimization algorithm for the generation of ONCV pseudopotentials,"** *Computer Physics Communications* **196**, 36–44 (2015). — Describes the SG15 ONCV pseudopotential set distributed for use with Qbox.

### Wannier-function methodology implemented in Qbox
10. Reference implementation of maximally localized Wannier functions in Qbox follows the algorithm described in *Computer Physics Communications* **155**, 1 (2003) (one-sided iterative Jacobi simultaneous-diagonalization approach), as cited in the Qbox User Guide.

### Selected applications demonstrating/validating Qbox's method and scalability
11. T. Ogitsu, E. Schwegler, F. Gygi, and G. Galli, **"Melting of Lithium Hydride Under Pressure,"** *Physical Review Letters* **91**, 175502 (2003).
12. M. Allesch, E. Schwegler, F. Gygi, and G. Galli, **"A First Principles Simulation of Rigid Water,"** *Journal of Chemical Physics* **120**(11), 5192–5198 (2004).
13. F. Gygi and G. Galli, **"Ab Initio Simulation in Extreme Conditions,"** *Materials Today* **8**(11), 26–32 (2005).
14. A. A. Correa, S. A. Bonev, and G. Galli, **"Carbon under Extreme Conditions: Phase Boundaries and Electronic Properties from First-Principles Theory,"** *Proceedings of the National Academy of Sciences* **103**(5), 1204–1208 (2006).
15. L. Dal Negro, S. Hamel, N. Zaitseva, J. H. Yi, A. Williamson, M. Stolfi, J. Michel, G. Galli, and L. C. Kimerling, **"Synthesis, Characterization, and Modeling of Nitrogen-Passivated Colloidal and Thin Film Silicon Nanocrystals,"** *IEEE Journal of Selected Topics in Quantum Electronics* **12**(6), 1151–1163 (2006).
16. C. Zhang, F. Gygi, and G. Galli, **"Charge state and entropic effects affecting the formation and dynamics of divacancies in 3C-SiC,"** *Physical Review Materials* **8**, 046201 (2024). — Demonstrates Qbox's client–server NEB capability.

### Interoperability / ecosystem context
17. F. Gygi, I. Duchemin, D. Donadio, et al., **"Code interoperability extends the scope of quantum simulations,"** *npj Computational Materials* **7**, 50 (2021). — Discusses Qbox in the context of a broader interoperable electronic-structure/molecular-dynamics software ecosystem.

---

*Note: for the primary architecture reference (item 1) and the Gordon Bell–winning benchmark paper (item 2), consult the original publications directly for full author lists, equations, and benchmark data, as only the abstracts and select passages were used here.*

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the Qbox 	Scalable plane-wave/pseudopotential code designed for large-scale parallel first-principles molecular dynamics on supercomputers. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
