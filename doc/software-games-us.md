# GAMESS (US): An Exhaustive Review

## 1. Overview

**GAMESS** (General Atomic and Molecular Electronic Structure System — US branch) is a general-purpose, ab initio and density-functional quantum chemistry software package. It originated on October 1, 1977 as part of a National Resources for Computations in Chemistry (NRCC) project, and in 1981 the original codebase split into two independently developed lines — **GAMESS (US)**, maintained at Iowa State University, and **GAMESS (UK)** — which have since diverged substantially. A separate derivative, **Firefly** (formerly PC GAMESS), also descends from the original US code.

GAMESS (US) is developed and maintained by the **Gordon Research Group** at Iowa State University (with Ames Laboratory / DOE affiliation), in collaboration with a large distributed network of academic contributors worldwide. It is distributed as **source-available freeware**: the full Fortran/C/C++ source code is provided free of charge to academic and industrial users, but it is *not* released under an OSI-approved open-source license, since redistribution is restricted.

| Attribute | Detail |
|---|---|
| Developer | Gordon Research Group, Iowa State University / Ames Laboratory |
| Initial release | October 1, 1977 |
| Latest stable release | ~2024 (rolling "R" point releases; development also tracked on GitHub) |
| Languages | Fortran (majority, legacy), with newer C/C++ libraries (libcchem, EXESS) |
| Platforms | Linux, Unix-like systems (FreeBSD, macOS), Windows |
| License | Proprietary/source-available freeware (not OSI open source) |
| Code size | ~2+ million lines of code |
| User base | Reported to exceed 140,000 researchers in 100+ countries |
| Homepage | msg.chem.iastate.edu/gamess |

## 2. Historical Development

- **1977–1981**: Original GAMESS code assembled by M. Dupuis, D. Spangler, and J. J. Wendoloski as part of the NRCC Software Catalog (Program QG01), University of California, Berkeley.
- **1981**: Codebase forks into GAMESS (US) and GAMESS (UK).
- **1993**: Publication of the canonical reference paper (Schmidt et al., *J. Comput. Chem.* 14, 1347), which remains the mandatory citation for any GAMESS-based publication.
- **1999**: Introduction of the first Distributed Data Interface (DDI), by Graham Fletcher and Mike Schmidt, enabling large-array storage across distributed memory for parallel execution.
- **2000s**: Growth to ~650,000 lines of Fortran by 2005; addition of analytic Hessians, direct SCF, direct AO integrals, effective core potentials, coupled-cluster methods, and solvation models.
- **2009**: Divide-and-conquer methods from Waseda University (M. Kobayashi, T. Akama, H. Nakai) incorporated.
- **2011**: **QuanPol**, the QM/MM polarizable force field program, completed at the University of Nebraska–Lincoln by Hui Li, Nandun Thellamurege, and Dejun Si, and merged into GAMESS.
- **2010s–2020s**: Major push toward exascale-ready architectures under the DOE **Exascale Computing Project (ECP)**, including GPU/accelerator support (CUDA, HIP/ROCm), new C++ libraries (libcchem, libcint, and later **EXESS**), OpenMP/OpenACC directive-based portability, and fragmentation-based scaling methods (FMO, EFMO).
- **2020**: "Recent Developments in GAMESS" (Barca et al., *J. Chem. Phys.* 152, 154102) documents a decade of new methods.
- **2023**: "Novel Methods on Novel Architectures" (Barca et al., *J. Chem. Theory Comput.* 19, 7031) documents GPU-era developments, quasi-atomic orbital (QUAO) bonding analysis, and density-fitted/fragmented coupled-cluster scaling.
- **2025–2026**: Continued collaboration (e.g., with Old Dominion University and AWS) to containerize GAMESS for CPU- and GPU-optimized cloud HPC deployment, improving reproducibility and accessibility for its global user base.

## 3. Core Electronic-Structure Capabilities

### 3.1 Self-Consistent-Field (Reference) Wavefunctions
GAMESS supports the following SCF reference types, each combinable (to varying extent) with conventional, direct, and parallel integral evaluation:

- **RHF** — restricted Hartree-Fock (closed shell)
- **ROHF** — restricted open-shell Hartree-Fock
- **UHF** — unrestricted Hartree-Fock
- **GVB** — generalized valence bond
- **MCSCF** — multi-configuration SCF (including determinant-based and CSF-based variants)

For each, GAMESS can compute energies, analytic gradients, and (for most) analytic or numerical Hessians, enabling geometry optimization, transition-state searches, and vibrational analysis.

### 3.2 Density Functional Theory (DFT)
- A broad library of exchange-correlation functionals: LDA, GGA (e.g., BLYP, PBE), hybrid (B3LYP, PBE0), meta-GGA, and range-separated/long-range-corrected functionals.
- Analytic gradients and Hessians for DFT energies.
- **Time-Dependent DFT (TD-DFT)** for excited-state energies, oscillator strengths, and (via QuanPol/EFP coupling) solvatochromic shifts.

### 3.3 Post-Hartree–Fock / Correlation Methods
- **Møller–Plesset perturbation theory**: MP2 (with analytic gradients, and parallel algorithms for both energy and gradient evaluation), and higher-order variants.
- **Configuration Interaction (CI)**: including unitary-group-based CI, determinant-based full and truncated CI.
- **Coupled-Cluster (CC) theory**: CCSD, CCSD(T), and equation-of-motion (EOM-CC) methods for excited, ionized, and electron-attached states (including active-space and multi-particle/hole treatments such as double-ionization-potential/double-electron-attachment EOM-CC).
- Excited states additionally accessible via CI, EOM, or TD-DFT.

### 3.4 Fragmentation and Multiscale Methods (Native to GAMESS)
- **Fragment Molecular Orbital method (FMO)**: enables all-electron ab initio treatment of very large systems (proteins, nanoparticles) by decomposing them into fragments, each solved with full QM machinery, with inter-fragment electrostatics handled self-consistently. Demonstrated on systems exceeding 17,000 atoms.
- **Effective Fragment Molecular Orbital (EFMO)**: a more approximate/scalable variant combining FMO with the Effective Fragment Potential method, targeted at GPU and exascale execution; used for systems such as mesoporous silica nanoparticles solvated by thousands of water molecules.
- **Effective Fragment Potential (EFP) method**: a first-principles-derived, non-empirical polarizable model potential for explicit solvent/environment effects, combinable with QM, PCM, and QuanPol regions in layered models.
- **Elongation Method**: models the stepwise formation of aperiodic polymers.

### 3.5 Implicit and Explicit Solvation Models
GAMESS implements multiple continuum solvation approaches:
- Onsager cavity model / self-consistent reaction field (SCRF)
- Polarizable Continuum Model (PCM)
- Conductor-like Screening Model (COSMO)
- Surface and Simulation of Volume Polarization for Electrostatics (SS(V)PE)
- SMx-family solvation models (temperature-dependent variants)

These can be layered — e.g., QM region + EFP explicit-solvent particles + PCM continuum — within a single calculation.

### 3.6 QM/MM Capabilities, Including QuanPol
GAMESS offers several routes to hybrid QM/MM modeling:

- **QuanPol** — the primary, fully integrated QM/MM/polarizable-force-field/continuum-solvation program in GAMESS. Developed by Hui Li, Nandun Thellamurege, and Dejun Si (University of Nebraska–Lincoln, initial implementation completed August 2011 under NSF support), QuanPol performs combined QM/MM calculations in which:
  - QM methods available include **HF, DFT, GVB, MCSCF, MP2, and TD-DFT**.
  - MM interactions use **induced-dipole polarizable force fields**, self-consistently and variationally coupled to the QM wavefunction.
  - Continuum solvation is handled via **induced surface charges** (e.g., FixSol/SphSol-type models), also solved self-consistently with the QM density.
  - Supported/compatible force fields include user-specified parameter sets as well as standard force fields such as **MMFF94, CHARMM, AMBER, and OPLS-AA**.
  - Analytic energy gradients are implemented for QM/MM/continuum combinations at the MP2 and TD-DFT levels (via modified Z-vector response treatments), enabling geometry optimization and vibrational/spectroscopic property prediction in solution or protein environments.
  - Applications include studies of hydrogen-bonding networks in proteins (e.g., photoactive yellow protein), UV/vis solvatochromic shifts, and pKa prediction for drug-like molecules.
- **SIMOMM / IMOMM** — Surface/embedded ONIOM-style QM/MM via the plug-in **TINKER** molecular mechanics program.
- **NEO (Nuclear-Electronic Orbital)** — a plug-in extension treating select nuclei (e.g., protons) quantum mechanically alongside electrons.
- Combined **QM/EFP/PCM** layered models as an alternative to explicit polarizable MM for solvent/environment effects.

### 3.7 Molecular Properties
- Analytic and numerical vibrational frequencies, IR and Raman intensities
- Radiative transition probabilities and spin-orbit coupling constants
- Static and frequency-dependent polarizabilities
- Reaction path following (IRC) and saddle-point (transition-state) searches
- Scalar relativistic corrections
- Molecular dynamics/trajectory-based "dynamical" runs

## 4. Parallelization and High-Performance Computing

- **Distributed Data Interface (DDI)**: GAMESS's core software layer for parallel execution, allowing distributed-memory systems to pool memory for the large arrays intrinsic to quantum-chemical methods. First introduced in 1999 (Fletcher/Schmidt); DDI2 introduced later by Ryan Olson.
- **Generalized Distributed Data Interface (GDDI)**: layers parallelism at the fragment level (e.g., for FMO/EFMO), relying on MPI as the communication substrate.
- Newer C++ libraries — **libcchem**, **libcint**, and **EXESS** — implement Rys-quadrature two-electron integrals, HF, MP2, and CCSD(T) on both CPU and GPU architectures (CUDA/HIP), depending on BLAS/diagonalization libraries and (for libcchem/EXESS) Global Arrays, Eigen, and CUDA/HIP.
- Portability strategies include OpenMP and OpenACC directive-based approaches to support diverse exascale-class architectures.
- Under the DOE **Exascale Computing Project**, GAMESS has been a flagship application, targeting heterogeneous-catalysis problems (e.g., mesoporous silica nanoparticle models with explicit solvent) that require many-thousand-atom treatments.
- Recent (2025–2026) collaboration between Old Dominion University, Iowa State University, and AWS has produced CPU- and GPU-optimized **containerized** GAMESS builds for cloud HPC (including ARM/Graviton-targeted builds), aimed at improving reproducibility and accessibility across heterogeneous computing environments, with reported energy-efficiency gains (~30%) at modest (~10%) speed cost.

## 5. Ecosystem and Tools

- **MacMolPlt**: companion visualization/input-preparation GUI (originally macOS-focused, later cross-platform for Linux/Windows/macOS).
- **GamessQ**: a lightweight batch job manager for desktop use (macOS, Windows, Linux).
- Development repository hosted on **GitHub**, with a separate public/stable distribution channel from msg.chem.iastate.edu.
- Widely deployed on academic HPC clusters (e.g., via EasyBuild modules) and teaching platforms (e.g., ChemCompute, used in undergraduate computational chemistry laboratories).

## 6. Licensing and Access

GAMESS (US) source code is distributed at no monetary cost to academic and (with registration) industrial users, but under a restrictive **proprietary freeware** license rather than an open-source license: redistribution of the source is not generally permitted, and users must register to obtain the code from Iowa State University. This distinguishes it from GAMESS-UK (which has separate commercial licensing) and from fully open-source alternatives such as NWChem or PSI4.

## 7. Notable Application Domains

- Drug design and structure-activity relationship (SAR) studies via FMO-based protein–ligand binding-energy estimation.
- Fragment- and structure-based drug design (FBDD/SBDD).
- pKa prediction for pharmaceutically relevant molecules.
- Heterogeneous catalysis modeling (e.g., mesoporous silica nanoparticles).
- Excited-state and photochemistry studies (e.g., photoactive yellow protein chromophore environments).
- Materials science and polymer/aperiodic-system modeling via the Elongation Method.
- Undergraduate and graduate computational chemistry education.

## 8. Strengths and Limitations

**Strengths**
- Very broad method coverage (SCF variants through CC/EOM-CC) in a single package.
- Native, tightly integrated fragmentation methods (FMO/EFMO) and native polarizable QM/MM (QuanPol) — not bolt-on external interfaces.
- Free access to full source code; strong presence in HPC/exascale research.
- Long-standing, well-documented citation culture tied to peer-reviewed theory papers.

**Limitations**
- Not OSI-licensed open source; redistribution restrictions.
- Legacy Fortran codebase (though undergoing active C++/GPU modernization) can complicate contribution and long-term maintainability compared to newer codes designed natively for GPUs.
- Input format and workflow are less modern/scriptable out-of-the-box than some competitors (e.g., Psi4's Python API), though this is mitigated by third-party wrappers.

---

# Publications Related to GAMESS Theory and Methods

## Foundational / Package Description
1. Schmidt, M. W.; Baldridge, K. K.; Boatz, J. A.; Elbert, S. T.; Gordon, M. S.; Jensen, J. H.; Koseki, S.; Matsunaga, N.; Nguyen, K. A.; Su, S.; Windus, T. L.; Dupuis, M.; Montgomery, J. A. "General Atomic and Molecular Electronic Structure System." *Journal of Computational Chemistry* **1993**, *14*, 1347–1363. https://doi.org/10.1002/jcc.540141112 — *the mandatory citation for any GAMESS-based publication.*

2. Gordon, M. S.; Schmidt, M. W. "Advances in electronic structure theory: GAMESS a decade later." In *Theory and Applications of Computational Chemistry: The First Forty Years*, Dykstra, C. E.; Frenking, G.; Kim, K. S.; Scuseria, G. E., Eds.; Elsevier, 2005; pp 1167–1189.

3. Barca, G. M. J.; Bertoni, C.; Carrington, L.; Datta, D.; De Silva, N.; Deustua, J. E.; Fedorov, D. G.; Gour, J. R.; Gunina, A. O.; Guidez, E.; Harville, T.; Irle, S.; Ivanic, J.; Kowalski, K.; Leang, S. S.; Li, H.; Li, W.; Lutz, J. J.; Magoulas, I.; Mato, J.; Mironov, V.; Nakata, H.; Pham, B. Q.; Piecuch, P.; Poole, D.; Pruitt, S. R.; Rendell, A. P.; Roskop, L. B.; Ruedenberg, K.; Sattasathuchana, T.; Schmidt, M. W.; Shen, J.; Slipchenko, L.; Sosonkina, M.; Sundriyal, V.; Tiwari, A.; Galvez Vallejo, J. L.; Westheimer, B.; Włoch, M.; Xu, P.; Zahariev, F.; Gordon, M. S. "Recent Developments in the General Atomic and Molecular Electronic Structure System." *The Journal of Chemical Physics* **2020**, *152*, 154102. https://doi.org/10.1063/5.0005188

4. Barca, G. M. J.; Bertoni, C.; Fedorov, D. G.; Ivanic, J.; Leang, S. S.; Pham, B. Q.; Piecuch, P.; Slipchenko, L. V.; Xu, P.; Gordon, M. S. *et al.* "The General Atomic and Molecular Electronic Structure System (GAMESS): Novel Methods on Novel Architectures." *Journal of Chemical Theory and Computation* **2023**, *19*, 7031–7055. https://doi.org/10.1021/acs.jctc.3c00379

## QM/MM and QuanPol
5. Thellamurege, N. M.; Si, D.; Cui, F.; Zhu, H.; Lai, R.; Li, H. "QuanPol: A Full Spectrum and Seamless QM/MM Program." *Journal of Computational Chemistry* **2013**, *34*, 2816–2833. https://doi.org/10.1002/jcc.23435

6. Li, H.; Gordon, M. S. "Polarization Energy Gradients in Combined Quantum Mechanics, Effective Fragment Potential, and Polarizable Continuum Model Calculations." *The Journal of Chemical Physics* **2007**, *126*, 124112. https://doi.org/10.1063/1.2711154

7. Li, H. "Analytic Energy Gradient in Combined Second-Order Møller–Plesset Perturbation Theory and Polarizable Force Field Calculation." *The Journal of Chemical Physics* (relevant QM/MM/C MP2 gradient paper; see also "A new parallel algorithm of MP2 energy calculations" and "New parallel algorithm for MP2 energy gradient calculations," describing the Z-vector-based response treatment used for QM/MM/continuum MP2 properties).

8. Li, H. "Analytic Energy Gradient in Combined Time-Dependent Density Functional Theory and Polarizable Force Field Calculation" — TD-DFT/QM/MM/polarizable force field gradient theory underlying QuanPol excited-state QM/MM calculations.

9. Li, H.; Pomelli, C. S.; Jensen, J. H. "Continuum Solvation Gradients Within the Polarizable Continuum Model Combined with a Discrete Solvation Model" / related PCM–EFP integration papers underlying GAMESS's layered QM/EFP/PCM solvation approach.

10. Nakano, H. *et al.* — Divide-and-conquer method papers (Waseda University group; Kobayashi, Akama, Nakai), incorporated into GAMESS January 2009, providing linear-scaling QM/MM-adjacent fragmentation.

## Fragment Molecular Orbital (FMO) and EFMO
11. Fedorov, D. G.; Kitaura, K. "Extending the Power of Quantum Chemistry to Large Systems with the Fragment Molecular Orbital Method." *The Journal of Physical Chemistry A* **2007**, *111*, 6904–6914.

12. Fedorov, D. G.; Nagata, T.; Kitaura, K. "Exploring Chemistry with the Fragment Molecular Orbital Method." *Physical Chemistry Chemical Physics* **2012**, *14*, 7562–7577.

13. Pruitt, S. R.; Bertoni, C.; Brorsen, K. R.; Gordon, M. S. "Efficient and Accurate Fragmentation Methods." *Accounts of Chemical Research* **2014**, *47*, 2786–2794.

## Effective Fragment Potential (EFP)
14. Gordon, M. S.; Freitag, M. A.; Bandyopadhyay, P.; Jensen, J. H.; Kairys, V.; Stevens, W. J. "The Effective Fragment Potential Method: A QM-Based MM Approach to Modeling Environmental Effects in Chemistry." *The Journal of Physical Chemistry A* **2001**, *105*, 293–307.

15. Gordon, M. S.; Slipchenko, L.; Li, H.; Jensen, J. H. "The Effective Fragment Potential: A General Method for Predicting Intermolecular Interactions." *Annual Reports in Computational Chemistry* **2007**, *3*, 177–193.

## Coupled-Cluster / EOM-CC Methods in GAMESS
16. Piecuch, P.; Kowalski, K.; Pimienta, I. S. O.; Fan, P.-D.; Lodriguito, M.; McGuire, M. J.; Kucharski, S. A.; Kuś, T.; Musiał, M. "Method of Moments of Coupled-Cluster Equations: A New Formalism for Designing Accurate Electronic Structure Methods for Ground and Excited States." *Theoretical Chemistry Accounts* **2004**, *112*, 349–393.

17. Deustua, J. E.; Shen, J.; Piecuch, P. "Converging High-Level Coupled-Cluster Energetics by Monte Carlo Sampling and Moment Expansions." *Physical Review Letters* **2017**, *119*, 223003.

## Distributed Data Interface / Parallel Algorithms
18. Fletcher, G. D.; Schmidt, M. W.; Bode, B. M.; Gordon, M. S. "The Distributed Data Interface in GAMESS." *Computer Physics Communications* **2000**, *128*, 190–200.

19. Olson, R. M.; Bentz, J. L.; Kendall, R. A.; Schmidt, M. W.; Gordon, M. S. "A Novel Approach to Parallel Coupled Cluster Calculations: Combining Distributed and Shared Memory Techniques for Modern Cluster Based Systems." *Journal of Chemical Theory and Computation* **2007**, *3*, 1312–1328.

20. Bentz, J. L.; Olson, R. M.; Gordon, M. S.; Schmidt, M. W.; Kendall, R. A. "Coupled Cluster Algorithms for Networks of Shared Memory Parallel Processors." *Computer Physics Communications* **2007**, *176*, 589–600.

## GPU / Exascale Modernization
21. Leang, S. S.; Rendell, A. P.; Gordon, M. S. "Quantum Chemical Calculations Using Accelerators: Redesigning Three-Center Two-Electron Integral Evaluation for a GPU." *Journal of Chemical Theory and Computation* **2014**, *10*, 908–912.

22. Sosonkina, M.; Mateescu, G.; Xu, P.; Sattasathuchana, T.; Pham, B.; Gordon, M. S.; Leang, S. S. "Runtime Performance of a GAMESS Quantum Chemistry Application Offloaded to GPUs." (OSTI technical report / conference paper), 2024.

## Solvation Models Referenced by GAMESS
23. Klamt, A.; Schüürmann, G. "COSMO: A New Approach to Dielectric Screening in Solvents with Explicit Expressions for the Screening Energy and Its Gradient." *Journal of the Chemical Society, Perkin Transactions 2* **1993**, 799–805.

24. Miertuš, S.; Scrocco, E.; Tomasi, J. "Electrostatic Interaction of a Solute with a Continuum: A Direct Utilization of ab initio Molecular Potentials for the Prevision of Solvent Effects." *Chemical Physics* **1981**, *55*, 117–129.

25. Chipman, D. M. "Reaction Field Treatment of Charge Penetration." *The Journal of Chemical Physics* **2000**, *112*, 5558–5565. (SS(V)PE method.)

## Review of QM/QM-MM Capability in Drug Design
26. Fedorov, D. G.; Nagata, T.; Li, H.; *et al.* "GAMESS as a Free Quantum-Mechanical Platform for Drug Research." *Current Topics in Medicinal Chemistry* **2012**, review article discussing FMO and QM/MM (including QuanPol) applications to biochemical and drug-design systems.

---

*Note: References 6–10 above summarize the theoretical lineage of QuanPol's gradient/response formalisms as documented in the GAMESS literature and user manual (INTRO.DOC/REFS.DOC); consult the official GAMESS documentation (`REFS.DOC`) distributed with the source code for the complete, version-specific bibliography, as GAMESS's own manual maintains the authoritative and continuously updated reference list for each implemented method.*

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the GAMESS (US) 	Free general quantum chemistry package supporting DFT, HF, and many post-HF methods, including QM/MM via the QuanPol module. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
