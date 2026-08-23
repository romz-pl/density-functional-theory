# ORCA — A Comprehensive Review

## 1. Overview

ORCA is a general-purpose quantum chemistry package featuring a variety of methods including semi-empirical, density functional theory, many-body perturbation, coupled cluster, and multireference methods, providing an easy-to-learn input structure and high accessibility of quantum chemical approaches and workflows. In addition to being efficient, user friendly, and, to the largest extent possible, platform independent, ORCA features a number of methods that are either unique to it or were first implemented in the course of its development — notably a range of spectroscopic and magnetic properties, and the linear- or low-order single- and multi-reference local correlation methods based on pair natural orbitals (domain-based local pair natural orbital methods).

ORCA features density functional theory, a range of wavefunction-based correlation methods, semi-empirical methods, and even force-field methods, together with a range of solvation and embedding models and a complete intrinsic QM/MM engine. A specialty of the program has always been a focus on transition metals and spectroscopy, as well as applicability to "real-life" chemical applications involving systems with a few hundred atoms.

## 2. Development History and Governance

- Development of ORCA started in 1997, while Frank Neese was a postdoc at Stanford University, and continued as he moved through the University of Bonn and eventually the Max Planck Institute for Coal Research (MPI KoFo), where the program is now primarily maintained.
- The program is mainly developed by Frank Neese, the Department of Molecular Theory and Spectroscopy at MPI KoFo, and FACCTs GmbH, which also manages commercial licensing to industry.
- FACCTs GmbH was co-founded in 2016 by Christoph Riplinger as a Max Planck Society spin-off, providing industrial users with tailored support while ensuring free distribution continues for non-commercial research.
- The core development team includes contributions from Ute Becker, Dmytro Bykov, Dmitry Ganyushin, Andreas Hansen, Robert Izsak, Dimitrios G. Liakos, Christian Kollmar, Simone Kossmann, Dimitrios A. Pantazis, Taras Petrenko, Christoph Reimann, Christoph Riplinger, Michael Roemelt, Barbara Sandhöfer, Igor Schapiro, Kantharuban Sivalingam, Frank Wennmohs, Boris Wezisla, and external collaborators including Mihály Kállay, Stefan Grimme, and Edward Valeev.
- Growth of the user base has been rapid: over 22,000 registered users as of January 2020; 80,000 by January 2025; and more than 100,000 registered academic users as of the current FACCTs listing, describing ORCA as the fastest-growing quantum chemistry package to date.

## 3. Version History and Major Milestones

| Version | Date | Key Notes |
|---|---|---|
| Early versions | 1997–2010s | Grew from Neese's postdoctoral project into a general-purpose package |
| 4.x series | ~2017–2021 | Established DLPNO methods, embedding, extensive spectroscopy |
| 5.0 | July 2021 | Introduced the SHARK integral package, described as highly performant, and a new programming paradigm — the "loop-kernel-consumer" (LKC) concept — enabling compact, streamlined code; at release, ORCA was already the second-most used quantum chemistry suite worldwide with roughly 35,000 academic users |
| 6.0 | July 2024 | Described as a major turning point, representing a near-complete rewrite of the code base, yielding major performance improvements, a clean and efficient code base facilitating future development, a large amount of new functionality, and new interface capabilities for interoperability with other quantum chemistry packages. Efficiency gains included a 2–4× faster chain-of-spheres exchange (COSX) algorithm, up to 20% acceleration in the split-RI-J method, and roughly 20% fewer SCF convergence cycles, enabling routine SCF calculations for systems exceeding 10,000 basis functions via features such as the LeanSCF module |
| 6.0.1 | November 2024 | Bug-fix / stability release |
| 6.1.0 | 17 June 2025 | Delivered optimizations for large-scale simulations via automated fragmentation and multiscale setups, alongside analytical computation of Raman intensities |
| 6.1.x line | 2025–2026 | Continued incremental refinement; introduction of the official Python interface (OPI), first supported from ORCA 6.1 |

The most recent major release event was accompanied by a newly published user manual of roughly 1,500 pages, reflecting the scope of the package, which is now developed by a team of around 30 developers and international collaborators.

## 4. Licensing Model

ORCA follows a "free-for-academia, paid-for-industry" dual-licensing model:

- ORCA is and will remain free for academic and personal use; qualifying users can download the latest Linux version directly, while a commercial license or free trial can be obtained through FACCTs.
- From the start, ORCA was freely available to academic researchers — a feature the developers have always considered essential and non-negotiable. Given rapidly growing industrial demand, FAccTs was founded to license ORCA to commercial users, while ensuring that academic access remains free in the long term.
- The academic End User License Agreement is restrictive about scope: the software may be used exclusively in ACADEMIA for ACADEMIC PURPOSES and for PRIVATE USE, and "ACADEMIA" is defined as schools, colleges, universities, and Max Planck Society institutes, explicitly excluding institutions connected to the military; "PRIVATE USE" explicitly excludes any for-profit, non-profit, or governmental organizational use, and excludes any use that results in or contributes to commercial activity.
- Users must register on the official ORCA Forum to obtain binaries; academic labs register on the ORCA website and non-profit users may contact the licensing team directly for institutional access.
- Auxiliary open-source tooling is licensed separately: the ORCA Python Interface (OPI) is released publicly under GPL-3.0, allowing anyone to use, modify, and distribute it freely, while FACCTs retains the ability to relicense contributions under proprietary terms for commercial distributions.

## 5. Core Methodology and Feature Set

### 5.1 Electronic-structure methods
ORCA features extensive capabilities in Hartree–Fock, DFT, single-reference correlation, and multi-reference correlation methods, plus semi-empirical quantum chemical methods (NDO methods such as MNDO, as well as the more recently developed XTB method from the Grimme group) and even classical force-field calculations. Multilevel and embedding approaches are available to bridge the domains of applicability of the many available methods, along with a large number of methods for spectroscopic property prediction, most of which can be combined with an extended set of supporting features to solve real chemical problems.

### 5.2 DFT and TD-DFT
DFT is one of ORCA's most heavily used engines, spanning the standard rungs of exchange-correlation functionals combined with dispersion corrections, and extending to time-dependent DFT (TD-DFT) for excited-state energies, absorption/emission spectra, and related response properties — an area closely tied to the package's spectroscopic focus.

### 5.3 Correlated wavefunction methods and DLPNO
A hallmark ORCA capability is the domain-based local pair natural orbital (DLPNO) family of local correlation methods, which make coupled-cluster-quality results accessible for much larger systems than canonical implementations allow:

- The severe scaling of coupled cluster methods has been a major research problem, motivating the development of local correlation methods; ORCA features DLPNO coupled-cluster methods that make CC calculations feasible on large systems, building on the earlier local pair natural orbital (LPNO) methodology.
- DLPNO-CCSD(T) allows single-point energy calculations for systems with hundreds of atoms while retaining essentially the accuracy of canonical CCSD(T), with errors typically smaller than 1 kcal/mol for relative energies; the dimension of the pair-natural-orbital space is controlled by the TCutPNO threshold.
- DLPNO-CCSD(T) exploits molecular sparsity via localized orbitals for a linear-scaling approximation to canonical CCSD(T); at present, per the PSI4 documentation comparing implementations, DLPNO-CCSD/(T) is only available for closed-shell RHF references and does not yet implement molecular point-group symmetry, running instead in C1 symmetry.
- Extensions include: implicit solvation coupling with C-PCM, where the overhead from added solvent terms amounts to less than 5% of the equivalent gas-phase job time; Local Energy Decomposition (LED) analysis, decomposing DLPNO-CCSD(T) interaction energies into geometric/electronic preparation, electrostatics, exchange, charge polarization, and London dispersion terms, later extended to open-shell systems calculated at the UHF-DLPNO-CCSD(T) level; and DLPNO-based tailored coupled-cluster methods (DLPNO-TCCSD/(T)) integrating DMRG references.
- Accuracy has been benchmarked extensively, e.g. a comprehensive study over the GMTKN55 benchmark superset examined closed- and open-shell systems and barrier-height subsets, with mean absolute deviations reported for "Tight" DLPNO-CCSD(T) settings, and a dedicated study on transition-metal complexes identified semicore correlation and dynamic-correlation-induced orbital relaxation as the two principal error sources, proposing a strategy to correct for both. Independent assessments have also flagged limitations: the DLPNO approximation's fidelity for non-covalent interaction energies in large supramolecular complexes has not been thoroughly vetted, and accuracy is known to deteriorate as molecular size increases, particularly at looser truncation settings.

### 5.4 Multireference and other methods
ORCA supports multireference wavefunction methods (e.g., CASSCF/NEVPT2-type approaches), semi-empirical methods, QM/MM, and a variety of solvation/embedding models, positioning it as a genuinely general-purpose package rather than a single-method specialist code.

### 5.5 Spectroscopy and magnetic properties
Consistent with its origin in spectroscopy-oriented groups, ORCA is noted for a specific emphasis on spectroscopic properties of open-shell molecules, useful not only to computational chemists but also to experimental chemists, physicists, and biologists interpreting experimental data with the help of calculations. EPR, NMR, and spin-orbit coupling property modules are integrated (in ORCA 6, several of these — formerly separate modules such as `orca_scfhess`, `orca_eprnmr`, and `orca_soc` — were rendered superfluous by the new unified computational philosophy introduced with the rewrite).

## 6. Performance and Efficiency

Efficiency has been a consistent design goal:

- Great effort has been invested in making accurate calculations as fast as possible; developments such as DLPNO and RIJCOSX can speed up calculations by orders of magnitude and even yield linear-scaling behavior with system size.
- The ORCA 6.0 rewrite delivered a 2–4× faster COSX algorithm and up to 20% acceleration in split-RI-J, alongside roughly 20% fewer SCF convergence cycles, enabling routine SCF treatment of systems with more than 10,000 basis functions.
- Standard density-fitting/RI approaches (RI-J, RI-JK, RIJCOSX) are used throughout to accelerate integral evaluation for HF, DFT, and post-HF methods.

## 7. Usability

ORCA is widely regarded as unusually accessible relative to peer packages, owing to:
- A plain-text, keyword-driven input format described repeatedly across sources as user-friendly, helpful not only for computational chemists but for experimentalists interested in interpreting spectroscopic data.
- A scripting layer, "Compound," for chaining multi-step workflows, and, more recently, the ORCA Python Interface (OPI), introduced because accessibility and automation of quantum chemical tasks have become increasingly important, driven by growing computational resources and demand for large-scale, high-quality data for machine-learning method development. OPI is compatible with ORCA ≥ 6.1.1 and Python ≥ 3.11, and is released as an open-source community project.
- FACCTs' commercial-facing workflow tool "Weasel," described as intended to bring quantum chemical simulations — which can otherwise seem accessible only to experts — closer to a broader base of users.

## 8. Platform Support and Deployment

ORCA is written in C++ and supports Linux, Microsoft Windows, and macOS. It is precompiled and distributed as binaries rather than requiring users to compile from source. It runs on essentially all major HPC systems used in computational chemistry (e.g., NERSC, OSC, HPC2N, and numerous university clusters), typically parallelized via MPI; several site documentation pages note ORCA's MPI usage is somewhat atypical — for example, despite using MPI, ORCA is generally only able to run on a single node due to the non-standard way it uses MPI, unlike conventional multi-node MPI programs.

## 9. Applications and Domains

Its main field of application is larger molecules, transition metal complexes, and their spectroscopic properties, using standard Gaussian basis functions with full parallelization. Illustrative published applications surfaced in this research include:
- Benchmarking DLPNO approaches for predicting phosphorescence energies in aromatic systems, using B3LYP-D3/def2-TZVPP-optimized geometries with DLPNO-CCSD(T) single-point refinement.
- Assessing DLPNO-CCSD(T) accuracy for noncovalent bond dissociation enthalpies in coinage-metal (Cu⁺, Ag⁺, Au⁺) cation complexes.
- Complete-basis-set DLPNO-CCSD(T) extrapolations to characterize vibrational properties of benzene adsorbed on ordered water ice surfaces.
- Use as a reference-level engine (alongside MRCC) for generating delta-learned coupled-cluster training data for condensed-phase liquid water simulations.
- Systematic conformational benchmarking of piperazine-based drug-design scaffolds against DLPNO-CCSD(T)/CBS reference energies.

## 10. Strengths

- Genuinely broad method coverage (semi-empirical → DFT/TD-DFT → correlated wavefunction → multireference) in one package.
- Best-in-class accessibility to near-CCSD(T) accuracy on large systems via DLPNO methods, at a fraction of canonical computational cost.
- Consistently competitive-to-leading raw performance due to continual algorithmic investment (SHARK integrals, RIJCOSX, LeanSCF, the 6.0 rewrite).
- Deep, long-standing specialization in spectroscopic and magnetic properties (EPR, NMR, optical spectra), valuable for interpreting experimental data.
- Free academic access with a large, active user community and forum-based support, plus a growing modern Python-native workflow ecosystem (OPI).
- Rapid release cadence and transparent versioned documentation (extensive user manual).

## 11. Limitations and Criticisms

- The academic license is legally restrictive: use is confined strictly to academic institutions and non-commercial private use, with explicit exclusions for any organization or activity connected to commercial, governmental, or for-profit work — meaning many real-world users must obtain a separate commercial license through FACCTs.
- DLPNO methods carry accuracy caveats that require user awareness: fidelity for large-system non-covalent interactions is not fully vetted and known to degrade with increasing molecular size, and truncation thresholds are, by design, opaque — the ORCA manual states that DLPNO-CCSD truncation parameters "should almost not be touched" and are intentionally undocumented, a practice that has drawn explicit methodological criticism from outside groups for discouraging user awareness of numerical settings that may materially affect results.
- Certain DLPNO-CC variants are currently restricted to closed-shell references and lack point-group symmetry exploitation, limiting applicability for some open-shell or highly symmetric systems relative to canonical implementations.
- As a large, continuously evolving C++ codebase undergoing a "near-complete rewrite" as recently as 2024, backward compatibility and reproducibility of older workflows across major versions can require care (module structure changes noted with ORCA 6.0's elimination of several legacy modules).
- Being closed-source (aside from the newer OPI interface component), the core ORCA codebase itself is not independently auditable/modifiable in the way fully open-source competitors (e.g., PSI4, NWChem) are, which can be a barrier for method-development groups needing to inspect or extend low-level code directly.

## 12. Comparison Context

Independent groups reproducing DLPNO-type methods elsewhere have used ORCA as the reference benchmark implementation, e.g. an open-source DLPNO-CCSD(T) implementation using a t1-transformed Hamiltonian was developed explicitly as an analogue to the ORCA implementation, reproducing its linear-scaling integral generation and contraction while simplifying the algorithm, and PSI4's own DLPNO-CCSD implementation is directly compared against ORCA's in its documentation, noting PSI4 typically recovers slightly more correlation energy at a given nominal PNO_CONVERGENCE setting due to additional cutoff criteria, though the two cannot be made to match exactly due to differing T1 formulations. This reflects ORCA's position as something close to a field-defining reference implementation for local correlation methods.

---

# Key Publications on ORCA's Theory and Implementation

**Primary software/methods overview papers:**

1. Neese, F. "The ORCA program system." *WIREs Computational Molecular Science* **2012**, *2*(1), 73–78. https://doi.org/10.1002/wcms.81 — original comprehensive program description.
2. Neese, F.; Wennmohs, F.; Becker, U.; Riplinger, C. "The ORCA quantum chemistry program package." *Journal of Chemical Physics* **2020**, *152*(22), 224108. https://doi.org/10.1063/5.0004608 — special software-issue overview covering features up to v4.2.
3. Neese, F. "Software Update: The ORCA Program System — Version 6.0." *WIREs Computational Molecular Science* **2025**, *15*(2), e70019. https://doi.org/10.1002/wcms.70019 — describes the philosophy and new features of the 6.0 rewrite.
4. Neese, F. "ORCA, an Ab Initio, Density Functional and Semi-Empirical Program Package." Version documentation, 2008 (early formal program description).

**DLPNO / local correlation theory:**

5. Riplinger, C.; Neese, F. "An efficient and near linear scaling pair natural orbital based local coupled cluster method." *Journal of Chemical Physics* **2013**, *138*, 034106.
6. Riplinger, C.; Sandhoefer, B.; Hansen, A.; Neese, F. "Natural triple excitations in local coupled cluster calculations with pair natural orbitals." *Journal of Chemical Physics* **2013**, *139*, 134101.
7. Riplinger, C.; Pinski, P.; Becker, U.; Valeev, E. F.; Neese, F. "Sparse maps — A systematic infrastructure for reduced-scaling electronic structure methods. II. Linear scaling domain based pair natural orbital coupled cluster theory." *Journal of Chemical Physics* **2016**, *144*, 024109.
8. Guo, Y.; Riplinger, C.; Becker, U.; Liakos, D. G.; Minenkov, Y.; Cavallo, L.; Neese, F. "Communication: An improved linear scaling perturbative triples correction for the domain based local pair-natural orbital based singles and doubles coupled cluster method [DLPNO-CCSD(T)]." *Journal of Chemical Physics* **2018**, *148*, 011101.
9. Altun, A.; Neese, F.; Bistoni, G. "Extrapolation to the limit of a complete pair natural orbital space in local coupled-cluster calculations." *Journal of Chemical Theory and Computation* **2020**, *16*, 6142–6149. https://doi.org/10.1021/acs.jctc.0c00344
10. Altun, A.; Riplinger, C.; Neese, F.; Bistoni, G. "Exploring the accuracy limits of PNO-based local coupled-cluster calculations for transition-metal complexes." *Journal of Chemical Theory and Computation* **2023**, *19*(5). https://doi.org/10.1021/acs.jctc.3c00087
11. Liakos, D. G.; Guo, Y.; Neese, F. "Comprehensive benchmark results for the domain based local pair natural orbital coupled cluster method (DLPNO-CCSD(T)) for closed- and open-shell systems." *Journal of Physical Chemistry A* **2020**, *124*(1). https://doi.org/10.1021/acs.jpca.9b05734

**Local energy decomposition (LED):**

12. Altun, A.; Neese, F.; Bistoni, G. "Local energy decomposition analysis of hydrogen-bonded dimers within a domain-based pair natural orbital coupled cluster study." *Beilstein Journal of Organic Chemistry* **2018**, *14*, 919–929. https://doi.org/10.3762/bjoc.14.79
13. Altun, A.; Saitow, M.; Neese, F.; Bistoni, G. "Local energy decomposition of open-shell molecular systems in the domain-based local pair natural orbital coupled cluster framework." *Journal of Chemical Theory and Computation* **2019**. https://doi.org/10.1021/acs.jctc.8b01145

**Implicit solvation with DLPNO:**

14. "Implicit solvation in domain based pair natural orbital coupled cluster (DLPNO-CCSD) theory." *Journal of Chemical Physics* (2021). PMID: 34347890.

**Integral technology and SCF:**

15. Neese, F.; Wennmohs, F.; Hansen, A.; Becker, U. "Efficient, approximate and parallel Hartree–Fock and hybrid DFT calculations. A 'chain-of-spheres' algorithm for the Hartree–Fock exchange." *Chemical Physics* **2009**, *356*, 98–109 (origin of the COSX/RIJCOSX approach).
16. Izsák, R.; Neese, F. "An overlap fitted chain of spheres exchange method." *Journal of Chemical Physics* **2011**, *135*, 144105.

**Basis sets and extrapolation:**

17. Neese, F.; Valeev, E. F. "Revisiting the atomic natural orbital approach for basis sets: robust systematic basis sets for explicitly correlated and conventional correlated ab initio methods?" *Journal of Chemical Theory and Computation* **2011**, *7*(1), 33–43 (underlies the CBS extrapolation exponent commonly cited as "Neese and Valeev").

**Python interface / modern accessibility:**

18. Tetenberg, T.; Neugebauer, H.; Plett, C.; Santhosh, N.; Bursch, M.; Riplinger, C. "ORCA Meets Python — The ORCA Python Interface OPI." *Journal of Chemical Theory and Computation* **2026**, *22*(10), 4951–4967. https://doi.org/10.1021/acs.jctc.5c02141

**Independent assessments/benchmarks referencing ORCA's DLPNO theory:**

19. Herbert, J. M. "Assessing the domain-based local pair natural orbital (DLPNO) approximation for non-covalent interactions in sizable supramolecular complexes." *Journal of Chemical Physics* **2024**, *161*, 054114. https://doi.org/10.1063/5.0206533

---

*Note: this review draws on ORCA versions through 6.1.x as of mid-2026; consult the official ORCA Forum (orcaforum.kofo.mpg.de) and FACCTs documentation (faccts.de/orca) for the current release and licensing terms, as both evolve continuously.*

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the ORCA 	Free-for-academic-use quantum chemistry package popular for DFT, TD-DFT, and correlated wavefunction calculations, known for its efficiency and usability. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
