# deMonNano: A Comprehensive Review

*An auxiliary-density-functional-theory-derived, density-functional tight-binding (DFTB) code for large molecular, biomolecular, and nanoscale systems*

---

## 1. Overview and Identity

**deMonNano** is a density-functional-theory-based tight-binding (DFTB) simulation package developed principally by the group of Mathias Rapacioli (Laboratoire de Chimie et Physique Quantiques, Université Paul Sabatier — Toulouse III), in a lineage that traces back to the original **deMon** ("densité de Montréal") DFT program family founded by Dennis Salahub and Andreas Köster. It is the "light," tight-binding-level counterpart of the full Gaussian-basis auxiliary-density-functional-theory (ADFT) code **deMon2k**, and the two are explicitly designed to interoperate in a multiscale modeling strategy: deMon2k provides reference electronic structure and can source/verify DFTB parameters, while deMonNano performs the fast, approximate calculations needed for large systems, long trajectories, and extensive configurational sampling.

The project is hosted at `demon-nano.ups-tlse.fr` and is distributed to the academic community as source code (Fortran), typically requested from the developers. Present and former developers listed by the project include M. Rapacioli, T. Heine, L. Dontot, M. Yusef Buey, E. Posenitskiy, N. Tarrat, P. Guibourg, F. Spiegelman, F. Louisnard, C. Marti, J. Cuny, M. Morinière, C. Dubosq, S. Patchkovskii, J. Frenzel, E. Michoulier, H. Duarte, L. Zchekhov, and D. Salahub.

---

## 2. Theoretical Foundation

### 2.1 From Kohn–Sham DFT to Tight Binding

deMonNano implements the standard hierarchy of **density-functional-based tight-binding (DFTB)** approximations, obtained by a systematic Taylor expansion of the Kohn–Sham DFT total energy around a reference density (typically a superposition of neutral atomic densities), combined with a minimal valence Slater-type/numerical atomic-orbital basis and pre-tabulated (Slater–Koster) two-center integrals:

- **Non-SCC (zeroth order) DFTB** — a standard, non-self-consistent tight-binding scheme with a DFT-derived band-structure (Hamiltonian) term and a short-range repulsive pair-potential fitted against reference DFT calculations.
- **SCC-DFTB (second order)** — the self-consistent-charge extension, in which Mulliken/gross atomic charge fluctuations are relaxed self-consistently through a monopole-monopole (γ-function) Coulomb interaction, substantially improving transferability, polarity description, and applicability to polar/charged and hydrogen-bonded systems. This is the DFTB flavor most commonly used for biomolecular and hydrogen-bonded (water, PAH–water, biomolecule) systems in deMonNano.
- **Third-order extension (DFTB3-like terms)** — improved charge-dependence of the Hubbard parameters/chemical hardness, better describing highly charged or deprotonated/protonated species (relevant to biomolecular acid–base chemistry).

### 2.2 Auxiliary Density Functional Theory (ADFT) Lineage

A distinguishing conceptual root of the deMon family (and, indirectly, of deMonNano's parameterization philosophy) is **Auxiliary Density Functional Theory (ADFT)**, developed by Köster and coworkers for deMon2k. In ADFT, the exchange-correlation energy and potential are evaluated from an auxiliary density obtained by variational Coulomb-fitting rather than from the full Kohn–Sham density, replacing costly four-center electron-repulsion integrals with three-center fitting integrals. This yields near-linear-scaling SCF cycles and underlies deMon2k's ability to handle QM regions of hundreds to over a thousand atoms and MM environments of hundreds of thousands of atoms in QM/MM contexts — a capability that motivates and parallels deMonNano's own large-system orientation at the DFTB level, and which is used as a reference/benchmarking and parameter-generation tool for deMonNano's DFTB Slater–Koster parameter sets.

### 2.3 Key Methodological Extensions Implemented in deMonNano

- **Dispersion corrections**: an *a posteriori* van der Waals/dispersion correction (a Rapacioli–Spiegelman–style D-type correction, later compatible with D3-like schemes) to remedy DFTB's (and DFT's) intrinsic difficulty with London dispersion — essential for stacked aromatic systems, molecular clusters, and non-covalent biomolecular interactions.
- **Improved / iteratively refined atomic charges** for hydrogen-bonded and liquid systems (e.g., improved SCC-DFTB water parameterizations using iterative Boltzmann inversion).
- **Configuration Interaction on top of DFTB (DFTB-CI)**: an extended DFTB-CI scheme to describe charge-transfer / charge-resonance excited states in cationic molecular clusters (e.g., PAH dimers), going beyond a single-determinant DFTB ground state.
- **Time-Dependent DFTB (TD-DFTB)** for electronic excitation energies and optical spectra, and non-adiabatic / surface-hopping extensions for excited-state dynamics (in conjunction with deep-learning-accelerated potential energy surfaces in more recent work).
- **QM/MM coupling with class-1 (biomolecular) force fields**: mechanical and polarizable-electrostatic embedding schemes coupling a DFTB QM region to an MM environment described by standard biomolecular force fields, enabling hybrid DFTB/MM simulations of solvated and embedded systems.
- **Periodic boundary conditions with k-point sampling**, for extended/crystalline and surface-supported systems (graphene, graphite, metal–organic frameworks, adsorbed clusters).
- **Metadynamics and enhanced sampling** coupling for free-energy surface reconstruction (e.g., alanine dipeptide conformational free energy, compared against ADFT and higher DFTB orders).
- **Automatic differentiation of the energy expression** within self-consistent tight-binding, easing implementation and validation of new energy terms/forces.
- **Sparse, parallelized SCF algorithms**: a sparse self-consistent-field solver exploiting the locality/sparsity of the DFTB density and Hamiltonian matrices, with an MPI-parallel implementation aimed explicitly at pushing DFTB calculations toward very large atom counts (reported informally as approaching million-atom scale test systems), directly targeting the large-biomolecular/nanoscale use case.
- **Global/energy-landscape exploration tools**: coupling of the DFTB potential energy surface with basin-hopping/threshold-algorithm and tree-based stochastic search methods for structural and conformational exploration of clusters and flexible molecules.
- **Dissipative/friction dynamics** extensions within the DFTB scheme for modeling energy exchange with an environment (e.g., surface-supported systems).

---

## 3. Core Capabilities

| Category | Capability |
|---|---|
| **Electronic structure** | Non-SCC and SCC (2nd/3rd order) DFTB; spin-polarized calculations; constrained/charged species |
| **Dynamics** | Born–Oppenheimer molecular dynamics; Monte Carlo sampling; metadynamics; dissipative/Langevin-type friction dynamics |
| **Excited states** | Time-dependent DFTB (TD-DFTB); DFTB-CI for charge-resonance/charge-transfer states; non-adiabatic dynamics with surface hopping |
| **Non-covalent interactions** | Empirical dispersion (D-type/D3-compatible) corrections; refined electrostatics/charge models for hydrogen bonding |
| **Multiscale coupling** | QM(DFTB)/MM with class-1 biomolecular force fields (mechanical and polarizable embedding) |
| **Periodicity** | k-point periodic boundary conditions for crystals, surfaces, 2D materials |
| **Exploration** | Global optimization / potential-energy-landscape exploration (threshold algorithm, tree-based stochastic search) |
| **Performance** | Sparse-matrix SCF solver with MPI parallelization for large-system scaling |
| **Interoperability** | ASE (Atomic Simulation Environment) calculator interface (`ase.calculators.demonnano`); parameter/Slater–Koster file compatibility with the broader DFTB ecosystem (e.g., dftb.org repositories) |

---

## 4. Application Domains

deMonNano's published application record (drawn from its own maintained bibliography and related literature) is centered on several overlapping domains:

- **Polycyclic aromatic hydrocarbons (PAHs) and carbonaceous nanograins** — structure, fragmentation, ionization, vibrational (IR) spectroscopy, and astrophysical/interstellar chemistry (PAH clusters, cations, dissociation dynamics, "very small grain" models of the interstellar medium).
- **Molecular and atomic clusters** — water clusters (neutral, protonated, hydrated ions), PAH–water complexes, metal clusters (silver, gold, including multiply charged and phosphine-ligated gold nanoclusters), and mixed clusters (e.g., silver–hydrocarbon, Fe–PAH, Si–PAH).
- **Nanoscale and surface-supported systems** — metal nanoparticles, benzene/PAH adsorption on graphene, gold electrode–water interfaces, metal–organic frameworks (MOFs).
- **Biomolecular and water/aqueous-phase systems** — liquid water parameterization and simulation, hydrated peptide/uracil/DNA-base-related clusters, conformational free-energy landscapes of peptides (e.g., alanine dipeptide via DFTB-metadynamics), and QM/MM DFTB studies of biomolecules in explicit or force-field-represented solvent.
- **Reactivity and materials chemistry** — organic contaminant transformation pathways (e.g., ibuprofen degradation products), phthalate conformational diversity, amorphous hydrogenated carbon structure generation, chiroptical/conformational studies of small organic/organometallic molecules.
- **Spectroscopy and dynamics under excitation** — anharmonic IR spectra via finite-temperature MD, non-adiabatic electronic relaxation in PAHs (armchair vs. zigzag edge effects), photodissociation and collision-induced dissociation of clusters and cluster cations.

*(Note: while deMonNano is frequently marketed and generically referenced — including in this review's framing — as suited to "large biomolecular and nanoscale systems" owing to its DFTB/O(N)-oriented design, water/aqueous-cluster and QM/MM-biomolecular work forms a smaller, more specialized fraction of its documented publication record relative to its dominant use in PAH/carbon-nanograin astrochemistry and metal-cluster/nanoscience contexts. Users seeking DFTB specifically pre-parameterized and extensively validated for proteins/nucleic acids at scale may also wish to compare against DFTB+, which has a larger dedicated biomolecular-parameter and biomolecular-application track record.)*

---

## 5. Relationship to the Broader DFTB Ecosystem

deMonNano is one of several independent DFTB implementations descended from the same theoretical lineage (Seifert; Elstner, Frauenheim, and coworkers), alongside codes such as **DFTB+**, **CP2K**'s DFTB module, **AMBER**'s DFTB interface, and others. It shares:

- Compatibility with standard **Slater–Koster parameter sets** (`.skf` files) distributed via the community repository at dftb.org, in principle allowing reuse of parameter sets developed by other groups (e.g., `mio`, `3ob`, `matsci`, `auorg`), alongside deMonNano/Toulouse-group-specific sets (e.g., for silver/gold clusters, water, carbon systems, biomolecular "BIO" parameter type as seen in its ASE interface examples).
- An **ASE calculator interface**, allowing deMonNano to be driven from Python workflows for geometry optimization, molecular dynamics, and property calculations, using `input_arguments` such as `DFTB`, `CHARGE`, and `PARAM` (e.g., `PARAM=PTYPE=BIO` for biomolecular-oriented parameterizations).
- A close working relationship with **deMon2k** (ADFT), which is used both as a high-level reference method for benchmarking/parameterizing deMonNano's DFTB Hamiltonians and repulsive potentials, and as a complementary tool in multiscale (DFTB↔DFT↔DFT/MM) simulation protocols.

---

## 6. Distribution, Documentation, and Access

- **Official site**: `http://demon-nano.ups-tlse.fr/` — provides an "About" page, Quickstart/Download instructions, a Tutorial, a Keywords reference (input specification), and a maintained "Selected Bibliography" page.
- **User's Guide**: *"The deMon-Nano User's Guide, Installation Guide and Reference Manual"* — a PDF manual documenting installation, file I/O conventions (e.g., `ioeri.scr`, `ioscf.scr`, `iocdf.scr`, `iogrd.scr` scratch files; `RHO.bin` visualization output), and keyword/input-block syntax.
- **License/access model**: as with deMon2k, deMonNano is distributed to the academic/research community typically upon request to the developers rather than as an open, unrestricted public download; it is not part of major package managers or the DFTB+/AMBER ecosystem's default distributions.
- **Training**: deMon2k/deMonNano joint tutorials have been run under CECAM auspices, covering ADFT/ADPT theory, DFT/MM approaches, DFTB basic theory and self-consistent-charge extensions, ab initio and Born–Oppenheimer/Car–Parrinello MD, and hands-on practical sessions.

---

## 7. Summary Assessment

**Strengths**
- Strong DFT grounding (via its historical and methodological ties to deMon2k/ADFT) gives improved transferability relative to purely empirical tight-binding schemes.
- Rich, actively developed feature set for cluster science: dispersion corrections, DFTB-CI charge-resonance states, TD-DFTB, periodic boundary conditions, and non-adiabatic/excited-state dynamics — features not uniformly available across all DFTB codes.
- Purpose-built performance features (sparse SCF, MPI parallelization, automatic differentiation) explicitly aimed at large-system and long-trajectory simulation.
- Interoperates with the broader Slater–Koster parameter ecosystem and with ASE for flexible scripting.
- Direct multiscale continuity with deMon2k/ADFT for validation and QM/MM biomolecular embedding.

**Limitations / considerations**
- The bulk of the *documented and self-cited* publication record is concentrated in astrochemistry/PAH cluster physics and metal-nanocluster science rather than large-scale protein/nucleic-acid biomolecular simulation per se; "large biomolecular" capability rests more on general DFTB/QM-MM methodology (shared with, and partly validated against, deMon2k) than on an extensive record of large-protein deMonNano production runs.
- Distribution is restrictive (request-based academic license), and documentation/tutorials are comparatively sparse relative to more widely adopted community DFTB codes.
- As with all DFTB methods, accuracy is contingent on the quality and transferability of the specific Slater–Koster parameter set used for a given chemical system; users must verify parameter availability/validation for their target biomolecular or nanomaterial system.

---

## 8. Bibliography

### 8.1 Foundational DFTB and ADFT Theory (general, underpinning deMonNano's methods)

1. Porezag, D.; Frauenheim, Th.; Köhler, Th.; Seifert, G.; Kaschner, R. *Construction of tight-binding-like potentials on the basis of density-functional theory: Application to carbon.* **Phys. Rev. B** 1995, 51, 12947–12957.
2. Elstner, M.; Porezag, D.; Jungnickel, G.; Elsner, J.; Haugk, M.; Frauenheim, Th.; Suhai, S.; Seifert, G. *Self-consistent-charge density-functional tight-binding method for simulations of complex materials properties.* **Phys. Rev. B** 1998, 58, 7260–7268.
3. Frauenheim, Th.; Seifert, G.; Elstner, M.; Hajnal, Z.; Jungnickel, G.; Porezag, D.; Suhai, S.; Scholz, R. *A self-consistent charge density-functional based tight-binding method for predictive materials simulations in physics, chemistry and biology.* **Phys. Status Solidi B** 2000, 217, 41–62.
4. Elstner, M.; Frauenheim, T.; Kaxiras, E.; Seifert, G.; Suhai, S. *A self-consistent charge density-functional based tight-binding scheme for large biomolecules.* **Phys. Status Solidi B** 2000, 217, 357–376.
5. Elstner, M.; Seifert, G.; Frauenheim, Th. *Tight-binding approach to time-dependent density-functional response theory.* **Phys. Rev. B** 2001, 63, 085108.
6. Köhler, C.; Seifert, G.; Gerstmann, U.; Elstner, M.; Overhof, H.; Frauenheim, Th. *Approximate density-functional calculations of spin densities in large molecular systems and complex solids.* **Phys. Chem. Chem. Phys.** 2001, 3, 5109.
7. Yang, Y.; Yu, H.; York, D.; Cui, Q.; Elstner, M. *Extension of the self-consistent-charge density-functional tight-binding method: third-order expansion of the density functional theory total energy and introduction of a modified effective Coulomb interaction.* **J. Phys. Chem. A** 2007, 111, 10861–10873.
8. Frauenheim, Th.; Seifert, G.; Elstner, M.; Niehaus, T.; Köhler, C.; Amkreutz, M.; Sternberg, M.; Hajnal, Z.; Di Carlo, A.; Suhai, S. *Atomistic simulations of complex materials: ground-state and excited-state properties.* **J. Phys.: Condens. Matter** 2002, 14, 3015–3047.
9. Koskinen, P.; Mäkinen, V. *Density-functional tight-binding for beginners.* **Comput. Mater. Sci.** 2009, 47, 237–253.
10. Gaus, M.; Cui, Q.; Elstner, M. *DFTB3: Extension of the self-consistent-charge density-functional tight-binding method (SCC-DFTB).* **J. Chem. Theory Comput.** 2011, 7, 931–948.
11. Elstner, M.; Seifert, G. *Density functional tight binding.* **Philos. Trans. R. Soc. A** 2014, 372, 20120483.
12. Gaus, M.; Cui, Q.; Elstner, M. *Density functional tight binding: application to organic and biological molecules.* **WIREs Comput. Mol. Sci.** 2014, 4, 49–61.
13. Cui, Q.; Elstner, M. *Density functional tight binding: values of semi-empirical methods in an ab initio era.* **Phys. Chem. Chem. Phys.** 2014, 16, 14368–14377.
14. Christensen, A. S.; Kubař, T.; Cui, Q.; Elstner, M. *Semiempirical quantum mechanical methods for noncovalent interactions for chemical and biochemical applications.* **Chem. Rev.** 2016, 116, 5301–5337.

### 8.2 deMon2k / ADFT (theoretical basis and multiscale partner code)

15. Köster, A. M.; et al. *deMon2k, The deMon developers.* Cinvestav, Mexico City, program suite.
16. Salahub, D. R.; Noskov, S. Yu.; Lev, B.; Zhang, R.; Ngo, V.; Goursot, A.; Calaminici, P.; Köster, A. M.; Alvarez-Ibarra, A.; Mejía-Rodríguez, D.; Řezáč, J.; Cailliez, F.; de la Lande, A. *QM/MM Calculations with deMon2k.* **Molecules** 2015, 20, 4780–4812.
17. Flores-Moreno, R.; Köster, A. M.; et al. *(ADFT method development papers — auxiliary density functional theory, variational Coulomb fitting).*
18. *Molecular Simulations with in-deMon2k QM/MM, a Tutorial-Review.* **Molecules** 2019, 24, 1653.
19. *QM/MM with Auxiliary DFT in deMon2k*, in *Multiscale Dynamics Simulations: Nano and Nano-bio Systems in Complex Environments*, Royal Society of Chemistry (book chapter).
20. *Current status of deMon2k for the investigation of the early stages of matter irradiation by time-dependent DFT approaches.* **Eur. Phys. J. Spec. Top.** 2023.

### 8.3 deMonNano-Specific Method Development

21. Rapacioli, M.; Spiegelman, F.; Talbi, D.; Mineva, T.; Goursot, A.; Heine, T.; Seifert, G. *Correction for dispersion and Coulombic interactions in molecular clusters with density-functional derived methods: Application to polycyclic aromatic hydrocarbons.* **J. Chem. Phys.** 2009, 130, 244304.
22. Rapacioli, M.; Spiegelman, F.; Scemama, A.; Mirtschink, A. *Modeling charge resonance in cationic molecular clusters: Combining DFT-tight binding with configuration interaction.* **J. Chem. Theory Comput.** 2011, 7, 44–55.
23. Rapacioli, M.; Simon, A.; Dontot, L.; Spiegelman, F. *Extensions of DFTB to investigate molecular complexes and clusters.* **Phys. Status Solidi B** 2012, 249, 245–258.
24. Gamboa, A.; Rapacioli, M.; Spiegelman, F. *Automatic differentiation of the energy within self-consistent tight-binding methods.* **J. Chem. Theory Comput.** 2013, 9, 3900–3907.
25. Scemama, A.; Rapacioli, M.; Renon, N. *A sparse self-consistent field algorithm and its parallel implementation: Application to density-functional-based tight binding.* **J. Chem. Theory Comput.** 2014, 10, 2344–2354.
26. Dontot, L.; Suaud, N.; Rapacioli, M.; Spiegelman, F. *An extended DFTB-CI model for charge-transfer excited states in cationic molecular clusters: model studies versus ab initio calculations in small PAH clusters.* **Phys. Chem. Chem. Phys.** 2016, 18, 3545–3557.
27. Cuny, J.; Tarrat, N.; Spiegelman, F.; Huguenot, A.; Rapacioli, M. *Density-functional tight-binding approach for metal clusters, nanoparticles, surfaces and bulk: application to silver and gold.* **J. Phys.: Condens. Matter** 2018, 30, 303001.
28. Simon, A.; Rapacioli, M.; Michoulier, E.; Zheng, L.; Korchagina, K.; Cuny, J. *Contribution of the density-functional based tight binding scheme to the description of water clusters: methods, applications and extension to bulk systems.* **Mol. Simul.** 2019, 45, 249–268.
29. Cuny, J.; Korchagina, K.; Menakbi, C.; Mineva, T. *Metadynamics combined with auxiliary density functional and density functional tight-binding methods: alanine dipeptide as a case study.* **J. Mol. Model.** 2017, 23, 72.
30. Spiegelman, F.; Tarrat, N.; Cuny, J.; Dontot, L.; Posenitskiy, E.; Martí, C.; Simon, A.; Rapacioli, M. *Density-functional tight-binding: basic concepts and applications to molecules and clusters.* **Adv. Phys.: X** 2020, 5, 1710252.
31. Rapacioli, M.; Tarrat, N. *Periodic DFTB for supported clusters: implementation and application on benzene dimers deposited on graphene.* **Computation** 2022, 10, 39.
32. Yusef Buey, M.; Mineva, T.; Rapacioli, M. *Coupling density functional based tight binding with class 1 force fields in a hybrid QM/MM scheme.* **Theor. Chem. Acc.** 2022, 141, 16.
33. Cinq, N.; Simon, A.; Louisnard, F.; Cuny, J. *Accurate SCC-DFTB parametrization of liquid water with improved atomic charges and iterative Boltzmann inversion.* **J. Phys. Chem. B** 2023, 127, 7590–7601.
34. Michoulier, E.; Lemoine, D.; Spiegelman, F.; Nave, S.; Rapacioli, M. *Dissipative friction dynamics within the density-functional based tight-binding scheme.* **Eur. Phys. J. Spec. Top.** 2023, 232, 1975–1983.
35. Rapacioli, M.; Yusef Buey, M.; Spiegelman, F. *Addressing electronic and dynamical evolution of molecules and molecular clusters: DFTB simulations of energy relaxation in polycyclic aromatic hydrocarbons.* **Phys. Chem. Chem. Phys.** 2024, 26, 1499–1515.
36. Milia, V.; Rapacioli, M.; Tarrat, N.; Zanon, C.; Cortés, J. *Exploring molecular energy landscapes by coupling the DFTB potential with a tree-based stochastic algorithm: investigation of the conformational diversity of phthalates.* **J. Chem. Inf. Model.** 2024, 64, 3290–3301.
37. Rapacioli, M.; Schön, J. C.; Tarrat, N. *Exploring energy landscapes at the DFTB quantum level using the threshold algorithm: the case of the anionic metal cluster Au₂₀⁻.* **Theor. Chem. Acc.** 2021, 140, 85.

### 8.4 Representative Biomolecular / Aqueous-Phase Applications

38. Korchagina, K.; Simon, A.; Rapacioli, M.; Spiegelman, F.; L'Hermite, J.-M.; Braud, I.; Zamith, S.; Cuny, J. *Theoretical investigation of the solid–liquid phase transition in protonated water clusters.* **Phys. Chem. Chem. Phys.** 2017, 19, 27288–27298.
39. Korchagina, K.; Simon, A.; Rapacioli, M.; Spiegelman, F.; Cuny, J. *Structural characterization of sulfur-containing water clusters using a density-functional based tight binding approach.* **J. Phys. Chem. A** 2016, 120, 9089–9100.
40. Cuny, J.; Calatayud, J. C.; Ansari, N.; Hassanali, A. A.; Rapacioli, M.; Simon, A. *Simulation of liquids with the tight-binding density-functional approach and improved atomic charges.* **J. Phys. Chem. B** 2020, 124, 7421–7432.
41. Zamith, S.; Zheng, L.; Cuny, J.; L'Hermite, J.-M.; Rapacioli, M. *Collision-induced dissociation of protonated uracil water clusters probed by molecular dynamics simulations.* **Phys. Chem. Chem. Phys.** 2021, 23, 27404–27416.
42. Braud, I.; Zamith, S.; Cuny, J.; Zheng, L.; L'Hermite, J.-M. *Size-dependent proton localization in hydrated uracil clusters: a joint experimental and theoretical study.* **J. Chem. Phys.** 2019, 150, 014303.
43. Zamith, S.; Kassem, A.; L'Hermite, J.-M.; Joblin, C.; Cuny, J. *Threshold collision induced dissociation of protonated water clusters.* **J. Chem. Phys.** 2023, 159, 184301.

### 8.5 deMonNano's Own Maintained "Selected Bibliography" (application-oriented; representative recent entries, non-exhaustive)

44. Singh, D. K.; Biennier, L.; Simon, A.; Chakraborty, S.; Dartois, E.; Georges, R.; Kassi, S.; Chandrasekaran, V.; Sabbah, H.; Joblin, C.; Ranjan, S.; Suwas, S.; Jagadeesh, G.; Arunan, E. *Shock-induced evolution: tracing the fate of coronene in astrophysical environments.* **Astron. Astrophys.** 2025, 704, A345.
45. Milia, V.; Rapacioli, M.; Zanon, C.; Cortés, J.; Tarrat, N. *The ARMAGNHAC database: a ratio-based molecular analyzer and generator of numerous hydrogenated amorphous carbons.* **J. Phys. Chem. A** 2025, 129, 10358–10367.
46. Salomon, G.; Rapacioli, M.; Schön, J. C.; Tarrat, N. *Predicting energy-dependent transformation products of environmental contaminants: the case of ibuprofen.* **Physics** 2025.
47. Alauzet, C.; Spiegelman, F.; Simon, A. *Modeling silver clusters–hydrocarbon interactions: a challenge for SCC-DFTB.* **Comput. Theor. Chem.** 2024, 1239, 114744.
48. Tarrat, N.; Rapacioli, M.; Poudel, B.; Spiegelman, F. *DFTB approach to multiply charged Auₙ^q+ clusters, Part I: Structural and electronic properties.* **J. Innov. Mater. Extreme Cond.** 2024, 5, 106.
49. Tarrat, N.; Rapacioli, M.; Poudel, B.; Spiegelman, F. *DFTB approach to multiply charged Auₙ^q+ clusters, Part II: Energetics, ionization and fragmentation.* **J. Innov. Mater. Extreme Cond.** 2024, 5, 115.
50. Vuong, V. Q.; Madridejos, J. M. L.; Aradi, B.; Sumpter, B. G.; Metha, G. F.; Irle, S. *Density-functional tight-binding for phosphine-stabilized nanoscale gold clusters.* **Chem. Sci.** 2020, 11, 13113–13128.
51. Dontot, L.; Spiegelman, F.; Rapacioli, M. *Structures and energetics of neutral and cationic pyrene clusters.* **J. Phys. Chem. A** 2019, 123, 9531–9543.
52. Simon, A.; Rapacioli, M.; Mascetti, J.; Spiegelman, F. *Vibrational spectroscopy and molecular dynamics of water monomers and dimers adsorbed on polycyclic aromatic hydrocarbons.* **Phys. Chem. Chem. Phys.** 2012, 14, 6771–6786.
53. Lukose, B.; et al. *Structure and electronic structure of metal-organic frameworks within the density-functional based tight-binding method.* (arXiv/journal preprint), relevant to periodic-DFTB deMonNano applications to MOFs.

*(Note: entries 44–52 are a representative, non-exhaustive sample from deMonNano's own continuously updated "Selected Bibliography" page, which as of retrieval lists on the order of 100 application papers spanning 2009–2025, concentrated in PAH/carbon-cluster astrochemistry, metal-cluster nanoscience, and water/aqueous-cluster studies. The full current list is maintained at `demon-nano.ups-tlse.fr/pages/demon_biblio.html`.)*

---

## 9. Primary Sources Consulted

- deMonNano official website — About, Quickstart/Download, Keywords, Selected Bibliography (`demon-nano.ups-tlse.fr`)
- *The deMon-Nano User's Guide, Installation Guide and Reference Manual* (PDF, hosted on the deMonNano site)
- ASE documentation, `ase.calculators.demonnano` interface page
- CECAM deMon-2k and deMonNano tutorial workshop description
- deMon2k program documentation and QM/MM ADFT review literature (Molecules, RSC book chapter)
- Peer-reviewed literature on SCC-DFTB/DFTB3 theory (Elstner et al., Frauenheim et al., Gaus/Cui/Elstner) as cited above
- Peer-reviewed deMonNano method and application papers as cited above (Rapacioli, Spiegelman, Simon, Tarrat, Cuny, and coworkers)

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the DeMonNano 	DFT-based tight-binding code (auxiliary density functional tight-binding) for large biomolecular and nanoscale systems. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
