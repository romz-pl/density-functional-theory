# DFTB+: An Exhaustive Review

## 1. Overview

**DFTB+** is a versatile, community-developed, open-source software package for fast, approximate quantum-mechanical (density-functional-theory-based) atomistic simulations. It implements the **Density Functional Tight Binding (DFTB)** family of methods and the **extended tight-binding (xTB)** family of methods (the GFN*n*-xTB models), both of which are derived as systematic simplifications of Kohn–Sham DFT. This heritage allows DFTB+ to reach accuracies that are qualitatively closer to DFT than classical force fields, while running roughly two to three orders of magnitude faster, enabling simulations of systems with thousands of atoms and access to nanosecond-scale molecular dynamics that are inaccessible to conventional ab initio methods.

DFTB+ can be used as:
- A **standalone command-line application**, driven by a human/machine-readable input format (Human-friendly Structured Data, HSD);
- A **library**, embedded into third-party simulation codes (e.g., ASE, GROMACS, AMBER, ChemShell, i-PI, PLUMED);
- A **calculation server**, accessed via socket communication (i-PI / IPI-protocol), enabling loose coupling to external drivers for path-integral MD, replica-exchange, or advanced sampling.

## 2. History and Governance

- Development began around **2004 at the University of Paderborn**, growing out of earlier tight-binding DFT codes developed by the groups of Gotthard Seifert, Thomas Frauenheim, and collaborators (the lineage traces back to tight-binding potentials constructed from DFT in the mid-1990s).
- Development subsequently continued as a **joint project between the University of Bremen (Germany)** and the **University of Strathclyde (Glasgow, UK)**.
- Since **2017**, DFTB+ has been a fully **open-source project**, released under the **GNU Lesser General Public License (LGPL v3)**.
- The project is coordinated by **Bálint Aradi** (University of Bremen) and **Ben Hourahine** (University of Strathclyde), with contributions from a large international community (institutions including the Max Planck Institute for the Structure and Dynamics of Matter, Oak Ridge National Laboratory, University of Costa Rica, University of Warwick, University of Bonn, RIKEN, UNIST, Argonne National Laboratory, University of Luxembourg, ETH Zurich, and others).
- Source hosted at **github.com/dftbplus/dftbplus**; auxiliary repositories cover Slater-Koster parameter generation (`skprogs`), the DFTB+ recipes/tutorial collection (`recipes`), HSD parsing tools (`hsd-python`), test parameter sets (`testparams`), and developer documentation (`develguide`).

## 3. Theoretical Foundation

### 3.1 DFTB as a Kohn–Sham Expansion

DFTB and xTB are both derived by expanding the Kohn–Sham total-energy functional in powers of a charge-density fluctuation δρ(**r**) around a reference density ρ₀(**r**), the latter typically constructed as a superposition of neutral, spherically symmetric atomic densities. Truncating this Taylor expansion at different orders yields the different levels of theory:

| Level | Order of expansion | Key characteristics |
|---|---|---|
| **DFTB1** (non-SCC) | 0th order | Equivalent to a standard, non-self-consistent tight-binding scheme; Hamiltonian/overlap matrices are two-center, pre-tabulated Slater–Koster integrals. |
| **DFTB2 / SCC-DFTB** | 2nd order | Adds a self-consistent redistribution of Mulliken (atomic) charges via a Coulombic charge-fluctuation ("γ") term; captures long-range electrostatics and partially cancels self-interaction error. This is the level introduced by Elstner *et al.* (1998). |
| **DFTB3** | 3rd order | Adds diagonal third-order charge-fluctuation terms and a modified short-range γ function (γʰ) to correct problematic hydrogen–heteroatom interactions; substantially improves proton affinities, hydrogen-bonding energetics, and reaction barriers. |
| **xTB (GFN0/GFN1/GFN2-xTB)** | Related expansion, different parametrization philosophy | Share the same overall DFTB-like Hamiltonian structure but use global/element-specific (largely pair-independent) parametrization, extended-Hückel-like off-diagonal elements, D3/D4 dispersion, and (in GFN2) anisotropic multipole electrostatics, targeting geometries, vibrational frequencies, and non-covalent interactions ("GFN" = Geometry, Frequency, Noncovalent) across nearly the whole periodic table (Z = 1–86). |

Because Hamiltonian and overlap matrix elements are strictly **two-center** and depend only on the interatomic distance (and orbital angular character), they can be pre-computed and tabulated as functions of distance in **Slater–Koster (SK) files**, avoiding on-the-fly evaluation of most integrals. The total energy is completed with a short-range, pairwise **repulsive energy** term fitted to reproduce reference DFT/experimental structures and energetics, compensating for double-counting and higher-order errors.

### 3.2 What DFTB+ Adds on Top of the Core Model

Building on this DFTB/xTB core, DFTB+ implements numerous extensions that broaden accuracy and applicability:

- **Spin polarization** (collinear and non-collinear), and treatment of **spin–orbit coupling**.
- **DFTB3+U / LDA+U**-like corrections for improved treatment of localized d/f electrons in transition metals and lanthanides/actinides.
- **Long-range-corrected (range-separated) DFTB**, addressing self-interaction error and improving charge-transfer excitations and orbital-energy gaps.
- **Approximate hybrid-functional DFTB**, incorporating a fraction of exact exchange.
- **Dispersion corrections**: empirical UFF-based dispersion, Grimme's DFT-D3, self-consistently charge-dependent DFT-D4, Tkatchenko–Scheffler and many-body dispersion (MBD) methods, important for non-covalently bound and layered/molecular-crystal systems.
- **Implicit solvation models**: Generalized Born (GB), the Analytical Linearized Poisson–Boltzmann (ALPB) model, and COSMO-type continuum solvation.
- **Time-Dependent DFTB (TD-DFTB)** for excited states: both a **Casida linear-response (LR-TDDFTB)** implementation and a **real-time (Ehrenfest / RT-TDDFTB)** implementation for simulating photoabsorption spectra, nonlinear optical response, and non-adiabatic excited-state dynamics; a velocity-gauge real-time TDDFTB implementation extends this to large periodic condensed-matter systems (thousands of atoms).
- **Non-adiabatic / excited-state molecular dynamics**, including surface-hopping-type schemes.
- **Electron transport** via the **non-equilibrium Green's function (NEGF)** formalism (through the `libNEGF` library), enabling calculation of coherent and dissipative charge transport in open quantum systems (molecular electronics, nanodevices), including a self-consistent treatment with the Poisson equation for electrostatics.
- **Extended Lagrangian Born–Oppenheimer molecular dynamics (XL-BOMD)**, allowing linear-scaling, stable long-time MD without a fully converged SCC cycle at every step.
- **Path-integral molecular dynamics (PIMD)** and coupling to external drivers via socket/i-PI communication, for treating nuclear quantum effects.
- **QM/MM embedding** interfaces (e.g., with ChemShell, Amber, GROMACS) and interfaces enabling **metadynamics/enhanced sampling** via PLUMED.
- **Machine-learned and neural-network-corrected repulsive potentials** (ΔTB-style corrections) to bridge accuracy gaps versus DFT/DFT-hybrid references, an area of active development discussed in the 2025 review.
- **Electrostatic embedding and periodic boundary conditions** with full k-point sampling, enabling calculations on molecules, clusters, surfaces, and 1D/2D/3D periodic solids on equal footing.
- **Force, stress tensor, and analytic/numerical Hessian (second-derivative) evaluation**, geometry and cell optimization, transition-state search utilities, and phonon/vibrational analysis (often through interfacing with external phonon codes such as Phonopy).

### 3.3 xTB / GFN*n* Methods within DFTB+

DFTB+ implements the **GFN0-xTB, GFN1-xTB, and GFN2-xTB** Hamiltonians natively (and can interoperate with the standalone `xtb` and `tblite` libraries). These extended tight-binding methods:
- Are parametrized broadly (elements up to Z = 86, i.e., through radon) using largely **element-specific rather than element-pair-specific** parameters, improving transferability to chemical space outside the fitting set.
- GFN1-xTB uses coordination-number-dependent energy levels, polarized Gaussian-type basis functions, a double-ζ hydrogen basis, diagonal third-order charge fluctuations, an explicit halogen-bonding correction, and D3 dispersion.
- GFN2-xTB replaces the simple monopole electrostatics with **anisotropic, multipole (up to quadrupole) electrostatics** via short-range-damped interactions of atomic multipoles, and incorporates the **self-consistent, charge-dependent D4 dispersion model**, removing the need for separate halogen/hydrogen-bond corrections.
- GFN0-xTB is a simpler, extremely fast tight-binding-only variant suited to very large systems and screening applications.

Because DFTB and xTB share the same underlying Hamiltonian architecture, DFTB+ treats them within a common numerical/software infrastructure, letting users switch between DFTB parametrizations (application-specific Slater–Koster sets) and general-purpose xTB parametrizations depending on the required trade-off between accuracy, transferability, and coverage of the periodic table.

## 4. Software Architecture and Numerical Methods

- **Language and build system**: DFTB+ is written primarily in modern **Fortran**, uses **CMake** (≥ 3.16) for configuration, and supports both serial (OpenMP-threaded) and **MPI-parallelized builds** (via MPICH or Open MPI), including hybrid MPI+OpenMP execution and **GPU-accelerated eigensolvers** (via MAGMA) for large systems.
- **Linear algebra**: Uses dense and **sparse-matrix** representations depending on system size and boundary conditions; large periodic/large-scale calculations exploit **ELSI** (the Electronic Structure Infrastructure library) to interface with multiple eigensolver/density-matrix-purification back-ends (e.g., ELPA, PEXSI, NTPoly), enabling access to **linear-scaling / O(N)** density-matrix-based solvers as an alternative to conventional diagonalization for very large systems.
- **Parallel scalability**: Distributed dense linear algebra uses **ScaLAPACK** (through the `scalapackfx` Fortran wrapper library maintained within the project) for MPI-parallel diagonalization.
- **Input/output**: Calculations are configured via the **HSD (Human-friendly Structured Data)** format, a JSON/XML-like hierarchical format with a dedicated Python parsing package (`hsd-python`, BSD-licensed) maintained by the project; geometries can be read/written in `gen`, XYZ, VASP-POSCAR, and other common formats.
- **Distribution**: The preferred installation route is via the **conda-forge** channel (using `conda`/`mamba`), which provides pre-built serial and MPI (MPICH/OpenMPI) variants; source builds and pre-built binaries are also available from the GitHub releases page. Optional license-restricted external components (e.g., certain Slater-Koster test parameter sets) are fetched via a dedicated `get_opt_externals` utility.
- **External library integrations**: DFTB+ interoperates with `libNEGF` (electron transport), `libMBD` (many-body dispersion), `DFTD3`/`DFTD4` (dispersion corrections), `PLUMED` (enhanced sampling/metadynamics), `ChIMES` (many-body correction potentials), Poisson-solver libraries for open-boundary electrostatics in transport calculations, `ELSI`, and `ScaLAPACK`. It can also be driven from or drive external codes through the **Atomic Simulation Environment (ASE)** Python interface and through socket/i-PI communication.

## 5. Core Simulation Capabilities

DFTB+ supports:
- **Single-point energy, force, and (in periodic systems) stress-tensor** evaluation.
- **Geometry optimization** (including lattice/cell optimization under constant pressure or fixed cell shape) and **transition-state search**.
- **Molecular dynamics** in various ensembles (NVE, NVT via multiple thermostats, NPT via barostats), including **Born–Oppenheimer MD** and **extended-Lagrangian (XL-BOMD)** schemes for improved long-time energy conservation at reduced cost.
- **Path-integral molecular dynamics** and coupling to external MD/sampling drivers.
- **Vibrational/phonon analysis** via numerical or analytic Hessians.
- **Electronic structure analysis**: band structures, density of states (total and projected/PDOS), Mulliken/other population analyses, and partial charges.
- **Excited-state properties**: absorption spectra (TD-DFTB/LR and RT), excited-state gradients/geometry optimization, non-adiabatic dynamics.
- **Charge/electron transport** in open, non-equilibrium two- or multi-terminal device geometries via the NEGF module.
- **Free-energy methods**: interfacing to PLUMED for metadynamics and related enhanced-sampling techniques to access reaction free-energy surfaces that are otherwise inaccessible on typical DFTB MD timescales.

## 6. Application Domains

Because of its favorable cost/accuracy trade-off, DFTB+ is used extensively across:
- **Biomolecular simulation**: proteins, enzymes, DNA/RNA base-pairing and stacking, QM/MM reaction mechanism studies.
- **Materials science**: defects in graphene and other 2D materials, amorphous and nanostructured solids, metal–organic frameworks (MOFs), semiconductor and oxide surfaces.
- **Nanoscience**: fullerenes, carbon nanotubes, clusters (including ligand-stabilized metal nanoclusters), nanowires.
- **Molecular electronics and quantum transport**: current–voltage characteristics of molecular junctions, nanodevices.
- **Photophysics/photochemistry**: UV-Vis absorption spectra, plasmonic nanostructures (via RT-TDDFTB), photocatalysis, non-adiabatic excited-state dynamics.
- **Astrochemistry/astrobiology**: formation of carbon clusters and complex molecules in space, mineral/schreibersite models.
- **Reactive chemistry**: thermal decomposition of energetic materials, rare-event sampling of chemical reactions (proton transfer, bond-breaking) via metadynamics.
- **Condensed-matter and periodic systems**: band structures, phonons, and thermal transport in crystals and low-dimensional materials.

## 7. Strengths and Limitations

**Strengths**
- Two to three orders of magnitude faster than comparable ab initio DFT, enabling system sizes (thousands of atoms) and timescales (ns-scale MD) inaccessible to first-principles methods.
- Retains meaningful quantum-mechanical character (bond breaking/forming, charge transfer, polarizability) largely absent from classical force fields.
- Broad feature set spanning ground- and excited-state properties, transport, and free-energy methods within a single, actively maintained, open-source code.
- Flexible deployment (standalone, library, socket server) and strong interoperability with the broader computational-chemistry/materials ecosystem (ASE, GROMACS, AMBER, PLUMED, i-PI).
- Community governance and continuous incorporation of state-of-the-art extensions (long-range correction, hybrid DFTB, D4 dispersion, machine-learned repulsive corrections).

**Limitations**
- Accuracy is fundamentally bounded by the semi-empirical parametrization: results are only as good as the underlying Slater–Koster/xTB parameter sets, which are element- and sometimes system-class-specific; transferability outside the fitting domain is not guaranteed.
- The minimal, confined atomic-orbital basis set and (for standard DFTB1–3) monopole approximation for charge fluctuations limit accuracy for systems with significant polarization, anisotropic charge distributions, or strong static correlation, though several of these are mitigated by extensions (multipole electrostatics in GFN2-xTB, hybrid/range-separated DFTB).
- Transition-metal and heavy-element chemistry, spin-state energetics, and strongly correlated systems remain more challenging than for light main-group organic/biomolecular chemistry, though targeted extensions (DFTB3+U, spin-polarized xTB) partially address this.
- As with all semi-empirical methods, systematic benchmarking against DFT/experiment is essential before applying a given parameter set to a new chemical problem.

## 8. Availability

- **License**: GNU Lesser General Public License v3 (LGPLv3) for the core code; auxiliary documentation/recipes repositories are typically CC-BY-SA; parametrization data typically CC-BY-SA (with per-set variation).
- **Source & releases**: `github.com/dftbplus/dftbplus`
- **Installation**: conda-forge (`mamba install 'dftbplus=*=nompi_*'` or MPI variants), pre-built binaries, or source build via CMake.
- **Website / documentation**: `dftbplus.org`; parametrization files available from `dftb.org`; developer guide at `dftbplus-develguide.readthedocs.io`.

---

# Publications Related to DFTB+ Theory and Methodology

## Foundational DFTB Theory

1. Porezag, D.; Frauenheim, T.; Köhler, T.; Seifert, G.; Kaschner, R. *Construction of tight-binding-like potentials on the basis of density-functional theory: Application to carbon.* **Phys. Rev. B** 1995, 51, 12947–12957.
2. Seifert, G.; Porezag, D.; Frauenheim, T. *Calculations of molecules, clusters, and solids with a simplified LCAO-DFT-LDA scheme.* **Int. J. Quantum Chem.** 1996, 58, 185–192.
3. Elstner, M.; Porezag, D.; Jungnickel, G.; Elsner, J.; Haugk, M.; Frauenheim, T.; Suhai, S.; Seifert, G. *Self-consistent-charge density-functional tight-binding method for simulations of complex materials properties.* **Phys. Rev. B** 1998, 58, 7260–7268. (The foundational SCC-DFTB / DFTB2 paper.)
4. Frauenheim, T.; Seifert, G.; Elstner, M.; Hajnal, Z.; Jungnickel, G.; Porezag, D.; Suhai, S.; Scholz, R. *A self-consistent charge density-functional based tight-binding method for predictive materials simulations in physics, chemistry and biology.* **Phys. Status Solidi B** 2000, 217, 41–62.
5. Frauenheim, T.; Seifert, G.; Elstner, M.; Niehaus, T.; Köhler, C.; Amkreutz, M.; Sternberg, M.; Hajnal, Z.; Di Carlo, A.; Suhai, S. *Atomistic simulations of complex materials: ground-state and excited-state properties.* **J. Phys.: Condens. Matter** 2002, 14, 3015–3047.
6. Seifert, G. *Tight-binding density functional theory: an approximate Kohn–Sham DFT scheme.* **J. Phys. Chem. A** 2007, 111, 5609–5613.
7. Gaus, M.; Cui, Q.; Elstner, M. *DFTB3: Extension of the Self-Consistent-Charge Density-Functional Tight-Binding Method (SCC-DFTB).* **J. Chem. Theory Comput.** 2011, 7, 931–948.
8. Yang, Y.; Yu, H.; York, D.; Cui, Q.; Elstner, M. *Extension of the self-consistent-charge density-functional tight-binding method: third-order expansion of the density functional theory total energy and introduction of a modified effective coulomb interaction.* **J. Phys. Chem. A** 2007, 111, 10861–10873.
9. Gaus, M.; Goez, A.; Elstner, M. *Parametrization and Benchmark of DFTB3 for Organic Molecules.* **J. Chem. Theory Comput.** 2013, 9, 338–354.
10. Seifert, G.; Joswig, J.-O. *Density-functional tight binding—an approximate density-functional theory method.* **WIREs Comput. Mol. Sci.** 2012, 2, 456–465.
11. Elstner, M.; Seifert, G. *Density functional tight binding.* **Philos. Trans. R. Soc. A** 2014, 372, 20120483.
12. Koskinen, P.; Mäkinen, V. *Density-functional tight-binding for beginners.* **Comput. Mater. Sci.** 2009, 47, 237–253.
13. Spiegelman, F.; Tarrat, N.; Cuny, J.; Dontot, L.; Posenitskiy, E.; Martí, C.; Simon, A.; Rapacioli, M. *Density-functional tight-binding: basic concepts and applications to molecules and clusters.* **Adv. Phys.: X** 2020, 5, 1710252.

## The DFTB+ Software Package

14. Aradi, B.; Hourahine, B.; Frauenheim, T. *DFTB+, a Sparse Matrix-Based Implementation of the DFTB Method.* **J. Phys. Chem. A** 2007, 111, 5678–5684. (First dedicated DFTB+ software paper.)
15. Hourahine, B.; Aradi, B.; Blum, V.; Bonafé, F.; Buccheri, A.; Camacho, C.; Cevallos, C.; Deshaye, M. Y.; Dumitrică, T.; Dominguez, A.; Ehlert, S.; Elstner, M.; van der Heide, T.; Hermann, J.; Irle, S.; Kranz, J. J.; Köhler, C.; Kowalczyk, T.; Kubař, T.; Lee, I. S.; Lutsker, V.; Maurer, R. J.; Min, S. K.; Mitchell, I.; Negre, C.; Niehaus, T. A.; Niklasson, A. M. N.; Page, A. J.; Pecchia, A.; Penazzi, G.; Persson, M. P.; Řezáč, J.; Sánchez, C. G.; Sternberg, M.; Stöhr, M.; Stuckenberg, F.; Tkatchenko, A.; Yu, V. W.-z.; Frauenheim, T. *DFTB+, a software package for efficient approximate density functional theory based atomistic simulations.* **J. Chem. Phys.** 2020, 152, 124101. (Primary, most-cited reference for the current code.)
16. *Erratum: "DFTB+, a software package for efficient approximate density functional theory based atomistic simulations"* [J. Chem. Phys. 152, 124101 (2020)]. **J. Chem. Phys.** 2020.
17. Hourahine, B.; Berdakin, M.; Bich, J. A.; Bonafé, F. P.; Camacho, C.; Cui, Q.; Deshaye, M. Y.; Díaz Mirón, G.; Ehlert, S.; Elstner, M.; Frauenheim, T.; Goldman, N.; González León, R. A.; van der Heide, T.; Irle, S.; Kowalczyk, T.; Kubař, T.; Lee, I. S.; Lien-Medrano, C. R.; Maryewski, A.; Melson, T.; Min, S. K.; Niehaus, T.; Niklasson, A. M. N.; Pecchia, A.; Reuter, K.; Sánchez, C. G.; Scheurer, C.; Sentef, M. A.; Stishenko, P. V.; Vuong, V. Q.; Aradi, B. *Recent Developments in DFTB+, a Software Package for Efficient Atomistic Quantum Mechanical Simulations.* **J. Phys. Chem. A** 2025, 129, 5373–5390. (Latest comprehensive review of code capabilities.)

## Extended Tight-Binding (xTB / GFN*n*) Methods

18. Grimme, S.; Bannwarth, C.; Shushkov, P. *A Robust and Accurate Tight-Binding Quantum Chemical Method for Structures, Vibrational Frequencies, and Noncovalent Interactions of Large Molecular Systems Parametrized for All spd-Block Elements (Z = 1–86).* **J. Chem. Theory Comput.** 2017, 13, 1989–2009. (GFN1-xTB.)
19. Bannwarth, C.; Ehlert, S.; Grimme, S. *GFN2-xTB—An Accurate and Broadly Parametrized Self-Consistent Tight-Binding Quantum Chemical Method with Multipole Electrostatics and Density-Dependent Dispersion Contributions.* **J. Chem. Theory Comput.** 2019, 15, 1652–1671.
20. Spicher, S.; Grimme, S. *Robust Atomistic Modeling of Materials, Organometallic, and Biochemical Systems.* **Angew. Chem. Int. Ed.** 2020, 59, 15665–15673. (GFN-FF.)
21. Bannwarth, C.; Caldeweyher, E.; Ehlert, S.; Hansen, A.; Pracht, P.; Seibert, J.; Spicher, S.; Grimme, S. *Extended tight-binding quantum chemistry methods.* **WIREs Comput. Mol. Sci.** 2021, 11, e01493. (General GFN-xTB review.)
22. Ehlert, S.; Stahn, M.; Spicher, S.; Grimme, S. *Robust and Efficient Implicit Solvation Model for Fast Semiempirical Methods.* **J. Chem. Theory Comput.** 2021, 17, 4250–4261. (ALPB solvation, used by both xTB and DFTB+.)
23. Grimme, S. *Exploration of Chemical Compound, Conformer, and Reaction Space with Meta-Dynamics Simulations Based on Tight-Binding Quantum Chemical Calculations.* **J. Chem. Theory Comput.** 2019, 15, 2847–2862.
24. Hourahine, B. *et al.* *Overview on Building Blocks and Applications of Efficient and Robust Extended Tight Binding.* **J. Phys. Chem. A** 2025, 129, 2667–2682.

## Dispersion and Non-Covalent Interaction Corrections

25. Grimme, S.; Antony, J.; Ehrlich, S.; Krieg, H. *A consistent and accurate ab initio parametrization of density functional dispersion correction (DFT-D) for the 94 elements H–Pu.* **J. Chem. Phys.** 2010, 132, 154104. (DFT-D3.)
26. Caldeweyher, E.; Ehlert, S.; Hansen, A.; Neugebauer, H.; Spicher, S.; Bannwarth, C.; Grimme, S. *A generally applicable atomic-charge dependent London dispersion correction.* **J. Chem. Phys.** 2019, 150, 154122. (DFT-D4.)

## Excited States and Time-Dependent DFTB

27. Niehaus, T. A.; Suhai, S.; Della Sala, F.; Lugli, P.; Elstner, M.; Seifert, G.; Frauenheim, T. *Tight-binding approach to time-dependent density-functional response theory.* **Phys. Rev. B** 2001, 63, 085108.
28. Niehaus, T. A. *Approximate time-dependent density functional theory.* **J. Mol. Struct.: THEOCHEM** 2009, 914, 38–49.
29. Xu, Q.; Del Ben, M.; Okyay, M. S.; Choi, M.; Ibrahim, K. Z.; Wong, B. M. *Velocity-Gauge Real-Time Time-Dependent Density Functional Tight-Binding for Large-Scale Condensed Matter Systems.* **J. Chem. Theory Comput.** 2023 (arXiv:2308.09782).

## Electronic-Structure Infrastructure and Solvers

30. Yu, V. W.-z.; Corsetti, F.; García, A.; Huhn, W. P.; Jacquelin, M.; Jia, W.; Lange, B.; Lin, L.; Lu, J.; Mi, W.; Seifitokaldani, A.; Vázquez-Mayagoitia, Á.; Yang, C.; Yang, H.; Blum, V. *ELSI: A unified software interface for Kohn–Sham electronic structure solvers.* **Comput. Phys. Commun.** 2018, 222, 267–285.

## Application/Extension Papers Illustrating Methodology Reach

31. Aradi, B.; Niklasson, A. M. N.; Frauenheim, T. *Extended Lagrangian Density Functional Tight-Binding Molecular Dynamics for Molecules and Solids.* **J. Chem. Theory Comput.** 2015, 11, 3357–3363.
32. Pecchia, A.; Di Carlo, A. *Atomistic theory of transport in organic and inorganic nanostructures.* **Rep. Prog. Phys.** 2004, 67, 1497–1561. (NEGF transport theory as used in DFTB+/libNEGF.)
33. Kubař, T.; Elstner, M. *et al. Density functional tight binding-based free energy simulations in the DFTB+ program.* (PLUMED/free-energy interface.)

---

*Note: Publication details above were compiled from bibliographic databases, journal publisher pages, and the DFTB+ project's own documentation and citation lists; consult the original sources for complete author lists, page numbers, and DOIs where required for formal citation.*

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the DFTB+ 	Fast approximate DFT code implementing density-functional tight-binding (DFTB) and extended tight-binding (xTB) methods for large systems. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
