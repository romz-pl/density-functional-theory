# Hotbit: An Open-Source Density-Functional Tight-Binding Code in Python/ASE — Exhaustive Review

## 1. Overview

**Hotbit** is an open-source **Density-Functional Tight-Binding (DFTB)** electronic-structure code implemented as a native Python calculator for the **Atomic Simulation Environment (ASE)**. It implements the self-consistent-charge (SCC) DFTB formalism and is distributed under the **GNU General Public License (GPL)**. Its stated design goals, per the project's own documentation, are to provide:

- an open-source DFTB implementation,
- a lightweight companion to full DFT — for fast electronic-structure analysis, quick access to dynamical properties, and exploratory ("playing around") calculations,
- a compact, accessible codebase that is easy to inspect and modify (parallelization is deliberately avoided, which trades off scalability to very large systems for code transparency),
- an intuitive, scriptable user interface suited to teaching and learning realistic electronic-structure simulation, and
- a DFTB **parametrization suite**, including an interface to the **libxc** exchange-correlation library.

| Attribute | Detail |
|---|---|
| **Method** | Self-consistent-charge density-functional tight-binding (SCC-DFTB) |
| **Language** | Python (with performance-critical routines historically supplemented by compiled extensions) |
| **Integration** | Native ASE `Calculator` — works with ASE `Atoms`, dynamics, optimizers, NEB, etc. |
| **License** | GNU GPL |
| **Primary authors** | Pekka Koskinen and Ville Mäkinen (University of Jyväskylä, Finland), with a Fortran-based tight-binding predecessor by Michael Moseler (Fraunhofer IWM, Freiburg) |
| **Other developers** | Lars Pastewka (Fraunhofer IWM), and community contributors |
| **Repository** | github.com/pekkosk/hotbit (originally hosted on a University of Jyväskylä Trac server; moved to GitHub in December 2015) |
| **Units (internal)** | Hartree and Bohr |
| **Units (ASE interface)** | eV and Ångström |
| **Project history** | Under version control from September 2008; source opened under GPL in April 2009 |

## 2. Theoretical Foundation

Hotbit implements the standard **SCC-DFTB** approach: total energy and Hamiltonian matrix elements are derived from a second-order expansion of the DFT total-energy functional in fluctuations of the electron density around a reference density built from superposed neutral-atom densities. This gives an energy expression with three principal contributions:

1. **Band-structure (non-self-consistent tight-binding) energy** — the sum of occupied Kohn–Sham-like eigenvalues obtained by diagonalizing a Hamiltonian built from two-center Slater–Koster integrals computed once and tabulated as a function of interatomic distance.
2. **Second-order (self-consistent-charge) correction** — a Coulombic term coupling fluctuations of atomic (Mulliken) partial charges via a distance- and chemical-hardness-dependent function ("gamma function") that interpolates between the point-charge limit at long range and an on-site Hubbard-like limit at short range; this term is iterated to self-consistency, analogous to solving a simplified Kohn–Sham problem.
3. **Short-range repulsive energy** — a pairwise, empirically or DFT-fitted potential that absorbs double-counting corrections, exchange-correlation contributions, and core-core repulsion, fitted separately for each element pair.

This framework traces back to the earlier tight-binding-from-DFT works of Foulkes and Haydock, Seifert and co-workers, and Porezag *et al.*, and was formalized as SCC-DFTB by Elstner *et al.* (1998). **Hotbit's own theoretical exposition and its self-contained derivation of this formalism, together with practical recipes for parametrization**, is given in the paper *"Density-functional tight-binding for beginners"* by Koskinen and Mäkinen (2009) — the code's canonical citation (see §6). Distinctive aspects of Hotbit's implementation and parametrization workflow include:

- **Confined pseudo-atomic orbital basis sets** generated numerically by solving the atomic Kohn–Sham problem for a compressed/confined potential (Hotbit includes its own all-electron atomic DFT solver for this purpose, exposed via the `KSAllElectron` class).
- **Systematic, largely automated fitting of the short-range repulsive potentials** to reference DFT energy curves (dimers, chains, bulk equations of state, etc.), implemented via a dedicated `RepulsiveFitting` class in the `hotbit.parametrization` sub-package.
- **A parametrization suite with a libxc interface**, allowing Slater–Koster tables to be generated self-consistently with a chosen LDA or GGA exchange-correlation functional from the libxc library (meta-GGA and hybrid functionals are explicitly noted as unsupported/pending).
- Support for **charge self-consistency and long-range electrostatics**, including options for how the Coulomb interaction between fluctuating charges is treated (e.g., via the "gamma-function"/Ewald-type treatment appropriate to molecules, clusters, wires, slabs, and periodic bulk systems).
- A dedicated treatment of **non-standard periodic boundary conditions**, developed and published by Hotbit-affiliated authors (Kit, Pastewka, Koskinen, 2011) for simulating distorted/curved/twisted materials such as nanoribbons and nanotubes without explicit large supercells.

## 3. Core Features and Capabilities

- **Ground-state total-energy and force calculations** for molecules, clusters, wires, surfaces, and periodic bulk systems, entirely through the standard ASE `Calculator` interface (`get_potential_energy`, `get_forces`, `get_stress`, etc.).
- **Self-consistent Mulliken-charge treatment** of electrostatics, with configurable convergence mixers (e.g., an Anderson mixer) and convergence thresholds.
- **k-point sampling** for periodic systems (bulk crystals, slabs, nanotubes) via the ASE-style interface.
- **Spin-unpolarized ground-state DFTB** — the developers explicitly restrict scope to spin-paired SCC-DFTB, reasoning that the approximations already inherent in the repulsive-potential tight-binding scheme make further refinements (spin-polarization, going beyond ground-state theory) less valuable than moving to full DFT when higher accuracy is required.
- **Access to derived electronic-structure quantities**: eigenvalues/eigenvectors, density of states, Mulliken populations/charges, and related analysis outputs, exposed through the calculator's internal objects (e.g., aliases such as `calc.rep` for the repulsion-potential object, accessible via the interactive "hotbit shell").
- **A command-line/interactive "hotbit shell"** (`hotbit -h`, `hotbit -l <class>`, `hotbit -e`, etc.) providing built-in help, source-code listing of internal classes, and example scripts that can be copied directly into a working directory — intended to support the code's teaching/learning goal.
- **A parametrization workflow** for generating new Slater–Koster tables and fitting repulsive potentials for arbitrary element pairs, rather than being restricted to a fixed, pre-supplied parameter set.
- **Full ASE ecosystem interoperability**: because Hotbit is implemented as a native ASE calculator, it can be dropped into any ASE-based workflow — structure relaxation (`ase.optimize`), molecular dynamics, nudged elastic band (NEB) transition-state searches, phonon/vibrational analysis, thermochemistry, genetic-algorithm structure search, and so on — using exactly the same `atoms.calc = Hotbit(...)` pattern as other ASE calculators (e.g., DFTB+, GPAW, VASP interfaces).

## 4. Typical Usage Pattern

Hotbit is used the same way as any ASE calculator: an ASE `Atoms` object is instantiated (or read from file), a `Hotbit` calculator object is created with the desired parameters (Slater–Koster tables, charge-mixing scheme, convergence criteria, k-points, etc.) and attached to the atoms, after which standard ASE methods retrieve energies, forces, and derived properties. The parametrization submodule (`hotbit.parametrization`) is used separately, in atomic units, to generate or refine the underlying Slater–Koster and repulsive-potential tables before they are used in eV/Å production calculations.

## 5. Installation, Maintenance, and Ecosystem Status

- Hotbit installs via the classic Python `setup.py install` workflow; a helper script (`env_exports`) sets the environment variables (parameter-file paths, etc.) needed at runtime.
- Building the full parametrization suite with GGA/LDA functionals requires a separate installation of **libxc** and its Python bindings.
- Slater–Koster parameter sets compatible with Hotbit are also available from the general **dftb.org** repository, though a separate license agreement with those maintainers is required to use them; Hotbit also ships some in-house ("unofficial") first-trial parametrizations explicitly marked as unsuitable for publication-quality results.
- **Development status**: Hotbit is a comparatively small, academically maintained project (on the order of tens of GitHub stars/forks) rather than an actively, continuously developed large-scale HPC code. Its explicit design decision to avoid parallelization limits its applicability to large systems (in contrast to codes such as DFTB+ or the Fortran/C++ DFTB implementations that target massively parallel execution). It is best understood as a research/teaching-oriented DFTB engine rather than a production HPC tool.
- **Relation to other DFTB software in the ASE universe**: Hotbit should be distinguished from the separately maintained **DFTB+** program, which is also usable as an ASE calculator (`ase.calculators.dftb.Dftb`) but is an independent, external Fortran/C++ code requiring its own Slater–Koster `.skf` files; Hotbit, by contrast, is a *native* Python/ASE implementation with its own internal parametrization machinery. Related, more actively maintained parametrization tooling in the same academic lineage includes **Hotcent** (for generating DFTB parameters) and **Tango** (for repulsive-potential fitting), used together with DFTB+ in more recent ASE-DFTB tutorials.

## 6. Strengths and Limitations

**Strengths**
- Fully transparent, pure-Python(+extensions) implementation — the entire method, from atomic-orbital generation to repulsive-potential fitting, is inspectable and modifiable by the user, which is valuable for teaching and for methodological experimentation.
- Seamless integration with the whole ASE ecosystem (optimizers, dynamics, NEB, I/O formats, other calculators for comparison).
- Includes its own parametrization suite, so it is not dependent on externally licensed Slater–Koster tables (though such tables can also be used where license terms permit).
- libxc interface allows exploration of different LDA/GGA functionals in the underlying parametrization.
- Good pedagogical documentation, anchored by a widely cited, self-contained "beginners" theory paper.

**Limitations**
- No built-in parallelization, limiting practical system sizes compared to production DFTB codes.
- Restricted to spin-unpolarized, ground-state DFTB — no time-dependent DFTB/excited-state treatment, and no spin-polarized treatment, is provided within Hotbit itself (these are covered by the broader DFTB literature and by other codes, e.g. DFTB+).
- Meta-GGA and hybrid functionals are not supported in the libxc-based parametrization pipeline (only LDA/GGA), per the project's own documentation.
- The project has a comparatively small maintenance community and has not seen the same level of continuous development as DFTB+ or other production DFTB packages; some parameter sets shipped with the code are explicitly flagged as not publication-quality.
- Documentation is spread between a GitHub Wiki, a project paper, and an interactive shell-based help system rather than a single comprehensive manual, which can raise the learning curve for newcomers unfamiliar with the ASE ecosystem.

## 7. Summary

Hotbit occupies a specific niche in the DFTB software landscape: a **small, transparent, education- and research-oriented, native-Python/ASE SCC-DFTB implementation** with its own integrated parametrization suite (including a libxc interface), created primarily by Pekka Koskinen and Ville Mäkinen at the University of Jyväskylä, building on a Fortran predecessor from Fraunhofer IWM. It trades raw scalability and breadth of method (no spin-polarization, no TD-DFTB, no built-in parallelism) for code clarity, scriptability, and ease of modification, making it well suited to teaching electronic-structure concepts, rapid prototyping, and the development/testing of new DFTB parametrizations — while large-scale production DFTB work is more commonly carried out with codes such as DFTB+.

---

# Publications Related to Hotbit's Theory and Development

## A. The Canonical Hotbit / SCC-DFTB Theory Papers

1. **Koskinen, P.; Mäkinen, V.** *"Density-functional tight-binding for beginners."* **Computational Materials Science** 47(1), 237–253 (2009). DOI: 10.1016/j.commatsci.2009.07.013. (Also available as arXiv:0910.5861.) — The primary citation for Hotbit; a self-contained pedagogical derivation of SCC-DFTB from DFT and a practical guide to parametrization (pseudo-atomic orbitals, matrix elements, repulsive-potential fitting), explicitly introducing Hotbit as the accompanying open-source software.

2. **Koskinen, P.; Häkkinen, H.; Seifert, G.; Sanna, S.; Frauenheim, Th.; Moseler, M.** *"Density-functional based tight-binding study of small gold clusters."* **New Journal of Physics** 8, 9 (2006). DOI: 10.1088/1367-2630/8/1/009. — Early application/validation paper by core Hotbit developers, assessing SCC-DFTB against DFT for gold-cluster geometries and electronic structure.

## B. Foundational SCC-DFTB and Tight-Binding-from-DFT Theory (cited by the Hotbit project as its theoretical lineage)

3. **Seifert, G.; Eschrig, H.; Bieger, W.** *"Eine approximative Variante des LCAO-Xα-Verfahrens."* **Zeitschrift für Physikalische Chemie** 267, 529 (1989).
4. **Foulkes, W. M. C.; Haydock, R.** *"Tight-binding models and density-functional theory."* **Physical Review B** 39, 12520 (1989).
5. **Porezag, D.; Frauenheim, Th.; Köhler, Th.; Seifert, G.; Kaschner, R.** *"Construction of tight-binding-like potentials on the basis of density-functional theory: Application to carbon."* **Physical Review B** 51, 12947 (1995).
6. **Elstner, M.; Porezag, D.; Jungnickel, G.; Elsner, J.; Haugk, M.; Frauenheim, Th.; Suhai, S.; Seifert, G.** *"Self-consistent-charge density-functional tight-binding method for simulations of complex materials properties."* **Physical Review B** 58, 7260–7268 (1998). — The foundational formulation of SCC-DFTB that Hotbit implements.
7. **Niehaus, T. A.** *"Entwicklung approximativer Methoden in der zeitabhängigen Dichtefunktionaltheorie."* PhD thesis, Universität Paderborn (2001).
8. **Niehaus, T. A.; Suhai, S.; Della Sala, F.; Lugli, P.; Elstner, M.; Seifert, G.; Frauenheim, Th.** *"Tight-binding approach to time-dependent density-functional response theory."* **Physical Review B** 63, 085108 (2001).
9. **Todorov, T. N.** *"Time-dependent tight binding."* **Journal of Physics: Condensed Matter** 13, 10125 (2001).
10. **Niehaus, T. A.; Heringer, D.; Torralva, B.; Frauenheim, Th.** *"Importance of electronic self-consistency in the TDDFT based treatment of nonadiabatic molecular dynamics."* **European Physical Journal D** 35, 467 (2005).

## C. Physics Papers Employing Hotbit (selected list, as compiled by the project's own Wiki)

11. Malola, S.; Häkkinen, H.; Koskinen, P. *"Gold in graphene: In-plane adsorption and diffusion."* **Applied Physics Letters** 94, 043106 (2009).
12. Koskinen, P.; Malola, S.; Häkkinen, H. *"Self-passivating edge reconstructions of graphene."* **Physical Review Letters** 101, 115502 (2008).
13. Malola, S.; Häkkinen, H.; Koskinen, P. *"Raman spectra of single-walled carbon nanotubes with vacancies."* **Physical Review B** 77, 155412 (2008).
14. Malola, S.; Häkkinen, H.; Koskinen, P. *"Comparison of Raman spectra and vibrational density of states between graphene nanoribbons."* **European Physical Journal D** 52, 71 (2009).
15. Rakyta, P.; Kormányos, A.; Cserti, J.; Koskinen, P. *"Exploring the graphene edges with coherent electron focusing."* **Physical Review B** 81, 115411 (2010).
16. Lin, X.; Nilius, N.; Sterrer, M.; Koskinen, P.; Häkkinen, H.; Freund, H.-J. *"Characterizing the periphery atoms of Au islands on MgO thin films."* **Physical Review B** (2010).
17. Malola, S.; Häkkinen, H.; Koskinen, P. *"Structural, chemical and dynamical trends in graphene grain boundaries."* **Physical Review B** 81, 165447 (2010).
18. Koskinen, P.; Kit, O. O. *"Efficient approach for simulating distorted materials."* **Physical Review Letters** 105, 106401 (2010).
19. Koskinen, P. *"Electronic and optical properties in carbon nanotubes under pure bending."* **Physical Review B** 81, 193409 (2010).
20. Koskinen, P.; Kit, O. O. *"Approximate modeling of spherical membranes."* **Physical Review B** 81, 235420 (2010).
21. Kit, O. O.; Pastewka, L.; Koskinen, P. *"Revised periodic boundary conditions: Fundamentals, electrostatics, and the tight-binding approximation."* **Physical Review B** 84, 155431 (2011).
22. Koskinen, P. *"Electromechanics of twisted graphene nanoribbons."* **Applied Physics Letters** 99, 013105 (2011).
23. Kit, O. O.; Tallinen, T.; Mahadevan, L.; Timonen, J.; Koskinen, P. *"Twisting graphene nanoribbons into carbon nanotubes."* **Physical Review B** 85, 085428 (2012).
24. Ramasubramaniam, A.; Koskinen, P.; Kit, O. O.; Shenoy, V. B. *"Edge-stress-induced spontaneous twisting of graphene nanoribbons."* **Journal of Applied Physics** 111, 054302 (2012).
25. Koskinen, P. *"Graphene nanoribbons subject to gentle bends."* **Physical Review B** 85, 205429 (2012).
26. Mäkinen, V.; Koskinen, P.; Häkkinen, H. *"Modeling thiolate-protected gold clusters with density-functional tight-binding."* **European Physical Journal D** 67, 38 (2012).
27. Korhonen, T.; Koskinen, P. *"Electronic structure trends of Möbius graphene nanoribbons from minimal-cell simulations."* (2012).
28. Stiehler, C.; Pan, Y.; Schneider, W.-D.; Koskinen, P.; Häkkinen, H.; Nilius, N.; Freund, H.-J. *"Electron quantization in arbitrarily shaped Au islands on MgO thin films."* **Physical Review B** 88, 115415 (2013).
29. Koskinen, P.; Fampiou, I.; Ramasubramaniam, A. *"Density-functional tight-binding simulations of curvature-controlled valley polarization and band-gap tuning in bilayer MoS₂."* **Physical Review Letters** 112, 186802 (2014).
30. Korhonen, T.; Koskinen, P. *"Electromechanics of graphene spirals."* **AIP Advances** 5, 127125 (2015).
31. Koskinen, P.; Korhonen, T. *"Plenty of motion at the bottom: Atomically thin liquid gold membrane."* **Nanoscale** 7, 10140 (2015).

## D. Supporting Infrastructure Referenced by Hotbit

32. **Larsen, A. H.; Mortensen, J. J.; Blomqvist, J.; Castelli, I. E.; Christensen, R.; Dułak, M.; Friis, J.; Groves, M. N.; Hammer, B.; Hargus, C.; Hermes, E. D.; Jennings, P. C.; Jensen, P. B.; Kermode, J.; Kitchin, J. R.; Kolsbjerg, E. L.; Kubal, J.; Kaasbjerg, K.; Lysgaard, S.; Maronsson, J. B.; Maxson, T.; Olsen, T.; Pastewka, L.; Peterson, A.; Rostgaard, C.; Schiøtz, J.; Schütt, O.; Strange, M.; Thygesen, K. S.; Vegge, T.; Vilhelmsen, L.; Walter, M.; Zeng, Z.; Jacobsen, K. W.** *"The Atomic Simulation Environment — a Python library for working with atoms."* **Journal of Physics: Condensed Matter** 29, 273002 (2017). — The ASE framework Hotbit integrates with as a native calculator.

---

*Note: References in Sections B and C reflect the theory and application-literature lists maintained on the Hotbit project's own GitHub Wiki ("About hotbit" page), supplemented by bibliographic verification of the primary Hotbit citation (Section A) via ScienceDirect, arXiv, and Semantic Scholar. Some entries in Section C (particularly volume/page details for items 16 and 27) were only partially specified in the source Wiki listing.*


---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the hotbit     Open-source density-functional tight-binding code implemented in Python, integrated with ASE. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
