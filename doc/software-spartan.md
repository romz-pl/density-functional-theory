# Spartan (Wavefunction, Inc.) — Exhaustive Software Review

## 1. Overview

Spartan is a commercial molecular modeling and computational chemistry package developed by **Wavefunction, Inc.** (Irvine, California), first released in 1991 and now in its latest major iteration, **Spartan'26**. It is built around a strong graphical user interface for building, running, and visualizing calculations, and it spans three classes of computational model:

- **Molecular mechanics** (force fields)
- **Quantum mechanics** (semi-empirical, Hartree–Fock, DFT, post-HF/correlated methods, and composite thermochemical recipes)
- **Machine learning models** (introduced in recent releases — corrected MMFF conformational energies, neural-network-predicted DFT geometries, and ML-based NMR protocols)

Since the mid-2000s, all quantum-mechanical methods beyond molecular mechanics and semi-empirical models have been powered by the **Q-Chem** computational engine, which Spartan wraps in its own interface, task/model dialogues, and visualization layer. This positions Spartan less as a "from-scratch" quantum chemistry code and more as an integrated modeling *environment*: build → choose task/method → compute → visualize/analyze, all without leaving the GUI or writing input decks by hand.

Spartan ships in several editions (**Spartan Parallel Suite/Spartan'26** for research/commercial/government use, and **Spartan Student Edition** for undergraduate teaching), plus companion products (**iSpartan** for iOS, and **Odyssey** for molecular dynamics teaching).

---

## 2. Availability, Platforms, and Licensing

| Aspect | Details |
|---|---|
| Platforms | Windows, macOS, Linux (x86-64) |
| Current release | Spartan'26 (previous: Spartan'24, '20, '18, '16, '14, '10, '08, '06, '04, '02, and earlier Unix/Mac/PC lineages back to 1991) |
| License type | Proprietary, commercial (single-seat and network/floating licenses) |
| Core scaling | Standard licenses configured for up to 16 CPU cores for select tasks; higher tiers available for >16-core use |
| Pricing tiers | Differentiated **Commercial**, **Government**, and **Academic** pricing, with academic pricing substantially discounted (illustrative Spartan'26 >16-core single license: ~$5,400 commercial / $3,600 government / $1,800 academic; network licenses roughly 50% higher) |
| Maintenance | Annual/3-year maintenance plans provide license transfers, priority support, and minor/major version updates; new purchases include one year free |
| Remote/server use | Can act as a **Computational Server**, accepting jobs submitted from other Spartan installations or from the iSpartan mobile app |
| Related products | iSpartan (iPad/iPhone/iPod Touch companion), Spartan Student Edition (teaching-focused subset), Odyssey (MD teaching tool) |

Because pricing is quoted per license tier and changes between releases, prospective purchasers should confirm current figures directly with Wavefunction rather than relying on any single cached quote.

---

## 3. Molecule Building and User Interface

Spartan's defining strength, repeatedly cited across its user base, is ease of molecular construction and calculation setup relative to script/input-file-driven quantum chemistry packages:

- **3D builders**: organic, inorganic, peptide/protein, nucleotide/oligonucleotide, and substituent builders, allowing construction of small organics through complex biomolecules and weakly bound aggregates/complexes.
- **2D sketching**: an integrated 2D sketch palette; on Windows, direct interoperability with ChemDraw (version 9.0+) is supported for 2D-to-3D conversion.
- **Dual 2D/3D display**: a dual-display window allows simultaneous 2D structure and 3D model views.
- **Calculations dialogue**: a single dialog specifies task (see §5) and computational model/method (see §4), removing the need to hand-write text input decks for standard jobs.
- **Real-time job monitor**: live visualization of geometry/transition-state optimizations as they progress.
- **Touch support**: touch-screen operation supported on compatible Windows devices.
- **Data analysis**: an internal spreadsheet with linear regression supports QSAR-style analysis of computed and experimental properties across sets of molecules.
- **Compact native files**: `.spartan` project files use compression and reduced-precision graphical data to keep file sizes manageable.

This design targets **broad usability**: bench chemists, students, and non-specialists can set up credible DFT/HF/semi-empirical jobs without deep familiarity with quantum chemistry input syntax, while experienced computational chemists retain access to a wide method/basis-set space and to text output for verification.

---

## 4. Computational Methods

### 4.1 Molecular Mechanics (Force Fields)
- **MMFF** (Merck Molecular Force Field, MMFF94), including an extended parameterization, plus a validation test suite
- **SYBYL** force field
- Applicable to systems of several thousand atoms; the only Spartan technique suited to biopolymers directly
- Used for rapid geometry pre-optimization, conformer generation, and (in recent releases) as the seed geometry for ML-accelerated DFT geometry prediction

### 4.2 Semi-Empirical Methods
- **MNDO** / MNDO(d)
- **AM1** (Austin Model 1)
- **PM3**, with parameter extensions to additional main-group and heavy elements
- **RM1** (Recife Model 1 — AM1 reparameterization for H, C, N, O, P, S, F, Cl, Br, I)
- **PM6**

### 4.3 Hartree–Fock and DFT
- **Hartree–Fock (HF/SCF)**, with implicit solvation via **SM8**
- **DFT**, also available with SM8 implicit solvation, spanning:
  - *Exchange functionals*: HF (exact exchange), Slater–Dirac (LDA), Becke88 (B88), Gill96, GG99, B(EDF1), PW91
  - *Correlation functionals*: VWN, Perdew86, LYP, PW91-correlation, Perdew–Zunger self-interaction correction, PBE
  - *Hybrid/composite and meta-GGA functionals*: **B3LYP**, **EDF1**, **EDF2** (parameterized specifically for vibrational frequency prediction), the **M06 suite** (M05, M06, M06-L, M06-2X, etc.), and the long-range-corrected, dispersion-corrected **ωB97X-D** — the latter is the workhorse functional behind Wavefunction's own precomputed spectral database (see §7)
  - *Basis sets*: the standard Pople family (3-21G, 6-31G\*, 6-311+G\*\*, etc.) is used throughout, consistent with the package's Pople-group lineage (co-founder Warren Hehre was a long-time collaborator of John Pople)

### 4.4 Post-Hartree–Fock / Correlated Methods
- **Møller–Plesset perturbation theory**: MP2 (including resolution-of-identity RI-MP2 and analytical-gradient variants), MP4
- **Coupled cluster methods**: CCD, CCSD-type approaches, quadratic CI-based variants, and related size-consistent doubles models
- **Configuration interaction**: CIS, CIS(D), QCISD, QCISD(T), RI-CIS(D)
- **Excited-state methods**: **TDDFT** (time-dependent DFT) and CI-singles-based approaches for vertical excitation energies, oscillator strengths, and excited-state geometries/gradients

### 4.5 Thermochemical Composite Methods
- **G3(MP2)** (Gaussian-3 reduced-order variant)
- **T1** (Wavefunction's own efficient heat-of-formation recipe, developed in-house and published in *J. Phys. Chem. A*)

### 4.6 Machine Learning Models (New in Recent Releases)
- **Corrected MMFF**: a neural-network correction layered on MMFF molecular mechanics to improve conformational energy differences in flexible organic molecules, approaching DFT-quality relative conformer energies at force-field cost
- **Est. Density Functional (ML geometry prediction)**: a neural network trained to predict ωB97X-D/6-31G\* equilibrium geometries starting from an MMFF geometry, bypassing iterative DFT optimization for applicable systems
- **MLHF-NMR protocol**: a machine-learning-assisted workflow (introduced with Spartan'26) for rapid, accurate ¹³C NMR chemical shifts of flexible organic/natural-product molecules, reported to run in under ~90 seconds for typical structures

This ML layer is Spartan's most significant recent architectural addition — it does not replace DFT/HF but sits alongside it as a faster surrogate for specific, well-validated tasks (conformer energetics, geometries, and NMR shift prediction), explicitly targeted at accuracy/speed tradeoffs relevant to bench and pharmaceutical chemists working with natural products and drug-like molecules.

---

## 5. Computational Tasks

Spartan structures nearly all work around a fixed menu of tasks, applicable across the model classes above:

- **Energy** — single-point energy and properties at a given geometry (computes the wavefunction for QM models)
- **Equilibrium geometry** — local energy minimization
- **Transition-state geometry** — first-order saddle-point location
- **Equilibrium conformer** — lowest-energy conformer search, often as a pre-step before a higher-level QM geometry optimization
- **Conformer distribution** — Boltzmann-weighted sampling of low-energy conformers for property averaging
- **Conformer library** — exhaustive conformer enumeration for similarity/pharmacophore work
- **Energy profile** — relaxed scans along a user-defined coordinate (e.g., reaction coordinate, dihedral scan)
- **Similarity analysis** — structural or pharmacophoric similarity scoring between molecules/conformers, including hydrogen-bond donor/acceptor, ionizable, hydrophobe, and aromatic feature matching

---

## 6. Visualization and Property Prediction

Graphical surfaces and property maps are a hallmark Spartan feature, widely used in chemistry education as well as research rationalization of reactivity/selectivity:

**Surfaces**
- Molecular orbitals (HOMO, LUMO, and others)
- Electron density (isodensity surfaces, size-defining)
- Spin density (radical/open-shell reactivity indicator)
- Van der Waals radius surfaces
- Solvent-accessible surface area
- Electrostatic potential surfaces

**Composite property maps**
- Electrostatic potential map (electrophilic-attack indicator)
- Local ionization potential map (relative ease of local ionization)
- LUMO map (nucleophilic-attack indicator)

**Spectra**
- IR/FT-IR and Raman spectra (with scaled vibrational frequencies)
- NMR: ¹H and ¹³C chemical shifts (empirically corrected) and coupling constants; ¹³C DEPT; 2D H–H and C–H correlation spectra
- UV/visible spectra (via TDDFT/CIS-type excited-state calculations)
- Import of experimental spectra for direct comparison: IR/UV-vis in JCAMP-DX (.dx) format, NMR in Chemical Markup Language (.cml)

---

## 7. Integrated Databases

Spartan distinguishes itself from many computational chemistry packages by bundling large precomputed and experimental reference databases directly into the interface:

| Database | Content | Approx. size |
|---|---|---|
| **Spartan Spectra & Properties Database (SSPD)** | Structures, energies, NMR/IR spectra, and wavefunctions computed at ωB97X-D/6-31G\* | >305,000–315,000 molecules |
| **Spartan Molecular Database (SMD)** | Structures/energies at HF (3-21G, 6-31G\*, 6-311+G\*\*), B3LYP (6-31G\*, 6-311+G\*\*), EDF1 (6-31G\*), MP2 (6-31G\*, 6-311+G\*\*), G3(MP2), and T1 | ~100,000 molecules |
| **NMRShiftDB** | Open-source experimental ¹H/¹³C shift database | — |
| **Cambridge Structural Database (CSD)** | Experimental small-molecule crystal structures | ~600,000 entries |
| **NIST Chemistry WebBook** | Experimental IR and UV-vis spectra | — |
| **Protein Data Bank (PDB)** | Biological macromolecular structures, on-line access | >95,000 entries |

These allow instant retrieval of precomputed structure/property/spectral data for known molecules, and on-the-fly regeneration of graphical surfaces from stored wavefunctions, without re-running a calculation.

---

## 8. Strengths (Broad-Usability Focus)

- **Low barrier to entry**: GUI-driven job setup removes the need to author quantum chemistry input files, making DFT/HF/semi-empirical methods accessible to organic/medicinal/industrial chemists without formal computational training.
- **Method breadth in one package**: molecular mechanics through post-HF correlated methods and composite thermochemistry are unified under one task/model dialogue, letting users escalate accuracy without switching software.
- **Strong visualization**: orbital/density/property maps are directly generated and manipulated in 3D, useful both for research rationalization and classroom instruction.
- **Precomputed data at scale**: the SSPD/SMD databases let routine questions ("what does the DFT-level HOMO of toluene look like?") be answered instantly rather than requiring a fresh calculation.
- **Practical spectral workflows**: empirically corrected NMR shift protocols and the newer ML-based MLHF-NMR workflow are aimed squarely at a common bench task (structure verification/elucidation via calculated vs. experimental NMR), with published accuracy benchmarks (~2 ppm for rigid molecules, ~3.5–4 ppm for flexible ones under earlier protocols).
- **Cross-platform and mobile-extendable**: native support on Windows/macOS/Linux, plus iSpartan for job submission from iOS devices via the Computational Server feature.
- **Backed by a modern QM engine**: since QM methods beyond semi-empirical route through Q-Chem, Spartan inherits a actively developed, published computational engine rather than a bespoke, less-vetted backend.

## 9. Limitations and Considerations

- **Commercial/proprietary licensing**: unlike open packages (e.g., ORCA, Psi4, NWChem), Spartan requires paid licenses (with academic discounts) and lacks source-level extensibility.
- **Ceiling on system size / scaling**: designed and marketed around desktop/workstation-scale problems (up to ~16 cores per standard license); it is not positioned as an HPC/cluster-scale package for very large periodic or biomolecular systems in the way that dedicated HPC codes are.
- **Interface-mediated method access**: the ease-of-use GUI trades off some of the fine-grained control (custom SCF convergence schemes, exotic method combinations, developer-level algorithmic options) available in script/keyword-driven packages.
- **No periodic/solid-state DFT**: Spartan's DFT implementation targets molecular (finite, non-periodic) systems; it is not a plane-wave/periodic solid-state DFT code.
- **ML models are narrow-scope surrogates**: the neural-network features (corrected MMFF, ML geometry, MLHF-NMR) are validated for specific tasks (conformational energetics, ωB97X-D geometry approximation, ¹³C shifts of flexible natural products) rather than general-purpose ML potentials across arbitrary chemistries.

---

## 10. Bottom-Line Assessment

Spartan occupies a distinctive niche: it is not the most powerful or most extensible quantum chemistry engine available, but it is arguably the most **usable** integrated environment for going from a sketched or built molecule to a DFT-, HF-, or semi-empirical-level structure, spectrum, or property map, with minimal setup friction. Its combination of (a) a mature, Q-Chem-backed method stack spanning semi-empirical through post-HF and composite thermochemistry, (b) strong 3D visualization of orbitals/density/property maps, (c) enormous precomputed reference databases, and (d) newer ML-accelerated workflows for conformer energetics and NMR prediction, makes it well suited to organic/medicinal chemistry structure elucidation, reaction/selectivity rationalization, and chemistry education — while being less suited to large-scale periodic, materials, or HPC-cluster-scale electronic structure work, where open, more extensible, cluster-native codes are typically preferred.

---

## 11. Key Publications Underpinning Spartan's Theory and Methods

### Software / engine itself
- Hehre, W. J.; Ohlinger, W. S. *Spartan'14 Tutorial and User's Guide*. Wavefunction, Inc.: Irvine, CA, 2013.
- Hehre, W. J. *A Guide to Molecular Mechanics and Quantum Chemical Calculations*. Wavefunction, Inc.: Irvine, CA, 2003. ISBN 1-890661-06-6.
- Krylov, A. I.; Gill, P. M. W. "Q-Chem: an engine for innovation." *WIREs Computational Molecular Science* **2013**, *3*(3), 317–326.

### Density functional theory — foundational
- Hohenberg, P.; Kohn, W. "Inhomogeneous Electron Gas." *Physical Review* **1964**, *136*(3B), B864–B871.

### Exchange functionals
- Dirac, P. A. M. "Note on Exchange Phenomena in the Thomas Atom." *Math. Proc. Cambridge Phil. Soc.* **1930**, *26*(3), 376–385.
- Becke, A. D. "Density-functional exchange-energy approximation with correct asymptotic behavior." *Physical Review A* **1988**, *38*(6), 3098–3100. (B88)
- Gill, P. M. W. "A new gradient-corrected exchange functional." *Molecular Physics* **1996**, *89*(2), 433–445.
- Gilbert, A. T. B.; Gill, P. M. W. "Decomposition of exchange-correlation energies." *Chemical Physics Letters* **1999**, *312*(5–6), 511–521. (GG99/Gill96 context)
- Perdew, J. P.; Wang, Y. "Accurate and simple analytic representation of the electron-gas correlation energy." *Physical Review B* **1992**, *45*(23), 13244–13249. (PW91)

### Correlation functionals and hybrids
- Perdew, J. P. "Density-functional approximation for the correlation energy of the inhomogeneous electron gas." *Physical Review B* **1986**, *33*(12), 8822–8824.
- Lee, C.; Yang, W.; Parr, R. G. "Development of the Colle-Salvetti correlation-energy formula into a functional of the electron density." *Physical Review B* **1988**, *37*(2), 785–789. (LYP)
- Vosko, S. H.; Wilk, L.; Nusair, M. "Accurate spin-dependent electron liquid correlation energies for local spin density calculations." *Canadian Journal of Physics* **1980**, *58*(8), 1200–1211. (VWN)
- Perdew, J. P.; Zunger, A. "Self-interaction correction to density-functional approximations for many-electron systems." *Physical Review B* **1986**, *33*(12), 8822–8824.
- Perdew, J. P.; Burke, K.; Ernzerhof, M. "Generalized Gradient Approximation Made Simple." *Physical Review Letters* **1996**, *77*(18), 3865–3868. (PBE)
- Stephens, P. J.; Devlin, F. J.; Chabalowski, C. F.; Frisch, M. J. "Ab Initio Calculation of Vibrational Absorption and Circular Dichroism Spectra Using Density Functional Force Fields." *Journal of Physical Chemistry* **1994**, *98*(45), 11623–11627. (B3LYP)
- Adamson, R. D.; Gill, P. M. W.; Pople, J. A. "Empirical density functionals." *Chemical Physics Letters* **1998**, *284*(5–6), 6–11. (EDF1)
- Gill, P. M. W.; Ching, Y. L.; George, M. W. "EDF2: A density functional for predicting molecular vibrational frequencies." *Australian Journal of Chemistry* **2004**, *57*(4), 365–370.

### Meta-GGA / Minnesota functionals
- Zhao, Y.; Schultz, N. E.; Truhlar, D. G. "Design of Density Functionals by Combining the Method of Constraint Satisfaction with Parameterization for Thermochemistry, Thermochemical Kinetics, and Noncovalent Interactions." *Journal of Chemical Theory and Computation* **2006**, *2*(2), 364–382. (M05)
- Zhao, Y.; Truhlar, D. G. "A new local density functional for main-group thermochemistry, transition metal bonding, thermochemical kinetics, and noncovalent interactions." *Journal of Chemical Physics* **2006**, *125*(19), 194101.
- Zhao, Y.; Truhlar, D. G. "Density Functional for Spectroscopy: No Long-Range Self-Interaction Error, Good Performance for Rydberg and Charge-Transfer States, and Better Performance on Average than B3LYP for Ground States." *Journal of Physical Chemistry A* **2006**, *110*(49), 13126–13130.
- Zhao, Y.; Truhlar, D. G. "The M06 suite of density functionals for main group thermochemistry, thermochemical kinetics, noncovalent interactions, excited states, and transition elements." *Theoretical Chemistry Accounts* **2008**, *120*(1–3), 215–241.

### Long-range-corrected / dispersion-corrected functionals
- Chai, J.-D.; Head-Gordon, M. "Systematic optimization of long-range corrected hybrid density functionals." *Journal of Chemical Physics* **2008**, *128*(8), 084106. (ωB97)
- Chai, J.-D.; Head-Gordon, M. "Long-range corrected hybrid density functionals with damped atom-atom dispersion corrections." *Physical Chemistry Chemical Physics* **2008**, *10*(44), 6615–6620. (ωB97X-D)

### Force fields
- Halgren, T. A. "Merck molecular force field. I. Basis, form, scope, parameterization, and performance of MMFF94." *Journal of Computational Chemistry* **1996**, *17*(5–6), 490–519.
- Clark, M.; Cramer, R. D., III; Van Opdenbosch, N. "Validation of the general purpose tripos 5.2 force field." *Journal of Computational Chemistry* **1989**, *10*(8), 982–1012. (SYBYL)

### Semi-empirical methods
- Dewar, M. J. S.; Thiel, W. "Ground states of molecules. 38. The MNDO method." *Journal of the American Chemical Society* **1977**, *99*(15), 4899–4907.
- Dewar, M. J. S.; Zoebisch, E. G.; Healy, E. F.; Stewart, J. J. P. "Development and use of quantum molecular models. 75." *Journal of the American Chemical Society* **1985**, *107*(13), 3902–3909. (AM1)
- Stewart, J. J. P. "Optimization of parameters for semiempirical methods I. Method." *Journal of Computational Chemistry* **1989**, *10*(2), 209–220. (PM3)
- Stewart, J. J. P. "Optimization of parameters for semiempirical methods II. Applications." *Journal of Computational Chemistry* **1989**, *10*(2), 221–264.
- Stewart, J. J. P. "Optimization of parameters for semiempirical methods. III." *Journal of Computational Chemistry* **1991**, *12*(3), 320–341.
- Stewart, J. J. P. "Optimization of parameters for semiempirical methods IV." *Journal of Molecular Modeling* **2004**, *10*(2), 155–164.
- Rocha, G. B.; Freire, R. O.; Simas, A. M.; Stewart, J. J. P. "RM1: A reparameterization of AM1 for H, C, N, O, P, S, F, Cl, Br, and I." *Journal of Computational Chemistry* **2006**, *27*(10), 1101–1111.
- Stewart, J. J. P. "Optimization of Parameters for Semiempirical Methods V." *Journal of Molecular Modeling* **2007**, *13*(12), 1173–1213. (PM6)

### Hartree–Fock, basis sets, and implicit solvation
- Ditchfield, R.; Hehre, W. J.; Pople, J. A. "Self-Consistent Molecular-Orbital Methods. IX. An Extended Gaussian-Type Basis for Molecular-Orbital Studies of Organic Molecules." *Journal of Chemical Physics* **1971**, *54*(2), 724–728.
- Marenich, A. V.; Olson, R. M.; Kelly, C. P.; Cramer, C. J.; Truhlar, D. G. "Self-Consistent Reaction Field Model for Aqueous and Nonaqueous Solutions Based on Accurate Polarized Partial Charges." *Journal of Chemical Theory and Computation* **2007**, *3*(6), 2011–2033. (SM8)

### Post-Hartree–Fock / correlated methods
- Møller, C.; Plesset, M. S. "Note on an Approximation Treatment for Many-Electron Systems." *Physical Review* **1934**, *46*(7), 618–622.
- Head-Gordon, M.; Pople, J. A.; Frisch, M. J. "MP2 energy evaluation by direct methods." *Chemical Physics Letters* **1988**, *153*(6), 503–506.
- Purvis, G. D.; Bartlett, R. J. "A full coupled-cluster singles and doubles model." *Journal of Chemical Physics* **1982**, *76*(4), 1910–1919.
- Raghavachari, K.; Trucks, G. W.; Pople, J. A.; Head-Gordon, M. "A fifth-order perturbation comparison of electron correlation theories." *Chemical Physics Letters* **1989**, *157*(6), 479–483.
- Pople, J. A.; Head-Gordon, M.; Raghavachari, K. "Quadratic configuration interaction." *Journal of Chemical Physics* **1987**, *87*(10), 5968–5975.
- Weigend, F.; Häser, M. "RI-MP2: first derivatives and global consistency." *Theoretical Chemistry Accounts* **1997**, *97*(1–4), 331–340.
- Distasio, R. A., Jr.; Steele, R. P.; Rhee, Y. M.; Shao, Y.; Head-Gordon, M. "An improved algorithm for analytical gradient evaluation in resolution-of-the-identity second-order Møller-Plesset perturbation theory." *Journal of Computational Chemistry* **2007**, *28*(5), 839–856.
- Krylov, A. I.; Sherrill, C. D.; Byrd, E. F. C.; Head-Gordon, M. "Size-consistent wave functions for nondynamical correlation energy." *Journal of Chemical Physics* **1998**, *109*(24), 10669–10678.

### Excited states / TDDFT / CIS
- Runge, E.; Gross, E. K. U. "Density-Functional Theory for Time-Dependent Systems." *Physical Review Letters* **1984**, *52*(12), 997–1000.
- Hirata, S.; Head-Gordon, M. "Time-dependent density functional theory for radicals." *Chemical Physics Letters* **1999**, *302*(5–6), 375–382.
- Maurice, D.; Head-Gordon, M. "Analytical second derivatives for excited electronic states using the single excitation configuration interaction method." *Molecular Physics* **1999**, *96*(10), 1533–1541.
- Head-Gordon, M.; Rico, R. J.; Oumi, M.; Lee, T. J. "A doubles correction to electronic excited states from configuration interaction in the space of single substitutions." *Chemical Physics Letters* **1994**, *219*(1–2), 21–29. (CIS(D))
- Rhee, Y. M.; Head-Gordon, M. "Scaled Second-Order Perturbation Corrections to Configuration Interaction Singles." *Journal of Physical Chemistry A* **2007**, *111*(24), 5314–5326. (RI-CIS(D))

### Composite/thermochemical methods
- Curtiss, L. A.; Raghavachari, K.; Trucks, G. W.; Pople, J. A. "Gaussian-2 theory for molecular energies of first- and second-row compounds." *Journal of Chemical Physics* **1991**, *94*(11), 7221–7231.
- Curtiss, L. A.; Raghavachari, K.; Redfern, P. C.; Rassolov, V.; Pople, J. A. "Gaussian-3 (G3) theory for molecules containing first and second-row atoms." *Journal of Chemical Physics* **1998**, *109*(18), 7764–7776.
- Curtiss, L. A.; Redfern, P. C.; Raghavachari, K.; Rassolov, V.; Pople, J. A. "Gaussian-3 theory using reduced Møller-Plesset order." *Journal of Chemical Physics* **1998**, *110*(10), 4703–4710. (G3(MP2))
- Ohlinger, W. S.; Klunzinger, P. E.; Deppmeier, B. J.; Hehre, W. J. "Efficient Calculation of Heats of Formation." *Journal of Physical Chemistry A* **2009**, *113*(10), 2165–2175. (T1)

### NMR / spectral prediction
- Wolinski, K.; Hinton, J. F.; Pulay, P. "Efficient implementation of the gauge-independent atomic orbital method for NMR chemical shift calculations." *Journal of the American Chemical Society* **1990**, *112*(23), 8251–8260. (GIAO)
- Kussmann, J.; Ochsenfeld, C. "Linear-scaling method for calculating nuclear magnetic resonance chemical shifts using gauge-including atomic orbitals within Hartree-Fock and density-functional theory." *Journal of Chemical Physics* **2007**, *127*(5), 054103.
- Hehre, W.; Klunzinger, P.; Deppmeier, B.; Driessen, A.; Uchida, N.; Hashimoto, M.; Fukushi, E.; Takata, Y. "Efficient Protocol for Accurately Calculating ¹³C Chemical Shifts of Conformationally Flexible Natural Products: Scope, Assessment, and Limitations." *Journal of Natural Products* **2019**, *82*(8), 2299–2306.
- Scott, A. P.; Radom, L. "Harmonic Vibrational Frequencies: An Evaluation of Hartree−Fock, Møller−Plesset, Quadratic Configuration Interaction, Density Functional Theory, and Semiempirical Scale Factors." *Journal of Physical Chemistry* **1996**, *100*(41), 16502–16513.
- Johnson, B. G.; Florián, J. "The prediction of Raman spectra by density functional theory. Preliminary findings." *Chemical Physics Letters* **1995**, *247*(1–2), 120–125.

### Machine learning models (recent Spartan-specific work)
- Hehre, T.; Klunzinger, P. E.; Deppmeier, B.; Ohlinger, W.; Hehre, W. J. "Practical Machine Learning Strategies. I. Correcting the MMFF Molecular Mechanics Model to More Accurately Provide Conformational Energy Differences in Flexible Organic Molecules." *Journal of Computational Chemistry* **2025**, *46*(1), e70016.
- Hehre, T.; Klunzinger, P. E.; Deppmeier, B.; Ohlinger, W.; Hehre, W. J. "Accurate Prediction of ωB97X-D/6-31G* Equilibrium Geometries from a Neural Net Starting from Merck Molecular Force Field (MMFF) Molecular Mechanics Geometries." *Journal of Chemical Information and Modeling* **2025**, *65*(5), 2314–2321.

### Education-focused literature using Spartan
- Shusterman, A. J.; Shusterman, G. P. "Teaching Chemistry with Electron Density Models." *Journal of Chemical Education* **1997**, *74*(7), 771–775.
- Hehre, W. J.; Shusterman, A.; Nelson, J. *Molecular Modeling Workbook for Organic Chemistry.* Wavefunction, Inc., 1998.
- Linenberger, K. J.; Cole, R. S.; Sarkar, S. "Looking Beyond Lewis Structures: A General Chemistry Modeling Experiment Focusing on Physical Properties and Geometry." *Journal of Chemical Education* **2011**, *88*(7), 962–965.
- Kim, H.; Sulaimon, S.; Menezes, S.; Son, A.; Menezes, W. J. C. "A Comparative Study of Successful Central Nervous System Drugs Using Molecular Modeling." *Journal of Chemical Education* **2011**, *88*(10), 1389–1393.

---

*Compiled from Wavefunction, Inc. product documentation (wavefun.com), the Spartan (chemistry software) Wikipedia entry and its cited primary literature, and institutional software listings (UCSD, University of Michigan).*


---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the Spartan 	Commercial molecular modeling package offering DFT alongside HF/semi-empirical methods, aimed at broad usability.. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
