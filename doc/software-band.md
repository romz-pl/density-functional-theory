# BAND: A Slater-Type-Orbital Periodic DFT Code for Solids, Surfaces, and Semiconductors

*A review of the BAND engine within the Amsterdam Modeling Suite (AMS), Software for Chemistry & Materials (SCM), Vrije Universiteit Amsterdam.*

---

## 1. Overview

BAND is the periodic (crystalline/slab/wire) density-functional theory engine of the Amsterdam Modeling Suite. It is the "sister code" of the molecular engine ADF, sharing the same numerical philosophy: an all-electron, atom-centered Slater-type-orbital (STO) basis, three-dimensional numerical integration of the Kohn–Sham matrix elements, and auxiliary-basis density fitting for the Coulomb potential. Where ADF treats finite (0-D) molecular systems, BAND extends this machinery to systems periodic in one, two, or three dimensions — polymers/nanotubes (1-D), slabs/surfaces/2-D materials (2-D), and bulk crystals (3-D) — using Bloch sums of local atomic orbitals rather than plane waves.

Because BAND is a Linear-Combination-of-Atomic-Orbitals (LCAO) code rather than a plane-wave/pseudopotential code, it does not need artificial vacuum-supercell periodicity to describe surfaces or low-dimensional systems: a 2-D slab is treated as *genuinely* two-dimensional, and a 1-D chain or nanotube as genuinely one-dimensional, without spurious images along the non-periodic direction(s). This is presented by SCM as one of BAND's principal differentiators relative to plane-wave codes such as VASP or Quantum ESPRESSO (both of which are also distributed with AMS as complementary plane-wave engines for dense 3-D bulk systems).

BAND is developed and maintained by SCM (Theoretical Chemistry, Vrije Universiteit Amsterdam) and is licensed commercially as part of the Amsterdam Modeling Suite, alongside ADF (molecular DFT), DFTB, ReaxFF, force fields, machine-learning potentials, and COSMO-RS.

---

## 2. Theoretical & Numerical Foundations

### 2.1 Basis set: numerical and Slater-type atomic orbitals
- BAND (like ADF) is an **all-electron** code: there are no pseudopotentials, norm-conserving or otherwise, and no frozen-core approximation is *required* (though frozen cores are available as an efficiency option).
- The basis functions are **Slater-type orbitals (STOs)**, which — unlike Gaussian-type orbitals — reproduce the correct cusp behavior of the wavefunction at the nucleus and the correct exponential decay far from it, giving compact, physically motivated basis sets.
- BAND additionally supports **numerical atomic orbitals (NAOs)** for the core region in some contexts, allowing a very accurate description of tightly bound core states without inflating the STO basis, which historically permitted accurate *total energies* (not just relative/formation energies) even for heavy elements.
- Standard STO basis-set families shared with ADF span the whole periodic table (Z = 1–118): SZ, DZ, DZP, TZP, TZ2P, QZ4P, plus even-tempered and diffuse-augmented sets, with matching frozen-core definitions when core freezing is used.
- Because the same basis-set and relativistic-treatment machinery is shared between ADF (molecular) and BAND (periodic), results are directly comparable across the cluster ↔ periodic-slab ↔ bulk-crystal continuum — useful for surface-science workflows that compare adsorbate/cluster and extended-slab pictures on equal footing.

### 2.2 Bloch functions, k-space, and periodicity of arbitrary dimension
- Atomic STOs are combined into symmetry-adapted **Bloch sums**; the Kohn–Sham problem is solved self-consistently at a discrete, user-controllable (or automatic) set of **k-points** sampling the Brillouin zone.
- BAND natively supports **0-D, 1-D, 2-D, and 3-D periodicity** in a single unified formalism — polymers and nanotubes (1-D), free-standing slabs and 2-D materials (2-D), and bulk crystals (3-D) — treating the non-periodic directions as truly open/vacuum rather than via artificial supercell replication.
- **k-space integration** uses analytic/semi-analytic quadratic-integration schemes over the 2-D and 3-D Brillouin zone developed specifically for this code (Wiesenekker/te Velde/Baerends), rather than simple linear tetrahedron or Monkhorst–Pack summation alone, improving convergence with respect to k-mesh density.

### 2.3 Density fitting and numerical integration
- The Coulomb (Hartree) potential and, where relevant, exchange integrals are evaluated via **auxiliary-basis density fitting**, avoiding explicit 4-index two-electron integral evaluation; the fit is handled with schemes tailored for periodic/molecular systems (Franchini/Philipsen/van Lenthe/Visscher density-fitting scheme).
- Exchange–correlation contributions are evaluated by **3-D numerical integration** on a fuzzy-cell (Becke-type) atomic grid, adapted to periodic systems.
- **Multi-resolution** integration/fitting strategies are available to balance cost against accuracy for large or heterogeneous systems.

### 2.4 Self-consistent field (SCF) and convergence
- Modern SCF convergence accelerators are provided (LISTi, EDIIS, ARH-type schemes), aimed at difficult cases such as small-gap or metallic/near-degenerate systems.
- Spin-unrestricted and non-collinear/spin-orbit self-consistent treatments are supported (see §2.5), including antiferromagnetic ordering.

### 2.5 Relativistic treatment
- BAND shares ADF's relativistic Hamiltonians, principally the **ZORA** (Zeroth-Order Regular Approximation) scheme, applicable scalar-relativistically or with **self-consistent spin–orbit coupling**, enabling efficient and accurate treatment of heavy elements throughout the periodic table without pseudopotentials.
- A **finite-nucleus model** is available in place of the point-charge nucleus for high-accuracy core-property work (hyperfine, NMR, EFG).

### 2.6 Exchange–correlation functionals
- BAND supports the standard hierarchy of density functionals: LDA, GGAs (e.g., PBE), meta-GGAs, and modern functionals such as **r²SCAN-D4**; **hybrid functionals** with a fraction of exact (Hartree–Fock-type) exchange are available via a periodic **Hartree–Fock RI (resolution-of-identity)** implementation.
- Full integration with the **LibXC** functional library extends the range of available XC functionals well beyond SCM's native set.
- **Dispersion corrections**: Grimme's **D3** (with Becke–Johnson damping) and **D4** schemes are available for van der Waals interactions in sparse/layered materials.
- Band-gap-correcting approaches beyond standard XC are available for semiconductors/insulators: **DFT+U** (Hubbard correction, e.g., for correlated oxides such as NiO), **model potentials** such as **TB-mBJ** (Tran–Blaha modified Becke–Johnson), and the **DFT-1/2** method — each aimed at correcting the well-known DFT band-gap underestimation without the cost of hybrids or many-body perturbation theory.
- A many-body perturbation theory (**MBPT**) route, including a **GW** approximation (eigenvalue-only self-consistent GW on top of, e.g., PBE), is available for quasiparticle band structures with STO or Gaussian-type basis input.

---

## 3. Ground-State and Structural Capabilities

- **Geometry optimization** and **lattice optimization** (relaxing both atomic positions and unit-cell parameters) with analytic energy gradients, formulated specifically for the periodic STO/NAO Bloch basis (Kadantsev, Klooster, de Boeij, Ziegler).
- **Analytic atomic forces** and **stress/elastic tensor** calculation for mechanical-property prediction (elastic constants, bulk/shear moduli).
- **Phonon calculations** (vibrational frequencies of periodic systems) via the AMS driver's numerical-differentiation "Properties" framework, usable for thermodynamic and lattice-dynamical analysis (e.g., graphene phonon dispersions).
- **Transition-state search** and **nudged elastic band (NEB)** workflows for periodic reaction pathways (e.g., surface reactions, diffusion barriers), inherited from the shared AMS driver.
- **Homogeneous external electric and magnetic fields**, useful for studying field-induced band-gap closing in 2-D semiconductors, Stark effects, and magnetic response.
- **Continuum solvation** at surfaces/interfaces via **COSMO** (Conductor-like Screening Model) and **SM12**, relevant to electrochemistry and catalysis in solution.

---

## 4. Electronic-Structure Analysis Tools

BAND provides an extensive suite of post-SCF analysis tools that are a particular strength of the atomic-orbital approach, since orbitals and density partitions map naturally onto chemically meaningful quantities:

- **Density of states (DOS)** and **projected/partial DOS (PDOS)**, decomposable by atom, orbital character, or user-defined fragment.
- **Band structures**, including **"fat bands"** (band structures weighted/colored by orbital or fragment character) along user-specified or automatically generated high-symmetry k-paths, with restart capability so DOS/band-structure post-processing can be redone without repeating the SCF.
- **COOP** (Crystal Orbital Overlap Population) and related bonding-indicator analyses for identifying bonding/antibonding/non-bonding character between atoms or fragments in extended solids — a signature BAND/ADF-family capability tracing to the atomic-orbital basis.
- **QTAIM** (Quantum Theory of Atoms in Molecules / Bader) and **ELF** (Electron Localization Function) analysis for real-space bonding and topological charge-density analysis.
- **Periodic Energy Decomposition Analysis (PEDA)**, including **PEDA-NOCV** (Natural Orbitals for Chemical Valence) — decomposing bonding/interaction energies (e.g., adsorbate–surface interactions) into physically interpretable orbital-interaction, electrostatic, Pauli-repulsion, and dispersion terms, extending ADF's well-known molecular EDA to periodic and semi-periodic (fragment-in-slab) systems.
- **Fragment analysis**: systems can be decomposed into interacting fragments (e.g., adsorbate + slab) with shared orbital/density integration, enabling chemisorption/physisorption analysis directly analogous to cluster calculations (e.g., CO on a Cu or MgO surface).
- **Local density of states (LDOS)** for simulated **STM** images.
- **Charge analysis**: Bader/QTAIM charges, Mulliken-type and other population schemes.
- **Effective mass** calculation at arbitrary k-points, directly relevant to charge-carrier mobility prediction in semiconductors.
- **3-D field visualization**: real-space plotting of densities, potentials, and orbitals on a user-specified grid.
- **Electronic transport (NEGF)**: Non-Equilibrium Green's Function methodology for two-probe/molecular-junction conductance calculations, including with applied bias, extending BAND (and ADF) to quantum-transport problems.

---

## 5. Spectroscopic Properties

- **NMR chemical shifts** in periodic systems, using a **Gauge-Including Atomic Orbitals (GIAO)** formulation adapted to Bloch-periodic boundary conditions.
- **EPR/ESR** parameters: **g-tensors** (Zeeman interaction) and **hyperfine A-tensors**, formulated for periodic systems using numerical and Slater-type orbitals — of particular use for paramagnetic defects in solids.
- **Electric field gradients (EFG)** and nuclear quadrupole interactions, plus **Q-tensor** properties.
- **Properties at nuclei** and **X-ray form factors**.
- **Optical/dielectric response** via **Time-Dependent Current-DFT (TD-CDFT)**, including treatment of both non-metallic (dielectric-function) and metallic response, with the **Vignale–Kohn** current-functional extension for improved optical spectra of metals and semiconductors; **EELS** (electron energy loss spectroscopy) is available from the same response-property machinery.
- **Relativistic (spin-orbit) TD-CDFT** for linear optical response of solids containing heavy elements.
- **Dipole moment and Berry-phase** polarization calculations for periodic systems.

---

## 6. Efficiency, Scaling, and Practical Features

- **Parallelization** across many CPU cores for both SCF and property calculations.
- **Linear-scaling techniques and distance cut-offs** for the atomic-orbital/density-fitting machinery, improving tractability for larger unit cells and low-symmetry, sparse structures.
- **Space-group and point-group symmetry** exploitation to reduce computational cost.
- **Restart infrastructure**: SCF, DOS, and band-structure post-processing can each be restarted independently from checkpoint (KF) files, useful for large screening or iterative-analysis workflows.
- Tight integration with the **AMS driver** and graphical user interface (GUI) for job setup, band-structure/Brillouin-zone path selection, and result visualization, plus scripting via **PLAMS** (Python).

---

## 7. Application Domains: Solids, Surfaces, and Semiconductors

Consistent with its atomic-orbital, all-electron design, BAND is marketed and used chiefly for:

- **Surfaces and interfaces**: true 2-D slab periodicity (no vacuum-image artifacts), adsorption energetics, work functions, surface reconstruction, heterogeneous catalysis mechanisms, and interface band alignment.
- **Low-dimensional and sparse materials**: nanotubes, 2-D materials (e.g., transition-metal dichalcogenide monolayers such as MoS₂), polymers, layered/van-der-Waals solids (via D3/D4 dispersion).
- **Semiconductors**: band-gap and band-structure prediction (with the DFT-gap-underestimation problem addressed via DFT+U, model potentials such as TB-mBJ, DFT-1/2, or GW), effective-mass/mobility prediction, optical/dielectric spectra via TD-CDFT, defect and paramagnetic-center EPR/hyperfine analysis, and field-effect studies (e.g., band-gap closing of 2-D semiconductors under an applied electric field).
- **Bulk crystalline solids**: total-energy, elastic, phonon, magnetic (including antiferromagnetic spin ordering), and core-level spectroscopic (NMR/EPR/EFG) properties, all-electron and relativistic where needed (heavy-element compounds, actinide/lanthanide solids).
- **Batteries, catalysis, and materials-science screening** more broadly, as part of AMS's integrated multiscale offering (kinetic Monte Carlo, microkinetics, machine-learning potentials for larger-scale/longer-time extensions beyond DFT).

---

## 8. Positioning Relative to Plane-Wave Codes

SCM explicitly markets BAND as the atomic-orbital complement to the plane-wave engines bundled with AMS (Quantum ESPRESSO, and via interface, VASP):

| | BAND (LCAO/STO) | Plane-wave codes (QE, VASP) |
|---|---|---|
| Basis | All-electron Slater-type (+ numerical) atomic orbitals | Plane waves + pseudopotentials |
| Core treatment | All electrons, no pseudopotentials | Pseudopotential/PAW required |
| Best suited for | Surfaces, 1-D/2-D systems, sparse matter, relativistic/heavy elements, detailed chemical/orbital analysis | Dense 3-D bulk periodic systems, large supercells |
| Chemical-bonding analysis | COOP, PEDA(-NOCV), fragment analysis, QTAIM, ELF | Generally more limited without post-processing |
| Vacuum handling for surfaces | True 2-D periodicity, no vacuum slab needed | Requires vacuum-gap supercell |

---

## 9. Publications Related to BAND's Theory

### 9.1 Core method and required general citations
1. G. te Velde and E. J. Baerends, *Precise density-functional method for periodic structures*, **Physical Review B 44, 7888 (1991)**. https://doi.org/10.1103/PhysRevB.44.7888
2. BAND 2026.1, SCM, Theoretical Chemistry, Vrije Universiteit, Amsterdam, The Netherlands, https://www.scm.com (program citation; contributor list per SCM's Required Citations page includes P.H.T. Philipsen, G. te Velde, E.J. Baerends, J.A. Berger, P.L. de Boeij, M. Franchini, J.A. Groeneveld, E.S. Kadantsev, R. Klooster, F. Kootstra, M.C.W.M. Pols, P. Romaniello, M. Raupach, D.G. Skachkov, J.G. Snijders, C.J.O. Verzijl, J.A. Celis Gil, J.M. Thijssen, G. Wiesenekker, C.A. Peeples, G. Schreckenbach, T. Ziegler).
3. E. J. Baerends, T. Ziegler, et al. (ADF program suite, foundational molecular-DFT methodology shared with BAND) — see the general Amsterdam Modeling Suite overview paper below.
4. Amsterdam Modeling Suite overview: *The Amsterdam Modeling Suite*, **The Journal of Chemical Physics 162, 162501 (2025)** (AIP Publishing) — describes the shared design principles (density fitting, numerical integration, Slater-type orbitals) underlying ADF/BAND.

### 9.2 Brillouin-zone integration
5. G. Wiesenekker, G. te Velde and E. J. Baerends, *Analytic quadratic integration over the two-dimensional Brillouin zone*, **Journal of Physics C: Solid State Physics 21, 4263 (1988)**. https://doi.org/10.1088/0022-3719/21/23/012
6. G. Wiesenekker and E. J. Baerends, *Quadratic integration over the three-dimensional Brillouin zone*, **Journal of Physics: Condensed Matter 3, 6721 (1991)**. https://doi.org/10.1088/0953-8984/3/35/005

### 9.3 Numerical integration and density fitting
7. M. Franchini, P. H. T. Philipsen, L. Visscher, *The Becke Fuzzy Cells Integration Scheme in the Amsterdam Density Functional Program Suite*, **Journal of Computational Chemistry 34, 1818 (2013)**. https://doi.org/10.1002/jcc.23323
8. M. Franchini, P. H. T. Philipsen, E. van Lenthe, L. Visscher, *Accurate Coulomb Potentials for Periodic and Molecular Systems through Density Fitting*, **Journal of Chemical Theory and Computation 10, 1994 (2014)**. https://doi.org/10.1021/ct500172n

### 9.4 Geometry optimization / analytic gradients
9. E. S. Kadantsev, R. Klooster, P. L. de Boeij and T. Ziegler, *The Formulation and Implementation of Analytic Energy Gradients for Periodic Density Functional Calculations with STO/NAO Bloch Basis Set*, **Molecular Physics 105, 2583 (2007)**. https://doi.org/10.1080/00268970701598063

### 9.5 Time-dependent (current-)DFT and optical response
10. F. Kootstra, P. L. de Boeij and J. G. Snijders, *Efficient real-space approach to time-dependent density functional theory for the dielectric response of nonmetallic crystals*, **Journal of Chemical Physics 112, 6517 (2000)**. https://doi.org/10.1063/1.481315
11. P. Romaniello and P. L. de Boeij, *Time-dependent current-density-functional theory for the metallic response of solids*, **Physical Review B 71, 155108 (2005)**. https://doi.org/10.1103/PhysRevB.71.155108
12. F. Kootstra, P. L. de Boeij, and J. G. Snijders, *Application of time-dependent density-functional theory to the dielectric function of various nonmetallic crystals*, **Physical Review B 62, 7071 (2000)**. https://doi.org/10.1103/PhysRevB.62.7071
13. P. Romaniello, P. L. de Boeij, F. Carbone, and D. van der Marel, *Optical properties of bcc transition metals in the range 0–40 eV*, **Physical Review B 73, 075115 (2006)**. https://doi.org/10.1103/PhysRevB.73.075115
14. P. Romaniello and P. L. de Boeij, *Relativistic two-component formulation of time-dependent current-density functional theory: Application to the linear response of solids*, **Journal of Chemical Physics 127, 174111 (2007)**. https://doi.org/10.1063/1.2780146

### 9.6 Vignale–Kohn current functional
15. J. A. Berger, P. L. de Boeij and R. van Leeuwen, *Analysis of the viscoelastic coefficients in the Vignale-Kohn functional: The cases of one- and three-dimensional polyacetylene*, **Physical Review B 71, 155104 (2005)**. https://doi.org/10.1103/PhysRevB.71.155104
16. J. A. Berger, P. Romaniello, R. van Leeuwen and P. L. de Boeij, *Performance of the Vignale-Kohn functional in the linear response of metals*, **Physical Review B 74, 245117 (2006)**. https://doi.org/10.1103/PhysRevB.74.245117
17. J. A. Berger, P. L. de Boeij, and R. van Leeuwen, *Analysis of the Vignale-Kohn current functional in the calculation of the optical spectra of semiconductors*, **Physical Review B 75, 035116 (2007)**. https://doi.org/10.1103/PhysRevB.75.035116

### 9.7 NMR
18. D. Skachkov, M. Krykunov, E. Kadantsev, and T. Ziegler, *The Calculation of NMR Chemical Shifts in Periodic Systems Based on Gauge Including Atomic Orbitals and Density Functional Theory*, **Journal of Chemical Theory and Computation 6, 1650 (2010)**. https://doi.org/10.1021/ct100046a
19. D. Skachkov, M. Krykunov, and T. Ziegler, *An improved scheme for the calculation of NMR chemical shifts in periodic systems based on gauge including atomic orbitals and density functional theory*, **Canadian Journal of Chemistry 89, 1150 (2011)**. https://doi.org/10.1139/v11-050

### 9.8 EPR/ESR (hyperfine and g-tensor)
20. E. S. Kadantsev and T. Ziegler, *Implementation of a Density Functional Theory-Based Method for the Calculation of the Hyperfine A-tensor in Periodic Systems with the Use of Numerical and Slater Type Atomic Orbitals: Application to Paramagnetic Defects*, **Journal of Physical Chemistry A 112, 4521 (2008)**. https://doi.org/10.1021/jp800494m
21. E. S. Kadantsev and T. Ziegler, *Implementation of a DFT Based Method for the Calculation of Zeeman g-tensor in Periodic Systems with the use of Numerical and Slater Type Atomic Orbitals*, **Journal of Physical Chemistry A 113, 1327 (2009)**. https://doi.org/10.1021/jp805466c

### 9.9 Electronic transport (NEGF)
22. C. J. O. Verzijl and J. M. Thijssen, *DFT-Based Molecular Transport Implementation in ADF/BAND*, **Journal of Physical Chemistry C 116, 24393 (2012)**. https://doi.org/10.1021/jp3044225

### 9.10 Electron energy density / in-situ electronegativity
23. S. Racioppi and M. Rahm, *In-Situ Electronegativity and the Bridging of Chemical Bonding Concepts*, **Chemistry – A European Journal 27, 18156 (2021)**. https://doi.org/10.1002/chem.202103477
24. S. Racioppi, P. Hyldgaard and M. Rahm, *Quantifying Atomic Volume, Partial Charge, and Electronegativity in Condensed Phases*, **Journal of Physical Chemistry C 128, 4009 (2024)**. https://doi.org/10.1021/acs.jpcc.3c07677

### 9.11 Doctoral theses documenting BAND's TD-(C)DFT methodology
25. F. Kootstra, Ph.D. thesis, Rijksuniversiteit Groningen (2001).
26. P. Romaniello, Ph.D. thesis, Rijksuniversiteit Groningen (2006).
27. J. A. Berger, Ph.D. thesis, Rijksuniversiteit Groningen (2006).

---

## 10. Summary Assessment

BAND occupies a distinctive niche among periodic DFT codes: rather than the plane-wave/pseudopotential paradigm dominant in solid-state physics software (VASP, Quantum ESPRESSO, CASTEP, CP2K's GPW mode, etc.), it applies the all-electron, Slater-type-orbital, numerical-integration methodology of the ADF molecular-DFT program to systems with 1-, 2-, and 3-D periodicity. This yields particular strengths in:

- **True low-dimensional periodicity** (surfaces, 2-D materials, nanotubes/wires) without vacuum-supercell artifacts;
- **All-electron relativistic treatment** (ZORA, spin-orbit) of heavy elements without pseudopotential approximations;
- **Chemically interpretable bonding analysis** (COOP, PEDA/PEDA-NOCV, fragment decomposition, QTAIM/ELF) that leverages the atomic-orbital basis directly;
- **A broad spectroscopic toolkit** (NMR, EPR/ESR, EFG, TD-CDFT optical/EELS response) implemented natively for periodic systems, which is comparatively rare in plane-wave codes;
- **Semiconductor-focused features** — DFT+U, TB-mBJ, DFT-1/2, and GW routes to improved band gaps, plus effective-mass and mobility analysis — addressing the well-known DFT band-gap problem from several angles within one code.

Its principal trade-off relative to plane-wave codes is efficiency on very large, dense, high-symmetry 3-D bulk supercells, which is why SCM bundles BAND alongside the plane-wave engine Quantum ESPRESSO (and a VASP interface) within AMS, letting users choose the LCAO or plane-wave approach according to system dimensionality, element composition, and the type of analysis required.


---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the BAND 	Periodic-system DFT code (Slater-type orbitals) within the Amsterdam Modeling Suite, for solids, surfaces, and semiconductors. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
