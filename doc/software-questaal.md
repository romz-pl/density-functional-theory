# Questaal (LMTO/ASA): A Review of the Electronic Structure Suite

## 1. Overview

Questaal is an open-source suite of first-principles electronic structure codes built on the linear muffin-tin orbital (LMTO) method. It descends from the LMTO formalism developed by O. K. Andersen's group in Stuttgart in the 1980s, and has been developed mainly by Mark van Schilfgaarde, with substantial contributions from Dimitar Pashov, Swagata Acharya, Walter R. L. Lambrecht, Jerome Jackson, Kirill D. Belashchenko, Athanasios Chantis, Francois Jamet, Takao Kotani, and others. The suite was designated a CCP9 flagship code in 2015, with EPSRC support for further development.

Questaal is distinctive among electronic structure packages in offering, within a single code base and a shared input/data infrastructure:

- Both **atomic-sphere-approximation (ASA)** and **full-potential (FP)** LMTO implementations of DFT.
- An all-electron implementation of **quasiparticle self-consistent GW (QSGW)**, including an extension with ladder diagrams in the screened interaction (QSGŴ).
- An interface allowing **QSGW (or DFT) to be combined with DMFT**, for strongly correlated systems.
- ASA-based Green's-function methods supporting the **Coherent Potential Approximation (CPA)** for disorder, magnetic exchange coupling, and layer/transport (Landauer–Büttiker) calculations.
- A relativistic (scalar-Dirac, and in parts full-Dirac) treatment of the electronic structure.

**Primary reference:**
Dimitar Pashov, Swagata Acharya, Walter R. L. Lambrecht, Jerome Jackson, Kirill D. Belashchenko, Athanasios Chantis, Francois Jamet, Mark van Schilfgaarde, *"Questaal: a package of electronic structure methods based on the linear muffin-tin orbital technique,"* Computer Physics Communications **249**, 107065 (2020). https://doi.org/10.1016/j.cpc.2019.107065

Users are asked to cite this paper for any results generated with Questaal codes.

---

## 2. Theoretical Foundations

### 2.1 The LMTO method and its two lineages in Questaal

Questaal implements LMTO theory in two related but distinct forms:

- **ASA-LMTO** (the historical, "second-generation" Andersen formalism): spheres are inflated to fill space (overlapping so the interstitial volume vanishes), which is a geometry violation but yields a very compact, efficient basis — ideal for close-packed metals and for cases requiring Green's-function machinery (exchange couplings, CPA disorder, transport).
- **Full-potential LMTO (lmf)**: generalizes the LMTO envelope functions to smooth (numerically nonsingular) Hankel functions and augments them without the shape approximation of the ASA, giving a short-ranged, tight-binding-like basis suitable for accurate total energies, forces, and all-electron GW/DMFT.

Both trace to a common tight-binding/screened-envelope-function representation: envelope (Hankel/Bessel) functions are screened via structure constants to become short-ranged, following Andersen and Jepsen's "Explicit, First-Principles Tight-Binding Theory" (Phys. Rev. Lett. 53, 2571, 1984). This screening is what lets Questaal's ASA codes (`lm`, `lmgf`, `lmpg`) share a common structure-constant generator (`lmstr`).

### 2.2 One-body vs. many-body pictures

The methods paper frames Questaal's design around the distinction between:

- **One-body (static, Hermitian) potentials** — DFT/LDA, Hartree–Fock, and the QSGW potential all fall in this class; self-consistency means the independent variable that generates the potential (density, wavefunctions, or Green's function respectively) converges to what the potential itself produces.
- **Many-body (dynamical, non-Hermitian) self-energies** — GW (in its non-quasiparticlized, dynamical form) and DMFT construct a frequency-dependent self-energy Σ(ω), whose imaginary part encodes quasiparticle lifetime/damping.

QSGW is positioned as Questaal's central methodological innovation: rather than using DFT/LDA as a fixed starting point for a one-shot GW correction (the traditional G₀W₀ approach), QSGW iterates the Green's function to self-consistency by constructing a static, Hermitian, Hartree–Fock-like potential Σ₀ that best approximates the dynamical GW self-energy at each cycle. This removes starting-point dependence, satisfies a variational principle (Ismail-Beigi, J. Phys.: Condens. Matter 29, 385501, 2017), and makes the resulting errors systematic — enabling controlled improvements such as the addition of ladder diagrams to the screened interaction W (QSGŴ; Cunningham, Grüning, Pashov, van Schilfgaarde, Phys. Rev. B 108, 165104, 2023), which corrects W's tendency to over-screen within plain RPA.

### 2.3 DMFT integration

DMFT is treated as a complementary, site-local, non-perturbative many-body theory: it captures dynamical local correlations (e.g., in narrow d/f manifolds) that RPA-based GW handles poorly, in exchange for tractability being restricted to a local (and usually diagonal, to avoid the fermion sign problem) subspace of correlated orbitals. Because QSGW already removes most of the nonlocal/charge-fluctuation error left by LDA, combining **QSGW + DMFT** is presented as particularly effective — QSGW supplies a good nonlocal, static reference; DMFT supplies the missing local spin-fluctuation/Hubbard-band physics.

### 2.4 Relativity

Questaal solves a scalar-relativistic approximation to the Dirac equation (Koelling & Harmon, J. Phys. C 10, 3107, 1977) throughout most of the suite, with spin-orbit coupling added perturbatively; parts of the code (notably in the ASA Green's-function branch, `lmgf`) implement a fully relativistic Dirac formulation.

---

## 3. Executable Codes

Questaal is not a single monolithic binary but a family of executables sharing common input/data formats, connected through driver scripts. The methods paper and online documentation organize them as follows.

### 3.1 Input-preparation utilities

| Program | Function |
|---|---|
| `blm` | Generates the main input file (`ctrl.ext`) from structural information; central to most tutorials |
| `cif2init`, `cif2site` | Convert CIF crystallographic files (via the external `cif2cell` tool) into Questaal input |
| `poscar2init`, `poscar2site` | Convert VASP POSCAR files into Questaal input |
| `site2init` | Generates an `init` file from a Questaal `site` file |
| `lmchk` | Crystal-structure diagnostics: neighbor tables, bond angles, augmentation-sphere overlap checks, automatic sphere-radius determination, empty-sphere placement |
| `lmscell` | Supercell/superlattice maker; includes a layered-system inspection mode (`lmscell --stack`) |
| `lmplan` | Legacy analyzer for layered-system geometry (largely superseded by `lmscell --stack`) |
| `rdcmd` | Shell-like command reader using Questaal's own parser/programming language |
| `lmxbs` | Generates input for the `xbs` crystal-structure visualization tool |

### 3.2 One-body (DFT-level) electronic structure codes

| Program | Method | Notes |
|---|---|---|
| `lmf` | Full-potential LMTO DFT | The standard band-structure/self-consistent-density program; the platform onto which QSGW and DMFT are layered |
| `lm` | ASA-LMTO DFT (band/Hamiltonian formulation) | Requires companion program `lmstr` to generate screened structure constants |
| `lmgf` | ASA-LMTO DFT, Green's-function formulation | Enables magnetic exchange interactions, CPA (chemical/spin disorder), static susceptibilities |
| `lmpg` | ASA-LMTO Green's function for layered/semi-infinite systems | 2D periodic boundary conditions + real-space principal-layer treatment in the third dimension; Landauer–Büttiker transport, non-equilibrium capability |
| `tbe` | Empirical tight-binding | Self-consistent (important for polar systems), Hubbard-parameter support, GPU-parallelized builds, path-integral treatment of quantum nuclear motion |
| `lmmc` | Fast LDA-based molecular code | Undocumented/lightly documented |

### 3.3 Beyond-DFT (many-body) Green's-function codes

| Program | Function |
|---|---|
| `lmgw.sh` | Driver script orchestrating a full QSGW calculation as a sequence of smaller executables (exchange potential, polarization/RPA operator, correlation self-energy, etc.), synchronized with `lmf` |
| `lmfgwd` | Interface/driver generating the information the GW machinery needs from `lmf` |
| `lmfgws` | Post-processing analysis of dynamical self-energies / spectral functions, used after GW or DMFT (not CPA) runs |
| `lmfdmft` | Interface to external DMFT impurity solvers: generates hybridization/Weiss-field information and reads back self-energies to feed into `lmf` |
| `lmdmft` | Constructs the nonlocal (site-diagonal, orbital-nonlocal) self-energy object used in DMFT embedding |

**Important architectural point:** the DMFT impurity solvers themselves are *not* distributed as part of Questaal. Questaal supplies the DFT/QSGW host, the projection to a correlated subspace, and the self-consistency bookkeeping, but the impurity problem must be solved by an external package (e.g., a CT-QMC solver such as those from the TRIQS ecosystem) compiled separately and interfaced via `lmfdmft`.

### 3.4 Post-processing and supporting utilities

| Program | Function |
|---|---|
| `lmdos` | Partial density-of-states generation from `lmf`, `lm`, or `tbe` output |
| `fplot` | Gnuplot-like general-purpose plotting utility |
| `plbnds` | Reformats band-structure data (from `lmf`/`lm`/`tbe`) or spectral-function data (from `lmfgws`) for plotting |
| `pldos` | Organizes `lmdos` output for plotting |
| `bandedge` | Script to locate band extrema and extract effective masses |
| `spectral` | Converts raw GW dynamical self-energy output into the `se` format read by `lmfgws` |
| `s2s5` | Converts a binary QSGW self-energy file to HDF5 |
| `mcx` | Matrix calculator |
| `pfit`, `dval` | Least-squares fitting / generic calculator utilities (undocumented) |
| `map-results-irr-to-fbz` | Post-processes magnetic-susceptibility results into exchange interactions |

### 3.5 Shared input infrastructure

Nearly all executables read a common, tree-structured, largely format-free input file `ctrl.ext`, parsed through a built-in preprocessor with variable declarations, conditionals, and algebraic expressions — effectively a small domain-specific input language. This allows a single `ctrl` file to range from a minimal tutorial example to a large database-like file spanning many materials and run modes. The `GWinput` file separately controls GW-specific parameters (k-mesh for polarizability, basis for the mixed product basis, energy contour, etc.).

---

## 4. DFT Layer: Full-Potential vs. ASA

### 4.1 Full-potential LMTO (`lmf`)

`lmf` is Questaal's principal DFT engine and the foundation for GW and DMFT. Its basis set is built from **smooth (nonsingular) Hankel functions** (Bott, Methfessel, Krabs, Schmid, J. Math. Phys. 39, 3393, 1998), which behave as ordinary Hankel functions outside a smoothing radius but remain finite and well-behaved at the origin, allowing a numerically robust, short-ranged, tight-binding-style envelope basis without the ASA's space-filling-sphere approximation.

Notable capabilities documented for `lmf`:
- Basis-set augmentation with local orbitals (semicore states) and, optionally, extra plane waves (`PWMODE`) to systematically converge difficult cases (e.g., open structures, or where LMTO alone leaves residual basis-incompleteness error) — the PMT (LMTO + APW) fusion method of Kotani & van Schilfgaarde (Phys. Rev. B81, 125117, 2010).
- Self-consistent total-energy and force calculations, enabling structural relaxation ("molecular statics").
- Noncollinear magnetism and perturbative or self-consistent spin-orbit coupling.
- Generation of energy bands, partial DOS, charge densities, and optical/dielectric properties.

### 4.2 ASA suite (`lm`, `lmgf`, `lmpg`)

The ASA codes trade the exact interstitial geometry for overlapping, space-filling spheres — a simplification pioneered by Andersen in the 1970s originally to treat transition metals. The methods paper and ASA-overview documentation describe:

- **Continuous principal quantum numbers** and floating linearization energies used to parameterize the energy-dependent phase shift into an energy-independent Hamiltonian (the defining "linearization" step of all LMTO variants).
- The **combined correction** term, which recovers some of the accuracy lost by ignoring the true (non-spherical, non-space-filling) potential shape.
- **Downfolding**, allowing high-lying or chemically inert orbitals to be eliminated from the active Hilbert space while retaining their effect via renormalized structure constants — useful for reducing basis size in complex structures.
- Automatic/adaptive **sphere-radius determination** and **empty-sphere insertion**, needed when treating open (non-close-packed) structures within the ASA's space-filling constraint.

Three ASA-level programs exist:

1. **`lm`** — the Hamiltonian (band) formulation; self-consistent densities and energy bands for crystals.
2. **`lmgf`** — a Green's-function formulation nearly equivalent to the `lm` Hamiltonian without the combined-correction term, built on an analytic approximation to KKR multiple-scattering theory. Distinct capabilities relative to `lm` include magnetic exchange-interaction calculations (mapping DFT total-energy differences onto a Heisenberg model via the magnetic force theorem / transverse susceptibility), spin-spin/spin-orbit/orbit-orbit susceptibilities, and the **Coherent Potential Approximation (CPA)** for chemically or spin-disordered (including finite-temperature paramagnetic) systems.
3. **`lmpg`** — the layer-Green's-function analogue of `lmgf`, using 2D periodicity in-plane and a real-space principal-layer technique along the third (surface-normal or transport) direction; supports Landauer–Büttiker transport calculations and non-equilibrium extensions (used, e.g., in superconducting trilayer/Andreev-level studies).

**Trade-off summary:** ASA codes are cheaper and give access to Green's-function-specific physics (CPA, exchange couplings, transport) not available (or not as naturally available) in `lmf`, but sacrifice full-potential accuracy; `lmf` is the higher-fidelity, general-purpose route and the required entry point for GW/DMFT.

---

## 5. Many-Body Perturbation Theory: GW and QSGW

### 5.1 Historical/algorithmic lineage

Questaal's GW implementation descends from an all-electron, mixed-product-basis GW method (Kotani & van Schilfgaarde, Solid State Commun. 121, 461, 2002), building on the foundational GW formalism of Hedin (Phys. Rev. 139, A796, 1965) and the first practical semiconductor/insulator GW calculations of Hybertsen and Louie (Phys. Rev. B34, 5390, 1986).

### 5.2 One-shot GW vs. QSGW

- **One-shot (G₀W₀) GW**: the traditional perturbative correction to an LDA/DFT starting point; available in Questaal but subject to well-documented starting-point sensitivity (van Schilfgaarde, Kotani, Faleev, Phys. Rev. B74, 245125, 2006).
- **Quasiparticle self-consistent GW (QSGW)**: iterates the Green's function G and screened interaction W to self-consistency by constructing an optimal static, Hermitian one-body potential from the dynamical self-energy at each cycle (Faleev, van Schilfgaarde, Kotani, Phys. Rev. Lett. 93, 126406, 2004; van Schilfgaarde, Kotani, Faleev, Phys. Rev. Lett. 96, 226402, 2006; algorithmic basis in Kotani, van Schilfgaarde, Faleev, Phys. Rev. B76, 165106, 2007). This is Questaal's flagship method and its principal point of differentiation from most other electronic-structure packages.
- **QSGŴ (with ladder diagrams)**: extends QSGW by including electron-hole ladder diagrams (via a Bethe–Salpeter-equation-derived vertex) in the construction of the screened interaction W, correcting RPA's tendency to over-screen — identified as the dominant remaining systematic error in plain QSGW (Cunningham, Grüning, Pashov, van Schilfgaarde, Phys. Rev. B 108, 165104, 2023).

### 5.3 Workflow and execution model

QSGW calculations are orchestrated by the shell driver `lmgw.sh`, which calls a sequence of smaller executables computing the bare exchange (Fock) potential, the RPA polarizability, the screened correlation self-energy, etc., synchronizing at each step with `lmf`, which continues to manage the underlying one-body (density/potential) machinery. `lmfgwd` serves as the driver/interface generating GW-specific data from an `lmf` calculation. Configuring `lmgw.sh` for HPC/cluster environments is documented separately, reflecting the significant computational cost of GW relative to DFT.

### 5.4 Derived spectral/response-function capabilities

Built on top of the GW machinery, Questaal implements:
- **RPA and RPA+BSE dielectric functions** — optical/dielectric response, with Bethe–Salpeter-equation corrections for excitonic effects.
- **Transverse magnetic susceptibility** calculations from GW.
- **Phonon-related properties** derived from the GW self-energy framework.
- **Photoemission modeling** using GW spectral functions.
- Analysis tools (`lmfgws`, `spectral`, `s2s5`) for post-processing the dynamical self-energy into spectral functions, density of states, and format conversions (including HDF5 export).

### 5.5 Total-energy status

An RPA-based total-energy functional (Luttinger–Ward) was implemented in Questaal's code antecedent as early as 2002 (Miyake, Aryasetiawan, Kotani, van Schilfgaarde, Usuda, Terakura, Phys. Rev. B66, 245103, 2002), but basis-set-incompleteness issues limited its practicality within the LMTO framework; RPA total energies have since become more tractable in large plane-wave codes (at high computational cost), and basis-set-corrected RPA total energies are described as an active area of ongoing Questaal development.

---

## 6. DMFT Integration (QSGW+DMFT / DFT+DMFT)

### 6.1 Conceptual framing

DMFT is treated in Questaal as a site-local, non-perturbative complement to the (nonlocal, but static/perturbative) QSGW self-energy. The DMFT self-energy Σ_DMFT is constructed to be local in real space but retains full frequency dependence and (in principle) full orbital-matrix structure within a chosen correlated subspace — typically the d- or f-orbital manifold at a given site — although in practice the matrix is usually restricted to its diagonal to avoid a severe fermion sign problem in the impurity solver.

### 6.2 Architecture

- **`lmfdmft`** acts as the interface layer: it projects the QSGW (or DFT) Hamiltonian/Green's function onto the correlated subspace, generates the input (hybridization function / Weiss field) required by an external impurity solver, and later reads the resulting self-energy back in to update the lattice Green's function.
- **`lmdmft`** constructs the nonlocal-in-site, orbital-local self-energy object used within the embedding scheme.
- **The impurity solver itself is external** to Questaal. Community practice (per external literature, e.g. the ComDMFT paper) references solvers from the **EDMFTF**/TRIQS ecosystem; Questaal's own DMFT tutorial series is built around **TRIQS**-based continuous-time quantum Monte Carlo (CT-QMC) solvers.

### 6.3 Documented workflow (per Questaal's DMFT tutorial series)

The DMFT tutorial pipeline (as organized on questaal.org) covers, in sequence: an introduction to the method; setting up and running the DMFT self-consistency loop; handling input/parameter subtleties; charge and static-magnetic self-consistency contributions; the density loop; the Maximum Entropy analytic-continuation method (needed to convert imaginary-axis CT-QMC results to real-frequency spectral functions); the self-energy loop and "dynamical double counting" (the correction needed to avoid double-counting correlation effects already partly captured by the GW/DFT host); and final spectral-function analysis. The documentation explicitly notes this pipeline is comparatively complex to set up and that parts of the tutorial have not been updated for some time — an important caveat for prospective users.

### 6.4 Downstream physics enabled

DMFT-derived response functions documented for the suite include dielectric functions, transverse magnetic susceptibility, and pairing susceptibility computed on top of DMFT self-energies — relevant to spectroscopy (photoemission modeling) and to magnetic/superconducting instability analysis in correlated materials (e.g., iron pnictides, transition-metal oxides, actinide compounds — the classic proving grounds for DFT+DMFT-type methods).

### 6.5 Position relative to the wider DFT+DMFT ecosystem

Questaal is one of a small number of packages offering a **GW+DMFT** (rather than plain LDA/DFT+DMFT) combination — placing it alongside packages such as ComDMFT (which implements ab initio LQSGW+DMFT and, in newer versions, fully self-consistent GW+EDMFT). The QSGW+DMFT combination is argued in the primary literature to be particularly well-matched: QSGW supplies an accurate, self-consistent, nonlocal (largely charge/exchange) starting point, while DMFT supplies the local dynamical spin-fluctuation and Hubbard-band physics that QSGW's RPA-level treatment misses.

---

## 7. Additional Distinguishing Capabilities

- **LDA+U**: available as a simpler, static correction to correlated-orbital physics, usable independently of GW/DMFT or as an auxiliary component.
- **Coherent Potential Approximation (CPA)** (via `lmgf`): mean-field treatment of substitutional chemical disorder and/or spin disorder (e.g., paramagnetic states above a magnetic ordering temperature), a capability essentially unique to the Green's-function branch of the ASA suite.
- **Magnetic exchange coupling** (via `lmgf`): extraction of Heisenberg-model exchange parameters (Jᵢⱼ) from the transverse magnetic susceptibility / magnetic force theorem, feeding into downstream spin-dynamics or Curie-temperature estimates.
- **Transport** (via `lmpg`): Landauer–Büttiker conductance/transmission calculations for layered and multi-terminal structures, including non-equilibrium and superconducting (Andreev-level) extensions.
- **Empirical tight-binding (`tbe`)**: a distinct, non-ab-initio pathway using Slater–Koster parameterization, notable for supporting self-consistency (important for polar/ionic systems), Hubbard-U-like corrections, GPU-accelerated builds, and a path-integral treatment of quantum nuclear motion (zero-point/quantum nuclear effects).
- **Wannier-90 interface**: documented support for exporting Questaal Hamiltonians to the Wannier90 ecosystem for maximally localized Wannier function construction.
- **Response-function breadth**: the developers state Questaal has (to their knowledge) a unique combination of spin, charge, and pairing response functions computable at both the GW and DMFT levels within a single, internally consistent framework.

---

## 8. Numerical/Software Design Notes

- **Basis functions**: smooth (nonsingular) Hankel functions underlie the full-potential basis, chosen for numerical robustness at the origin while retaining short-rangedness in real space.
- **Structure constants**: precomputed and stored (via `lmstr` for the ASA codes) as screened, short-ranged tight-binding structure constants, following the Andersen–Jepsen screening transformation.
- **Relativity**: scalar-relativistic Dirac treatment (Koelling–Harmon) as the default across the suite, with a fully relativistic Dirac formulation available in parts of the ASA/Green's-function branch, and perturbative or self-consistent spin-orbit coupling elsewhere.
- **Brillouin-zone integration, spin/spin-orbit handling, and rotation conventions** are documented as distinct numerics topics, reflecting the suite's long development history and the need to support noncollinear magnetism consistently across DFT, GW, and DMFT layers.
- **Editors**: a family of embedded binary/data-file editors is built into the codes for manipulating restart files, self-energy files, and structural data directly, rather than requiring separate tooling.

---

## 9. Practical Considerations

- **Licensing/distribution**: Questaal is open source; the package is obtained via a download link on questaal.org, with the CPC paper serving as the mandatory citation for any published results.
- **Build complexity**: as a suite of many interdependent executables (DFT, GW driver scripts, DMFT interfaces, and numerous pre-/post-processing utilities) sharing common data formats, Questaal has a nontrivial build and configuration process, particularly for GW on HPC systems (a dedicated tutorial addresses this) and for DMFT, which additionally requires compiling and interfacing an external impurity solver not distributed with Questaal itself.
- **Documentation state**: core DFT, ASA, and GW documentation and tutorials are actively maintained; the DMFT tutorial pipeline is explicitly flagged by the developers as more complex to execute and, in places, not recently updated — worth budgeting extra time for when setting up DMFT or QSGW+DMFT workflows.
- **Community/training**: the suite has been the subject of dedicated CECAM/Daresbury training schools ("Questaal schools") covering DFT/ASA basics, QSGW, and DMFT (including hands-on TRIQS-based CT-QMC exercises), reflecting an active user/developer community centered on King's College London, Case Western Reserve University, University of Nebraska–Lincoln, and STFC Daresbury Laboratory.

---

## 10. Summary Assessment

Questaal occupies a distinctive niche among electronic-structure packages: it is simultaneously (a) a mature, all-electron ASA/full-potential LMTO DFT code with capabilities (CPA disorder, layer transport, magnetic exchange couplings) more commonly associated with dedicated Green's-function/KKR packages, and (b) one of the few packages offering production-grade **quasiparticle self-consistent GW**, with the further, comparatively rare extension to **ladder-diagram-corrected QSGŴ** and to **QSGW+DMFT** for strongly correlated materials. Its main practical costs are the complexity inherent in coordinating many separate executables/interfaces (particularly for GW and DMFT workflows) and a reliance on external, separately compiled DMFT impurity solvers. For users specifically interested in correlated-electron spectroscopy — combining accurate nonlocal exchange/correlation physics (via QSGW) with local dynamical correlation (via DMFT) — Questaal is one of a very small number of open-source codes offering this combination natively.

---

## Key References

1. D. Pashov, S. Acharya, W. R. L. Lambrecht, J. Jackson, K. D. Belashchenko, A. Chantis, F. Jamet, M. van Schilfgaarde, *Questaal: a package of electronic structure methods based on the linear muffin-tin orbital technique*, Comput. Phys. Commun. **249**, 107065 (2020).
2. O. K. Andersen, *Linear methods in band theory*, Phys. Rev. B **12**, 3060 (1975).
3. O. K. Andersen, O. Jepsen, *Explicit, First-Principles Tight-Binding Theory*, Phys. Rev. Lett. **53**, 2571 (1984).
4. M. Methfessel, M. van Schilfgaarde, R. A. Casali, *A full-potential LMTO method based on smooth Hankel functions*, in *Electronic Structure and Physical Properties of Solids*, Lecture Notes in Physics **535**, 114 (Springer, 2000).
5. T. Kotani, M. van Schilfgaarde, *All-electron GW approximation with the mixed basis expansion based on the full-potential LMTO method*, Solid State Commun. **121**, 461 (2002).
6. S. V. Faleev, M. van Schilfgaarde, T. Kotani, *All-electron self-consistent GW approximation: Application to Si, MnO, and NiO*, Phys. Rev. Lett. **93**, 126406 (2004).
7. M. van Schilfgaarde, T. Kotani, S. V. Faleev, *Quasiparticle self-consistent GW theory*, Phys. Rev. Lett. **96**, 226402 (2006).
8. T. Kotani, M. van Schilfgaarde, S. V. Faleev, *Quasiparticle self-consistent GW method: a basis for the independent-particle approximation*, Phys. Rev. B **76**, 165106 (2007).
9. T. Kotani, M. van Schilfgaarde, *A fusion of the LAPW and the LMTO methods: the augmented plane wave plus muffin-tin orbital (PMT) method*, Phys. Rev. B **81**, 125117 (2010).
10. B. Cunningham, M. Grüning, D. Pashov, M. van Schilfgaarde, *QSGŴ: Quasiparticle Self-consistent GW with ladder diagrams in W*, Phys. Rev. B **108**, 165104 (2023).
11. S. Ismail-Beigi, *Justifying quasiparticle self-consistent schemes via gradient optimization in Baym–Kadanoff theory*, J. Phys.: Condens. Matter **29**, 385501 (2017).
12. I. Turek et al., *Electronic structure of disordered alloys, surfaces and interfaces* (Kluwer, Boston, 1996).
13. D. D. Koelling, B. N. Harmon, *A technique for relativistic spin-polarised calculations*, J. Phys. C: Solid State Phys. **10**, 3107 (1977).
14. Official documentation and tutorials: https://www.questaal.org/docs/ and https://www.questaal.org/tutorial/

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the Questaal (LMTO/ASA) 	Suite of full-potential and atomic-sphere-approximation LMTO codes for electronic structure, including DFT+DMFT and GW. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
