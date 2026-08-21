# SPR-KKR: A Comprehensive Review

**Fully Relativistic Spin-Polarized KKR Green's-Function DFT Code for Solids, Surfaces, and Disordered Systems**

---

## 1. Overview

SPR-KKR (Spin-Polarized Relativistic Korringa–Kohn–Rostoker) is an all-electron, Green's-function-based *ab initio* electronic-structure package developed and maintained principally by the group of Hubert Ebert (Ludwig-Maximilians-Universität München), together with long-standing collaborators (D. Ködderitzsch, J. Minár, S. Mankovsky, S. Bornemann, S. Polesya, and others). Rather than solving for Bloch eigenstates directly, SPR-KKR formulates the electronic structure problem in terms of the single-particle Green's function via multiple-scattering theory (MST). This Green's-function foundation is the source of nearly all of the package's distinguishing capabilities.

- **Method family:** Korringa–Kohn–Rostoker (KKR) multiple-scattering theory, in its full-potential and atomic-sphere-approximation (ASA) variants.
- **Relativity:** Supports non-relativistic, scalar-relativistic, and fully relativistic (four-component Dirac) modes, with spin–orbit coupling treated on an equal footing with exchange splitting in the fully relativistic mode.
- **Disorder treatment:** Native, non-perturbative treatment of substitutional chemical disorder via the Coherent Potential Approximation (CPA), avoiding the need for large supercells.
- **Primary language/platform:** Fortran core with Tcl-based control scripts; historically Linux-only; a graphical front end (XBAND) and, more recently, Python tooling have been developed around it.
- **License model:** Distributed under a group/site licence upon registration (not fully open-source in the conventional sense), from the Munich group's website.

---

## 2. Theoretical Foundation

### 2.1 Multiple-Scattering / KKR-Green's-Function Formalism

SPR-KKR builds the single-particle retarded Green's function of the solid,

$$G(\vec r,\vec r',E) = \sum_{\Lambda_1\Lambda_2} Z^n_{\Lambda_1}(\vec r,E)\,\tau^{nn'}_{\Lambda_1\Lambda_2}(E)\,Z^{n'\times}_{\Lambda_2}(\vec r',E) - \left[Z^n_{\Lambda_1}(\vec r,E) J^{n\times}_{\Lambda_1}(\vec r',E)\Theta(r'-r) + J^n_{\Lambda_1}(\vec r,E) Z^{n\times}_{\Lambda_1}(\vec r',E)\Theta(r-r')\right]\delta_{nn'}\, ,$$

where $Z^n_\Lambda$ and $J^n_\Lambda$ are the regular and irregular solutions of the single-site Dirac equation on site $n$, and $\tau^{nn'}_{\Lambda\Lambda'}(E)$ is the scattering-path operator connecting sites $n$ and $n'$, obtained from the site $t$-matrices and the KKR structure constants via the multiple-scattering (Dyson-like) matrix equation. In the fully relativistic mode the composite index $\Lambda=(\kappa,\mu)$ combines the relativistic spin–orbit quantum number $\kappa$ and the magnetic quantum number $\mu$.

The underlying single-site problem is the Dirac equation of relativistic spin-density-functional theory,

$$\mathcal H_D = -ic\,\vec{\boldsymbol\alpha}\cdot\vec\nabla + \tfrac12 c^2(\boldsymbol\beta-1) + V(\vec r) + \beta\,\vec{\boldsymbol\sigma}\cdot\vec B(\vec r)\, ,$$

with $V(\vec{r})$ and $\vec{B}(\vec{r})$ the effective electrostatic and exchange–correlation magnetic-field potentials of (spin-)density functional theory.

Because the formalism is built directly around $G$ rather than around eigenstates, SPR-KKR gains several structural advantages:

- **Direct access to spectroscopic and response quantities** (DOS, Bloch spectral functions, susceptibilities, conductivities) without an explicit sum over bands.
- **Embedding via the Dyson equation:** local perturbations of a host system (impurities, clusters, surfaces, interfaces) can be treated by updating only the relevant scattering matrices, without artificial periodic supercells.
- **Natural treatment of substitutional disorder** through configurational (ensemble) averaging of $\tau$, most commonly via the CPA.

### 2.2 Levels of Relativistic Treatment

- **Non-relativistic** (Schrödinger-based).
- **Scalar-relativistic**: mass–velocity and Darwin terms retained, spin–orbit coupling dropped; compatible with non-collinear magnetism and spin spirals.
- **Fully relativistic**: four-component Dirac formalism with spin–orbit coupling and exchange splitting treated on equal footing, breaking rotational spin symmetry so that only the total angular momentum projection $J_z$ remains a good quantum number in general magnetic configurations. This mode is essential for magnetocrystalline anisotropy, Dzyaloshinskii–Moriya interactions (DMI), X-ray magnetic circular/linear dichroism, spin–orbit torques, and topological/altermagnetic band-structure effects.

### 2.3 Potential Treatment

- **Atomic Sphere Approximation (ASA):** space-filling, (optionally overlapping) spherical potentials — computationally efficient, historically the default mode.
- **Full-Potential (FP) mode:** shape-unrestricted cell potentials for higher accuracy, at higher computational cost, important for open/low-symmetry structures and total-energy-sensitive properties.

### 2.4 Exchange–Correlation and Beyond-LSDA Treatments

- Local (spin) density approximation (LSDA), e.g., the Vosko–Wilk–Nusair parametrization.
- Generalized Gradient Approximation (GGA), e.g., Perdew–Burke–Ernzerhof (PBE).
- **LSDA+U / rotationally invariant LDA+U**, incorporated into the fully relativistic SPR-KKR formalism to treat on-site Coulomb correlation (e.g., in f-electron systems), combinable with CPA for correlated disordered alloys.
- **KKR+DMFT (Dynamical Mean-Field Theory):** the local, site-projected Green's function from KKR is used to construct the DMFT self-energy problem (e.g., via the Spin-polarized T-matrix Fluctuation-Exchange, SPTF, solver, or perturbative FLEX-type solvers), giving a relativistic LSDA+DMFT scheme for strongly correlated 3d/4f/5f systems.
- **Bogoliubov–de Gennes (BdG) extension:** self-consistent treatment of superconducting pairing within the KKR-GF framework (layered/heterostructure geometries), enabling access to gap anisotropy and proximity effects from first principles.

### 2.5 Disorder: Coherent Potential Approximation (CPA)

Chemical (substitutional) disorder on one or more sublattices is treated by the single-site CPA (Soven; Győrffy formulation for KKR), self-consistently determining an effective medium such that a single "average" scatterer causes zero further scattering on average. This is combined seamlessly with the relativistic and LSDA+U extensions, and is central to alloy, dilute-magnetic-semiconductor, and finite-temperature-magnetism (via the "disordered local moment", DLM, picture) calculations.

### 2.6 Two-Dimensional / Layered (Screened / TB-KKR) Techniques

A tight-binding/screened-KKR reformulation (building on the SPR-TB-KKR developments, incorporating work originally due to R. Zeller and the Jülich group, later merged into the Munich package) allows genuinely semi-infinite 2-D-periodic geometries — true surfaces, interfaces, and multilayers — to be treated without artificial 3-D periodic repetition, alongside the layer-KKR (matrix-inversion in reciprocal 2-D space) approach used for photoemission and diffraction calculations.

---

## 3. Core Capabilities and Property Modules

SPR-KKR is organized as a suite of task-driven executables (e.g., `kkrscf`, `kkrgen`, `kkrspec`, etc., historically driven by Tcl scripting and an input "task" file) sharing a common self-consistent potential. Principal capabilities include:

### 3.1 Ground-State Electronic Structure
- Self-consistent charge and spin densities for 3-D periodic crystals, disordered alloys (via CPA), and low-dimensional/embedded systems.
- Band structures, densities of states (total, site-, and orbital-resolved), Fermi surfaces.
- Total energies, magnetic moments (spin and orbital), and hyperfine fields.
- Non-collinear and spin-spiral magnetic structures (scalar-relativistic and fully relativistic modes).

### 3.2 Disordered and Low-Dimensional Systems
- Random substitutional alloys and partial/fractional site occupations via CPA.
- Finite-temperature magnetism through the disordered local moment (DLM) approach combined with CPA.
- Surfaces, interfaces, thin films, and embedded clusters via Green's-function embedding / screened KKR, circumventing large supercells.

### 3.3 Magnetic Interactions and Micromagnetics
- Extraction of interatomic (isotropic) exchange coupling constants $J_{ij}$ (Lichtenstein-type mapping onto a classical Heisenberg model).
- Relativistic extensions yielding the Dzyaloshinskii–Moriya interaction (DMI) tensor, exchange-stiffness tensor, and spin-wave stiffness directly from the fully relativistic KKR-GF technique.
- Magnetocrystalline anisotropy energies (via the magnetic force theorem or total-energy differences).
- Interfacing with Monte Carlo and atomistic spin-dynamics frameworks (e.g., estimation of ordering/Curie temperatures) using the computed exchange parameters.

### 3.4 Spectroscopies
- X-ray absorption spectroscopy (XAS) and X-ray magnetic circular/linear dichroism (XMCD/XMLD).
- Photoemission: one-step-model angle-resolved photoemission spectroscopy (ARPES), including matrix-element, final-state, and surface effects; extensions to hard X-ray ARPES (HARPES) and photoelectron diffraction (PED/XPD) via a layer-KKR (k-space) formulation valid over a wide energy range.
- Compton profiles and magnetic Compton scattering.
- Core-level and valence-band photoemission with disorder ("alloy analogy" model) and finite-temperature (vibrational/spin-disorder) effects.

### 3.5 Response Functions and Transport
- Linear-response quantities: susceptibilities, conductivity tensors.
- Electronic and spin transport, including anomalous and spin Hall effects, within the Kubo–Greenwood/Kubo–Středa formalism implemented on top of the KKR-GF.
- Bloch spectral functions $A(\vec k,E)$ as the CPA-averaged, $\vec k$-resolved analogue of the DOS, used to visualize band "smearing" from disorder.

### 3.6 Correlated-Electron Extensions
- LSDA+U (rotationally invariant, relativistic) for localized correlated orbitals (notably actinide/lanthanide f-states).
- LSDA+DMFT with perturbative (SPTF/FLEX-type) impurity solvers for dynamical (energy-dependent, complex) self-energies, giving access to correlated spectral functions, quasiparticle renormalization, and satellite structures.

### 3.7 Superconductivity
- Fully relativistic Bogoliubov–de Gennes extension of the KKR-GF method for layered superconducting heterostructures, giving access to gap anisotropy, proximity-induced pairing, and (in principle) Yu–Shiba–Rusinov-type impurity states.

---

## 4. Numerical / Technical Characteristics

| Aspect | Description |
|---|---|
| Basis / representation | Multiple-scattering (Green's function), not a conventional wavefunction basis; angular-momentum expansion truncated at $l_\text{max}$ (commonly 2–4) |
| Potential shape | ASA (fast) or full potential (FP, shape-unrestricted) |
| Relativity | Non-, scalar-, or fully relativistic (Dirac) |
| Disorder | CPA (single-site), combinable with LSDA+U and DMFT |
| BZ integration | Regular $k$-mesh (e.g., Monkhorst–Pack-like grids); mesh density is a user-set convergence parameter, often tens of $k$-points per axis for exchange-parameter or spectroscopy runs |
| Lattice sums | Structure constants typically evaluated via Ewald summation |
| Exchange–correlation | LSDA (Vosko–Wilk–Nusair) and GGA (PBE), among others |
| Parallelism | MPI support for computationally demanding (e.g., high-throughput, DMFT, transport) calculations |
| Platform | Historically Linux; core written in Fortran with Tcl-based job/task control |
| Distribution | Registration-based site/group licence from the Munich group (not a conventional open-source repository) |

---

## 5. Ecosystem, Interfaces, and Tooling

- **XBAND:** the package's traditional graphical user interface for setting up and analysing calculations.
- **ASE2SPRKKR:** a modern, actively developed Python interface integrating SPR-KKR into the Atomic Simulation Environment (ASE). It extends ASE's `Atoms` object to represent fractional/CPA site occupations, automates input generation and validation, parses SPR-KKR's output files (`.bsf`, `.dos`, `.spc`, `.jxc`, `.xas`, etc.) with built-in visualization, and interoperates with external structure-relaxation codes (e.g., GPAW, VASP) so that a relaxed geometry can be handed directly to SPR-KKR for spectroscopic post-processing — enabling high-throughput and multi-code workflows.
- **CECAM/community activity:** SPR-KKR is regularly featured alongside other KKR-Green's-function codes (e.g., JuKKR from Jülich, KKRnano, EMTO) at CECAM workshops on Green's-function methods for spectroscopy, transport, and magnetism, reflecting its role as one of the principal general-purpose KKR-GF platforms in the community.
- **Version history reflected in the literature:** citations trace versions from ~6.3 through 7.7 and 8.x (e.g., version 8.5 cited in 2022-era work), indicating continuous development.

---

## 6. Typical Application Domains

- Bulk metals, semiconductors, and compounds (band structures, DOS, moments).
- Random substitutional alloys (e.g., Heusler alloys, permanent-magnet compounds, dilute magnetic semiconductors).
- Magnetic thin films, surfaces, and multilayers (e.g., ultrathin films on oxide or metal substrates, spin-orbit-torque heterostructures).
- Finite-temperature magnetism and magnetic phase diagrams via DLM + CPA and derived exchange parameters feeding Monte Carlo/atomistic spin dynamics.
- Topological materials and altermagnets, where fully relativistic spin–orbit treatment is essential to resolve Kramers-degeneracy lifting mechanisms.
- Actinide/lanthanide systems (via LSDA+U and DMFT) — historically including early relativistic spin-polarized studies of δ-plutonium.
- Interpretation of experimental XAS/XMCD, ARPES/HARPES, photoelectron diffraction, and Compton-scattering data via direct one-step/Green's-function spectral simulations.
- Superconducting heterostructures and proximity effects.

---

## 7. Strengths and Limitations

**Strengths**
- Native, non-supercell treatment of chemical disorder (CPA) integrated with relativistic and correlated-electron extensions.
- Full four-component relativistic (Dirac) treatment on equal footing with exchange splitting — well suited to spin–orbit-driven phenomena (MAE, DMI, XMCD, spin Hall/anomalous Hall effects, altermagnetism).
- Direct Green's-function access to a very broad palette of spectroscopic and response properties without recomputing eigenstates.
- Mature embedding techniques for surfaces, interfaces, and impurities.
- Long, continuously maintained development history with tight coupling between method papers and code releases.

**Limitations**
- All-electron, multiple-scattering formalism is less naturally suited to structural relaxation/large-scale total-energy minimization than plane-wave or LAPW codes; structural optimization is often outsourced to other codes (e.g., VASP, GPAW) with SPR-KKR used for the electronic/spectroscopic post-processing (as facilitated by ASE2SPRKKR).
- CPA is a single-site (mean-field-like) approximation to disorder and does not capture short-range order/clustering effects without extensions (e.g., non-local CPA variants elsewhere in the KKR literature).
- Historically Linux/Fortran/Tcl-centric distribution model with registration-gated access, rather than a fully open, pip/conda-installable package (though the Python ecosystem around it is improving via ASE2SPRKKR).
- Angular-momentum truncation ($l_\text{max}$) and ASA-vs-FP choice require convergence testing, as with other KKR-family codes.

---

## 8. Key Publications on the Theory

### 8.1 Foundational and Method-Overview Papers

1. J. Korringa, *On the calculation of the energy of a Bloch wave in a metal*, **Physica** 13, 392 (1947).
2. W. Kohn and N. Rostoker, *Solution of the Schrödinger equation in periodic lattices with an application to metallic lithium*, **Phys. Rev.** 94, 1111 (1954).
3. P. Soven, *Coherent-Potential Model of Substitutional Disordered Alloys*, **Phys. Rev.** 156, 809 (1967).
4. G. M. Stocks, W. M. Temmerman, and B. L. Győrffy, *Complete solution of the Korringa–Kohn–Rostoker coherent-potential-approximation equations: Cu–Ni alloys*, **Phys. Rev. Lett.** 41, 339 (1978).
5. B. L. Győrffy, *Coherent-potential approximation for a nonoverlapping-muffin-tin-potential model of random substitutional alloys*, **Phys. Rev. B** 5, 2382 (1972).
6. H. Ebert, H. Freyer, A. Vernes, and G.-Y. Guo, *Correlated four-component band structure calculations of transition metal compounds*, and related early relativistic-KKR development papers by H. Ebert and coworkers (1980s–1990s).
7. G. Schadler, P. Weinberger, A. M. Boring, and R. C. Albers, *The relativistic spin-polarized KKR-Green's function — applications to the band structure of plutonium*, early relativistic spin-polarized KKR-GF papers (mid-1980s), foundational to the SPR-KKR approach.

### 8.2 The Central SPR-KKR Method/Review Article

8. **H. Ebert, D. Ködderitzsch, and J. Minár**, *Calculating condensed matter properties using the KKR-Green's function method — recent developments and applications*, **Reports on Progress in Physics** 74, 096501 (2011). DOI: 10.1088/0034-4885/74/9/096501.
   — This is the principal, most widely cited review of the modern KKR-GF method as implemented in SPR-KKR: it outlines the multiple-scattering/Green's-function formalism, relativistic and CPA extensions, and surveys applications across spectroscopy, transport, magnetism, and correlated-electron physics.

### 8.3 Extensions to the SPR-KKR Formalism

9. H. Ebert and S. Mankovsky, *Anisotropic exchange coupling in diluted magnetic semiconductors: Ab initio spin-density functional theory*, **Phys. Rev. B** 79, 045209 (2009).
10. A. B. Shick, V. Drchal, and J. Kudrnovský/H. Ebert et al., *Incorporation of the rotationally invariant LDA+U scheme into the SPR-KKR formalism: application to disordered alloys*, **J. Phys. Chem. Solids** (2003) / related SPR-KKR + LDA+U papers.
11. J. Minár, L. Chioncel, A. Perlov, H. Ebert, M. I. Katsnelson, and A. I. Lichtenstein, *Multiple-scattering formalism for correlated systems: A KKR-DMFT approach*, **Phys. Rev. B** 72, 045125 (2005).
12. J. Minár, H. Ebert, and L. Chioncel, *Correlation effects in ground-state properties probed by DFT+DMFT*, review-type article on the combined LSDA+DMFT method as implemented within KKR-based codes.
13. J. Braun, J. Minár, and H. Ebert, *Correlation, temperature and disorder: Recent developments in the one-step description of angle-resolved photoemission*, **Physics Reports** 740, 1–34 (2018).
14. S. Mankovsky, S. Polesya, and H. Ebert, *Spin-wave stiffness and micromagnetic exchange interactions expressed by means of the KKR Green function approach*, **Phys. Rev. B** (2018/2019); arXiv:1810.13175.
15. S. Mankovsky and H. Ebert, *First-principles calculation of the parameters used by atomistic magnetic simulations*, **Electron. Struct.** 4, 034004 (2022).
16. A. I. Lichtenstein, M. I. Katsnelson, V. P. Antropov, and V. A. Gubanov, *Local spin density functional approach to the theory of exchange interactions in ferromagnetic metals and alloys*, **J. Magn. Magn. Mater.** 67, 65 (1987) — the exchange-parameter mapping scheme extended relativistically within SPR-KKR.
17. L. Udvardi, L. Szunyogh, K. Palotás, and P. Weinberger, *First-principles relativistic study of spin waves in thin magnetic films*, **Phys. Rev. B** 68, 104436 (2003) — DMI mapping scheme, extended within the SPR-KKR relativistic KKR-GF framework.
18. Layer-KKR / photoelectron-diffraction extension: *Layered multiple scattering approach to Hard X-ray photoelectron diffraction: theory and application*, **npj Computational Materials** (2025), implementing PED/XPD within SPRKKR via a k-space layer-KKR formulation.
19. Relativistic Bogoliubov–de Gennes extension: *Relativistic spin-polarized KKR theory for superconducting heterostructures: Oscillating order parameter in the Au layer of Nb/Au/Fe trilayers* (and related BdG-KKR papers building on scalar- and fully relativistic implementations).
20. R. Eddhib et al., *ASE2SPRKKR: a unified Python framework integrating the Spin-Polarized Relativistic Korringa–Kohn–Rostoker method into the Atomic Simulation Environment*, arXiv:2608.05957 (2026) — the modern Python/ASE interface paper.

### 8.4 Complementary Textbooks / Broader Method References

21. J. Zabloudil, R. Hammerling, L. Szunyogh, and P. Weinberger, *Electron Scattering in Solid Matter: A Theoretical and Computational Treatise*, Springer Series in Solid-State Sciences 147, Springer (2005) — a standard reference text on KKR multiple-scattering theory underlying codes such as SPR-KKR.
22. J. S. Faulkner, G. M. Stocks, and Y. Wang, *Multiple Scattering Theory: Electronic Structure of Solids*, IOP Publishing, Bristol (2018).
23. H. Ebert, *Fully relativistic band structure calculations for magnetic solids — formalism and application*, in **Electronic Structure and Physical Properties of Solids: The Uses of the LMTO Method**, Lecture Notes in Physics 535, Springer (2000) — an early, detailed exposition of the fully relativistic spin-polarized KKR (i.e., SPR-KKR) formalism by the lead developer.

### 8.5 Software / Version Citations

24. H. Ebert et al., *The Munich SPR-KKR package*, versions 6.3–8.x, Ludwig-Maximilians-Universität München. (Cited variably as, e.g., "H. Ebert, A spin polarized relativistic Korringa–Kohn–Rostoker (SPR-KKR) code for calculating solid state properties," with version numbers 6.3, 7.7, 8.5, etc., available via the group's software distribution page.)

---

## 9. Suggested Citation Practice

Users of SPR-KKR conventionally cite, at minimum:
1. The **Rep. Prog. Phys. 74, 096501 (2011)** method paper (Ebert, Ködderitzsch, Minár), and
2. The **specific SPR-KKR software version** used (e.g., "The Munich SPR-KKR package, version X.Y"),

supplemented by the specific extension papers (CPA, LDA+U, DMFT, exchange-parameter/DMI theory, photoemission one-step model, BdG, ASE2SPRKKR, etc.) relevant to the property being calculated.

---

*Note: item numbering under Section 8 groups publications by theoretical role rather than strict chronology; several entries (especially in 8.1 and 8.3) represent the seminal papers most consistently referenced across SPR-KKR-based studies rather than an exhaustive bibliography of the entire KKR-GF literature.*

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the SPR-KKR 	Fully relativistic spin-polarized KKR Green's-function DFT code for solids, surfaces, and disordered systems. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
