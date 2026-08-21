# RSPt (Relativistic Spin-Polarized toolkit): A Comprehensive Review

## 1. Overview and Identity

**RSPt** — the acronym stands for **Relativistic Spin-Polarized toolkit** — is an all-electron, full-potential electronic structure code based on the **Full-Potential Linear Muffin-Tin Orbital (FP-LMTO)** method. It is described by its documentation as a code for electronic structure calculations based on the Full-Potential Linear Muffin-Tin Orbital method. The code traces its lineage to work by **John M. Wills**, whose original FP-LMTO implementation dates to the mid-1980s, and it has since been developed collaboratively, primarily by groups at **Uppsala University** (Sweden), **Los Alamos National Laboratory** (USA), and collaborators in France. RSPt offers a robust and flexible set of tools to calculate total energies, magnetic moments, band structures, Fermi surfaces and densities of states for all elements, and combinations thereof, over a wide range of volumes and structures.

The canonical reference for the method is the 2010 Springer book:

> J.M. Wills, M. Alouani, P. Andersson, A. Delin, O. Eriksson, O. Grechnyev, *Full-Potential Electronic Structure Method: Energy and Force Calculations with Density Functional and Dynamical Mean Field Theory*, Springer Series in Solid-State Sciences, Vol. 167 (Springer, Berlin, 2010).

The book notes that in 1986, when the code that was to become RSPt was developed enough to be useful, it was one of the first full-potential, all-electron, relativistic implementations of DFT. The code was initially developed as a personal research tool rather than a collaborative effort or product, and as it was meant to be minimally restrictive/flexible, it required experience to use effectively — attributes that inhibited but did not prevent its spread as a research tool. The official project website is `fplmto-rspt.org`.

---

## 2. Theoretical / Methodological Foundation

### 2.1 The Full-Potential LMTO Method

RSPt implements the **FP-LMTO** approach, which is distinguished from earlier "muffin-tin" and atomic-sphere-approximation (ASA) LMTO methods by making **no shape approximation** to the potential or charge density:

- RSPt is based on the Full-Potential Linear Muffin-Tin Orbital method, which allows for very small basis sets and fast calculations, without any restriction on the symmetry of the potential.
- Due to the full-potential character, the code does not have limitations dictated by the geometry of the problem under consideration.
- No restrictions are imposed on charge densities or potentials of the studied systems, which is especially important for anisotropic and layered structures.

Space is partitioned, as in all muffin-tin-derived methods, into non-overlapping **muffin-tin (MT) spheres** around each atomic site plus an **interstitial region**, but the potential and density inside and outside the spheres are expanded without spherical/shape truncation:

- The full-potential methodology involves expanding the non-spherical electron density and potential both in the muffin-tin regions and in the interstitial region, with a general definition of the mathematical functions used and their symmetry properties.
- The method offers minimal, double, and triple basis sets, whose convergence properties are documented numerically.

### 2.2 Multiple Energy Sets / "Dual" Basis

A defining methodological feature is the ability to treat **semi-core and valence states with the same angular quantum number but different principal quantum number** as separate basis "energy sets," an approach sometimes called a **dual basis**:

- RSPt accommodates multiple energy sets (i.e., valence and semi-core states) with the same angular quantum numbers but different principal quantum number, which is ideal to provide an accurate description of semi-core states.
- In both FP-LMTO and comparable FP-LAPW implementations, a dual basis set is employed to incorporate two valence wave functions with the same orbital quantum number (e.g., 3p and 4p functions for Fe).

In practical calculations, basis functions are constructed with several **kinetic-energy tails** per orbital character (commonly two or three), an approach illustrated repeatedly in the literature, e.g. three kinetic energy tails used for 4s and 4p states, with values -0.3, -2.8, and -1.6 Ry, or three kinetic energy tails for valence sp states at -0.1, -2.3, and -1.5 Ry, with only the first two tails used for other basis functions.

### 2.3 Relativistic Treatment

Consistent with its name, RSPt supports a hierarchy of relativistic treatments:

- **Scalar-relativistic** calculations (mass-velocity and Darwin terms, no spin–orbit), used in many production studies, e.g. scalar relativistic electronic structure calculations performed with FP-LMTO RSPt.
- **Scalar-relativistic + spin–orbit coupling (SR+SO)** as a perturbative correction.
- **Fully relativistic (four-component, Dirac-based)** treatment is available through closely related theoretical work; the underlying relativistic spin-polarized LMTO formalism (Ebert-type radial Dirac equation solutions expanded in the LMTO basis) has been described in detail in the literature associated with the RSPt lineage, where the Bloch state is expanded in terms of single-site solutions of the radial Dirac equation and their energy derivatives for the relativistic spin-polarized case.
- A dedicated fully-relativistic, four-component **Dirac-Kohn-Sham** derivative code, **dirac-fp**, was built directly on the FP-LMTO approach used in RSPt to benchmark scalar-relativistic, SR+SO, and full-Dirac treatments against one another, using the full potential linear muffin-tin orbital approach to electronic structure, comparing scalar relativistic, scalar relativistic plus spin-orbit, and Dirac treatments on FCC materials.

### 2.4 Exchange-Correlation Functionals

RSPt supports the standard hierarchy of semi-local functionals:

- Spin-polarized and/or spin-orbit calculations with several LDA and GGA functionals implemented.
- Commonly used forms in the literature include **LSDA (Perdew–Wang / PW92)** and **GGA-PBE** — e.g. LSDA and GGA treated respectively via the Perdew-Wang and Perdew-Burke-Ernzerhof parametrizations.

---

## 3. Beyond-DFT Capabilities

A major distinguishing strength of RSPt relative to many other DFT packages is its deep integration of **beyond-DFT electronic-structure methods**, developed largely by the Uppsala theory group.

### 3.1 DFT+U / LSDA+U

- Up-to-date implementations of beyond-DFT methods such as SIC, DFT+U, or DFT+DMFT are available.
- The +U correction is implemented in its most general, **rotationally invariant, full 4-index Slater-parametrized form**: the interaction vertex is taken as a full spin- and orbital-rotationally-invariant 4-index U-matrix, parametrized in terms of Slater parameters F⁰, F², and F⁴ for 3d states.
- Double-counting corrections such as the **fully localized limit (FLL)** are supported, e.g. the double-counting correction based on the fully localized limit, using a 4-index U-matrix constructed to be equivalent to formulations used in VASP in the absence of spin-orbit coupling.
- Orbital polarization correction (OPC) is also implemented as an approximate, computationally cheap analogue: RSPt includes an implementation of the orbital polarization correction, whose main idea is to include an approximate description of the second Hund's rule into the DFT problem; from a multipolar decomposition, the OPC term can be shown to be contained within LDA+U.

### 3.2 Self-Interaction Correction (SIC)

SIC is listed among the beyond-DFT methods implemented in RSPt. This has historically been applied to strongly localized f-electron systems (e.g., rare-earth-based magnetic scattering studies referenced in the literature on self-interaction-corrected relativistic magnetic scattering theory).

### 3.3 DFT+DMFT (Dynamical Mean-Field Theory)

This is arguably RSPt's flagship advanced capability and the reason the code is a "tool of choice" for the Uppsala materials-theory / DMFT group:

- RSPt includes a fully general DFT+DMFT implementation, with self-consistent cycles over self-energy and charge density.
- **Charge self-consistency (CSC-DMFT)** was formally implemented and documented in Grånäs, Di Marco, Thunström, Nordström, Eriksson, Björkman & Wills, *Comput. Mater. Sci.* **55**, 295 (2012): full charge self-consistence over the electron density was implemented into the LDA+DMFT scheme based on FP-LMTO, with computational details on the construction of the electron density from the density matrix provided; the method was tested on the charge-transfer insulator NiO with a static Hartree-Fock impurity solver, on bcc Fe with the spin-polarized T-matrix fluctuation-exchange (SPTF) solver, and on the permanent magnet SmCo₅ using both SPTF and Hubbard-I solvers because the Coulomb-interaction strength differs drastically between the Sm and Co sites.
- Multiple **impurity solvers** are available within the RSPt DMFT framework, including:
  - **Spin-dependent T-matrix Fluctuation Exchange (SPTF)** — developed by the Uppsala group as part of RSPt's DMFT solver suite.
  - **Hubbard-I approximation (HIA)** — used for strongly localized f-shells, shown to give excellent results for strongly correlated f-orbitals and chosen as the impurity solver for localized 4f orbitals such as Sm in SmCo₅.
  - **Exact diagonalization (ED)** solvers, including with hybridization/bath states — RSPt was used to perform LDA+DMFT calculations both with and without hybridization effects, using an exact-diagonalization solver, with bath states constructed from the muffin-tin "head" projection of the LMTOs.
  - A more recent **configuration-interaction (CI) solver** has also been built on top of RSPt's DMFT infrastructure for efficient treatment of correlated f-electron problems on modest hardware.
- **Local-orbital projection schemes**: RSPt provides at least two distinct schemes for constructing the correlated local orbitals used in DMFT — the **Löwdin-orthonormalized LMTO ("ORT")** basis and the **muffin-tin-head-projected ("MTH")** orbitals: the non-orthonormal LMTO basis set is used in DFT, the Löwdin-orthonormalized LMTO ("ORT") is used as a first projection scheme, and muffin-tin head-projected orbitals ("MTH") are used as a second projection scheme.

### 3.4 Interatomic Exchange Interactions (Jij) and Spin Dynamics Interfacing

RSPt has strong, well-developed post-processing capability for extracting **Heisenberg exchange parameters** $J_{ij}$, both at the DFT and DFT+DMFT levels, using the **magnetic force theorem / infinitesimal rotation method**:

- Once the electronic structure is converged, magnetic excitations are mapped onto a Heisenberg Hamiltonian, with pairwise exchange interactions extracted via the method of infinitesimal rotation of the spins, using the local magnetic force approach.
- At the correlated level: the crucial ingredient is a dynamical on-site exchange potential built from the difference of spin-dependent Hamiltonians and self-energies, traced over Matsubara frequencies and angular-momentum states together with the inter-site Green's functions.
- These $J_{ij}$ values feed naturally into **atomistic spin-dynamics** codes (e.g. UppASD) for finite-temperature magnetism studies (Curie temperatures, spin-wave spectra, etc.).

---

## 4. Ground-State / Standard DFT Capabilities

Beyond its "beyond-DFT" specialization, RSPt is a full production DFT engine offering the standard suite of ground-state properties:

- Total energies, equilibrium volumes/lattice constants, bulk moduli, equations of state
- Spin and orbital magnetic moments (collinear and non-collinear-adjacent settings, though note the caveat in §6 below)
- Band structures and Fermi surfaces
- Densities of states (total, atom-, orbital-, and spin-resolved)
- Forces (for structural relaxation), as reflected in the book's subtitle *"Energy and Force Calculations..."*
- **Open-core / spin-polarized core treatment** for strongly localized shells (notably rare-earth 4f electrons) as an alternative to explicit valence treatment: 4f electrons treated as spin-polarized core states (open-core approximation), a setup shown to avoid artificial hybridization effects seen in simple DFT treatments of hcp Gd.
- **Brillouin-zone integration** via the **tetrahedron method with Blöchl corrections**, as commonly cited in application papers: integration over the Brillouin zone using the tetrahedron method with Blöchl's correction.

---

## 5. Software Characteristics

| Aspect | Details |
|---|---|
| **Language / architecture** | Fortran-based, MPI-parallel self-consistent-field code |
| **Distribution** | Described in some third-party listings as an "Open Source" project, though access has historically been gated through registration/collaboration with the developer group rather than an anonymous public repository — the book's chapter *"Obtaining RSPt from the Web"* documents this acquisition process |
| **License cost** | Listed as free of charge by HPC-center documentation. |
| **HPC deployment** | Distributed via module systems on national HPC infrastructure — e.g. on PDC's Dardel supercomputer, RSPt is built with EasyBuild and run via MPI batch submission (`srun -n <ntasks> rspt`), with example input files shipped in the code's testsuite directory. |
| **Documentation** | Official book (Wills et al. 2010), an FAQ hosted at Uppsala University's physics department, and community-produced video tutorials |
| **Community tutorials** | A publicly available video-tutorial series and companion GitHub repository walk through crystal-structure setup, basis/energy-set configuration, self-consistency, adding semi-core states, spin-polarized core calculations, LDA+U vs LSDA+U setup, band-structure plotting, and Jij calculations, using example systems such as Co₂FeAl, NiO, and hcp Gd/Tb. |

---

## 6. Notable Limitations / Caveats Reported in the Literature

- **Fixed-axis collinearity in some contexts**: at least one cited application explicitly notes the code used is a collinear spin code with a fixed spin quantization axis, requiring different spin-axis orientations to be considered as separate calculations — a constraint relevant to magnetocrystalline-anisotropy-type studies, though this reflects the setup used in that particular study rather than necessarily an absolute code-wide restriction across all RSPt modes.
- **Steep learning curve**: as the developers themselves acknowledge, because the code required knowledge of its inner workings and was built to be maximally flexible rather than maximally convenient, effective use required real user experience — a recurring theme distinguishing RSPt from more "black-box" DFT codes (e.g., VASP, Quantum ESPRESSO).
- **Access model**: unlike fully open, anonymously downloadable codes, obtaining and building RSPt has historically involved more direct engagement with the developer group (registration, EasyBuild recipes on specific HPC systems, etc.) rather than a simple public package-manager install.
- **Niche computational core**: as an all-electron full-potential method with explicit muffin-tin/interstitial geometric partitioning, RSPt is computationally heavier per atom than pseudopotential plane-wave codes, though the LMTO minimal-basis philosophy is specifically designed to offset this (small basis sets, fast individual calculations).

---

## 7. Position in the DFT Code Landscape

RSPt sits within the **all-electron, full-potential** family of DFT codes, alongside methods such as **FP-LAPW** (e.g., WIEN2k, Elk) and other full-potential LMTO/screened-KKR-type efforts (e.g., **Questaal**, **EMTO**). Comparative usage in the literature is common:

- Studies have used two independent all-electron approaches side-by-side — FP-LAPW as implemented in WIEN2k, and FP-LMTO as implemented in RSPt — to cross-validate DMFT-derived exchange parameters.
- RSPt (FP-LMTO) and Elk (FP-LAPW) have been directly compared for magnetic-oxide systems such as FeTe, BiFeO₃, and hexaferrites, exploiting the shared full-potential, dual-basis philosophy of both methods.
- RSPt has also served as an all-electron reference/benchmark for validating pseudopotential-based approaches, such as in-situ pseudopotential construction compared against Quantum ESPRESSO.

RSPt's particular niche — **small, physically transparent LMTO basis sets + genuinely full-potential treatment + tightly integrated, charge-self-consistent DFT+DMFT** — makes it especially attractive for the **strongly correlated electron materials community** (transition-metal oxides, rare-earth/actinide intermetallics, Heusler alloys, permanent-magnet materials, and magnetic thin films/heterostructures), rather than for high-throughput screening or very large supercell simulations, which are typically better served by linear-scaling or plane-wave/pseudopotential codes.

---

## 8. Representative Application Domains (from the literature)

- Magnetic Heusler alloys and half-metallic heterostructures (e.g., Co₂MnAl/CoMnVAl)
- Rare-earth magnetism (hcp Gd under pressure, Sm-Co permanent magnets)
- Transition-metal-oxide Mott/charge-transfer insulators (NiO)
- Surface and thin-film magnetism, layer-resolved exchange interactions (late 3d-metal surfaces)
- Actinide and lanthanide correlated-electron physics (Ce compounds, Pu)
- 2D-material-supported magnetic adatoms (e.g., transition-metal adatoms on monolayer NbSe₂)
- Multiferroics and hexaferrites (BiFeO₃, SrFe₁₂O₁₉)
- X-ray magnetic circular dichroism and magnetic scattering of X-rays

---

## 9. Key Primary References

1. J.M. Wills, M. Alouani, P. Andersson, A. Delin, O. Eriksson, O. Grechnyev, *Full-Potential Electronic Structure Method: Energy and Force Calculations with Density Functional and Dynamical Mean Field Theory*, Springer Series in Solid-State Sciences 167 (Springer, Berlin, 2010).
2. J.M. Wills, O. Eriksson, M. Alouani, D.L. Price, "Full-Potential LMTO Total Energy and Force Calculations," in *Electronic Structure and Physical Properties of Solids*, ed. H. Dreysse (Springer, Berlin, 2000).
3. O. Grånäs, I. Di Marco, P. Thunström, L. Nordström, O. Eriksson, T. Björkman, J.M. Wills, "Charge self-consistent dynamical mean-field theory based on the full-potential linear muffin-tin orbital method: Methodology and applications," *Comput. Mater. Sci.* **55**, 295 (2012).
4. Official project site: fplmto-rspt.org (Uppsala University materials-theory group)

---

*Note: This review draws on official project/documentation summaries and citation-format method descriptions from peer-reviewed papers that used RSPt, since the primary technical monograph and project site are not fully text-accessible for direct excerpting.*

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the RSPt (Relativistic Spin Polarized toolkit) 	Full-potential linear muffin-tin orbital (FP-LMTO) all-electron DFT code for solids. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
