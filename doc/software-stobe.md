# StoBe: A Gaussian-Basis DFT Code for X-ray Absorption Spectroscopy — Exhaustive Review

## 1. Overview and Identity

**StoBe** (also written **StoBe-deMon**, name derived from **Sto**ckholm–**Be**rlin) is a Linear-Combination-of-Gaussian-Type-Orbitals (LCGTO) Kohn–Sham Density Functional Theory (DFT) program for molecules and finite atomic/molecular clusters. It is historically one of the most widely used codes for simulating **near-edge X-ray absorption fine structure (NEXAFS/XANES)**, **X-ray emission spectra (XES)**, **X-ray photoelectron spectra (XPS)**, and related core-level spectroscopies via all-electron Gaussian-basis DFT.

- **Full name / branding:** StoBe-deMon
- **Lineage:** A specialized fork of the **deMon** ("**de**nsité **Mon**tréal") DFT program family, originally developed under the direction of **Dennis R. Salahub** at the Université de Montréal.
- **Distinguishing focus:** Core-level electronic structure and core-excitation spectroscopy — this is what separates StoBe from sibling deMon branches (deMon-KS, deMon2k).
- **License model:** Closed source; distributed under an academic/commercial license from the Fritz-Haber-Institut (Berlin) and Stockholm University groups. Source code available only with a purchased license.
- **Official/reference hosts:**
  - Fritz Haber Institute of the Max Planck Society (Berlin) — K. Hermann's group
  - Stockholm University — L. G. M. Pettersson's group (FYSIKUM)

## 2. Historical Development

### 2.1 Origins in deMon
- The deMon program was created for DFT calculations on atoms, molecules, and clusters, with a first widely available version around 1992, originating from the Salahub group in Montréal.
- Two largely independent lines of deMon development emerged in the 1990s:
  - The **Montpellier** branch (deMon-KS), oriented toward general molecular DFT and gradient-based methods.
  - The **Stockholm** branch, where core-level spectroscopy capability was implemented for the first time.
- In 1997 the Montpellier and Stockholm lines were merged into a unified **deMon-KS3** release.

### 2.2 Birth of StoBe
- Following the deMon-KS3 merger, the Stockholm branch was further and independently modified — in collaboration with **K. Hermann** at the **Fritz-Haber-Institut** in Berlin — specifically to strengthen and extend core-level spectrum calculations.
- This Berlin-modified Stockholm branch was subsequently distributed under its own name, **StoBe**, reflecting the Stockholm–Berlin collaboration, and became particularly well known and well suited for core-level (X-ray) spectroscopy work.
- Development and maintenance have since been jointly carried by the Hermann group (Berlin) and the Pettersson group (Stockholm), with numerous named contributing authors over the years (see citation block below).

### 2.3 Version history (as documented in the literature)
| Version | Approx. year cited in literature |
|---|---|
| StoBe-deMon (early distributed versions) | late 1990s–2000s |
| StoBe-deMon v3.0 | 2007 |
| StoBe-deMon v3.1 | 2011 |
| StoBe-deMon v3.3 | referenced in more recent (2020s) publications |

### 2.4 Relationship to other deMon-family codes
- **deMon-KS** — the general-purpose molecular DFT branch that eventually evolved separately.
- **deMon2k** — a modern, actively maintained, MPI-parallelized descendant of deMon-KS (led by A. M. Köster, P. Calaminici, M. E. Casida, D. R. Salahub, and the "International deMon Developers Community"), which also implements its own X-ray absorption/emission spectrum module, but is architecturally and developmentally distinct from StoBe.
- **ALLCHEM** — an independent, from-scratch DFT project initiated in Hannover (A. M. Köster group, 1995) exploring auxiliary-density methods; conceptually related in lineage but not a StoBe fork.
- StoBe and deMon2k should not be treated as interchangeable: StoBe remains the code most specifically engineered and most consistently cited for transition-potential core-excitation spectrum work, while deMon2k is the general-purpose, actively developed successor line.

## 3. Theoretical / Methodological Foundations

### 3.1 Electronic structure framework
- Self-consistent Kohn–Sham DFT.
- **LCGTO-DFT**: molecular orbitals expanded in Gaussian-type orbital (GTO) basis sets (all-electron, or with effective core potentials, ECPs, for heavier/non-excited centers).
- **Auxiliary (fitting) basis sets** are used to expand the electron density (and thus the Coulomb/exchange-correlation potential) in an approach descending from the Sambe–Felton / Dunlap variational Coulomb-fitting scheme central to the deMon lineage, avoiding four-center two-electron integral evaluation.
- Standard GGA/hybrid exchange–correlation functionals available (e.g., Becke '88 exchange, Perdew '86 correlation — "BP" or "PW86/PP" combinations — are the most frequently cited defaults in XAS work using StoBe).

### 3.2 Core-excitation methodology: the Transition-Potential (TP-DFT) approach
This is the theoretical core of StoBe's relevance to XAS:

- **Slater transition-state / transition-potential (TP) method**: A core-excited final state is approximated by performing a self-consistent DFT calculation with a **fractional occupation** (typically **half**, i.e., a "half core hole," HCH) removed from the relevant core orbital of the absorbing atom. Excitation energies are then obtained as orbital-energy differences (Kohn–Sham eigenvalue differences) from this transition-potential state, and oscillator strengths from the corresponding transition dipole matrix elements.
- This was systematically formalized and validated for X-ray absorption spectrum simulation by **Triguero, Pettersson, and Ågren (1998)**, who established the half-core-hole TP-DFT protocol as implemented in StoBe and benchmarked it against static-exchange, ΔSCF, and correlated wavefunction approaches.
- **ΔKohn–Sham (ΔKS) / ΔIP-TP hybrid corrections**: because TP eigenvalue differences describe *relative* transition energies well but can be less accurate for *absolute* binding/transition energies, a combined approach uses a full core-hole ΔSCF calculation for the absolute core-level binding energy and the TP calculation for relative (chemical-shift-resolved) transition energies — sometimes labeled ΔIP-TP in the later literature discussing this StoBe-style methodology.
- **Double/augmented (diffuse) basis technique**: to properly describe Rydberg and continuum-like final states above the ionization threshold, the basis set on the excited (core-hole) center is augmented with a large set of diffuse s, p, d Gaussian exponents (commonly cited as "19s19p19d" or "21 s,p,d diffuse exponents") after the ground/transition-potential density has converged; the Kohn–Sham matrix is rediagonalized in this extended basis to generate the additional virtual states needed for the spectral region far above threshold.
- **Stieltjes imaging**: because the augmented-basis diagonalization yields a discrete pseudospectrum rather than a true continuum, the portion of the computed spectrum above the ionization potential is reconstructed from these discrete states using the **Stieltjes imaging** procedure (moment-theory-based continuum reconstruction originally due to Langhoff and co-workers), converting discrete oscillator strengths into a smooth photoabsorption/photoionization cross-section profile.
- **Basis sets for the excited (core-hole) atom**: commonly the **IGLO-III** basis (originally developed for NMR shielding calculations, chosen here for its flexible, well-balanced core/valence description important for core-hole relaxation) fully uncontracted, sometimes supplemented by extra polarization/diffuse functions; non-excited atoms typically use smaller valence-only bases combined with effective core potentials (ECPs).
- **Effective core potentials (ECPs)**: used on spectator/non-excited atoms (and sometimes on heavy absorbing centers, e.g., transition metals) both to reduce cost and to aid SCF convergence of the core-hole state; model potentials such as those of Andzelm et al. are commonly cited for metal centers (e.g., Ru).

### 3.3 Related capabilities built on the same LCGTO-DFT/TP framework
- **NEXAFS/XANES** (near-edge X-ray absorption fine structure)
- **XES** (X-ray emission spectroscopy) — computed from ground-state Kohn–Sham eigenstates/transition moments rather than the core-hole TP state
- **XPS** (core-level binding energies / chemical shifts) via ΔKS and the "unrestricted generalized transition state" (UGTS) method
- **X-ray Raman scattering (XRS)** spectra (TP-based)
- **Auger decay** transition rate calculations (built on core-hole excited-state dynamics)
- Standard ground-state DFT properties: geometry optimization, vibrational analysis, population analyses, etc., inherited from the deMon lineage

## 4. Typical Computational Workflow in StoBe XAS Studies

A representative, frequently reproduced recipe drawn from the literature surveyed:

1. Build a finite cluster model of the system of interest (molecule, chemisorbed species + substrate slab/cluster, solvation shell, etc.).
2. Optimize geometry (or import experimental/MD-derived structure) using a standard GGA functional (commonly BP/PW86).
3. Select the absorbing atom(s); assign it a large, flexible, often fully uncontracted core-quality basis (e.g., IGLO-III).
4. Assign smaller valence-electron bases + ECPs to spectator atoms to control cost.
5. Perform a **half-core-hole transition-potential SCF** calculation on the absorbing atom's core orbital.
6. Extract discrete transition energies (orbital-energy differences) and oscillator strengths (dipole transition moments) for bound/valence-like final states.
7. Augment the absorbing atom's basis with diffuse s/p/d functions; rediagonalize to obtain a pseudospectrum of higher-lying (Rydberg/continuum-like) states.
8. Apply **Stieltjes imaging** to reconstruct the smooth continuum absorption profile above the ionization potential.
9. Optionally correct absolute transition energies using a **ΔSCF** (full core hole) calculation of the core ionization potential, combined multiplicatively/additively with the TP relative energies (ΔIP-TP).
10. Convolute discrete-state stick spectra with a Gaussian/Lorentzian broadening function (lifetime + experimental resolution) to produce the final simulated spectrum for comparison with experiment.

## 5. Domains of Application (illustrative, not exhaustive)

- Gas-phase small molecules (benchmark systems: H₂O, CO, N₂, NH₃, formaldehyde, etc.)
- Liquid water and ice — cluster-based modeling of hydrogen-bonding structure via O K-edge XAS (notably work by Pettersson, Nilsson, Leetmaa, Ljungberg and coworkers)
- Chemisorption/surface science — molecules on metal surfaces (e.g., CO on Ru(0001), N on Cu(100)) using cluster models of the surface
- Organic/bio-molecular systems — peptide bonds, amino acids, proteins (site-specific core-level electronic structure)
- Solid-state/molecular-crystal systems modeled via finite cluster extraction (e.g., C₆₀)
- Metal–organic interfaces and thin-film heterostructures (e.g., organic semiconductors on ferromagnetic metal substrates)
- Strained carbon-cage/diamondoid molecules
- Comparative benchmarking studies against other core-spectroscopy codes/methods (FEFF multiple-scattering, CPMD plane-wave core-hole pseudopotential method, TDDFT-based approaches, OCEAN/BSE, etc.)

## 6. Strengths and Recognized Limitations

**Strengths**
- Mature, extensively validated half-core-hole TP-DFT protocol specifically tuned for K-edge (and other) XAS/XES/XPS simulation.
- All-electron (or ECP-assisted) Gaussian-basis treatment gives good near-edge accuracy and access to XES/XPS/Auger properties within one consistent framework.
- Long track record (25+ years) with a large, cross-validated body of published spectra across gas, liquid, surface, and solid-cluster systems.
- Explicit machinery (diffuse-basis augmentation + Stieltjes imaging) for extending discrete bound-state calculations into the ionization continuum — a capability not universal among core-spectroscopy codes.

**Limitations**
- Closed-source/licensed distribution limits independent code auditing and broad community-driven development compared to open-source successors (e.g., some ERKALE, PSIXAS, or DFT/TDDFT-based open implementations that have since adopted or extended the same TP philosophy).
- Cluster-based (finite, non-periodic) representation of extended/periodic systems (surfaces, bulk solids) introduces finite-size/embedding approximations; comparisons with periodic plane-wave or multiple-scattering methods (CPMD, FEFF, OCEAN) are a recurring theme in the literature precisely because of this limitation.
- The (Δ)TP method is a mean-field, single-particle-excitation-energy approximation; it lacks explicit treatment of multi-electron correlation/relaxation effects captured by higher-level wavefunction methods (e.g., coupled-cluster core-valence-separated linear response), though it is computationally far cheaper and scales to much larger systems.
- Absolute transition energies from the pure half-core-hole TP eigenvalue differences can require empirical/ΔSCF correction for best agreement with experiment.
- Development pace/community size is smaller than for actively maintained modern periodic or open-source codes; deMon2k remains the more actively developed general-purpose successor within the same family.

## 7. Notable People / Institutions Associated with StoBe

| Name | Association |
|---|---|
| Dennis R. Salahub | Original deMon program (Université de Montréal); principal author on StoBe distribution |
| Lars G. M. Pettersson | Stockholm University (FYSIKUM); co-lead developer, transition-potential XAS methodology |
| Klaus Hermann | Fritz-Haber-Institut, Berlin; co-lead developer, principal maintainer of StoBe distribution/documentation |
| Hans Ågren | KTH Royal Institute of Technology, Stockholm; co-developer of TP-DFT XAS methodology |
| Luis Triguero | Lead author of the foundational TP-DFT XAS validation study |
| Mathias Leetmaa, Mattias P. Ljungberg, Anders Nilsson | Stockholm group; application/extension to condensed-phase (water/ice) XAS |
| M. E. Casida, C. Daul, A. Goursot, A. Koester, E. Proynov, A. St-Amant | Additional named authors on official StoBe-deMon version citations |
| Contributing authors (v3.0/3.1 releases) | V. Carravetta, H. Duarte, C. Friedrich, N. Godbout, J. Guan, C. Jamorski, M. Leboeuf, M. Nyberg, S. Patchkovskii, L. Pedocchi, F. Sim, A. Vela |

## 8. Access

- Official information/contact page: Fritz Haber Institute of the Max Planck Society — StoBe page (Klaus Hermann, Department of Inorganic Chemistry).
- Historically also referenced at fhi-berlin.mpg.de/KHsoftware/StoBe/ and via a "StoBe Software" license entity (Stockholm, Sweden).
- Distributed with the BALSAC visualization/structure package documentation bundle in some historical releases.
- Full source requires a purchased license; documentation (manual, basis-set database) has historically been available more openly.

---

## 9. Key Publications on StoBe's Theory and Methodology

### 9.1 Core code citation / official manual references
- K. Hermann, L. G. M. Pettersson, M. E. Casida, C. Daul, A. Goursot, A. Koester, E. Proynov, A. St-Amant, D. R. Salahub (contributing authors: V. Carravetta, H. Duarte, C. Friedrich, N. Godbout, J. Guan, C. Jamorski, M. Leboeuf, M. Leetmaa, M. Nyberg, S. Patchkovskii, L. Pedocchi, F. Sim, L. Triguero, A. Vela), **"StoBe-deMon"**, Version 3.0, StoBe Software, Stockholm, Sweden (2007).
- K. Hermann, L. G. M. Pettersson, M. E. Casida, C. Daul, A. Goursot, A. Koester, E. Proynov, A. St-Amant, D. R. Salahub, et al., **"StoBe-deMon"**, Version 3.1 (2011).
- StoBe-deMon Version 3.3 (as cited in more recent, 2020s literature), StoBe Software.

### 9.2 Foundational transition-potential DFT (TP-DFT) methodology for XAS
- L. Triguero, L. G. M. Pettersson, H. Ågren, **"Calculations of near-edge x-ray-absorption spectra of gas-phase and chemisorbed molecules by means of density-functional and transition-potential theory"**, *Phys. Rev. B* **58**, 8097 (1998). — The principal foundational paper establishing the half-core-hole TP-DFT protocol used in StoBe.
- L. Triguero, L. G. M. Pettersson, H. Ågren, *J. Phys. Chem. A* **102**, 10599 (1998). — Companion/extension study on the TP-DFT method for core-level spectra.
- L. Triguero, O. Plashkevych, L. G. M. Pettersson, H. Ågren, **"Separate state vs. transition state Kohn–Sham calculations for XPS core level binding energies and chemical shifts"**, *J. Electron Spectrosc. Relat. Phenom.* (1999). — ΔKS vs. UGTS benchmarking for XPS binding energies within the same framework.

### 9.3 Core-hole potentials, screening, and related theoretical review
- M. Leetmaa, M. P. Ljungberg, A. Lyubartsev, A. Nilsson, L. G. M. Pettersson, **"Theoretical approximations to X-ray absorption spectroscopy of liquid water and ice"**, *J. Electron Spectrosc. Relat. Phenom.* (2010). — Critical review of the TP-DFT approach as applied via StoBe-deMon to condensed-phase water/ice XAS.
- E. L. Shirley, L. G. M. Pettersson, D. Prendergast, **"Core-hole potentials and related effects"**, *International Tables for Crystallography*, Vol. I (2021). — Review contextualizing the (half) core-hole transition-potential approximation (including StoBe-based results) against other constrained-occupancy/self-consistent-field screening methods.

### 9.4 Underlying deMon program lineage and auxiliary-density DFT methodology
- Overview history: **"The Game of the Name"**, deMon Software User Guide/documentation (demon-software.com) — historical account of the deMon-KS1 → deMon-KS3 → StoBe branching.
- Foundational deMon DFT/auxiliary-density fitting methodology descends from the Sambe–Felton/Dunlap variational Coulomb-fitting approach used throughout the deMon program family (as implemented and documented across deMon-KS/StoBe/deMon2k releases).
- A. M. Köster, P. Calaminici, M. E. Casida, R. Flores-Moreno, G. Geudtner, A. Goursot, T. Heine, A. Ipatov, F. Janetzko, J. M. del Campo, S. Patchkovskii, J. U. Reveles, D. R. Salahub, A. Vela, **"deMon2k"**, The International deMon Developers Community, Cinvestav-IPN, Mexico (2006). — Sibling/successor code citation, useful for tracing shared theoretical ancestry.

### 9.5 Basis sets and effective core potentials used in the StoBe XAS protocol
- N. Godbout, D. R. Salahub, J. Andzelm, E. Wimmer, **"Optimization of Gaussian-type basis sets for local spin density functional calculations. Part I."**, *Can. J. Chem.* **70**, 560 (1992). — Source of the (DZVP/effective-core-potential-compatible) valence basis sets used for spectator atoms.
- W. Kutzelnigg, U. Fleischer, M. Schindler, **"The IGLO-Method: Ab initio Calculation and Interpretation of NMR Chemical Shifts and Magnetic Susceptibilities"**, in *NMR Basic Principles and Progress*, Vol. 23, Springer-Verlag, Berlin (1991), pp. 165–262. — Source of the IGLO-III basis set family adopted for the core-excited absorbing atom in StoBe XAS calculations.
- L. G. M. Pettersson, U. Wahlgren, O. Gropen, *J. Chem. Phys.* **86**, 2176 (1987). — Related basis-set/ECP methodology referenced alongside IGLO-III usage.
- K. Hermann, **"Orbital basis set, auxiliary functions and pseudopotentials for StoBe (Stockholm-Berlin version of deMon, a DFT molecule/cluster package)"**, online database, Fritz Haber Institute of the Max Planck Society, Berlin (Dec. 2000). — The canonical basis-set/ECP/auxiliary-function reference database for StoBe calculations.

### 9.6 Stieltjes imaging (continuum reconstruction) theory
- P. W. Langhoff, C. T. Corcoran, **"Stieltjes imaging of photoabsorption and dispersion profiles"**, *J. Chem. Phys.* **61**, 146 (1974). — Foundational Stieltjes-imaging/moment-theory paper underlying the above-threshold continuum reconstruction step used in StoBe XAS spectra.
- Additional Langhoff and co-workers papers on Stieltjes–Chebyshev moment theory and molecular photoionization cross sections (various dates, 1970s–1980s) — general theoretical basis for discretized-continuum reconstruction techniques as adapted within the LCGTO-DFT/StoBe framework.

### 9.7 Exchange-correlation functionals commonly cited in StoBe XAS work
- A. D. Becke, **"Density-functional exchange-energy approximation with correct asymptotic behavior"**, *Phys. Rev. A* **38**, 3098 (1988).
- J. P. Perdew, **"Density-functional approximation for the correlation energy of the inhomogeneous electron gas"**, *Phys. Rev. B* **33**, 8822 (1986); Erratum, *Phys. Rev. B* **34**, 7406 (1986).

### 9.8 Application/validation studies frequently cited as StoBe methodology touchstones
- L. Triguero, O. Plashkevych, L. G. M. Pettersson, H. Ågren — series of papers benchmarking Kohn–Sham DFT predictions of core photoelectron binding energies and chemical shifts (UGTS vs. ΔKS) against wavefunction-based references.
- Studies applying and cross-validating the TP-DFT/StoBe method against FEFF (multiple scattering), CPMD (plane-wave, core-hole pseudopotential), and TDDFT-based approaches for both molecular and condensed-phase (water/ice, metal-adsorbate) XAS, forming the broader empirical validation literature for the method (e.g., comparative XAS studies of CO/Ru(0001), water/ice O K-edge, and C₆₀ core excitation).

---

*Note: Several of the above entries (particularly §9.6–9.8) represent methodological building blocks broadly cited alongside StoBe applications rather than works published by the StoBe developer group itself; they are included because they constitute the theoretical backbone (Stieltjes imaging, exchange-correlation functionals, basis sets) that the code's core-excitation methodology depends on.*

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the StoBe 	Gaussian-basis DFT codes historically used for X-ray absorption spectroscopy simulations. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
