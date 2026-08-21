# AkaiKKR (MACHIKANEYAMA): An Exhaustive Review

## 1. Overview

AkaiKKR — also distributed under the name **MACHIKANEYAMA** — is a first-principles electronic-structure package built around the **Korringa–Kohn–Rostoker (KKR) Green's-function method**, developed and maintained primarily by **Hisazumi Akai** (Osaka University / ISSP, University of Tokyo) as a program package used for first-principles calculation of electronic structures of metals, semiconductors and compounds, in the framework of the local density approximation or generalized gradient approximation (LDA/GGA) of density functional theory. The package has been continuously developed since the late 1970s and is still being developed by various authors; each program in the package is written in FORTRAN 77, is completely self-contained without external library dependencies, and runs on anything from a laptop to a supercomputer under UNIX, Linux, macOS or Windows.

The code's defining scientific niche is its combination of an **all-electron KKR Green's-function formalism** with the **coherent potential approximation (CPA)** for chemical/magnetic disorder — a pairing that makes it one of the most widely used tools worldwide for first-principles modeling of random substitutional alloys.

---

## 2. Theoretical and Methodological Foundations

### 2.1 The KKR Green's-function approach

Unlike wavefunction/eigenvalue-diagonalization methods (plane-wave pseudopotential codes, LAPW, etc.), KKR is a **multiple-scattering theory**. Space is partitioned into non-overlapping (or space-filling) cells around each atomic site; within each cell the Kohn–Sham/Schrödinger equation is solved as a scattering problem, producing single-site regular and irregular wavefunctions and a **t-matrix**. The site t-matrices are combined via the **scattering path operator (τ-matrix)**,

$$
\underline{\tau}^{nm}(\epsilon) = \underline{t}^{n}(\epsilon)\,\delta_{nm} + \underline{t}^{n}(\epsilon)\sum_{k\neq n}\underline{g}^{nk}(\epsilon)\,\underline{\tau}^{km}(\epsilon),
$$

which directly yields the **single-particle Green's function** of the full system This is done using multiple scattering theory: each atomic cell is treated as an electron scattering center, for which the local ("single-site") solutions and the scattering t-matrix can be calculated, and the t-matrices are used to construct the multiple scattering path matrix τ. Physical observables (charge density, DOS, magnetic moments) follow directly from energy integrals/contour integrations of the Green's function rather than from summing over occupied eigenstates, avoiding explicit diagonalization of a secular matrix for each k-point/band.

Because the method works with the Green's function rather than Bloch wavefunctions, it naturally generalizes to non-periodic and statistically disordered environments — impurities, semi-infinite surfaces, and substitutional alloys — without supercells.

**Key advantages inherited from the KKR-GF formalism:**
- **All-electron treatment**: no pseudopotential/PAW approximation; core and valence states solved self-consistently, giving direct access to core-level and electric-field-gradient information.
- **No plane-wave cutoff error**: it does not suffer from any serious truncation errors such as those of the plane-wave cutoff. Basis completeness is instead controlled by the angular-momentum cutoff ℓ_max.
- **Direct access to the Green's function** as a good starting point for first-principles linear response theory, many-body theory, and so on — relevant for transport (Kubo–Greenwood), spin-fluctuation theory, and exchange-coupling calculations.

### 2.2 Potential shape: ASA and muffin-tin construction

AkaiKKR's traditional and most widely used mode employs the **muffin-tin / atomic-sphere approximation (ASA)** for the crystal potential — spherically symmetric potentials inside touching or slightly overlapping atomic spheres, with a constant (muffin-tin zero) interstitial potential. This is the origin of much of its speed: single-site scattering problems reduce to radial Schrödinger/Dirac equations, and the whole multi-site problem reduces to inverting a KKR structural matrix of moderate size, indexed by $(\ell,m)$ per site with $L_{\max}=(\ell_{\max}+1)^2$. Typical production calculations use $\ell_{\max}=2$–$3$ (s,p,d,[f]).

The trade-off is that ASA/muffin-tin geometry is less accurate than full-potential treatments for open or strongly non-spherical structures (e.g., layered/molecular crystals, strongly covalent semiconductors, surfaces with large corrugation) — a limitation shared with early implementations of essentially all KKR codes and discussed further in §7.

### 2.3 Exchange–correlation functionals and relativistic treatment

AkaiKKR supports standard **LDA** parametrizations (von Barth–Hedin, Moruzzi–Williams–Janak, Vosko–Wilk–Nusair) and **GGA** functionals (Perdew–Wang, PBE), selectable at input time. Spin-polarized (collinear) calculations are standard; the potential/Hamiltonian is typically solved in the **scalar-relativistic approximation**, with some interfaces/derivative implementations extending to full spin-orbit/relativistic (Dirac) treatments for magnetocrystalline-anisotropy work (see §5).

### 2.4 The Fast KKR-CPA Algorithm — the central technical contribution

The single most consequential technical innovation associated with AkaiKKR is Akai's 1989 **fast KKR-CPA algorithm**: A fast KKR CPA method is explained and its convergence properties are examined numerically. It is shown that a step number of N ≃ 300, which determines the number of k-points used for the numerical integration in the k-space as well as the number of iteration steps in determining the coherent t-matrix, is sufficient for most purposes, including total-energy calculations.

Conceptually, CPA replaces a chemically/magnetically random lattice by an effective, translationally-invariant **coherent medium**, whose single-site t-matrix $\bar{t}$ is determined self-consistently by requiring that embedding any real alloy constituent into the effective medium produces, on average, zero additional scattering:

$$
c_A\, \tau^{A} + c_B\, \tau^{B} + \cdots = \bar{\tau},
$$

with concentrations $c_i$ and component-resolved scattering path operators $\tau^i$ built from the coherent medium via the usual Dyson-type embedding. This is the **Soven equation** generalized to the multi-site KKR structure [3] P. Soven, Phys. Rev., 156 (1967) 809, and it is solved self-consistently alongside the usual DFT charge-density loop. Akai's contribution was to make the k-space Brillouin-zone integration required for the coherent τ-matrix numerically efficient and to demonstrate that a modest number of iteration/integration steps suffices for converged total energies — this is what earns the "fast KKR-CPA" label and what underlies AkaiKKR's practical efficiency for high-throughput alloy screening. This algorithm and its extension to spin-polarized systems were validated on ferromagnetic Ni–Fe systems in the framework of the KKR CPA combined with the local spin density approximation.

Because CPA is applied at the **single-site (mean-field) level**, chemical/spin configurations enter only through site-averaged occupation probabilities $c_i$; the method is exact in the dilute limit and in the ordered limit, but by construction neglects **short-range order (SRO)** and **local lattice relaxation/environment effects** around a disordered site — a well-known limitation of single-site CPA generally (see §7), addressed elsewhere in the KKR literature by cluster extensions such as **KKR-NLCPA** The Korringa-Kohn-Rostoker non-local coherent potential approximation systematically provides a hierarchy of improvements upon the widely used KKR-CPA approach and includes non-local correlations in the disorder configurations by means of a self-consistently embedded cluster; it remains fully causal, becomes exact in the limit of large cluster sizes, and reduces to the KKR-CPA for a single-site cluster. — which is not itself part of the standard AkaiKKR distribution.

---

## 3. Software Architecture and Distribution

| Aspect | Detail |
|---|---|
| **Language** | FORTRAN 77 |
| **Dependencies** | Completely self-contained; no additional libraries required |
| **Portability** | UNIX, Linux, Mac OS, Windows — any platform with a Fortran compiler |
| **Distribution** | Source-code registration/download via the official ISSP-hosted site; historically required user account registration |
| **Development lineage** | Machikaneyama package → *Machikaneyama2000* → *Machikaneyama2002* → *AkaiKKR CPA20xx* release series (e.g. `akaikkr_cpa2021v01`) |
| **Primary author/maintainer** | Hisazumi Akai (Osaka University; ISSP, University of Tokyo) |
| **Hosting institution** | Institute for Solid State Physics (ISSP), University of Tokyo (`kkr.issp.u-tokyo.ac.jp`) |

The package is organized as a set of driver programs invoked through short input decks with a "first key" selecting the run mode — e.g. **`go`** for a self-consistent-field ground-state calculation, **`dos`** for density-of-states evaluation, and further modes for Bloch spectral functions and related post-processing you can get the local density of state as well as the total dos using "dos" instead of "go" as the first key of your input; all the dos data are output as numerical data that can be plotted with any plotting software, and the total dos can be plotted using the "gpd" routine provided with the package. The pattern (`go`, `dos`, `fsm`, `gofmg`, `spc31`, …) observed in the package's own regression-test suite indicates a modular family of executables sharing a common input format, each specialized to a class of physical output (ground state, DOS, fixed-spin-moment, Bloch spectral function, etc.).

A companion **`AkaiKKRPythonUtil`** package (Apache-2.0 licensed, `AkaiKKRteam` on GitHub) provides Python-side helpers for input generation (including CIF-to-input conversion) and rudimentary output parsing/plotting, and ships with a self-test suite exercising representative systems (elemental Fe, Ni, Cu; compounds such as Co₂MnSi, GaAs; disordered/magnetic test cases such as NiFe, FeRh₀.₅Pt₀.₅, SmCo₅, and a local-moment-disorder case `Fe_lmd`) Copyright (c) 2021-2023 AkaiKKRteam, distributed under the Apache License, Version 2.0; it tests the AkaiKKR package and generates AkaiKKR input files from cif files.

More recently, an independent, actively developed community package, **`akaitools`**, addresses a long-standing usability gap: AkaiKKR has been in use for decades, but structured programmatic access to its output remains largely manual. The official AkaiKKRPythonUtil package parses output through individual methods that return a mix of scalars, lists, tuples, dicts, and DataFrames with no unified data model; its higher-level plotting layer couples parsing and figure generation together; the package is not on PyPI and is incompatible with current Python environments and ships no documentation site. `akaitools` is a Python package that parses AkaiKKR output files into structured, type-annotated Python objects, covering three output types: SCF results (convergence history and per-atom electronic/magnetic properties), spin-resolved orbital-projected DOS for each CPA component, and Bloch spectral functions on a user-defined k-point path. Practically, it offers frozen dataclass models backed by NumPy arrays with eV-conversion helpers, DOS/SCF export to pandas via `.to_dataframe()`, Matplotlib-based plotting, a terminal CLI (`akaitools go|dos|spc`), and programmatic input-file generation via an `InputFile` class supporting CPA alloys, multi-site structures, and k-path definitions for spectral-function calculations. Its emergence is itself informative: it signals that AkaiKKR's native text-based I/O had, over decades of use, become a genuine bottleneck for the kind of high-throughput/data-driven workflows (materials screening, ML-model training sets) increasingly demanded of legacy scientific codes.

---

## 4. Core Capabilities

- **Self-consistent electronic-structure calculation** (`go`-type runs): charge density, total energy, magnetic moments, Fermi energy, for ordered crystals, impurities embedded in a host, and CPA-disordered alloys/mixed crystals.
- **Density of states**: total and site/orbital-projected DOS, spin-resolved, including CPA-component-resolved DOS for each alloy species on a given sublattice.
- **Bloch spectral functions (BSF)**: the CPA analogue of a band structure — the disorder-broadened $A(\mathbf{k},E)$ obtained from the imaginary part of the configurationally averaged Green's function along a chosen k-path, used to visualize how alloying broadens and shifts otherwise sharp bands.
- **Magnetic structure options**: ferromagnetic, antiferromagnetic, and **disordered local moment (DLM)** treatments (paramagnetic-state modeling via CPA over up/down spin orientations on otherwise equivalent sites) — used extensively for finite-temperature magnetism and local-moment-disorder studies (cf. the `Fe_lmd` test case above and the Akai–Dederichs local-moment-disorder formalism referenced in the magnetism literature).
- **Fixed-spin-moment (FSM) calculations**, useful for mapping total-energy-vs-moment curves in itinerant magnets.
- **Coherent potential approximation** for substitutional site disorder on any sublattice, with an arbitrary number of components and concentrations, allowing the effects of atomic disorder and substitution doping to be treated, with the key advantage that doping is generated automatically — without the need to create a supercell — since it is generated directly by CPA.
- **Total-energy-based structural/thermodynamic screening**: equilibrium lattice constants (via Birch–Murnaghan-type equation-of-state fits to total energy vs. volume), formation energies of random alloys and high-entropy alloys, elastic constants (via full-potential extensions), used heavily in high-entropy-alloy (HEA) screening pipelines.

---

## 5. Application Domains

AkaiKKR's use in the literature clusters around several recurring problem classes:

1. **Random substitutional alloys and solid solutions** — the archetypal use case, from binary transition-metal alloys (Ni–Fe, Cu–Ni, Fe–Cr) to complex multi-component **high-entropy alloys (HEAs)**, where CPA sidesteps the combinatorial supercell-enumeration problem entirely. A representative large-scale application combined a high-throughput phase-prediction algorithm for PtPd-based high-entropy alloys via combined KKR-CPA and artificial-neural-network (ANN) techniques, using KKR-CPA to generate 2,720 data points of formation energy and lattice parameters, followed by ANN training on 36,556 data points across 9,139 HEA systems with 137,085 features, verified by R² close to unity.
2. **Magnetism in itinerant and dilute-magnetic systems**: ferromagnetic transition-metal alloys, magnetocrystalline anisotropy of disordered alloys computed via the spin-polarized relativistic KKR-CPA a first-principles theory of magnetocrystalline anisotropy of disordered alloys within the framework of the spin-polarized fully relativistic Korringa-Kohn-Rostoker coherent-potential approximation, in which the method calculates the MAE directly rather than by subtracting two large total energies, since MAE (~μeV) is orders of magnitude smaller than the total energy; dilute magnetic semiconductors; half-metallic ferromagnets in doped semiconductor hosts (CdTe, SnC, HfO₂, MgSe) reported across numerous applied-physics studies.
3. **Rare-earth permanent-magnet design**, where the AkaiKKR/CPA methodology has been used alongside related tools for understanding and optimization of hard magnetic compounds from first principles.
4. **Finite-temperature magnetism / Curie-temperature estimation**: DLM-CPA treatments of the paramagnetic state combined with exchange-parameter extraction and Monte Carlo/mean-field estimates of $T_C$.
5. **Transport and residual-resistivity calculations** in disordered alloys, exploiting the Green's-function/Kubo formalism, e.g. impurity resistivity of fcc and hcp Fe-based alloys relevant to thermal stratification at the top of the core of super-Earths, where the Kohn-Sham equation was solved via the KKR Green's-function method implemented in AkaiKKR, using the atomic sphere approximation, lmax = 3, scalar-relativistic effects, and CPA for substitutional chemical disorder.
6. **Superconductivity screening in disordered/HEA superconductors**, where free energy and electronic-structure calculations were performed using the KKR method, implemented together with CPA for disordered systems, in which a random arrangement of all elements is replaced by an ordered lattice representing an average over all possible disordered-lattice configurations within the simple unit cell (bcc or fcc), using the AkaiKKR (MACHIKANEYAMA) package with a PBE exchange-correlation functional constructed in the semirelativistic muffin-tin approximation. A related Matthias-rule study of bcc HEA superconductors used the same lineage of methods and citations.
7. **Full-potential KKR-CPA extensions**, used for elastic-constant prediction pipelines: computational materials design based on first-principles electronic-structure calculations, demonstrated for construction of an elastic-constant prediction model via machine learning applied to a database of elastic constants of 2,555 BCC HEAs generated by the full-potential Korringa-Kohn-Rostoker coherent-potential-approximation (FPKKR-CPA) method.
8. **Mineral/oxide and thin-film semiconductor property prediction** (spray-pyrolysis thin films, doped oxides), typically as one leg of a combined experiment+DFT study.

---

## 6. Typical Computational Settings and Workflow

Representative published settings converge on a fairly narrow, well-tested parameter regime:

- **Potential shape**: muffin-tin or ASA, semirelativistic (scalar-relativistic).
- **Exchange-correlation**: PBE-GGA is now most common; VWN/Moruzzi–Williams–Janak LDA remains used, with several studies noting that the specific choice of exchange-correlation functional may not significantly affect derived properties such as impurity resistivity.
- **Angular-momentum cutoff**: $\ell_{\max}=2$–$3$ is standard; larger $\ell_{\max}$ is occasionally used for full-potential refinements.
- **Brillouin-zone sampling**: on the order of a few hundred irreducible k-points for SCF/DOS convergence (consistent with Akai's own N≈300 step-number guidance for the CPA k-integration).
- **Workflow**: crystal structure (often from a CIF) → input-deck generation (manually or via `AkaiKKRPythonUtil`/`akaitools`) → `go` SCF run to self-consistency → `dos`/spectral post-processing → downstream analysis (equation-of-state fitting, magnetic-moment/Curie-temperature extraction, elastic-constant fitting, ML-feature generation) via community Python tooling.

---

## 7. Strengths

- **Speed and compactness for CPA-alloy screening.** The fast KKR-CPA algorithm makes AkaiKKR genuinely practical for scanning large alloy composition spaces (hundreds to thousands of compositions), which is precisely why it underlies most of the HEA/high-throughput-screening literature above.
- **No supercell needed for disorder.** CPA handles arbitrary concentrations natively — including physically inaccessible fractional compositions used to trace trends — generating doping automatically without the need to create a supercell.
- **All-electron accuracy** without a plane-wave cutoff, giving direct, unbiased access to core-level and near-nucleus quantities.
- **Self-contained, dependency-free FORTRAN 77**, trivially portable and easy to build on essentially any platform with a Fortran compiler — a real practical advantage over codes with heavy external-library stacks.
- **Direct Green's-function output**, giving a natural entry point to transport, spectral-function, and linear-response extensions beyond the ground state.
- **Decades of validation** across a very wide range of physical problems (magnetism, transport, thermodynamics, superconductivity screening), with the foundational 1989 Akai paper remaining one of the most consistently cited references in the KKR-CPA literature.

## 8. Limitations

- **Single-site (mean-field) CPA neglects short-range order and local relaxation.** As with all standard CPA implementations, correlations between neighboring sites' occupations, local lattice distortions around unlike neighbors, and clustering/ordering tendencies are not captured; more sophisticated variants such as KKR-NLCPA exist in the broader literature to include non-local correlations in the disorder configurations by means of a self-consistently embedded cluster but are not part of the standard AkaiKKR distribution.
- **ASA/muffin-tin geometry** is less well suited to open structures, strongly directional bonding, and surfaces/interfaces with large potential non-sphericity, relative to genuinely full-potential all-electron methods (though full-potential KKR-CPA extensions are used in some derivative studies, e.g. elastic-constant work).
- **Legacy Fortran 77 codebase and text-based I/O.** The lack of a structured output format has been an explicit, recently documented pain point: its output is unstructured plain text with no programmatic interface, leaving data extraction entirely to the user and making systematic or high-throughput studies impractical — a limitation only now being addressed by third-party tooling (`akaitools`) rather than the core package itself.
- **Access friction.** Distribution has historically required registration on the ISSP-hosted site rather than open publication on a code-hosting platform, and a website security renewal invalidated old user accounts, requiring re-registration — a minor but real barrier to reproducibility and onboarding relative to fully open-source competitors.
- **Documentation is comparatively sparse in English** and largely tutorial-oriented (e.g., the widely cited Sato/Ogura/Akai guide), rather than a comprehensive reference manual; much community knowledge is distributed across a web-forum/bulletin-board archive rather than centralized documentation.
- **Officially supported Python tooling was, until recently, unmaintained and difficult to install**, per the `akaitools` authors' assessment noted above.

## 9. Comparison with Related KKR Codes

| Code | Distinguishing focus | Disorder treatment | Potential shape |
|---|---|---|---|
| **AkaiKKR (MACHIKANEYAMA)** | Speed/compactness, CPA alloys, magnetism | Fast single-site KKR-CPA | ASA/muffin-tin (primary); FP extensions in derivative work |
| **JuKKR** (Jülich) | Relativistic effects, scattering, impurities, large systems | Impurity embedding (KKRimp), nano/disordered solids (KKRnano) | Full-potential The Korringa-Kohn-Rostoker (KKR) Green's function method is a highly accurate all-electron method for DFT, and the most important features of the Jülich KKR codes include relativistic calculations, scattering effects, and finite-sized clusters or very large systems. |
| **SPR-KKR** (Munich) | Spin-polarized relativistic (Dirac) KKR, spectroscopy | KKR-CPA | ASA and full-potential modes |
| **MECCA** | High-pressure/high-temperature, warm dense matter | k-space CPA | ASA, scalar-relativistic MECCA is applicable to the whole pressure and temperature range of interest, beyond what is available from pseudopotential methods, though as presently implemented it is a static DFT code that does not sample ionic degrees of freedom explicitly. |
| **GreenALM** | Disordered-alloy-focused KKR-DFT | CPA / DLM | Green's-function-based charge density GreenALM implements a usual DFT Kohn-Sham iteration loop that eventually produces the total energy and charge density; the only difference from a more widely adopted Hamiltonian approach is that the charge density itself is calculated from the Green's function instead of a set of wavefunctions. |

Within this landscape, AkaiKKR occupies the niche of the **fast, lightweight, alloy-screening workhorse**, trading some of the geometric flexibility of full-potential codes (JuKKR, SPR-KKR full-potential mode) for computational efficiency and CPA maturity — which explains its dominant footprint in the HEA/high-throughput alloy-design literature specifically, versus the surface/spectroscopy/topological-materials work more associated with JuKKR and SPR-KKR.

---

## 10. Summary Assessment

AkaiKKR (MACHIKANEYAMA) is a mature, scientifically well-validated, and computationally efficient realization of the KKR-CPA method, whose core algorithmic contribution — Akai's 1989 fast KKR-CPA scheme — remains foundational to the field of first-principles random-alloy theory nearly four decades later. Its principal comparative advantage is the combination of **all-electron accuracy, dependency-free portability, and genuinely fast CPA convergence**, which has made it the tool of choice for large-scale compositional screening of disordered alloys, magnetic materials, and — increasingly — high-entropy alloys and superconductors. Its principal weaknesses are structural rather than physical: an aging text-based I/O layer only now being modernized by third-party community tools, comparatively sparse English documentation, some access friction relative to fully open repositories, and the intrinsic single-site-mean-field limitations shared by all standard CPA implementations (no short-range order, no explicit local relaxation). For problems that are fundamentally about **substitutional disorder and its effect on electronic/magnetic/thermodynamic properties**, AkaiKKR remains one of the most capable and heavily used tools available; for problems demanding fine structural relaxation, strongly non-spherical potentials, or beyond-mean-field treatment of correlated disorder, full-potential or cluster-CPA codes are the more natural choice.

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the AkaiKKR (MACHIKANEYAMA) 	All-electron Korringa-Kohn-Rostoker (KKR) Green's-function DFT code, well suited to disordered alloys via the coherent potential approximation. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
