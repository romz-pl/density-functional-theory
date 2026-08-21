# *Exciting*: A Full-Potential All-Electron LAPW+lo Code for Excited-State Properties

**An exhaustive technical review, with emphasis on GW, the Bethe–Salpeter Equation (BSE), and TDDFT**

---

## 1. Overview

*exciting* is an open-source, all-electron, full-potential electronic-structure package built on the (linearized) augmented planewave plus local-orbital [(L)APW+lo] method. This basis set is widely regarded as the most numerically precise scheme for solving the Kohn–Sham (KS) equations of density-functional theory (DFT), and *exciting* has demonstrated micro-Hartree-level total-energy precision. Because it treats **all electrons** — core, semi-core, and valence — on an equal footing without pseudopotentials or shape approximations to the potential, density, or wavefunctions, the code is applicable to any element of the periodic table and gives direct, unbiased access to core-electron physics.

As the name suggests, *exciting*'s defining scientific focus is **excited-state and spectroscopic properties**. Beyond ground-state DFT, it implements:

- **Time-dependent DFT (TDDFT)** — both linear-response (LR-TDDFT) and real-time (RT-TDDFT), including Ehrenfest dynamics
- **The GW approximation** of many-body perturbation theory (MBPT) — one-shot G₀W₀ and quasiparticle self-consistent GW (QSGW) — for quasiparticle band structures
- **The Bethe–Salpeter equation (BSE)** — for optical and core-level absorption spectra with excitonic (electron–hole) effects, including non-equilibrium and dynamically screened variants
- Resonant inelastic X-ray scattering (RIXS), pump–probe spectroscopy, exciton–phonon coupling, density-functional perturbation theory (DFPT) for phonons and electron–phonon coupling, Wannier interpolation, constrained DFT, and more

The code is developed primarily by the group of Claudia Draxl (Humboldt-Universität zu Berlin) together with an international collaboration, is written in modern Fortran (2018), and is distributed under the GNU GPL at [exciting-code.org](https://exciting-code.org/) and on [GitHub](https://github.com/exciting/exciting). Releases are named after chemical elements; the current version is *exciting*-sodium.

---

## 2. The LAPW+lo Basis and Ground-State Framework

### 2.1 Basis-set construction

*exciting* partitions the unit cell into non-overlapping **muffin-tin (MT) spheres** centered on the atoms and an **interstitial region**. In the interstitial region, basis functions are planewaves; inside each MT sphere, they are expanded as atomic-like radial functions times spherical harmonics, smoothly matched to the planewaves at the sphere boundary:

$$
\phi_{\mathbf{G}}^{\mathbf{k}}(\mathbf{r}) =
\begin{cases}
\dfrac{1}{\sqrt{\Omega}} e^{i(\mathbf{G}+\mathbf{k})\cdot\mathbf{r}}, & \mathbf{r}\in I,\\[4pt]
\sum_{\ell m \xi} \mathcal{A}_{\ell m \alpha \xi}^{\mathbf{G}+\mathbf{k}}\, u_{\ell\alpha\xi}(r_\alpha;\epsilon_{\ell\alpha})\, Y_{\ell m}(\hat{\mathbf r}_\alpha), & \mathbf{r}\in \mathrm{MT}_\alpha .
\end{cases}
$$

The radial functions $u_{\ell\alpha\xi}$ solve a radial Schrödinger-like equation at fixed *linearization energies* $\epsilon_{\ell\alpha}$; energy derivatives of these functions (the LAPW ingredient) or additional **local orbitals (LOs)** confined entirely to one MT sphere (the "+lo" ingredient) provide the variational flexibility needed to (i) avoid the nonlinear eigenvalue problem of the original APW method, (ii) accurately describe semi-core states, and (iii) systematically improve high-angular-momentum and high-energy components essential for excited-state theory. *exciting* automatically places linearization energies via a Wigner–Seitz prescription, avoiding material-specific tuning.

Basis-set completeness has two independent, systematically improvable knobs — the planewave cutoff (interstitial region) and the LO set (MT region) — and *exciting* has recently introduced a **dual basis self-validation (DBSV)** procedure and automated, transferable cutoff/MT-radius selection to remove much of the historical reliance on user experience for basis construction.

Core states are obtained separately, by solving the **radial Dirac equation** inside the MT spheres self-consistently with the same KS potential, giving relativistic effects (including spin–orbit coupling for core levels) automatically and making core-level spectroscopy a native, first-class feature rather than an add-on.

### 2.2 Exchange-correlation functionals

*exciting* supports the full hierarchy of Jacob's ladder relevant to solids:

- **LDA / GGA** (direct or via the bundled Libxc library)
- **meta-GGAs** (SCAN, r2SCAN, TASK, TPSS, HLE17, and any other Libxc mGGA), handled via a self-consistent GGA-seeded scheme because the kinetic-energy-density term makes the radial-solver problem non-multiplicative
- **DFT-1/2**, a Slater-transition-state-inspired correction that reaches near-GW/hybrid accuracy for band gaps at semilocal cost, and which also serves as a low-cost, well-behaved starting point for G₀W₀
- **Hybrid functionals** (PBE0, HSE06) solved via the generalized Kohn–Sham (gKS) equations, with **two independent implementations** of the Fock-exchange operator: a mixed-product-basis projection and an **adaptively compressed exchange (ACE)** scheme that removes spurious empty-state dependence and, combined with a hybrid-consistent radial solver, recovers micro-Hartree total-energy precision
- **Optimized effective potential (OEP)** for orbital-dependent local potentials
- Interfaces for van-der-Waals corrections and machine-learned functionals via community libraries

### 2.3 Spin–orbit coupling: SVLO

For heavy elements with strong SOC, *exciting* implements **second variation with local orbitals (SVLO)**, which augments the conventional second-variational SOC treatment with explicit (including Dirac-type $p_{1/2}$) local orbitals near the nucleus. This slashes the number of unoccupied scalar-relativistic states needed for convergence (e.g., ~500 vs. several thousand for γ-CsPbI₃), giving a demonstrated ~3.6× speedup at equal precision. The same idea has been extended to accelerate SOC-BSE calculations (Section 4.4 below).

### 2.4 Constrained DFT, Wannier interpolation, symmetry

- **Constrained DFT (cDFT)** fixes non-equilibrium occupation numbers to mimic excitonic or pump-probe-generated carrier distributions at a fraction of BSE's cost, especially useful for spatially confined electron–hole pairs (e.g., the supercell core-hole approach) and for late-time-delay pump–probe spectra.
- **Maximally localized Wannier functions (MLWFs)** can be built from KS/gKS or GW input (isolated and entangled band subspaces, with automatic initial-guess generation from LOs — no manual projector input required), enabling efficient interpolation of bands, DOS, group velocities, and electron–phonon matrix elements.
- Crystal symmetry is exploited throughout: potential/density are expanded in symmetry-adapted lattice harmonics, and inversion symmetry allows a real (rather than complex) eigenvalue solver, giving roughly a 4× diagonalization speedup.

### 2.5 Lattice dynamics (DFPT)

*exciting* implements phonons via **density-functional perturbation theory** (Sternheimer-equation approach, fully symmetry-reduced over **k**/**q**) as well as finite-difference/supercell methods, including Born effective charges and the dielectric tensor for the LO–TO splitting in polar materials. Electron–phonon coupling (EPC) matrix elements feed Fan–Migdal and Debye–Waller self-energies for temperature-dependent quasiparticle renormalization (zero-point renormalization, band-gap temperature dependence), Eliashberg-function-based superconductivity parameters, and — going beyond standard perturbation theory — polaron and exciton-polaron physics.

---

## 3. GW: Quasiparticle Band Structures

### 3.1 Motivation and formalism

KS eigenvalues are Lagrange multipliers, not quasiparticle (QP) energies, and can deviate substantially from the true many-body excitation spectrum — most visibly in the severe band-gap underestimation of semilocal DFT. The **GW approximation** to Hedin's self-energy, $\Sigma = iGW$ (vertex $\Gamma \approx 1$), is the standard MBPT route to QP band structures, connecting directly to photoemission (XPS/ARPES) experiments.

*exciting* implements the standard **one-shot G₀W₀** scheme perturbatively on top of a DFT/gKS reference, using the linearized QP equation

$$
\mathrm{Re}[\epsilon_{n\mathbf{k}}^{\rm QP}] \approx \epsilon_{n\mathbf{k}} + Z_{n\mathbf{k}}\,
\langle \psi_{n\mathbf{k}} | \mathrm{Re}[\Sigma(\mathbf{r},\mathbf{r}',\epsilon_{n\mathbf{k}})] - v_{\rm xc}(\mathbf{r})\,\delta(\mathbf{r}-\mathbf{r}') | \psi_{n\mathbf{k}} \rangle .
$$

Any GS/gKS functional (LDA, GGA, mGGA, DFT-1/2, hybrids) can serve as the G₀W₀ starting point.

### 3.2 The mixed-product basis

Two-particle quantities (polarizability, dielectric matrix, screened interaction) are expanded not in the LAPW+lo basis itself but in an **auxiliary mixed-product basis** (spherical harmonics inside MT spheres, planewaves in the interstitial region) specifically tailored to represent *products* of KS wavefunctions,

$$
\psi_{n\mathbf{k}}(\mathbf{r})\,\psi^*_{m\mathbf{k}-\mathbf{q}}(\mathbf{r}) = \sum_i M^i_{nm}(\mathbf{k},\mathbf{q})\,\mathcal{B}_i^{\mathbf{q}}(\mathbf{r}).
$$

The RPA dielectric matrix, screened Coulomb interaction $W = \varepsilon^{-1}v$, and the exchange/correlation parts of the self-energy are all built from this basis, following an approach that traces back to Aryasetiawan/Gunnarsson-style all-electron mixed-basis GW and its refinement for LAPW+lo.

### 3.3 High-performance computing and numerical precision

- A **task-based GW workflow** parallelizes over **k**/**q** points, unoccupied-state bands, and frequencies, with GPU offloading (vendor BLAS + OpenMP-offload kernels) giving **4–10× speedups**.
- Because excited-state properties are especially sensitive to basis-set incompleteness, *exciting* employs extensive **high-energy local orbitals (HELOs)** for GW and an **incomplete-basis-set correction (IBC)**.
- Accurate treatment of the long-range $\mathbf{q}\to0$ limit of the dielectric matrix/screened interaction is essential for anisotropic and low-dimensional systems; *exciting* interfaces with the code-agnostic **IDieL** library (spherical-harmonic expansion of head/wings/body over a scaled Brillouin zone) for a parameter-free, fully anisotropic evaluation, extending earlier 2D-specific Coulomb-cutoff treatments (validated for MoS₂, black phosphorus).
- The **expansion-and-addition-screening (EAS)** method exploits the additive nature of polarizability in weakly bound van der Waals heterostructures (e.g., organic-molecule/MoS₂ interfaces), cutting G₀W₀ polarizability cost by ~56% (25% total speedup) and BSE polarizability cost by ~69% (46% total speedup).

### 3.4 Beyond G₀W₀: quasiparticle self-consistent GW (QSGW)

To remove G₀W₀'s well-known starting-point dependence, *exciting* implements **QSGW**, in which the dynamical self-energy is replaced by an optimal *static, non-local* exchange-correlation potential

$$
v_{\rm xc}^{\rm opt}(\mathbf{k}) = \tfrac12\sum_{nl}|\psi_{n\mathbf{k}}\rangle\left[\mathrm{Re}(\Sigma_{nl\mathbf{k}}(\epsilon_{n\mathbf{k}})) + \mathrm{Re}(\Sigma_{nl\mathbf{k}}(\epsilon_{l\mathbf{k}}))\right]\langle\psi_{l\mathbf{k}}|,
$$

iterated to self-consistency. This eliminates starting-point dependence and produces improved electron densities that systematically reduce the delocalization error of semilocal functionals (demonstrated for Ar and CaO). Being all-electron, *exciting*'s QSGW captures core–valence exchange-correlation contributions to the self-energy that pseudopotential approaches must neglect — an error that self-consistent iteration can otherwise amplify.

### 3.5 Cross-code validation

Benchmark studies comparing G₀W₀ band gaps across Abinit, *exciting*, FHI-aims, and GPAW for a set of solids find agreement within ~0.1–0.3 eV overall, with the two **all-electron codes (*exciting* and FHI-aims) showing the closest mutual agreement** (≲15 meV at the DFT level, ~0.1 eV at G₀W₀), underlining the reference-quality role of all-electron LAPW+lo GW.

---

## 4. Bethe–Salpeter Equation (BSE): Optical and Core Excitations

### 4.1 Formalism

The BSE is *exciting*'s workhorse for **neutral (electron–hole) excitations** and their signatures in absorption, EELS, and X-ray spectra. It is recast as an effective two-particle eigenvalue problem in a transition basis of occupied↔unoccupied products, with the Bethe–Salpeter Hamiltonian

$$
H^{\rm BSE}(\mathbf{q}) = \begin{pmatrix} A(\mathbf{q}) & B(\mathbf{q}) \\ -B^*(\mathbf{q}) & -A^*(\mathbf{q}) \end{pmatrix},
\qquad
A = D + \gamma_x V - W, \quad B = \gamma_x V - W,
$$

where $D$ is the diagonal QP energy-difference term, $V$ the bare exchange interaction (responsible for singlet–triplet splitting, controlled by $\gamma_x$), and $W$ the **statically screened direct interaction**. *exciting* supports the **Tamm–Dancoff approximation (TDA)**, where the anti-resonant block $B$ is dropped, *and* the **full BSE** including resonant–anti-resonant coupling — a capability not universal among BSE codes. Special cases available include the independent-particle approximation (IPA, no kernel) and RPA (W neglected).

### 4.2 Low-scaling BSE implementation

Conventional BSE matrix construction/diagonalization scales as $\mathcal{O}(N_o^2N_u^2N_{\mathbf k}^2)$ / $\mathcal{O}(N_o^3N_u^3N_{\mathbf k}^3)$, quickly becoming prohibitive. *exciting* implements **interpolative separable density fitting (ISDF)** to approximate wavefunction pair products via a small set of real-space interpolation points (found by centroidal Voronoi tessellation), combined with an **iterative Lanczos solver**. This reduces the overall scaling to $\Theta(N_o N_u N_{\mathbf k}\log N_{\mathbf k})$, reproducing direct-BSE spectra exactly (demonstrated for LiF) with speedups from ~6× up to two orders of magnitude for dense **k**-meshes or large cells.

### 4.3 Core-level BSE (the "same-footing" advantage)

Because core states are obtained from the same all-electron Dirac-equation solver as valence states, *exciting* solves the BSE for **core→conduction transitions exactly as for optical (valence→conduction) transitions**, simply by substituting core initial states into the same matrix-element machinery. This gives XANES/NEXAFS/XAS spectra spanning any absorption edge for any element, on equal footing with optical absorption — a hallmark, code-defining capability described in detail in the dedicated core-level-spectroscopy literature (Draxl & Cocchi).

### 4.4 SVLO in BSE

The SVLO basis (Section 2.3) has been extended consistently to BSE planewave and momentum matrix elements, giving substantially faster convergence of exciton binding energies with respect to the number of unoccupied states in strongly spin–orbit-coupled materials (demonstrated for γ-CsPbI₃).

### 4.5 Non-equilibrium and dynamically coupled BSE

A **non-equilibrium BSE** variant incorporates photoinduced carrier populations (from cDFT or RT-TDDFT) into the momentum matrix elements and screened interaction, enabling transient excitonic spectra for pump–probe simulations.

### 4.6 Resonant inelastic X-ray scattering (RIXS)

Built on the Kramers–Heisenberg formula, the standalone **BRIXS** package (tightly integrated with *exciting*'s workflow tools) combines a core-edge BSE calculation and an optical BSE calculation to compute RIXS cross-sections, including arbitrary incident/scattered photon polarizations, and extends naturally to *pumped* RIXS via non-equilibrium BSE.

### 4.7 Lattice screening and exciton–phonon coupling

Beyond purely electronic screening, *exciting* has implemented **phonon-assisted (Fröhlich-type) screening** of the electron–hole interaction within BSE, giving exciton-binding-energy renormalization and absorption-onset red-shifts (demonstrated for ZnS, MgO, GaN — effects on the order of tens of meV, dominated by long-range coupling to polar LO phonons). This connects to a broader exciton–phonon coupling (EXPC) module addressing self-trapped excitons and excitonic polarons.

---

## 5. TDDFT (Linear-Response and Real-Time)

### 5.1 Linear-response TDDFT (LR-TDDFT)

*exciting* implements frequency-domain LR-TDDFT with a range of exchange-correlation kernels (RPA and beyond) and full inclusion of finite-momentum-transfer effects, making it suitable for inelastic X-ray and electron-energy-loss spectroscopy at general **q**. TDDFT is generally preferred over BSE for systems with weak electron–hole interaction (small molecules, metals) because it is computationally far cheaper than MBPT while remaining formally exact in principle (subject to the xc-kernel approximation used).

### 5.2 Real-time TDDFT (RT-TDDFT)

A full **all-electron, full-potential RT-TDDFT** implementation propagates the time-dependent KS equations directly in the LAPW+lo basis, benchmarked against both *exciting*'s own LR-TDDFT and against the (pseudopotential/real-space) Octopus code, with satisfactory agreement in all cases. Demonstrated applications include:

- Laser-pulse-driven excitation dynamics in monolayer MoS₂
- Third-harmonic generation in silicon (a genuinely nonlinear, non-perturbative response only accessible in real time)
- Pump–probe simulation of the dielectric function of diamond

RT-TDDFT can also be coupled to **Ehrenfest dynamics** for combined electron–ion non-adiabatic propagation, and its photoexcited-carrier output feeds directly into non-equilibrium BSE (Section 4.5) for post-hoc excitonic/core-spectroscopy analysis of driven systems.

### 5.3 Comparative role of TDDFT vs. BSE in *exciting*

The code's documentation explicitly frames the trade-off: TDDFT (LR or RT) is the natural, cost-effective choice for weakly correlated electron–hole systems and finite-momentum-transfer spectroscopy, while BSE is reserved for systems with strong excitonic/core-hole binding — with both routes converging into shared machinery for pump–probe and non-equilibrium spectroscopy.

---

## 6. Pump–Probe Spectroscopy and Non-Equilibrium Methods

*exciting* combines **cDFT**, **RT-TDDFT**, and **non-equilibrium BSE** into a coherent pump–probe toolchain: cDFT models the picosecond-timescale, thermalized carrier distribution after a pump pulse; RT-TDDFT tracks femtosecond-timescale, non-thermal carrier occupations directly; either route supplies the excited-state occupations that non-equilibrium BSE uses to compute transient optical or (via non-equilibrium RIXS) X-ray absorption spectra. This machinery has been applied, e.g., to ultrafast dynamic Coulomb screening of X-ray core excitons in photoexcited semiconductors.

---

## 7. Computational Infrastructure

| Aspect | Details |
|---|---|
| Language | Fortran 2018 (performance-critical parts HPC-optimized) |
| Build system | CMake; GNU Fortran (≥ GCC 11), Intel oneAPI IFX, Cray, AMD Flang supported |
| Required libraries | FFTW3-compatible FFT library, BLAS/LAPACK (OpenBLAS, MKL, Cray LibSci, …) |
| Bundled | Libxc, FoX, BSPLINE-FORTRAN |
| Optional acceleration | ELPA, ScaLAPACK, HDF5 (I/O), SIRIUS interface |
| GPU support | Intel oneAPI, CUDA/cuBLAS, ROCm; unified-memory-optimized builds available |
| Parallelism | MPI + OpenMP, with task-based GPU-offloaded workflows for GW |
| License | GNU GPL |
| Distribution | [exciting-code.org](https://exciting-code.org/), [github.com/exciting/exciting](https://github.com/exciting/exciting) |
| Testing | Three-tier CI: unit tests, regression tests, full-tutorial workflow tests |
| Tutorials | Jupyter-notebook-based, browser-executable, CI-integrated |
| Ancillary tools | excitingtools, excitingworkflow, excitingscripts (Python); BRIXS/pyBRIXS (RIXS); ElaStic (elastic constants); CELL (cluster expansion); interfaces to SIRIUS and to the coupled-cluster code Cc4s; interface to elphbolt for Boltzmann transport beyond the relaxation-time approximation |

---

## 8. Publications on the Underlying Theory and Implementation

The following papers document the theoretical foundations and code implementations referenced throughout this review. The *exciting* GitHub repository explicitly requests that users cite the general package paper [1] plus any feature-specific papers relevant to their calculation.

### Foundational / general package

1. Gulans, A.; Kontur, S.; Meisenbichler, C.; Nabok, D.; Pavone, P.; Rigamonti, S.; Sagmeister, S.; Werner, U.; Draxl, C. **"exciting — a full-potential all-electron package implementing density-functional theory and many-body perturbation theory."** *J. Phys.: Condens. Matter* **26**, 363202 (2014). — The primary code-reference paper.
2. Raya-Moreno, M.; Buccheri, A.; Dasch, N. A.; et al.; Draxl, C. **"An exciting approach to theoretical spectroscopy."** arXiv:2601.11388 (2026). — The comprehensive, up-to-date methods review underlying most of the present document.

### Ground-state methodology, basis sets, and precision

3. Jansen, H. J. F.; Freeman, A. J. **"Total-energy full-potential linearized augmented-plane-wave method for bulk solids: Electronic and structural properties of tungsten."** *Phys. Rev. B* **30**, 561 (1984). — Foundational LAPW reference.
4. Sjöstedt, E.; Nordström, L.; Singh, D. J. **"An alternative way of linearizing the augmented plane-wave method."** *Solid State Commun.* **114**, 15 (2000). — APW+lo method.
5. Gulans, A.; Kontur, S.; Meisenbichler, C.; et al. (see [1]).
6. Huhn, W. P.; Blum, V. (and related benchmarking works) — precision benchmarking of all-electron LAPW+lo vs. other basis types.

### Hybrid functionals, mGGA, DFT-1/2

7. Betzinger, M.; Friedrich, C.; Blügel, S. **"Hybrid functionals within the FLAPW method: Implementation and applications of PBE0."** *Phys. Rev. B* **81**, 195117 (2010).
8. Lin, L. **"Adaptively compressed exchange operator."** *J. Chem. Theory Comput.* **12**, 2242 (2016). — ACE method adopted for hybrid-functional exchange in *exciting*.
9. Ferreira, L. G.; Marques, M.; Teles, L. K. **"Approximation to density functional theory for the calculation of band gaps of semiconductors."** *Phys. Rev. B* **78**, 125116 (2008). — DFT-1/2.
10. Fischer, D.; et al. **"Probing the LDA-1/2 method as a starting point for G0W0 calculations."** *Phys. Rev. B* **94**, 235141 (2016).

### Spin–orbit coupling (SVLO)

11. Recent *exciting* SVLO implementation paper describing the second-variation-with-local-orbitals scheme and its extension to BSE (see review, Ref. [30] therein).

### GW implementation

12. Jiang, H.; Gómez-Abal, R. I.; Rinke, P.; Scheffler, M. **"Localized and itinerant states in lanthanide oxides united by GW@LDA+U."** *Phys. Rev. Lett.* **102**, 126403 (2009); and related mixed-product-basis GW methodology.
13. Nabok, D.; Gulans, A.; Draxl, C. **"Accurate all-electron G0W0 quasiparticle energies employing the full-potential augmented plane-wave method."** *Phys. Rev. B* **94**, 035418 (2016). — Core *exciting* G₀W₀ implementation paper.
14. Kotani, T.; van Schilfgaarde, M. **"All-electron GW approximation with the mixed basis expansion based on the full-potential LMTO method."** *Solid State Commun.* **121**, 461 (2002). — Mixed-product-basis GW precursor.
15. Friedrich, C.; Blügel, S.; Schindlmayr, A. **"Efficient implementation of the GW approximation within the all-electron FLAPW method."** *Phys. Rev. B* **81**, 125102 (2010).
16. QSGW implementation in *exciting*: van Schilfgaarde, M.; Kotani, T.; Faleev, S. **"Quasiparticle self-consistent GW theory."** *Phys. Rev. Lett.* **96**, 226402 (2006), as adapted and implemented in *exciting* (Ref. [34] of the 2026 review).
17. Azizi, M.; et al. **"Precision benchmarks for solids: G0W0 calculations with different basis sets."** arXiv:2411.19701 / *Comput. Mater. Sci.* (2025). — Cross-code (Abinit/*exciting*/FHI-aims/GPAW) G₀W₀ benchmark.
18. Rasmussen, F. A.; Schmidt, P. S.; Winther, K. T.; Thygesen, K. S. **"Efficient many-body calculations for two-dimensional materials using exact limits for the screened potential: Band gaps of MoS2, h-BN, and phosphorene."** *Phys. Rev. B* **94**, 155406 (2016). — 2D **q**→0 methodology adopted/extended in *exciting*.
19. IDieL library papers for anisotropic **q**→0 treatment (GreenX ecosystem).
20. Kim, M.; et al. **"Efficient G0W0 and BSE calculations of heterostructures within an all-electron framework."** arXiv:2507.17960 (2025). — Expansion-and-addition-screening (EAS) implementation.

### Bethe–Salpeter equation

21. Sagmeister, S.; Ambrosch-Draxl, C. **"Time-dependent density functional theory versus Bethe–Salpeter equation: an all-electron study."** *Phys. Chem. Chem. Phys.* **11**, 4451 (2009). — Cited requirement for BSE/TDDFT usage in *exciting*.
22. Vorwerk, C.; Aurich, B.; Cocchi, C.; Draxl, C. **"Bethe–Salpeter equation for absorption and scattering spectroscopy: implementation in the exciting code."** *Electron. Struct.* **1**, 037001 (2019) [arXiv:1904.05575]. — Core BSE implementation paper, including finite-momentum-transfer and core-spectroscopy extensions.
23. Draxl, C.; Cocchi, C. **"exciting core-level spectroscopy."** arXiv:1709.02288 — Core-level BSE methodology and applications.
24. Low-scaling (ISDF + Lanczos) BSE implementation paper (Ref. [36], [223], [227]–[228] of the 2026 review).

### Real-time and linear-response TDDFT

25. Pela, R. R.; Draxl, C. **"All-electron full-potential implementation of real-time TDDFT in exciting."** *Electron. Struct.* **3**, 037001 (2021) [arXiv:2102.02630]. — RT-TDDFT implementation, benchmarked against Octopus.
26. Sagmeister & Ambrosch-Draxl (2009), see [21], for LR-TDDFT vs. BSE.

### RIXS

27. Vorwerk, C.; Cocchi, C.; Draxl, C.; et al. — BRIXS/RIXS theoretical framework papers (Kramers–Heisenberg implementation on top of BSE), corresponding to Refs. [37]–[38], [248]–[250] of the 2026 review.

### Non-equilibrium spectroscopy / pump–probe / cDFT

28. Non-equilibrium BSE and cDFT pump–probe methodology (Refs. [115]–[116] of the 2026 review).
29. Ultrafast dynamic Coulomb screening of X-ray core excitons in photoexcited semiconductors, arXiv:2412.01945 — combined cDFT/RT-TDDFT + non-equilibrium BSE application.

### Lattice dynamics, electron–phonon coupling, exciton–phonon coupling

30. DFPT implementation in *exciting* (Ref. [22] of the 2026 review).
31. Fröhlich-type phonon-assisted BSE screening implementation (Ref. [241] of the 2026 review).

### Wannier interpolation

32. Tillack, S.; et al. **"Maximally localized Wannier functions within the (L)APW+LO method."** arXiv:1908.11156.

### Symmetry, dual basis self-validation, automated basis construction

33. Dual basis self-validation (DBSV) methodology paper (Ref. [28] of the 2026 review; "to be published").

---

## 9. Summary

*exciting* occupies a distinctive niche among electronic-structure packages: it combines the **highest-precision all-electron LAPW+lo ground-state framework** with essentially the **complete methodological hierarchy for excited-state and spectroscopic properties** — TDDFT (linear-response and real-time), G₀W₀ and QSGW quasiparticle theory, BSE (valence and core, static and non-equilibrium, with and without the Tamm–Dancoff approximation), RIXS, and exciton–phonon/lattice-screening physics — all built on the same all-electron foundation that lets valence and core states be treated on strictly equal footing. Recent development has focused heavily on removing the traditional computational bottlenecks of these methods (low-scaling ISDF-Lanczos BSE, task-based GPU-accelerated GW, expansion-and-addition screening for heterostructures, automated basis-set construction) without sacrificing the code's benchmark-grade numerical precision, positioning it as a reference tool for cross-validating other electronic-structure and MBPT implementations.

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the Exciting 	Full-potential all-electron LAPW+lo code with a strong focus on excited-state properties (GW, BSE, TDDFT).Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
