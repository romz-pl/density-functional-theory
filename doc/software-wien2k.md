# WIEN2k: A Comprehensive Review

## 1. Overview

**WIEN2k** is an all-electron, full-potential electronic structure code for solids based on the **(linearized) augmented plane wave plus local orbitals [(L)APW+lo]** method within the framework of density functional theory (DFT). It is developed and maintained at the **Institute of Materials Chemistry, TU Wien** (Vienna University of Technology), Austria, under the leadership of Peter Blaha and Karlheinz Schwarz, with major contributions from G. K. H. Madsen, D. Kvasnicka, J. Luitz, R. Laskowski, F. Tran, and L. D. Marks, among others.

WIEN2k solves the Kohn–Sham equations self-consistently for periodic solids, treating **all electrons** (core and valence) explicitly and relativistically, without recourse to pseudopotentials. It is widely regarded as one of the most accurate general-purpose DFT codes available for solids, and is frequently used as a benchmark/reference code against which pseudopotential and other all-electron methods are validated.

- **Initial release:** 1990 (as the original WIEN code); successive versions WIEN93, WIEN97, WIEN2k
- **Current name origin:** "2k" refers to the year 2000, when the Vienna-developed WIEN code was renamed WIEN2k
- **Latest stable release:** WIEN2k_24.1 (August 2024)
- **Language:** Fortran 90 (with C/C++ utility components and Perl/shell scripting for the driver layer)
- **Operating systems:** Linux/Unix (no native Windows support; runs under WSL)
- **User base:** Licensed to well over 3,000–3,400 research groups worldwide; associated publications number in the tens of thousands of citations

## 2. Licensing Model

WIEN2k is **not free/open-source software**; it is distributed under a **proprietary, source-available license** administered by TU Wien.

| License tier | Approximate fee (one-time) |
|---|---|
| Academic institutions | € 400 |
| Governmental laboratories | € 1,000 |
| Commercial / industrial users | € 4,000 |

Key licensing characteristics:
- The license is a **permanent (perpetual), one-time purchase** rather than a subscription, and includes the **full source code**, compiled distribution, and free access to subsequent point-release updates within a major version line.
- Licenses are typically issued as **group/site licenses**, valid for use by members of a single research group (PI, postdocs, and doctoral students working under them), with no formal limit on the number of machines used, though each independent group/PI generally needs its own license.
- Because the source is distributed, local compilation against optimized numerical libraries (Intel MKL, ScaLAPACK, FFTW, etc.) is expected and necessary to achieve good performance.
- Redistribution of the code or derivative binaries outside the licensed group is not permitted under the license terms.
- Some HPC centers/clusters require users to submit proof of an individual or group license before granting access to a shared installation.

## 3. Theoretical and Methodological Foundation

### 3.1 The APW lineage
WIEN2k's basis-set methodology traces a well-documented historical lineage:

1. **APW (Augmented Plane Wave), J. C. Slater (1937):** the unit cell is partitioned into non-overlapping atomic (muffin-tin) spheres and an interstitial region; plane waves are used in the interstitial region and numerical atomic-like radial functions (solutions of the radial Schrödinger equation at a given energy) are used inside the spheres. The original formulation is **energy-dependent**, leading to a nonlinear eigenvalue problem that is computationally expensive.
2. **LAPW (Linearized APW), O. K. Andersen (1975):** introduces a linearization of the energy dependence by using the radial function and its energy derivative at a fixed linearization energy, converting the problem into a standard linear eigenvalue problem. Koelling and Arbman (1975) turned this into a practical computational scheme.
3. **Local orbitals (LO), D. J. Singh (1991):** additional basis functions confined to the muffin-tin spheres, used to correctly treat semicore states and avoid linearization errors ("ghost bands") without energy-dependent basis functions.
4. **APW+lo, E. Sjöstedt, L. Nordström, D. J. Singh (2000); G. K. H. Madsen et al.:** a variant that restores the APW's superior plane-wave convergence (roughly halving the number of plane waves needed relative to LAPW) by combining APW at a fixed linearization energy with a small additional local orbital ("lo") to restore variational freedom, at the cost of extra local orbitals.
5. **Higher-derivative local orbitals (HDLOs):** address residual linearization inaccuracies in high-precision calculations by adding local orbitals built from second energy derivatives.

WIEN2k uses a **mixed basis**: APW+lo is typically applied to chemically important, hard-to-converge low-*l* channels (s, p, d, sometimes f depending on the atom), while standard LAPW+LO is used for higher angular momenta, combining the convergence benefits of APW+lo with the simplicity of LAPW for less critical channels. The authors describe the full scheme informally as "(L)APW+lo+LO+HDLO," abbreviated **APW+lo**.

### 3.2 Full-potential, all-electron, relativistic treatment
- **Full potential:** unlike shape-approximated (muffin-tin or ASA) methods, WIEN2k expands both the electron density and the Kohn–Sham potential without shape approximation — as Fourier series in the interstitial region and as lattice harmonics (symmetry-adapted combinations of spherical harmonics) times radial functions inside the spheres.
- **All-electron:** core states (from 1s upward) are treated together with valence states; no pseudopotential approximation is used. Core states are solved fully relativistically via the radial Dirac equation; valence states are treated either scalar-relativistically (mass–velocity correction and Darwin term) or, optionally, fully relativistically with spin-orbit coupling included via a second-variational treatment.
- **Basis-set convergence** is governed primarily by two parameters: $R_{MT}K_{max}$ (product of the smallest muffin-tin radius and the plane-wave cutoff) and $l_{max}$ (angular-momentum cutoff for the spherical-harmonic expansion inside the spheres), plus $G_{max}$ for the potential/density Fourier expansion.

### 3.3 Exchange-correlation functionals
WIEN2k supports a broad range of exchange-correlation treatments, including:
- **LDA** (Local Density Approximation), e.g. Perdew–Wang or Ceperley–Alder–Perdew–Zunger parametrizations
- **GGA** (Generalized Gradient Approximation), notably **PBE** (Perdew–Burke–Ernzerhof) and **PBEsol**, as well as the Engel–Vosko (EV-GGA) functional optimized for exchange potentials
- **meta-GGA** functionals (e.g. SCAN and related)
- The **Tran–Blaha modified Becke–Johnson potential (TB-mBJ)**, introduced in 2009 and incorporated into WIEN2k around 2010–2011, a semi-local exchange potential (not derived from an energy functional, so it cannot be used for force calculations) that yields band gaps for semiconductors and insulators approaching GW-level accuracy at a fraction of the computational cost. Reparametrizations (e.g. by Koller et al.) are also supported.
- **DFT+U** (Hubbard correction) for correlated *d*/*f* electron systems
- **Hybrid functionals** (e.g. PBE0, HSE-type screened hybrids), at higher computational cost
- **Nonlocal van der Waals (vdW) correlation functionals**, implemented via an FFT-based scheme (Román-Pérez–Soler approach) adapted for the all-electron density near the nucleus

### 3.4 Solving strategy
The self-consistency (SCF) cycle in WIEN2k proceeds through discrete executable programs orchestrated by shell/Perl scripts (`lapw0`, `lapw1`, `lapw2`, `lcore`, `mixer`, etc.), each responsible for one stage: potential generation, secular-equation setup and diagonalization, valence density generation, core-state solution, and charge-density mixing. This modular, file-based pipeline architecture (rather than a monolithic executable) is a distinguishing structural feature of the code and facilitates task-level and k-point parallelization.

## 4. Computational Capabilities and Properties

WIEN2k can compute a very broad range of ground-state and derived properties, including:

- **Electronic structure:** self-consistent total energy, band structure, density of states (total and orbital/atom-projected), Fermi surfaces
- **Structural properties:** equilibrium lattice parameters, atomic forces and structural (geometry) optimization, bulk modulus and elastic constants (via companion packages), equations of state
- **Magnetic properties:** collinear and non-collinear spin-polarized calculations, spin-orbit coupling, orbital moments, hyperfine fields
- **Spectroscopic properties:** X-ray absorption/emission spectra, electric field gradients (EFG), NMR chemical shielding tensors, Mössbauer isomer shifts
- **Optical properties:** linear (and some nonlinear) optical response, dielectric functions, absorption/reflectivity spectra
- **Electric polarization:** Berry-phase-based calculation of macroscopic electric polarization (ferroelectrics)
- **Transport properties:** via interfacing with external Boltzmann transport codes (e.g., BoltzTraP/BoltzTraP2) for Seebeck coefficient, electrical/thermal conductivity, thermoelectric figure of merit
- **Vibrational/phonon properties:** via companion packages implementing the direct (frozen-phonon/supercell) method
- **Wannier functions:** via the WIEN2WANNIER interface to Wannier90, enabling band interpolation, transport calculations, and input generation for DMFT (dynamical mean-field theory)
- **Many-body extensions:** GW approximation and Bethe–Salpeter equation (BSE) implementations for quasiparticle band structures and optical excitations (with higher computational cost than semilocal functionals)

## 5. Software Architecture and Usability

- **Structural input:** the crystal structure is described in a `case.struct` file (space group, lattice parameters, atomic positions, muffin-tin radii); auxiliary converters (e.g., `cif2struct`) allow import from CIF and other formats.
- **Symmetry handling:** automated detection and use of all 230 space groups, for both centrosymmetric and non-centrosymmetric structures, reducing computational cost via symmetry-adapted basis functions and irreducible Brillouin-zone sampling.
- **w2web:** a browser-based graphical user interface (GUI) that allows users to set up, run, and monitor calculations remotely, including SCF cycle configuration, convergence criteria, and parallelization settings.
- **Parallelization:** WIEN2k supports multiple parallelization strategies — MPI-based k-point parallelization, parallel matrix diagonalization via ScaLAPACK, and OpenMP shared-memory parallelization within nodes — allowing scaling from workstations to large HPC clusters.
- **External library integration:** performance-critical linear algebra can leverage optimized vendor libraries (Intel MKL, ScaLAPACK, FFTW, ELPA, etc.).
- **Companion/auxiliary codes** commonly used alongside WIEN2k include:
  - **BoltzTraP / BoltzTraP2** — semiclassical Boltzmann transport properties
  - **Gibbs2** — quasi-harmonic thermodynamic properties
  - **XCrySDen** — structure and orbital/Wannier density visualization
  - **WIEN2WANNIER** — interface to Wannier90
  - **PHONON / phonon-related packages** — lattice-dynamical properties
  - **ELAST / ortho-elastic** and related third-party packages — elastic constants for various crystal symmetries
  - **AIM (Atoms in Molecules)** — Bader-type charge-density topological analysis

## 6. Typical Application Domains

WIEN2k is used extensively across condensed-matter physics, materials science, mineralogy/geophysics, and chemistry, including:

- Semiconductors and insulators (band-gap prediction, especially via TB-mBJ or hybrid functionals)
- Transition-metal oxides and strongly correlated systems (with DFT+U or hybrid corrections)
- Actinide and lanthanide compounds, requiring relativistic and spin-orbit treatments
- Half-Heusler alloys and thermoelectric materials
- Perovskites and halide perovskites (photovoltaics, optoelectronics)
- Magnetic materials and spintronics (hyperfine fields, magnetic anisotropy)
- Surfaces, thin films, and 2D materials (via slab geometries)
- High-pressure mineral physics and geophysics
- NMR and Mössbauer spectroscopy interpretation via first-principles EFG and hyperfine parameter calculations

## 7. Strengths and Limitations

**Strengths**
- Among the most accurate all-electron, full-potential DFT methods available for solids; frequently used as an accuracy benchmark against pseudopotential codes in cross-code validation studies
- Handles heavy elements (lanthanides, actinides) and core-level/spectroscopic properties naturally, since all electrons are treated explicitly
- Rich and mature ecosystem of validated exchange-correlation functionals and post-processing tools accumulated over more than three decades of development
- Robust, automated symmetry handling across all space groups

**Limitations**
- Proprietary licensing (cost and source-distribution restrictions) contrasts with fully open-source alternatives (e.g., Quantum ESPRESSO, ABINIT)
- Fortran/script-driven, file-based architecture has a steeper learning curve than some modern, more monolithic packages
- No native Windows support
- All-electron LAPW-type methods scale less favorably with system size than pseudopotential plane-wave methods for very large unit cells, limiting routine applicability to systems of many hundreds to low thousands of atoms without significant computational resources
- mBJ and related semi-local potentials, while good for band gaps, are not exchange-correlation *energy* functionals, so atomic forces are not available in that mode (structural relaxation must be done with LDA/GGA first)

## 8. Publications on the Underlying Theory and Method

The following publications document the theoretical foundations, methodological developments, and reference implementations underlying WIEN2k, in roughly chronological/developmental order.

### 8.1 Foundational quantum-mechanical and DFT background
- P. Hohenberg and W. Kohn, "Inhomogeneous Electron Gas," *Phys. Rev.* **136**, B864 (1964).
- W. Kohn and L. J. Sham, "Self-Consistent Equations Including Exchange and Correlation Effects," *Phys. Rev.* **140**, A1133 (1965).

### 8.2 The APW / LAPW / APW+lo method lineage
- J. C. Slater, "Wave Functions in a Periodic Potential," *Phys. Rev.* **51**, 846 (1937). — original APW method.
- O. K. Andersen, "Linear methods in band theory," *Phys. Rev. B* **12**, 3060 (1975). — linearization of the APW method (LAPW).
- D. D. Koelling and G. O. Arbman, "Use of energy derivative of the radial solution in an augmented plane wave method: application to copper," *J. Phys. F: Met. Phys.* **5**, 2041 (1975). — practical LAPW implementation.
- E. Wimmer, H. Krakauer, M. Weinert, and A. J. Freeman, "Full-potential self-consistent linearized-augmented-plane-wave method for calculating the electronic structure of molecules and surfaces: O2 molecule," *Phys. Rev. B* **24**, 864 (1981). — full-potential LAPW formulation.
- M. Weinert, "Solution of Poisson's equation: Beyond Ewald-type methods," *J. Math. Phys.* **22**, 2433 (1981). — full-potential Poisson-equation solution technique used in FP-LAPW.
- D. J. Singh, "Ground-state properties of lanthanum: Treatment of extended-core states," *Phys. Rev. B* **43**, 6388 (1991); D. J. Singh and H. Krakauer, "H-point phonon in molybdenum: Superlinear convergence with the LAPW local orbital method," *Phys. Rev. B* **43**, 1441 (1991). — introduction of local orbitals (LO) to the LAPW basis.
- E. Sjöstedt, L. Nordström, and D. J. Singh, "An alternative way of linearizing the augmented plane-wave method," *Solid State Commun.* **114**, 15 (2000). — the APW+lo method.
- G. K. H. Madsen, P. Blaha, K. Schwarz, E. Sjöstedt, and L. Nordström, "Efficient linearization of the augmented plane-wave method," *Phys. Rev. B* **64**, 195134 (2001). — mixed LAPW/APW+lo basis efficient implementation.

### 8.3 Textbook and pedagogical references
- D. J. Singh and L. Nordström, *Planewaves, Pseudopotentials, and the LAPW Method*, 2nd ed. (Springer, New York, 2006).

### 8.4 Reviews and the WIEN2k reference/user-guide publications
- P. Blaha, K. Schwarz, P. Sorantin, and S. B. Trickey, "Full-potential, linearized augmented plane wave programs for crystalline systems," *Comput. Phys. Commun.* **59**, 399 (1990). — original WIEN code description.
- K. Schwarz, "DFT calculations of solids with LAPW and WIEN2k," *J. Solid State Chem.* **176**, 319–328 (2003).
- K. Schwarz, P. Blaha, and S. B. Trickey, "Electronic structure of solids with WIEN2k," *Mol. Phys.* **108**, 3147–3166 (2010).
- P. Blaha, K. Schwarz, G. K. H. Madsen, D. Kvasnicka, J. Luitz, R. Laskowski, F. Tran, and L. D. Marks, *WIEN2k: An Augmented Plane Wave Plus Local Orbitals Program for Calculating Crystal Properties* (Vienna University of Technology, Austria, 2018). — the official WIEN2k user manual/reference.
- P. Blaha, K. Schwarz, F. Tran, R. Laskowski, G. K. H. Madsen, and L. D. Marks, "WIEN2k: An APW+lo program for calculating the properties of solids," *J. Chem. Phys.* **152**, 074101 (2020). — the primary, most current, and most cited overview/reference paper for WIEN2k (version 19), covering methodology, features, and available external libraries/functionals in detail.

### 8.5 Exchange-correlation functionals used in / with WIEN2k
- J. P. Perdew, K. Burke, and M. Ernzerhof, "Generalized Gradient Approximation Made Simple," *Phys. Rev. Lett.* **77**, 3865 (1996). — the PBE-GGA functional.
- E. Engel and S. H. Vosko, "Exact exchange-based quasiparticle calculations," *Phys. Rev. B* **47**, 13164 (1993). — EV-GGA exchange potential.
- A. D. Becke and E. R. Johnson, "A simple effective potential for exchange," *J. Chem. Phys.* **124**, 221101 (2006). — original Becke–Johnson (BJ) exchange potential.
- F. Tran and P. Blaha, "Accurate Band Gaps of Semiconductors and Insulators with a Semilocal Exchange-Correlation Potential," *Phys. Rev. Lett.* **102**, 226401 (2009). — introduction of the Tran–Blaha modified Becke–Johnson (TB-mBJ) potential.
- D. Koller, F. Tran, and P. Blaha, "Merits and limits of the modified Becke-Johnson exchange potential," *Phys. Rev. B* **83**, 195134 (2011).
- D. Koller, F. Tran, and P. Blaha, "Improving the modified Becke-Johnson exchange potential," *Phys. Rev. B* **85**, 155109 (2012). — reparametrized mBJ coefficients.
- Y.-S. Kim, M. Marsman, G. Kresse, F. Tran, and P. Blaha, "Towards efficient band structure and effective mass calculations for III-V direct band-gap semiconductors," *Phys. Rev. B* **82**, 205212 (2010).
- J. A. Camargo-Martínez and R. Baquero, "Performance of the modified Becke-Johnson potential for semiconductors," *Phys. Rev. B* **86**, 195106 (2012).
- H. Jiang, "Band gaps from the Tran-Blaha modified Becke-Johnson approach: A systematic investigation," *J. Chem. Phys.* **138**, 134115 (2013).

### 8.6 Relativistic treatment and spin-orbit coupling
- D. D. Koelling and B. N. Harmon, "A technique for relativistic spin-polarised calculations," *J. Phys. C: Solid State Phys.* **10**, 3107 (1977). — scalar-relativistic approximation.
- P. Novák, "Spin-orbit coupling in the LAPW method" (internal WIEN2k documentation/report describing the second-variational spin-orbit implementation used in the code).

### 8.7 Related post-processing and extended-method publications
- G. K. H. Madsen and D. J. Singh, "BoltzTraP. A code for calculating band-structure dependent quantities," *Comput. Phys. Commun.* **175**, 67 (2006). — Boltzmann transport post-processor commonly interfaced with WIEN2k.
- J. Kuneš, R. Arita, P. Wissgott, A. Toschi, H. Ikeda, and K. Held, "Wien2wannier: From linearized augmented plane waves to maximally localized Wannier functions," *Comput. Phys. Commun.* **181**, 1888 (2010). — Wannier-function interface.
- G. Román-Pérez and J. M. Soler, "Efficient Implementation of a van der Waals Density Functional: Application to Double-Wall Carbon Nanotubes," *Phys. Rev. Lett.* **103**, 096102 (2009). — FFT-based vdW functional method adapted in WIEN2k.

### 8.8 Independent benchmarking / accuracy-verification studies
- K. Lejaeghere, G. Bihlmayer, T. Björkman, *et al.*, "Reproducibility in density functional theory calculations of solids," *Science* **351**, aad3000 (2016). — large cross-code benchmark including WIEN2k as a reference all-electron method.
- E. Ambrosch-Draxl and J. O. Sofo, "Linear optical properties of solids within the full-potential linearized augmented planewave method," *Comput. Phys. Commun.* **175**, 1 (2006). — optical-property implementation validation.

---

*Note: Publication reference details (volume/page numbers) are drawn from the standard literature record as commonly cited in WIEN2k-related papers and the official WIEN2k documentation; users are encouraged to verify exact bibliographic details against the official WIEN2k reference page (susi.theochem.tuwien.ac.at) or the publisher record before formal citation.*

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the WIEN2k 	Widely used all-electron full-potential linearized augmented plane-wave plus local orbitals (FP-LAPW+lo) code for solids; commercial/academic license. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
