# FHI-aims: All-Electron, Numeric Atom-Centered Orbital DFT Code — A Comprehensive Review

## 1. Overview

FHI-aims (**F**ritz **H**aber **I**nstitute **a**b **i**nitio **m**olecular **s**imulations) is an all-electron, full-potential electronic structure code based on **numeric atom-centered orbitals (NAOs)**. Since its foundation in 2004, the code was designed with a clear set of goals: to be numerically precise across the periodic table, to be "all-electron" (not pseudopotential-based), and to handle periodic systems (extended solids, surfaces, nanostructures) as well as non-periodic systems (molecules and clusters) on an equal footing. It was first released to the community in 2009.

FHI-aims enables first-principles simulations with very high numerical accuracy for production calculations, with excellent scalability up to very large system sizes (thousands of atoms) and up to very large, massively parallel supercomputers (tens of thousands of CPU cores). Its accuracy is benchmarked against the best available all-electron reference methods (e.g., the Δ-test comparisons and follow-up studies), and its scalability has been demonstrated in dedicated large-system studies.

Today, FHI-aims is a worldwide community project with well over 150 individual contributors, spanning institutions such as FHI, Duke University, TU Munich, USTC Hefei, Aalto University, University of Luxembourg, TU Graz, Cardiff University, and many others.

---

## 2. Core Methodology

### 2.1 Numeric Atom-Centered Orbital (NAO) Basis

FHI-aims expands Kohn–Sham orbitals in numerically tabulated atom-centered basis functions of the form

$$\varphi_i(\mathbf{r}) = \frac{u_i(r)}{r}\, Y_{lm}(\Omega)$$

where $Y_{lm}(\Omega)$ denotes the real spherical harmonics, and the radial function $u_i(r)$ is obtained by solving a radial Schrödinger-like equation containing an effective potential $v_i(r)$ and a confining potential $v_{cut}(r)$ that switches on beyond an onset radius and forces $u_i(r)$ to strictly vanish at $r_{onset}+w$. This confinement ensures strict spatial localization of every basis function while retaining full numerical flexibility of the radial shape — unlike Gaussian- or Slater-type orbitals, the NAO radial part is not restricted to any particular analytic form, though GTOs and STOs can be represented as special (numerically tabulated) cases within the same framework.

FHI-aims uses atomic orbitals solved from a set of spherically symmetric free-atom Kohn–Sham equations as its basis functions, computed on a numerically precise one-dimensional logarithmic grid and then orthonormalized before being interpolated onto a three-dimensional overlapping atom-centered integration grid consisting of spherical shells around each atom, with points on each shell distributed following Lebedev/Delley quadrature prescriptions.

### 2.2 Hierarchical Basis Sets

FHI-aims provides preconstructed, hierarchical numeric atom-centered basis sets for elements 1–102, enabling systematic convergence from fast qualitative calculations to meV-level total-energy convergence. Preconstructed default settings ("light," "intermediate," "tight," "really_tight") bundle consistent choices for integration grids, the Hartree potential treatment, and basis cutoffs. For correlated wavefunction methods, numeric atom-centered valence correlation-consistent basis sets (NAO-VCC-*n*Z) are available for systematic convergence of many-body perturbation methods (currently for H–Ar), alongside support for standard Gaussian basis sets from quantum chemistry.

### 2.3 All-Electron, Full-Potential Treatment

Because all electrons (core and valence) are treated explicitly with no shape approximation to potentials or wavefunctions, FHI-aims avoids the approximations inherent to pseudopotential methods, while achieving all-electron accuracy at a computational cost comparable to plane-wave/pseudopotential implementations. A dedicated, validated **frozen-core approximation** is also available for cases where core relaxation is not required, improving efficiency without sacrificing controlled accuracy.

### 2.4 Hartree Potential and Resolution-of-Identity (RI) Techniques

The electrostatic (Hartree) potential is evaluated efficiently as a correction relative to the sum of free-atom potentials via multipole-expanded atom-centered components, which is central to the code's favorable scaling. For exact-exchange and correlated methods, FHI-aims employs a **resolution-of-identity (RI)** approach using auxiliary NAO basis functions (generated on-the-fly by default, via products of primary basis functions with Gram–Schmidt orthogonalization) to reduce four-center Coulomb integrals to two- and three-center quantities, enabling Hartree–Fock, hybrid DFT, MP2, RPA, and *GW* within the same NAO infrastructure.

---

## 3. Electronic-Structure Methods Implemented

| Category | Capabilities |
|---|---|
| **DFT functionals** | LDA (Perdew–Wang 1992, Perdew–Zunger 1981, Vosko–Wilk–Nusair), GGA (PBE, PBEsol, BLYP, AM05, RPBE, revPBE, PBEint), meta-GGAs, hybrid functionals (PBEh/PBE0, HSE, B3LYP-type), and range-separated/screened hybrids |
| **Exact exchange / wavefunction methods** | Hartree–Fock, MP2, coupled-cluster theory, double hybrids |
| **Many-body perturbation theory** | *GW* (G₀W₀ and self-consistent variants), Random Phase Approximation (RPA) total energies, Bethe–Salpeter Equation (BSE) for neutral excitations |
| **DFT+U** | Spherically-averaged Hubbard correction, using NAOs directly as projector functions |
| **Excited-state / dynamics methods** | Linear-response TDDFT, real-time TDDFT (RT-TDDFT), imaginary-time TDDFT (it-TDDFT) |
| **Response properties** | Real-space, all-electron perturbation theory (PT) for homogeneous electric fields — polarizabilities, dielectric constants, Raman intensities (harmonic and anharmonic) |
| **Relativistic treatments** | Scalar-relativistic approaches (ZORA-type), with development toward the full four-component Dirac equation |
| **Structure & dynamics** | Analytical gradients and stress tensors for structure relaxation and *ab initio* molecular dynamics (AIMD); Born effective charges, polarization, and topological invariants (ℤ₂) via a Berry-phase approach |
| **Embedding & environment** | QM/MM embedding with norm-conserving pseudopotentials; implicit solvation formalisms for molecules |

The DFT+U implementation in particular exploits the fact that FHI-aims's highly localized NAOs straightforwardly serve as projector functions for the Hubbard correction, yielding the necessary occupations of the correlated subspace at no additional cost, validated on bulk NiO and polaron formation at rutile TiO₂(110).

---

## 4. Numerical Performance and Scalability

- System size range up to thousand(s) of atoms, with O(N)-like scaling for the most expensive operations — the limiting factor beyond this range being the conventional O(N³) eigensolver.
- The code is seamlessly parallel in both time and memory, from desktop machines up to tens of thousands of CPUs, using the specifically optimized, massively parallel eigensolver ELPA to minimize scalability barriers, alongside the ELSI infrastructure that unifies access to eigensolvers and density-matrix solvers.
- A memory-efficient Hartree–Fock/exact-exchange implementation exploits MPI-3 intra-node shared-memory arrays.
- Recent algorithmic advances have enabled **efficient all-electron hybrid DFT for systems beyond 10,000 atoms**: restructured Hartree-potential and exact-exchange (EXX) evaluation avoids branching in inner loops and reduces subroutine calls by two orders of magnitude, benefiting both semilocal and hybrid DFT calculations, without introducing new approximations.

---

## 5. Outputs and Interoperability

Standard outputs include electron densities, Kohn–Sham orbitals, band structures, densities of states, and the electron localization function (ELF), all ready for visualization with standard tools. The code connects to a broad software ecosystem:

- **i-PI** — for advanced (path-integral) molecular dynamics
- **ASE** (Atomic Simulation Environment) — scripting/workflow interface
- **FHI-vibes** — vibrational/phonon and thermodynamics workflows
- **GIMS** (Graphical Interface for Materials Simulations) — a browser-based GUI for building inputs and running/analyzing calculations, also supporting the *exciting* code
- **ELPA / ELSI / CECAM Electronic Structure Library** — shared HPC solver infrastructure
- **NOMAD** — data repository and FAIR-data integration (FHI-aims is part of the NOMAD Laboratory ecosystem)
- **TDEP** — temperature-dependent effective potentials for anharmonic thermal transport and Raman spectroscopy

---

## 6. Licensing and Access

FHI-aims is **not open source**; it is distributed under an academic/commercial license model administered by the non-profit association **MS1P e.V.** (Molecular Simulations from First Principles), to which the Max Planck Society has transferred the rights to use and exploit the code.

- **Academic license**: covers use by a single research group at a university or non-profit research lab; payment is voluntary, with suggested academic license fees around €2,000 (up to 5 users) or €4,000 (up to 20 users), plus VAT, used to help sustain code maintenance.
- **Commercial license**: required for use in a non-academic (e.g., company) setting; individual quotes are provided on request.
- Source code access, once licensed, is via a private GitLab server; community support runs through a dedicated **Slack workspace** and the **FHI-aims Club** registration portal.
- The most recent stable release at the time of writing is **250822** (22 August 2025).

---

## 7. Typical Application Domains

- Molecular and cluster chemistry (organic and inorganic molecules, nanoclusters)
- Periodic solids: bulk crystals, surfaces, interfaces, 2D materials
- Defects and polarons in oxides (e.g., TiO₂, NiO)
- Excited-state spectroscopy: optical absorption, Raman, core-level/X-ray spectroscopy via BSE and TDDFT
- Vibrational and thermal transport properties (phonons, TDEP-based anharmonic effects)
- Ferroelectrics and topological materials (Berry-phase polarization, Z₂ invariants)
- Large biomolecular and nanostructure systems enabled by scalable hybrid DFT
- Machine-learning and multiscale workflows coupling DFT data to interatomic potentials

---

## 8. Key Publications on FHI-aims Theory and Implementation

The following list covers the primary methodological and theoretical papers underlying FHI-aims, organized by subject area.

### Foundational Method and Code Papers
- V. Blum, R. Gehrke, F. Hanke, P. Havu, V. Havu, X. Ren, K. Reuter, M. Scheffler, *"Ab initio molecular simulations with numeric atom-centered orbitals,"* **Computer Physics Communications** 180, 2175–2196 (2009).
- V. Havu, V. Blum, P. Havu, M. Scheffler, *"Efficient O(N) integration for all-electron electronic structure calculation using numeric basis functions,"* **Journal of Computational Physics** 228, 8367–8379 (2009).
- V. Blum, M. Rossi, S. Kokott, M. Scheffler, *"The FHI-aims Code: All-electron, ab initio materials simulations towards the exascale,"* in *Roadmap on Electronic Structure Codes in the Exascale Era*, **Electronic Structure** (2022/2023); arXiv:2209.12747.

### Resolution-of-Identity, Exact Exchange, and Correlated Methods
- X. Ren, P. Rinke, V. Blum, J. Wieferink, A. Tkatchenko, A. Sanfilippo, K. Reuter, M. Scheffler, *"Resolution-of-identity approach to Hartree-Fock, hybrid density functionals, RPA, MP2, and GW with numeric atom-centered orbital basis functions,"* **New Journal of Physics** 14, 053020 (2012).
- A. C. Ihrig, J. Wieferink, I. Y. Zhang, M. Ropo, X. Ren, P. Rinke, M. Scheffler, V. Blum, *"Accurate localized resolution of identity approach for linear-scaling hybrid density functionals and for many-body perturbation theory,"* **New Journal of Physics** 17, 093020 (2015).
- S. Kokott, F. Merz, Y. Yao, C. Carbogno, M. Rossi, V. Havu, M. Rampp, M. Scheffler, V. Blum, *"Efficient all-electron hybrid density functionals for atomistic simulations beyond 10,000 atoms,"* **The Journal of Chemical Physics** 161, 024112 (2024).

### Large-Scale / Linear-Scaling Solvers
- A. Marek, V. Blum, R. Johanni, V. Havu, B. Lang, T. Auckenthaler, A. Heinecke, H.-J. Bungartz, H. Lederer, *"The ELPA library — scalable parallel eigenvalue solutions for electronic structure theory and computational science,"* **Journal of Physics: Condensed Matter** 26, 213201 (2014).
- V. W.-z. Yu, F. Corsetti, A. García, W. P. Huhn, M. Jacquelin, W. Jia, B. Lange, L. Lin, J. Lu, W. Mi, A. Seifitokaldani, Á. Vázquez-Mayagoitia, C. Yang, H. Yang, V. Blum, *"ELSI: A unified software interface for Kohn–Sham electronic structure solvers,"* **Computer Physics Communications** 222, 267–285 (2018).
- V. Blum et al., linear-scaling / scalability demonstration references, e.g. **Physical Review Letters** 111, 065502 (2013) and **Computer Physics Communications** 192, 60–69 (2015).

### DFT+U and Beyond-DFT Corrections
- M. Behler, N. Bogdanov, W. R. L. Lambrecht, et al. (implementation team), *"Intricacies of DFT+U, Not Only in a Numeric Atom Centered Orbital Framework,"* **Journal of Chemical Theory and Computation** 15, 1743–1755 (2019).

### Response Properties, Perturbation Theory, and Spectroscopy
- H. Shang, N. Raimbault, P. Rinke, M. Scheffler, M. Rossi, C. Carbogno, *"All-electron, real-space perturbation theory for homogeneous electric fields: theory, implementation, and application within DFT,"* **New Journal of Physics** 20, 073040 (2018).
- C. Carbogno, N. Rybin, S. Panahian Jand, A. Akkoush, C. Mera Acosta, Z. Yuan, M. Rossi, *"Polarisation, Born Effective Charges, and Topological Invariants via a Berry-Phase Approach,"* arXiv:2501.02550 (2025).

### Time-Dependent DFT
- X.-P. Li et al., FHI-aims RT-TDDFT overview, *"Real-Time Time-Dependent Density Functional Theory within FHI-aims"* (conference/short overview paper).
- W. Huang, F. Zhu, X. Ren, et al., *"All-electron real-time and imaginary-time time-dependent density functional theory within a numeric atom-centered basis function framework,"* **The Journal of Chemical Physics** 155, 154801 (2021).

### GW, BSE, and Excited-State Benchmarks
- Numeric atom-centered orbital *GW*/BSE benchmark and basis-set works, e.g. *"Precision benchmarks for solids: G0W0 calculations with different basis sets"* (arXiv:2411.19701), and BSE-in-FHI-aims implementation papers on molecular excitation benchmarks against Thiel's set.
- Correlation-consistent NAO basis-set development papers, e.g. *"Developing correlation-consistent numeric atom-centered orbital basis sets for Krypton: Applications in RPA-based correlated calculations"* (arXiv:2309.06145).

### Basis Sets and Frozen-Core Approximation
- *"Accurate Frozen Core Approximation for All-Electron Density-Functional Theory,"* arXiv:2106.06412.
- Separable RI-approach benchmarking work, *"Benchmarking the accuracy of the separable resolution of the identity approach for correlated methods in the numeric atom-centered orbitals framework,"* arXiv:2310.11058.

### QM/MM and Solvation
- QM/MM embedding infrastructure with norm-conserving pseudopotentials: **Journal of Chemical Physics** 141, 024105 (2014).
- Implicit solvation formalisms: **Journal of Chemical Theory and Computation** 12, 4052–4066 (2016) and 13, 5582–5603 (2017).

### Accuracy / Reproducibility Benchmarks Involving FHI-aims
- K. Lejaeghere et al., *"Reproducibility in density functional theory calculations of solids,"* **Science** 351, aad3000 (2016) — the Δ-test cross-code benchmark including FHI-aims.
- Additional accuracy demonstrations: **Journal of Physical Chemistry Letters** 8, 1449–1457 (2017); **Physical Review Materials** 1, 033803 (2017).

---

*Note: Several entries above (especially recent arXiv preprints and short conference proceedings) may since have been published in final peer-reviewed form; author lists and exact page/volume numbers for a few secondary items should be cross-checked against the official FHI-aims publication list at fhi-aims.org for citation purposes.*

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the FHI-aims 	All-electron, numeric atom-centered orbital DFT code offering high accuracy across molecules, clusters, and periodic solids. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
