# SALMON: An Exhaustive Review of the Open-Source Real-Time TDDFT Code for Light–Matter Interaction Simulations

## 1. Overview

**SALMON** (Scalable Ab-initio Light–Matter simulator for Optics and Nanoscience) is an open-source software package for first-principles simulation of electron dynamics and optical properties in molecules, nanostructures, and crystalline solids. It is built on **real-time, real-space time-dependent density functional theory (TDDFT)**, solving the time-dependent Kohn–Sham (TDKS) equation directly in the time domain rather than through frequency-domain linear-response formalisms based on matrix diagonalization.

- **Homepage:** http://salmon-tddft.jp
- **Source code:** https://github.com/SALMON-TDDFT/SALMON2 (development version, "SALMON2")
- **License:** Apache License, Version 2.0
- **Primary language:** Fortran 2003, with hybrid MPI + OpenMP (and OpenACC/GPU support in recent releases)
- **Build system:** CMake (≥3.14 for SALMON2)

SALMON emerged from the **merger of two predecessor codes**:
- **GCEED** (Grid-based Coupled Electron and Electromagnetic field Dynamics), developed at the Institute for Molecular Science, targeting molecules and nanostructures with large-scale spatial/orbital parallelization.
- **ARTED** (Ab-initio Real-Time Electron Dynamics simulator), developed mainly at the University of Tsukuba, targeting crystalline solids and featuring coupled Maxwell–TDDFT multiscale propagation of light pulses through bulk media.

By unifying these, SALMON provides a single, consistent framework applicable to **isolated systems** (molecules, clusters, nanoparticles) and **periodic systems** (crystalline solids) alike, including the coupled propagation of electromagnetic fields through matter.

---

## 2. Scientific and Numerical Foundations

### 2.1 Governing equation

SALMON propagates the time-dependent Kohn–Sham equation

$$
i\hbar \frac{\partial}{\partial t} \psi_p(\mathbf{r}, t) = \left[ -\frac{\hbar^2}{2m}\nabla^2 + V_{\text{ion}} + e^2\!\int \frac{\rho(\mathbf{r}',t)}{|\mathbf{r}-\mathbf{r}'|}\,d\mathbf{r}' + V_{\text{xc}} + V_{\text{ext}} \right]\psi_p(\mathbf{r},t)
$$

directly on a **uniform real-space Cartesian grid**, rather than in a plane-wave or localized-basis representation. This real-space, real-time approach was pioneered by Yabana and Bertsch (1996) and later extended to crystalline solids by Bertsch, Iwata, Rubio, and Yabana (2000).

Key numerical/methodological ingredients:

- **Pseudopotentials:** Norm-conserving pseudopotentials (Troullier–Martins type) with the Kleinman–Bylander separable form for the nonlocal part.
- **Exchange–correlation treatment:** Adiabatic approximation — the ground-state XC functional form is reused at each time step, evaluated on the instantaneous time-dependent density. Supported functionals include LDA (Perdew–Zunger), the Becke–Johnson (BJ) potential and its Tran–Blaha modification (useful for reproducing insulator band gaps), and meta-GGA/hybrid functionals in later extensions.
- **Boundary conditions:** Vacuum (isolated systems) or fully periodic (crystalline solids); a cuboid-shaped simulation cell is used, so non-cuboid primitive cells require an appropriate supercell.
- **Gauge choice for periodic systems:** Velocity gauge, preserving lattice periodicity of the Hamiltonian at all times; Bloch orbitals $w_{n\mathbf{k}}(\mathbf{r},t)$ are evolved with $\psi_p(\mathbf{r},t) = e^{i\mathbf{k}\cdot\mathbf{r}}w_{n\mathbf{k}}(\mathbf{r},t)$.
- **Time propagation:** Explicit integration via truncated Taylor expansion of the time-evolution operator, optionally combined with a predictor–corrector scheme.
- **External field treatment:**
  - **Impulsive (delta-function) kick** — for linear-response quantities (polarizability, photoabsorption spectrum, dielectric function) extracted from a single time-propagation run via Fourier transform of the induced dipole/current.
  - **Finite pulsed fields** (single- or two-pulse, arbitrary CEP, or user-supplied time-domain tables) — for simulating nonlinear, non-perturbative response to intense ultrashort laser pulses, including pump–probe-type setups.
- **Ground-state initialization:** A static Kohn–Sham (DFT) SCF calculation precedes the real-time propagation to generate the initial orbitals, densities, band structures/orbital energies, density of states, and interatomic forces.

### 2.2 Multiscale Maxwell + TDDFT coupling

A distinguishing capability of SALMON (inherited from ARTED) is the **simultaneous, coupled solution of Maxwell's equations for the macroscopic electromagnetic field and the microscopic TDKS equation** in each unit cell of a spatial grid representing the bulk medium. The 1D wave equation for the vector potential,

$$
\frac{1}{c^2}\frac{\partial^2 A}{\partial t^2} - \frac{\partial^2 A}{\partial X^2} = \frac{4\pi}{c} I(X,t),
$$

is solved on a coarse macroscopic grid $X$, while at every macroscopic grid point an independent microscopic TDKS calculation supplies the local induced current $I(X,t)$ that feeds back into the field propagation. This enables first-principles simulation of:

- Laser-pulse propagation, reflection, and transmission at a bulk material surface
- Nonlinear/irreversible energy transfer from an intense pulse into the electronic system as a function of penetration depth
- High-harmonic generation during propagation through solids
- Coherent phonon generation, dielectric breakdown, and ultrafast (attosecond-scale) changes in dielectric response

At present, multiscale propagation is supported along a single spatial direction (e.g., a linearly polarized pulse normally incident on a flat bulk surface), since full 2D/3D Maxwell–TDDFT coupling is computationally far more demanding.

### 2.3 Observables and outputs

- **Ground state:** Kohn–Sham orbitals and energies/bands, density of states and projected DOS, electron localization function (ELF), forces
- **Linear response:** Oscillator-strength distribution / photoabsorption spectrum (molecules, nanostructures), frequency-dependent dielectric function and optical conductivity (solids)
- **Nonlinear real-time dynamics:** Time-dependent electron density and current, induced dipole/polarization, electronic excitation energy, number density of excited electron–hole pairs, ionic forces (enabling coupled electron–ion/Ehrenfest-type dynamics and structural optimization/MD)
- **Field propagation:** Spatiotemporal electromagnetic field profiles, Drude–Lorentz model outputs, metasurface optical response, energy transfer from light to electrons as a function of depth

---

## 3. Software Architecture and Parallelization

- SALMON exposes its computational modes through Fortran namelists (`&calculation`, `&system`, `&rgrid`, `&kgrid`, `&tgrid`, `&propagation`, `&emfield`, `&multiscale`, `&analysis`, `&hartree`, `&opt`, `&md`, etc.), providing a single unified input-file syntax across ground-state, linear-response, real-time, and multiscale calculation modes.
- **Parallelization strategy:** Hybrid MPI + OpenMP, with domain decomposition along two axes — the **spatial grid** (division of the simulation cuboid) and the **orbital index** (and, for periodic systems, **k-points**). Spatial-domain parallelization is favored during the SCF ground-state stage (due to sequential Gram–Schmidt orthogonalization), while orbital-based parallelization dominates during real-time propagation, since each orbital can be evolved largely independently. An automatic process-assignment heuristic is available for process counts factorizable into 2, 3, and 5.
- **Scalability:** The code has been demonstrated on large parallel machines — e.g., the K computer (SPARC64 VIIIfx), Oakforest-PACS (Intel Xeon Phi Knights Landing), and Fujitsu FX100/A64FX systems — with electron-dynamics calculations reported for systems of up to a few thousand atoms and near-linear parallel efficiency (>99% in some benchmark core-shell nanoparticle calculations up to ~2,000 MPI processes).
- **GPU acceleration:** Recent SALMON2 releases (v2.2.x) add OpenACC-based GPU acceleration for DFT/TDDFT computation (validated with the NVIDIA HPC SDK), along with preconditioned conjugate-gradient SCF solvers and general performance improvements.
- **I/O and restart:** Supports restart from checkpoint data; pseudopotential file formats include `.cpi` (fhi98pp-style) and `.fhi` (ABINIT-compatible) formats.

---

## 4. Version History and Development Status

| Version | Release | Notable changes |
|---|---|---|
| v1.0–v1.2.x | 2017–2019 | Initial public releases merging GCEED and ARTED functionality |
| v1.2.2 | 2019-05-21 | Stable release of the "version 1" code line |
| **Code paper published** | 2019-02-01 | SALMON code paper published in *Computer Physics Communications* |
| v2.0.0 | 2020-07-22 | Major rewrite ("SALMON2"): common spatial grid for Maxwell–TDDFT simultaneous calculation, two-temperature-model coupling with FDTD, revised parallelization scheme, revised restart-file format, new/changed input keywords (some v1 functionality, e.g. orbital projection, temporarily unavailable) |
| v2.0.1 | 2021-01-29 | Bug fixes |
| v2.0.2 | 2021-11-06 | Bug fixes (total-energy noise in weak-pulse TDDFT, cube-file density output, external-field printing, restart/projection fixes) |
| v2.1.0 | 2022-03-31 | Feature and stability updates |
| v2.2.0 | 2023-04-21 | Feature and stability updates |
| v2.2.1 | 2024-05-31 | OpenACC GPU acceleration (NVIDIA HPC SDK v23.11 verified), CG preconditioning |
| v2.2.2 | current | Further GPU/OpenACC fixes, corrected pseudopotential/LAPACK links, index-checking and compiler-compatibility fixes |

Development is ongoing on GitHub under the `SALMON-TDDFT` organization, which hosts:
- **SALMON2** — the actively developed/maintained code
- **SALMON** — the legacy version-1 code line
- **SALMON-DOCS** — manual/documentation source
- **SALMON-inputs** — a database of input files reproducing published papers
- **SALMON-VR** and **SALMON2-evaluation-scripts** — auxiliary repositories

The project maintains a active user community, periodic international schools/workshops (e.g., "Ab Initio Electron Dynamics Simulations," most recently scheduled for March 2026 in Nara, Japan), and video tutorials.

---

## 5. Comparison with Related Codes

SALMON belongs to the family of **real-space, real-time TDDFT** codes and shares numerical heritage with several related packages:

- **Octopus** — shares a common methodological root with SALMON (a collaborative real-space/real-time TDDFT effort circa 2000); widely used, broad feature set for molecules and solids.
- **ARTED / GCEED** — SALMON's direct predecessors and now superseded by it.
- **Qbox, FPSID** — plane-wave-basis real-time TDDFT implementations, some with coupled ion dynamics.
- **SIESTA-based and Elk FP-LAPW implementations** — real-time TDDFT using numerical atomic orbitals or full-potential LAPW bases, respectively.
- **NWChem** — real-time TDDFT capability within a standard quantum-chemistry package.

SALMON differentiates itself primarily through its **native, tightly coupled multiscale Maxwell–TDDFT propagation** for simulating pulsed-light propagation through bulk/thin-film media — a capability inherited from ARTED and not commonly found as a built-in feature elsewhere — combined with strong HPC scalability for large isolated nanostructures inherited from GCEED.

---

## 6. Typical Application Domains

Based on the developers' own example calculations and the broader literature using SALMON:

- **Molecular photoabsorption spectroscopy** (e.g., small molecules such as acetylene) and few-cycle laser-pulse-driven electronic excitation dynamics
- **Metal/semiconductor core–shell nanostructures** (e.g., Ag@Si nanoparticles) — size- and composition-dependent optical response
- **Dielectric function calculations for bulk semiconductors** (e.g., silicon), including comparison of LDA vs. Becke–Johnson potentials against experimental optical data
- **Strong-field / high-intensity laser–solid interaction physics**, including:
  - Nonlinear electron-hole pair generation and multiphoton absorption (e.g., in TiS₂ under intense few-cycle pulses)
  - Coherent phonon generation
  - Ultrafast/attosecond dynamical Franz–Keldysh effect
  - Dielectric breakdown of insulators under intense pulses
  - High-harmonic generation in bulk solids and during pulse propagation through dielectrics
  - Attosecond band-gap dynamics and light-matter energy transfer (relevant to attosecond metrology experiments)
- **Metasurface and near-field nano-optics**, where non-dipole, spatially nonuniform optical near-fields drive otherwise-forbidden transitions and second-harmonic generation

---

## 7. Requirements, Availability, and Licensing

- **Compilers/libraries:** MPI-capable C and Fortran (2003+) compilers, a LAPACK implementation (Intel MKL or Fujitsu SSL-II have been most thoroughly tested)
- **Build:** CMake-based (`configure.py --arch=<ARCH> --prefix=<DIR>`, then `make && make install`); an `--disable-mpi` flag permits single-process builds
- **Execution:** `salmon.cpu < inputfile.inp > fileout.out` (serial/single node) or `mpiexec -n NPROC salmon.cpu < inputfile.inp > fileout.out` (parallel); a `salmon.mic` binary variant exists for Xeon Phi many-core systems
- **License:** Apache License 2.0 — permissive, allows free reuse, modification, and redistribution, including in derivative/commercial works, provided attribution and license terms are preserved
- **Cost:** Free and fully open source

---

## 8. Strengths and Limitations

**Strengths**
- Unified treatment of isolated and periodic systems within one code base and input syntax
- Full nonlinear, non-perturbative real-time propagation — not restricted to the linear-response regime
- Native coupled Maxwell–TDDFT multiscale capability for realistic pulse-propagation simulations, largely unique among general-purpose TDDFT packages
- Demonstrated strong parallel scalability to large systems (thousands of atoms) on leadership-class HPC systems, with increasing GPU support
- Permissive open-source licensing and an actively maintained public input-file database tied to published research

**Limitations**
- Currently restricted to **spin-unpolarized (spin-saturated)** systems and **nonmagnetic** materials
- Cuboid-only simulation cell — non-cuboid crystal primitive cells require supercells, increasing computational cost
- k-point sampling does not yet exploit crystal symmetry to reduce cost
- Multiscale Maxwell–TDDFT coupling is currently limited to 1D (normal-incidence, planar) light propagation geometries
- Real-space time-dependent exchange-correlation is treated adiabatically; current/vector-potential-dependent XC functionals (e.g., Vignale–Kohn) proposed for periodic systems are not yet implemented
- As with all real-time propagation approaches, obtaining high-resolution spectra (especially at low energy/long times) can require long propagation times and correspondingly higher computational cost compared to some frequency-domain linear-response methods

---

## 9. How to Cite SALMON

The developers request that the core code paper always be cited when SALMON contributes meaningfully to a publication, with additional citations depending on the specific functionality used (see the theory publication list below, items 1–7).

---

## 10. Publications Related to SALMON's Theoretical Foundations

The following is the developer-designated list of foundational theory and methodology papers underlying SALMON, as maintained on the official "About SALMON" citation page, supplemented with the primary code-description paper.

1. **Code paper (cite for any significant use of SALMON):**
   M. Noda, S. A. Sato, Y. Hirokawa, M. Uemoto, T. Takeuchi, S. Yamada, A. Yamada, Y. Shinohara, M. Yamaguchi, K. Iida, I. Floss, T. Otobe, K.-M. Lee, K. Ishimura, T. Boku, G. F. Bertsch, K. Nobusada, K. Yabana.
   *SALMON: Scalable ab-initio light-matter simulator for optics and nanoscience.*
   **Computer Physics Communications**, 235, 356–365 (2019).

2. **Massively parallel implementation (large isolated systems):**
   M. Noda, K. Ishimura, K. Nobusada, K. Yabana, T. Boku.
   *Massively-parallel electron dynamics calculations in real-time and real-space: toward applications to nanostructures of more than ten-nanometers in size.*
   **Journal of Computational Physics**, 265, 145–155 (2014).

3. **Formalism for periodic systems / dielectric function:**
   G. F. Bertsch, J.-I. Iwata, A. Rubio, K. Yabana.
   *Real-space, real-time method for the dielectric function.*
   **Physical Review B**, 62(12), 7998 (2000).

4. **Original real-time TDDFT linear-response (impulsive-kick) implementation:**
   K. Yabana, G. F. Bertsch.
   *Time-dependent local-density approximation in real time.*
   **Physical Review B**, 54, 4484–4487 (1996).

5. **Multiscale Maxwell–TDDFT formalism for strong fields in solids:**
   K. Yabana, T. Sugiyama, Y. Shinohara, T. Otobe, G. F. Bertsch.
   *Time-dependent density functional theory for strong electromagnetic fields in crystalline solids.*
   **Physical Review B**, 85(4), 045134 (2012).

6. **Parallelization of coupled Maxwell–TDDFT calculations:**
   S. A. Sato, K. Yabana.
   *Maxwell + TDDFT multi-scale simulation for laser-matter interactions.*
   **Journal of Advanced Simulation in Science and Engineering**, 1(1), 98–110 (2014). doi:10.15748/jasse.1.98

7. **Many-core-processor computational aspects (periodic systems):**
   Y. Hirokawa, T. Boku, S. A. Sato, K. Yabana.
   *Electron dynamics simulation with time-dependent density functional theory on large scale symmetric mode Xeon Phi cluster.*
   In: **Parallel and Distributed Processing Symposium Workshops, 2016 IEEE International**, pp. 1202–1211 (2016).

### Additional foundational and closely related methodology references (cited within the SALMON code paper)

- E. Runge, E. K. U. Gross. *Density-functional theory for time-dependent systems.* **Physical Review Letters**, 52(12), 997 (1984). — Foundational TDDFT theorem underlying the entire approach.
- G. F. Bertsch, J.-I. Iwata, A. Rubio, K. Yabana. *(See item 3 above; foundational periodic real-space/real-time formalism.)*
- K. Yabana, G. Bertsch. *Time-dependent local-density approximation in real time: Application to conjugated molecules.* **International Journal of Quantum Chemistry**, 75(1), 55–66 (1999).
- K. Yabana, T. Nakatsukasa, J. Iwata, G. Bertsch. *Real-time, real-space implementation of the linear response time-dependent density-functional theory.* **Physica Status Solidi B**, 243(5), 1121–1138 (2006).
- N. Troullier, J. L. Martins. *Efficient pseudopotentials for plane-wave calculations.* **Physical Review B**, 43(3), 1993 (1991). — Norm-conserving pseudopotential construction used by SALMON.
- L. Kleinman, D. Bylander. *Efficacious form for model pseudopotentials.* **Physical Review Letters**, 48(20), 1425 (1982). — Separable nonlocal pseudopotential form used in SALMON.
- J. P. Perdew, A. Zunger. *Self-interaction correction to density-functional approximations for many-electron systems.* **Physical Review B**, 23(10), 5048 (1981). — LDA functional used as SALMON's default XC treatment.
- A. D. Becke, E. R. Johnson. *A simple effective potential for exchange.* **Journal of Chemical Physics**, 124(22), 221101 (2006). — Becke–Johnson potential option in SALMON.
- F. Tran, P. Blaha. *Accurate band gaps of semiconductors and insulators with a semilocal exchange-correlation potential.* **Physical Review Letters**, 102, 226401 (2009). — TB-mBJ potential option in SALMON.
- M. Noda, K. Ishimura, K. Nobusada. *Program package of photoinduced electron dynamics: GCEED (grid-based coupled electron and electromagnetic field dynamics).* Proceedings of Computational Science Workshop 2014 (CSW2014), p. 011010 (2015). — Predecessor code to SALMON.
- Y. Hirokawa, T. Boku, S. A. Sato, K. Yabana. *Performance evaluation of large scale electron dynamics simulation under many-core cluster based on Knights Landing.* Proceedings of HPC Asia 2018, pp. 183–191 (2018).
- G. Vignale, W. Kohn. *Current-dependent exchange-correlation potential for dynamical linear response theory.* **Physical Review Letters**, 77, 2037–2040 (1996). — Discusses the vector-potential XC extension noted as not yet implemented in SALMON.

---

## 11. Key References for Further Reading

- SALMON official website: http://salmon-tddft.jp
- Online manual: https://salmon-tddft.jp/webmanual/current/html/index.html
- GitHub organization: https://github.com/salmon-tddft
- Noda, M. et al. *SALMON: Scalable Ab-initio Light-Matter simulator for Optics and Nanoscience.* **Computer Physics Communications** 235, 356–365 (2019); preprint arXiv:1804.01404.
- List of "Publications using SALMON" (application papers by the broader user community): https://salmon-tddft.jp/publications_using_SALMON.html

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the SALMON 	Open-source real-time TDDFT code for light-matter interaction simulations in molecules and solids. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
