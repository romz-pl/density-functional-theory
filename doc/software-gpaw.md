# GPAW: A Real-Space / Plane-Wave / LCAO DFT Code Implementing the PAW Method

## 1. Overview

GPAW is an open-source, Python-based electronic structure code that implements Kohn–Sham density-functional theory (DFT) using the **projector augmented-wave (PAW) method**. It is developed primarily at the Technical University of Denmark (DTU), the University of Iceland, Aalto University, and a broad international collaborator base, and is distributed under the **GNU General Public License v3.0**. The package is tightly integrated with the **Atomic Simulation Environment (ASE)**, which supplies the Python interface for building atomic structures, driving structure optimization, molecular dynamics, and post-processing — GPAW itself functions as an ASE `Calculator` backend.

A defining architectural feature of GPAW is that it supports **three interchangeable representations of the Kohn–Sham wavefunctions**, all built on the same underlying PAW formalism:

- **Real-space uniform grids** (finite-difference, FD mode) — the code's original and namesake representation
- **Plane waves** (PW mode)
- **Linear combination of atomic orbitals** (LCAO mode) — a numerical, localized atomic-orbital basis

These three modes are mutually independent but connected through transformations via the real-space grid, so a calculation can be initialized cheaply in LCAO and refined in grid or plane-wave mode. This multi-basis flexibility is considered largely unique among mainstream periodic/molecular DFT codes, most of which commit to a single basis type.

## 2. The PAW Method in GPAW

GPAW's PAW implementation follows Blöchl's original 1994 formulation. The core idea:

- A linear transformation operator maps rapidly-oscillating **all-electron** valence wavefunctions near the nucleus onto **smooth pseudo-wavefunctions** that can be represented economically on a numerical grid, in plane waves, or in a localized basis.
- Inside atom-centered augmentation spheres, the pseudo-quantities are re-expanded in terms of atomic partial waves and matched to projector functions; outside the spheres, pseudo and all-electron quantities coincide.
- Unlike norm-conserving or ultrasoft pseudopotentials, PAW retains a **frozen-core, all-electron character**: the full nodal structure of the valence wavefunctions can be reconstructed, giving access to properties (e.g., core-level shifts, hyperfine parameters, X-ray spectra) that are inaccessible to standard pseudopotential methods.
- PAW becomes formally exact in the limit of a complete projector/partial-wave set, and reduces to Vanderbilt ultrasoft pseudopotentials in an appropriate limit — giving smooth, non-normalized pseudo-wavefunctions that are cheap for elements such as first-row atoms (O, N) and transition metals with localized 3d orbitals.
- Atomic PAW **datasets/"setups"** (analogous to pseudopotential files) are generated per element and distributed as an official "GPAW setups" bundle, encoding core densities, partial waves, projectors, and compensation charges.

GPAW's original contribution (Mortensen, Hansen & Jacobsen, 2005) was to implement this PAW formalism directly on **uniform real-space grids** combined with multigrid (Jacobi/Gauss–Seidel-type) iterative solvers and multigrid preconditioning — an alternative to the plane-wave PAW implementations found in codes such as VASP and ABINIT. Real-space grids offer systematic, single-parameter convergence (grid spacing), good parallel scalability via domain decomposition, and natural treatment of non-periodic (cluster/molecule) boundary conditions without artificial periodic images. The later addition of a genuine plane-wave mode and an LCAO mode extended GPAW's applicability to larger systems and different accuracy/cost trade-offs while sharing the same PAW core.

## 3. Core Ground-State DFT Capabilities

- **Exchange–correlation functionals:** LDA, GGA (PBE, revPBE, RPBE, PW91, and others), meta-GGAs, and hybrid functionals; extensive additional functionals available through the **Libxc** library.
- **GLLB-SC** potential for improved band-gap estimates via the derivative discontinuity.
- **DFT+U** (Hubbard-corrected DFT) for localized d/f electron systems.
- Self-consistent field solvers with density mixing, various eigensolvers (RMM-DIIS, conjugate gradient, Davidson), and multigrid preconditioning for the grid mode.
- **Spin-polarized and non-collinear magnetism**, including external magnetic fields, spin–orbit coupling, orbital magnetization, magnetic anisotropy, spin-spiral calculations via the generalized Bloch theorem, and adiabatic magnon dispersions from the magnetic force theorem.
- **Symmetry handling**, k-point sampling, and Berry-phase-based properties: spontaneous polarization, Born effective charges, and piezoelectric response tensors for solids.
- **ΔSCF** (delta self-consistent field) method for approximate excited-state total energies.
- Continuum solvent (implicit solvation) models and **jellium** background-charge boundary conditions.
- Point-group symmetry analysis of vibrational modes.

## 4. Structure, Dynamics, and Optimization (via ASE)

Because GPAW plugs into ASE as a calculator, it inherits ASE's full toolkit for:

- Geometry/structure relaxation (BFGS, FIRE, and other optimizers)
- Nudged elastic band (NEB) calculations for reaction paths and barriers
- Global structure optimization algorithms
- Ab initio molecular dynamics, and QM/MM embedding for hybrid simulations
- Maximally localized Wannier function construction
- Simulated scanning tunneling microscopy (STM) images
- Phonon and electron–phonon coupling calculations

## 5. Excited States and Response Properties

GPAW implements a substantial range of many-body and excited-state methods layered on top of the ground-state PAW machinery:

- **Time-dependent DFT (TDDFT):** both linear-response (Casida-type, for molecules, and via the dielectric-response/density-response function for extended systems) and real-time propagation formulations.
- **Many-body GW approximation** for quasiparticle band structures (G₀W₀ implemented in plane-wave mode).
- **Bethe–Salpeter Equation (BSE)** solver for optical absorption spectra including excitonic effects.
- **Random-phase approximation (RPA)** correlation energies and the linear dielectric response (RPA and ALDA kernels).
- **Variational excited-state methods** via direct orbital optimization (e.g., generalized-mode following/SCF approaches for excited states, including treatment beyond ΔSCF).
- Non-equilibrium Green's function (NEGF) electron transport calculations under finite bias, using the LCAO basis.
- X-ray absorption spectroscopy (XAS).
- Non-linear optical response tensors for solids.
- Magnetic excitation spectra and dynamic magnetic response from TDDFT.
- Calculation of **charged crystal point defects** (formation energies, corrections for periodic charged-defect supercells).

## 6. Parallelization and Performance

- Domain decomposition (real-space grid), plane-wave/G-vector parallelization, band parallelization, k-point parallelization, and spin parallelization, combinable for large-scale MPI runs on HPC clusters.
- Use of **BLAS/ScaLAPACK/BLACS**, **FFTW**, and libxc as core numerical dependencies; core routines are written in C for performance with a Python front end/driver.
- **GPU acceleration** has recently been added, achieved with comparatively modest code modifications by leveraging the **CuPy** library for array operations, enabling offloading of key kernels to GPUs.
- Demonstrated scalability to large processor counts and (per project announcements) multi-million-CPU-hour allocations for GPAW-based research.

## 7. Software Architecture and Development Model

- Written primarily in **Python**, with performance-critical kernels in **C**.
- Hosted on **GitLab** (migrated from earlier hosting) under GPL v3.0, encouraging community contributions, transparency, and reuse.
- Modular design: the separation between PAW core routines and the three wavefunction representations allows new methods to be added to one mode without necessarily reimplementing them for all three; the review paper notes, however, that feature parity across FD, PW, and LCAO modes is incomplete — for example, some newer features (e.g., RPA total energies, stress tensor) currently work only in PW mode, while others (like NEGF transport) require LCAO. The GPAW documentation maintains an up-to-date feature-availability table by mode.
- Installable via `pip install gpaw` (PyPI) with dependencies on ASE, NumPy, SciPy, Libxc, a C compiler, and a BLAS library; MPI, BLACS/ScaLAPACK, and FFTW are optional but strongly recommended for production/parallel work.
- Active release cadence (versioned by year/month, e.g., 24.1, 24.6, 25.1, 25.7, 26.7) with regular developer meetings, mailing lists, and a Matrix chat channel for community support.

## 8. Typical Usage Pattern

A minimal GPAW calculation is expressed as an ASE `Atoms` object combined with a `GPAW` calculator specifying the mode (`PW`, `LCAO`, or default finite-difference grid), exchange-correlation functional, and other parameters — e.g., building an `Atoms` object for a molecule, attaching `GPAW(xc='PBE', mode=PW(300))` as its calculator, and then calling ASE methods such as `get_potential_energy()` or `get_forces()`. This close coupling with ASE's Pythonic API is central to GPAW's usability and to its integration into higher-level materials-discovery and workflow tools.

## 9. Position in the DFT Software Landscape

- GPAW is one of relatively few actively maintained, general-purpose electronic-structure packages offering a **genuine choice between real-space grid, plane-wave, and localized-basis representations within a single, consistent PAW formalism** — most competing PAW/pseudopotential codes (VASP, Quantum ESPRESSO, ABINIT, CP2K, etc.) commit primarily to one representation family.
- Its real-space-grid PAW implementation was among the first of its kind, distinguishing it historically from the plane-wave PAW implementations that dominated the field after Kresse & Joubert's PAW-in-plane-waves formulation.
- The LCAO mode provides a fast, lower-cost alternative for pre-relaxation, screening, and large/molecular-dynamics-scale simulations, with the option to switch to the more accurate grid or plane-wave representation for final, publication-quality results.
- The package's open, Python-native, ASE-integrated design has made it a common backend for high-throughput materials screening pipelines and for coupling DFT with machine-learning interatomic-potential frameworks.

## 10. Notable Limitations / Ongoing Development

- Feature coverage is not uniform across the three wavefunction modes; users should consult GPAW's mode-feature compatibility table before designing calculations relying on newer or more specialized functionality.
- As with plane-wave/grid PAW methods generally, results depend on the quality and transferability of the atomic PAW datasets (setups) used, which must be chosen (and where relevant, validated) for the elements and properties of interest.
- The project's public roadmap and 2024 review outline ongoing work including expanded GPU support, further unification of feature availability across modes, and continued extension of excited-state and magnetic-response capabilities.

---

# Key Publications on GPAW Theory and Implementation

These are the primary references documenting GPAW's theoretical foundations, method development, and the underlying PAW formalism, as recommended by the GPAW developers themselves for citation.

1. **P. E. Blöchl**, "Projector augmented-wave method," *Physical Review B* **50**, 17953 (1994). — The original formulation of the PAW method underlying GPAW.

2. **P. E. Blöchl, C. J. Först, and J. Schimpl**, "Projector augmented wave method: ab initio molecular dynamics with full wave functions," *Bulletin of Materials Science* **26**, 33 (2003). — Extended PAW methodology reference cited in GPAW documentation.

3. **J. J. Mortensen, L. B. Hansen, and K. W. Jacobsen**, "Real-space grid implementation of the projector augmented wave method," *Physical Review B* **71**, 035109 (2005). — The foundational GPAW paper, introducing the real-space, finite-difference/multigrid PAW implementation.

4. **J. Enkovaara, C. Rostgaard, J. J. Mortensen, J. Chen, M. Dułak, L. Ferrighi, J. Gavnholt, C. Glinsvad, V. Haikola, H. A. Hansen, H. H. Kristoffersen, M. Kuisma, A. H. Larsen, L. Lehtovaara, M. Ljungberg, O. Lopez-Acevedo, P. G. Moses, J. Ojanen, T. Olsen, V. Petzold, N. A. Romero, J. Stausholm-Møller, M. Strange, G. A. Tritsaris, M. Vanin, M. Walter, B. Hammer, H. Häkkinen, G. K. H. Madsen, R. M. Nieminen, J. K. Nørskov, M. Puska, T. T. Rantala, J. Schiøtz, K. S. Thygesen, and K. W. Jacobsen**, "Electronic structure calculations with GPAW: a real-space implementation of the projector augmented-wave method," *Journal of Physics: Condensed Matter* **22**, 253202 (2010). — The comprehensive first-generation review paper, covering TDDFT, ΔSCF, XAS, Wannier orbitals, exchange-correlation functionals, and parallelization.

5. **A. H. Larsen, M. Vanin, J. J. Mortensen, K. S. Thygesen, and K. W. Jacobsen**, "Localized atomic basis set in the projector augmented wave method," *Physical Review B* **80**, 195112 (2009). — Introduces the LCAO (localized atomic-orbital basis) mode within the PAW/GPAW framework.

6. **J. J. Mortensen, A. H. Larsen, M. Kuisma, A. V. Ivanov, A. Taghizadeh, A. Peterson, A. Haldar, A. O. Dohn, C. Schäfer, E. Ö. Jónsson, E. D. Hermes, F. A. Nilsson, G. Kastlunger, G. Levi, H. Jónsson, H. Häkkinen, J. Fojt, J. Kangsabanik, J. Sødequist, J. Lehtomäki, J. Heske, J. Enkovaara, K. T. Winther, M. Dulak, M. M. Melander, M. Ovesen, M. Louhivuori, M. Walter, M. Gjerding, O. Lopez-Acevedo, P. Erhart, R. Warmbier, R. Würdemann, S. Kaappa, S. Latini, T. M. Boland, T. Bligaard, T. Skovhus, T. Susi, T. Maxson, T. Rossi, X. Chen, Y. L. A. Schmerwitz, J. Schiøtz, T. Olsen, K. W. Jacobsen, and K. S. Thygesen**, "GPAW: An open Python package for electronic structure calculations," *The Journal of Chemical Physics* **160**, 092503 (2024), doi: 10.1063/5.0182685. — The most recent, comprehensive review, covering the plane-wave and LCAO modes, GW/BSE, magnetism/spin–orbit, non-linear optics, defect calculations, GPU acceleration via CuPy, and the project's future roadmap.

### Supporting / Method-Specific References Commonly Cited Alongside GPAW

- **S. Lehtola, C. Steigemann, M. J. T. Oliveira, and M. A. L. Marques**, "Recent developments in libxc — a comprehensive library of functionals for density functional theory," *SoftwareX* **7**, 1 (2018). — Reference for the Libxc exchange-correlation functional library used by GPAW.
- **A. H. Larsen, J. J. Mortensen, J. Blomqvist, I. E. Castelli, R. Christensen, M. Dułak, J. Friis, M. N. Groves, B. Hammer, C. Hargus, E. D. Hermes, P. C. Jennings, P. B. Jensen, J. Kermode, J. R. Kitchin, E. L. Kolsbjerg, J. Kubal, K. Kaasbjerg, S. Lysgaard, J. B. Maronsson, T. Maxson, T. Olsen, L. Pastewka, A. Peterson, C. Rostgaard, J. Schiøtz, O. Schütt, M. Strange, K. S. Thygesen, T. Vegge, L. Vilhelmsen, M. Walter, Z. Zeng, and K. W. Jacobsen**, "The atomic simulation environment—a Python library for working with atoms," *Journal of Physics: Condensed Matter* **29**, 273002 (2017). — The companion reference for ASE, the framework GPAW is built on and which developers ask users to cite alongside GPAW papers.

*Note: for calculation-specific features (e.g., GW, BSE, RPA correlation, spin-spiral/magnon methods, non-linear optics, DFT+U), the GPAW documentation and the 2024 review paper each cite dozens of additional focused method papers; consult the "How should I cite GPAW?" section of the official documentation (gpaw.readthedocs.io/faq.html) for the full, continuously updated list relevant to a specific feature used.*

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the GPAW 	Real-space/plane-wave/LCAO DFT code implementing the projector augmented-wave (PAW) method, built on the Python ASE framework. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
