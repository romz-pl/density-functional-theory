# Elk: An All-Electron Full-Potential LAPW Code — Exhaustive Review

## 1. Overview

Elk is an all-electron full-potential linearised augmented-plane wave (LAPW) code with many advanced features, written originally at Karl-Franzens-Universität Graz as a milestone of the EXCITING EU Research and Training Network. The code is designed to be as simple as possible so that new developments in electronic structure theory can be added quickly and reliably, and it is freely available under the GNU General Public License (GPLv3). It is presently developed primarily at the Max Planck Institute of Microstructure Physics (Halle) and the Max Born Institute (Berlin), and is described by its authors as a testbed for new methods, while also being extensively used for production work, especially for materials that are particularly sensitive to the type of approximation used or for which pseudopotential methods are inappropriate.

A defining characteristic of the code is that almost all features can be used in combination with each other, and its FAQ notes that Elk is specifically designed to make the implementation of new developments in DFT as easy as possible, with a cleanly written, fully documented, and open-source codebase. As of this review, the latest version is 11.0.2, and the package has been under continuous development for well over 15 years.

---

## 2. Design Philosophy and Usability

- **Single input file** — just one input file is required, with all input parameters optional; the file elk.in should be easier to understand than a GUI, and no GUI is provided.
- **Documentation** — every subroutine and function has a LaTeX header, readable by the Protex Perl script and compiled into a single reference document (elk.pdf) via `make doc`.
- **Units** — Elk exclusively uses Hartree for energy and Bohr for length; Rydbergs, eV, and Ångströms are never used internally.
- **Language and standards** — the code follows strict Fortran 2008 compliance, uses only one programming language, and has a free-form style input file.
- **Community stewardship** — Elk has received both a SourceForge Community Choice award and a Community Leader Award. Notably, no "AI"-generated code is accepted for inclusion in Elk.

---

## 3. Core Methodology: The LAPW Basis

Elk uses the all-electron LAPW approach in which the basis functions consist of a plane wave augmented by a linear combination of atomic-like functions in the muffin-tins, chosen so the resulting function is continuous across the muffin-tin boundary. Efficient evaluation of the Hartree potential is central to the method: considerable effort has gone into a highly optimised implementation of Weinert's method for solving Poisson's equation in the LAPW basis, which is also used with complex charge densities to compute Coulomb matrix elements needed for Hartree-Fock, the optimised effective potential (OEP) method, and reduced density matrix functional theory (RDMFT).

Key basis-related features include:
- LAPW basis with local orbitals, and APW radial derivative matching to arbitrary orders at the muffin-tin surface (super-LAPW, etc.)
- An arbitrary number of local orbitals is allowed, so that all core states can be made valence if desired
- Every element in the periodic table is available, with core states treated using the radial Dirac equation
- Total energies are resolved into individual components

---

## 4. Feature Set

### 4.1 Exchange-Correlation Functionals
LSDA and GGA functionals are available, and variational meta-GGA is available through Libxc, giving access, per the CECAM workshop description, to an interface to the ETSF exchange-correlation library Libxc, making available almost every LDA and GGA functional ever invented.

### 4.2 Structure and Symmetry
Elk determines lattice and crystal symmetry groups from the input lattice and atomic coordinates, can determine atomic coordinates from space-group data via the Spacegroup utility, outputs XCrysDen and V_Sim files, automatically reduces conventional to primitive cells, automatically determines muffin-tin radii, fully symmetrises the density/magnetisation and their conjugate fields, and automatically determines and reduces the k-point set.

### 4.3 Magnetism
This is one of Elk's most distinctive strengths:
- Spin-polarised calculations are performed in the most general way, referring in the code only to (n(r); m(r)) and (vs(r); Bs(r))
- Spin-orbit coupling (SOC) is included in the second-variational scheme; non-collinear magnetism (NCM) is supported with arbitrary on-site magnetic fields
- Fixed spin-moment and fixed tensor moment (experimental) calculations are supported, including with SOC and NCM
- Spin spirals for any q-vector, spin-polarised cores, automatic determination of the magnetic anisotropy energy (MAE), and classical spin/orbital dipole interaction contributions to the Kohn–Sham field are all available

### 4.4 Plotting and Visualization
Band structure with angular momentum character; total/partial DOS with irreducible representation projection; 1/2/3D charge density, exchange-correlation and Coulomb potential, ELF, magnetisation, B_xc and ∇·B_xc plots; 3D Fermi surfaces; wavefunction plots; electric field plots; experimental STM imaging from LDOS; and current density plots are all built in.

### 4.5 Forces, Structural Relaxation, and Phonons
Forces (including incomplete basis set corrections and core corrections) work with spin-orbit coupling, non-collinear magnetism, and LDA+U. Structural capabilities include structural optimisation of atomic positions and lattice vectors, iso-volumetric unit-cell optimisation, and molecular dynamics within the Born–Oppenheimer approximation.

For lattice dynamics: phonons for arbitrary q-vectors are computed with density functional perturbation theory (DFPT), and also via the supercell method; phonon dispersion and DOS, non-analytic LO-TO splitting for polar semiconductors (via the dielectric tensor and Born effective charges), and thermodynamic quantities (free energy, entropy, heat capacity) from the phonon DOS are all available. Superconductivity-related quantities include electron-phonon coupling matrices, phonon linewidths, the Eliashberg function α²F(ω), the electron-phonon coupling constant λ, the McMillan-Allen-Dynes critical temperature Tc, and self-consistent solution of the Eliashberg equations. Phonon calculations can also be distributed across networked computers.

### 4.6 Advanced and Beyond-DFT Methods
This is where Elk most distinguishes itself among LAPW packages:

| Method | Capability |
|---|---|
| Exact exchange (EXX/OEP) | Exact exchange optimised effective potential method and EXX energies, both compatible with SOC and NCM |
| Hartree-Fock | Hartree-Fock for solids, with SOC and NCM support |
| DFT+U | Fully localised limit (FLL), around mean field (AFM), and interpolation between the two; compatible with SOC, NCM, and spin-spirals. According to the developers, Elk has the most sophisticated implementation of DFT+U available, usable together with spin-orbit coupling, non-collinear magnetism and spin-spirals, with the ability to interpolate between around-mean-field (AMF) and the fully localised limit (FLL) |
| RDMFT | Reduced density matrix functional theory for solids — Elk is the only solid-state code capable of calculating ground-state properties and spectra using one-body RDMFT |
| BSE | Bethe-Salpeter equation, including beyond the Tamm-Dancoff approximation, compatible with SOC and NCM |
| TDDFT | Linear-response TDDFT for optical response, femtosecond-scale real-time evolution of the electronic state, and Ehrenfest dynamics via TDDFT molecular dynamics. Elk is the only code capable of real-time propagation for extended systems, and the only solid-state code capable of real-time dynamics of spins and charge under strong electromagnetic fields |
| GW | GW approximation spectral functions and GW spectral function band structures, both experimental and compatible with SOC/NCM |
| Ultra long-range (ULR) | An experimental method for calculations using a generalisation of the Bloch ansatz, enabling solid-state DFT calculations on systems of almost unlimited size, up to the micron length scale, while retaining microscopic electronic-structure detail, based on an additional sum over a finer reciprocal-space grid around each k-point |
| Born effective charges | Static and experimental dynamical Born effective charges |

### 4.7 Miscellaneous Specialized Outputs
Mössbauer hyperfine parameters (isomer shift, EFG, hyperfine contact fields); first-order optical response; Kerr angle/MOKE output (experimental); ELNES; non-linear optical second-harmonic generation; L, S and J expectation values; effective mass tensors; equation-of-state fitting; iterative diagonalisation with fine-grained parallelisation round out the feature list.

### 4.8 Parallelism and Programming
Elk supports OpenMP parallelisation, MPI parallelisation, and efficient hybrid OpenMP+MPI parallelism.

---

## 5. Recent Developments (Selected Highlights from Release History)

- **v11.0.2** — further improved Wannier90 interface, allowing individual angular-momentum values per species; added tesseral tensor moments
- **v10.9.5** — added calculated Wannier90 s/p/d/f projectors and subspace diagonalisation as an experimental alternative to iterative diagonalisation
- **v9.1.15** — added piezoelectric and magnetoelectric tensors (experimental); enabled exchange-correlation functionals with Laplacian terms, including the deorbitalised functionals of Mejía-Rodríguez and Trickey
- **v8.6.7** — added dynamical Born effective charges, published in Phys. Rev. B 106, L180303 (2022)
- **v7.0.12** — added electron-phonon mean-field theory (highly experimental), and an interface to the DMFT code TRIQS maintained as a separate branch; also the ultra long-range method was published in Phys. Rev. Lett. 125, 256402 (2020)
- **v6.2.8** — Wannier90 interface added, usable for Hartree-Fock band structures and with non-collinear spin-polarised calculations; self-consistent density GW became available (experimental)

Elk interoperates with several external tools, including PyProcar for pre/post-processing, the Elk Optics Analyzer for optics output, DGrid, ASE, Libxc, NOMAD, and Phonopy.

---

## 6. Licensing, Governance, and Special Versions

- **License**: GNU General Public License (v3 or later)
- **Main authors/contributors**: a large list led by Kay Dewhurst, Sangeeta Sharma, Lars Nordström, Francesco Cricchio, Oscar Grånäs, Hardy Gross, and dozens of others
- **Deriving codes**: contributors are encouraged to fork and specialise, provided they use a different name for their version — e.g. "Elk-Uppsala". Known specialised branches include Elk-TRIQS (DMFT interface) and elk-w90-improved (full Wannier90 interface)
- **Sister code**: the exciting code, which shares Elk's origins in the EXCITING EU Research and Training Network

---

## 7. Distribution and Deployment

Elk is packaged and available on numerous HPC systems (e.g., modules maintained at DIPC and PDC/KTH) and is downloadable via SourceForge, with the standard build using its own bundled or system BLAS/LAPACK/FFT libraries. Deployment commonly uses hybrid MPI+OpenMP execution for best performance.

---

## 8. Citing Elk

Formal citation is not mandatory but is appreciated; a reference to the Elk website via the provided BibTeX entry will suffice for general use of the code.

---

## 9. Publications Related to Elk's Theory and Methods

The following list gathers the principal theoretical/methodological publications underlying features implemented in Elk, based on citations from the developers and related literature:

### Overview / General Code
- J. K. Dewhurst, S. Sharma, *Development of the Elk LAPW Code* (MPI Halle overview article/talk)
- Weinert, M., *Solution of Poisson's equation: beyond Ewald-type methods*, J. Math. Phys. **22**, 2433 (1981) — underlying method for Elk's Poisson-equation solver

### LDA+U
- Around-mean-field / fully-localised-limit interpolation scheme: Phys. Rev. B **67**, 153106 (2003)

### Reduced Density Matrix Functional Theory (RDMFT)
- S. Sharma, J. K. Dewhurst, N. N. Lathiotakis, E. K. U. Gross, *Reduced Density Matrix Functional for Many-Electron Systems*, Phys. Rev. B **78**, 201103(R) (2008) [arXiv:0801.3787]
- S. Sharma, J. K. Dewhurst, S. Shallcross, E. K. U. Gross, *Spectral Density and Metal-Insulator Phase Transition in Mott Insulators within Reduced Density Matrix Functional Theory*, Phys. Rev. Lett. (2013)
- N. N. Lathiotakis, N. Helbig, E. K. U. Gross, Phys. Rev. A **72**, 030501 (2005)
- N. N. Lathiotakis, N. Helbig, E. K. U. Gross, Phys. Rev. B **75**, 195120 (2007)

### Ultra Long-Range (ULR) Method
- T. Müller, S. Sharma, E. K. U. Gross, J. K. Dewhurst, *Extending Solid-State Calculations to Ultra Long-Range Length Scales*, Phys. Rev. Lett. **125**, 256402 (2020) [arXiv:2008.12573]

### Coupled Electron-Phonon Bogoliubov Method
- Editors' Suggestion, Phys. Rev. B **105**, 174509 (2022)

### Dynamical Born Effective Charges
- Phys. Rev. B **106**, L180303 (2022)

### Meta-GGA / Deorbitalisation
- Bartók, A. P., Yates, J. R., *Regularized SCAN functional*, J. Chem. Phys. **150**, 161101 (2019)
- Mejía-Rodríguez, D., Trickey, S. B., *Deorbitalised meta-GGA exchange-correlation functionals with Laplacian terms*, Phys. Rev. B **98**, 115161 (2018)
- Meta-GGA partial deorbitalisation extension covering the full Libxc kinetic-energy functional range: arXiv:2304.02363

### Born Effective Charges (King–Smith–Vanderbilt Method)
- King-Smith, R. D., Vanderbilt, D., *Theory of Polarization of Crystalline Solids*, Phys. Rev. B **47**, 1651(R) (1993)

### Non-linear Optics (formalism correction)
- Original formalism: Phys. Rev. B **67**, 165332 (2003) (subsequently corrected and re-implemented in Elk after issues identified by Y. Liang and X. Gonze)

### Positron / Two-Component DFT (used alongside Elk in Fermi-surface studies)
- Two-component DFT for electron-positron systems, cited in ELK-based 2D-ACAR studies (arXiv:2104.12563)

---

**Note:** Full up-to-date reference lists, the complete BibTeX citation entry, and the detailed LaTeX-compiled subroutine documentation (`elk.pdf`) are maintained directly on the [Elk homepage](https://elk.sourceforge.io/) and should be consulted for the authoritative and continuously updated bibliography.

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the Elk 	Open-source, full-potential all-electron LAPW code, notable for extensive support of advanced DFT and beyond-DFT methods. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
