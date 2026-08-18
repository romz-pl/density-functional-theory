# ONETEP: A Comprehensive Technical Review

**Order-N Electronic Total Energy Package — Linear-Scaling, Plane-Wave-Accuracy Density Functional Theory**

---

## 1. Overview

ONETEP (Order-N Electronic Total Energy Package) is a linear-scaling density functional theory (DFT) code designed to perform large-basis, plane-wave-accuracy quantum-mechanical calculations on systems containing thousands to tens of thousands of atoms, on parallel high-performance computers. It was conceived to close the gap between the accuracy of conventional plane-wave pseudopotential DFT codes (such as CASTEP or VASP) — which scale as $\mathcal{O}(N^3)$ with the number of atoms $N$ — and linear-scaling ($\mathcal{O}(N)$) methods, which historically sacrificed basis-set completeness and transferability to achieve favourable scaling.

ONETEP achieves this by reformulating the plane-wave pseudopotential method in terms of the single-particle density matrix, expressed using strictly localized, in-situ-optimized local orbitals — Non-orthogonal Generalized Wannier Functions (NGWFs) — themselves expanded in a systematic, plane-wave-equivalent basis of periodic sinc (psinc) functions. The result is a code that recovers full plane-wave accuracy and variational convergence with respect to basis-set size, while the computational cost of the calculation scales linearly with system size once localization constraints are imposed.

ONETEP has been under continuous development since the early 2000s by a consortium of UK academic groups (originally Cambridge, Southampton, and Imperial College London; later joined by Warwick and, more recently, Gdańsk University of Technology), and is distributed under an academic licence via Cambridge Enterprise, with a commercial licence route available (historically marketed through Accelrys/BIOVIA Materials Studio).

---

## 2. Theoretical and Methodological Foundations

### 2.1 The linear-scaling density-matrix reformulation

Conventional Kohn–Sham DFT constructs the electron density from $\mathcal{O}(N)$ orthonormal eigenstates $\psi_i(\mathbf{r})$, each expanded in an $\mathcal{O}(N)$-sized basis set, with an enforced mutual orthogonality that costs $\mathcal{O}(N^3)$ operations. ONETEP instead works directly with the one-particle density matrix,

$$\rho(\mathbf{r}, \mathbf{r}') = \sum_i \psi_i(\mathbf{r}) f_i \psi_i^{*}(\mathbf{r}') = \phi_\alpha(\mathbf{r})\, K^{\alpha\beta}\, \phi_\beta^{*}(\mathbf{r}')$$

where $\phi_\alpha(\mathbf{r})$ are the NGWFs (Greek indices summed implicitly), and $K^{\alpha\beta}$ is the **density kernel** — a generalized, non-orthogonal analogue of the occupation-number matrix that simultaneously encodes orbital occupancy and the mutual non-orthogonality of the underlying local orbitals.

Because electronic structure exhibits "nearsightedness" in systems with a non-zero band gap (the Kohn/Prodan–Kohn principle), both the NGWFs and the density kernel can be truncated: NGWFs are strictly confined within spherical localization regions of radius $R_\phi$ centred on atoms, and the density kernel elements $K^{\alpha\beta}$ are set to zero beyond a chosen cutoff distance between orbital centres. This yields sparse matrix representations of the Hamiltonian, overlap, and density-kernel matrices, and sparse linear algebra can then be used throughout, recovering strictly linear-scaling computational effort with system size — asymptotically, once the kernel cutoff is finite.

### 2.2 NGWFs and the psinc basis set

The NGWFs are not fixed, pre-tabulated atomic orbitals (as in typical linear-scaling codes using numerical atomic orbitals or Gaussians); instead they are **variationally optimized in situ** during the calculation, adapting to the specific chemical environment of each atom. This "self-consistent, optimized minimal basis" strategy allows ONETEP to use a genuinely minimal number of orbitals per atom while retaining accuracy equivalent to a large, converged plane-wave basis — a combination that is difficult to achieve with fixed local-orbital sets.

To allow this in-situ optimization while retaining full basis-set systematic improvability, each NGWF is expanded in periodic cardinal sine (**psinc**) functions, $D(\mathbf{r}-\mathbf{r}_m)$, centred on the points of a regular real-space Cartesian grid. Psinc functions are the real-space, band-width-limited analogue of a plane-wave basis truncated at a kinetic-energy cutoff: they are mathematically equivalent, via Fourier transform, to a plane-wave basis set restricted to reciprocal-space vectors below the cutoff. This equivalence means that:

- The psinc/NGWF basis inherits the systematic improvability of plane waves (a single cutoff-energy parameter controls convergence, exactly as in conventional plane-wave codes).
- The basis is inherently orthogonal at the level of the underlying psinc grid, unlike Gaussian or Slater-type bases, avoiding basis-set superposition error (BSSE) in the underlying representation.
- Kinetic energies and short-range physics are captured with the same fidelity as full plane-wave calculations, in contrast to atom-centred bases which can struggle to represent regions far from any atom (e.g. transition states, or interstitial/vacuum regions).

An **FFT-box** technique is used to keep this scheme's computational cost independent of overall system size: each NGWF is periodic only within a small box (typically of the order of 20×20×20 Å) surrounding its localization region, rather than the full simulation cell, so that FFT operations (used for kinetic-energy evaluation, interpolation to avoid aliasing, and Hamiltonian application) do not grow with system size.

### 2.3 Nested-loop energy minimization

A single-point ONETEP calculation proceeds via **nested-loop conjugate-gradients minimization** of the total-energy functional with respect to two independent sets of variables:

1. **Inner loop** — optimization of the density kernel $K^{\alpha\beta}$ for fixed NGWFs, subject to the constraints that the kernel be idempotent (ensuring the NGWFs represent physically valid, mutually orthogonalizable orbitals) and normalized (trace equal to the correct electron count). This is most commonly performed via the **LNV (Li–Nunes–Vanderbilt) purification method**, a polynomial expansion technique that enforces idempotency implicitly through an auxiliary, non-idempotent kernel and iterative purification, avoiding costly explicit orthogonalization. Alternative or complementary approaches implemented include density/kernel mixing (Kernel-DIIS, restricted to insulating systems and formally cubic-scaling) and, for metals, ensemble/finite-temperature methods (see below).
2. **Outer loop** — optimization of the NGWF expansion coefficients (i.e., the shapes of the local orbitals themselves) for a fixed density kernel, also via conjugate gradients.

Because both the kernel and the NGWFs are variationally optimized, the total energy is variational with respect to *both* basis-set completeness (psinc cutoff) and localization-region size, giving ONETEP well-defined, systematically improvable convergence parameters analogous to conventional plane-wave DFT, in contrast to many other local-orbital linear-scaling schemes.

### 2.4 Metallic and finite-temperature systems: EDFT and AQuA-FOE

The LNV kernel-optimization scheme intrinsically assumes a non-zero HOMO–LUMO gap (integer occupancies), which is inappropriate for metallic or near-metallic systems. ONETEP therefore also implements **Ensemble DFT (EDFT)**, in which Kohn–Sham states are allowed fractional occupancies determined self-consistently via a Fermi–Dirac distribution at a finite electronic temperature, and the Helmholtz free energy (rather than the total energy) is minimized. Because EDFT in its direct-diagonalization form requires an explicit Hamiltonian diagonalization each SCF cycle — reintroducing cubic scaling — ONETEP additionally provides the **AQuA-FOE** (Andersen-Quasi/Approximate, Fermi Operator Expansion) approach, which uses Chebyshev/hyperbolic-tangent expansions of the Fermi operator to construct finite-temperature density kernels for large metallic systems without diagonalization, restoring near-linear scaling for such systems.

### 2.5 Pseudopotentials and the PAW-in-linear-scaling formalism

ONETEP uses norm-conserving and, more importantly, its own adaptation of the **projector augmented-wave (PAW)** method, reformulated for compatibility with a non-orthogonal, localized-orbital, density-matrix-based framework (rather than the eigenstate-based PAW of conventional plane-wave codes). This allows access to quantities that depend sensitively on the wavefunction/density near the nucleus (e.g., core-level spectra, hyperfine parameters) while retaining the smooth, "soft" pseudo-density-matrix representation compatible with the psinc basis.

### 2.6 Forces, geometry optimization and molecular dynamics

Because psinc/NGWF localization regions are fixed spheres centred on atoms that must (in principle) move during geometry relaxation or dynamics, care is required to maintain smooth, translationally consistent potential-energy surfaces. ONETEP was explicitly designed so that on-the-fly NGWF optimization and use of the psinc basis produce smooth potential-energy surfaces whose ionic (Hellmann–Feynman) forces are consistent with the energy gradient, enabling reliable geometry optimization and both Born–Oppenheimer and (with appropriate propagation schemes) ab initio molecular dynamics on large systems.

### 2.7 Parallelization

ONETEP was designed from inception as a parallel code, originally within the MPI (Message Passing Interface) paradigm and now supporting hybrid MPI+OpenMP parallelism: MPI divides computational resources into processes (each typically holding a subset of atoms/NGWFs and the associated FFT boxes), which can be further subdivided into OpenMP threads occupying individual CPU cores. This decomposition, combined with the linear-scaling algorithm and the fixed-size FFT box, allows ONETEP to scale efficiently to very large processor counts on modern HPC systems. Recent releases have also introduced experimental GPU acceleration via OpenACC for several of the most computationally intensive kernels (structure factor evaluation, density calculation, local potential integrals, ion–ion summation, and partial Hartree–Fock exchange).

---

## 3. Key Capabilities

ONETEP has grown from a proof-of-concept linear-scaling total-energy code into a broad electronic-structure platform. Major capabilities include:

| Category | Features |
|---|---|
| **Electronic structure core** | Linear-scaling DFT via NGWFs/psinc basis; LNV density-kernel optimization; density-kernel truncation for asymptotic $\mathcal{O}(N)$ scaling; PAW and norm-conserving pseudopotentials |
| **Boundary conditions** | Fully periodic (bulk solids), 2D slab/surface, 1D wire/nanotube, and fully open (isolated molecule/cluster) boundary conditions within a single formalism |
| **Exchange–correlation** | LDA, GGA (e.g. PBE), meta-GGA, and hybrid functionals with (approximate, linear-scaling) exact/Hartree–Fock exchange; empirical dispersion corrections (e.g. DFT-D2/D3-type) |
| **Metals & finite temperature** | Ensemble DFT (EDFT) with Fermi–Dirac, Gaussian, and Methfessel–Paxton smearing; AQuA-FOE Fermi-operator-expansion method for large-scale metallic systems |
| **Strongly correlated electrons** | DFT+U (Dudarev and generalized/orbital-resolved formulations) with Hubbard projectors defined via NGWFs or Wannier-like functions; forces and phonons including +U contributions; ONETEP+TOSCAM coupling to Dynamical Mean-Field Theory (DMFT) for strongly correlated d/f-electron systems |
| **Excited states & spectroscopy** | Linear-response time-dependent DFT (LR-TDDFT) for optical excitations of large systems, using dedicated optimization of unoccupied-state (conduction) NGWFs; core-loss (EELS/XAS-type) spectroscopy; density of states and band-structure-type analysis (via projections, since no explicit Bloch bands are computed) |
| **Molecular dynamics & vibrations** | Born–Oppenheimer ab initio molecular dynamics; vibrational/phonon calculations |
| **Environment & embedding** | Minimal-parameter implicit solvation model (polarizable continuum-type, self-consistent with the electron density); QM/MM coupling (including with polarizable force fields such as AMOEBA); QM-in-QM embedding / embedded mean-field theory (EMFT) |
| **Charge/interaction partitioning & analysis** | Distributed multipole analysis (DMA); electron localization function (ELF) and localized orbital locator (LOL) descriptors; natural bond orbital (NBO)-type analysis; density-derived electrostatic and chemical (DDEC-type) charge partitioning; methods for partitioning energies/interactions between molecular fragments |
| **Constrained & tight-binding methods** | Constrained DFT (cDFT), including charge- and spin-constrained variants, for charge-transfer and diabatic-state studies; a self-consistent-charge Density Functional Tight-Binding (DFTB) mode for very rapid approximate calculations |
| **Electronic transport** | Methods for computing electronic transport properties of extended/large-scale systems |

---

## 4. Typical Workflow and Practical Characteristics

A ONETEP calculation is controlled via a plain-text `.dat` keyword input file (species definitions, pseudopotential files, NGWF radii and counts per species, kinetic-energy cutoff, kernel cutoff, and task-specific blocks). Unlike conventional plane-wave DFT, ONETEP calculations require simultaneous convergence of **two** interlinked quantities — the density kernel and the NGWFs — and results are only reliable once both the RMS NGWF gradient and the kernel/Hamiltonian commutator have converged; the ONETEP developers explicitly describe the code as "not yet quite... 'black box'" compared with conventional plane-wave codes, and recommend verbose output and convergence testing (particularly with respect to psinc cutoff energy, NGWF localization radius, and density-kernel cutoff) when approaching new types of system.

ONETEP has demonstrated calculations on systems including DNA fragments (~2600 atoms), carbon nanotubes (~4000 atoms), bulk silicon (~4096 atoms), point defects in Al₂O₃, full proteins and protein–ligand complexes (e.g. the ~2600-atom T4 lysozyme L99A/M102Q system), semiconductor nanorods and nanoparticles, and metalloproteins (e.g. ~1000-atom myoglobin models); reported applications extend to tens of thousands of atoms on suitable HPC resources, particularly for systems where regions of vacuum or low electronic complexity (surfaces, nanostructures, biomolecules) reduce the effective cost relative to dense, fully occupied bulk solids. The favourable crossover point (system size beyond which ONETEP outperforms cubic-scaling plane-wave codes) is lower for low-dimensional systems (molecules, surfaces, nanotubes) than for dense three-dimensional bulk solids.

---

## 5. Licensing and Availability

- ONETEP is distributed under an **academic user licence** administered by **Cambridge Enterprise Ltd** (University of Cambridge), free of charge to researchers with a permanent or fixed-term fellowship position at UK-based (and, at various points, other) academic institutions and UK national laboratories; the licence covers the licence-holder's research group at the same institution for research-only use.
- Non-academic/commercial use has historically been available via **BIOVIA/Accelrys Materials Studio**, and directly from the developers on a commercial basis.
- Official resources: project website `onetep.org`; documentation at `docs.onetep.org`; licensing portal via Cambridge Enterprise (`licensing.enterprise.cam.ac.uk/product/onetep`).
- Primary developing institutions: University of Cambridge, University of Southampton, Imperial College London, University of Warwick, and (more recently) Gdańsk University of Technology.

---

## 6. Summary Assessment

**Strengths**
- Genuine plane-wave-equivalent accuracy combined with linear ($\mathcal{O}(N)$) scaling, via a systematically improvable, variationally optimized minimal local-orbital basis — a combination rarely achieved simultaneously by competing linear-scaling codes.
- Very broad and continuously expanding feature set (DFT+U/DMFT, TDDFT, implicit solvent, QM/MM/QM-in-QM, cDFT, transport, spectroscopy), making it applicable across materials science, chemistry, and structural biology.
- Smooth potential energy surfaces enabling reliable geometry optimization and dynamics despite the use of localized orbitals.
- Long-standing, actively maintained academic development community with an openly documented codebase (`docs.onetep.org`) and tutorials.

**Limitations / considerations**
- Requires more user expertise and careful convergence testing than routine plane-wave DFT, owing to the coupled convergence of NGWFs and the density kernel and the sensitivity of results to localization-region and kernel-cutoff choices.
- True asymptotic linear scaling requires finite density-kernel truncation, which introduces an additional convergence parameter not present in conventional DFT; achieving both accuracy and favourable scaling for metallic or highly delocalized systems needs finite-temperature/EDFT or AQuA-FOE methods rather than the default LNV scheme.
- Academic licensing (rather than fully open-source distribution) limits unrestricted redistribution and community-driven code contribution relative to some open-source DFT packages.
- The crossover system size beyond which ONETEP is more efficient than cubic-scaling plane-wave codes depends strongly on dimensionality and density of the system, so it does not uniformly outperform conventional DFT codes for small-to-moderate, dense three-dimensional systems.

---

## 7. Key Publications on ONETEP Theory and Methodology

### 7.1 Foundational method and code overview papers

- P. D. Haynes, C.-K. Skylaris, A. A. Mostofi, M. C. Payne, "ONETEP: linear-scaling density-functional theory with plane waves," *Psi-k Newsletter* **72**, 78 (2005).
- C.-K. Skylaris, P. D. Haynes, A. A. Mostofi, M. C. Payne, "Introducing ONETEP: Linear-scaling density functional simulations on parallel computers," *Journal of Chemical Physics* **122**, 084119 (2005). doi:10.1063/1.1839852
- C.-K. Skylaris, P. D. Haynes, A. A. Mostofi, M. C. Payne, "Using ONETEP for accurate and efficient O(N) density functional calculations," *Journal of Physics: Condensed Matter* **17**, 5757 (2005).
- N. D. M. Hine, P. D. Haynes, A. A. Mostofi, C.-K. Skylaris, M. C. Payne, "Linear-scaling density-functional theory with tens of thousands of atoms: Expanding the scope and scale of calculations with ONETEP," *Computer Physics Communications* **180**, 1041 (2009). doi:10.1016/j.cpc.2008.12.023
- J. C. A. Prentice, J. Aarons, J. C. Womack, A. E. A. Allen, L. Andrinopoulos, L. Anton, R. A. Bell, A. Bhandari, G. A. Bramley, R. J. Charlton, R. J. Clements, D. J. Cole, G. Constantinescu, F. Corsetti, S. M.-M. Dubois, K. K. B. Duff, J. M. Escartín, A. Greco, Q. Hill, L. P. Lee, E. Linscott, D. D. O'Regan, M. J. S. Phipps, L. E. Ratcliff, Á. Ruiz Serrano, E. W. Tait, G. Teobaldi, V. Vitale, N. Yeung, T. J. Zuehlsdorff, J. Dziedzic, P. D. Haynes, N. D. M. Hine, A. A. Mostofi, M. C. Payne, C.-K. Skylaris, "The ONETEP linear-scaling density functional theory program," *Journal of Chemical Physics* **152**, 174111 (2020). doi:10.1063/5.0004445 — the current, most comprehensive overview/review paper for the code.

### 7.2 NGWFs, basis set, and density-kernel optimization

- C.-K. Skylaris, A. A. Mostofi, P. D. Haynes, O. Diéguez, M. C. Payne, "Nonorthogonal generalized Wannier function pseudopotential plane-wave method," *Physical Review B* **66**, 035119 (2002).
- A. A. Mostofi, C.-K. Skylaris, P. D. Haynes, M. C. Payne, "Total-energy calculations on a real space grid with localized functions and a plane-wave basis," *Computer Physics Communications* **147**, 788 (2002).
- C.-K. Skylaris, A. A. Mostofi, P. D. Haynes, C. J. Pickard, M. C. Payne, "Accurate kinetic energy evaluation in electronic structure calculations with localized functions on real space grids," *Computer Physics Communications* **140**, 315 (2001).
- Á. Ruiz-Serrano, C.-K. Skylaris, "A variational method for density functional theory calculations on metallic systems with thousands of atoms," *Journal of Chemical Physics* **139**, 054107 (2013). (Density-kernel/EDFT methodology.)
- C.-K. Skylaris, P. D. Haynes, "Achieving plane wave accuracy in linear-scaling density functional theory applied to periodic systems: A case study on crystalline silicon," *Journal of Chemical Physics* **127**, 164712 (2007).
- L. Anton, C.-K. Skylaris, "Density kernel optimization in the ONETEP code," *Journal of Physics: Condensed Matter* **20**, 294207 (2008).

### 7.3 Ensemble DFT, finite temperature, and metallic systems

- Á. Ruiz-Serrano, C.-K. Skylaris, "A variational method for density functional theory calculations on metallic systems with thousands of atoms," *Journal of Chemical Physics* **139**, 054107 (2013).
- J. Aarons, M. Sarwar, D. Thompsett, C.-K. Skylaris, "Perspective: Methods for large-scale density functional calculations on metallic systems," *Journal of Chemical Physics* **145**, 220901 (2016). doi:10.1063/1.4972007

### 7.4 PAW formalism in linear-scaling DFT

- N. D. M. Hine, "Linear-scaling density functional theory using the projector augmented wave method," *Journal of Physics: Condensed Matter* **29**, 024001 (2017).

### 7.5 DFT+U, Hubbard corrections, and DMFT coupling

- D. D. O'Regan, M. C. Payne, A. A. Mostofi, "Subspace representations in ab initio methods for strongly correlated systems," *Physical Review B* **83**, 245124 (2011).
- D. D. O'Regan, N. D. M. Hine, M. C. Payne, A. A. Mostofi, "Projector self-consistent DFT+U using nonorthogonal generalized Wannier functions," *Physical Review B* **82**, 081102(R) (2010).
- D. J. Cole, D. D. O'Regan, M. C. Payne, "Ligand Discrimination in Myoglobin from Linear-Scaling DFT+U," *Journal of Physical Chemistry Letters* **3**, 1448 (2012).
- E. B. Linscott, D. J. Cole, N. D. M. Hine, M. C. Payne, C. Weber, "ONETEP + TOSCAM: Uniting Dynamical Mean Field Theory and Linear-Scaling Density Functional Theory," *Journal of Chemical Theory and Computation* **16**, 4899 (2020). doi:10.1021/acs.jctc.0c00162

### 7.6 Time-dependent DFT and optical/core-level spectroscopy

- L. E. Ratcliff, N. D. M. Hine, P. D. Haynes, "Calculating optical absorption spectra for large systems using linear-scaling density functional theory," *Physical Review B* **84**, 165131 (2011).
- T. J. Zuehlsdorff, N. D. M. Hine, J. S. Spencer, N. M. Harrison, D. J. Riley, P. D. Haynes, "Linear-scaling time-dependent density-functional theory in the linear response formalism," *Journal of Chemical Physics* **139**, 064104 (2013).

### 7.7 Implicit solvation and QM/MM/embedding

- J. Dziedzic, H. H. Helal, C.-K. Skylaris, A. A. Mostofi, M. C. Payne, "Minimal parameter implicit solvent model for ab initio electronic-structure calculations," *Europhysics Letters* **95**, 43001 (2011).
- J. Dziedzic, S. J. Fox, T. Fox, C. S. Tautermann, C.-K. Skylaris, "Large-scale DFT calculations in implicit solvent — a case study on the T4 lysozyme L99A/M102Q protein," *International Journal of Quantum Chemistry* **113**, 771 (2013).
- S. J. Fox, J. Dziedzic, T. Fox, C. S. Tautermann, C.-K. Skylaris, "Density functional theory calculations on entire proteins for free energies of binding: Application to a model polar binding site," *Proteins: Structure, Function, and Bioinformatics* **82**, 3335 (2014).

### 7.8 Analysis and charge-partitioning tools

- R. J. Clements, J. C. Womack, C.-K. Skylaris, "Electron localisation descriptors in ONETEP: a tool for interpreting localisation and bonding in large-scale DFT calculations," *Electronic Structure* **2**, 027001 (2020). doi:10.1088/2516-1075/ab8d19

---

*Compiled from the ONETEP developers' publications, the official ONETEP documentation (`docs.onetep.org`), the ONETEP project website (`onetep.org`), and associated peer-reviewed literature.*

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the ONETEP 	Linear-scaling plane-wave-based DFT code using local orbitals ("psinc" functions), designed for very large systems. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
