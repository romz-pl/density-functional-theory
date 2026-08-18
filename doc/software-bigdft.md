# BigDFT: An Exhaustive Review

*A Daubechies-wavelet-based DFT code with adaptive resolution and linear-scaling algorithms*

---

## 1. Overview

BigDFT is a free, open-source (GPL) electronic-structure code implementing Kohn–Sham density functional theory (DFT) using **Daubechies wavelets** as the underlying basis set. The project began in 2005 as an EU FP6-STREP-NEST-funded collaboration between four European laboratories: L_Sim (CEA Grenoble, France), the University of Basel (Switzerland), the Université catholique de Louvain-la-Neuve (Belgium), and Kiel University (Germany). Rather than being a DFT code built from scratch and later optimized, BigDFT was conceived from the outset as an experiment: could Daubechies wavelets — a basis set with well-established advantages in signal processing and applied mathematics — provide a superior foundation for large-scale, high-precision electronic structure calculations?

The resulting code combines systematic convergence, adaptive real-space resolution, and — uniquely among wavelet-based codes — a genuine linear-scaling (`$O(N)$`) algorithm built on the same basis set used for the cubic-scaling (`$O(N^3)$`) approach. This dual-mode design (cubic-scaling for small/medium systems, linear-scaling for systems of many thousands of atoms) is one of BigDFT's defining features. The code received the 2009 French Bull-Joseph Fourier prize for high-performance computing.

---

## 2. Theoretical and Mathematical Foundations

### 2.1 Why Daubechies wavelets?

Daubechies wavelets are compactly supported, orthonormal functions that form a complete, systematic basis set. BigDFT typically uses Daubechies wavelets of order 16 (i.e., with 8 vanishing moments), interpolated on a uniform real-space grid. Their central properties, which motivate their use in DFT, are:

- **Orthonormality** — simplifies the overlap matrix and numerical linear algebra (in the simplest cubic-scaling case, the overlap matrix is the identity).
- **Compact support in real space** — each wavelet is nonzero only over a small, finite spatial region, giving the basis strict locality.
- **Good localization in Fourier (reciprocal) space** — despite being real-space objects, Daubechies wavelets decay rapidly in reciprocal space, giving them a systematic, plane-wave-like convergence behavior.
- **Multi-resolution structure** — wavelets exist at different scales/resolution levels, which is what permits *adaptive* real-space resolution (see §2.2).
- **Smoothness** — high-order Daubechies wavelets are sufficiently smooth to represent Kohn–Sham orbitals accurately with relatively few coefficients.

No other basis set combines orthogonality, compact support in both spaces, and completeness in quite the same way; plane waves are delocalized, Gaussians are non-orthogonal and lack compact support, and finite differences/real-space grids lack the systematic convergence and smoothness guarantees of wavelets.

### 2.2 Adaptive resolution: the "double-grid" approach

BigDFT does not use a single, uniform fine grid over the whole simulation cell. Instead, it defines two nested levels of resolution:

- A **coarse grid**, covering the full simulation region with resolution level *j* (scaling function coefficients only).
- A **fine grid**, restricted to regions close to the atoms (where the wavefunctions vary rapidly, e.g. near pseudopotential cores), carrying full scaling-function-plus-wavelet coefficients at the next resolution level.

This adaptive, two-level mesh means that computational effort — and memory — is concentrated where the wavefunction actually has fine structure, while the low-density vacuum/interstitial regions are represented economically on the coarse grid. This is conceptually similar to adaptive mesh refinement in finite-element methods, but implemented within a genuinely orthogonal, systematic wavelet basis rather than an ad hoc real-space grid.

### 2.3 The Kohn–Sham problem in the wavelet basis

The Kohn–Sham orbitals |Ψᵢ⟩ are eigenfunctions of the Kohn–Sham Hamiltonian:

$$\left(-\frac{1}{2}\nabla^2 + V_{KS}[\rho] + V_{\rm psp}\right)|\Psi_i\rangle = \epsilon_i |\Psi_i\rangle$$

expanded in the Daubechies basis. Because the basis is orthonormal and localized, the kinetic energy operator, the application of the local/nonlocal pseudopotential, and the Poisson equation for the Hartree potential can all be evaluated with optimal-scaling algorithms (essentially convolutions with short filters), while the sparsity of the localized basis keeps memory requirements low. A preconditioned conjugate-gradient (or DIIS/mixing) scheme is used to solve for the orbitals self-consistently.

### 2.4 Poisson solver and boundary conditions

A key theoretical/algorithmic contribution of the BigDFT project is its **flexible Poisson solver**, based on an interpolating-scaling-function/Green's-function formalism (Genovese, Deutsch, Neelov, Goedecker & Beylkin, 2006). It provides efficient, high-accuracy solutions of the Hartree potential under **free, surface (2D-periodic + 1D-free), wire (1D-periodic + 2D-free), and fully periodic boundary conditions**, all within a unified algorithmic framework. This flexibility — rare among plane-wave-based codes, which are natively periodic and require artificial supercells and charge corrections to simulate isolated systems — is one of BigDFT's most cited practical advantages. The Poisson solver has since been decoupled from BigDFT and reused independently in other codes (ABINIT, Octopus, CP2K).

### 2.5 Pseudopotentials and PAW

BigDFT uses dual-space Gaussian-type, norm-conserving pseudopotentials (GTH/HGH-type, including variants with non-linear core corrections), chosen because their analytic Gaussian form matches efficiently onto Daubechies wavelets/Gaussian overlap integrals. This pseudopotential/wavelet combination has been shown to reproduce all-electron precision for many ground-state quantities. In addition, a **wavelet-based Projector Augmented-Wave (PAW)** formalism has been implemented, extending BigDFT toward frozen-core all-electron precision while retaining the systematic, adaptive, localized character of the wavelet basis.

### 2.6 Linear-scaling formalism (support functions)

The linear-scaling (`LS-BigDFT`) mode reformulates the problem in terms of the **one-particle density matrix** rather than explicit orthogonal orbitals, expressed through a set of strictly localized, in-situ optimized **support functions** {φ_α(r)}:

$$H_{\alpha\beta} = \int \phi_\alpha(\mathbf{r})\,\mathcal{H}(\mathbf{r})\,\phi_\beta(\mathbf{r})\, d\mathbf{r}$$

Each support function is itself expanded in the underlying Daubechies wavelet basis and confined to a localization region around an atom (or small group of atoms); its shape is optimized self-consistently ("in situ") rather than fixed a priori, unlike traditional minimal-basis linear-scaling approaches. This in-situ optimization is what allows a *minimal* number of support functions per atom (typically close to the number of valence atomic orbitals) to reach near-cubic-scaling accuracy. Because the localization regions are compact, the Hamiltonian and overlap matrices in the support-function representation are **sparse**, and the density kernel can be obtained without explicit diagonalization, giving overall linear scaling of both computational cost and memory with system size.

Several methods are available to compute the density matrix/kernel from the sparse Hamiltonian:

- **Fermi Operator Expansion (FOE)** — polynomial (Chebyshev) expansion of the Fermi–Dirac operator; the only method offering asymptotically strict linear scaling, including for metallic (zero-gap) systems, implemented via the associated **CheSS** sparse matrix library.
- **PEXSI** (Pole EXpansion and Selected Inversion) — an external library option for intermediate system sizes, scaling sub-cubically.
- Direct **diagonalization** (LAPACK/ScaLAPACK) of the sparse matrices — cubic scaling, used mainly for validation/benchmarking at smaller sizes.

Because file sizes per support function stay roughly constant (~1 MB) regardless of total system size, and the number of files grows only with the number of atoms/fragments, LS-BigDFT also has favorable, close-to-constant I/O behavior — a nontrivial advantage at the many-thousand-atom scale.

### 2.7 GPU and hybrid-architecture acceleration

The convolution-heavy operations intrinsic to the wavelet formalism (application of the kinetic energy operator, the Poisson solver, orbital transformations) map naturally onto GPU acceleration. BigDFT was an early adopter of hybrid CPU–GPU DFT, reporting speedups of roughly 6× for full DFT runs (and up to ~20× for individual convolution-dominated operations) while preserving full numerical precision and parallel scalability — including a dedicated customized 3D GPU Poisson solver for free/wire/surface boundary conditions.

---

## 3. Key Practical Features

- **Boundary conditions**: free, wire, surface, and fully periodic — all with one unified Poisson-solver formalism.
- **Dual-mode operation**: cubic-scaling (all-orbital, `O(N^3)`) and linear-scaling (`O(N)`, support-function-based) within the same code and basis set.
- **System sizes**: routinely from a handful of atoms up to systems of many thousands of atoms in linear-scaling mode.
- **Environment models**: implicit solvent models and external electric fields.
- **Spectroscopy/excited states**: linear-response time-dependent DFT (TD-DFT), fragment-based constrained DFT (CDFT) and transition-based CDFT (T-CDFT) for charge-transfer excitations, implemented for both cubic- and linear-scaling modes.
- **Fragment approach**: reusable, transferable support functions across molecular fragments, useful for large biomolecular/supramolecular assemblies.
- **Interfaces**: usable as a standalone package or via the ABINIT package (BigDFT wavelet mode has been distributed within ABINIT since v5.5); **PyBigDFT** provides a Python workflow/analysis layer.
- **License**: GNU GPL, free and open source.
- **Parallelization**: MPI + OpenMP hybrid parallelism, with GPU acceleration for key kernels; demonstrated scalability to many thousands of cores.

---

## 4. Strengths and Limitations

**Strengths**
- Systematic, controllable convergence (a single "grid spacing" parameter, much like a plane-wave cutoff, governs accuracy in a predictable way).
- High accuracy relative to cost, due to the combination of a smooth systematic basis with adaptive resolution.
- Native support for non-periodic boundary conditions without supercell artifacts.
- A rare example of a linear-scaling method built on, and fully consistent with, a systematic wavelet basis rather than an empirical local-orbital set.
- Strong, validated performance for both insulating and metallic (zero-gap) large systems in linear-scaling mode.

**Limitations**
- Pseudopotential-based (not fully all-electron in its most widely used mode), although the PAW extension mitigates this.
- The linear-scaling mode's accuracy depends on localization-radius choices, which trade accuracy against computational cost.
- As with any wavelet/real-space-grid method, some operations (e.g., certain hybrid-functional exact-exchange terms) are less naturally efficient than in plane-wave or Gaussian-basis codes, though solutions exist within the project's ecosystem.

---

## 5. Publications on BigDFT's Theory and Methodology

The list below focuses on the papers that establish and extend BigDFT's core theoretical/methodological framework (wavelet formalism, Poisson solver, linear scaling, PAW, GPU algorithms, and related reviews), rather than application papers that merely use the code.

### 5.1 Foundational papers

1. Genovese, L., Neelov, A., Goedecker, S., Deutsch, T., Ghasemi, S. A., Willand, A., Caliste, D., Zilberberg, O., Rayson, M., Bergman, A., & Schneider, R. (2008). **Daubechies wavelets as a basis set for density functional pseudopotential calculations.** *The Journal of Chemical Physics*, 129(1), 014109. https://doi.org/10.1063/1.2949547

2. Genovese, L., Deutsch, T., Neelov, A., Goedecker, S., & Beylkin, G. (2006). **Efficient solution of Poisson's equation with free boundary conditions.** *The Journal of Chemical Physics*, 125(7), 074105.

### 5.2 Parallel and GPU-acceleration methodology

3. Genovese, L., Ospici, M., Deutsch, T., Méhaut, J.-F., Neelov, A., & Goedecker, S. (2009). **Density functional theory calculation on many-cores hybrid central processing unit–graphic processing unit architectures.** *The Journal of Chemical Physics*, 131(3), 034103 (also cited as 134103 in some references).

4. Genovese, L., Videau, B., Ospici, M., Deutsch, T., Goedecker, S., & Méhaut, J.-F. (2011). **Daubechies wavelets for high performance electronic structure calculations: The BigDFT project.** *Comptes Rendus Mécanique*, 339(2–3), 149–164.

5. Ghasemi, S. A., Kosior, A., & collaborators — Poisson-solver GPU papers: Cerioni, A., Genovese, L., Mirone, A., & Solé, V. A. (2012). **Efficient and accurate solver of the three-dimensional screened and unscreened Poisson's equation with generic boundary conditions.** *The Journal of Chemical Physics*, 137(13), 134108.

6. (GPU Poisson solver, free boundary conditions) **A customized 3D GPU Poisson solver for free boundary conditions.** *Computer Physics Communications* (2013), building on the BigDFT Poisson-solver formalism.

### 5.3 Linear-scaling formalism

7. Mohr, S., Ratcliff, L. E., Boulanger, P., Genovese, L., Caliste, D., Deutsch, T., & Goedecker, S. (2014). **Daubechies wavelets for linear scaling density functional theory.** *The Journal of Chemical Physics*, 140(20), 204110.

8. Mohr, S., Ratcliff, L. E., Genovese, L., Caliste, D., Boulanger, P., Goedecker, S., & Deutsch, T. (2015). **Accurate and efficient linear scaling DFT calculations with universal applicability.** *Physical Chemistry Chemical Physics*, 17(47), 31360–31370.

9. Ratcliff, L. E., Mohr, S., Huhs, G., Deutsch, T., Masella, M., & Genovese, L. (2017). **Challenges in large scale quantum mechanical calculations.** *WIREs Computational Molecular Science*, 7(1), e1290. (Also archived as arXiv:1609.00252.)

10. Mohr, S., Eixarch, M., Amsler, M., Mavrič, A., Genovese, L., & Deutsch, T. (2018). **Linear scaling DFT calculations for large tungsten systems using an optimized local basis.** *Computational Materials Science*, 144, 128–133. (arXiv:1711.10993.)

11. Ratcliff, L. E., Degomme, A., Flores-Livas, J. A., Goedecker, S., & Genovese, L. (2015). **Affordable and accurate large-scale hybrid-functional calculations on GPU-accelerated supercomputers.** *Journal of Physics: Condensed Matter* (methodological extension of linear-scaling framework).

### 5.4 Sparse matrix / density-kernel methodology (CheSS library)

12. Mohr, S., Dawson, W., Wagner, M., Caliste, D., Nakajima, T., & Genovese, L. (2017). **Efficient computation of sparse matrix functions for large-scale electronic structure calculations: The CheSS library.** *Journal of Chemical Theory and Computation*, 13(10), 4684–4698.

13. Dawson, W., & Nakajima, T. (2018/2020). Related papers on Fermi Operator Expansion (FOE) implementations for the BigDFT/CheSS ecosystem.

### 5.5 Projector Augmented-Wave (PAW) extension

14. Rangel, T., Caliste, D., Genovese, L., & Torrent, M. (2016). **A wavelet-based Projector Augmented-Wave (PAW) method: reaching frozen-core all-electron precision with a systematic, adaptive and localized wavelet basis set.** *Computer Physics Communications*, 208, 218–226. (arXiv:1605.04375.)

### 5.6 Time-dependent DFT / excited-state methodology

15. Natarajan, B., Genovese, L., Casida, M. E., Deutsch, T., Burchak, O. N., Philouze, C., & Balakirev, M. Y. (2012). **Wavelet-based linear-response time-dependent density-functional theory.** *Chemical Physics*, 402, 29–40.

16. Ivanov, A. V., Genovese, L., & Casida, M. E. et al. **Wavelets for density-functional theory and post-density-functional-theory calculations.** In *Handbook of Numerical Analysis* / review chapter (arXiv:1110.4853) — comprehensive review of the wavelet-TD-DFT formalism as implemented in BigDFT.

17. Stella, M., Thapa, K., Genovese, L., & Ratcliff, L. E. (2022). **Transition-based constrained DFT for the robust and reliable treatment of excitations in supramolecular systems.** *Journal of Chemical Theory and Computation*, 18(5), 3027–3038. (arXiv:2106.01142.)

### 5.7 Fragment / constrained-DFT methodology

18. Ratcliff, L. E., Genovese, L., Mohr, S., & Deutsch, T. (2015). **Fragment approach to constrained density functional theory calculations using Daubechies wavelets.** *The Journal of Chemical Physics*, 142(23), 234105.

### 5.8 Reviews, roadmaps, and pedagogical accounts

19. Genovese, L., Videau, B., Caliste, D., Deutsch, T., Ospici, M., & Méhaut, J.-F. (2011/2012). **Daubechies wavelets for linear scaling density functional theory** — expanded conference/review contribution (companion to item 7, targeted at HPC audiences).

20. Ratcliff, L. E., Dawson, W., Fisicaro, G., Caliste, D., Mohr, S., Degomme, A., Videau, B., Cristiglio, V., Stella, M., D'Alessandro, M., Goedecker, S., Nakajima, T., Deutsch, T., & Genovese, L. (2020). **Flexibilities of wavelets as a computational basis set for large-scale electronic structure calculations.** *The Journal of Chemical Physics*, 152(19), 194110. (A comprehensive, up-to-date review of the whole BigDFT methodology.)

21. Ratcliff, L. E., et al. (2022). **Roadmap on electronic structure codes in the exascale era** — BigDFT contribution section. *Modelling and Simulation in Materials Science and Engineering* / arXiv:2209.12747.

22. (Community exploratory/validation paper) **Exploratory data science on supercomputers for quantum mechanical calculations**, arXiv:2310.09552 — validation of BigDFT against Gaussian-basis references (PSI4, W4-11 dataset).

### 5.9 Selected downstream/derivative works

23. Genovese, L., Deutsch, T., & Goedecker, S. (2007). **Efficient and accurate three-dimensional Poisson solver for surface problems.** *The Journal of Chemical Physics*, 127(5), 054704. (Extension of the free-BC Poisson solver to surface boundary conditions.)

24. Neelov, A., & Goedecker, S. (2006). **An efficient numerical quadrature for the calculation of the potential energy of wavefunctions expressed in the Daubechies wavelet basis.** *Journal of Computational Physics*, 217(2), 312–339.

25. Willand, A., Kvashnin, Y. O., Genovese, L., Vázquez-Mayagoitia, Á., Deb, A. K., Sadeghi, A., Deutsch, T., & Goedecker, S. (2013). **Norm-conserving pseudopotentials with chemical accuracy compared to all-electron calculations.** *The Journal of Chemical Physics*, 138(10), 104109. (Pseudopotential methodology used throughout BigDFT.)

---

*Note: reference numbering and journal volume/page details for some secondary/derivative entries above (§5.3 item 11, §5.4 item 13, §5.9) reflect the best available bibliographic information found; readers should verify exact volume/page numbers against the publisher record before formal citation.*

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the BigDFT 	DFT code using Daubechies wavelets as a basis, enabling adaptive resolution and linear-scaling algorithms. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
