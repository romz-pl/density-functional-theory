# PROFESS: An Exhaustive Review of the Orbital-Free DFT Code

## 1. Overview

**PROFESS** (**PR**inceton **O**rbital-**F**ree **E**lectronic **S**tructure **S**oftware) is a first-principles electronic structure package that implements **orbital-free density functional theory (OF-DFT)**. Unlike conventional Kohn–Sham DFT (KS-DFT), which represents the non-interacting kinetic energy through a set of single-particle orbitals that must be computed, orthogonalized, and diagonalized self-consistently, OF-DFT approximates the kinetic energy as an *explicit functional of the electron density itself*, $T_s[n(\mathbf{r})]$, via a **kinetic energy density functional (KEDF)**. Because no orbitals are ever constructed, the cubic-scaling bottlenecks of KS-DFT (orthonormalization, diagonalization) are eliminated, and the total energy, forces, and stresses can be evaluated with (quasi-)linear scaling in system size. This makes PROFESS one of the primary tools for mesoscale, "large but still quantum mechanical," simulations — historically pushing electronic-structure calculations up to millions of atoms on modest computing resources.

PROFESS originated in Prof. Emily A. Carter's research group (UCLA, then Princeton University), building on earlier OF-DFT algorithm work by Stuart C. Watson. It has since evolved through several major versions and is now maintained as an open-source project independent of any single institution.

---

## 2. Theoretical Foundation

### 2.1 Orbital-free vs. Kohn–Sham DFT

Both approaches rest on the Hohenberg–Kohn theorems, which guarantee that the ground-state energy of an interacting electron system is a unique functional of the density $n(\mathbf{r})$ alone:

$$E[n] = T_s[n] + E_H[n] + E_{xc}[n] + \int v_{ext}(\mathbf{r})\, n(\mathbf{r})\, d\mathbf{r}$$

Kohn–Sham DFT sidesteps the unknown non-interacting kinetic energy functional $T_s[n]$ by reintroducing a fictitious set of single-particle orbitals $\{\psi_i\}$ such that $T_s = -\tfrac{1}{2}\sum_i \int \psi_i^\ast \nabla^2 \psi_i\, d\mathbf{r}$, which is exact but requires solving $N_{\text{electron}}$ coupled eigenvalue equations self-consistently — the source of cubic ($O(N^3)$) computational scaling for the orthogonalization/diagonalization steps as system size $N$ grows.

OF-DFT instead evaluates $T_s[n]$ directly as an explicit approximate functional of the density, bypassing orbitals entirely. The ground-state density is obtained by direct unconstrained (or singly-constrained, for particle-number conservation) minimization of the total energy functional with respect to $n(\mathbf{r})$ (or, in PROFESS's actual implementation, with respect to $\phi(\mathbf{r}) = \sqrt{n(\mathbf{r})}$, which automatically enforces non-negativity of the density). This removes the wavefunction-orthogonality bottleneck and allows the electronic energy and potential terms to scale essentially linearly with system size.

### 2.2 Kinetic Energy Density Functionals (KEDFs)

The accuracy of any OF-DFT calculation is dominated by the quality of the approximate KEDF, since the kinetic energy is typically the largest energy term. PROFESS implements a hierarchy of KEDFs of increasing sophistication:

- **Thomas–Fermi (TF) functional** — the local-density approximation to the kinetic energy, exact in the limit of a slowly varying, uniform electron gas:
$$T_{TF}[n] = C_F \int n(\mathbf{r})^{5/3}\, d\mathbf{r}, \qquad C_F = \frac{3}{10}(3\pi^2)^{2/3}$$

- **von Weizsäcker (vW) functional** — exact for a single-orbital (or two-electron closed-shell) system, and dominant in regions of rapidly varying density (e.g., near nuclei or in the vacuum tail):
$$T_{vW}[n] = \frac{1}{8}\int \frac{|\nabla n(\mathbf{r})|^2}{n(\mathbf{r})}\, d\mathbf{r}$$

- **TF$\lambda$vW / TFvW hybrids** — linear combinations of the above two limiting forms.

- **Nonlocal, linear-response-based KEDFs** — the workhorse functionals for metallic systems, constructed so that the KEDF reproduces the exact linear response (Lindhard) function of the homogeneous electron gas:
  - **Wang–Teter (WT) functional** — the first widely used nonlocal KEDF of this class, using a density-independent convolution kernel.
  - **Wang–Govind–Carter (WGC) functional** — generalizes WT by introducing a **density-dependent response kernel**, giving substantially improved accuracy for nearly-free-electron metals such as Al, Mg, and their alloys. An analytic closed form of the WGC kernel (amenable to efficient reciprocal-space evaluation under both periodic and Dirichlet/aperiodic boundary conditions) was later derived, enabling its use in PROFESS's plane-wave/FFT framework.
  - **Perrot** and **Smargiassi–Madden** functionals — related nonlocal forms sharing the same linear-response ansatz.
  - **Huang–Carter (HC) functional** — extends the nonlocal approach to **semiconductors**, whose more inhomogeneous, gapped electron densities are poorly described by the metal-oriented WT/WGC kernels; the kernel depends on both the density and the reduced density gradient.
  - **Enhanced von Weizsäcker WGC (EvW-WGC)** — a further refinement combining WGC with a locally enhanced vW term, improving transferability across the metal-to-semiconductor range.

- **Single-point (semilocal) generalized-gradient/Laplacian-level KEDFs** — including gradient-expansion and Pauli-potential-based forms suited for finite-temperature (free-energy) simulations at warm-dense-matter conditions, where nonlocal kernels become computationally unwieldy at high electron temperature.

All of these nonlocal kernels rely fundamentally on the **linear response function of the homogeneous electron gas (Lindhard function)**, which is why they perform best for "nearly-free-electron" metals whose valence electron density resembles a perturbed homogeneous gas, and progressively worse for more localized, covalently bonded, or strongly inhomogeneous systems (transition metals, molecules, semiconductors away from their bulk-like environment).

### 2.3 Local Pseudopotentials

A defining constraint of OF-DFT is that, because no orbitals exist, the **angular-momentum-dependent projection operators used in standard nonlocal pseudopotentials (NLPPs)** cannot be applied — there is nothing for them to project onto. PROFESS therefore requires **local pseudopotentials (LPPs)**: simple multiplicative operators, identical for all angular momentum channels, that regularize the electron–nucleus Coulomb singularity while removing chemically inert core electrons. Local pseudopotentials are generally more accurate for main-group, nearly-free-electron metals (Al, Mg, Li, Na) with isotropic electron distributions, and less accurate for transition metals and semiconductors, whose density has genuine angular structure that an NLPP would capture but an LPP cannot. Advances such as **bulk-derived local pseudopotentials (BLPS)**, constructed by inverting Kohn–Sham equations on bulk valence densities, and optimized effective local pseudopotentials (OEPP), have substantially improved LPP quality and expanded PROFESS's useful element coverage (with published, tested LPPs available for elements including Al, Mg, Si, Li, Na, and select others).

### 2.4 Numerical Method

PROFESS represents the electron density (as $\phi = \sqrt{n}$) on a **real-space grid under periodic boundary conditions**, with reciprocal-space (plane-wave/FFT-based) evaluation of the kinetic, Hartree, and nonlocal-KEDF kernel terms. The ground-state density is found by direct minimization of the OF-DFT total-energy functional using truncated-Newton, conjugate-gradient, or related optimization algorithms, rather than a self-consistent-field eigenvalue procedure. Energies, atomic forces, and the stress tensor are all implemented analytically, enabling geometry optimization, unit-cell relaxation, and (in later versions) *ab initio* molecular dynamics.

---

## 3. Capabilities and Computational Scaling

- **(Quasi-)linear scaling**: all purely electronic energy and potential terms (kinetic, Hartree, exchange-correlation, electron–ion) scale linearly, or close to it ($O(N\log N)$ for the FFT-based nonlocal KEDF convolutions), with the number of atoms/grid points.
- **Quadratic-scaling ion terms**: the direct real-space Ewald-type ion–ion and ion–electron interaction terms exhibit quadratic scaling with the number of ions in the originally published implementation, which becomes the eventual bottleneck at very large atom counts (later versions introduced linear-scaling treatments of these terms as well).
- **Single-processor million-atom demonstrations**: the original 2008 release could treat tens of thousands of atoms with quantum mechanics on a single processor; subsequent parallelization (PROFESS 2.0) enabled the celebrated **million-atom bulk aluminum** calculation.
- **Molecular dynamics**: PROFESS 3.0 added *ab initio* MD capability driven by OF-DFT forces, including thermostatted (canonical-ensemble) dynamics, enabling finite-temperature simulations (e.g., of liquid metals) that would be prohibitively expensive with KS-DFT MD.
- **Finite-temperature / warm dense matter**: coupling with Quantum Espresso enabled consistent OF-DFT-driven and KS-DFT-driven *ab initio* MD "on the same footing" across a huge temperature range (thousands to millions of Kelvin), with reported speedups of a few times to several hundred times relative to KS-DFT MD depending on temperature and system size.

### Practical limitations

- **Accuracy is bounded by KEDF quality.** OF-DFT is formally exact (per Hohenberg–Kohn), but no currently known KEDF is exact in practice; PROFESS's nonlocal functionals are most reliable for **main-group, nearly-free-electron metals** (Al, Mg, Li, Na, K and their alloys) and progressively less reliable for **transition metals** and **semiconductors**, though functionals such as Huang–Carter partially close this gap for semiconductors.
- **Pseudopotential restriction.** The requirement for local pseudopotentials limits the set of elements with well-validated potentials and generally excludes straightforward, highly accurate treatment of transition metals.
- **Not (originally) suited to molecules/finite systems** in the same way as extended solids, since the linear-response KEDFs are derived from the homogeneous electron gas; molecular/finite-system OF-DFT is an active but harder research frontier (addressed partly through PROFESS-based studies and separate tools).

---

## 4. Version History

| Version | Year | Key contributors | Principal advances |
|---|---|---|---|
| **Precursor algorithms** | 2000 | S. C. Watson, E. A. Carter | Linear-scaling parallel algorithms for first-principles treatment of metals; foundational algorithmic groundwork developed during Watson's postdoc in the Carter group (UCLA), following his PhD with Paul Madden's group at Oxford. |
| **PROFESS 1.0** | 2008 | G. S. Ho, V. L. Lignères, E. A. Carter | First official public release (Princeton). Implemented energy, force, and stress functionals under periodic boundary conditions; TF, vW, TFvW$\alpha$, Wang–Teter, and Wang–Govind–Carter KEDFs; local pseudopotential support; demonstrated tens of thousands of atoms on a single processor. Written in Fortran 90; distributed via the Computer Physics Communications (CPC) Program Library. |
| **PROFESS 2.0** | 2010 | L. Hung, C. Huang, I. Shin, G. S. Ho, V. L. Lignères, E. A. Carter | Parallelization and fully linear-scaling algorithms, enabling treatment of over one million atoms of bulk aluminum with explicit quantum mechanics — a landmark demonstration of OF-DFT scalability (building on Hung & Carter's separate million-atom demonstration paper). |
| **PROFESS 3.0** | 2015 | M. Chen, J. Xia, C. Huang, J. M. Dieterich, L. Hung, I. Shin, E. A. Carter | Added *ab initio* molecular dynamics capability (OF-DFT-driven MD), thermostats, and further algorithmic/functional advances, positioning PROFESS as a tool for dynamical, finite-temperature simulations of metals (e.g., liquid tin, liquid lithium melting studies). |
| **PROFESS@QE (coupled)** | 2014 | V. V. Karasiev, T. Sjostrom, S. B. Trickey | Not a PROFESS version per se, but a significant fork/coupling: implementation of finite-temperature (free-energy) OF-DFT functionals (finite-T Thomas–Fermi, second-order gradient, finite-T GGA) within PROFESS, coupled to Quantum Espresso so that OF-DFT and KS-DFT MD could be run with identical algorithms/thermostats/convergence settings — used extensively for warm-dense-matter and equation-of-state studies. |
| **libKEDF** | 2017 | J. M. Dieterich, W. C. Witt, E. A. Carter | An accelerated, modular C++ library of nonlocal KEDFs, developed as a stepping-stone toward re-architecting PROFESS; demonstrated the benefits of moving off the original Fortran codebase. |
| **PROFESS 4 (`profess`)** | ongoing | Broader open-source community (maintainer W. C. Witt et al.) | Full rewrite: a C++ core library, a `profess` executable for simple calculations, and Python bindings for complex workflows, prioritizing interoperability with other materials-modeling tools. The project is now developed independently of any single institution, hosted with documentation at profess.dev (source at the `EACcodes/PROFESS` GitHub organization, which also hosts the legacy Fortran PROFESS 3 code). |
| **PROFESS-AD** | 2023 | (Carter-group-adjacent development) | A prototyping/research variant built on automatic differentiation, enabling direct evaluation of higher-order derivatives (bulk moduli, elastic constants, force constants) and facilitating rapid development/testing of new KEDFs, with features intended to migrate into the mainstream performance-optimized PROFESS code. |

---

## 5. Notable Applications and Landmark Studies

- **Million-atom bulk aluminum**: explicit quantum-mechanical (OF-DFT) treatment of over one million Al atoms, a scale essentially unreachable by conventional KS-DFT, illustrating PROFESS's core value proposition for mesoscale materials science.
- **Large lattice defects and nanostructures**: dislocations, grain boundaries, and large nanostructures (nanowires, quantum dots) — systems too large for standard first-principles methods but too subtle for empirical potentials — are natural targets for OF-DFT/PROFESS.
- **Melting and liquid-metal dynamics**: molecular dynamics studies of liquid lithium melting and liquid tin dynamics using PROFESS 3.0's OF-DFT MD.
- **Warm dense matter / equation-of-state studies**: PROFESS coupled with Quantum Espresso has been used for finite-temperature simulations up to millions of Kelvin, and for constructing wide-range equation-of-state tables (e.g., for CHON compounds relevant to inertial confinement fusion), combining KS-DFT, OF-DFT, and extrapolation.
- **Semiconductor and alloy properties**: with the Huang–Carter and related KEDFs, PROFESS has been used to compute bulk moduli, equilibrium volumes, phase-transition pressures, point-defect energies, and surface energies for Si and other semiconductors, as well as Al–Mg intermetallics and Li–Si amorphous alloys.

---

## 6. Relationship to Other OF-DFT Codes

PROFESS is one of several mature, actively used OF-DFT packages, alongside:

- **ATLAS** — a real-space finite-difference OF-DFT implementation (demonstrated single-point energy calculations for over 100 million aluminum atoms).
- **DFTpy** — an object-oriented, Python-based OF-DFT platform.
- **CONUNDrum** and **DRAGON** — additional OF-DFT codes.
- **GPAW** — primarily a projector-augmented-wave KS-DFT code that also includes an orbital-free DFT module.
- **PROFESS-AD** — the automatic-differentiation-based prototyping tool discussed above, designed to feed new methods back into mainstream PROFESS.
- Some general-purpose KS-DFT codes (e.g., **ABINIT**, **ABACUS**) have also incorporated OF-DFT capabilities.

`libKEDF` — the standalone accelerated KEDF library spun out of PROFESS development — has itself been used as an interoperable nonlocal-KEDF engine referenced by other OF-DFT efforts.

---

## 7. Distribution and Access

- **Legacy Fortran PROFESS 3**: available via the `EACcodes/PROFESS` GitHub repository, and historically distributed through the Computer Physics Communications (CPC) Program Library (Catalogue ID AEBN_v1_0 for the original 2008 release), requiring FFTW and MINPACK-2 as external dependencies; designed for Linux with Intel/AMD compilers.
- **Modern PROFESS 4 (`profess`)**: a C++ core library plus a `profess` executable and Python bindings, documented at profess.dev, developed as an open-source, cross-institutional project.

---

## 8. Summary Assessment

PROFESS occupies a distinctive niche in first-principles materials simulation: by replacing Kohn–Sham orbitals with an explicit kinetic energy density functional, it achieves scaling and system sizes far beyond conventional DFT, at the cost of accuracy that is fundamentally tied to the quality of the KEDF and the availability of good local pseudopotentials. It is best suited to — and has been most successfully applied to — **large-scale simulations of main-group, nearly-free-electron metals and their alloys** (Al, Mg, Li, Na, Si and related systems), including defects, nanostructures, liquids, and warm-dense-matter states, where million-atom or high-temperature regimes place these problems beyond the reach of orbital-based methods. Its two-decade development history — from Watson & Carter's original linear-scaling algorithms, through the Ho/Lignères/Carter, Hung/Huang/Shin, and Chen/Xia/Dieterich-era Fortran releases, to the current C++/Python `profess` rewrite and the automatic-differentiation PROFESS-AD prototype — mirrors the broader maturation of OF-DFT and nonlocal KEDF theory (Wang–Teter → Wang–Govind–Carter → Huang–Carter and successors) as a distinct, complementary branch of density functional methodology alongside Kohn–Sham DFT.

---

## 9. Bibliography — Publications on PROFESS Theory and Development

### 9.1 Core PROFESS code papers

1. S. C. Watson, E. A. Carter, "Linear-scaling parallel algorithms for the first principles treatment of metals," *Computer Physics Communications* **128**, 67 (2000). DOI: 10.1016/S0010-4655(00)00064-3
2. G. S. Ho, V. L. Lignères, E. A. Carter, "Introducing PROFESS: A new program for orbital-free density functional theory calculations," *Computer Physics Communications* **179**, 839–854 (2008). DOI: 10.1016/j.cpc.2008.07.002
3. L. Hung, C. Huang, I. Shin, G. S. Ho, V. L. Lignères, E. A. Carter, "Introducing PROFESS 2.0: A parallelized, fully linear scaling program for orbital-free density functional theory calculations," *Computer Physics Communications* **181**, 2208–2209 (2010). DOI: 10.1016/j.cpc.2010.09.001
4. L. Hung, E. A. Carter, "Accurate simulations of metals at the mesoscale: Explicit treatment of 1 million atoms with quantum mechanics," *Chemical Physics Letters* **475**, 163–170 (2009). DOI: 10.1016/j.cplett.2009.04.059
5. M. Chen, J. Xia, C. Huang, J. M. Dieterich, L. Hung, I. Shin, E. A. Carter, "Introducing PROFESS 3.0: An advanced program for orbital-free density functional theory molecular dynamics simulations," *Computer Physics Communications* **190**, 228–230 (2015). DOI: 10.1016/j.cpc.2014.12.021
6. V. V. Karasiev, T. Sjostrom, S. B. Trickey, "Finite-temperature orbital-free DFT molecular dynamics: Coupling PROFESS and Quantum Espresso," *Computer Physics Communications* **185**, 3240–3249 (2014). DOI: 10.1016/j.cpc.2014.06.001; arXiv:1406.0835
7. J. M. Dieterich, W. C. Witt, E. A. Carter, "libKEDF: An accelerated library of kinetic energy density functionals," *Journal of Computational Chemistry* **38**, 1552–1560 (2017). DOI: 10.1002/jcc.24806

### 9.2 Kinetic energy density functional (KEDF) theory used in PROFESS

8. L.-W. Wang, M. P. Teter, "Kinetic-energy functional of the electron density," *Physical Review B* **45**, 13196 (1992). [Wang–Teter functional]
9. Y. A. Wang, N. Govind, E. A. Carter, "Orbital-free kinetic-energy functionals for the nearly free electron gas," *Physical Review B* **58**, 13465 (1998).
10. Y. A. Wang, N. Govind, E. A. Carter, "Orbital-free kinetic-energy density functionals with a density-dependent kernel," *Physical Review B* **60**, 16350 (1999); erratum *Physical Review B* **64**, 129901 (2001). [Wang–Govind–Carter (WGC) functional]
11. B. Zhou, Y. A. Wang, E. A. Carter, "Failure of the Wang–Teter and related approximate kinetic energy density functionals for silicon," *Physical Review B* **69**, 125109 (2004).
12. C. Huang, E. A. Carter, "Nonlocal orbital-free kinetic energy density functional for semiconductors," *Physical Review B* **81**, 045206 (2010). DOI: 10.1103/PhysRevB.81.045206 [Huang–Carter (HC) functional]
13. I. Shin, E. A. Carter, "Enhanced von Weizsäcker Wang-Govind-Carter kinetic energy density functional for semiconductors," *Journal of Chemical Physics* **140**, 18A531 (2014). DOI: 10.1063/1.4869867
14. W. Mi, A. Genova, M. Pavanello, "Nonlocal kinetic energy functionals by functional integration," *Journal of Chemical Physics* **148**, 184107 (2018). DOI: 10.1063/1.5023926
15. J. Xia, E. A. Carter, "Single-point kinetic energy density functionals: A pointwise kinetic energy density analysis and numerical convergence investigation," *Physical Review B* **91**, 045124 (2015).
16. W. C. Witt, E. A. Carter, "Kinetic energy density of nearly free electrons. I. Response functionals of the external potential," *Physical Review B* **100**, 122106 (2019).
17. W. C. Witt, E. A. Carter, "Kinetic energy density of nearly free electrons. II. Response functionals of the electron density," *Physical Review B* **100**, 125107 (2019).
18. W. C. Witt, K. Jiang, E. A. Carter, "Upper bound to the gradient-based kinetic energy density of noninteracting electrons in an external potential," *Journal of Chemical Physics* **151**, 064113 (2019).
19. (Derivation of the analytic WGC kernel form used computationally in PROFESS) — "Analytic form for a nonlocal kinetic energy functional with a density-dependent kernel for orbital-free density functional theory under periodic and Dirichlet boundary conditions," *Physical Review B* **78**, 045105 (2008). DOI: 10.1103/PhysRevB.78.045105

### 9.3 Local pseudopotentials for OF-DFT

20. J.-D. Chai, J. D. Weeks, "Orbital-free density functional theory: Kinetic potentials and ab initio local pseudopotentials," *Physical Review B* **75**, 205122 (2007).
21. J.-D. Chai, V. L. Lignères, G. Ho, E. A. Carter, J. D. Weeks, "Orbital-free density functional theory: Linear scaling methods for kinetic potentials, and applications to solid Al and Si," *Chemical Physics Letters* **473**, 263–267 (2009).
22. C. Huang, E. A. Carter, "Transferable local pseudopotentials for magnesium, aluminum and silicon," *Physical Chemistry Chemical Physics* **10**, 7109–7120 (2008). [Bulk-derived local pseudopotentials, BLPS]

### 9.4 Reviews, methodological surveys, and PROFESS-adjacent studies

23. Y. A. Wang, E. A. Carter, "Orbital-free kinetic-energy density functional theory," in *Theoretical Methods in Condensed Phase Chemistry*, S. D. Schwartz (Ed.), Kluwer, pp. 117–184 (2002).
24. V. L. Lignères, E. A. Carter, "An introduction to orbital-free density functional theory," in *Handbook of Materials Modeling: Methods*, S. Yip (Ed.), Springer, pp. 137–148 (2005). DOI: 10.1007/978-1-4020-3286-8_9
25. V. V. Karasiev, S. B. Trickey, "Issues and challenges in orbital-free density functional calculations," *Computer Physics Communications* **183**, 2519–2527 (2012); arXiv:1109.6602
26. J. Xia, C. Huang, I. Shin, E. A. Carter, "Can orbital-free density functional theory simulate molecules?," *Journal of Chemical Physics* **136**, 084102 (2012). DOI: 10.1063/1.3685604
27. W. C. Witt, B. G. del Rio, J. M. Dieterich, E. A. Carter, "Orbital-free density functional theory for materials research," *Journal of Materials Research* **33**, 777–795 (2018). DOI: 10.1557/jmr.2017.462
28. W. Mi, K. Luo, S. B. Trickey, M. Pavanello, "Orbital-free density functional theory: An attractive electronic structure method for large-scale first-principles simulations," *Chemical Reviews* **123**, 12039–12104 (2023).
29. Q. Xu, C. Ma, W. Mi, Y. Wang, Y. Ma, "Recent advancements and challenges in orbital-free density functional theory," *WIREs Computational Molecular Science* **14**, e1724 (2024). DOI: 10.1002/wcms.1724
30. Y. Ke, F. Libisch, J. Xia, L.-W. Wang, E. A. Carter, "Angular-momentum-dependent orbital-free density functional theory," *Physical Review Letters* **111**, 066402 (2013).
31. Y. Ke, F. Libisch, J. Xia, E. A. Carter, "Angular momentum dependent orbital-free density functional theory: Formulation and implementation," *Physical Review B* **89**, 155112 (2014).

### 9.5 Related/competing codes and automatic-differentiation extension (cited for context)

32. W. Mi, X. Shao, C. Su, Y. Zhou, S. Zhang, Q. Li, H. Wang, L. Zhang, M. Miao, Y. Wang, et al., "ATLAS: A real-space finite-difference implementation of orbital-free density functional theory," *Computer Physics Communications* **200**, 87–95 (2016). arXiv:1507.07373
33. X. Shao, K. Jiang, W. Mi, A. Genova, M. Pavanello, "DFTpy: An efficient and object-oriented platform for orbital-free DFT simulations," *WIREs Computational Molecular Science* **11**, e1482 (2021).
34. (PROFESS-AD) "Automatic differentiation for orbital-free density functional theory," *Journal of Chemical Physics* **158**, 124801 (2023). DOI: 10.1063/5.0142294; arXiv:2212.03231
35. W. C. Witt, B. W. Shires, C. W. Tan, W. J. Jankowski, C. J. Pickard, "Random structure searching with orbital-free density functional theory," *Journal of Physical Chemistry A* **125**, 1650–1660 (2021).

---

*Compiled from the PROFESS project documentation (profess.dev), the Computer Physics Communications publication record, and the broader orbital-free DFT literature. Version and authorship details for the core code releases are drawn directly from the official PROFESS project history page.*

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt:  Create an exhaustive review of the PROFESS 	Orbital-free DFT code for large-scale simulations using kinetic-energy density functionals instead of Kohn-Sham orbitals. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
