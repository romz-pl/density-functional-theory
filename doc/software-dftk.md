# DFTK.jl (Density-Functional Toolkit): An Exhaustive Technical Review

## 1. Overview and Positioning

DFTK.jl (Density-Functional ToolKit) is a Julia library implementing plane-wave, pseudopotential Kohn–Sham density-functional theory (DFT) for periodic systems. It is developed under the JuliaMolSim GitHub organization, led primarily by Michael F. Herbst and Antoine Levitt, with substantial contributions from the applied-mathematics groups of Eric Cancès and collaborators (École des Ponts ParisTech, Inria Paris, Sorbonne Université) and, more recently, EPFL's Mathematics for Materials Modelling (MatMat) group.

Unlike production codes such as **Abinit**, **Quantum ESPRESSO**, or **VASP**, DFTK is explicitly *not* optimized primarily for raw throughput on the largest possible systems. Its stated design goal is to be a **research vehicle at the intersection of numerical analysis, applied mathematics, and electronic-structure physics** — a codebase small and clean enough (~10,000 lines of pure Julia, as of the current stable release) that novel algorithms can be prototyped, verified, and pushed to production-scale problems (systems up to roughly 1000 electrons) within the *same* code, rather than being confined to disposable single-purpose research scripts.

Consequently, DFTK's core value proposition is threefold:

1. **Algorithmic transparency** — every routine (SCF, eigensolvers, response theory, symmetry handling) is written in high-level, inspectable Julia, avoiding hidden Fortran/C++ black boxes.
2. **Mathematical rigor** — many features exist because they are the direct output of peer-reviewed numerical-analysis results, and the code is used as the reference implementation for those papers.
3. **Multiple dispatch and genericity** — DFTK exploits Julia's type system to make numerical precision (`Float32`, `Float64`, `Double64`, interval arithmetic), automatic differentiation, and GPU execution largely orthogonal, "bolt-on" concerns rather than separate code paths.

---

## 2. Mathematical and Physical Model

### 2.1 The periodic Kohn–Sham problem

DFTK solves the Kohn–Sham equations for a crystal with periodic lattice $\Omega$, discretized via Bloch's theorem into a family of $k$-point problems over the Brillouin zone. For each $k$, one seeks orbitals $\{\psi_{k,n}\}$ solving

$$
H_k(\rho) \, \psi_{k,n} = \varepsilon_{k,n}\, \psi_{k,n}, \qquad
H_k(\rho) = \tfrac{1}{2}(-i\nabla + k)^2 + V(\rho),
$$

with the density self-consistently defined as

$$
\rho(r) = \sum_k w_k \sum_n f(\varepsilon_{k,n}) \, |\psi_{k,n}(r)|^2 .
$$

The potential $V(\rho) = V_\text{nuc} + V_\text{H}^\rho + V_\text{xc}(\rho)$ combines:

- the **nuclear/pseudopotential attraction** term (local + nonlocal, Kleinman–Bylander projector form),
- the **Hartree potential**, obtained as the unique zero-mean solution of the periodic Poisson equation
  $$-\Delta V_\text{H}^\rho(r) = 4\pi\left(\rho(r) - \frac{1}{|\Omega|}\int_\Omega \rho\right),$$
- the **exchange–correlation potential** $V_\text{xc}(\rho)$, delegated to **libxc** for LDA/GGA/meta-GGA functionals.

The occupation function $f$ (Fermi–Dirac, Gaussian, Methfessel–Paxton, or Marzari–Vanderbilt cold smearing) handles finite-temperature and metallic occupation, chosen so $\int \rho = N_\text{electrons}$.

### 2.2 Discretization

Orbitals are expanded in plane waves truncated by a kinetic-energy cutoff $E_\text{cut}$, giving a finite-dimensional nonlinear eigenvalue problem per $k$-point. DFTK documents this discretization explicitly as a mathematical object (see the *Periodic problems and plane-wave discretizations* and *Comparing discretization techniques* chapters), treating basis-set truncation error, $k$-point sampling, and the choice between plane waves and other bases (finite differences, finite elements) as first-class, quantifiable sources of numerical error rather than incidental implementation details.

### 2.3 Beyond standard Kohn–Sham

Because the `Model`/`Hamiltonian` abstraction is generic over the terms composing $V(\rho)$, DFTK is used to implement and study equations well outside conventional solid-state DFT, including:

- **Cohen–Bergstresser** empirical pseudopotential (linear) models,
- the **Gross–Pitaevskii equation** (nonlinear Schrödinger equation for Bose–Einstein condensates), in 1D and with an external magnetic field,
- **anyonic** models (2D magnetic Thomas–Fermi-type problems),
- fully custom analytic potentials.

This generality is deliberate: it lets applied mathematicians validate numerical methods on *simplified* model problems using the exact same infrastructure that scales to *real* crystalline silicon or graphene calculations — directly supporting the "toy problem → production scale in one code" philosophy.

---

## 3. Algorithmic Content

### 3.1 Self-consistent field (SCF) algorithms

SCF convergence in metals and other difficult systems is one of DFTK's signature research areas.

- **Density/potential mixing** with Anderson/Pulay acceleration, and specifically **LDOS (local density of states) mixing**, a black-box preconditioner that auto-detects whether a system is metallic or insulating and adapts accordingly (Herbst & Levitt, *J. Phys.: Condens. Matter*, 2021).
- **Adaptive/robust line search (damping)** for SCF iterations, replacing heuristic fixed-damping schemes with a rigorously justified adaptive procedure (Herbst & Levitt, *J. Comput. Phys.*, 2022).
- **Direct minimization** methods (as an alternative to fixed-point SCF), with convergence properties analyzed rigorously by Cancès, Kemlin & Levitt (*SIAM J. Matrix Anal. Appl.*, 2021), including a Rayleigh–Ritz correction step.
- Multiple solver back-ends are directly comparable within the same framework (`Comparison of DFT solvers`, `Analysing SCF convergence` examples), which is itself a research tool: it lets users empirically benchmark preconditioners and mixing schemes on identical Hamiltonians.

### 3.2 Eigensolvers

Diagonalization of $H_k(\rho)$ at fixed density is performed with iterative Krylov-subspace eigensolvers (via `KrylovKit.jl`), including locally optimal block preconditioned conjugate gradient (LOBPCG)-style approaches typical of large-scale plane-wave DFT, with dynamic band-count adaptation to guarantee occupation of the lowest-lying orbitals up to a specified occupation threshold.

### 3.3 Density-functional perturbation theory (DFPT) and response

DFTK implements a full **DFPT** stack for linear-response properties (polarizability, elastic constants, phonons), and — notably — couples DFPT with **forward-mode algorithmic differentiation (AD)** to build "AD-DFPT": derivatives of *any* output quantity with respect to *any* input parameter (atomic geometry, functional parameters, pseudopotential parameters) can be obtained without hand-deriving adjoint/perturbation equations (Schmitz, Ploumhans & Herbst, *npj Comput. Mater.*, 2026). This is used for inverse materials design (e.g., band-gap targeting), learning of exchange-correlation functional parameters, and uncertainty propagation from DFT inputs to relaxed structures — an explicitly *algorithmic-research* application enabled by the code's genericity.

Numerical stability of response-property calculations in *metals* (where perturbation theory is complicated by the absence of a spectral gap) was separately analyzed and implemented following Cancès, Herbst, Kemlin, Levitt & Stamm (*Lett. Math. Phys.*, 2023). Efficient Krylov-subspace strategies specifically for the linear-response (Sternheimer-type) equations have also been developed and integrated (Herbst & Sun, 2025).

### 3.4 Band-structure and spectral algorithms

The **modified-operator method** for computing band diagrams of crystalline materials — a technique designed to remove spurious features/improve conditioning near band crossings — was developed and validated using DFTK (Cancès, Hassan & Vidal, *Math. Comp.*, 2024).

### 3.5 Symmetry and Brillouin-zone sampling

Brillouin-zone symmetry reduction for $k$-point sampling is delegated to **spglib**, with an internal search for irreducible $k$-point subsets, reducing the effective computational cost of Monkhorst–Pack-type grids.

### 3.6 Parallelism and hardware

- **MPI-based distributed parallelism** over $k$-points.
- **GPU support**: NVIDIA GPUs are considered largely mature; AMD GPU support is preliminary.
- **Multi-level threading** across $k$-points, eigenvector blocks, FFTs, and linear algebra.

The separation of the numerical algorithm from the hardware back-end is again a consequence of Julia's array abstractions (broadcasting, generic linear algebra) rather than hand-written CUDA/HIP kernels, which keeps the "one algorithm, many back-ends" property intact.

### 3.7 Genericity as an algorithmic tool

A structurally important, and somewhat unusual, feature is **arbitrary floating-point-type support**: the entire SCF/eigensolver/DFPT stack is written generically enough to run in `Float32`, `Float64`, or extended precision (`Double64` via DoubleFloats.jl). This is not a convenience feature — it is used as a *research instrument* for:

- rigorously separating **discretization error** from **floating-point/round-off error** in convergence studies,
- validating error-estimator theory (Section 4) against near-exact reference solutions computed in extended precision,
- studying low-precision (`Float32`) feasibility for large-scale or GPU-bound calculations.

---

## 4. Verification, Error Control, and Mathematical Rigor

This is arguably DFTK's most distinctive contribution relative to conventional plane-wave codes, and the area with the deepest publication record.

### 4.1 A posteriori and a priori error estimation

A substantial and coherent line of work, largely emerging from the Cancès/Levitt/Dusson/Kemlin collaboration and implemented directly in DFTK, targets **rigorous, computable error bounds** on DFT observables:

- **A posteriori error estimation** for the *non-self-consistent* Kohn–Sham equations (Herbst, Levitt & Cancès, *Faraday Discussions*, 2020) — the first step toward quantifying how discretization error (basis-set truncation, $k$-point sampling) propagates into physical observables.
- **Practical error bounds for properties** (notably atomic forces) in plane-wave calculations (Cancès, Dusson, Kemlin & Levitt, *SIAM J. Sci. Comput.*, 2022), with a worked example distributed as `Practical error bounds for the forces` in the documentation itself — i.e., the theory is directly executable, not just cited.
- **Fully guaranteed and computable error bounds on the total energy** for periodic Kohn–Sham equations with convex density functionals (Bordignon, Dusson, Cancès, Kemlin, Reyes & Stamm, 2024) — extending the guarantees from *properties* to the *energy* itself, under precise convexity assumptions on the functional.
- **A priori error analysis** of linear and nonlinear periodic Schrödinger equations with analytic potentials (Cancès, Kemlin & Levitt, *J. Sci. Comput.*, 2024), establishing convergence rates for the plane-wave discretization under explicit smoothness assumptions — theoretical underpinning for *why* plane-wave cutoffs converge the way they empirically do.

Together these results move DFT error control from "increase $E_\text{cut}$ until the answer looks converged" toward certified, quantitative bounds — a research program essentially unique to DFTK among actively maintained plane-wave codes.

### 4.2 Spectral-theoretic analysis

- **Feshbach–Schur method** analysis for Fourier spectral discretizations of Schrödinger operators (Dusson, Sigal & Stamm, *Math. Comp.*, 2023), providing rigorous convergence theory for eigenvalue problems of the type solved at each SCF iteration.
- **Kohn–Sham inversion with mathematical guarantees** (Herbst, Bakkestuen & Laestadius, *Phys. Rev. B*, 2025) — placing the (numerically delicate) inverse-DFT problem (recovering a potential from a target density) on rigorous footing, using DFTK as the computational testbed.

### 4.3 Empirical cross-validation against production codes

Beyond internal mathematical guarantees, DFTK is benchmarked empirically:

- The documentation explicitly references the **Bosoni et al. verification study** (acwf-verification.materialscloud.org), an independent large-scale cross-code comparison of equations of state for unary compounds and oxides across the periodic table, and states that DFTK agrees closely with the other participating codes.
- Early releases explicitly reported **close agreement with Abinit** on thoroughly tested reference systems (silicon, graphite, manganese) as a baseline correctness check before broader adoption.
- The library maintains an extensive **automated unit-test suite** (documented under *Developer resources → Unit test system*) exercised via continuous integration on every change, covering physical correctness (total energies, forces, band structures against reference data) as well as internal numerical consistency (SCF convergence, symmetry handling, response-property consistency).

### 4.4 Verification-by-construction: examples as regression tests

Several research papers ship **companion computational scripts** that are either included directly in DFTK's example suite or hosted in linked repositories (see Section 6 and the publication list), meaning the numerical results reported in the mathematical literature are directly reproducible against the live codebase — an unusually strong form of open, ongoing verification, since regressions in the numerics would surface as reproducibility failures against these published benchmarks.

### 4.5 Algorithmic differentiation as a verification tool

The AD-DFPT framework (Section 3.3) is also explicitly framed as an **error-control** device: because forward-mode AD gives exact derivatives of DFT outputs with respect to any input, it enables direct propagation of parameter uncertainty (e.g., pseudopotential or functional parameters) into uncertainty on relaxed structures and computed properties — turning "error control" from a purely discretization-focused concept into one that also spans *model*-parameter uncertainty.

---

## 5. Software Engineering and Reproducibility Infrastructure

- **Language and size**: pure Julia, ≈10,000 lines (growing from ≈5k in 2019–2020 to ≈7k in 2021–2022 to ≈10k currently), reflecting steady, incremental growth rather than a large monolithic rewrite.
- **Dependencies of note**: `libxc` (exchange-correlation functionals), `FFTW` (Fourier transforms), `spglib` (crystal symmetry), `Optim.jl` (optimization algorithms), `KrylovKit.jl` (Krylov eigen/linear solvers) — each independently citable, and DFTK's documentation explicitly asks users to cite the *dependency* papers alongside DFTK itself, an unusually careful attribution practice.
- **Interoperability**: AtomsBase.jl and AtomsCalculators.jl integration (Julia materials-science ecosystem interop), ASE (Atomistic Simulation Environment) support via `asedftk`, ETSF Nanoquanta I/O, and Wannierization via Wannier.jl or Wannier90.
- **Precision and hardware genericity**: `Float32`/`Float64`/`Double64` support and GPU execution are implemented as generic code paths rather than special-cased branches — a direct benefit of Julia's multiple dispatch and generic array programming model.
- **Archival and versioning**: Zenodo-archived releases (DOI `10.5281/zenodo.3541724` covers all versions) provide long-term citable, versioned snapshots of the code, supporting exact reproducibility of results tied to a specific release.
- **Governance**: developed under the JuliaMolSim GitHub organization; contributions and issues are openly tracked; a dedicated *Developer's style guide*, *Notation and conventions*, and *Data structures* documentation aim to keep the mathematically dense codebase approachable to new contributors.
- **Funding**: ISCD (Sorbonne Université), École des Ponts ParisTech, Inria Paris, RWTH Aachen University, the Swiss National Science Foundation, the European Research Council (Horizon 2020, grant agreement No. 810367), and the Simons Foundation — a funding profile dominated by mathematics/numerical-analysis and basic-research grants rather than industrial/production computing grants, consistent with its stated mission.

---

## 6. Illustrative List of Research Enabled by / Conducted Through DFTK

Beyond the algorithm papers in Section 4, DFTK has served as the computational backbone for a range of independent applied-mathematics and physics research, including (non-exhaustively):

- Low-temperature behavior of DFT for metals via Sommerfeld expansion and DFPT (Gonze, Tantardini & Levitt, *Phys. Rev. B*, 2026).
- Stochastic DFT via multilevel Monte Carlo (Quan & Chen, 2025).
- Neural-network acceleration of iterative solvers for nonlinear Schrödinger eigenvalue problems (Petersheim, Pietschmann, Püschel & Ruess, 2025).
- Magnetic Thomas–Fermi theory for 2D abelian anyons (Levitt, Lundholm & Rougerie, 2025).
- Energy-adaptive Riemannian conjugate-gradient methods for DFT (Petersheim, Püschel & Stykel, 2025).
- Dirac cones in mean-field models of graphene (Cazalis, *Pure Appl. Anal.*, 2024).
- Continuous moiré-scale models for twisted bilayer graphene (Cancès, Garrigue & Gontier, *Phys. Rev. B*, 2023).

This body of work illustrates DFTK's role less as an end-user production DFT package and more as **shared infrastructure for the numerical-analysis-of-electronic-structure community** — a niche few other actively maintained DFT codes occupy.

---

## 7. Strengths and Limitations

**Strengths**

- Tight, verifiable link between published numerical-analysis theorems and running code — theory and implementation co-evolve.
- Genericity (precision, hardware, differentiability) obtained largely "for free" via Julia's type system, rather than through duplicated code paths.
- Transparent, inspectable algorithms lower the barrier for methodological innovation (new preconditioners, solvers, error estimators) relative to legacy Fortran/C++ codebases.
- Strong emphasis on reproducibility: papers ship companion scripts runnable against the live library; Zenodo-archived, DOI-citable releases.
- Rigorous, *guaranteed* (not just empirical) error bounds are available for a growing class of problems (currently requiring, e.g., convexity of the functional) — a capability essentially absent from mainstream production codes.

**Limitations**

- Explicitly scoped to systems of up to roughly 1000 electrons; not intended to compete with highly optimized production codes (VASP, Quantum ESPRESSO, Abinit) at the largest system sizes.
- Some advanced features are marked **preliminary**: exact exchange/hybrid DFT, phonon computations, and AMD GPU support.
- As a research-first codebase, feature completeness (e.g., breadth of exchange-correlation approximations beyond what libxc offers, exotic pseudopotential formats) trails the most mature production packages.
- Rigorous error bounds on the *energy* currently require convexity assumptions on the density functional, limiting immediate applicability to some standard GGA/meta-GGA functionals without further theoretical extension.

---

## 8. Publications Related to DFTK's Theory

### 8.1 Reference paper (cite this for the code itself)

- M. F. Herbst, A. Levitt, E. Cancès. **"DFTK: A Julian approach for simulating electrons in solids."** *Proc. JuliaCon Conf.* **3**, 69 (2021). DOI: [10.21105/jcon.00069](https://doi.org/10.21105/jcon.00069)

### 8.2 Papers describing DFTK's algorithms

- N. F. Schmitz, B. Ploumhans, M. F. Herbst. **"Algorithmic differentiation for plane-wave DFT: materials design, error control and learning model parameters."** *npj Computational Materials* **12**, 6 (2026). DOI: [10.1038/s41524-025-01880-3](https://doi.org/10.1038/s41524-025-01880-3)
- M. F. Herbst, B. Sun. **"Efficient Krylov methods for linear response in plane-wave electronic structure calculations."** (2025). [ArXiv:2505.02319](https://arxiv.org/abs/2505.02319)
- E. Cancès, M. Hassan, L. Vidal. **"Modified-Operator Method for the Calculation of Band Diagrams of Crystalline Materials."** *Math. Comp.* **93**, 1203 (2024). DOI: [10.1090/mcom/3897](https://doi.org/10.1090/mcom/3897); [ArXiv:2210.00442](https://arxiv.org/abs/2210.00442)
- E. Cancès, M. F. Herbst, G. Kemlin, A. Levitt, B. Stamm. **"Numerical stability and efficiency of response property calculations in density functional theory."** *Letters in Mathematical Physics* **113**, 21 (2023). [ArXiv:2210.04512](https://arxiv.org/abs/2210.04512)
- M. F. Herbst, A. Levitt. **"A robust and efficient line search for self-consistent field iterations."** *Journal of Computational Physics* **459**, 111127 (2022). [ArXiv:2109.14018](https://arxiv.org/abs/2109.14018)
- M. F. Herbst, A. Levitt. **"Black-box inhomogeneous preconditioning for self-consistent field iterations in density functional theory."** *Journal of Physics: Condensed Matter* **33**, 085503 (2021). DOI: [10.1088/1361-648X/abcbdb](https://doi.org/10.1088/1361-648X/abcbdb); [ArXiv:2009.01665](https://arxiv.org/abs/2009.01665)

### 8.3 Verification, error control, and mathematical-analysis papers (research conducted with DFTK)

- X. Gonze, C. Tantardini, A. Levitt. **"Low-temperature behavior of density-functional theory for metals based on density-functional perturbation theory and Sommerfeld expansion."** *Physical Review B* **113**, 035125 (2026). DOI: [10.1103/yj83-j9p1](https://doi.org/10.1103/yj83-j9p1)
- X. Quan, H. Chen. **"Stochastic Density Functional Theory Through the Lens of Multilevel Monte Carlo Method."** (2025). [ArXiv:2512.04860](https://arxiv.org/abs/2512.04860v2)
- D. Petersheim, J.-F. Pietschmann, J. Püschel, K. Ruess. **"Neural Network Acceleration of Iterative Methods for Nonlinear Schrödinger Eigenvalue Problems."** (2025). [ArXiv:2507.16349](https://arxiv.org/abs/2507.16349)
- A. Levitt, D. Lundholm, N. Rougerie. **"Magnetic Thomas-Fermi theory for 2D abelian anyons."** (2025). [ArXiv:2504.13481](https://arxiv.org/abs/2504.13481)
- D. Petersheim, J. Püschel, T. Stykel. **"Energy-adaptive Riemannian Conjugate Gradient method for density functional theory."** (2025). [ArXiv:2503.16225](https://arxiv.org/abs/2503.16225)
- A. Bordignon, G. Dusson, E. Cancès, G. Kemlin, R. A. L. Reyes, B. Stamm. **"Fully guaranteed and computable error bounds on the energy for periodic Kohn-Sham equations with convex density functionals."** (2024). [ArXiv:2409.11769](https://arxiv.org/abs/2409.11769)
- M. F. Herbst, V. H. Bakkestuen, A. Laestadius. **"Kohn-Sham inversion with mathematical guarantees."** *Phys. Rev. B* **111**, 205143 (2025). DOI: [10.1103/PhysRevB.111.205143](https://doi.org/10.1103/PhysRevB.111.205143); [ArXiv:2409.04372](https://arxiv.org/abs/2409.04372)
- J. Cazalis. **"Dirac cones for a mean-field model of graphene."** *Pure and Applied Analysis* **6**, 129 (2024). DOI: [10.2140/paa.2024.6.129](https://doi.org/10.2140/paa.2024.6.129); [ArXiv:2207.09893](https://arxiv.org/abs/2207.09893)
- E. Cancès, G. Kemlin, A. Levitt. **"A Priori Error Analysis of Linear and Nonlinear Periodic Schrödinger Equations with Analytic Potentials."** *J. Sci. Comp.* **98**, 25 (2024). DOI: [10.1007/s10915-023-02421-0](https://doi.org/10.1007/s10915-023-02421-0); [ArXiv:2206.04954](https://arxiv.org/abs/2206.04954)
- E. Cancès, L. Garrigue, D. Gontier. **"A simple derivation of moiré-scale continuous models for twisted bilayer graphene."** *Physical Review B* **107**, 155403 (2023). DOI: [10.1103/PhysRevB.107.155403](https://doi.org/10.1103/PhysRevB.107.155403); [ArXiv:2206.05685](https://arxiv.org/abs/2206.05685)
- G. Dusson, I. Sigal, B. Stamm. **"Analysis of the Feshbach-Schur method for the Fourier spectral discretizations of Schrödinger operators."** *Mathematics of Computation* **92**, 217 (2023). DOI: [10.1090/mcom/3774](https://doi.org/10.1090/mcom/3774); [ArXiv:2008.10871](https://arxiv.org/abs/2008.10871)
- E. Cancès, G. Dusson, G. Kemlin, A. Levitt. **"Practical error bounds for properties in plane-wave electronic structure calculations."** *SIAM Journal on Scientific Computing* **44**, B1312 (2022). DOI: [10.1137/21M1456224](https://doi.org/10.1137/21M1456224); [ArXiv:2111.01470](https://arxiv.org/abs/2111.01470)
- E. Cancès, G. Kemlin, A. Levitt. **"Convergence analysis of direct minimization and self-consistent iterations."** *SIAM Journal on Matrix Analysis and Applications* **42**, 243 (2021). DOI: [10.1137/20M1332864](https://doi.org/10.1137/20M1332864); [ArXiv:2004.09088](https://arxiv.org/abs/2004.09088)
- M. F. Herbst, A. Levitt, E. Cancès. **"A posteriori error estimation for the non-self-consistent Kohn-Sham equations."** *Faraday Discussions* **224**, 227 (2020). DOI: [10.1039/D0FD00048E](https://doi.org/10.1039/D0FD00048E); [ArXiv:2004.13549](https://arxiv.org/abs/2004.13549)

### 8.4 Key dependency papers (cite alongside DFTK when relevant)

- S. Lehtola, C. Steigemann, M. J. T. Oliveira, M. A. L. Marques. **"Recent developments in libxc — A comprehensive library of functionals for density functional theory."** *SoftwareX* **7**, 1 (2018). DOI: [10.1016/j.softx.2017.11.002](https://doi.org/10.1016/j.softx.2017.11.002)
- M. Frigo, S. G. Johnson. **"The design and implementation of FFTW3."** *Proceedings of the IEEE* **93**, 216 (2005). DOI: [10.1109/JPROC.2004.840301](https://doi.org/10.1109/JPROC.2004.840301)
- A. Togo, I. Tanaka. **"Spglib: a software library for crystal symmetry search."** [ArXiv:1808.01590](https://arxiv.org/abs/1808.01590) (2018).
- P. K. Mogensen, A. N. Riseth. **"Optim: A mathematical optimization package for Julia."** *Journal of Open Source Software* **3**, 615 (2018). DOI: [10.21105/joss.00615](https://doi.org/10.21105/joss.00615)

---

*Sources: DFTK.jl official documentation (docs.dftk.org / JuliaMolSim.github.io/DFTK.jl), including the Publications and Research-conducted-with-DFTK pages, DFTK Features page, and DFTK GitHub repository (JuliaMolSim/DFTK.jl), all accessed August 2026.*

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the DFTK.jl (Density-Functional Toolkit) 	Julia-based, plane-wave pseudopotential DFT code emphasizing algorithmic research, verification, and mathematical rigor. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
