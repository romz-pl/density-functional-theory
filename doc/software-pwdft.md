# PWDFT and PWDFT.jl: An Exhaustive Review of Plane-Wave DFT Codes for Large-Scale and GPU-Accelerated Electronic Structure Calculations

> **Note on naming.** "PWDFT" is not a single project — it is a name independently adopted by at least **three unrelated codebases**, each with distinct authorship, language, scope, and hardware targets. This review treats them separately to avoid the confusion that pervades much of the secondary literature:
>
> 1. **PWDFT.jl** — a pedagogical/research Julia package by Fadjar Fathurrahman and coworkers (Bandung Institute of Technology).
> 2. **PWDFT (C++, NWChemEx)** — Eric Bylaska's GPU-accelerated plane-wave engine, the modernized descendant of NWChem's PSPW/Band modules, developed for the NWChemEx project.
> 3. **PWDFT (C/C++, LBNL)** — the plane-wave DFT sub-module of the DGDFT package developed by the Lin Lin / Wei Hu / Chao Yang group, published in *Computer Physics Communications*, and extended by **PWDFT-SW** for the Sunway supercomputer.
>
> All three are reviewed below, with primary emphasis on (1) and (2) as requested, and (3) included because it is frequently cross-cited and shares the "PWDFT" name and large-scale/GPU orientation.

---

## 1. PWDFT.jl (Julia)

### 1.1 Overview

PWDFT.jl is a plane-wave, pseudopotential Kohn–Sham DFT package written in pure Julia. It was developed as both a research tool and a transparent, pedagogically accessible alternative to established Fortran/C++ plane-wave codes (Quantum ESPRESSO, ABINIT, VASP), aiming to demonstrate that a high-level, JIT-compiled language could deliver competitive performance for electronic-structure workloads while remaining easy to read, extend, and modify.

- **Repository:** `f-fathurrahman/PWDFT.jl` on GitHub
- **Primary author:** Fadjar Fathurrahman, with Mohammad Kemal Agusta, Adhitya Gandaryus Saputro, and Hermawan Kresno Dipojono (Institut Teknologi Bandung)
- **Language:** Julia (≥ 0.7 originally; actively tracks modern Julia)
- **License/registration:** Openly available on GitHub; historically *not* registered in Julia's General package registry, so installation requires cloning/`develop`-ing the package rather than `Pkg.add("PWDFT")`
- **License for use:** MIT-style permissive (per repository); check the repository's `LICENSE`/README for current terms

### 1.2 Design and Implementation

PWDFT.jl decomposes a typical KS-DFT calculation into three explicit stages, each backed by transparent, user-inspectable Julia data types:

1. **Structure initialization** — an `Atoms` type describing a molecular or crystalline structure (Cartesian or fractional coordinates, lattice vectors, species, periodic images).
2. **Hamiltonian construction** — a `Hamiltonian` type assembling the kinetic operator, local/non-local pseudopotential terms, Hartree potential, and exchange–correlation potential on a plane-wave/FFT grid.
3. **Self-consistent solution** — SCF or direct energy-minimization drivers that iterate the Hamiltonian to self-consistency.

Core low-level operations — wavefunction orthogonalization, application of kinetic and potential operators, and iterative diagonalization — are implemented in pure Julia rather than wrapping external Fortran/C kernels, which the authors present as a demonstration of Julia's suitability for numerically intensive scientific software without sacrificing code clarity.

**Basis set and grid.** Kohn–Sham orbitals are expanded in a plane-wave basis defined by a kinetic-energy cutoff (`ecutwfc`, in Hartree), the same basis-set philosophy used by Quantum ESPRESSO, ABINIT, and VASP. Real-space and reciprocal-space grids are connected via FFTs (via `FFTW.jl`).

**Pseudopotentials.** PWDFT.jl uses **Goedecker–Teter–Hutter (GTH)** norm-conserving, separable, analytic pseudopotentials, which are bundled directly in the repository. GTH pseudopotentials were chosen because their fully analytic real/reciprocal-space forms simplify implementation and are also used by reference codes such as ABINIT and CP2K/Quickstep, easing benchmarking.

**Exchange–correlation functionals.** Only two XC options are supported natively at present:
- `xcfunc="VWN"` — LDA (Vosko–Wilk–Nusair), the default
- `xcfunc="PBE"` — GGA (Perdew–Burke–Ernzerhof)

Both are accessed through **Libxc.jl**, a Julia wrapper around the Libxc exchange-correlation library, rather than being hand-coded, which keeps the door open for extending to the full Libxc functional set in the future.

**Solvers for the Kohn–Sham problem.** Several algorithms are implemented:
- `KS_solve_SCF!` — standard SCF with electron-density mixing (for semiconducting and metallic systems)
- `KS_solve_SCF_potmix!` — SCF with Hartree/XC *potential* mixing rather than density mixing
- `KS_solve_Emin_PCG!` — direct total-energy minimization via preconditioned conjugate gradient (following the approach of Arias and coworkers), applicable to systems with a band gap

Mixing schemes available for SCF include: simple (linear) mixing, linear-adaptive mixing, Anderson mixing, Broyden mixing, Pulay (DIIS-like) mixing, periodic Pulay mixing ("ppulay"), and restarted Pulay mixing ("rpulay").

**Metallic/finite-temperature systems.** Fermi–Dirac smearing of orbital occupations is supported via `use_smearing=true` with a smearing temperature `kT` (Hartree, default 0.001), combined with `extra_states` to include unoccupied bands near the Fermi level. As of the documented feature set, spin-polarized calculations (`Nspin=2`) are only supported in combination with smearing (no fixed-magnetization scheme yet).

**k-point sampling.** Standard Monkhorst–Pack k-point meshes (`meshk=[nx,ny,nz]`) are supported for crystalline/periodic calculations, alongside Γ-point-only molecular/cluster calculations.

**Capabilities.**
- Total-energy calculations for molecules, surfaces, and 3-D periodic crystalline systems within a periodic unit cell (no isolated/non-periodic corrections such as truncated-Coulomb methods are implemented, so care is needed for charged or strongly non-periodic systems)
- Band-structure calculations along high-symmetry k-paths (documented with a worked Si (fcc) example producing a band diagram)
- Both SCF and direct energy-minimization ground-state solvers

**Units.** Internally, PWDFT.jl uses Hartree atomic units throughout (energies in Hartree, lengths in bohr).

### 1.3 GPU Acceleration Status

PWDFT.jl as released and documented by its original authors is a **CPU-oriented, pure-Julia implementation**; the CPC (2020) paper and the GitHub README describe CPU benchmarking only (single-core and, implicitly, Julia's native threading/BLAS backends), with no native CUDA/GPU backend documented in the mainline package at the time of the paper. Its relevance to GPU-accelerated electronic structure work is chiefly:
- as a **reference/pedagogical implementation** used by later GPU-oriented and alternative-language plane-wave codes (e.g., the Python packages `eminus` and educational DFT++-style codes) for cross-validation of total energies;
- as an example of how Julia's array/GPU ecosystem (`CUDA.jl`, broadcasting, generic array types) *could* in principle be leveraged for GPU dispatch, in the same spirit as GPU work carried out in the related Julia DFT package DFTK.jl (a separate, independently developed plane-wave DFT package also written in Julia, which has had explicit GSoC-funded GPU-porting efforts). PWDFT.jl and DFTK.jl are distinct projects; they should not be conflated.

Users interested in GPU-accelerated plane-wave DFT specifically through PWDFT.jl should check the current state of the repository directly, as GPU support is an area that could evolve independently of the original 2020 publication.

### 1.4 Validation and Benchmarking

The CPC (2020) paper benchmarks PWDFT.jl (Julia 1.3.1) for both **accuracy** and **CPU performance** against three established plane-wave/pseudopotential codes chosen specifically because they can use the same GTH pseudopotentials and plane-wave basis:
- **ABINIT** (v8.6.1)
- **Quantum ESPRESSO / PWscf** (from QE v6.4)
- **GPAW** (v19.8.1)

Total energies from PWDFT.jl are validated against ABINIT as the primary reference, and the paper's stated conclusion is that PWDFT.jl reproduces ABINIT's total energies to good agreement, establishing it as a correctness-verified, if not performance-leading, implementation. Independent later work (e.g., the `eminus` Python package paper, 2024) has also cross-validated total energies against PWDFT.jl, JDFTx, and Quantum ESPRESSO for a common set of molecules/solids using SVWN5 XC, finding close agreement across all four independently implemented codes despite very different programming languages (Julia, C++, Fortran, Python) — a useful independent confirmation of PWDFT.jl's correctness.

### 1.5 Example Usage Pattern

```julia
# Build a structure
atoms = Atoms(xyz_file="CH4.xyz", LatVecs=gen_lattice_sc(16.0))

# Build the Hamiltonian
ecutwfc = 15.0   # Hartree
pspfiles = ["../pseudopotentials/pade_gth/C-q4.gth",
            "../pseudopotentials/pade_gth/H-q1.gth"]
Ham = Hamiltonian(atoms, pspfiles, ecutwfc)

# Solve
KS_solve_SCF!(Ham, betamix=0.2)
# or
KS_solve_Emin_PCG!(Ham)
```

### 1.6 Position in the Ecosystem

PWDFT.jl is one of several Julia-language entrants into electronic structure theory that appeared roughly contemporaneously (alongside Fermi.jl for quantum chemistry and, separately, DFTK.jl for plane-wave DFT), reflecting a broader circa-2019–2021 interest in evaluating Julia as a language for scientific HPC software — offering near-C/Fortran performance with much higher code transparency and easier extensibility than legacy Fortran plane-wave codes.

---

## 2. PWDFT (C++, NWChemEx / Bylaska)

### 2.1 Overview

This PWDFT is a **from-scratch, modern C++ rewrite** of the pseudopotential plane-wave (PSPW) and band-structure (BAND) capabilities historically found in **NWChem's NWPW module**, developed under the **NWChemEx** project umbrella (the next-generation successor to NWChem) with a specific mandate to target GPU-accelerated exascale-class supercomputers.

- **Repository:** `ebylaska/PWDFT` on GitHub, with project site at `https://ebylaska.github.io/PWDFT/`
- **Primary author/maintainer:** Eric J. Bylaska (Pacific Northwest National Laboratory, PNNL), with the NWChemEx development team
- **Tagline:** "Autonomous Plane-Wave DFT Engine for Exascale Materials Discovery"
- **Language:** C++ (built via CMake), interoperating with MPI, OpenMP, and vendor GPU toolchains
- **License:** GPL-3.0
- **Association:** Explicitly tied to the NWChemEx electronic-structure software ecosystem

### 2.2 Design and Implementation

The code is organized primarily under an `Nwpw/` source tree (mirroring NWChem's historical NWPW module naming) and is built with standard out-of-source CMake workflows:

```bash
cmake -S Nwpw/ -B build
cd build
make
```

It can also be built as a **shared library** (`libpwdft.so`/`.dylib`) exposing callable entry points for embedding in other applications (e.g., `pwdft::cpsd`, `pwdft::cpmd`, `pwdft::pspw_minimizer`, `pwdft::pspw_geovib`), reflecting its intended role as a reusable computational engine within the broader NWChemEx software stack rather than only a standalone executable.

**Physical/algorithmic capabilities** (inherited conceptually from and extending NWChem's NWPW/PSPW/BAND modules):
- **PSPW**: pseudopotential, Γ-point plane-wave DFT for molecules, liquids, crystals, and surfaces
- **BAND**: pseudopotential plane-wave band-structure code for crystals/surfaces with small band gaps (semiconductors, metals), using full Brillouin-zone k-point sampling
- **PAW**: projector augmented-wave capability, historically integrated into the PSPW module
- Car–Parrinello and Born–Oppenheimer (extended-Lagrangian) ab initio molecular dynamics (AIMD), including NVE/NVT ensembles
- Steepest-descent and conjugate-gradient wavefunction/geometry minimization
- Norm-conserving pseudopotentials (Hamann, Troullier–Martins, Hartwigsen–Goedecker–Hutter/HGH), with optional semicore corrections, read from CPI/TETER-format files and a bundled pseudopotential library
- LDA and PBE96 exchange–correlation functionals (spin-restricted and unrestricted), self-interaction correction (SIC), perturbative OEP, Hartree–Fock, and hybrid (exact-exchange) functionals
- Support for input decks in NWChem's native `.nw` format (the code reads and processes NWChem-style input, e.g. `pwdft cco-cu_surf30.nw`)

### 2.3 GPU Acceleration

This is the PWDFT codebase most directly and explicitly built for **multi-vendor GPU acceleration**, targeting the three major current GPU ecosystems via CMake build flags:

| Target | CMake flag | Toolchain | Example DOE system |
|---|---|---|---|
| NVIDIA CUDA | `-DNWPW_CUDA=ON` | `nvcc`/CUDA, cuBLAS, cuFFT, cuSOLVER | NERSC Perlmutter, ALCF Polaris |
| AMD HIP/ROCm | `-DNWPW_HIP=ON` | `hipcc`, `amd-mixed` module, `GPU_TARGETS=gfx90a` | OLCF Frontier |
| Intel SYCL | `-DNWPW_SYCL=On` | `icpx`/`icx`/`ifx` | ALCF Aurora/Sunspot |

The repository provides working build-and-submission recipes for each of these DOE leadership-class machines, including MPI-rank-to-GPU affinity binding scripts (e.g., Polaris's reverse-order GPU assignment script, Frontier's `srun --gpus-per-node=8 --gpu-bind=closest`, and Aurora's tile-compact binding), reflecting a mature, production-oriented HPC deployment rather than a research prototype.

**Documented performance characteristics.** The repository's own timing tables (for a CCO/Cu-surface test case) show:
- On a CPU-only workstation ("WE45090"), CPU time drops from ~27.9 s (1 core) to ~9.5 s (8 cores)
- On the same workstation with GPU acceleration ("WE45090-GPU"), CPU time drops to ~4.2–4.8 s, i.e., roughly a further ~2× speed-up over the best 8-core CPU-only result, dominated by non-local pseudopotential and FFM/FMF (wavefunction–potential contraction) operations, while the 3-D FFT becomes comparatively more GPU-costly (attributed to being memory/data-movement bound rather than compute bound)
- On NERSC Perlmutter, GPU runs at 1–4 GPU ranks (~1.1–2.9 s) substantially outperform CPU-only Perlmutter runs even at 8–16 CPU ranks (~1.8–2.7 s), while CPU-only scaling continues to help up to 128 ranks (~1.1–1.3 s), illustrating the classic plane-wave-DFT tension between GPU throughput on dense linear algebra/FFT kernels and MPI-parallel-FFT communication bottlenecks at scale

The authors explicitly note, based on these timings, that **parallel FFTs are memory-bound rather than compute-bound**, motivating hybrid MPI+OpenMP FFT strategies and exploration of alternative FFT pipelining/kernel strategies (including references to custom "Stockham" FFT kernel experiments in the repository's `Miscellaneous/programfft` directory) as an ongoing area of GPU-era algorithmic work.

### 2.4 Relationship to NWChem

This PWDFT is best understood as the **spiritual and architectural successor** to NWChem's mature PSPW/BAND/PAW plane-wave modules (which have existed within NWChem for over two decades and are documented extensively in the NWChem user manual and wiki), rewritten in modern C++ specifically to be GPU-portable across CUDA/HIP/SYCL and to serve as a standalone/embeddable library for the NWChemEx next-generation software ecosystem, rather than being locked into NWChem's legacy Fortran/Global-Arrays runtime.

---

## 3. PWDFT (C/C++, LBNL / DGDFT Group) and PWDFT-SW

### 3.1 Overview

A third, independently developed "PWDFT" is a C/C++ plane-wave DFT sub-module developed within the **DGDFT** (Discontinuous Galerkin DFT) software suite, associated with the group of Lin Lin, Chao Yang, and collaborators (originally Lawrence Berkeley National Laboratory), and more recently extended for the domestic Chinese "new Sunway" exascale supercomputer by Wei Hu, Jinlong Yang, and collaborators at USTC.

- **Nature of problem (per its CPC catalog entry):** electronic-structure calculations based on Kohn–Sham DFT, formulated as a constrained energy-minimization / nonlinear eigenvalue problem
- **Solution method:** SCF iterations and direct constrained-minimization algorithms with a variety of eigensolvers
- **Distribution:** CPC Library (Computer Physics Communications Program Library); developer repository historically hosted at `bitbucket.org/berkeleylab/scales`
- **Parallelization:** MPI/OpenMP multi-core, with a multi-level parallel strategy enabling calculations on systems with tens of thousands of atoms

### 3.2 PWDFT-SW: Sunway Extension

**PWDFT-SW** (Jiang, Cao, Chen, Qin, Hu, An, and Yang, 2024) is a substantially refactored parallel implementation of this PWDFT line specifically targeting the **new Sunway supercomputer**, whose processes are constrained to only 16 GB of memory each — a hard limit that makes conventional O(N²)-memory, O(N³)-compute plane-wave DFT workflows impractical at the many-thousand-atom scale without significant algorithmic and systems-level redesign.

Key contributions of PWDFT-SW:
- Extensive refactoring and recalibration of core plane-wave DFT kernels (including pseudopotential evaluation) to match the Sunway architecture's heterogeneous management-processing-element (MPE) / compute-processing-element (CPE) design
- Substantially reduced computational cost and memory footprint relative to the baseline implementation
- A demonstrated **64.8× speed-up** for a 4,096-silicon-atom system
- Scaling of plane-wave DFT calculations up to systems of **16,384 carbon atoms** — presented as extending the practical size limit of plane-wave DFT well beyond typical few-thousand-atom ceilings

This PWDFT line (and DGDFT more broadly) has also been run at extreme scale on the Sunway TaihuLight system (predecessor to the "new Sunway"), where DGDFT — of which PWDFT is a constituent sub-module — was scaled to over 8.5 million processor cores for systems with tens of thousands of atoms, a result that was highlighted as a case study in extreme-scale DFT HPC.

### 3.3 Distinguishing PWDFT (LBNL/DGDFT) from the Other Two

This PWDFT should not be confused with:
- **PWDFT.jl**, which is a small-to-medium-scale, single/multi-core Julia package aimed at pedagogical clarity and moderate-size systems, with no documented 16k-atom or Sunway-scale capability;
- **PWDFT (C++/NWChemEx, Bylaska)**, which targets DOE leadership-class GPU machines (CUDA/HIP/SYCL) rather than the Sunway (SW26010-class) many-core heterogeneous architecture, and is tied to the NWChem/NWChemEx software lineage rather than DGDFT/SCALES.

---

## 4. Comparative Summary

| Aspect | PWDFT.jl | PWDFT (C++, NWChemEx) | PWDFT (C/C++, DGDFT/LBNL) + PWDFT-SW |
|---|---|---|---|
| Language | Julia | C++ | C/C++ |
| Lead author(s) | F. Fathurrahman et al. (ITB) | E. J. Bylaska (PNNL) | Lin Lin / C. Yang lineage; PWDFT-SW: Q. Jiang, W. Hu, J. Yang et al. (USTC) |
| Lineage | New, independent implementation | Successor to NWChem NWPW (PSPW/BAND/PAW) | Sub-module of DGDFT/SCALES |
| Primary goal | Transparent, pedagogical, correctness-validated plane-wave DFT in a high-level language | GPU-portable, exascale-ready plane-wave engine for NWChemEx | Extreme-scale (10⁴+ atom) plane-wave DFT on leadership supercomputers |
| Pseudopotentials | GTH (norm-conserving, analytic) | Hamann, Troullier–Martins, HGH, PAW | Norm-conserving (implementation-specific) |
| XC functionals | LDA-VWN, GGA-PBE (via Libxc.jl) | LDA, PBE96, SIC, OEP, Hartree–Fock/hybrid | LDA/GGA class functionals |
| k-points | Γ-point and Monkhorst–Pack meshes | Γ-point (PSPW) and full BZ (BAND) | Primarily large-Γ-point/large-supercell regime |
| AIMD | Not a stated focus | Car–Parrinello, Born–Oppenheimer AIMD | Not the primary emphasis (ground-state HPC scaling) |
| GPU support | Not natively documented in mainline (CPU-focused) | Yes — CUDA, HIP/ROCm, SYCL (multi-vendor) | Not GPU-specific; instead targets Sunway SW26010-class many-core CPUs |
| Demonstrated scale | Molecules, surfaces, small-to-moderate crystals (validation-scale) | Surfaces/interfaces/nanomaterials on DOE GPU machines (Perlmutter, Frontier, Aurora, Polaris) | Up to 16,384 atoms (PWDFT-SW); DGDFT to 8.5M cores / tens of thousands of atoms |
| Reference publication | *Comput. Phys. Commun.* **256**, 107372 (2020) | No single unifying CPC-style paper identified; documented via NWChem/NWChemEx literature and repository | CPC Program Library entry ("PWDFT") + PWDFT-SW (2024, arXiv/IEEE) |
| License | Open on GitHub (see repo) | GPL-3.0 | CPC Program Library terms |

---

## 5. Publications Related to the Underlying Theory and Software

### 5.1 PWDFT.jl — primary and closely related references

1. **F. Fathurrahman, M. K. Agusta, A. G. Saputro, H. K. Dipojono**, "PWDFT.jl: A Julia package for electronic structure calculation using density functional theory and plane wave basis," *Computer Physics Communications* **256**, 107372 (2020). https://doi.org/10.1016/j.cpc.2020.107372 — *the primary citation for PWDFT.jl.*
2. **M. Bockstedte, A. Kley, J. Neugebauer, M. Scheffler**, "Density-functional theory calculations for polyatomic systems: Electronic structure, static and elastic properties and ab initio molecular dynamics," *Comput. Phys. Commun.* **107**, 187 (1997). — cited by PWDFT.jl's authors as methodological background for plane-wave/pseudopotential DFT implementation.
3. **S. Ismail-Beigi, T. A. Arias**, "New algebraic formulation of density functional calculation," *Comput. Phys. Commun.* **128**, 1–45 (2000). — theoretical basis for the direct energy-minimization (preconditioned conjugate-gradient) approach used in `KS_solve_Emin_PCG!`.
4. **C. Yang, J. C. Meza, B. Lee, L.-W. Wang**, "KSSOLV — a MATLAB toolbox for solving the Kohn–Sham equations," *ACM Trans. Math. Softw.* **36**, 1–35 (2009). — related solver-design reference.
5. **R. M. Martin**, *Electronic Structure: Basic Theory and Practical Methods*, Cambridge University Press (2004). — general plane-wave DFT theory textbook cited as background.
6. **J. Kohanoff**, *Electronic Structure Calculations for Solids and Molecules: Theory and Computational Methods*, Cambridge University Press (2006). — general theory reference cited by the project.
7. **D. Marx, J. Hutter**, *Ab Initio Molecular Dynamics: Basic Theory and Advanced Methods*, Cambridge University Press (2009). — cited as background for AIMD-adjacent theory (Car–Parrinello/Born–Oppenheimer MD), though AIMD is not itself a core PWDFT.jl feature.
8. **S. Goedecker, M. Teter, J. Hutter**, "Separable dual-space Gaussian pseudopotentials," *Phys. Rev. B* **54**, 1703 (1996) — the original GTH pseudopotential formalism used throughout PWDFT.jl (not explicitly listed in the repository's reference list surfaced here, but foundational to its pseudopotential implementation; recommended for completeness).
9. **J. P. Perdew, K. Burke, M. Ernzerhof**, "Generalized Gradient Approximation Made Simple," *Phys. Rev. Lett.* **77**, 3865 (1996) — foundational PBE GGA functional reference (accessed via Libxc.jl).
10. **S. H. Vosko, L. Wilk, M. Nusair**, "Accurate spin-dependent electron liquid correlation energies for local spin density calculations: a critical analysis," *Can. J. Phys.* **58**, 1200 (1980) — foundational VWN LDA correlation functional reference (accessed via Libxc.jl).

### 5.2 PWDFT (C++, NWChemEx / Bylaska) and its NWChem/PSPW theoretical lineage

11. **E. J. Bylaska**, "Plane-Wave DFT Methods for Chemistry," in *Annual Reports in Computational Chemistry*, Vol. 13, Elsevier (2017), Chapter on PSPW methodology — a detailed methodological description of the pseudopotential plane-wave (PSPW) method as implemented in NWChem/PWDFT, covering geometry optimization, AIMD, exact exchange, and AIMD/MM.
12. **E. J. Bylaska, K. Tsemekhman, S. B. Baden, J. H. Weare, H. Jónsson**, "Parallel implementation of γ-point pseudopotential plane-wave DFT with exact exchange," *J. Comput. Chem.* **32**, 54–69 (2011). https://doi.org/10.1002/jcc.21598 — key methods paper for exact-exchange/hybrid-DFT capability in the PSPW/PWDFT line.
13. **E. Aprà, E. J. Bylaska, W. A. de Jong, N. Govind, K. Kowalski, T. P. Straatsma, M. Valiev, H. J. J. van Dam, Y. Alexeev, J. Anchell, et al.**, "NWChem: Past, present, and future," *J. Chem. Phys.* **152**, 184102 (2020). https://arxiv.org/abs/2004.12023 — comprehensive overview of NWChem including the NWPW (PSPW/BAND) plane-wave module architecture that PWDFT (C++) modernizes.
14. **E. J. Bylaska, E. Aprà, K. Kowalski, M. Jacquelin, W. A. de Jong, A. Vishnu, B. Palmer, J. Daily, T. P. Straatsma, J. R. Hammond, M. Klemm**, "Transitioning NWChem to the Next Generation of Manycore Machines," report/chapter, ORNL/PNNL/OSTI (2017). https://doi.org/10.2172/1422349 — describes early many-core/thread-level parallelization of the NWPW/AIMD methods, directly relevant to the motivation for the GPU-native C++ rewrite.
15. **E. J. Bylaska, K. M. Rosso**, background reference (2018), cited in periodic one/two-electron integral formulations for the pseudopotential plane-wave method (see item 17 below) — Hamiltonian construction methodology underlying PWDFT (C++).
16. **R. Car, M. Parrinello**, "Unified Approach for Molecular Dynamics and Density-Functional Theory," *Phys. Rev. Lett.* **55**, 2471 (1985). — the foundational Car–Parrinello AIMD method implemented in the PSPW/PWDFT Car–Parrinello (`cpmd`) driver.
17. **W. Song, E. J. Bylaska, et al.**, "Periodic One-Electron and Two-Electron Integrals using the Pseudopotential Plane-Wave Method," *Materials Theory* **7**, 2 (2023). https://doi.org/10.1186/s41313-022-00049-5 — recent methodological paper detailing integral evaluation (Filon integration for exact-exchange Brillouin-zone folding, Bessel-transform techniques) directly underlying PWDFT (C++)'s Hamiltonian and exact-exchange machinery.
18. NWChem/NWPW user documentation: "Plane-Wave Density Functional Theory (NWPW)," NWChem User Manual / Wiki (`nwchemgit.github.io`) — the operational/theoretical reference manual describing PSPW, BAND, and PAW task input, algorithms, and capabilities that PWDFT (C++) reimplements and extends for GPU platforms.

### 5.3 PWDFT (C/C++, DGDFT/LBNL) and PWDFT-SW

19. CPC Program Library catalog entry, "Plane wave density functional theory (PWDFT)," program summary describing the C/C++ PWDFT module (SCF and constrained-minimization solvers for the Kohn–Sham nonlinear eigenvalue problem), associated with the `berkeleylab/scales` repository.
20. **W. Hu, X. Qin, Q. Jiang, J. Chen, H. An, W. Jia, F. Li, X. Liu, D. Chen, J. Yang**, "Extreme-Scale Density Functional Theory High Performance Computing of DGDFT for Tens of Thousands of Atoms using Millions of Cores on Sunway TaihuLight," (2020). https://arxiv.org/abs/2003.00407 — describes the DGDFT method (of which PWDFT is a sub-module) and its extreme-scale HPC deployment.
21. **Q. Jiang, Z. Cao, J. Chen, X. Qin, W. Hu, H. An, J. Yang**, "PWDFT-SW: Extending the Limit of Plane-Wave DFT Calculations to 16K Atoms on the New Sunway Supercomputer," *IEEE Transactions on Parallel and Distributed Systems* (2024/2025); preprint: arXiv:2406.10765. https://arxiv.org/abs/2406.10765 — the primary PWDFT-SW methods and performance paper.
22. **Q. Jiang, Z. Cao, X. Cui, L. Wan, X. Qin, H. Cao, H. An, J. Chen, J. Liu, W. Hu, et al.**, "Extending the limit of LR-TDDFT on two different approaches: Numerical algorithms and new Sunway heterogeneous supercomputer," *Parallel Computing* **120**, 103085 (2024) — related work by the same group extending time-dependent DFT capability on the same Sunway hardware/software stack.

### 5.4 Broader context and comparative literature (cross-cited alongside "PWDFT")

23. **M. F. Herbst, A. Levitt, E. Cancès**, "DFTK: A Julian approach for simulating electrons in solids," *Proc. JuliaCon Conf.* **3**, 69 (2021). https://doi.org/10.21105/jcon.00069 — a separate Julia plane-wave DFT package (DFTK.jl) frequently discussed alongside PWDFT.jl in the Julia scientific-computing ecosystem; not the same project, but relevant comparative/contextual reading, including its own independent GPU-porting efforts (GSoC 2022).
24. **X. Andrade et al.**, INQ framework — a GPU-native, from-scratch real-time/ground-state plane-wave DFT code, cited comparatively in GPU-DFT literature alongside PWDFT/PWDFT-SW.
25. **J. Gao et al.**, "PyPWDFT: A Lightweight Python Software for Single-Node 10K Atom Plane-Wave Density Functional Theory Calculations," *J. Chem. Theory Comput.* **21**, 2353 (2025) — a related, differently authored Python plane-wave DFT code explicitly benchmarked against/discussed alongside the "PWDFT" family for large-atom-count single-node calculations.
26. **T. Schmidt et al.** (eminus development team), "eminus – Pythonic electronic structure theory," (2024), https://arxiv.org/abs/2410.19438 — independent Python plane-wave DFT package that cross-validates total energies against PWDFT.jl, JDFTx, and Quantum ESPRESSO.

---

## 6. Practical Notes for Users

- **If your goal is to learn plane-wave DFT internals or prototype algorithms in a high-level, hackable language:** PWDFT.jl is the natural choice — its explicit `Atoms`/`Hamiltonian` data types and pure-Julia SCF/CG solvers are unusually transparent compared to legacy Fortran plane-wave codes, at the cost of GPU support not being a documented, first-class feature.
- **If your goal is production GPU-accelerated plane-wave DFT on DOE leadership machines (Perlmutter, Frontier, Aurora, Polaris) with AIMD/Car–Parrinello capability and a PAW/exact-exchange feature set inherited from NWChem:** the C++ PWDFT under NWChemEx (`ebylaska/PWDFT`) is the relevant code, with explicit CUDA/HIP/SYCL build paths.
- **If your goal is pushing plane-wave DFT to many-thousand-atom system sizes on extreme-scale (especially Sunway-class) hardware:** the DGDFT-associated C/C++ PWDFT and its PWDFT-SW extension are the relevant references, though this line is less GPU-oriented and more oriented toward Sunway's heterogeneous many-core-CPU architecture.
- **Always verify current repository state before citing feature claims**, especially for GPU support in PWDFT.jl and for version-specific pseudopotential/XC-functional coverage in all three codes, as all remain under active development and this review reflects information gathered as of August 2026.

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the PWDFT / PWDFT.jl 	Plane-wave DFT codes (C++/Julia implementations) developed for large-scale and GPU-accelerated electronic structure calculations. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
