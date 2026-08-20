# ABACUS — An Exhaustive Review

**Atomic-orbital Based Ab-initio Computation at USTC**

---

## 1. Overview

ABACUS is an open-source, first-principles electronic structure and molecular dynamics package built around Kohn-Sham Density Functional Theory (DFT). Its defining architectural feature is **dual basis-set support**: calculations can be performed using either a **plane-wave (PW) basis** or a **numerical atomic orbital (NAO/LCAO) basis**, with a third hybrid mode (LCAO-in-PW) also available. This lets users trade off between the systematic completeness of plane waves and the computational efficiency/locality of atomic orbitals within a single, consistent codebase — allowing convenient cross-validation of results between the two representations.

The project was initiated in 2007 by **Prof. Lixin He's group** at the Key Laboratory of Quantum Information, Chinese Academy of Sciences (CAS), University of Science and Technology of China (USTC). Development has since expanded into a multi-institutional collaboration involving the **Institute of Physics (CAS)**, **Peking University**, **AI for Science Institute (AISI) Beijing**, and the **Hefei Institute of Artificial Intelligence**, and the code is now maintained as part of the **DeepModeling** open-source ecosystem.

- **License:** GNU GPLv3
- **Primary repository:** `github.com/deepmodeling/abacus-develop` (mirrored at `github.com/abacusmodeling/abacus-develop`)
- **Official site:** abacus.ustc.edu.cn (documentation portal: abacus.deepmodeling.com)
- **Language:** C++ (C++11+), MPI-parallelized, with GPU backends
- **Distribution:** Source build (CMake or Makefile), conda-forge package, Docker containers, and precompiled releases

---

## 2. Design Philosophy and Basis Sets

### 2.1 Plane-Wave (PW) Basis
The PW mode expands the Kohn-Sham orbitals in a plane-wave basis truncated by an energy cutoff (`ecutwfc`), combined with pseudopotentials to treat core electrons. This is the standard approach used by codes such as Quantum ESPRESSO and VASP, offering systematic, easily-converged accuracy and straightforward force/stress evaluation.

### 2.2 Numerical Atomic Orbitals (NAO/LCAO)
ABACUS's namesake strength lies in its **numerically tabulated, atom-centered basis functions**, generated via an optimization scheme originally developed by **Chen, Guo, and He (CGH)**. These orbitals are organized into hierarchical, systematically improvable basis sets (e.g., single-ζ, double-ζ, double-ζ + polarization "DZP", up to more complete sets), allowing users to dial in the accuracy/cost trade-off. A **second-generation NAO construction strategy** (Lin, Ren, He, 2021) improved basis compactness by incorporating gradients of reference wavefunctions, and more recent work has explored NAOs built from **contracted truncated spherical waves** for further systematic improvability.

Because LCAO Hamiltonian and overlap matrices are sparse (orbitals are localized), ABACUS can exploit this sparsity for large-scale and linear/sub-cubic-scaling algorithms, in contrast to the dense-matrix operations inherent to PW methods — this is central to ABACUS's positioning as a **large-scale / high-throughput** code.

### 2.3 Pseudopotentials
ABACUS uses **norm-conserving pseudopotentials** to describe nuclear ion–valence electron interactions, compatible with standard pseudopotential libraries (e.g., SG15, PD03, and other UPF-format sets), alongside a curated library of matching numerical orbitals.

---

## 3. Core Functionality

### 3.1 Electronic Structure Methods
- Kohn-Sham DFT (ground-state SCF total energy, band structure, DOS/PDOS)
- **Stochastic DFT** (PW only) — enables finite-temperature/warm-dense-matter calculations at reduced cost by stochastic trace estimation
- **Orbital-free DFT (OF-DFT)**, including machine-learning-based non-local kinetic energy density functionals
- **Real-time time-dependent DFT (rt-TDDFT)**, including a hybrid-gauge approach for improved accuracy with NAOs
- DFT + Dynamical Mean-Field Theory (DFT+DMFT) within the LCAO framework (via external interfacing)

### 3.2 Exchange-Correlation Functionals
- LDA
- GGA (PBE and related)
- meta-GGA (via LibXC; includes SCAN, rSCAN, r²SCAN implementations in NAO basis)
- Hybrid functionals: **PBE0** and **HSE06** (LCAO basis, using localized resolution-of-identity techniques for efficiency; PW support is also present/under continued development)

### 3.3 Structural and Dynamical Simulation
- Geometry optimization: Conjugate Gradient (CG), BFGS, FIRE
- Cell relaxation and stress-tensor calculations (including accurate NAO-basis stress implementations)
- Ab initio molecular dynamics (AIMD) with NVE/NVT ensembles
- Interfaces to machine-learned interatomic potentials: **DeePMD**, Lennard-Jones potentials
- Interface to **Phonopy** for phonon calculations
- Interface to **DP-GEN** for automated training-data generation workflows

### 3.4 Response Properties and Advanced Analysis
- Electric polarization via the **Berry phase** method
- Berry curvature calculation using non-orthogonal atomic orbitals
- Interface to **Wannier90** for maximally localized Wannier functions
- Mulliken population/charge analysis (LCAO)
- Projected density of states (PDOS) and band unfolding (for disordered/alloyed/doped systems)
- Electrostatic potential output
- Electric field and dipole correction (for slab/surface calculations)
- DFT+U (LCAO)
- Semi-empirical van der Waals corrections (Grimme DFT-D2/D3)
- Implicit solvation model
- Constant-potential / grand-canonical electrochemistry simulation (added in v3.9)

### 3.5 Beyond-DFT / Correlated Methods
- **Random Phase Approximation (RPA)** correlation energy via the companion **LibRPA** package — subquadratic- and low-scaling real-space RPA implementations using NAOs, applicable to periodic systems
- G0W0 calculations (via interfacing, for benchmark/cluster studies)

### 3.6 Machine-Learning / AI Integration
ABACUS is explicitly positioned (per its own 2025 flagship publication) as an **"Electronic Structure Analysis Package for the AI Era."** It functions as a data-generation and integration platform for:
- **DeePKS** — a neural-network correction scheme bridging DFT and machine-learned potentials, bringing near-hybrid-functional accuracy at LCAO cost
- **DeePMD-kit / DP-GEN** — machine-learned interatomic potential training and active-learning workflows, including generation of **DPA (Deep Potential with Attention)** large atomic models
- **DeepH** — deep-learning prediction of DFT Hamiltonians for large-scale electronic-structure inference
- **DeePTB** — deep-learning tight-binding models
- **HamGNN** — graph neural network Hamiltonian prediction

---

## 4. High-Performance Computing and Scalability

ABACUS is architected for **large-scale and high-throughput** first-principles calculations:

- **MPI parallelization** scaling to O(10³) CPU cores across both PW and LCAO code paths
- **GPU acceleration**: CUDA support for PW calculations, and more recently dedicated **GPU acceleration of NAO/LCAO-based algorithms** (2024 work demonstrates substantial speedups for LCAO Hamiltonian construction, density-matrix operations, and force/stress evaluation on GPUs, addressing a gap relative to other LCAO codes such as SIESTA and FHI-aims)
- Support for domestic Chinese accelerator hardware (e.g., DCU) alongside standard CUDA GPUs
- **ABACUS 3.9** (Dec 2024) delivered a major optimization of **Davidson diagonalization** for PW calculations, more than halving CPU-side PW computation time and improving DCU efficiency by ~30%, alongside new constant-potential simulation capability
- **ABACUS LTS (Long-Term Support) v3.10** (April 2025) — a stability-focused release aimed at supporting a broader research ecosystem
- Pole-expansion-and-selected-inversion (PEXSI) techniques have been explored to accelerate atomic-orbital-based electronic structure calculations
- High-throughput workflow tooling: the **APNS (ABACUS-Pseudopotential-Numerical atomic orbital-Square)** framework supports systematic, automated benchmarking across pseudopotential/orbital combinations for high-throughput screening applications

---

## 5. Ecosystem and Interfaces

| Component | Role |
|---|---|
| **LibRPA** | Companion package for low-scaling RPA correlation-energy calculations using NAOs |
| **PYATB** | Python package for ab initio tight-binding-model-based electronic structure analysis, built on ABACUS outputs |
| **DeePKS-kit** | Neural-network-based DFT accuracy correction framework |
| **DeePMD-kit / DP-GEN** | Machine-learned potential training and active learning |
| **DeepH / DeePTB / HamGNN** | Deep-learning Hamiltonian/tight-binding prediction |
| **Wannier90** | Maximally localized Wannier function construction |
| **Phonopy** | Phonon property calculation |
| **ShengBTE** | Boltzmann transport equation solver for thermal transport (compatible via ABACUS force constants) |

ABACUS is distributed and versioned as part of the broader **DeepModeling** open-source community, which also hosts DeePMD-kit and related AI-for-science tooling — reflecting the project's strategic emphasis on bridging first-principles accuracy with machine-learning-accelerated simulation at scale.

---

## 6. Installation and Access

- **Build systems:** CMake (≥3.16, recommended) or traditional Makefiles
- **Requirements:** C++11-compatible compiler (GCC ≥5 or Intel C++), MPI (Intel MPI / MPICH / Open MPI), BLAS/LAPACK (OpenBLAS or Intel MKL), FFTW3; optional ScaLAPACK, ELPA, LibXC, LibRI, Wannier90, LibTorch (for ML interfaces), and CUDA/ROCm for GPU builds
- **Distribution channels:** GitHub source releases, conda-forge (`conda install -c conda-forge abacus`), Docker images, and precompiled binaries from the official site
- **Documentation:** Full manual, input-parameter reference, and tutorials hosted at abacus.deepmodeling.com, with example inputs in the repository

---

## 7. Positioning Relative to Other DFT Codes

- Versus **pure PW codes** (Quantum ESPRESSO, VASP, ABINIT): ABACUS offers comparable PW functionality plus a native, systematically improvable NAO option in the same codebase, avoiding the need to switch software when moving between basis-set philosophies.
- Versus **other NAO/LCAO codes** (SIESTA, OpenMX, FHI-aims, CP2K): ABACUS distinguishes itself through its deep native integration with machine-learning workflows (DeePKS, DeePMD, DP-GEN) and its dual-basis architecture, and — per its developers — is explicitly engineered toward large-scale/high-throughput data generation for training general-purpose ML interatomic potentials, an emphasis less central to some legacy LCAO codes.
- Active GPU-acceleration efforts for LCAO specifically address a historically underdeveloped area relative to GPU-accelerated PW codes.

---

## 8. Applications Reported in the Literature

Published applications spanning ABACUS's user base illustrate its breadth: 2D materials and moiré systems (twisted bilayers, kagome monolayers), topological and magnetic materials (MnBi2Te4, CrTe2 films, bilayer nickelate spin-density waves), battery/electrochemistry research (Li-metal dendrite formation, sodium-ion battery electrodes, solid electrolytes), warm dense matter and shock physics (via stochastic DFT), semiconductor and perovskite solar-cell materials, catalysis and surface chemistry, nuclear materials (deuterium/tritium retention in liquid metals), and machine-learned potential training pipelines for metals, alloys, and oxides.

---

## 9. Summary Assessment

**Strengths**
- Genuine dual-basis (PW + NAO) architecture within one consistent code, enabling internal cross-validation
- Systematically improvable, well-documented NAO basis-set hierarchy
- Strong, first-class integration with the machine-learning/AI-for-science ecosystem (DeePMD, DeePKS, DP-GEN, DeepH)
- Active large-scale HPC development, including growing GPU support for both PW and LCAO paths
- Fully open source (GPLv3) with an active multi-institutional developer base and frequent releases

**Considerations**
- Some advanced features (e.g., hybrid functionals in PW mode) are noted as being under continued development/testing
- Documentation and orbital/pseudopotential library coverage, while extensive, is less exhaustive than some longer-established codes (e.g., VASP's ecosystem)
- As a rapidly-evolving research code with frequent point releases, users doing production/high-throughput work may need to pin specific versions for reproducibility

---

## 10. Key Citation (General Reference)

> Z.-K. Cao et al. (ABACUS Developer Team), *ABACUS: An Electronic Structure Analysis Package for the AI Era*, arXiv:2501.08697 (2025); published in *J. Chem. Phys.* (2025).

---

# Publications Related to ABACUS Theory and Methods

The following list covers the package's **methods and theoretical development** literature — the algorithms, basis-set constructions, and formalisms underlying ABACUS — as compiled from the official ABACUS publications page and cross-referenced against primary sources. (A separate, much larger "Applications" literature exists — illustrative examples are cited in Section 8 above but are not exhaustively reproduced here, as the request specifies package-theory publications.)

## Foundational / General-Purpose Papers

1. Li, P.; Liu, X.; Chen, M.; Lin, P.; Ren, X.; Lin, L.; Yang, C.; He, L. *Large-scale ab initio simulations based on systematically improvable atomic basis.* **Comput. Mater. Sci.** 112, 503–517 (2016).
2. Lin, P.; Ren, X.; Liu, X.; He, L. *Ab initio electronic structure calculations based on numerical atomic orbitals: Basic formalisms and recent progresses.* **WIREs Comput. Mol. Sci.** 14, e1687 (2024).
3. Cao, Z.-K. et al. (ABACUS Developer Team). *ABACUS: An Electronic Structure Analysis Package for the AI Era.* **arXiv:2501.08697** (2025).
4. Liu, X.; Chen, M.; Li, P.; Shen, Y.; Ren, X.; Guo, G.-C.; He, L. *Introduce first-principles simulation package ABACUS based on systematically improvable atomic orbitals.* **Acta Phys. Sin.** 64, 187104 (2015).

## Numerical Atomic Orbital (NAO) Basis-Set Construction

5. Chen, M.; Fang, W.; Sun, G.-Z.; Guo, G.-C.; He, L. *Method to construct transferable minimal basis sets for ab initio calculations.* **Phys. Rev. B** 80, 165121 (2009).
6. Chen, M.; Guo, G.-C.; He, L. *Systematically improvable optimized atomic basis sets for ab initio calculations.* **J. Phys.: Condens. Matter** 22, 445501 (2010). *(First-generation NAO basis)*
7. Chen, M.; Guo, G.-C.; He, L. *Electronic structure interpolation via atomic orbitals.* **J. Phys.: Condens. Matter** 23, 325501 (2011).
8. Lin, P.; Ren, X.; He, L. *Strategy for constructing compact numerical atomic orbital basis sets by incorporating the gradients of reference wavefunctions.* **Phys. Rev. B** 103, 235131 (2021). *(Second-generation NAO basis)*
9. Huang, Y. et al. *Systematically Improvable Numerical Atomic Orbital Basis Using Contracted Truncated Spherical Waves.* **arXiv:2603.13995** (2026 preprint).

## Hybrid Functionals and Exchange-Correlation Methods

10. Lin, P.; Ren, X.; He, L. *Accuracy of Localized Resolution of the Identity in Periodic Hybrid Functional Calculations with Numerical Atomic Orbitals.* **J. Phys. Chem. Lett.** 11, 3082–3088 (2020).
11. Lin, P.; Ren, X.; He, L. *Efficient Hybrid Density Functional Calculations for Large Periodic Systems Using Numerical Atomic Orbitals.* **J. Chem. Theory Comput.** 17(1), 222–239 (2021).
12. Lin, P.; Ji, Y.; He, L.; Ren, X. *Efficient Hybrid-Functional-Based Force and Stress Calculations for Periodic Systems with Thousands of Atoms.* **J. Chem. Theory Comput.** 21, 3394 (2025).
13. Liu, R.; Zheng, D.; Liang, X.; Ren, X.; Chen, M.; Li, W. *Implementation of the meta-GGA exchange-correlation functional in numerical atomic orbital basis: With systematic testing on SCAN, rSCAN, and r²SCAN functionals.* **J. Chem. Phys.** 159(7), 074109 (2023).

## Forces, Stress, and Structural Optimization

14. Zheng, D.; Ren, X.; He, L. *Accurate stress calculations based on numerical atomic orbital bases: Implementation and benchmarks.* **Comput. Phys. Commun.** 267, 108043 (2021).

## Algorithmic / Numerical Methods

15. Lin, L.; Chen, M.; Yang, C.; He, L. *Accelerating atomic orbital-based electronic structure calculation via pole expansion and selected inversion.* **J. Phys.: Condens. Matter** 25, 295501 (2013).
16. Luo, K.; Wang, T.; Ren, X. *Direct minimization on the complex Stiefel manifold in Kohn-Sham density functional theory for finite and extended systems.* **Comput. Phys. Commun.** 312, 109596 (2025).
17. Liu, Q.; Chen, M. *Plane-wave-based stochastic-deterministic density functional theory for extended systems.* **Phys. Rev. B** 106, 125132 (2022).

## GPU Acceleration and HPC

18. (2024 preprint/study) *GPU Acceleration of Numerical Atomic Orbitals-Based Density Functional Theory Algorithms within the ABACUS package.* **arXiv:2409.09399** (2024).

## Response Properties, Berry Phase, and Band-Structure Methods

19. Jin, G.; Zheng, D.; He, L. *Calculation of Berry curvature using non-orthogonal atomic orbitals.* **J. Phys.: Condens. Matter** 33, 325503 (2021).
20. Dai, Z.; Jin, G.; He, L. *First-principles calculations of the surface states of doped and alloyed topological materials via band unfolding method.* **Comput. Mater. Sci.** 213, 111656 (2022).
21. Jin, G.; Pang, H.; Ji, Y.; Dai, Z.; He, L. *PYATB: an efficient Python package for electronic structure calculations using ab initio tight-binding model.* **Comput. Phys. Commun.** 291, 108844 (2023).

## DFT+U and Correlated-Electron Methods

22. Qu, X.; Jiang, H.; He, L.; Ren, X. *DFT+U within the framework of linear combination of numerical atomic orbitals.* **J. Chem. Phys.** 156, 234104 (2022).
23. Qu, X.; Xu, P.; Li, R.; Li, G.; He, L.; Ren, X. *Density Functional Theory Plus Dynamical Mean Field Theory within the Framework of Linear Combination of Numerical Atomic Orbitals Formulation and Benchmarks.* **J. Chem. Theory Comput.** 18(9), 5589–5606 (2022).

## Random Phase Approximation (RPA) / Beyond-DFT Correlation

24. Shi, R.; Lin, P.; Zhang, M.-Y.; He, L.; Ren, X. *Subquadratic-scaling real-space random phase approximation correlation energy calculations for periodic systems with numerical atomic orbitals.* **Phys. Rev. B** 109, 035103 (2024).
25. Shi, R.; Zhang, M.-Y.; Lin, P.; He, L.; Ren, X. *LibRPA: A software package for low-scaling first-principles calculations of random phase approximation electron correlation energy based on numerical atomic orbitals.* **Comput. Phys. Commun.** 309, 109496 (2025).

## Real-Time TDDFT

26. He, F.; Ren, X.; Jiang, J.; Zhang, G.; He, L. *Real-time, time-dependent density functional theory study on photoinduced isomerizations of azobenzene under a light field.* **J. Phys. Chem. Lett.** 13, 427 (2022).
27. Zhao, H.; He, L. *Hybrid Gauge Approach for Accurate Real-Time TDDFT Simulations with Numerical Atomic Orbitals.* **J. Chem. Theory Comput.** 21, 3335–3341 (2025).

## Machine-Learning / AI Integration Methods (DFT–ML Interfaces)

28. Li, W.; Ou, Q.; Chen, Y.; Cao, Y.; Liu, R.; Zhang, C.; Zheng, D.; Cai, C.; Wu, X.; Wang, H.; Chen, M.; Zhang, L. *DeePKS+ABACUS as a bridge between expensive quantum mechanical models and machine learning potentials.* **J. Phys. Chem. A** 126, 9154 (2022).
29. Ou, Q.; Tuo, P.; Li, W.; Wang, X.; Chen, Y.; Zhang, L. *DeePKS model for halide perovskites with the accuracy of a hybrid functional.* **J. Phys. Chem. C** 127, 18755 (2023).
30. Sun, L.; Chen, M. *Truncated non-local kinetic energy density functionals for simple metals and Si.* **Phys. Rev. B** 108, 075158 (2023).
31. Sun, L.; Chen, M. *Machine-learning-based non-local kinetic energy density functional for simple metals and alloys.* **Phys. Rev. B** 109, 115135 (2024).
32. Sun, L.; Chen, M. *Multi-channel machine learning based nonlocal kinetic energy density functional for semiconductors.* **Electronic Structure** 6, 045006 (2024).
33. Tang, Z.; Li, H.; Lin, P.; Gong, X.; Jin, G.; He, L.; Jiang, H.; Ren, X.; Duan, W.; Xu, Y. *A deep equivariant neural network approach for efficient hybrid density functional calculations.* **Nat. Commun.** 15, 8815 (2024).
34. Gu, Q.; Zhouyin, Z.; Pandey, S. K.; Zhang, P.; Zhang, L.; E, W. *Deep learning tight-binding approach for large-scale electronic simulations at finite temperatures with ab initio accuracy.* **Nat. Commun.** 15, 6772 (2024).
35. Zhang, D.; Liu, X.; Zhang, X. et al. *DPA-2: a large atomic model as a multitask learner.* **npj Comput. Mater.** 10, 293 (2024).
36. Liang, X.; Liu, R.; Chen, M. *A Deep Learning Framework for the Electronic Structure of Water: Toward a Universal Model.* **J. Chem. Theory Comput.** 21, 14 (2025).
37. Zeng, J.; Zhang, D.; Peng, A.; et al. *DeePMD-kit v3: A Multiple-Backend Framework for Machine Learning Potentials.* **J. Chem. Theory Comput.** 21, 4375 (2025).
38. Zeng, J.; Giese, T. J.; Zhang, D.; Wang, H.; York, D. M. *DeePMD-GNN: A DeePMD-kit Plugin for External Graph Neural Network Potentials.* **J. Chem. Inf. Model.** 65, 3154 (2025).

## Chinese-Language Methods Papers

39. 陈默涵. *密度泛函理论软件ABACUS进展及其与深度学习算法的融合及应用* (Progress of the DFT software ABACUS and its integration with deep learning algorithms). **金属学报 (Acta Metall. Sin.)**, Issue 10 (2024).
40. 沈瑜, 李会民, 刘晓辉. *第一性原理计算软件包ABACUS中格点积分的优化* (Optimization of grid integration in ABACUS). **科研信息化技术与应用** 6, 12 (2015).
41. 赵慰, 赵永华, 刘晓辉, 何力新. *第一性原理计算软件包在GPU集群上的加速* (Acceleration of the first-principles software package on GPU clusters). **计算机科学与探索** 8, 897 (2014).

---

### Sources
- Official ABACUS website — Publications page: https://abacus.ustc.edu.cn/publication/list.htm
- Official ABACUS website — Home/Citation page: https://abacus.ustc.edu.cn/main.htm
- Official ABACUS website — Features page: https://abacus.ustc.edu.cn/features/list.htm
- GitHub repository: https://github.com/deepmodeling/abacus-develop
- ABACUS flagship 2025 paper: arXiv:2501.08697 / PubMed: 41263655
- GPU-LCAO acceleration study: arXiv:2409.09399
- NAO contracted spherical wave basis study: arXiv:2603.13995


---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the Abacus 	Open-source DFT package supporting both plane-wave and numerical atomic orbital bases, developed for large-scale and high-throughput calculations. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
