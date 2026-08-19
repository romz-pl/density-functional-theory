# CONQUEST: A Review of the Linear-Scaling DFT Code for Million-Atom Simulations

## 1. Overview

CONQUEST is an open-source, massively parallel density functional theory (DFT) code designed from its inception to enable electronic-structure calculations on very large systems — from a few hundred atoms up to several million. It is developed jointly by University College London (UCL, London Centre for Nanotechnology), the National Institute for Materials Science (NIMS, Japan), and the Institut des Sciences Moléculaires (ISM, Université de Bordeaux). The project traces its origin to work by Mike Gillan, Emilio Hernández, and Chris Goringe in the mid-1990s, and has been continuously developed since, most recently under the stewardship of David Bowler, Tsuyoshi Miyazaki, Ayako Nakata, and Lionel Truflandier.

CONQUEST distinguishes itself from most other DFT packages by offering **two solvers within a single code**:

- An **exact diagonalisation** solver, giving conventional, cubic-scaling but fully accurate DFT, usable for systems up to roughly 1,000–10,000 atoms depending on basis choice.
- A **linear-scaling (O(N)) solver**, based on direct optimisation of the density matrix, which allows the computational cost (both time and memory) to grow linearly with the number of atoms — enabling simulations of **millions of atoms** on large HPC platforms.

This dual capability lets users move seamlessly along an accuracy/cost hierarchy — from empirical or ab initio tight-binding, through self-consistent tight-binding DFT, to full plane-wave-accuracy DFT — within one consistent framework.

CONQUEST is released as free, open-source software under the **MIT licence**, with source code hosted on GitHub (`OrderN/CONQUEST-release`), full documentation on ReadTheDocs, and integration with the Atomic Simulation Environment (ASE). As of late 2024 the code stood at version 1.4, with 1.5 documentation in preparation.

---

## 2. Theoretical Foundations

### 2.1 Support functions and the density matrix

Rather than solving directly for extended Kohn–Sham eigenstates, CONQUEST represents the Kohn–Sham orbitals $\Psi_\nu(\mathbf{r})$ as linear combinations of strictly localised **support functions** $\phi_{i\alpha}(\mathbf{r})$ centred on atoms:

$$\Psi_\nu(\mathbf{r}) = \sum_{i\alpha} c_{i\alpha}(\nu)\,\phi_{i\alpha}(\mathbf{r})$$

The one-particle density matrix is written in a separable form,

$$\rho(\mathbf{r},\mathbf{r}') = \sum_{i\alpha,j\beta} \phi_{i\alpha}(\mathbf{r})\, L_{i\alpha,j\beta}\, \phi_{j\beta}(\mathbf{r}')$$

where $L$ is the density matrix expressed in the support-function basis. Because $\phi_{i\alpha}$ is strictly zero outside a chosen cutoff radius around atom $i$, and because $L_{i\alpha,j\beta}$ is forced to zero once atoms $i$ and $j$ are separated by more than a chosen range, both the basis and the density matrix become **sparse**. Exploiting this sparsity in the storage and multiplication of all matrices is what allows the memory and CPU cost to grow linearly, rather than cubically, with system size $N$.

### 2.2 Three nested optimisation loops

The self-consistent ground state is found through three nested loops:

1. **Inner loop** — with support functions and charge density fixed, the density matrix (or Kohn–Sham eigenstates, if diagonalising) is optimised to minimise the energy.
2. **Middle loop** — self-consistency between the charge density and the effective (Hartree + exchange-correlation) potential is achieved by iteratively reducing the density residual.
3. **Outer loop** — the support functions themselves are variationally optimised (relevant for the blip-function basis, where the support functions are not fixed a priori).

### 2.3 Linear-scaling density matrix optimisation (O(N))

For the O(N) solver, CONQUEST uses the **density matrix minimisation (DMM)** method originally proposed by Li, Nunes and Vanderbilt, combined with an auxiliary density matrix and McWeeny purification to enforce idempotency and the correct electron count while avoiding diagonalisation entirely. This machinery was established in the code's founding papers by Hernández, Gillan and Goringe (1995–1997) and refined by Bowler and Gillan (1998–1999), who identified and solved several forms of numerical **ill-conditioning** that arise when minimising energy with respect to strictly localised functions (particularly for blip bases). A later refinement, generalized canonical purification (Truflandier, Dianzinga & Bowler, 2016), improved the robustness of density-matrix construction. This LNV/McWeeny approach remains the linear-scaling solver used in the code today.

### 2.4 Basis sets

CONQUEST supports two fundamentally different types of support-function basis:

- **Blip functions (B-splines)**, defined on a regular real-space grid. These form a strictly localised, systematically improvable basis: as the blip grid spacing is reduced, the basis becomes arbitrarily complete, and full plane-wave accuracy can be recovered. This basis was described in the code's founding theoretical papers (Hernández, Gillan & Goringe, *Phys. Rev. B* 1997) and is the natural choice when maximum, systematically convergeable accuracy is required, at higher computational cost, since the support functions themselves must be optimised numerically.
- **Pseudo-atomic orbitals (PAOs)**, numerical atomic-orbital-like functions generated by an accompanying basis-generation utility, conceptually similar to the bases used by SIESTA and PLATO. PAOs are fixed, analytic functions (not variationally reoptimised in the outer loop in the same way blips are), giving much faster, cheaper calculations while still spanning a hierarchy from minimal (single-zeta) ab initio tight-binding bases up to large, near-converged multiple-zeta-polarised sets. CONQUEST can also read `.ion` PAO/pseudopotential files produced by SIESTA.

If a restricted PAO basis is used non-self-consistently, CONQUEST effectively behaves as an ab initio tight-binding code; as the PAO basis is enlarged, or when blips are used, full self-consistent DFT accuracy is recovered.

### 2.5 Multi-site support functions (MSSF)

A major methodological development, introduced by Nakata, Bowler and Miyazaki (2014, 2015), is the **multi-site support function (MSSF)** approach: linear combinations of PAOs on a group of neighbouring atoms are contracted into a smaller number of optimised local orbitals per atom. This drastically reduces the number of basis functions carried through the calculation while largely preserving PAO-level accuracy, and it enables **exact diagonalisation** calculations to be pushed from ~1,000 atoms up to roughly **10,000 atoms** on modest computational resources, without invoking the O(N) solver at all. MSSF accuracy has been validated for bulk Si, Al, Fe, NiO and hydrated DNA fragments, and it has since been reviewed in detail (Nakata, Bowler & Miyazaki, *J. Phys. Soc. Jpn.* 2022).

### 2.6 Pseudopotentials, exchange-correlation, forces and stress

- CONQUEST uses norm-conserving pseudopotentials, and is fully compatible with Don Hamann's ONCVPSP-generated pseudopotentials, including the **PseudoDojo** and **SG15** libraries — the same pseudopotential formalism used by plane-wave codes such as ABINIT and Quantum ESPRESSO/PWSCF, facilitating direct cross-code comparison.
- Exchange-correlation functionals are provided via the **LibXC** library (LDA and GGA fully supported; metaGGA and hybrid/exact-exchange functionals were under active development as of the 2022 code paper).
- **Atomic forces** are computed via a unified algorithmic framework spanning the full hierarchy from empirical tight-binding to full ab initio, self-consistent and non-self-consistent GGA forces (Miyazaki, Bowler, Choudhury & Gillan, 2004; Torralba et al., 2009).
- **Stress and cell optimisation** are implemented for both the diagonalisation and O(N) solvers, allowing full structural relaxation of both ionic positions and the simulation cell.
- **DFT+U** and **constrained DFT (cDFT)**, the latter using a Becke-weight population constraint with analytic forces (Sena, Miyazaki & Bowler, 2011), are available for modelling localised electron states and electron-transfer processes.
- **Dispersion corrections** (semi-empirical DFT-D2/D3, Tkatchenko–Scheffler) and van der Waals density functionals (vdW-DF) are supported, alongside Resta's approach to electric polarisation.

### 2.7 Molecular dynamics

CONQUEST supports NVE, NVT and NPT molecular dynamics with both the diagonalisation and linear-scaling solvers. Reliable MD with O(N) density matrices is more demanding than structural relaxation, because the density matrix must be accurate at *every* MD step to preserve a correct trajectory, not just in a final, well-converged geometry. To make this affordable, CONQUEST implements **extended Lagrangian Born–Oppenheimer molecular dynamics (XLBOMD)** together with the DMM approach (Hirakawa, Suzuki, Bowler & Miyazaki, 2017), avoiding full reconvergence of the density matrix at every step while maintaining long-term energy conservation. This has enabled stable, linear-scaling first-principles MD of complex biological systems with 10,000+ atoms (Arita, Bowler & Miyazaki, 2014; Otsuka, Taiji, Bowler & Miyazaki, 2016).

---

## 3. Parallelisation and Performance

CONQUEST was engineered from the start for **massively parallel** execution:

- **Spatial domain decomposition**: the simulation cell is partitioned into regions, each managed by an MPI process, with atoms and integration-grid points distributed accordingly.
- **Load balancing** is achieved via a space-filling-curve partitioning scheme (Brazdova & Bowler, 2008), ensuring even distribution of atoms/work across processes for irregular or heterogeneous systems.
- **Sparse matrix multiplication** is central to the O(N) solver's performance; a dedicated parallel algorithm for multiplying the sparse Hamiltonian, overlap and density matrices was developed and benchmarked early in the project (Bowler, Miyazaki & Gillan, 2001) and remains foundational to the code's scalability.
- **Demonstrated scaling**: the O(N) solver shows essentially **ideal weak scaling** (fixed number of atoms per MPI process, cores increased in proportion to system size). On the Japanese K computer, near-ideal parallel efficiency was retained using more than **200,000 cores**, and the code has since also been run at large scale on Fugaku. Applications have reached **over 1,000,000 silicon atoms**, and the underlying feasibility of million-atom DFT was first explicitly demonstrated by Bowler & Miyazaki (2010) on the UK's HECToR supercomputer.
- Using the exact-diagonalisation solver with a full PAO basis, systems up to roughly 1,000 atoms are tractable on modest resources (200–500 cores); with MSSF contraction, diagonalisation extends to roughly 10,000 atoms with similar resources.

---

## 4. Practical Scope and Applications

CONQUEST has been applied across a broad range of materials and molecular systems, including:

- **Semiconductor surfaces and nanostructures**: Ge/Si(001) hut-cluster self-assembly (systems up to ~23,000 atoms), Ge(105) reconstructions, Si and SiGe core–shell nanowires.
- **Ferroelectric thin films**: polarisation textures and vortex/domain-wall formation in PbTiO₃ films on SrTiO₃, and interactions between domain walls and surface trenches, using cells of up to 5,000 atoms.
- **Complex oxides**: topologically protected vortices and ferroelectric domain walls in hexagonal YGaO₃.
- **Biomolecular systems**: DNA fragments, the gramicidin A ion channel embedded in a hydrated lipid bilayer (17,000+ atom self-consistent calculations), and platinum-based anticancer drug aquation free-energy profiles.
- **Nanocatalysis**: gold and cobalt-phosphide nanoparticle/nanorod catalysts, metal–organic framework (MOF) electrochromic and gas-storage materials.
- **Interfaces and molecular electronics**: graphene–metal interfaces, single-molecule charge-transport junctions, on-surface coordination chemistry.
- **Materials informatics**: unsupervised-learning analysis of local atomic structure (e.g., in amorphous silica) and machine-learning force-field hyperparameter optimisation, using CONQUEST-generated reference data.

The `MateriApps` materials-simulation portal and the code's own documentation summarise CONQUEST's practical operating regime as follows: for systems of a few hundred atoms, conventional plane-wave codes are often more efficient (CONQUEST's parallel design carries communication overhead that is not worthwhile at small scale); but for systems from several hundred to many millions of atoms, and particularly for highly parallel HPC workloads, CONQUEST is presented as one of the few production DFT codes capable of reaching that regime with full ab initio accuracy.

---

## 5. Comparison with Related Linear-Scaling Codes

CONQUEST belongs to a small family of production O(N) DFT codes that exploit the nearsightedness (locality) of the density matrix in insulating/semiconducting systems:

- **ONETEP** — uses variationally optimised non-orthogonal generalised Wannier functions (NGWFs) on a psinc (plane-wave-derived) basis.
- **BigDFT** — uses a Daubechies wavelet basis to achieve adaptive real-space locality.
- **CONQUEST** — uses strictly localised support functions (blips or PAOs) with direct density-matrix minimisation and sparse-matrix parallelism.

A shared limitation across this class of methods is noted in the literature: because they rely on density-matrix locality, they are best suited to systems with a non-vanishing electronic gap (insulators, semiconductors, large-gap molecular systems); metals and small-gap systems remain challenging for the linear-scaling solver, and the CONQUEST developers have identified this as an important direction for future methodological work.

---

## 6. Availability and Access

- **Source code**: `github.com/OrderN/CONQUEST-release` (MIT licence)
- **Documentation/manual**: `conquest.readthedocs.io` (also distributed with each release)
- **Project pages**: `ordern.github.io`
- **Interoperability**: ASE (Atomic Simulation Environment) integration; reads SIESTA `.ion` files; uses ONCVPSP/PseudoDojo/SG15 pseudopotentials compatible with ABINIT and PWSCF.
- **Current version** (as of the most recent public release referenced): v1.4 (December 2024), with v1.5 documentation subsequently issued.
- CONQUEST issues a machine-generated **BibTeX file** of the specific papers relevant to each run, so users can cite precisely the methodology invoked in a given calculation.

---

## 7. Publications on CONQUEST Theory and Methodology

The list below is restricted to papers describing the **underlying theory, algorithms, and code architecture** of CONQUEST (as opposed to the much larger body of papers reporting scientific *applications* using the code). It is arranged chronologically.

### Foundational papers (1995–2000)

- E. Hernández and M. J. Gillan, "Self-consistent first-principles technique with linear scaling," *Phys. Rev. B* **51**, 10157 (1995).
- E. Hernández, M. J. Gillan and C. M. Goringe, "Linear-scaling density-functional-theory technique: The density-matrix approach," *Phys. Rev. B* **53**, 7147 (1996).
- E. Hernández, M. J. Gillan and C. M. Goringe, "Basis functions for linear-scaling first-principles calculations," *Phys. Rev. B* **55**, 13485 (1997) — description of the B-spline (blip) basis set.
- C. M. Goringe, E. Hernández, M. J. Gillan and I. J. Bush, "Linear-scaling DFT-pseudopotential calculations on parallel computers," *Comput. Phys. Commun.* **102**, 1 (1997).
- M. J. Gillan, D. R. Bowler, C. M. Goringe and E. Hernández, "First Principles Order N Calculations on Very Large Systems," in *The Physics of Complex Liquids*, Proc. Int. Symposium, Nagoya, 1997 (World Scientific, 1998).
- D. R. Bowler and M. J. Gillan, "Length-scale ill conditioning in linear-scaling DFT," *Comput. Phys. Commun.* **112**, 103 (1998).
- D. R. Bowler and M. J. Gillan, "Density matrices in O(N) electronic structure calculations: theory and applications," *Comput. Phys. Commun.* **120**, 95 (1999) — the McWeeny/auxiliary-density-matrix linear-scaling solver still used today.
- D. R. Bowler, I. J. Bush and M. J. Gillan, "Practical methods for ab initio calculations on thousands of atoms," *Int. J. Quantum Chem.* **77**, 831 (2000).

### 2001–2010: parallelisation, forces, PAOs

- D. R. Bowler, T. Miyazaki and M. J. Gillan, "Parallel Sparse Matrix Multiplication for Linear Scaling Electronic Structure Calculations," *Comput. Phys. Commun.* **137**, 255 (2001).
- D. R. Bowler and M. J. Gillan, "An embedding scheme based on quantum linear-scaling methods," *Chem. Phys. Lett.* **355**, 306 (2002).
- D. R. Bowler, T. Miyazaki and M. J. Gillan, "Recent progress in linear scaling ab initio electronic structure techniques," *J. Phys.: Condens. Matter* **14**, 2781 (2002).
- T. Miyazaki, D. R. Bowler, R. Choudhury and M. J. Gillan, "Atomic force algorithms in density functional theory electronic-structure techniques based on local orbitals," *J. Chem. Phys.* **121**, 6186 (2004).
- T. Miyazaki, R. Choudhury, D. R. Bowler and M. J. Gillan, "Large-scale ab-initio calculations," Proc. 3rd Int. Conf. Comput. Model. Simul. Mater. (Techna Group, 2005).
- D. R. Bowler, R. Choudhury, M. J. Gillan and T. Miyazaki, "Recent progress with large-scale ab initio calculations: the CONQUEST code," *Phys. Status Solidi B* **243**, 989 (2006).
- T. Miyazaki, D. R. Bowler, R. Choudhury and M. J. Gillan, "Density functional calculations of Ge(105): Local basis sets and O(N) methods," *Phys. Rev. B* **76**, 115327 (2007).
- M. J. Gillan, D. R. Bowler, A. S. Torralba and T. Miyazaki, "Order-N first-principles calculations with the CONQUEST code," *Comput. Phys. Commun.* **177**, 14 (2007).
- V. Brazdova and D. R. Bowler, "Automatic data distribution and load balancing with space-filling curves: implementation in CONQUEST," *J. Phys.: Condens. Matter* **20**, 275223 (2008).
- A. S. Torralba, M. Todorović, V. Brázdová, R. Choudhury, T. Miyazaki, M. J. Gillan and D. R. Bowler, "Pseudo-atomic orbitals as basis sets for the O(N) DFT code CONQUEST," *J. Phys.: Condens. Matter* **20**, 294206 (2008).
- T. Miyazaki, D. R. Bowler, M. J. Gillan and T. Ohno, "The energetics of hut-cluster self-assembly in Ge/Si(001) from linear-scaling DFT calculations," *J. Phys. Soc. Jpn.* **77**, 123706 (2008).
- T. Otsuka, T. Miyazaki, T. Ohno, D. R. Bowler and M. J. Gillan, "Accuracy of order-N density-functional theory calculations on DNA systems using CONQUEST," *J. Phys.: Condens. Matter* **20**, 294201 (2008).
- A. S. Torralba, D. R. Bowler, T. Miyazaki and M. J. Gillan, "Non-self-consistent Density-Functional Theory Exchange-Correlation Forces for GGA Functionals," *J. Chem. Theory Comput.* **5**, 1499 (2009).
- D. R. Bowler and T. Miyazaki, "Calculations for millions of atoms with density functional theory: linear scaling shows its potential," *J. Phys.: Condens. Matter* **22**, 074207 (2010) — first demonstration of practical million-atom DFT with CONQUEST.

### 2011–2020: constrained DFT, O(N) solver refinements, MSSF, MD, code-paper era

- A. M. P. Sena, T. Miyazaki and D. R. Bowler, "Linear Scaling Constrained Density Functional Theory in CONQUEST," *J. Chem. Theory Comput.* **7**, 884 (2011).
- D. R. Bowler and T. Miyazaki, "O(N) methods in electronic structure calculations," *Rep. Prog. Phys.* **75**, 036503 (2012) — comprehensive review of linear-scaling methods generally, including CONQUEST.
- M. Todorović, D. R. Bowler, M. J. Gillan and T. Miyazaki, "Density-functional theory study of gramicidin A ion channel geometry and electronic properties," *J. R. Soc. Interface* **10**, 20130547 (2013).
- A. Nakata, D. R. Bowler and T. Miyazaki, "Efficient Calculations with Multisite Local Orbitals in a Large-Scale DFT Code CONQUEST," *J. Chem. Theory Comput.* **10**, 4813 (2014) — introduction of the MSSF method.
- M. Arita, D. R. Bowler and T. Miyazaki, "Stable and Efficient Linear Scaling First-Principles Molecular Dynamics for 10,000+ atoms," *J. Chem. Theory Comput.* **10**, 5419 (2014).
- M. Arita, S. Arapan, D. R. Bowler and T. Miyazaki, "Large-scale DFT simulations with a linear-scaling DFT code CONQUEST on K-computer," *J. Adv. Simul. Sci. Eng.* **1**, 87 (2014).
- A. Nakata, D. Bowler and T. Miyazaki, "Optimized multi-site local orbitals in the large-scale DFT program CONQUEST," *Phys. Chem. Chem. Phys.* **17**, 31427 (2015).
- C. O'Rourke and D. R. Bowler, "Linear scaling density matrix real time TDDFT: Propagator unitarity and matrix truncation," *J. Chem. Phys.* **143**, 102801 (2015).
- L. A. Truflandier, R. M. Dianzinga and D. R. Bowler, "Communication: Generalized canonical purification for density matrix minimization," *J. Chem. Phys.* **144**, 091102 (2016).
- T. Otsuka, M. Taiji, D. R. Bowler and T. Miyazaki, "Linear-scaling first-principles molecular dynamics of complex biological systems with the Conquest code," *Jpn. J. Appl. Phys.* **55**, 1102B1 (2016).
- A. Nakata, Y. Futamura, T. Sakurai, D. R. Bowler and T. Miyazaki, "Efficient Calculation of Electronic Structure Using O(N) Density Functional Theory," *J. Chem. Theory Comput.* **13**, 4146 (2017).
- T. Hirakawa, T. Suzuki, D. R. Bowler and T. Miyazaki, "Canonical-ensemble extended Lagrangian Born–Oppenheimer molecular dynamics for the linear scaling density functional theory," *J. Phys.: Condens. Matter* **29**, 405901 (2017).
- C. O'Rourke, S. Y. Mujahed, C. Kumarasinghe, T. Miyazaki and D. R. Bowler, "Structural properties of silicon–germanium and germanium–silicon core–shell nanowires," *J. Phys.: Condens. Matter* **30**, 465303 (2018).
- C. Romero-Muñiz, A. Nakata, P. Pou, D. R. Bowler, T. Miyazaki and R. Pérez, "High-accuracy large-scale DFT calculations using localized orbitals in complex electronic systems: the case of graphene–metal interfaces," *J. Phys.: Condens. Matter* **30**, 505901 (2018).
- D. R. Bowler, J. S. Baker, J. T. L. Poulton, S. Y. Mujahed, J. Lin, S. Yadav, Z. Raza and T. Miyazaki, "Highly accurate local basis sets for large-scale DFT calculations in CONQUEST," *Jpn. J. Appl. Phys.* **58**, 100503 (2019).
- J. S. Baker, T. Miyazaki and D. R. Bowler, "The pseudoatomic orbital basis: electronic accuracy and soft-mode distortions in ABO₃ perovskites," *Electron. Struct.* **2**, 025002 (2020).
- L. A. Truflandier, R. M. Dianzinga and D. R. Bowler, "Notes on density matrix perturbation theory," *J. Chem. Phys.* **153**, 164105 (2020).
- A. Nakata, J. S. Baker, S. Y. Mujahed, J. T. L. Poulton, S. Arapan, J. Lin, Z. Raza, S. Yadav, L. Truflandier, T. Miyazaki and D. R. Bowler, "Large scale and linear scaling DFT with the CONQUEST code," *J. Chem. Phys.* **152**, 164112 (2020) — the primary modern CONQUEST code/methodology reference (arXiv:2002.07704).

### 2021–2025: code paper, MSSF review, exascale roadmap

- D. R. Bowler, T. Miyazaki, A. Nakata and L. Truflandier, "The CONQUEST code: large scale and linear scaling DFT," *Modelling Simul. Mater. Sci. Eng.* (2022) (arXiv:2205.08941) — concise, up-to-date description of code capabilities, solvers, and roadmap.
- A. Nakata, D. R. Bowler and T. Miyazaki, "Large-Scale DFT Methods for Calculations of Materials with Complex Structures," *J. Phys. Soc. Jpn.* **91**, 091011 (2022) — detailed review of the MSSF method.
- V. Gavini *et al.* (including CONQUEST developers), "Roadmap on electronic structure codes in the exascale era," *Modelling Simul. Mater. Sci. Eng.* **31**, 063301 (2023).
- (2025/2026) A stress-calculation methodology paper — "Stress calculation in linear scaling DFT: convergence and dynamics" — extends the exact-diagonalisation and linear-scaling stress formalism within CONQUEST.

---

## 8. Summary

CONQUEST occupies a distinctive niche among DFT packages: it is one of very few production codes able to combine **fully self-consistent, forces-and-stress-capable DFT** with genuine **linear-scaling to millions of atoms**, while also offering a conventional, exact-diagonalisation mode competitive at the thousand-to-ten-thousand-atom scale via multi-site support functions. Its architecture — strictly localised support functions (blip or PAO), sparse-matrix density-matrix minimisation, and a spatially decomposed, near-ideally weak-scaling parallel implementation — has been validated on some of the world's largest supercomputers (K computer, Fugaku, ARCHER), and applied to problems spanning semiconductor nanostructures, ferroelectric thin films, complex oxides, catalytic nanoparticles, and large biomolecular assemblies. Its open MIT licence, ASE integration, and compatibility with standard pseudopotential libraries (PseudoDojo, SG15) make it accessible to the broader electronic-structure community, while its ongoing development — metaGGA/hybrid functionals, improved metallic-system solvers, and continued exascale optimisation — points to an active and expanding role in large-scale materials simulation.

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the CONQUEST 	Linear-scaling DFT code for simulations of very large systems (up to millions of atoms) using numerical atomic orbitals or blip basis functions. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
