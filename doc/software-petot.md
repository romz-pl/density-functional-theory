# PEtot: A Review

## 1. Overview

**PEtot** ("Parallel total Energy") is a portable, parallel plane‑wave pseudopotential (PWP) density functional theory (DFT) code for atomistic total‑energy and electronic‑structure calculations, developed primarily by **Lin‑Wang Wang** and collaborators at the Computational Research Division / Materials Sciences Division of **Lawrence Berkeley National Laboratory (LBNL)**, with support from the U.S. Department of Energy (DOE). It is explicitly designed for **large‑system simulations on large parallel computers** — historically IBM SP machines at NERSC, and subsequently Linux/Cray/Blue Gene clusters, and GPU‑accelerated clusters. The code is distributed as free, open source under an **LBNL BSD‑style license**.

PEtot solves the Kohn–Sham equations self‑consistently using a plane‑wave basis set combined with norm‑conserving and ultrasoft pseudopotentials, following the standard total‑energy pseudopotential methodology pioneered in the Car–Parrinello/iterative‑minimization tradition. Its distinguishing contribution is not novel DFT physics per se, but rather a set of **parallelization and numerical‑implementation innovations** that let PWP‑DFT scale to unusually large atom counts and processor counts for its era.

## 2. Development History and Versions

PEtot exists in three released versions, reflecting an evolution in scale and parallelization sophistication:

| Version | Approx. size | Notes |
|---|---|---|
| Version 1 | ~0.5 MB source | Original release |
| Version 2 | ~4 MB source | Documented on the public PEtot homepage; single/simple parallelization scheme |
| Version 3 | ~50 MB source | Current recommended version; adds three levels of parallelization and all‑band algorithms |

Per the code's own documentation, Version 3 introduced **three simultaneous levels of parallelization** — over G‑vectors (plane‑wave coefficients), over the band (state) index, and over k‑points — together with **all‑band conjugate‑gradient algorithms**, making it markedly faster and more scalable for large systems than the earlier versions. The last public update to the homepage tarballs is dated **June 18, 2010** ("to reduce memory in Hpsi").

## 3. Underlying Physical/Numerical Methodology

- **Basis set:** Plane waves (reciprocal‑space representation of Kohn–Sham orbitals), the standard choice for periodic solids and supercells.
- **Pseudopotentials:** Both **norm‑conserving** pseudopotentials (Martins‑type libraries) and **ultrasoft** pseudopotentials (Vanderbilt‑type libraries) are supported; the homepage links directly to J. L. Martins's and D. Vanderbilt's pseudopotential distribution sites.
- **Nonlocal projector evaluation:** PEtot implements Lin‑Wang Wang's **mask‑function real‑space method** for evaluating the nonlocal (projector) part of the pseudopotential in real space rather than reciprocal space, which is significantly cheaper for large systems while preserving high accuracy.
- **Eigensolver / SCF minimization:** Iterative, conjugate‑gradient‑based minimization of the Kohn–Sham energy functional (in the spirit of Payne–Teter–Allan–Arias‑type "iterative minimization techniques" widely used in PWP codes), including all‑band (blocked) conjugate‑gradient variants in version 3 that better exploit parallel dense linear algebra and communication patterns.
- **Charge‑density mixing:** Uses Thomas–Fermi‑based charge‑mixing schemes for accelerating self‑consistency (see Raczkowski, Canning & Wang below).
- **Parallel FFTs:** Custom parallel 3‑D FFT implementations (work led by Andrew Canning) are central to PEtot's scalability, since FFTs between real and reciprocal space dominate communication cost in PWP‑DFT.

## 4. Parallelization Architecture

PEtot's core technical identity is its **multi‑level MPI parallelization**, designed to let plane‑wave DFT — normally hard to scale past a few thousand cores — run efficiently on tens to hundreds of thousands of processors:

1. **G‑space (plane‑wave/FFT‑grid) parallelization** — distributing reciprocal‑space coefficients and real‑space FFT grids across processors.
2. **Band‑index (state) parallelization** — distributing the set of occupied/computed Kohn–Sham orbitals ("bands") across processor groups.
3. **k‑point parallelization** — distributing independent Brillouin‑zone sample points across processor groups (available from version 3 onward).

This hybrid G‑vector/band‑index/k‑point scheme allows PEtot to calculate systems of up to roughly **a thousand atoms on thousands of processors**, and — in derivative divide‑and‑conquer schemes built on top of PEtot (see §6) — to be used as the per‑fragment DFT "engine" for simulations of tens of thousands of atoms.

## 5. Performance Engineering and Hardware Ports

Because PEtot is a DOE/NERSC‑lineage code, it has repeatedly served as a testbed for HPC performance‑engineering research:

- **Parallel FFT communication studies**: Canning, Shalf, Wang, Wasserman & Gajbe (2009) compared different communication structures for scalable parallel 3‑D FFTs in first‑principles codes including PEtot.
- **GPU porting (2011–2013)**: Wang and Jia and collaborators ported PEtot to multi‑node GPU clusters:
  - An initial GPU implementation (SC'11 poster/paper, 2011) achieved roughly a **10× speed‑up** over the CPU version using a hybrid reciprocal‑space/band‑index parallelization scheme, benchmarked on up to 256 CPU/GPU units, and was described as the first GPU DFT‑PWP code demonstrated scalable to large numbers of CPU/GPU computing units.
  - A follow‑up, more deeply optimized GPU implementation (Jia et al., *Computer Physics Communications*, 2013) moved essentially the entire calculation onto the GPU, introduced new algorithms to reduce MPI communication data volume, and adopted new GPU/CPU numerical libraries, achieving **13×–22× speed‑ups** over CPU PEtot for a representative 512‑atom system, with detailed scaling analysis versus GPU/CPU unit count.
- **Commercial/derivative code**: The GPU‑era PEtot lineage evolved into **PWmat**, a commercial plane‑wave GPU materials‑simulation code explicitly described in the literature as "developed based on a CPU DFT plane‑wave pseudopotential code: PEtot," with reported speed‑ups of 18×–30× over CPU PEtot.
- PEtot is regularly cited as a **reference/comparison code** in the DFT software literature — e.g. in surveys of plane‑wave DFT packages alongside VASP, CASTEP, CPMD, ABINIT, Quantum ESPRESSO (PWSCF), DACAPO, SOCORRO, DFT++, PARATEC, CP2K, SPHINX, and Qbox — and as a large‑system baseline in recent (2024) work extending plane‑wave DFT to Sunway‑class exascale supercomputers (PWDFT‑SW), which explicitly groups PEtot with PWDFT as one of the few PWP‑DFT codes "designed for large‑size atoms (thousands of atoms)" as opposed to codes tuned for hundreds of atoms.

## 6. Role as a Building Block: LS3DF and Divide‑and‑Conquer Simulations

PEtot's most consequential downstream impact is as the **direct-DFT solver kernel inside the Linearly Scaling Three‑Dimensional Fragment (LS3DF) method**, developed by Lin‑Wang Wang, Zhengji Zhao, Juan Meza, Byounghak Lee, Hongzhang Shan and collaborators at LBNL:

- LS3DF is an O(N) (linear‑scaling), divide‑and‑conquer electronic‑structure method that decomposes a large system into overlapping small fragments, solves each fragment with a **direct LDA/DFT code — explicitly identified as PEtot** — and then patches the fragment results together with a scheme that cancels artificial boundary (surface) effects, reproducing essentially the same accuracy as a full direct DFT calculation.
- The LS3DF paper (Wang, Lee, Shan, Zhao, Meza, Strohmaier & Bailey, *Proceedings of SC08*, 2008) won the **2008 ACM Gordon Bell Prize** ("Special" category, for algorithmic innovation), reporting up to 60.3 Tflop/s (23.4% of peak) on 30,720 Cray XT4 cores and 107.5 Tflop/s (24.2% of peak) on 131,072 BlueGene/P cores.
- Later LS3DF performance papers report scaling to **442 Tflop/s on 147,456 processors** (Cray XT5 "Jaguar," OLCF) and runs on up to 163,840 BlueGene/P processors (ALCF), applied to systems of up to 36,000 atoms — with the per‑fragment solves still performed by "PEtot_F," a fragment‑adapted variant of PEtot.
- A 2024 paper ("10‑Million Atoms Simulation of First‑Principle Package LS3DF," *JCST*) reports further-optimized LS3DF/PEtot‑lineage code scaling to a **10‑million‑silicon‑atom system** at 34.8 PFLOPS (21.2% of peak) on the Sugon supercomputer, underscoring the continued relevance of the PEtot computational kernel decades after its initial release.
- PEtot/LS3DF‑based simulations have also been used for INCITE‑allocated production science (e.g., large CdSe/CdS nanorod and ZnTeO alloy electronic‑structure studies).

## 7. Typical Use Cases

Based on its design goals and the literature that cites/uses it, PEtot (directly or via LS3DF/PWmat) has been applied to:

- Large semiconductor and nanostructure supercells (quantum dots, nanorods, alloy systems such as ZnTe₁₋ₓOₓ)
- Total‑energy, self‑consistent electronic‑structure calculations of systems from tens to thousands of atoms directly, and up to tens of thousands–millions of atoms via the LS3DF divide‑and‑conquer extension
- Benchmarking and methodology papers on parallel FFTs, GPU acceleration, and exascale‑class plane‑wave DFT performance engineering

## 8. Availability

- **Homepage:** the PEtot homepage (hosted historically at LBNL's CMSN/CRD site, e.g. `cmsn.lbl.gov/html/PEtot/PEtot.html`, and mirrored at `hpcrd.lbl.gov/~linwang/PEtot/PEtot.html`) hosts source tarballs for all three versions, input/output file documentation (`etot.input`, `atom.config`, `vwr.atom`, `kpt.file`, `symm.file`, `kpgen.input`, `moment.input`, `report`, `pmatrix`, `maskr`, `graph.j`), links to norm‑conserving and ultrasoft pseudopotential libraries, and a downloadable mask‑function package for the real‑space nonlocal‑pseudopotential implementation.
- **License:** LBNL BSD license — free to use and modify for any purpose, with no warranty from the developers.
- **Primary developer/contact listed:** Lin‑Wang Wang, Computational Research Division (CRD), LBNL.

## 9. Summary Assessment

PEtot occupies a specific niche in the plane‑wave DFT software ecosystem: rather than competing head‑to‑head with feature‑rich, widely‑used packages such as VASP, Quantum ESPRESSO, or CASTEP on breadth of functionality, its primary value proposition has been **raw parallel scalability for large atom counts on large processor counts**, achieved through (a) real‑space mask‑function evaluation of nonlocal pseudopotentials, (b) custom parallel FFTs, (c) a three‑tier (G‑space/band/k‑point) MPI parallelization scheme with all‑band algorithms, and (d) early and aggressive GPU porting. Its greatest scientific and technical legacy is arguably indirect: serving as the per‑fragment DFT engine inside the Gordon Bell Prize‑winning LS3DF divide‑and‑conquer method, and as the algorithmic ancestor of the commercial GPU code PWmat — both of which have pushed plane‑wave DFT simulation sizes from thousands toward millions of atoms on leadership‑class supercomputers.

---

# Publications Related to PEtot's Theory and Methodology

### Core PEtot method / numerical implementation
1. Wang, L.-W. **"Mask‑function real‑space implementations of nonlocal pseudopotentials."** *Phys. Rev. B* **64**, 201107(R) (2001). — The foundational real‑space nonlocal‑pseudopotential (mask‑function) formalism underlying PEtot's evaluation of the nonlocal projector terms.
2. Canning, A., Wang, L.-W., Williamson, A., Zunger, A. **"Parallel empirical pseudopotential electronic structure calculations for million atom systems."** *J. Comput. Phys.* **160**, 29–41 (2000). https://doi.org/10.1006/jcph.2000.6440 — Parallel folded‑spectrum empirical‑pseudopotential precursor establishing the parallelization philosophy later carried into PEtot's ab initio (DFT) framework.
3. Raczkowski, D., Canning, A., Wang, L.-W. **"Thomas–Fermi charge mixing for obtaining self‑consistency in density functional calculations."** *Phys. Rev. B* **64**, R121101 (2001). — Charge‑density mixing scheme used for SCF convergence in PEtot.
4. Canning, A., Shalf, J., Wang, L.-W., Wasserman, H., Gajbe, M. **"A comparison of different communication structures for scalable parallel three dimensional FFTs in first principle codes."** *Proc. ParCo09*, Lyon, France (2009). — Analysis of the parallel FFT communication schemes central to PEtot's scalability.
5. Payne, M. C., Teter, M. P., Allan, D. C., Arias, T. A., Joannopoulos, J. D. **"Iterative minimization techniques for ab initio total‑energy calculations: molecular dynamics and conjugate gradients."** *Rev. Mod. Phys.* **64**, 1045–1097 (1992). — Foundational iterative‑diagonalization/conjugate‑gradient methodology for plane‑wave pseudopotential total‑energy calculations, of which PEtot's SCF solver is a parallel implementation.

### GPU acceleration of PEtot
6. Jia, W. et al. **"Large scale plane wave pseudopotential density functional theory calculations on GPU clusters."** *Proceedings of SC11 (2011 International Conference for High Performance Computing, Networking, Storage and Analysis)*, Article No. 71 (2011). https://doi.org/10.1145/2063384.2063479
7. Jia, W., Fu, J., Cao, Z., Wang, L., Chi, X., Gao, W., Wang, L.-W. **"The analysis of a plane wave pseudopotential density functional theory code on a GPU machine."** *Comput. Phys. Commun.* **184**, 9–18 (2013). https://doi.org/10.1016/j.cpc.2012.08.002 — 13×–22× GPU speed‑up analysis of PEtot.
8. Related: Jia, W. et al. **"Fast plane wave density functional theory molecular dynamics calculations on multi‑GPU machines."** *J. Comput. Phys.* (cited as a follow‑on GPU‑MD performance paper in the LS3DF/PEtot lineage).

### LS3DF (divide‑and‑conquer method built on PEtot)
9. Wang, L.-W., Zhao, Z., Meza, J. **"Linear scaling three‑dimensional fragment method for large‑scale electronic structure calculations."** *Phys. Rev. B* **77**, 165113 (2008).
10. Zhao, Z., Meza, J., Wang, L.-W. **"A divide and conquer linear scaling three dimensional fragment method for large scale electronic structure calculations."** *J. Phys.: Condens. Matter* **20**, 294203 (2008).
11. Wang, L.-W., Lee, B., Shan, H., Zhao, Z., Meza, J., Strohmaier, E., Bailey, D. H. **"Linearly scaling 3D fragment method for large‑scale electronic structure calculations."** *Proceedings of SC08 (ACM/IEEE Conference on Supercomputing)* (2008); LBNL‑959E. **2008 ACM Gordon Bell Prize (Special Category).** https://doi.org/10.1145/1413370.1413437
12. Wang, L.-W. et al. **"The linearly scaling 3D fragment method for large scale electronic structure calculations."** *J. Phys.: Conf. Ser.* **180**, 012079 (2009). https://doi.org/10.1088/1742-6596/180/1/012079 — Extended performance results (442 Tflop/s on 147,456 cores) and CdSe/CdS nanorod application.
13. Yan, Y.-J. et al. **"10‑Million Atoms Simulation of First‑Principle Package LS3DF on Sugon Supercomputer."** *J. Comput. Sci. Technol.* (2024). — Modern re‑optimization of the LS3DF/PEtot kernel reaching 10‑million‑atom silicon simulations.

### Related plane‑wave/pseudopotential DFT context and comparisons
14. Bylaska, E. J., Tsemekhman, K., Baden, S. B., Weare, J. H., Jónsson, H. **"Parallel implementation of γ‑point pseudopotential plane‑wave DFT with exact exchange."** *J. Comput. Chem.* **32**, 54–69 (2011). https://doi.org/10.1002/jcc.21598 — Related parallel PWP‑DFT algorithm work (NWChem), frequently discussed alongside PEtot in the PWP‑DFT literature.
15. Genovese, L., Neelov, A., Goedecker, S., Deutsch, T., Ghasemi, S. A., Willand, A., Galiste, D., Zilberberg, O., Rayson, M., Bergman, A., Schneider, R. **"Daubechies wavelets as a basis set for density functional pseudopotential calculations."** *J. Chem. Phys.* **129**, 014109 (2008). — Cited comparison point in PEtot‑GPU papers among alternative DFT bases.
9. (Cross‑reference) Chen, Y. et al., **"PWDFT‑SW: Extending the Limit of Plane‑Wave DFT Calculations to 16K Atoms on the New Sunway Supercomputer,"** arXiv:2406.10765 (2024) — modern survey explicitly situating PEtot among large‑system PWP‑DFT codes.

### PEtot software homepage / primary documentation
16. Wang, L.-W. **PEtot Homepage and source distribution** (versions 1–3), Lawrence Berkeley National Laboratory, Computational Research Division. Last updated June 18, 2010. Available (historically) at `cmsn.lbl.gov/html/PEtot/PEtot.html` and `hpcrd.lbl.gov/~linwang/PEtot/PEtot.html`.

---

*Note: Several early source URLs for the PEtot homepage (LBNL CMSN/HPCRD web space) are no longer reliably reachable; content above was reconstructed from an archived snapshot of the homepage plus the peer‑reviewed literature that describes, uses, or cites PEtot. Readers seeking the current download location should search LBNL's Computational Research Division / Lin‑Wang Wang's group pages, as institutional URLs for legacy codes are prone to being moved or retired.*

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the PEtot 	Portable, plane-wave pseudopotential density functional theory (DFT) computer program designed for large-scale electronic structure calculations on parallel computer architectures. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
