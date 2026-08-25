# MOPAC (Molecular Orbital PACkage) — An Exhaustive Review

## 1. Overview

MOPAC is a semi-empirical quantum chemistry program that computes the electronic structure, geometry, energetics, and derived physical properties of molecules, crystals, and nanostructures. It was first released in 1983 and is one of the longest-lived, most widely used codes in computational chemistry. MOPAC occupies the "middle ground" of atomistic simulation, between fully classical molecular-mechanics force fields and expensive *ab initio* electronic-structure methods (Hartree–Fock, post-HF, and DFT), by using quantum-mechanical models whose Hamiltonian matrix elements are simplified and fit ("parameterized") against experimental and high-level theoretical reference data rather than evaluated exactly.

- **Type:** Computational chemistry / electronic-structure software
- **Primary paradigm:** Semi-empirical molecular orbital theory based on the Neglect of Diatomic Differential Overlap (NDDO) approximation
- **Original author:** James J. P. Stewart
- **Current steward:** Molecular Sciences Software Institute (MolSSI), Virginia Tech
- **Written in:** Fortran (98.8% of the codebase), with CMake build scripts, a thin C API layer, and Python bindings/tooling
- **Platforms:** Linux, macOS, Microsoft Windows
- **License:** Apache License 2.0 (open source since 2022)
- **Initial release:** 1983
- **Website:** openmopac.net
- **Repository:** github.com/openmopac/mopac

## 2. What Problem MOPAC Solves

Quantum-chemical simulation is used as a partner to experiment: to interpret results, extend their reach, and sometimes substitute for experiments outright. The obstacle is cost — *ab initio* methods (Hartree–Fock, DFT, post-HF) scale steeply with system size and basis-set quality. Semi-empirical methods were developed specifically to preserve enough quantum-mechanical structure (and thus transferability across chemical space) while cutting computational cost dramatically.

Relative to routine modern DFT calculations, MOPAC-class semi-empirical calculations are reported to be roughly **1000 times faster** while being **about half as accurate**. This trade-off makes MOPAC well suited to:

- High-throughput virtual screening of candidate molecules/materials
- Interactive/teaching applications where near-real-time feedback matters
- Providing a fast pre-check or starting geometry/energy estimate before committing to a much more expensive *ab initio* or DFT calculation
- Large biomolecular systems (proteins, enzymes) where full *ab initio* or DFT treatment is computationally out of reach for most users
- Situations where only modest computing resources (no supercomputer access) are available

It continues to be used, per its current maintainers, as a genuine "middle ground" between *ab initio* quantum mechanics and classical force fields — including in workflows where cheap semi-empirical results are used to accelerate or complement DFT-based pipelines on large systems.

## 3. History and Licensing Timeline

| Period | Status |
|---|---|
| 1981 | Development consolidates prior semi-empirical thermochemistry work from Michael Dewar's group at UT Austin into a unified package |
| 1983 | First public release of MOPAC, distributed as public-domain software via the Quantum Chemistry Program Exchange (QCPE), QCPE program #455 |
| 1990 | Major software paper published describing architecture and features in detail (*J. Computer-Aided Molecular Design*) |
| 1993 | MOPAC acquired by Fujitsu Limited; becomes commercial, closed-source software |
| 1993–2016 | ~30 years as commercial software; final commercial line is "MOPAC 2016" |
| 2022 | Open-source release under the Apache 2.0 license, stewarded by MolSSI as a direct continuation of the codebase |
| 2026 (May) | Most recent stable release at time of writing: v23.2.5 |

The open-source `openmopac/mopac` repository is explicitly described by its maintainers as "a direct continuation of the commercial development and distribution of MOPAC that ended at MOPAC 2016," rather than a rewrite or fork.

## 4. Theoretical Foundations

### 4.1 NDDO and the semi-empirical MO framework
MOPAC's core models share the basic mathematical structure of Hartree–Fock theory evaluated in a minimal valence atomic-orbital basis, but with most of the computationally expensive two-electron integrals either neglected or replaced by simplified, parameterized expressions — the Neglect of Diatomic Differential Overlap (NDDO) approximation. This traces back to early semi-empirical MO theory (Pople, 1953; Pople, Santry & Segal, 1965), refined into the Modified Neglect of Diatomic Overlap (MNDO) model form by Dewar and Thiel (1977), and further developed through the Austin Model 1 (AM1) parameterization (Dewar, Zoebisch, Healy & Stewart, 1985).

### 4.2 The MNDO/AM1/PM-family lineage implemented in MOPAC
MOPAC has historically served as the primary development and reference platform for the MNDO family of Hamiltonians. Successive re-parameterizations and functional-form refinements produced the PMx series (PM3, PM6, PM7), with PM7 (Stewart, 2012) being described in the current JOSS software paper as the flagship model used to illustrate MOPAC's accuracy/cost balance. A dedicated protein/enzyme-optimized model was released more recently (Stewart & Stewart, 2023).

### 4.3 Extension to the periodic table and beyond organic gas-phase chemistry
Originally built for thermochemistry of organic molecules in vacuum, restricted to elements without valence d orbitals, MOPAC's scope has since been extended to:
- Most of the periodic table, via re-optimized NDDO parameters covering roughly 70 elements (Stewart, 2007)
- Solids/crystals, via a practical method for modeling periodic systems with semi-empirical Hamiltonians (Stewart, 2000)
- Molecules in solution, via the COSMO continuum solvation model (Klamt & Schüürmann, 1993)
- Electronic spectroscopy, via a re-release incorporating the INDO/s semi-empirical model with configuration-interaction (CI) excited states (Gieseking, 2021)

### 4.4 Linear-scaling protein/enzyme methodology — MOZYME
Proteins and enzymes commonly contain thousands of atoms, which is prohibitive for standard quantum-chemical self-consistent-field (SCF) procedures whose cost grows steeply with system size. MOPAC addresses this with the **MOZYME** solver, which reformulates the semi-empirical SCF problem in a localized molecular orbital (LMO) basis, giving costs that scale only *linearly* with molecular size (Stewart, 1996). MOZYME can also locate transition states near an enzyme active site embedded within a much larger protein structure (Stewart, 2017), and MOPAC includes practical support for working directly with Protein Data Bank (PDB) structure files, including adding hydrogen atoms that X-ray crystallography does not reliably resolve.

## 5. Software Architecture and Features

- **Primary interface:** command-line program driven by an input file specifying approximate atomic coordinates plus keywords controlling the calculation; produces an output file with heat of formation, optimized geometry, and other computed properties.
- **Application programming interface (API):** since the 23.0.0 release, MOPAC's most commonly used functionality is exposed via an API that avoids disk I/O entirely, addressing a specific bottleneck for high-throughput workloads on systems with many CPU cores but slow filesystems. The API has a C-language binding layer to sidestep Fortran's application binary interface (ABI) portability problems, which in turn underlies a Python wrapper package (`mopactools`) that brings MOPAC into the Python scientific-software ecosystem.
- **Build system:** CMake-based; requires a Fortran compiler, BLAS/LAPACK, Python 3, and NumPy. Optional integration with the MolSSI Driver Interface (MDI) for coupling MOPAC as an engine within multi-code simulation workflows.
- **Distribution:** pre-built binaries for Linux, macOS, and Windows attached to every GitHub Release; also distributed via the `conda-forge` channel (`conda install -c conda-forge mopac`).
- **Third-party integration:** commonly used through graphical front ends such as WebMO, and increasingly as a computational engine driven from Python via its API.
- **Community/governance:** actively maintained by MolSSI, which is supported by U.S. National Science Foundation grant CHE-2136142; developed openly on GitHub with continuous-integration testing and code-coverage tracking.

## 6. Position Relative to Other Software

The JOSS software paper for the open-source release situates MOPAC among several related tools:

- **Other implementations of MNDO-family models:** Gaussian, CP2K, Sparrow, and ULYSSES also implement semi-empirical MNDO-family Hamiltonians, though MOPAC is historically the model-development platform for this family.
- **Ab initio comparators:** Gaussian is cited as the representative example of full *ab initio*/DFT software against which MOPAC's ~1000x speed and ~2x lower accuracy trade-off is framed.
- **High-throughput-focused competitors:** Sparrow is specifically noted as targeting ultra-fast semi-empirical calculations for high-throughput computational campaigns, a niche that overlaps with one of MOPAC's core use cases.

## 7. Practical Use Cases Highlighted by Developers

1. **Interactive chemistry education and exploration** — the low cost per calculation supports near-real-time, interactive quantum-mechanical demonstrations for students.
2. **High-throughput virtual screening** — of molecules/materials, where evaluating very large candidate libraries at *ab initio* cost is infeasible.
3. **Pre-screening before expensive calculations** — using MOPAC results to estimate outcomes and catch problems (bad geometries, convergence issues) before committing to costly DFT or *ab initio* runs, sometimes as an explicit precursor step within mixed semi-empirical/DFT workflows on large systems.
4. **Protein and enzyme modeling** — via MOZYME's linear-scaling LMO-based SCF solver, including active-site transition-state searches within full protein structures, for cases where only the semi-empirical level is computationally affordable to users without supercomputer access.

## 8. Summary Assessment

**Strengths**
- Extremely mature, continuously developed lineage of semi-empirical Hamiltonians (MNDO → AM1 → PM3/PM6/PM7) with MOPAC as their reference implementation
- Orders-of-magnitude speed advantage over *ab initio*/DFT methods, enabling problem sizes (proteins, high-throughput screens) otherwise inaccessible
- Linear-scaling MOZYME solver specifically targets large biomolecular systems
- Broad periodic-table coverage and solid-state/solution-phase extensions beyond its original organic-gas-phase scope
- Now genuinely open source (Apache 2.0) with modern engineering practices (CMake, CI, conda-forge distribution, a disk-free API, Python bindings) after three decades as commercial software
- Actively governed and funded through MolSSI/NSF, reducing single-maintainer risk

**Limitations**
- Roughly half the accuracy of routine modern DFT calculations, per its own developers' characterization — not a substitute for *ab initio*/DFT where high accuracy is required
- Accuracy and transferability are bounded by the quality and scope of the empirical parameterization for a given element/property combination
- Historically closed-source for ~30 years (1993–2022), which limited independent community contribution and auditability during that period
- Best understood as a complementary, cost-efficient screening/pre-processing tool alongside DFT in large-system workflows, rather than a full replacement for it

## 9. Basic Citation

> J. E. Moussa and J. J. P. Stewart, "MOPAC: An open-source semiempirical molecular orbital program," *Journal of Open Source Software*, 2026, 11(119), 8025. DOI: [10.21105/joss.08025](https://doi.org/10.21105/joss.08025)

Software archival citation: DOI [10.5281/zenodo.6511958](https://doi.org/10.5281/zenodo.6511958)

---

# Publications Related to MOPAC's Theory

The following list draws primarily from the reference bibliography of the current open-source MOPAC software paper (Moussa & Stewart, 2026), supplemented with foundational and directly related theory papers it cites. It is organized to trace the theoretical lineage from foundational semi-empirical MO theory through the MNDO/AM1/PM-family Hamiltonians to MOPAC-specific methodological extensions.

### Foundational semi-empirical molecular orbital theory
1. Pople, J. A. (1953). *Electron interaction in unsaturated hydrocarbons.* Transactions of the Faraday Society, 49, 1375–1385. DOI: [10.1039/TF9534901375](https://doi.org/10.1039/TF9534901375)
2. Pople, J. A., Santry, D. P., & Segal, G. A. (1965). *Approximate self-consistent molecular orbital theory. I. Invariant procedures.* The Journal of Chemical Physics, 43, S129–S135. DOI: [10.1063/1.1701475](https://doi.org/10.1063/1.1701475)
3. Dewar, M. J. S. (1975). *Quantum organic chemistry.* Science, 187, 1037–1044. DOI: [10.1126/science.187.4181.1037](https://doi.org/10.1126/science.187.4181.1037)

### The MNDO / AM1 model forms
4. Dewar, M. J. S., & Thiel, W. (1977). *Ground states of molecules. 38. The MNDO method. Approximations and parameters.* Journal of the American Chemical Society, 99, 4899–4907. DOI: [10.1021/ja00457a004](https://doi.org/10.1021/ja00457a004)
5. Dewar, M. J. S., Zoebisch, E. G., Healy, E. F., & Stewart, J. J. P. (1985). *Development and use of quantum mechanical molecular models. 76. AM1: A new general purpose quantum mechanical molecular model.* Journal of the American Chemical Society, 107, 3902–3909. DOI: [10.1021/ja00299a024](https://doi.org/10.1021/ja00299a024)

### MOPAC-specific model development and core software papers
6. Stewart, J. J. P. (1983). *MOPAC, QCPE program #455.* Quantum Chemistry Program Exchange Bulletin, 3, 43.
7. Stewart, J. J. P. (1990). *MOPAC: A semiempirical molecular orbital program.* Journal of Computer-Aided Molecular Design, 4, 1–103. DOI: [10.1007/BF00128336](https://doi.org/10.1007/BF00128336)
8. Stewart, J. J. P. (1996). *Application of localized molecular orbitals to the solution of semiempirical self-consistent field equations* [MOZYME]. International Journal of Quantum Chemistry, 58, 133–146. DOI: [10.1002/(SICI)1097-461X(1996)58:2%3C133::AID-QUA2%3E3.0.CO;2-Z](https://doi.org/10.1002/(SICI)1097-461X(1996)58:2%3C133::AID-QUA2%3E3.0.CO;2-Z)
9. Stewart, J. J. P. (2000). *A practical method for modeling solids using semiempirical methods.* Journal of Molecular Structure, 556, 59–67. DOI: [10.1016/S0022-2860(00)00651-7](https://doi.org/10.1016/S0022-2860(00)00651-7)
10. Stewart, J. J. P. (2007). *Optimization of parameters for semiempirical methods V: Modification of NDDO approximations and application to 70 elements.* Journal of Molecular Modeling, 13, 1173–1213. DOI: [10.1007/s00894-007-0233-4](https://doi.org/10.1007/s00894-007-0233-4)
11. Stewart, J. J. P. (2012). *Optimization of parameters for semiempirical methods VI: More modifications to the NDDO approximations and re-optimization of parameters* [PM7]. Journal of Molecular Modeling, 19, 1–32. DOI: [10.1007/s00894-012-1667-x](https://doi.org/10.1007/s00894-012-1667-x)
12. Stewart, J. J. P. (2017). *An investigation into the applicability of the semiempirical method PM7 for modeling the catalytic mechanism in the enzyme chymotrypsin.* Journal of Molecular Modeling, 23, 154. DOI: [10.1007/s00894-017-3326-8](https://doi.org/10.1007/s00894-017-3326-8)
13. Gieseking, R. L. M. (2021). *A new release of MOPAC incorporating the INDO/s semiempirical model with CI excited states.* Journal of Computational Chemistry, 42, 365–378. DOI: [10.1002/jcc.26455](https://doi.org/10.1002/jcc.26455)
14. Stewart, J. J. P., & Stewart, A. C. (2023). *A semiempirical method optimized for modeling proteins.* Journal of Molecular Modeling, 29, 284. DOI: [10.1007/s00894-023-05695-1](https://doi.org/10.1007/s00894-023-05695-1)
15. Moussa, J. E., & Stewart, J. J. P. (2026). *MOPAC: An open-source semiempirical molecular orbital program.* Journal of Open Source Software, 11(119), 8025. DOI: [10.21105/joss.08025](https://doi.org/10.21105/joss.08025)

### Solvation modeling used within MOPAC
16. Klamt, A., & Schüürmann, G. (1993). *COSMO: A new approach to dielectric screening in solvents with explicit expressions for the screening energy and its gradient.* Journal of the Chemical Society, Perkin Transactions 2, 799–805. DOI: [10.1039/P29930000799](https://doi.org/10.1039/P29930000799)

### Context, comparison, and application literature cited alongside MOPAC's theory
17. Frisch, M. J. et al. (2016). *Gaussian 16, Revision C.01.* Gaussian, Inc. (comparator *ab initio*/DFT package also implementing MNDO-family models)
18. Kühne, T. D. et al. (2020). *CP2K: An electronic structure and molecular dynamics software package – Quickstep: Efficient and accurate electronic structure calculations.* The Journal of Chemical Physics, 152(19), 194103. DOI: [10.1063/5.0007045](https://doi.org/10.1063/5.0007045)
19. Bosia, F., Zheng, P., Vaucher, A., Weymuth, T., Dral, P. O., & Reiher, M. (2023). *Ultra-fast semi-empirical quantum chemistry for high-throughput computational campaigns with Sparrow.* The Journal of Chemical Physics, 158(5), 054118. DOI: [10.1063/5.0136404](https://doi.org/10.1063/5.0136404)
20. Menezes, F., & Popowicz, G. M. (2022). *ULYSSES: An efficient and easy to use semiempirical library for C++.* Journal of Chemical Information and Modeling, 62(16), 3685–3694. DOI: [10.1021/acs.jcim.2c00757](https://doi.org/10.1021/acs.jcim.2c00757)
21. Weymuth, T., & Reiher, M. (2021). *Immersive interactive quantum mechanics for teaching and learning chemistry.* Chimia, 75, 45–49. DOI: [10.2533/chimia.2021.45](https://doi.org/10.2533/chimia.2021.45)
22. Pyzer-Knapp, E. O., Suh, C., Gómez-Bombarelli, R., Aguilera-Iparraguirre, J., & Aspuru-Guzik, A. (2015). *What is high-throughput virtual screening? A perspective from organic materials discovery.* Annual Review of Materials Research, 45, 195–216. DOI: [10.1146/annurev-matsci-070214-020823](https://doi.org/10.1146/annurev-matsci-070214-020823)
23. Moussa, J. E. (2025). *The enduring relevance of semiempirical quantum mechanics.* The Journal of Physical Chemistry A, 129, 8465–8477. DOI: [10.1021/acs.jpca.5c03425](https://doi.org/10.1021/acs.jpca.5c03425)
24. Burley, S. K., Berman, H. M., Kleywegt, G. J., Markley, J. L., Nakamura, H., & Velankar, S. (2017). *Protein Data Bank (PDB): The single global macromolecular structure archive.* In *Protein Crystallography: Methods and Protocols* (pp. 627–641). Springer New York. DOI: [10.1007/978-1-4939-7000-1_26](https://doi.org/10.1007/978-1-4939-7000-1_26)
25. Boyd, D. B. (2013). *Quantum Chemistry Program Exchange, facilitator of theoretical and computational chemistry in pre-internet history.* In *Pioneers of Quantum Chemistry* (pp. 221–273). American Chemical Society. DOI: [10.1021/bk-2013-1122.ch008](https://doi.org/10.1021/bk-2013-1122.ch008)

---

*Sources: MOPAC JOSS software paper (Moussa & Stewart, 2026, DOI 10.21105/joss.08025) and its bibliography; the official `openmopac/mopac` GitHub repository and README; MolSSI's MOPAC project page; Wikipedia's MOPAC entry (infobox facts cross-checked). Compiled August 25, 2026.*


---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the MOPAC 	Primarily semi-empirical quantum chemistry package; re-released as open source, sometimes used alongside DFT workflows for large systems. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
