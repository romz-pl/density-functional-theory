# Firefly (PC GAMESS): An Exhaustive Review

## 1. Overview

Firefly — known for most of its history as **PC GAMESS** — is a freely available ab initio and density functional theory (DFT) computational chemistry package. It was developed to deliver high performance specifically on Intel-compatible x86, AMD64, and EM64T processors, and it descends from the GAMESS (US) source code base while having evolved into an independently developed program.

The project was launched in 1994 at the Chemistry Department of Moscow State University (MSU) under the coordination of **Dr. Alex A. Granovsky**, with significant contributions over the years from Dr. Anastasia V. Bochenkova and Dr. James W. Kress, among others. The first public release (version 4.0, build #1080) appeared on March 18, 1997.

| Attribute | Detail |
|---|---|
| Original name | PC GAMESS |
| Current name | Firefly (since December 4, 2009) |
| Lead developer | Alex A. Granovsky |
| Institution | Moscow State University, Chemistry Department |
| Initial release | March 18, 1997 |
| Latest stable release | 8.2.0 (September 19, 2016) |
| Beta / RC line | 8.0.0 RC (July 25, 2012) |
| Languages | Fortran, C |
| Operating systems | Windows, Linux, OS X |
| Target platforms | x86, x86-64 (Intel-compatible, AMD64, EM64T) |
| License | Proprietary freeware (registration required) |
| Website | classic.chem.msu.su/gran/gamess/ |

## 2. History and Relationship to GAMESS (US)

- **1994–1999**: The project began as "PC GAMESS," built directly on the GAMESS (US) source code from the Iowa State University Quantum Chemistry Group (ISUQCG), incorporating versions up to October 25, 1999.
- **Post-1999**: No further GAMESS (US) code was incorporated. From this point on, development proceeded completely independently, even though input/output-level compatibility with GAMESS (US) was deliberately preserved as much as possible.
- **2008**: The Firefly Project Team formally dissociated from ISUQCG and GAMESS (US).
- **2009**: The name "PC GAMESS" was retired. From December 4, 2009 onward, only the name **Firefly** is valid, and licenses for pre-Firefly PC GAMESS builds (7.1.B and earlier) were revoked.
- **2012**: Firefly 8.0.0 RC was released for public beta testing, a substantial rewrite offering better speed and reliability than the 7.1.G line.
- **2016**: Firefly 8.2.0, the most recent stable release, was published.
- **2019**: Alex A. Granovsky, the sole principal developer, passed away (born January 15, 1971 – died November 9, 2019), which has effectively frozen further first-party development of the package.

Despite its shared ancestry with GAMESS (US), by the developers' own account roughly **60–70% of the codebase was rewritten**, particularly:
- Platform-specific subsystems (memory allocation, disk I/O, networking)
- Core mathematical routines (e.g., matrix algebra)
- Quantum-chemical method implementations (Hartree–Fock, Møller–Plesset perturbation theory, and DFT)

This rewrite is the primary reason Firefly is significantly faster than the original GAMESS (US) code on the same hardware.

## 3. Core Scientific Capabilities

Firefly implements a broad range of standard and advanced electronic-structure methods:

- **SCF wavefunctions**: RHF, ROHF, UHF, GVB, MCSCF
- **Post-HF correlation**: MP2 (fast analytic energy/gradient), MP3, MP4 (SDQ and full)
- **Multireference methods**: MCSCF (with semi-numerical gradients from state-averaged orbitals and state-tracking), MRMP2, MCQDPT, and Firefly's signature **XMCQDPT2** method
- **Density Functional Theory (DFT)** and **Time-Dependent DFT (TDDFT)**, using a DFT code and functional library that is completely distinct from that of GAMESS (US)
- **Excited-state / conical-intersection tools**: MCSCF-based location of conical intersections and interstate crossings; XMCQDPT2-based potential energy surface mapping
- **Geometry optimization and relaxed surface scans**, with dedicated optimization engines different from GAMESS (US)
- **Vibrational analysis** (Hessian evaluation) for stationary-point characterization

Firefly notably **lacks** some capabilities found in modern GAMESS (US) releases, such as native coupled-cluster (CC) methods and the Fragment Molecular Orbital (FMO) method.

## 4. x86-Oriented Performance Engineering

Firefly's defining engineering goal was maximum throughput on Intel-compatible x86/x86-64 hardware rather than breadth of portability. Techniques cited by the developers and documentation include:

- **Real-time data compression/decompression** for I/O-bound calculation stages
- **Efficient, modern two-electron integral evaluation algorithms** for direct (integral-not-stored) SCF/DFT
- A **Quantum Fast Multipole Method (QFMM)** implementation for accelerating Coulomb-term evaluation, with Firefly-specific routines (e.g., the `NEARJ=2` "fastints"-based near-field routine) distinct from the GAMESS (US)-derived HONDO-based routine
- Highly optimized **MP2 energy and gradient modules**
- A fast **RHF MP3/MP4** energy code
- **Parallel execution** on SMP systems and clusters (or both simultaneously), with attention to scalability on large clusters and many-core systems; later versions (from 6.3 onward) introduced a hierarchical "eXtreme Parallel" (xp) model and peer-to-peer communication layer intended to scale DFT energies/gradients and post-HF corrections (MP2, MCQDPT2) across many nodes, including hybrid MPI/OpenMP execution modes

### Independent speed benchmarking

In the community-run **"Quantum Chemistry Speed Test"** (a public, single-CPU HF/B3LYP benchmark comparing GAMESS (US), Firefly, and Q-Chem, among others), Firefly's DFT implementation placed second overall, behind only the commercial Q-Chem package, while outperforming other free/open packages by a substantial margin. Representative single-CPU timings from that benchmark (B3LYP, default grids per program) showed Firefly completing the test system in well under half the wall-clock time of GAMESS (US) using its own default grid, illustrating the practical effect of Firefly's grid pruning and integral-evaluation optimizations. Because grid defaults differ between programs (Firefly's default DFT grid is a pruned Lebedev grid using 63 radial points per atom with 302 angular points per shell, versus GAMESS (US)'s unpruned 96-radial-point default), direct energy/time comparisons require normalizing the integration grid — which the benchmark authors did as part of the test protocol.

## 5. DFT-Specific Notes

- Firefly's exchange–correlation functional library and numerical-integration (quadrature) machinery are implemented independently of GAMESS (US); the two programs do not necessarily support an identical set of functionals.
- Analytic gradients are available for a range of DFT functionals, including hybrid GGAs, supporting efficient geometry optimization at the DFT level.
- TDDFT is supported for excited-state energies, complementing the multireference (XMCQDPT2/MCQDPT) route for systems where a single-reference DFT description is inadequate (e.g., near conical intersections).
- Firefly is commonly used in the literature for standard workflows such as ground-state geometry optimization at B3LYP with a modest basis set (e.g., 6-31G(d,p) or cc-pVDZ), followed by single-point property or spectral calculations (IR/Raman, UV/Vis via TDDFT) — sometimes cross-validated against XMCQDPT2 results for charge-transfer states where DFT/TDDFT is known to be less reliable.

## 6. XMCQDPT2: Firefly's Signature Multireference Method

The most scientifically distinctive contribution associated with Firefly is **XMCQDPT2** (Extended Multi-Configuration Quasi-Degenerate Perturbation Theory, second order), developed by Granovsky as a reformulation of Nakano's multi-state multi-configuration quasi-degenerate perturbation theory (MCQDPT).

Key characteristics:
- It is a **multi-state, multireference, second-order perturbation theory**, comparable in spirit to multistate CASPT2/MS-CASPT2 and to multireference MP2 (MRMP2).
- It specifically targets known weaknesses of standard MCQDPT and other multistate multireference perturbation theories, such as **erratic behavior of energies near conical intersections and avoided crossings** (the "intruder state" problem), through a state-specific partitioning and an extended model-space formulation.
- A later **resolvent-fitting approximation** was introduced by Granovsky to accelerate MCQDPT2/XMCQDPT2 calculations.
- **Analytical gradient theory** for XMCQDPT2 (including the resolvent-fitted variant) was subsequently developed (partly by external groups, e.g., Park and co-workers), enabling direct geometry optimization and characterization of minimum-energy conical intersections (MECIs) rather than relying solely on numerical gradients.
- An extension, **CAP-XMCQDPT2**, combines XMCQDPT2 with the complex absorbing potential (CAP) technique to compute positions and widths of metastable electronic resonances (shape and Feshbach resonances).

XMCQDPT2 has been widely applied to photochemical and photobiological problems — e.g., retinal chromophore photoisomerization, green fluorescent protein (GFP) chromophore models, photoactive yellow protein (PYP) chromophore isomerization, and stilbene photoisomerization — where multireference treatment of excited-state potential energy surfaces and conical intersections is essential and where the method has generally compared favorably with experiment and with other high-level theoretical approaches.

## 7. Licensing and Distribution

- Firefly is **freeware but not open-source**; the source code is not publicly distributed.
- Obtaining the software requires reading and agreeing to the Firefly License Agreement and submitting basic user/system information through the official registration process.
- Upon registration, users receive a password needed to extract the binaries from a password-protected RAR archive (itself packaged in a ZIP file) available from the official download area.
- Windows, Linux, and (historically) OS X binary distributions have been provided.
- Users are required to cite Firefly explicitly in publications using a specific reference format (e.g., "Alex A. Granovsky, Firefly version 8, http://classic.chem.msu.su"), rather than only mentioning it informally in the text.

## 8. Ecosystem and Related Tools

- **Compatibility layer**: Firefly maintains a high degree of input/output file-format compatibility with GAMESS (US), which allows many existing GAMESS-style workflows, tutorials, and third-party pre/post-processing tools to be reused with only minor adaptation.
- **Graphical front-ends and utilities**: third-party tools such as **Winmostar** and **FiCo (Firefly Commander)**, a batch-job processing utility, provide GUI-based or batch-oriented access to Firefly.
- **Visualization**: Firefly output (geometries, vibrational modes, orbitals) is commonly visualized with external packages such as **MOLDEN** and **MOLEKEL**.
- **Structure conversion**: tools like Open Babel are used to convert common structure formats (e.g., MOL files) into Firefly/PC GAMESS input format.
- **Deployment on HPC clusters**: Firefly has been packaged as an environment module on academic HPC systems (e.g., DESY's Maxwell cluster), typically as a statically linked MPICH1-based parallel build.

## 9. Strengths and Limitations

**Strengths**
- Strong single-node and cluster performance on x86/x86-64 hardware, validated by independent benchmarking against GAMESS (US) and commercial packages.
- A genuinely distinctive, actively cited multireference method (XMCQDPT2) not available in most competing free packages.
- High input/output compatibility with GAMESS (US), easing migration and reuse of existing know-how, tutorials, and scripts.
- No-cost availability for academic and individual use.

**Limitations**
- **Not open source**; source code is unavailable, restricting community contribution, auditability, and modern packaging (e.g., no conda/apt packages).
- Development has been effectively dormant since the death of the sole principal developer, Alex A. Granovsky, in November 2019; the last stable release (8.2.0) dates to 2016.
- Missing some modern method categories present in current GAMESS (US) and other major packages, notably coupled-cluster theory and fragment molecular orbital (FMO) methods.
- Restricted to x86/x86-64 architectures by design, unlike more portable codes.
- Registration/licensing friction (manual approval and password-protected archives) compared to packages with direct or open-source downloads.

## 10. Summary Assessment

Firefly (PC GAMESS) occupies a specific niche in the computational chemistry software landscape: a historically GAMESS (US)-derived, but since 1999 fully independently developed, freeware package whose main value propositions are (1) markedly optimized performance for HF/DFT/MP2–MP4 workloads on ordinary x86 workstations and clusters, and (2) the XMCQDPT2 multistate multireference perturbation theory, a method of real and continuing significance in photochemistry and photobiology research. Its practical relevance today is best understood as twofold: it remains in active use as a computational engine in many published multireference/excited-state studies (via XMCQDPT2 and CASSCF workflows), even as the package itself is no longer under active first-party development following its lead developer's death.

---

# Bibliography: Publications on Firefly's Underlying Theory

The following references document the theoretical and methodological foundations most closely associated with Firefly, particularly its signature multireference perturbation theory developments. Where formal journal citation details were not confirmed in available sources, the software self-citation is given as documented by the developers.

### Primary software citation
1. Granovsky, A. A. *Firefly version 8*, http://classic.chem.msu.su — the official citation required by the Firefly License Agreement for any publication using the program.

### XMCQDPT2 core theory
2. Granovsky, A. A. "Extended multi-configuration quasi-degenerate perturbation theory: The new approach to multi-state multi-reference perturbation theory." *Journal of Chemical Physics*, **134**, 214113 (2011). — The primary methodological paper introducing XMCQDPT2; the Firefly documentation explicitly directs users to read this paper before performing XMCQDPT2 calculations.

### XMCQDPT2 gradient theory
3. Park, J. W. "Analytical Gradient Theory for Resolvent-Fitted Second-Order Extended Multiconfiguration Perturbation Theory (XMCQDPT2)." *Journal of Chemical Theory and Computation*, **17** (10), 6122–6133 (2021).
4. Park, J. W. "Analytical First-Order Derivatives of Second-Order Extended Multiconfiguration Quasi-Degenerate Perturbation Theory (XMCQDPT2): Implementation and Application." *Journal of Chemical Theory and Computation* (2020).

### Extensions of XMCQDPT2
5. Kunitsa, A. A.; Granovsky, A. A.; Bravaya, K. B. "CAP-XMCQDPT2 method for molecular electronic resonances." *Journal of Chemical Physics*, **146** (18), 184107 (2017).

### MCSCF gradient methodology (related Firefly theory)
6. Granovsky, A. A. "Communication: An efficient approach to compute state-specific nuclear gradients for a generic state-averaged multi-configuration self consistent field wavefunction." *Journal of Chemical Physics*, **143** (23), 231101 (2015).

### Illustrative applications establishing/validating XMCQDPT2 theory
7. Granovsky, A. A. "Photoisomerization of Stilbene: The Detailed XMCQDPT2 Treatment." *Journal of Chemical Theory and Computation*, **9** (11), 4973–4990 (2013).
8. Bochenkova, A. V. "Multiconfigurational Methods Including XMCQDPT2 Theory for Excited States of Light-Sensitive Biosystems." Book chapter (2024), https://doi.org/10.1016/B978-0-12-821978-2.00133-1.

### Software/version reference commonly cited alongside theory papers
9. Granovsky, A. A. "Firefly version 8," http://classic.chem.msu.su (accessed via various dates in the literature) — cited in conjunction with the XMCQDPT2 theory papers whenever Firefly is used as the computational engine.

*Note: Item 2 (the 2011 J. Chem. Phys. paper) is the paper Firefly's own documentation identifies as the essential theoretical reference for XMCQDPT2 and should be treated as the primary citation for the method's theory. Full author lists, volume/page details, and additional co-authors for some entries should be verified against the original journal record before use in formal citations, as source snippets did not always confirm complete bibliographic metadata.*


---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the Firefly (PC GAMESS) 	Quantum chemistry package derived from GAMESS(US), optimized for x86 architectures, supporting DFT calculations. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
