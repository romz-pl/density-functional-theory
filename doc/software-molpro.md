# MOLPRO — Exhaustive Software Review

*Commercial quantum chemistry package specializing in high-accuracy correlated wavefunction methods, with integrated DFT/TD-DFT capability*

---

## 1. Overview

MOLPRO is a general-purpose, commercially licensed *ab initio* quantum chemistry program package developed and maintained by **Hans-Joachim Werner** (Universität Stuttgart) and **Peter J. Knowles** (Cardiff University), together with a large community of contributing authors. Its development history stretches back over 50 years, originating in the late 1960s from Wilfried Meyer's pseudo-natural-orbital configuration-interaction (PNO-CI) and coupled electron-pair (PNO-CEPA) methods, and Peter Pulay's pioneering analytical Hartree–Fock gradient program.

The program's defining identity is **accuracy-first electron correlation treatment**: MOLPRO is one of the reference implementations for internally contracted multireference configuration interaction (IC-MRCI) and for explicitly correlated (F12) coupled-cluster and multireference methods, and it is widely regarded as a gold-standard tool for benchmark-quality thermochemistry, spectroscopy, and potential energy surface (PES) generation. Efficient DFT and TD-DFT functionality is included as a complementary, production-quality capability rather than the program's primary focus.

MOLPRO is closed-source and requires a paid licence (with discounted academic rates and free trial licences); it is distributed by **TTI GmbH** (Stuttgart, Germany).

---

## 2. Development History and Governance

| Era | Milestone |
|---|---|
| Late 1960s | Meyer develops PNO-CI/PNO-CEPA; Pulay develops the first analytical HF energy gradient program |
| 1970s | Methods applied to unprecedented-accuracy molecular property calculations for small molecules |
| 1980 | Werner and Meyer develop a new state-averaged MCSCF method, laying groundwork for later multireference developments |
| 1980s–1990s | Development of internally contracted MRCI (IC-MRCI), enabling global potential energy surfaces for ground and excited electronic states |
| 2000s | Local correlation methods, density fitting, explicitly correlated (F12) methods introduced |
| 2010s | PNO-based local coupled-cluster methods extend explicitly correlated CCSD(T)-F12 accuracy to systems of ~100 atoms |
| 2020 | Comprehensive package paper published in *J. Chem. Phys.* (the current primary citation) |
| 2020s | Python scripting framework and new graphical front end **iMolpro** introduced; PNO-LCCSD(T)-F12 refined for noncovalent interactions and conformational/isomerization energies, extending routine applicability to 100–200 atom systems |
| June 2026 | Version **2026.1** released (current release at time of writing) |

Governance and IP rest with TTI GmbH; the scientific/technical leadership remains with Werner and Knowles, supported by a large, named developer community (full author list maintained at molpro.net).

---

## 3. Licensing and Distribution

- **Model:** Commercial, closed-source, node/group/site/service licence tiers, typically sold in 1–4-year terms.
- **Academic pricing (indicative, per year, EUR, non-commercial only):**
  - Single-workstation (≤16 cores): ~€550 first year / ~€500 renewal
  - Single-workstation (≤32 cores): ~€850 first year / ~€750 renewal
  - Research-group licence: ~€1,700 first year / ~€1,500 renewal
  - Site licence (whole institution): ~€5,100 first year / ~€4,500 renewal
  - Service licence (computer-centre, all registered users, source code included): ~€6,800 first year / ~€6,000 renewal
- **Commercial-use pricing:** provided on request (not publicly listed).
- **Trials:** 1-month free trial licences available for Linux/macOS.
- **Developing-country programme:** discounted/subsidised access scheme for researchers in selected countries.
- **Source code:** binaries are standard; source access is automatic for service-licence holders and available on request for group/site licences, primarily to allow optimized builds on specialist HPC architectures.
- **Export restrictions:** because MOLPRO bundles Intel-licensed libraries, it cannot legally be distributed to Russia, Cuba, Iran, North Korea, Sudan, Syria, or parts of Ukraine under U.S. Export Administration Regulations.
- **Platforms:** Linux and macOS (per current stable release metadata); historically also supported on a range of parallel/HPC architectures.

---

## 4. Core Scientific Capabilities

### 4.1 Wavefunction-based electron correlation methods (flagship strength)
- **Multireference methods:** State-averaged CASSCF/MCSCF (second-order convergence algorithms originating from Werner–Knowles and Knowles–Werner developments), internally contracted MRCI (IC-MRCI), CASPT2, and their explicitly correlated (F12) analogues (CASPT2-F12, MRCI-F12) — these are historically MOLPRO's most distinctive and widely cited capability, enabling global, spectroscopically accurate potential energy surfaces for ground and excited electronic states.
- **Single-reference coupled-cluster methods:** Highly efficient CCSD, CCSD(T), and explicitly correlated CCSD(T)-F12 implementations, including density-fitted/integral-direct variants.
- **Explicit correlation (F12) methods:** A central MOLPRO specialty — F12-corrected wavefunction methods closely approach the complete basis set (CBS) limit already at triple-zeta basis quality, substantially reducing basis-set error and computational cost relative to conventional CBS extrapolation.
- **Local correlation / PNO methods:** Local coupled-cluster approaches (e.g., PNO-LCCSD(T)-F12) reduce the formal scaling of correlated calculations with molecular size, enabling near-CCSD(T)/CBS-quality energies for molecules of roughly 100–200 atoms — a scale historically inaccessible to canonical coupled-cluster methods.
- **Perturbation theory:** MP2 (including local and explicitly correlated variants), CASPT2, and related multireference perturbation approaches.
- **Composite/embedding techniques:** Projection-based embedding methods for treating a chemically active subsystem at high level within a larger environment.

### 4.2 Density Functional Theory (secondary but production-grade)
- Efficient Kohn–Sham DFT with a large library of exchange-correlation functionals.
- TD-DFT for excited-state properties.
- DFT is positioned in MOLPRO's own documentation as a complementary capability ("efficient DFT and TD-DFT methods with a very large number of density functionals are also available") rather than the program's differentiating strength — that role is filled by the correlated wavefunction methods above.

### 4.3 Properties, gradients, and relativistic treatments
- Analytical energy gradients for a wide range of closed- and open-shell methods, for both ground and excited electronic states.
- A broad range of molecular property calculations (e.g., NMR shielding tensors, magnetizabilities, rotational g-tensors via density-fitted local MP2/GIAO approaches).
- Anharmonic vibrational spectra calculations.
- Relativistic effects via Douglas–Kroll–Hess (DKH) scalar-relativistic Hamiltonians, Breit–Pauli spin–orbit treatments, and effective core potentials (ECPs).
- Basis sets: full support for correlation-consistent (cc-pVXZ, aug-cc-pVXZ, cc-pwCVXZ, etc.) and other standard Gaussian basis families, plus complementary auxiliary basis sets for density fitting/RI.

### 4.4 Interfaces and usability
- **Input:** traditional MOLPRO command-based input language.
- **iMolpro:** a newer graphical front end (now available on Windows, macOS, and Linux) providing a complete environment for building, submitting (locally or to remote servers), and analysing calculations; open source and freely distributed (requires a licensed MOLPRO installation to run).
- **gmolpro:** the older graphical user interface, still functional but considered obsolescent in favor of iMolpro.
- **Python framework:** calculations can be specified and analysed programmatically via a Python interface, in addition to conventional input files.
- **Extensibility:** general interfaces exist for coupling MOLPRO to external codes (e.g., dynamics packages such as SHARC use MOLPRO as an electronic-structure back end for non-adiabatic molecular dynamics).

---

## 5. Typical Application Domains

- Benchmark-quality thermochemistry (heats of formation, reaction/isomerization/conformational energetics)
- High-resolution rovibrational spectroscopy and global PES construction for small-to-medium molecules
- Non-adiabatic/excited-state dynamics (via multireference methods and dynamics-code interfaces such as SHARC)
- Non-covalent interaction energies (a specific, recently emphasized use case for PNO-LCCSD(T)-F12)
- Astrochemistry and molecular collision studies (e.g., ro-vibrational quenching cross-sections, interstellar reaction kinetics)
- Cold-molecule and ultracold-collision physics (e.g., alkali-dimer and alkaline-earth-monofluoride interaction potentials)
- Transition-metal and larger-molecule electronic structure via local/PNO correlation methods and DFT

---

## 6. Strengths and Limitations

### Strengths
- Best-in-class, extensively validated implementations of IC-MRCI, CASPT2, and explicitly correlated (F12) coupled-cluster/multireference methods.
- Local/PNO correlation methods extend genuinely high-accuracy (near-CBS CCSD(T)) treatment to molecules of ~100–200 atoms, competitive with or exceeding other local-correlation codes (e.g., ORCA's DLPNO family) for the specific niche of high-accuracy energetics on larger systems.
- Long, continuous development lineage with deep expertise in multireference and explicitly correlated theory.
- Mature analytical gradient and property infrastructure across many correlated methods, not just DFT.
- Flexible interfacing (Python framework, iMolpro, external dynamics-code coupling).

### Limitations
- Closed-source and commercially licensed — a barrier relative to free/open packages (e.g., Psi4, PySCF) for cost-sensitive groups, though academic and trial pricing is comparatively modest for smaller licence tiers.
- DFT functionality, while efficient and broad in functional coverage, is not MOLPRO's differentiating feature and is comparatively less emphasized in its own documentation than the correlated wavefunction methods.
- Export-controlled distribution excludes several countries.
- Source-code access is gated by licence tier, limiting deep customization for many academic users relative to fully open-source alternatives.
- Learning curve for the native input language and advanced multireference/local-correlation options can be steep for new users, though iMolpro substantially lowers this barrier.

---

## 7. Version and Release Notes (Recent)

- **2026.1** (released June 2026) — current stable release; includes new features and improvements summarized in the official "recent changes" documentation.
- **2025.2** — prior stable release cited in general reference sources.
- Historical versions (2012.1, 2015, 2019.2, 2022.1, 2022.2, etc.) are frequently cited in the primary literature, reflecting the package's long, incrementally documented release history.

---

## 8. Summary Assessment

MOLPRO occupies a distinctive niche among quantum chemistry packages: rather than competing primarily on DFT breadth or general-purpose convenience, it is the tool of choice when the scientific requirement is **highest attainable accuracy** in correlated wavefunction theory — particularly multireference problems (bond breaking, excited states, near-degenerate electronic structure) and explicitly correlated coupled-cluster energetics at or near the complete-basis-set limit. Its DFT/TD-DFT module makes it usable as a general-purpose package for routine work, but the software's reputation and core value proposition rest on IC-MRCI, CASPT2-F12, MRCI-F12, CCSD(T)-F12, and PNO-based local coupled-cluster methods. It is best suited to research groups whose work depends on benchmark-level accuracy (spectroscopy, PES generation, thermochemical kinetics, cold-molecule physics) and who can accommodate its commercial licensing model.

---

## 9. Key Publications — MOLPRO Package and Theory

### 9.1 Primary package-description papers
1. Werner, H.-J.; Knowles, P. J.; Manby, F. R.; Black, J. A.; Doll, K.; Heßelmann, A.; Kats, D.; Köhn, A.; Korona, T.; Kreplin, D. A.; Ma, Q.; Miller, T. F.; Mitrushchenkov, A.; Peterson, K. A.; Polyak, I.; Rauhut, G.; Sibaev, M. **"The Molpro quantum chemistry package."** *J. Chem. Phys.* **152**, 144107 (2020). https://doi.org/10.1063/5.0005081 — *the current primary citation for the package.*
2. Werner, H.-J.; Knowles, P. J.; Knizia, G.; Manby, F. R.; Schütz, M. **"Molpro: a general-purpose quantum chemistry program package."** *WIREs Comput. Mol. Sci.* **2**, 242–253 (2012). https://doi.org/10.1002/wcms.82

### 9.2 Recent methodological papers (explicitly correlated local coupled cluster)
3. Hansen, A.; Knowles, P. J.; Werner, H.-J. **"Accurate Calculation of Noncovalent Interactions Using PNO-LCCSD(T)-F12 in Molpro."** *J. Phys. Chem. A* **129**, 4812–4833 (2025). https://doi.org/10.1021/acs.jpca.5c02316
4. Werner, H.-J.; Hansen, A. **"Local Wave Function Embedding: Correlation Regions in PNO-LCCSD(T)-F12 Calculations."** *J. Phys. Chem. A* **128**, 10936–10947 (2024). https://doi.org/10.1021/acs.jpca.4c06852
5. Werner, H.-J.; Hansen, A. **"Accurate Calculation of Isomerization and Conformational Energies of Larger Molecules Using Explicitly Correlated Local Coupled Cluster Methods in Molpro and ORCA."** *J. Chem. Theory Comput.* **19**, 7007–7030 (2023). https://doi.org/10.1021/acs.jctc.3c00270

### 9.3 Foundational multireference and correlation-method theory
6. Werner, H.-J.; Knowles, P. J. **"A second order multiconfiguration SCF procedure with optimum convergence."** *J. Chem. Phys.* **82**, 5053 (1985).
7. Knowles, P. J.; Werner, H.-J. **"An efficient second-order MC SCF method for long configuration expansions."** *Chem. Phys. Lett.* **115**, 259 (1985).
8. Shamasundar, K. R.; Knizia, G.; Werner, H.-J. **"A new internally contracted multi-reference configuration interaction method."** *J. Chem. Phys.* **135**, 053101 (2011).
9. Knowles, P. J.; Hampel, C.; Werner, H.-J. **"Coupled cluster theory for high spin, open shell reference wave functions."** *J. Chem. Phys.* **99**, 5219–5227 (1993); erratum *J. Chem. Phys.* **112**, 3106–3107 (2000).

### 9.4 Explicit correlation (F12) theory
10. Adler, T. B.; Knizia, G.; Werner, H.-J. **"A simple and efficient CCSD(T)-F12 approximation."** *J. Chem. Phys.* **127**, 221106 (2007).
11. Knizia, G.; Adler, T. B.; Werner, H.-J. **"Simplified CCSD(T)-F12 methods: Theory and benchmarks."** *J. Chem. Phys.* **130**, 054104 (2009).
12. Peterson, K. A.; Adler, T. B.; Werner, H.-J. **"Systematically convergent basis sets for explicitly correlated wavefunctions: The atoms H, He, B–Ne, and Al–Ar."** *J. Chem. Phys.* **128**, 084102 (2008).

### 9.5 Local correlation and NMR/magnetic property methods
13. Loibl, S.; Manby, F. R.; Schütz, M. **"Density fitted, local Hartree–Fock treatment of NMR chemical shifts using London atomic orbitals."** *Mol. Phys.* **108**, 477–485 (2010).
14. Loibl, S.; Schütz, M. **"NMR shielding tensors for density fitted local second-order Møller–Plesset perturbation theory using gauge including atomic orbitals."** *J. Chem. Phys.* **137**, 084107 (2012).
15. Loibl, S.; Schütz, M. **"Magnetizability and rotational g tensors for density fitted local second-order Møller-Plesset perturbation theory using gauge-including atomic orbitals."** *J. Chem. Phys.* **141**, 024108 (2014).

### 9.6 Historical foundational methods (pre-MOLPRO-branding origins)
16. Meyer, W. Original PNO-CI (pseudo-natural orbital configuration interaction) method papers, late 1960s–1970s (*Int. J. Quantum Chem.*, *J. Chem. Phys.*, *Theor. Chim. Acta* series).
17. Meyer, W. PNO-CEPA (pseudo-natural orbital coupled electron pair approximation) method papers, early 1970s.
18. Pulay, P. First analytical Hartree–Fock energy gradient method papers, late 1960s–1970s (*Mol. Phys.*).

---

## 10. Primary Sources Consulted

- Official MOLPRO website: https://www.molpro.net/
- Official product/licensing catalogue: https://www.molpro.net/info/products
- *J. Chem. Phys.* **152**, 144107 (2020), AIP Publishing: https://doi.org/10.1063/5.0005081
- *WIREs Comput. Mol. Sci.* **2**, 242 (2012): https://doi.org/10.1002/wcms.82
- Wikipedia, "MOLPRO": https://en.wikipedia.org/wiki/MOLPRO
- Citation/usage context drawn from multiple recent primary-literature papers citing MOLPRO (2024–2026 arXiv preprints and journal articles) for application-domain and version-history corroboration.

*Note: pricing figures reflect the current published academic non-commercial rate card and may change; commercial-use pricing is quoted directly by the vendor on request.*

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the MOLPRO 	Commercial quantum chemistry package strong in high-accuracy correlated methods, also offering DFT functionality. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
