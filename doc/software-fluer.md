# FLEUR: An Exhaustive Review of the Open-Source FLAPW Code

## 1. Overview

**FLEUR** (Full-potential Linearized augmented plane wave in EURope) is an open-source, all-electron density functional theory (DFT) code that implements the **full-potential linearized augmented plane-wave (FLAPW) method**. It is developed primarily at the **Peter Grünberg Institute / Institute for Advanced Simulation, Forschungszentrum Jülich, Germany**, under the leadership of Stefan Blügel and a large team of contributors (Daniel Wortmann, Gregor Michalicek, Gustav Bihlmayer, and many others).

| Attribute | Detail |
|---|---|
| Method | All-electron FLAPW (DFT) |
| Language | Fortran (with MPI/OpenMP parallelization) |
| License | MIT License (open source) |
| Operating system | Linux (HPC-oriented) |
| Website | [www.flapw.de](https://www.flapw.de) |
| Source repository | GitLab (iffgit.fz-juelich.de/fleur/fleur), mirrored on GitHub (`JuDFTteam/FLEUR`) |
| Latest stable line | MaX-R7.x series (e.g., MaX-R7.1, MaX-R7.2), with a rolling `v26`/master documentation track |
| Primary developer institution | Forschungszentrum Jülich (PGI-1/IAS-1), Germany |

FLAPW is widely regarded as a **"gold-standard" reference method** in DFT because it makes no shape approximation to the potential (true full-potential) and treats *all* electrons — core and valence — explicitly, rather than relying on pseudopotentials. This makes FLEUR one of the benchmark codes used in large-scale DFT-code verification studies (e.g., the 2016 *Science* reproducibility study and its follow-up "How to verify the precision of DFT implementations" work).

---

## 2. The FLAPW Method in FLEUR

### 2.1 Core methodology
FLAPW partitions space into two regions:
- **Muffin-tin (MT) spheres** centered on each atomic nucleus, where wavefunctions are expanded in atomic-like radial functions times spherical harmonics.
- **Interstitial region**, where wavefunctions are expanded in plane waves.

The "linearized" part (LAPW) comes from using energy-independent radial basis functions plus their energy derivatives (Andersen 1975; Koelling & Arbman 1975), avoiding the need to solve secular equations at every energy. The "full-potential" extension (Wimmer, Freeman, Krakauer, Weinert, 1981–1982) removes the muffin-tin/shape approximation on the potential itself, making the method applicable to open, low-symmetry, and complex structures — not just close-packed metals.

Key methodological refinements implemented in FLEUR:
- **Local orbitals (LOs)** for precise treatment of semicore states and elimination of the linearization error for valence and unoccupied states.
- **APW+lo** basis option alongside conventional LAPW+LO.
- Systematic, convergable basis set (unlike pseudopotential plane-wave cutoffs alone, FLAPW convergence can be checked against $l_{max}$, $K_{max}$, and MT radii).

### 2.2 Dimensionality and geometries
FLEUR is unusual among FLAPW codes in natively supporting:
- **3D bulk crystals**
- **2D film geometry** (true thin films/surfaces without artificial vacuum-slab periodic repetition)
- **1D wire/nanotube geometry** (an extension unique to FLEUR, supporting a wide range of chiral symmetries)

This film-geometry capability is a long-standing FLEUR specialty, historically important for surface magnetism and thin-film physics.

---

## 3. Magnetism and Spin-Orbit Coupling: FLEUR's Core Strength

FLEUR's reputation rests heavily on its treatment of complex magnetism, making it a code of choice in the spintronics and magnetism community.

### 3.1 Non-collinear magnetism
FLEUR implements **ab initio non-collinear magnetism** within the FLAPW framework (Kurz, Förster, Nordström, Bihlmayer, Blügel, 2004), allowing local magnetic moments to point in arbitrary directions. This underlies calculations of:
- Spin spirals and generalized Bloch theorem approaches (allowing spin-spiral calculations without large supercells)
- Frustrated and complex magnetic structures
- **Magnetic Skyrmion lattices**

### 3.2 Spin-orbit coupling (SOC)
SOC can be included either perturbatively (second-variation) or fully self-consistently, enabling:
- **Magnetocrystalline anisotropy energy (MAE)** calculations
- **Dzyaloshinskii–Moriya interaction (DMI)** from first principles — a FLEUR specialty (Heide, Bihlmayer, Blügel, 2009; Zimmermann et al., 2014), directly relevant to chiral spin textures and skyrmion stability
- **Rashba and Dresselhaus effects**
- **Topological insulators and Chern insulators**
- Spin-dependent transport: anomalous, spin, and inverse spin Hall effects; spin-orbit torque; anomalous Nernst effect — accessed via linear-response Kubo formalism, often in combination with Wannier90-generated tight-binding models

### 3.3 Exchange parameters and spin excitations
FLEUR can extract:
- Interatomic exchange (Heisenberg) parameters for atomistic spin-dynamics simulations, providing a bridge to multiscale magnetic modeling
- Magnon (spin-wave) dispersions, in part via coupling to the SPEX GW code with ladder-diagram corrections
- Hyperfine fields and field gradients (a natural consequence of the all-electron treatment)

### 3.4 Representative application domains
- Ultrathin magnetic films and surfaces
- Chiral magnetic order and skyrmions at surfaces (e.g., Cr/W(110), Fe/Ir(111)-type systems)
- Graphene and 2D-material spin-orbit physics
- Topological insulators
- Transition-metal, lanthanide, and actinide compounds (elements FLEUR handles natively and accurately due to its all-electron nature)

---

## 4. Beyond Ground-State DFT: Extended Capabilities

| Capability | Notes |
|---|---|
| Exchange-correlation functionals | LDA, a wide range of GGAs, LDA+U |
| Hybrid functionals | PBE0, HSE — implemented directly in the all-electron FLAPW basis |
| Optimized Effective Potential (OEP) / exact exchange | Local exact-exchange potentials within FLAPW |
| GW approximation | Via the companion **SPEX** code, sharing the FLAPW basis for consistent quasiparticle band structures |
| Wannier functions | Maximally localized Wannier functions (MLWFs) constructed within the FLAPW formalism, with a Wannier90 interface |
| DFT+U+V | Intersite Coulomb interaction extension of LDA+U, with U/V parameters obtainable from constrained RPA (cRPA) |
| Constrained RPA (cRPA) | Effective Coulomb interaction (Hubbard U, Hund's J) calculations, often paired with SPEX and Wannier90 |
| Phonons | Density-functional perturbation theory (DFPT) implementation solving the Sternheimer equation within the muffin-tin-centered LAPW basis (a relatively recent, 2023–2024 addition) |
| Green-function embedding / transport | Via the companion **G-Fleur** code for semi-infinite and ballistic transport geometries |
| Forces | All-electron force calculations including core-state and MT-boundary discontinuity contributions, enabling structural relaxation |
| External electric fields | Supported for film/surface geometries |

FLEUR is thus better described as a **code family**: the core FLEUR DFT engine, the **inpgen** input generator, **SPEX** (GW/many-body perturbation theory), and **G-Fleur** (Green-function transport/embedding), all sharing the FLAPW basis-set philosophy.

---

## 5. Software Engineering, Ecosystem, and Distribution

- **License & governance**: FLEUR moved to a fully open-source **MIT license**, hosted on a GitLab instance at Forschungszentrum Jülich (iffgit.fz-juelich.de/fleur/fleur), with a public read-only mirror on GitHub (`JuDFTteam/FLEUR`). Development follows a `develop` → `stable` → `release` branch model.
- **Parallelization**: MPI + OpenMP hybrid parallelization; can leverage optimized linear-algebra and I/O libraries such as **ScaLAPACK, ELPA, Elemental, MAGMA, HDF5, LibXC**, and MPI, targeting HPC clusters.
- **Input/output**: Modern XML-based input/output (`inp.xml`, `out.xml`), which enabled tighter integration with external workflow tools.
- **Ecosystem/interfaces**:
  - **AiiDA-FLEUR** — an AiiDA plugin providing full provenance-tracked workflows (SCF, structure relaxation, DMI, magnetic anisotropy convergence, spin-spiral dispersion, etc.), MIT-licensed, developed by JuDFTteam.
  - **masci-tools** — shared Python utilities for parsing/writing FLEUR I/O, used across the ecosystem.
  - **pymatgen-io-fleur** and **ase-fleur** — I/O plugins integrating FLEUR with the Pymatgen and ASE materials-informatics ecosystems.
  - **GFleur** — Green's-function embedding extension for transport/semi-infinite geometries.
- **European materials infrastructure**: FLEUR is one of the flagship codes of the **MaX (Materials Design at the Exascale)** EU Centre of Excellence, which funds much of its HPC-oriented development (hence the "MaX-Rx.y" release naming).
- **Verification**: FLEUR is one of the reference codes used in cross-code DFT precision/reproducibility benchmarking efforts (e.g., the *Science* 2016 study and its 2024 successor).

---

## 6. Citing FLEUR

The FLEUR team requests that users cite (a) the website, (b) the version-independent Zenodo software entry, and (c) the specific FLEUR version used:

- Software citation (version-independent): D. Wortmann *et al.*, **FLEUR**, Zenodo, DOI: [10.5281/zenodo.7576163](https://doi.org/10.5281/zenodo.7576163)
- Website: [www.flapw.de](https://www.flapw.de)

Beyond the software citation, FLEUR's own documentation explicitly asks users to cite the relevant *methodological* papers behind whichever features they used (see list below).

---

## 7. Publications Related to FLEUR's Theory and Methodology

### 7.1 Foundational LAPW / FLAPW method papers
1. O. K. Andersen, "Linear methods in band theory," *Phys. Rev. B* **12**, 3060 (1975).
2. D. D. Koelling and G. O. Arbman, "Use of energy derivative of the radial solution in an augmented plane wave method: application to copper," *J. Phys. F: Metal Phys.* **5**, 2041 (1975).
3. H. Krakauer, M. Posternak, and A. J. Freeman, "Linearized augmented plane-wave method for the electronic band structure of thin films," *Phys. Rev. B* **19**, 1706 (1979).
4. E. Wimmer, A. J. Freeman, H. Krakauer, and M. Weinert, "Full-potential self-consistent linearized-augmented-plane-wave method for calculating the electronic structure of molecules and surfaces: O₂ molecule," *Phys. Rev. B* **24**, 864 (1981).
5. M. Weinert, E. Wimmer, and A. J. Freeman, "Total-energy all-electron density functional method for bulk solids and surfaces," *Phys. Rev. B* **26**, 4571 (1982).
6. D. J. Singh and L. Nordström, *Planewaves, Pseudopotentials, and the LAPW Method*, Springer (2005). [Comprehensive textbook treatment]
7. S. Blügel and G. Bihlmayer, "Full-Potential Linearized Augmented Planewave Method," in *Computational Nanoscience: Do It Yourself!*, NIC Series Vol. 31, p. 85, John von Neumann Institute for Computing, Jülich (2006).

### 7.2 Basis-set extensions and linearization-error elimination
8. D. Singh, "Ground-state properties of lanthanum: Treatment of extended-core states," *Phys. Rev. B* **43**, 6388 (1991). [Local orbitals]
9. C. Friedrich, A. Schindlmayr, S. Blügel, and T. Kotani, "Elimination of the linearization error in GW calculations based on the linearized augmented-plane-wave method," *Phys. Rev. B* **74**, 045104 (2006).
10. G. Michalicek, M. Betzinger, C. Friedrich, and S. Blügel, "Elimination of the linearization error and improved basis-set convergence within the FLAPW method," *Comp. Phys. Commun.* **184**, 2670 (2013).

### 7.3 Low-dimensional geometries (films, wires)
11. Y. Mokrousov, G. Bihlmayer, and S. Blügel, "A full-potential linearized augmented planewave method for one-dimensional systems: gold nanowire and iron monowires in a gold tube," *Phys. Rev. B* **72**, 045402 (2005).

### 7.4 Magnetism
12. Ph. Kurz, F. Förster, L. Nordström, G. Bihlmayer, and S. Blügel, "Ab initio treatment of non-collinear magnets with the full-potential linearized augmented planewave method," *Phys. Rev. B* **69**, 024415 (2004).
13. M. Heide, G. Bihlmayer, and S. Blügel, "Describing Dzyaloshinskii–Moriya spirals from first principles," *Physica B* **404**, 2678 (2009).
14. B. Zimmermann, M. Heide, G. Bihlmayer, and S. Blügel, "First-principles analysis of a homochiral cycloidal magnetic structure in a monolayer Cr on W(110)," *Phys. Rev. B* **90**, 115427 (2014).
15. E. Şaşıoğlu, A. Schindlmayr, Ch. Friedrich, F. Freimuth, and S. Blügel, "Wannier-function approach to spin excitations in solids," *Phys. Rev. B* **81**, 054434 (2010).

### 7.5 Forces and structural relaxation
16. R. Yu, D. Singh, and H. Krakauer, "All-electron and pseudopotential force calculations using the linearized-augmented-plane-wave method," *Phys. Rev. B* **43**, 6411 (1991).
17. D. A. Klüppelberg, M. Betzinger, and S. Blügel, "Atomic force calculations within the all-electron FLAPW method: Treatment of core states and discontinuities at the muffin-tin sphere boundary," *Phys. Rev. B* **91**, 035105 (2015).

### 7.6 Transport and Green-function embedding
18. D. Wortmann, H. Ishida, and S. Blügel, "Ab initio Green-function formulation of the transfer matrix: Application to complex bandstructures," *Phys. Rev. B* **65**, 165103 (2002).
19. D. Wortmann, H. Ishida, and S. Blügel, "Embedded Green-function approach to the ballistic electron transport through an interface," *Phys. Rev. B* **66**, 075113 (2002).

### 7.7 Hybrid functionals and exact exchange
20. M. Betzinger, C. Friedrich, S. Blügel, and A. Görling, "Local exact exchange potentials within the all-electron FLAPW method and a comparison with pseudopotential results," *Phys. Rev. B* **83**, 045105 (2011).
21. M. Betzinger, C. Friedrich, and S. Blügel, "Hybrid functionals within the all-electron FLAPW method: implementation and applications of PBE0," *Phys. Rev. B* **81**, 195117 (2010).
22. M. Schlipf, M. Betzinger, C. Friedrich, M. Ležaić, and S. Blügel, "HSE hybrid functional within the FLAPW method and its application to GdN," *Phys. Rev. B* **84**, 125142 (2011).
23. M. Betzinger, C. Friedrich, A. Görling, and S. Blügel, "Precise response functions in all-electron methods: Application to the optimized-effective-potential approach," *Phys. Rev. B* **85**, 245124 (2012).

### 7.8 Wannier functions
24. F. Freimuth, Y. Mokrousov, D. Wortmann, S. Heinze, and S. Blügel, "Maximally Localized Wannier Functions within the FLAPW formalism," *Phys. Rev. B* **78**, 035120 (2008).

### 7.9 GW approximation and many-body perturbation theory (SPEX)
25. C. Friedrich, S. Blügel, and A. Schindlmayr, "Efficient calculation of the Coulomb matrix and its expansion around k=0 within the FLAPW method," *Comp. Phys. Comm.* **180**, 347 (2009).
26. C. Friedrich, S. Blügel, and A. Schindlmayr, "Efficient implementation of the GW approximation within the all-electron FLAPW method," *Phys. Rev. B* **81**, 125102 (2010).

### 7.10 Constrained RPA / Hubbard U
27. E. Şaşıoğlu, C. Friedrich, and S. Blügel, "Effective Coulomb interaction in transition metals from constrained random-phase approximation," *Phys. Rev. B* **83**, 121101(R) (2011).

### 7.11 DFT+U+V (recent extension)
28. W. Beida, G. Bihlmayer, C. Friedrich, G. Michalicek, D. Wortmann, and S. Blügel, "Implementation and application of a DFT+U+V approach within the all-electron FLAPW method," *Phys. Rev. B* **114**, 045129 (2026).

### 7.12 Phonons (DFPT within FLAPW)
29. C.-R. Gerhorst, A. Neukirchen, D. A. Klüppelberg, G. Bihlmayer, M. Betzinger, G. Michalicek, D. Wortmann, and S. Blügel, "Phonons from Density-Functional Perturbation Theory using the All-Electron Full-Potential Linearized Augmented Plane-Wave Method FLEUR," *Electron. Struct.* **6**, 017001 (2024).

### 7.13 High-performance computing / numerical implementation
30. A. Lichtenstein *et al.*, "High-performance generation of the Hamiltonian and Overlap matrices in FLAPW methods" (arXiv:1602.06589).
31. "Accelerating the computation of FLAPW methods on heterogeneous architectures" (arXiv:1712.07206).

### 7.14 Cross-code verification / benchmarking studies featuring FLEUR
32. K. Lejaeghere, G. Bihlmayer, T. Björkman, P. Blaha, S. Blügel, V. Blum, D. Caliste, I. E. Castelli, *et al.*, "Reproducibility in density functional theory calculations of solids," *Science* **351**, aad3000 (2016).
33. E. Bosoni, L. Beal, M. Bercx, P. Blaha, S. Blügel, J. Bröder, M. Callsen, S. Cottenier, *et al.*, "How to verify the precision of density-functional-theory implementations via reproducible and universal workflows," *Nat. Rev. Phys.* **6**, 45 (2024).

### 7.15 Primary software citation
34. D. Wortmann, G. Michalicek, N. Baadji, W. Beida, M. Betzinger, G. Bihlmayer, T. Bornhake, J. Bröder, T. Burnus, J. Enkovaara, F. Freimuth, C. Friedrich, C.-R. Gerhorst, S. Granberg Cauchi, U. Grytsiuk, A. Hanke, J.-P. Hanke, M. Heide, S. Heinze, R. Hilgers, H. Janssen, D. A. Klüppelberg, R. Kovacik, P. Kurz, M. Ležaić, G. K. H. Madsen, Y. Mokrousov, A. Neukirchen, M. Redies, S. Rost, M. Schlipf, A. Schindlmayr, M. Winkelmann, and S. Blügel, **FLEUR** [Computer software], Zenodo (2023–), DOI: [10.5281/zenodo.7576163](https://doi.org/10.5281/zenodo.7576163).

---

## 8. Summary Assessment

**Strengths**
- Reference-grade accuracy via the all-electron, full-potential FLAPW method — used as a benchmark by other DFT codes.
- Best-in-class support for non-collinear magnetism, spin-orbit coupling, DMI, and magnetic-anisotropy calculations.
- Native 1D/2D/3D geometry support (unique film/wire treatment without artificial periodic images).
- Broad periodic-table applicability, especially strong for lanthanides, actinides, and transition metals where pseudopotential transferability is often problematic.
- Rich beyond-DFT ecosystem: GW (via SPEX), hybrid functionals, DFT+U(+V), cRPA, Wannier functions, phonons via DFPT, and transport via G-Fleur.
- Fully open source (MIT license) with an active, well-integrated Python/AiiDA/ASE/Pymatgen ecosystem for high-throughput and automated workflows.

**Limitations**
- All-electron FLAPW is computationally more expensive than pseudopotential plane-wave methods, limiting routine system sizes compared to codes like VASP or Quantum ESPRESSO.
- Steeper learning curve; input/setup (muffin-tin radii, basis cutoffs, local-orbital choices) requires more physical/method expertise than typical pseudopotential codes.
- Primarily Linux/HPC-oriented, Fortran-based codebase, with a smaller (though growing) non-specialist user base relative to more turnkey DFT packages.

FLEUR remains one of the principal all-electron FLAPW implementations worldwide (alongside WIEN2k and Elk), and is particularly the code of choice within the community when non-collinear magnetism, spin-orbit-driven phenomena (DMI, topological states, spin-orbit torque), or precise all-electron benchmarking are the primary research objectives.

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the FLEUR 	Open-source full-potential linearized augmented plane-wave (FLAPW) code, strong in magnetism and spin-orbit coupling studies. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
