# VASP (Vienna Ab initio Simulation Package) — An Exhaustive Review

## 1. Overview

The Vienna Ab initio Simulation Package (VASP) is a Fortran-based software package for performing ab initio quantum-mechanical calculations of periodic systems, including electronic structure calculations and ab initio molecular dynamics, from first principles. It is developed and maintained by the Computational Materials Physics group at the University of Vienna (Georg Kresse and collaborators), and since 2018 has been commercially distributed and licensed through VASP Software GmbH.

At its core, VASP solves the many-body Schrödinger problem approximately using a plane-wave basis set combined with either:

- **Ultrasoft pseudopotentials (US-PP)**, following the Vanderbilt scheme, or
- **The projector-augmented-wave (PAW) method**, which is the recommended and most widely used approach today.

VASP is one of the most heavily used codes in condensed-matter physics, computational materials science, surface science, and computational chemistry, and it underlies a very large share of published DFT literature (it is regularly among the most-cited software packages in physics and materials science).

## 2. Core Theoretical Framework

VASP computes the electronic ground state via one of several levels of theory:

- **Density Functional Theory (DFT)** — solving the Kohn–Sham equations self-consistently; this is the default and most-used mode.
- **Hartree–Fock (HF)** — solving the Roothaan-type equations directly.
- **Hybrid functionals** — mixing a fraction of exact (Hartree–Fock) exchange with a GGA/meta-GGA functional (e.g., PBE0, HSE06, B3LYP).
- **Many-body perturbation theory (MBPT)**, including:
  - The **GW approximation** (G₀W₀, GW₀, self-consistent GW, quasiparticle self-consistent variants) for quasiparticle band structures.
  - The **Bethe–Salpeter Equation (BSE)** for optical absorption spectra including excitonic effects.
  - **Random Phase Approximation (RPA)** total energies via the Adiabatic-Connection Fluctuation-Dissipation Theorem (ACFDT).
  - **Second-order Møller–Plesset perturbation theory (MP2)** for periodic systems.
  - Time-dependent DFT / time-dependent Hartree–Fock (Casida equation).
- **DFT+U** (Dudarev and Liechtenstein formulations) for strongly correlated d/f-electron systems.
- **Machine-learned force fields (MLFF)** — an on-the-fly Bayesian active-learning scheme that trains interatomic potentials during ab initio molecular dynamics, dramatically reducing the number of explicit DFT evaluations needed for long/large MD trajectories.

Electronic wavefunctions and the local potential are expanded in a **plane-wave basis**, with core-valence interactions represented via PAW or US-PP. The self-consistency cycle is handled by efficient **iterative matrix diagonalization** (blocked Davidson, or the residual minimization method with direct inversion of the iterative subspace, RMM-DIIS), combined with density mixing schemes (typically Broyden/Pulay mixing).

## 3. Core Capabilities

| Category | Capabilities |
|---|---|
| **Ground-state DFT** | Total energy, forces, stress tensor; structural relaxation; equations of state |
| **Exchange-correlation** | LDA, GGA (PBE, PW91, RPBE, etc.), meta-GGA (SCAN, r2SCAN, TPSS), hybrids (HSE06, PBE0, B3LYP) |
| **van der Waals corrections** | DFT-D2/D3/D4, Tkatchenko–Scheffler, many-body dispersion (MBD), vdW-DF/optPBE/optB88/rev-vdW-DF2 non-local functionals |
| **Magnetism** | Collinear and non-collinear spin-polarized calculations, spin-orbit coupling, DFT+U |
| **Excited states / spectroscopy** | GW quasiparticle energies, BSE optical spectra, RPA correlation energies, linear-response dielectric functions, IPES/XPS-related core-level shifts |
| **Molecular dynamics** | Born-Oppenheimer AIMD (NVE, NVT via Nosé–Hoover/Langevin thermostats, NPT via Parrinello–Rahman), path-integral MD |
| **Machine learning** | On-the-fly Bayesian MLFF generation and prediction-only MD with trained force fields |
| **Vibrational/phonon properties** | Finite-difference and DFPT (density-functional perturbation theory) phonons, IR/Raman intensities |
| **Electronic structure analysis** | Band structures, densities of states (DOS, PDOS), Bader/Wannier-compatible outputs, Berry-phase polarization, Wannier90 interface |
| **Transport & response** | Elastic tensors, piezoelectric tensors, dielectric/optical tensors, effective masses |
| **NEB / transition-state search** | Nudged elastic band and variants for reaction pathways and diffusion barriers |
| **Parallelization** | MPI (k-point, band, plane-wave, orbital), OpenMP hybrid, GPU acceleration via OpenACC (NVIDIA GPUs) |

## 4. Numerical/Algorithmic Backbone

- **Iterative diagonalization:** blocked-Davidson and RMM-DIIS algorithms avoid full-matrix diagonalization, which is essential given the huge basis sets involved.
- **Charge/density mixing:** Broyden/Pulay-type mixers accelerate SCF convergence.
- **k-point sampling:** Monkhorst–Pack and Γ-centered grids, symmetry reduction.
- **PAW method:** Reconstructs the full all-electron wavefunction from pseudo-wavefunctions plus atom-centered augmentation, giving near all-electron (FLAPW-comparable) accuracy at pseudopotential-like cost.
- **Cubic-scaling RPA/GW algorithms:** later developments (Kaltak, Klimeš, Kresse) reduced the formal scaling of RPA/GW calculations from O(N⁴) to O(N³), substantially extending accessible system sizes.

## 5. Licensing and Distribution

VASP is **proprietary, closed-source software**; it is not free or open-source. A paid academic or commercial license (through VASP Software GmbH, or historically the University of Vienna) is required, and source code access is restricted to license holders. This is a notable contrast with several competing DFT codes (Quantum ESPRESSO, ABINIT, CP2K, GPAW) that are open source.

## 6. Version History Highlights

- **VASP 4.x** — early widely-used versions based on US-PP and PAW; established VASP as a mainstream DFT workhorse.
- **VASP 5.x** — expanded exchange-correlation functional library, hybrid functionals, GW/BSE modules, van der Waals functionals.
- **VASP 6.x** (2019–present) — major rewrite adding: an integrated one-step GW workflow, improved MP2/RPA implementations, machine-learned force fields, GPU (OpenACC) acceleration, improved non-collinear magnetism and spin-orbit handling, and continued performance/scalability improvements for exascale-class HPC systems.
- **VASP 6.5.x** (2024–2025) — further MLFF refinements, expanded meta-GGA and dielectric-response capability.
- **VASP 6.6.x** (2025–2026) — most recent stable line; 6.6.1 is a patch release specifically correcting a force-calculation bug introduced for collinear spin-polarized meta-GGA calculations in 6.6.0.

## 7. Strengths

- **Accuracy and maturity:** PAW implementation gives near all-electron accuracy; extremely well-validated across decades of literature.
- **Breadth of methods:** Few codes offer as unified a stack spanning DFT → hybrids → GW/BSE → RPA → MP2 → MLFF-accelerated MD within a single package.
- **Performance and scalability:** Multi-level parallelization (MPI/OpenMP/GPU) makes it competitive on modern HPC clusters; long-standing, heavily optimized Fortran codebase.
- **Ecosystem:** Extensive pseudopotential (POTCAR) library, mature documentation/wiki, large user community, interfaces to companion tools (Wannier90, Phonopy, ASE, pymatgen, VASPKIT, VTST tools for NEB/transition states).
- **Machine-learning acceleration:** The built-in on-the-fly Bayesian MLFF module is one of the more mature "self-training" force-field implementations bundled directly into a mainstream DFT code, letting users get large-scale/long-time MD without needing an external ML pipeline.

## 8. Limitations / Criticisms

- **Cost and access:** Commercial licensing is a real barrier for smaller groups, some countries, and reproducibility/open-science efforts, especially compared to open-source alternatives.
- **Closed source:** Source-level scrutiny and community-driven core development are restricted to license holders, limiting independent code auditing.
- **Plane-wave/periodic-only formalism:** Like other plane-wave PAW codes, VASP is fundamentally built around periodic boundary conditions; isolated molecules/clusters require supercell approaches, which can be less natural than in Gaussian-basis quantum chemistry codes.
- **Steep learning curve:** Correct/robust use (choice of pseudopotentials, k-point/energy cutoff convergence, functional choice, GW/BSE workflow parameters) requires substantial domain expertise; incorrect settings are a common source of unreliable published results in the wider literature.
- **GW/BSE/RPA cost:** Despite cubic-scaling advances, many-body perturbation theory calculations remain computationally expensive relative to standard DFT, limiting system sizes.
- **GPU support scope:** GPU acceleration (OpenACC, NVIDIA-focused) does not yet cover every feature to the same degree as the CPU code path, and build/portability across different GPU vendors/HPC software stacks can be nontrivial.

## 9. Typical Use Cases

- Bulk and surface electronic structure (band structures, DOS) of metals, semiconductors, insulators
- Defect and dopant energetics in semiconductors
- Catalysis and surface reaction studies (adsorption energies, reaction pathways via NEB)
- Battery and energy-materials modeling (ion diffusion, voltage profiles, interfaces)
- Magnetic materials and spintronics (non-collinear magnetism, spin-orbit coupling, DFT+U for correlated oxides)
- Optical and excitonic properties of 2D materials and semiconductors (GW+BSE)
- High-temperature/high-pressure phase behavior via AIMD and MLFF-accelerated MD
- Phonon and thermal transport properties

## 10. Summary Assessment

VASP remains one of the most capable, best-validated, and most widely adopted plane-wave DFT codes in computational materials science, distinguished by its combination of numerical robustness, methodological breadth (DFT through GW/BSE/RPA/MP2 and now MLFF-driven MD), and strong HPC performance. Its principal drawbacks are non-technical: proprietary licensing and closed-source distribution constrain access and independent verification relative to open-source competitors, and its many advanced features demand real expertise to use correctly. For groups with access to a license and the requisite expertise, it remains a de facto standard reference code against which many other electronic-structure methods and packages are benchmarked.

---

# Publications Related to VASP's Underlying Theory

These are the foundational and key extension papers describing the theoretical methods implemented in VASP. (Citation style: author(s), journal, volume, page/article number, year.)

### Foundational VASP method papers
- G. Kresse and J. Hafner, *Ab initio molecular dynamics for liquid metals*, Phys. Rev. B **47**, 558 (1993)
- G. Kresse and J. Hafner, *Ab initio molecular-dynamics simulation of the liquid-metal–amorphous-semiconductor transition in germanium*, Phys. Rev. B **49**, 14251 (1994)
- G. Kresse and J. Furthmüller, *Efficiency of ab-initio total energy calculations for metals and semiconductors using a plane-wave basis set*, Comput. Mater. Sci. **6**, 15 (1996)
- G. Kresse and J. Furthmüller, *Efficient iterative schemes for ab initio total-energy calculations using a plane-wave basis set*, Phys. Rev. B **54**, 11169 (1996)
- J. Hafner, *Ab-initio simulations of materials using VASP: Density-functional theory and beyond*, J. Comput. Chem. **29**, 2044 (2008)

### PAW method and pseudopotentials
- P. E. Blöchl, *Projector augmented-wave method*, Phys. Rev. B **50**, 17953 (1994)
- G. Kresse and D. Joubert, *From ultrasoft pseudopotentials to the projector augmented-wave method*, Phys. Rev. B **59**, 1758 (1999)
- D. Vanderbilt, *Soft self-consistent pseudopotentials in a generalized eigenvalue formalism*, Phys. Rev. B **41**, 7892 (1990)

### Exchange–correlation functionals used within VASP
- J. P. Perdew, K. Burke, and M. Ernzerhof, *Generalized Gradient Approximation Made Simple*, Phys. Rev. Lett. **77**, 3865 (1996)
- J. P. Perdew and A. Zunger, *Self-interaction correction to density-functional approximations for many-electron systems*, Phys. Rev. B **23**, 5048 (1981)
- J. Heyd, G. E. Scuseria, and M. Ernzerhof, *Hybrid functionals based on a screened Coulomb potential*, J. Chem. Phys. **118**, 8207 (2003) (HSE hybrid functional)
- A. V. Krukau, O. A. Vydrov, A. F. Izmaylov, and G. E. Scuseria, *Influence of the exchange screening parameter on the performance of screened hybrid functionals*, J. Chem. Phys. **125**, 224106 (2006)
- J. Sun, A. Ruzsinszky, and J. P. Perdew, *Strongly Constrained and Appropriately Normed Semilocal Density Functional (SCAN)*, Phys. Rev. Lett. **115**, 036402 (2015)

### DFT+U for correlated electron systems
- S. L. Dudarev, G. A. Botton, S. Y. Savrasov, C. J. Humphreys, and A. P. Sutton, *Electron-energy-loss spectra and the structural stability of nickel oxide: An LSDA+U study*, Phys. Rev. B **57**, 1505 (1998)
- A. I. Liechtenstein, V. I. Anisimov, and J. Zaanen, *Density-functional theory and strong interactions: Orbital ordering in Mott-Hubbard insulators*, Phys. Rev. B **52**, R5467 (1995)

### Many-body perturbation theory: GW, BSE, RPA/ACFDT, MP2
- L. Hedin, *New Method for Calculating the One-Particle Green's Function with Application to the Electron-Gas Problem*, Phys. Rev. **139**, A796 (1965) (origin of the GW approximation)
- M. Shishkin and G. Kresse, *Implementation and performance of the frequency-dependent GW method within the PAW framework*, Phys. Rev. B **74**, 035101 (2006)
- M. Shishkin and G. Kresse, *Self-consistent GW calculations for semiconductors and insulators*, Phys. Rev. B **75**, 235102 (2007)
- F. Fuchs, J. Furthmüller, F. Bechstedt, M. Shishkin, and G. Kresse, *Quasiparticle band structure based on a generalized Kohn-Sham scheme*, Phys. Rev. B **76**, 115109 (2007)
- T. Sander, E. Maggio, and G. Kresse, *Beyond the Tamm-Dancoff approximation for extended systems using exact diagonalization*, Phys. Rev. B **92**, 045209 (2015) (BSE implementation)
- S. Albrecht, L. Reining, R. Del Sole, and G. Onida, *Ab Initio Calculation of Excitonic Effects in the Optical Spectra of Semiconductors*, Phys. Rev. Lett. **80**, 4510 (1998)
- M. Rohlfing and S. G. Louie, *Electron-Hole Excitations in Semiconductors and Insulators*, Phys. Rev. Lett. **81**, 2312 (1998)
- J. Paier, M. Marsman, and G. Kresse, *Dielectric properties and excitons for extended systems from hybrid functionals*, Phys. Rev. B **78**, 121201 (2008)
- M. Kaltak, J. Klimeš, and G. Kresse, *Cubic scaling algorithm for the random phase approximation: Implementation and application to water*, Phys. Rev. B **90**, 054115 (2014)
- J. Klimeš, M. Kaltak, and G. Kresse, *Predictive GW calculations using plane waves and pseudopotentials*, J. Chem. Phys. **143**, 102816 (2015)
- P. Liu, M. Kaltak, J. Klimeš, and G. Kresse, *Cubic scaling GW: Towards fast quasiparticle calculations*, Phys. Rev. B **94**, 165109 (2016)
- B. Ramberger, T. Schäfer, and G. Kresse, *Analytic Interatomic Forces in the Random Phase Approximation*, Phys. Rev. Lett. **118**, 106403 (2017)
- M. Marsman, A. Grüneis, J. Paier, and G. Kresse, *Second-order Møller–Plesset perturbation theory applied to extended systems*, J. Chem. Phys. **130**, 184103 (2009)
- H. N. Rojas, R. W. Godby, and R. J. Needs, *Space-Time Method for Ab Initio Calculations of Self-Energies and Dielectric Response Functions of Solids*, Phys. Rev. Lett. **74**, 1827 (1995)

### van der Waals corrections
- S. Grimme, J. Antony, S. Ehrlich, and H. Krieg, *A consistent and accurate ab initio parametrization of density functional dispersion correction (DFT-D) for the 94 elements H-Pu*, J. Chem. Phys. **132**, 154104 (2010)
- A. Tkatchenko and M. Scheffler, *Accurate Molecular Van Der Waals Interactions from Ground-State Electron Density and Free-Atom Reference Data*, Phys. Rev. Lett. **102**, 073005 (2009)

### Machine-learned force fields (on-the-fly, Bayesian)
- R. Jinnouchi, F. Karsai, and G. Kresse, *On-the-fly machine learning force field generation: Application to melting points*, Phys. Rev. B **100**, 014105 (2019)
- R. Jinnouchi, J. Lahnsteiner, F. Karsai, G. Kresse, and M. Bokdam, *Phase Transitions of Hybrid Perovskites Simulated by Machine-Learning Force Fields Trained on the Fly with Bayesian Inference*, Phys. Rev. Lett. **122**, 225701 (2019)
- R. Jinnouchi, F. Karsai, C. Verdi, R. Asahi, and G. Kresse, *Descriptors representing two- and three-body atomic distributions and their effects on the accuracy of machine-learned interatomic potentials*, J. Chem. Phys. **152**, 234102 (2020)
- R. Jinnouchi, K. Miwa, F. Karsai, G. Kresse, and R. Asahi, *On-the-Fly Active Learning of Interatomic Potentials for Large-Scale Atomistic Simulations*, J. Phys. Chem. Lett. **11**, 6946 (2020)
- R. Jinnouchi, F. Karsai, and G. Kresse, *Making free-energy calculations routine: Combining first principles with machine learning*, Phys. Rev. B **101**, 060201 (2020)

### k-point sampling and smearing methods
- H. J. Monkhorst and J. D. Pack, *Special points for Brillouin-zone integrations*, Phys. Rev. B **13**, 5188 (1976)
- M. Methfessel and A. T. Paxton, *High-precision sampling for Brillouin-zone integration in metals*, Phys. Rev. B **40**, 3616 (1989)
- P. E. Blöchl, O. Jepsen, and O. K. Andersen, *Improved tetrahedron method for Brillouin-zone integrations*, Phys. Rev. B **49**, 16223 (1994)

### Foundational density-functional theory
- P. Hohenberg and W. Kohn, *Inhomogeneous Electron Gas*, Phys. Rev. **136**, B864 (1964)
- W. Kohn and L. J. Sham, *Self-Consistent Equations Including Exchange and Correlation Effects*, Phys. Rev. **140**, A1133 (1965)

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the Vienna Ab Initio Simulation Package (VASP). Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
