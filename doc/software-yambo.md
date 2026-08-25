# Yambo: A Review of the Many-Body Perturbation Theory Code for Excited-State Properties

## 1. Overview

**Yambo** is an open-source, plane-wave, pseudopotential-based *ab initio* code for computing excited-state electronic and optical properties of condensed matter systems (solids, surfaces, nanostructures, and molecules) from first principles. It implements the **Green's function formulation of many-body perturbation theory (MBPT)** — most notably the **GW approximation** for quasiparticle (charged) excitations and the **Bethe–Salpeter equation (BSE)** for neutral (optical) excitations — alongside **time-dependent density functional theory (TDDFT)** and **non-equilibrium Green's function (NEGF)** real-time methods.

Yambo does not perform ground-state DFT calculations itself. Instead, it is designed to **post-process the Kohn–Sham (KS) electronic structure produced by a DFT code**, most commonly **Quantum ESPRESSO**, but also **Abinit**, using these as a starting point (wavefunctions, KS eigenvalues, pseudopotentials) for the many-body corrections. This two-step workflow (DFT ground state → MBPT excited state) is standard in the field and is shared by comparable codes such as BerkeleyGW, VASP's GW/BSE module, Abinit's own MBPT routines, West, Fiesta, MolGW, and Exciting.

The project originated from a legacy code named "SELF," was released under the GNU GPL, and is now developed and maintained collaboratively by roughly 20 active contributors, hosted publicly on GitHub, and supported in part by the European MaX Centre of Excellence for HPC applications in materials science.

---

## 2. Theoretical Scope and Capabilities

### 2.1 Electronic (quasiparticle) properties — the GW approximation
- Quasiparticle energy corrections to DFT/KS eigenvalues via the **GW self-energy**, solving the quasiparticle equation derived from Hedin's equations.
- Frequency treatment of the screened Coulomb interaction W via:
  - **Plasmon-pole models** (Godby–Needs and related approximations) for computational efficiency.
  - **Full-frequency GW**, integrating the self-energy over real or imaginary frequency axes.
  - A **multipole approximation** to the frequency dependence of the dielectric screening, offering a middle ground between plasmon-pole speed and full-frequency accuracy (notably efficient for metals and 2D materials).
- **COHSEX approximation** as a static, simplified limit of GW.
- Quasiparticle **lifetimes and renormalization (Z) factors**, not just energy shifts.
- Support for calculations on **metals**, including specialized full-frequency and stochastic-integration screening approaches.
- **Efficient/accelerated GW** algorithms: interpolation of the screened interaction in momentum and frequency space, stochastic integration of the screened potential, and extrapolar-energy correction schemes to accelerate convergence with respect to empty states.
- **GW from occupied states only** approaches, reducing the need for expensive unoccupied-state summations.

### 2.2 Optical and neutral excitations
- **Bethe–Salpeter equation (BSE)** solved on top of GW (or scissor-corrected DFT) quasiparticle energies, capturing electron–hole (excitonic) interactions.
- Support for the **Tamm–Dancoff approximation** and options to go beyond it (capturing exciton–plasmon mixing/coupling effects).
- **RPA (random phase approximation)** independent-particle response as a baseline.
- **TDDFT** within the **adiabatic LDA (TD-LDA/ALDA)** and **long-range-corrected (LRC)** kernels as computationally cheaper alternatives to full BSE.
- Numerical efficiency improvements enabling treatment of systems with large numbers of electron–hole pairs (beyond 10⁵).
- **Electron energy loss spectroscopy (EELS)** and dynamical polarizability calculations.
- Post-processing tools for **exciton wavefunction analysis and localization**, in both real and reciprocal space, including visualization of exciton character.
- **Magneto-optical properties**, e.g., the magneto-optical Kerr effect in magnetic materials.
- **Surface spectroscopy** capabilities (e.g., reflectance anisotropy/differential reflectivity for surfaces and 2D systems).

### 2.3 Electron–phonon coupling and temperature effects
- Static and dynamic (frequency-dependent) many-body perturbation theory treatment of the **electron–phonon interaction**.
- Calculation of **temperature-dependent electronic structure and optical spectra**, including zero-point renormalization of band gaps.
- Excitation–phonon interaction, relevant to phonon-assisted optical processes and exciton–phonon coupling.

### 2.4 Real-time and non-equilibrium simulations (Yambo-RT / NEGF)
- **Real-time propagation** of the reduced single-particle density matrix (Bloch states) to simulate:
  - **Linear and non-linear optical response** (e.g., harmonic generation).
  - **Pump–probe experiments** with finite laser pulses, capturing coherent and incoherent carrier dynamics.
  - **Non-equilibrium Green's function (NEGF)** dynamics for dissipative processes and carrier/exciton relaxation.
- Simulation support relevant to **time-resolved ARPES** spectroscopy and non-equilibrium quasiparticle/excitonic energies.
- Various time propagators/integrators (Euler, RK2, RK4, Heun, exponential/implicit schemes) selectable for stability/accuracy trade-offs.

### 2.5 Additional / auxiliary features
- **Exact Coulomb truncation** techniques for accurate treatment of reduced-dimensionality systems (slabs, wires, isolated molecules) in a periodic supercell, avoiding spurious image interactions.
- Support for **norm-conserving pseudopotentials**.
- A range of exchange-correlation functionals inherited from the DFT starting point (LDA, PBE, PW91, and other GGAs).
- Advanced post-processing and data-analysis tooling for inspecting the simulation flow and intermediate databases.

---

## 3. Software Architecture and Design

- **Languages**: primarily **Fortran**, with components in **C** for I/O and system-level operations.
- **Interfacing with DFT codes**: Yambo requires a converted, pre-computed ground-state electronic structure. It ships with interface utilities (`p2y` for Quantum ESPRESSO, `a2y` for Abinit, and related tools) that convert the DFT output into Yambo's internal database format.
- **Data management**: Yambo organizes intermediate and output data into a set of **NetCDF/HDF5-format databases** (often numbering in the tens of files, with the largest reaching several gigabytes for large systems), chosen to optimize I/O performance and data portability across HPC platforms.
- **User interface**: a friendly, self-documenting **command-line and runlevel-based interface**, where a `yambo -h`-style flag system generates human-readable input files listing all relevant variables for a given run type. Yambo reports estimated elapsed and remaining time during execution.
- **Parallelization**: hybrid **MPI + OpenMP** parallelism across multiple levels (k-points, bands, frequencies, matrix elements, etc.), designed for **massively parallel, distributed-memory HPC architectures**.
- **GPU acceleration**: support for **NVIDIA GPU-accelerated architectures** (e.g., DGX systems), demonstrated to scale efficiently to over 1000 GPUs in benchmark studies (e.g., ~1440 GPUs / 360 nodes on JUWELS-Booster).
- **Portability**: runs on large HPC systems (HPE Cray EX, BullSequana XH, Fujitsu, IBM Power, NVIDIA DGX, Lenovo ThinkSystem), workstations (x86, ARM), and personal computers under Linux, macOS, and Windows (via Docker container).
- **Automation / workflow integration**: Python-based workflow tools (e.g., `YamboCalculation`, `YamboConvergence` within the AiiDA plugin ecosystem) automate the complex, multi-parameter convergence procedures characteristic of GW-BSE calculations, supporting reproducibility and high-throughput studies.
- **Licensing**: released under the **GNU General Public License (GPL)**, open source.

---

## 4. Typical Simulation Workflow

1. **Ground-state DFT calculation** in Quantum ESPRESSO or Abinit to obtain converged Kohn–Sham wavefunctions and eigenvalues (including empty/unoccupied states, which must be well converged for GW).
2. **Interface/conversion step** (`p2y`/`a2y`) to translate the DFT output into Yambo's internal database format.
3. **Setup run** in Yambo to initialize databases and report system information (k-points, symmetries, bands, etc.).
4. **Response function calculation**: independent-particle or RPA dielectric response, needed as the screening input to GW.
5. **GW quasiparticle calculation**: self-energy evaluation (plasmon-pole, multipole, or full-frequency) to obtain corrected quasiparticle energies, typically reported as a "scissor-like" or state-dependent correction to the DFT band structure.
6. **Optional BSE/TDDFT step**: construction and diagonalization (direct or iterative, e.g., Haydock/Lanczos) of the excitonic Hamiltonian using the GW-corrected energies, yielding optical absorption spectra, exciton binding energies, and exciton wavefunctions.
7. **Optional real-time/NEGF step**: time propagation for nonlinear optics or pump–probe simulations, if non-equilibrium properties are of interest.
8. **Post-processing**: extraction and plotting of spectra, band structures, exciton character, and other observables using Yambo's auxiliary tools.

Because of the many convergence parameters involved (empty bands, dielectric matrix cutoff, k-point sampling, frequency integration parameters), convergence studies are a substantial and well-recognized part of any Yambo-based project, and the code and its associated Python tooling include specific features to streamline this process.

---

## 5. Typical Applications

- Accurate **quasiparticle band gaps** and band alignments in semiconductors, insulators, and 2D materials.
- **Defect quasiparticle levels** in wide-gap materials (relevant to qubit and defect-photophysics studies).
- **Optical absorption spectra** capturing excitonic effects, including in low-dimensional and van der Waals materials.
- **Exciton, plasmon, and magnon physics**, including localization and spatial character of excitons.
- **ARPES-relevant quasiparticle spectral functions**, including time-resolved (pump–probe) ARPES simulations.
- **Temperature-dependent electronic and optical properties** via electron–phonon coupling (e.g., zero-point renormalization, temperature-dependent band gaps).
- **Non-linear optical properties** (e.g., harmonic generation) via real-time propagation.
- Benchmark and reproducibility studies of GW methodology across codes (e.g., cross-code G₀W₀ comparison efforts).

---

## 6. Strengths and Limitations (as generally characterized in the literature and documentation)

**Strengths**
- Comprehensive, actively maintained implementation spanning GW, BSE, TDDFT, electron–phonon coupling, and real-time/NEGF methods within a single package.
- Strong HPC scalability, including GPU support, aimed at large-scale, high-throughput MBPT calculations.
- Flexible interfacing with widely used, freely available ground-state codes (Quantum ESPRESSO, Abinit).
- Rich set of numerical approximations (plasmon-pole, multipole, full-frequency, stochastic) allowing users to trade off accuracy against computational cost.
- Active community, dedicated schools/workshops, and integration with high-throughput/automation frameworks (e.g., AiiDA).

**Limitations / considerations**
- As with all plane-wave pseudopotential GW-BSE codes, calculations are computationally demanding and convergence with respect to empty bands, dielectric cutoff, and k-point sampling requires careful, often costly, testing.
- Being a post-processing MBPT code, results are contingent on the quality and convergence of the underlying DFT ground state supplied by the interfaced code.
- Norm-conserving pseudopotential support is the standard route (though this has been extended over the project's history); ultrasoft/PAW compatibility considerations depend on the interface version and DFT code used.

---

## 7. Key Publications on Yambo's Theory and Methodology

### Foundational code papers
- **A. Marini, C. Hogan, M. Grüning, D. Varsano**, *"yambo: An ab initio tool for excited state calculations,"* Computer Physics Communications **180**, 1392–1403 (2009). https://doi.org/10.1016/j.cpc.2009.02.003 — the original code paper describing the GW and BSE implementation.
- **D. Sangalli, A. Ferretti, H. Miranda, C. Attaccalite, I. Marri, E. Cannuccia, P. Melo, M. Marsili, F. Paleari, A. Marrazzo, G. Prandini, P. Bonfà, M. O. Atambo, F. Affinito, M. Palummo, A. Molina-Sánchez, C. Hogan, M. Grüning, D. Varsano, A. Marini**, *"Many-body perturbation theory calculations using the yambo code,"* Journal of Physics: Condensed Matter **31**, 325902 (2019). https://doi.org/10.1088/1361-648X/ab15d0 (arXiv:1902.03837) — the comprehensive update describing electron–phonon coupling, real-time propagation, BSE efficiency improvements, exciton analysis tools, and HPC parallelization strategy.

### Methodological/theoretical underpinnings
- **L. Hedin**, *"On correlation effects in electron spectroscopies and the GW approximation,"* Journal of Physics: Condensed Matter **11**, R489 (1999) — foundational GW theory.
- **G. Onida, L. Reining, A. Rubio**, *"Electronic excitations: density-functional versus many-body Green's-function approaches,"* Reviews of Modern Physics **74**, 601 (2002) — key review of the GW-BSE formalism underlying Yambo's approach.
- **S. L. Adler**, *"Quantum theory of the dielectric constant in real solids,"* Physical Review **126**, 413 (1962).
- **N. Wiser**, *"Dielectric constant with local field effects included,"* Physical Review **129**, 62 (1963).

### Coulomb truncation and reduced-dimensionality systems
- **C. A. Rozzi, D. Varsano, A. Marini, E. K. U. Gross, A. Rubio**, *"Exact Coulomb cutoff technique for supercell calculations,"* Physical Review B **73**, 205119 (2006).

### Beyond Tamm–Dancoff / exciton–plasmon coupling
- **M. Grüning, A. Marini, X. Gonze**, *"Exciton-Plasmon States in Nanoscale Materials: Breakdown of the Tamm–Dancoff Approximation,"* Nano Letters **9**, 2820–2824 (2009).

### GW methodology and efficiency developments
- **P. Umari, G. Stenuit, S. Baroni**, *"GW quasiparticle spectra from occupied states only,"* Physical Review B **81**, 115104 (2010).
- **D. Alfè, D. A. Leon (et al., various author orderings across related works), A. Ferretti, D. Varsano, E. Molinari, C. Cardoso**, *"Efficient full frequency GW for metals using a multipole approach for the dielectric screening,"* Physical Review B **107**, 155130 (2023).
- **A. Guandalini, D. A. Leon, P. D'Amico, C. Cardoso, A. Ferretti, M. Rontani, D. Varsano**, *"Efficient GW calculations via interpolation of the screened interaction in momentum and frequency space: The case of graphene,"* Physical Review B **109**, 075120 (2024).
- **A. Guandalini, P. D'Amico, A. Ferretti, D. Varsano**, *"Efficient GW calculations in two dimensional materials through a stochastic integration of the screened potential,"* npj Computational Materials **9**, 44 (2023).
- Related work on frequency dependence via multipole approximation ("Frequency dependence in GW made simple using a multipole approximation").

### High-throughput and reproducibility
- **M. Bonacci, J. Qiao, N. Spallanzani, A. Marrazzo, G. Pizzi, E. Molinari, D. Varsano, A. Ferretti, D. Prezzi**, *"Towards high-throughput many-body perturbation theory: efficient algorithms and automated workflows,"* npj Computational Materials **9**, 74 (2023).
- **T. Rangel, M. Del Ben, D. Varsano, G. Antonius, F. Bruneval, F. H. da Jornada, M. J. van Setten, O. K. Orhan, D. D. O'Regan, A. Canning, A. Ferretti, A. Marini, G.-M. Rignanese, J. Deslippe, S. G. Louie, et al.**, *"Reproducibility in G₀W₀ Calculations for Solids,"* Computer Physics Communications **255**, 107242 (2020) — cross-code benchmarking effort including Yambo.

### Related/comparable code papers cited for context in Yambo's own documentation
- **C. Faber, I. Duchemin, T. Deutsch, C. Attaccalite, V. Olevano, X. Blase**, *"Electron–phonon coupling and charge-transfer excitations in organic systems from many-body perturbation theory: The Fiesta code,"* Journal of Materials Science (2013), describing the Fiesta GW/BSE Gaussian-basis code.
- **F. Bruneval, T. Rangel, S. M. Hamed, M. Shao, C. Yang, J. B. Neaton**, *"molgw 1: Many-body perturbation theory software for atoms, molecules, and clusters,"* Computer Physics Communications **208**, 149–161 (2016).

*Note: For the complete and continuously updated bibliography, consult the official Yambo publications page (yambo-code.eu / yambo-code.org), which lists papers describing each major methodological extension (e.g., electron–phonon, NEGF real-time, exciton analysis tools, GPU implementation) as they are published.*

---

## 8. Access and Resources

- **Official website**: yambo-code.eu / yambo-code.org
- **Source code**: hosted publicly on GitHub
- **License**: GNU GPL (open source)
- **Interfaced ground-state codes**: Quantum ESPRESSO, Abinit
- **Community resources**: an active user wiki/forum, periodic "Yambo Schools" (hands-on training events), and support through the MaX Centre of Excellence for HPC in materials science.



---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the Yambo 	Many-body perturbation theory code (GW, BSE) built to work alongside ground-state DFT codes (e.g., Quantum ESPRESSO, Abinit) for excited states. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
