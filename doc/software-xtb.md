# xtb (Extended Tight-Binding) — An Exhaustive Review

**Developer:** Grimme group, Mulliken Center for Theoretical Chemistry, University of Bonn
**Repository:** https://github.com/grimme-lab/xtb
**License:** GNU Lesser General Public License v3 (or later)
**Language:** Fortran (with a C API and Python/other-language bindings)
**Documentation:** https://xtb-docs.readthedocs.io

---

## 1. Overview

`xtb` is an open-source program package implementing the **semi-empirical extended tight-binding (xTB)** family of quantum-chemical methods developed by Stefan Grimme's group. It is designed to occupy the middle ground between classical force fields and Kohn–Sham DFT: it retains an explicit, self-consistent electronic-structure treatment (so it produces orbitals, charges, dipole moments, and a proper wavefunction/density) while running orders of magnitude faster than DFT, enabling routine treatment of systems from small molecules to large biomolecular or materials systems (up to several thousand atoms, essentially the whole periodic table up to Z = 86, radon).

The methods implemented are collectively branded **GFN**, an acronym originally standing for **G**eometries, **F**requencies, and **N**on-covalent interactions — the three classes of properties the methods were first parametrized to reproduce reliably. Over time the GFN family has broadened well beyond these original targets to cover thermochemistry, spectroscopy, reactivity, and dynamics.

### 1.1 Design philosophy

- **Minimal, mostly element-specific parametrization** (as opposed to atom-pair-specific parametrization used in classical DFTB), which is what allows near-complete periodic-table coverage without an explosion of fitted parameters.
- **Physically motivated energy decomposition** built from a self-consistent-charge tight-binding (SCC-TB) electronic term plus explicit classical/semiclassical correction terms for repulsion, dispersion, electrostatics, and (in some variants) halogen bonding.
- **"Black-box" robustness** — intended to run without manual tuning across essentially arbitrary input geometries and elements, including transition metals, radicals, and ions, which is unusual for semi-empirical methods.
- **Speed as a first-class design goal** — used for large-scale conformer/rotamer sampling, molecular dynamics, and high-throughput screening where DFT is computationally prohibitive.

---

## 2. The GFN Method Family

`xtb` is best understood as a *suite* of increasingly refined Hamiltonians rather than a single method. Each generation trades off cost, generality, and accuracy differently.

### 2.1 GFN0-xTB
- A **non-iterative** (non-self-consistent) zeroth-order tight-binding method.
- Total energy: repulsion + D4-type dispersion + short-range bond (SRB) correction + electronegativity-equilibration (EEQ) charge model + extended-Hückel-type (EHT) band-structure energy.
- Because there is no SCF cycle, it avoids the convergence failures that can plague self-consistent semi-empirical methods on difficult systems (e.g., polar/charged proteins) and is roughly **2–20× faster** than GFN1/GFN2-xTB.
- Primarily aimed at very fast geometry pre-optimization, periodic/solid-state structures, and large-scale screening where some accuracy is sacrificed for robustness and speed.

### 2.2 GFN1-xTB
- The first broadly parametrized member of the family, using an effective **self-consistent-charge (SCC) tight-binding Hamiltonian** related to DFTB methodology but built from predominantly element-specific (not atom-pair-specific) parameters.
- Employs a Slater-type minimal valence basis (angular momenta up to *l* = 2) with **coordination-number-dependent** on-site energy levels — a distinctive feature that lets the same atomic parameters adapt to different bonding environments.
- Total energy: electronic (band-structure) + repulsion + D3(BJ) dispersion + explicit halogen-bonding term.
- Parametrized specifically to reproduce structures, vibrational frequencies, and non-covalent interaction energies across the periodic table (Z ≤ 86).

### 2.3 GFN2-xTB
- The most widely used member of the family and, as of the mid-2020s, the **default general-purpose method** in `xtb`.
- Improves on GFN1-xTB with:
  - **Anisotropic (multipole) electrostatics** — atomic monopole *and* dipole/quadrupole terms (isotropic electrostatics + exchange-correlation, IES+IXC, plus anisotropic electrostatics and exchange-correlation, AES+AXC).
  - A **self-consistent, charge-dependent D4 dispersion** correction (rather than a post-hoc D3 correction), obtained naturally from second-order density fluctuations in the tight-binding picture.
  - Removal of the explicit pairwise halogen-bond term (subsumed into the improved electrostatics/dispersion treatment) and a Fermi-smearing electronic entropy term for finite-temperature (e.g., high-T MD) simulations.
  - Fully **analytical nuclear gradients**.
- Benchmarked extensively against structures, non-covalent interaction energies, conformational energies, barrier heights, and dipole moments; shown to improve not only on the original "GFN" target properties but also on several "off-target" electronic-structure properties relative to GFN1-xTB, GFN-xTB, PM6-D3H4X, and DFTB3-D3(BJ).

### 2.4 GFN-FF
- A **fully classical, non-electronic** atomistic force field derived from and consistent with the GFN philosophy: it replaces the extended-Hückel electronic term with classical bond, angle, and torsion terms while retaining the same EEQ electrostatics and D4-type dispersion used in GFN0-xTB.
- Total energy split into covalent (bond/angle/torsion) and non-covalent (electrostatics, dispersion, hydrogen bonding, halogen bonding) contributions, with charge-dependent corrections for H-bonds and halogen bonds.
- Automatically generates its own topology/connectivity and force-field parameters "on the fly" from the input geometry, requiring **no manual atom-typing**, unlike most classical force fields (e.g., AMBER, CHARMM, UFF).
- Intended for very large systems (proteins, host–guest complexes, organometallic assemblies, materials) where GFN1/2-xTB would already be too slow, while remaining more physically grounded than generic FFs like UFF.

### 2.5 g-xTB (next-generation, general-purpose)
- A newer (2025) method from the same group, explicitly designed as a prospective **replacement for the whole GFNn-xTB line**, addressing systematic weaknesses of GFN2-xTB such as:
  - No Hartree–Fock exchange (GFN2-xTB behaves like a GGA-level functional, which is fine for geometries but weak for thermochemistry/barriers).
  - Minimal, environment-independent basis sets.
- Incorporates range-separated approximate Fock exchange and higher-order charge-fluctuation terms.
- Reported to roughly **halve mean absolute errors** relative to GFN2-xTB across many properties, reaching a WTMAD-2 of about 9.3 kcal/mol on the GMTKN55 benchmark set (competitive with low-cost DFT), with particularly large gains for transition-metal complexes, spin-state energetics, and orbital-energy gaps — while remaining only modestly slower (≲30%) than GFN2-xTB for single-point energies. As of its introduction, it lacks analytical gradients, making geometry optimization/frequencies noticeably more expensive than with GFN1/2-xTB.

### 2.6 Other tight-binding-based methods bundled with xtb
- **PTB (density-matrix / "polarized" tight-binding)** — targets more accurate molecular *properties* (atomic charges, bond orders, dipole moments, vibrational intensities) rather than geometries per se; also used together with GFN2-xTB for Raman spectrum calculations.
- **CEH (Charge-Extended Hückel model)** — a fast partial-charge model.
- **sTDA-xTB** — a simplified Tamm–Dancoff approximation built on the xTB Hamiltonian for cheap UV/vis excited-state calculations.
- **xTB-IFF / aISS** — an intermolecular force-field-based algorithm for automatically generating low-energy non-covalent complex/dimer/oligomer structures from monomer coordinates, including reactivity exploration and site-directed docking.

---

## 3. Core Functionality

The `xtb` binary/library exposes a broad set of calculation types built on top of any of the GFN Hamiltonians:

| Category | Capabilities |
|---|---|
| **Single-point properties** | Total energy, orbitals/eigenvalues, Mulliken/CM5-type and EEQ atomic charges, bond orders, dipole/multipole moments, HOMO–LUMO gap |
| **Structure** | Geometry optimization (internal-coordinate or Cartesian, with an approximate normal-coordinate rational function optimizer, ANCopt) |
| **Vibrational analysis** | Analytical/numerical Hessians, harmonic frequencies, IR intensities, thermochemistry (ZPE, enthalpy, entropy, free energy) |
| **Spectroscopy** | IR spectra, Raman activity (via GFN2-xTB/PTB combination), UV/vis excitation energies (sTDA-xTB) |
| **Dynamics** | Born–Oppenheimer molecular dynamics (MD); **metadynamics** for enhanced conformational sampling |
| **Multi-scale modeling** | ONIOM-style QM/QM' or QM/MM-like layered calculations |
| **Periodic systems** | Periodic boundary conditions (initially via GFN-xTB/GFN1-xTB extensions, later including GFN-FF) for crystals, surfaces, and frameworks such as MOFs |
| **Non-covalent complex generation** | Automated docking/oligomer building via the aISS algorithm |
| **Implicit solvation** | GBSA, ALPB, ddCOSMO, and CPCM-X (see §4) |
| **External-driver mode** | Can act as an energy/gradient engine driven by external optimizers or other QC codes |

### 3.1 Implicit solvation models

`xtb` supports several increasingly refined continuum solvation approaches:
- **GBSA** (Generalized Born + Solvent-Accessible Surface Area) — the original implicit-solvent option.
- **ALPB** (Analytical Linearized Poisson–Boltzmann) — introduced as the improved default from v6.3.3 onward, offering better accuracy than GBSA, with internally parametrized solvents and adjustable reference-state corrections for solvation free energies (e.g., a 1 mol/L "bar1mol" state, or a reference state matching COSMO-RS conventions).
- **ddCOSMO** (domain-decomposition COSMO) — a more scalable conductor-like screening model.
- **CPCM-X** — an extended conductor-like polarizable-continuum model parametrized against the Minnesota Solvation Database, added more recently; not yet fully compatible with all calculation types (e.g., has limitations for geometry optimization).

### 3.2 Companion tools in the Grimme-lab ecosystem

`xtb` is the computational engine behind, and closely integrated with, several sibling tools from the same group:
- **CREST** (Conformer–Rotamer Ensemble Sampling Tool) — automated conformer/rotamer/tautomer/protonation-state search built on xTB metadynamics.
- **CENSO** — conformer ensemble refinement/ranking with higher-level (DFT) single points.
- **QCG** — automated construction of explicit solvent shells/clusters around a solute, with conformational sampling, built on CREST.
- **QCxMS** — calculation of electron-impact mass spectra using xTB-driven "chemical" molecular dynamics.
- Third-party integration: xTB Hamiltonians are also available as engines inside other QC packages (e.g., an I/O-based interface has been available in ORCA since version 4.2).

---

## 4. Technical/Software Characteristics

- **Core implementation:** Fortran, with a documented **C API** enabling use as a library from other languages; community bindings exist for Python and other ecosystems.
- **Build systems:** supports both **Meson** and **CMake**.
- **Distribution:** available as source (GitHub), as statically linked precompiled binaries (Linux/Windows) on the GitHub releases page, via **conda-forge** (Linux x86_64/aarch64/ppc64le and macOS x86_64/arm64), and via a Homebrew tap (`grimme-lab/homebrew-qc`) for macOS.
- **Continuous releases:** a "bleeding edge" continuous-release tag tracks the latest `main`-branch builds.
- **License:** LGPL v3+ — free to redistribute and modify, distributed without warranty of merchantability or fitness for a particular purpose.
- **Community/maintainer base:** actively maintained by a substantial list of contributors from the Grimme group and collaborators (documented in the repository's contributor list), with ongoing feature releases (recent additions include COSMO file output, sandwich-potential models, periodic boundary conditions for GFN-FF, GFN-FF dipole moments, the DIPRO dimer-projection method, and QCSchema/QChem/FHI-aims geometry I/O support).
- **Documentation:** a dedicated Sphinx-based user guide (`xtb-docs`, hosted on Read the Docs) covers command-line usage, theory background, and worked examples for each major feature area.

---

## 5. Strengths

1. **Cost–accuracy balance.** Routinely 2–3+ orders of magnitude faster than DFT while retaining an explicit, self-consistent electronic-structure description (unlike pure force fields), making it suitable for conformational sampling, high-throughput screening, and MD of systems DFT cannot practically reach.
2. **Broad periodic-table coverage** (Z ≤ 86) using primarily element-specific (not pairwise) parameters — a key practical advantage over classical DFTB, which typically requires costly pairwise reparametrization for every new element combination.
3. **Robustness/"black-box" character.** Designed to run without manual re-tuning across diverse chemistries, including many transition-metal and organometallic systems, radicals, and charged species, and to remain numerically stable even where SCF convergence is difficult.
4. **Native gradients and Hessians.** Analytical gradients (GFN1/GFN2-xTB, GFN-FF) enable efficient geometry optimization, frequency calculations, and MD without finite-difference overhead.
5. **Integrated implicit solvation**, spanning GBSA through the more modern ALPB and CPCM-X models, with explicit solvation-free-energy reference-state handling.
6. **Rich ecosystem.** Tight integration with CREST, CENSO, QCG, and QCxMS extends xtb into conformer searching, ensemble refinement, explicit microsolvation, and mass-spectrometry prediction — use cases well beyond a typical single-point/optimization engine.
7. **Free, open-source, actively developed**, with multiple installation routes (conda-forge, Homebrew, prebuilt binaries, source builds) and a documented C API for embedding in other software.
8. **Increasing sophistication over time.** The GFN0 → GFN1 → GFN2 → GFN-FF → g-xTB progression shows a clear, benchmarked trajectory of addressing earlier methods' known weaknesses (e.g., g-xTB explicitly targets GFN2-xTB's underestimated barrier heights, compressed orbital gaps, and unreliable transition-metal geometries).

## 6. Limitations and Known Weaknesses

1. **Accuracy ceiling relative to DFT/wavefunction methods.** Despite the "DFT-level accuracy" branding, GFN-xTB methods remain semi-empirical approximations; GFN2-xTB in particular behaves similarly to a GGA-level functional (no exact exchange), which is known to systematically underestimate reaction barrier heights and compress HOMO–LUMO/orbital-energy gaps.
2. **Transition-metal and spin-state reliability.** While parametrized across most of the periodic table, GFN1/GFN2-xTB geometries and relative spin-state energetics for transition-metal complexes can be less reliable than for main-group chemistry, sometimes yielding distorted or unphysical structures — a primary motivation for developing g-xTB.
3. **GFN0-xTB accuracy trade-off.** Its non-self-consistent, non-iterative treatment (which gives it speed and robustness) comes with reduced accuracy, since higher-order energy terms present in GFN1/GFN2-xTB are omitted.
4. **g-xTB immaturity (as of introduction).** Lacks analytical gradients at release, making geometry optimization and frequency calculations substantially slower than with GFN1/GFN2-xTB; it is not yet the default production method for all workflows.
5. **Basis-set/methodological rigidity.** The Slater-type minimal-basis, element-specific parametrization approach — central to xtb's efficiency and periodic-table coverage — inherently limits systematic improvability compared to ab initio/DFT approaches with adjustable basis sets.
6. **Solvation model maturity varies.** Newer solvation models (e.g., CPCM-X) have documented limitations for certain calculation types (e.g., geometry optimization), and not all solvation models are available for every calculator (e.g., not yet implemented for the PTB and `tblite` external calculators as of recent releases).
7. **Periodic/solid-state support is a later addition** relative to the molecular-focus core methods, and periodic-system accuracy (e.g., for metal–organic frameworks) is generally good for structural parameters (~75% of cell parameters within 5% of experiment in one benchmark) but not uniformly at DFT accuracy, especially for less common metal environments.
8. **Learning curve/format-specific I/O quirks**, typical of specialized computational-chemistry command-line tools, including solvent-name case sensitivity, distinct old (`--gbsa`) vs. new (`--alpb`) solvation flags across versions, and the need to track which features are compatible with which calculators/methods.

---

## 7. Typical Use Cases

- Fast geometry pre-optimization ahead of higher-level DFT/wavefunction refinement.
- Conformer/rotamer/protonation-state ensemble generation (via CREST) for downstream free-energy or NMR/UV-vis property prediction.
- Molecular dynamics and metadynamics-based enhanced sampling of chemical reactivity, especially where DFT-MD would be too expensive.
- Vibrational (IR/Raman) spectrum prediction and thermochemistry for moderate-to-large molecules.
- High-throughput virtual screening of molecular libraries, materials candidates, or host–guest/binding-affinity studies.
- Structure prediction and screening for metal–organic frameworks and other periodic materials.
- As an engine inside larger workflows (e.g., ORCA's xtb interface, machine-learning potential training-data generation, automated reaction-network exploration).

---

## 8. Key Publications on xtb Theory and Methods

The following references cover the theoretical foundations, individual method developments, and companion-model publications underlying the `xtb` program package. DOIs are included where established.

### 8.1 Core GFNn-xTB methods

- Grimme, S.; Bannwarth, C.; Shushkov, P. "A Robust and Accurate Tight-Binding Quantum Chemical Method for Structures, Vibrational Frequencies, and Noncovalent Interactions of Large Molecular Systems Parametrized for All spd-Block Elements (Z = 1–86)." *J. Chem. Theory Comput.* **2017**, *13*, 1989–2009. (GFN1-xTB) DOI: 10.1021/acs.jctc.7b00118
- Bannwarth, C.; Ehlert, S.; Grimme, S. "GFN2-xTB — An Accurate and Broadly Parametrized Self-Consistent Tight-Binding Quantum Chemical Method with Multipole Electrostatics and Density-Dependent Dispersion Contributions." *J. Chem. Theory Comput.* **2019**, *15*, 1652–1671. DOI: 10.1021/acs.jctc.8b01176
- Pracht, P.; Caldeweyher, E.; Ehlert, S.; Grimme, S. "A Robust Non-Self-Consistent Tight-Binding Quantum Chemistry Method for large Molecules." *ChemRxiv* preprint, **2019**. (GFN0-xTB) DOI: 10.26434/chemrxiv.8326202
- Spicher, S.; Grimme, S. "Robust Atomistic Modeling of Materials, Organometallic, and Biochemical Systems." *Angew. Chem. Int. Ed.* **2020**, *59*, 15665–15673. (GFN-FF) DOI: 10.1002/anie.202004239

### 8.2 Review articles

- Bannwarth, C.; Caldeweyher, E.; Ehlert, S.; Hansen, A.; Pracht, P.; Seibert, J.; Spicher, S.; Grimme, S. "Extended Tight-Binding Quantum Chemistry Methods." *WIREs Comput. Mol. Sci.* **2021**, *11*, e1493. DOI: 10.1002/wcms.1493
- Katbashev, A.; Stahn, M.; Rose, T.; Alizadeh, V.; Friede, M.; Plett, C.; Steinbach, P.; Ehlert, S. "Overview on Building Blocks and Applications of Efficient and Robust Extended Tight-Binding Methods." (Extended review of the xtb software ecosystem; check publisher listing for final journal/year.)

### 8.3 Solvation models

- Ehlert, S.; Stahn, M.; Spicher, S.; Grimme, S. "Robust and Efficient Implicit Solvation Model for Fast Semiempirical Methods." *J. Chem. Theory Comput.* **2021**, *17*, 4250–4261. (GBSA and ALPB) DOI: 10.1021/acs.jctc.1c00471
- Stahn, M.; Ehlert, S.; Grimme, S. "Extended Conductor-like Polarizable Continuum Solvation Model (CPCM-X) for Semiempirical Methods." *J. Phys. Chem. A* **2023**, *127*, 7036–7043. DOI: 10.1021/acs.jpca.3c04382

### 8.4 Dispersion correction (DFT-D4 / D4 model, used within GFN2-xTB and GFN-FF)

- Caldeweyher, E.; Bannwarth, C.; Grimme, S. "Extension of the D3 Dispersion Coefficient Model." *J. Chem. Phys.* **2017**, *147*, 034112. DOI: 10.1063/1.4993215
- Caldeweyher, E.; Ehlert, S.; Hansen, A.; Neugebauer, H.; Spicher, S.; Bannwarth, C.; Grimme, S. "A Generally Applicable Atomic-Charge Dependent London Dispersion Correction." *J. Chem. Phys.* **2019**, *150*, 154122. DOI: 10.1063/1.5090222
- Caldeweyher, E.; Mewes, J.-M.; Ehlert, S.; Grimme, S. "Extension and Evaluation of the D4 London-Dispersion Model for Periodic Systems." *Phys. Chem. Chem. Phys.* **2020**, *22*, 8499–8512. DOI: 10.1039/D0CP00502A

### 8.5 Excited states / spectroscopy

- Grimme, S.; Bannwarth, C. "Ultra-Fast Computation of Electronic Spectra for Large Systems by Tight-Binding Based Simplified Tamm–Dancoff Approximation (sTDA-xTB)." *J. Chem. Phys.* **2016**, *145*, 054103. (sTDA-xTB) DOI: 10.1063/1.4959605

### 8.6 Metadynamics / conformer sampling / mass spectrometry (CREST/QCxMS ecosystem)

- Grimme, S. "Exploration of Chemical Compound, Conformer, and Reaction Space with Meta-Dynamics Simulations Based on Tight-Binding Quantum Chemical Calculations." *J. Chem. Theory Comput.* **2019**, *15*, 2847–2862. DOI: 10.1021/acs.jctc.9b00143
- Pracht, P.; Bohle, F.; Grimme, S. "Automated Exploration of the Low-Energy Chemical Space with Fast Quantum Chemical Methods." *Phys. Chem. Chem. Phys.* **2020**, *22*, 7169–7192. (CREST) DOI: 10.1039/C9CP06869D
- Asgeirsson, V.; Bauer, C. A.; Grimme, S. "Quantum Chemical Calculation of Electron Ionization Mass Spectra for General Organic and Inorganic Molecules." *Chem. Sci.* **2017**, *8*, 4879–4895. (QCxMS) DOI: 10.1039/C7SC00601B

### 8.7 Next-generation method

- Grimme, S. et al. "g-xTB: A General-Purpose Extended Tight-Binding Electronic Structure Method." *ChemRxiv* preprint, **2025**. DOI: 10.26434/chemrxiv-2025-bjxvt

### 8.8 Foundational precursor (GFN-xTB, the original 2017 model this family builds on)

- Grimme, S.; Bannwarth, C.; Shushkov, P. (see §8.1 — this is the same 2017 paper commonly cited as "GFN-xTB"/"GFN1-xTB"; the naming "GFN1-xTB" became standard once GFN2-xTB was introduced to disambiguate versions.)

> **Note on citation practice:** the `xtb` program itself prints/maintains a version-specific list of "please cite" references (available in BibTeX form from the repository) that should be consulted for the exact set of papers relevant to whichever GFN Hamiltonian, dispersion correction, and solvation model are used in a given calculation, since the appropriate citation set depends on which features are invoked.

---

## 9. Summary Assessment

`xtb` is a mature, actively developed, and widely adopted open-source semi-empirical quantum chemistry package that has become a de facto standard for fast, broadly applicable tight-binding calculations across the periodic table. Its central value proposition — an explicit, physically motivated electronic-structure method that is dramatically cheaper than DFT yet more transferable and robust than classical force fields — is well supported by a large and growing benchmark literature. Its main limitations mirror those of any semi-empirical approach: an intrinsic accuracy ceiling relative to correlated wavefunction/DFT methods, particular caution warranted for transition-metal and spin-state problems, and periodic uneven maturity across newer features (periodic boundary conditions, newer solvation models, g-xTB). The ongoing development trajectory (GFN0 → GFN1 → GFN2 → GFN-FF → g-xTB, alongside the CREST/CENSO/QCG/QCxMS ecosystem) indicates a deliberate, benchmark-driven effort to close the remaining accuracy gaps while preserving the method family's defining speed and generality advantages.


---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the xtb (Grimme group) 	Open-source semi-empirical extended tight-binding package built to approximate DFT-level accuracy at very low computational cost. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
