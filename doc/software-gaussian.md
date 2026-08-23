# Gaussian: An Exhaustive Review of the Commercial Quantum Chemistry Package

## 1. Overview

Gaussian is a commercial electronic structure program originally developed by John Pople and coworkers at Carnegie Mellon University, first released as **Gaussian 70**. It has been continuously developed for over five decades and is now produced and licensed by Gaussian, Inc. (Wallingford, CT). It is one of the most widely used quantum chemistry packages in academia and industry, spanning chemistry, biochemistry, chemical engineering, physics, and materials science. The current production release is **Gaussian 16**, distributed with the companion GUI **GaussView 6** for job setup and visualization. No successor ("Gaussian 17") has been publicly announced as of 2026; Gaussian 16 continues to receive periodic revisions.

Gaussian starts from the fundamental equations of quantum mechanics to predict energies, molecular structures, vibrational frequencies, and a broad range of molecular properties for stable species, reactive intermediates, and transition states, across essentially the whole periodic table.

**Key facts**
- **Original author / origin:** John Pople and coworkers, Carnegie Mellon University (1970)
- **Current publisher:** Gaussian, Inc., Wallingford, Connecticut, USA
- **Current version:** Gaussian 16 (companion GUI: GaussView 6)
- **License model:** Commercial/proprietary; academic and commercial pricing tiers; parallel execution via the Linda interconnect for distributed-memory clusters
- **Platforms:** Linux/Unix, Windows, macOS
- **Notable ecosystem integration:** OpenEye Orion® cloud platform (Gaussian Module) for cloud-scale quantum chemistry workflows

---

## 2. Core Architecture

Gaussian is organized internally as a series of modular executables called **"links"** (e.g., L101, L502, L913…) that are chained together by an overlay-driven control program. Each link handles a discrete task: reading the route section, generating molecular integrals, running SCF/DFT iterations, computing correlation energies, evaluating analytic derivatives, and so on. This link-based architecture underlies Gaussian's flexibility in combining job types (single point, optimization, frequency, IRC, dynamics) with virtually any supported method/basis-set combination — a pairing Gaussian calls a **"model chemistry."**

Selected components:
- **L502/L503/L508:** SCF solvers (conventional, direct, quadratically convergent)
- **L510/L405:** MCSCF/CASSCF machinery
- **L601/L607:** Population analysis, NBO
- **L913/L914/L923:** Post-SCF energies, CI-Singles/RPA/ZINDO excited states, SAC-CI
- **L716/L1101–L1112:** Derivative and second-derivative (Hessian/frequency) evaluation
- **L120:** ONIOM control

---

## 3. Density Functional Theory (DFT) Capabilities

Gaussian offers one of the broadest DFT functional libraries of any general-purpose package, spanning every rung of "Jacob's Ladder" — from LDA through GGA, meta-GGA, hybrid, range-separated/long-range-corrected, and double-hybrid functionals — plus empirical dispersion corrections and user-configurable hybrids.

### 3.1 Exchange functionals
Slater (Xα/LSDA), Becke 88 (B), Perdew–Wang 91 (PW91), Barone-modified PW91 (mPW), Gill 96 (G96), PBE, OPTX, TPSS, revised TPSS, Becke–Roussel (BRx), PKZB, ωPBEh/HSE, PBEh.

### 3.2 Correlation functionals
VWN (VWN3/VWN5), LYP, Perdew 81 (P81), Perdew 86 (P86), Perdew–Wang 91 (PW91), PBE, Becke95 (B95), TPSS, revised TPSS, KCIS, BRC, PKZB, VP86, V5LYP.

### 3.3 Other "pure" (non-hybrid) functionals
VSXC, the HCTH functional family, τ-HCTH, B97-D, M06-L, SOGGA11, M11-L, MN12-L, N12, MN15-L.

### 3.4 Hybrid functionals
B3LYP, B3P86, B3PW91, "B1" and variants, B98, B97-1, B97-2, PBE1PBE (PBE0), HSEh1PBE and variants (HSE06-type screened hybrids), O3LYP, TPSSh, τ-HCTH-hybrid, BMK, APF/APFD, M05, M05-2X, M06, M06-HF, M06-2X, M08-HX, PW6B95, PW6B95D3, M11, SOGGA11-X, N12-SX, MN12-SX, MN15, HISS/HISSbPBE, X3LYP, BHandHLYP, plus **fully user-configurable hybrid functionals** (arbitrary linear combinations of HF exchange with exchange/correlation functionals).

### 3.5 Long-range-corrected / range-separated functionals
CAM-B3LYP (Coulomb-attenuating method), ωB97, ωB97X, ωB97X-D (with empirical dispersion), LC-ωPBE, and the general **LC-** prefix, which can apply Hirao-style long-range correction to most pure functionals (e.g., LC-BLYP).

### 3.6 Double-hybrid functionals
B2PLYP, mPW2PLYP (and dispersion-corrected variants), DSD-PBEP86, PBE0-DH, PBE-QIDH.

### 3.7 Empirical dispersion corrections
PFD, Grimme's GD2 (DFT-D2), GD3 (DFT-D3, zero-damping), GD3BJ (DFT-D3 with Becke–Johnson damping); functionals with dispersion built in by default include APFD, B97D3, ωB97X-D, and B2PLYP-D3. The standalone `EmpiricalDispersion` keyword allows any supported dispersion scheme to be layered onto compatible functionals.

### 3.8 DFT numerical-integration and fitting infrastructure
Multiple pruned/unpruned integration grid qualities; Lebedev angular quadratures; auxiliary **density-fitting sets** (DGA1, DGA2, W06, and legacy sets tuned for SVP/TZVP basis sets) with optional automatic generation and default resolution-of-identity-style acceleration.

### 3.9 Time-Dependent DFT (TD-DFT)
Vertical excitation energies, oscillator strengths, and analytic gradients for excited-state geometry optimization; state-specific and linear-response solvation coupling (see §6); vibronic/Franck–Condon and Herzberg–Teller band-shape simulation for absorption, fluorescence, and resonance-Raman spectra.

---

## 4. Hartree–Fock and Post-Hartree-Fock (Wavefunction) Methods

Gaussian supports the full traditional ab initio hierarchy, with restricted (R), unrestricted (U), and restricted-open-shell (RO) variants available for most methods.

### 4.1 Mean-field
- **HF** (default method if none is specified): RHF, UHF, ROHF
- **GVB-PP** (generalized valence bond, perfect pairing)

### 4.2 Perturbation theory
- **MP2, MP3, MP4** (SDQ and SDTQ variants), **MP5**
- Direct and semi-direct MP2 algorithms; analytic MP2 gradients and frequencies
- Spin-restricted (RMP) and spin-unrestricted post-HF perturbation theory for open-shell references

### 4.3 Configuration interaction
- **CIS** (configuration interaction singles) — with analytic gradients and second derivatives, used for excited states
- **CISD**, **QCISD** (quadratic CI, size-consistency corrected), **QCISD(T)**
- **CASSCF/RASSCF** (complete/restricted active-space SCF) multi-reference methods, optionally combined with MP2 for dynamic correlation on top of a CAS reference

### 4.4 Coupled cluster
- **CCD, CCSD, CCSD(T)** ("gold standard" single-reference method), **BD** (Brueckner Doubles), **BD(T)**
- Analytic gradients for geometry optimization at the CCSD level
- **EOM-CCSD** (equation-of-motion CC) for excited states, including solvated (PCM-coupled) variants

### 4.5 Electron propagator / Green's function methods
- **EPT** (electron propagator theory), including outer valence Green's function (OVGF)-type approaches and higher-order **ADC(3)**-related models for ionization potentials and electron affinities

### 4.6 Excited-state / spectroscopy-specific wavefunction methods
- **ZINDO** (semi-empirical INDO/S) for UV/Vis spectra
- **SAC-CI** (symmetry-adapted-cluster configuration interaction) for excited, ionized, and electron-attached states, with analytic gradients

### 4.7 Composite ("thermochemical recipe") methods
Combine multiple levels of theory to approach CCSD(T)/complete-basis-set-quality thermochemistry at greatly reduced cost:
- **G1, G2, G3, G4** and reduced-order variants (**G2(MP2)**, **G3(MP2)**, **G4(MP2)**), plus **G3B3/G4** using DFT geometries and zero-point energies
- **CBS-4, CBS-Q, CBS-QB3, CBS-APNO** (Petersson's complete-basis-set extrapolation family)
- **W1U, W1BD** (Martin-style Weizmann-1 composite methods, unrestricted and Brueckner-doubles based)

### 4.8 Semi-empirical methods
AM1, PM3, PM3MM, PM6, PDDG, MNDO, MINDO/3, and related parameterizations — useful for large systems or as initial guesses/geometries preceding higher-level single points.

### 4.9 Molecular mechanics and hybrid QM methods
- Built-in molecular mechanics (e.g., AMBER-type force fields) for MM and mixed QM/MM regions
- **ONIOM** (two- or three-layer "our own N-layer integrated molecular orbital and molecular mechanics") for treating large systems (including systems with thousands of atoms) by partitioning them into regions of differing levels of theory (e.g., QM:QM or QM:MM), with full support for energies, analytic gradients, frequencies, and electric-field derivatives

### 4.10 Relativistic treatments
Douglas–Kroll–Hess (DKH) scalar relativistic Hamiltonians; relativistic effective core potentials (pseudopotentials) for heavy elements (transition metals, lanthanides, actinides).

---

## 5. Basis Sets

Gaussian ships with an extensive internal basis-set library and can read externally supplied sets:
- **Pople-style split-valence sets**: STO-3G, 3-21G, 6-21G, 4-31G, 6-31G, 6-311G, with the standard diffuse (`+`/`++`) and polarization (`*`/`**`, or explicit `(d,p)`, `(3df,2p)`, etc.) augmentations
- **Dunning correlation-consistent sets**: cc-pVDZ/TZ/QZ/5Z and augmented (aug-cc-pV*Z) variants, designed for systematic basis-set extrapolation to the complete basis set (CBS) limit
- **Effective core potentials (ECPs) / pseudopotentials**: LanL2DZ, LanL2MB, SDD (Stuttgart–Dresden), CEP-4G/31G/121G, and energy-consistent pseudopotentials for transition metals, lanthanides, and actinides
- **Density-fitting/auxiliary basis sets** for accelerated Coulomb evaluation in DFT
- Support for basis functions and ECPs of **arbitrary angular momentum**, and for user-defined/general basis input (`Gen`/`GenECP`)
- Program limits: up to 250,000 atoms, 750,000 primitive shells, and 250,000 contracted shells (integral engine limits)

---

## 6. Solvation and Environmental Models

- **PCM (Polarizable Continuum Model)** family: standard/dielectric PCM, **C-PCM** (conductor-like), **IEF-PCM** (integral equation formalism), IPCM
- **SMD** (universal solvation model based on solute electron density, Truhlar/Cramer/Marenich) for solvation free energies across arbitrary solvents
- Analytic energy, gradient, and (for many methods) Hessian derivatives in solution, including at correlated (MP2, CCSD, EOM-CCSD) and TD-DFT levels
- State-specific and linear-response TD-DFT/PCM coupling for excited states in solution, including non-equilibrium solvation effects

---

## 7. Job Types and Molecular Properties

### 7.1 Job types
Single-point energy (SP), geometry optimization (Opt: minima, transition structures, higher-order saddle points), frequency/thermochemistry analysis (Freq, including anharmonic VPT2 treatment), intrinsic reaction coordinate following (IRC), IRC-max, potential energy surface scans (Scan), polarizability/hyperpolarizability (Polar), direct ab initio molecular dynamics (**ADMP** and **BOMD**), excitation energy transfer (EET), force evaluation, wavefunction stability testing (Stable), molecular volume, population-analysis-only reruns, and combined `method2/basis2 // method1/basis1` optimize-then-single-point workflows.

### 7.2 Geometry optimization algorithms
Berny algorithm, GDIIS, Newton–Raphson, Fletcher–Powell, Murtaugh–Sargent, eigenvector-following (EF); optimizations in redundant internal, Z-matrix (internal), Cartesian, or mixed coordinates. Transition-state searches via Synchronous Transit-Guided Quasi-Newton (**QST2/QST3**) and linear-synchronous-transit methods.

### 7.3 Representative molecular properties accessible
- IR, Raman (including pre-resonance and resonance Raman), VCD (vibrational circular dichroism), and ROA (Raman optical activity) spectra, both harmonic and anharmonic
- UV/Vis absorption and electronic circular dichroism (via CIS, ZINDO, TD-DFT, EOM-CC, SAC-CI)
- Vibronic band shapes via Franck–Condon and Herzberg–Teller analysis
- NMR shielding/chemical shifts and spin–spin coupling constants (GIAO-based)
- Atomic charges (Mulliken, natural population analysis/NBO, CHelpG, Merz–Kollman electrostatic-potential-derived, Hirshfeld)
- Dipole/multipole moments, polarizabilities and hyperpolarizabilities, optical rotation
- Thermochemistry (enthalpies, entropies, Gibbs free energies) from harmonic and anharmonic vibrational analysis
- Electron density and molecular orbitals on cube-file grids (`cubegen`) for visualization
- Electron affinities, ionization potentials (via composite methods, CCSD, or electron propagator theory)
- Hyperfine coupling constants, antiferromagnetic coupling (via fragment-guess/stability analysis)

---

## 8. Software Ecosystem

- **GaussView 6**: graphical front end for building molecules, generating input files, and visualizing results (orbitals, densities, spectra, vibrations, NBOs).
- **Linda**: distributed-memory parallel interconnect enabling Gaussian jobs to scale across multiple networked nodes.
- **GMMX**: conformational searching add-on.
- Native **NBO (Natural Bond Orbital)** analysis engine bundled with the program.
- Third-party interoperability: output is text-based and widely parsed by external tools (e.g., **cclib**), and Gaussian can serve as the QM engine within third-party QM/MM or cloud workflows — notably the **OpenEye Orion® Gaussian Module**, which exposes Gaussian's methods as scalable cloud workflows (Floes) on AWS infrastructure for pharmaceutical and materials-science R&D.

---

## 9. Strengths, Limitations, and Positioning

**Strengths**
- Breadth of methods and functionals is arguably unmatched among general-purpose commercial packages, letting a single program cover routine DFT screening, high-accuracy composite thermochemistry, multi-reference problems, and large hybrid QM/MM systems.
- Mature, heavily validated implementations with an enormous, decades-deep literature base (see §10) underpinning nearly every feature.
- Extensive spectroscopic property support (vibrational, NMR, chiroptical) with both harmonic and anharmonic treatments.
- Robust solvation modeling coupled to correlated wavefunction and TD-DFT methods.

**Limitations / considerations**
- Commercial license cost can be a barrier relative to open-source alternatives (e.g., Psi4, NWChem, ORCA — the latter free for academics).
- Some accuracy caveats are inherent to the methods rather than the software: e.g., B3LYP is well known to be strong for geometries but not always optimal for main-group thermochemistry, where double hybrids or other functionals may perform better; standard (non-dispersion-corrected) DFT does not capture long-range van der Waals interactions without an explicit `EmpiricalDispersion` correction.
- CCSD(T) and other high-level post-HF methods remain computationally expensive and scale poorly to very large systems relative to local-correlation approaches (e.g., DLPNO-CCSD(T)) found in some competing codes.
- Periodic/solid-state electronic structure is not a primary focus (Gaussian is molecule-centric); users needing full periodic DFT typically pair Gaussian with, or switch to, dedicated solid-state codes.

**Notable competing/complementary packages**: ORCA, Psi4, NWChem, GAMESS, Q-Chem, Molcas/OpenMolcas, Jaguar, Spartan, TURBOMOLE, ADF, CRYSTAL (periodic).

---

## 10. Key Publications on Gaussian's Underlying Theory

The following list highlights the primary literature underpinning major methods and features implemented in Gaussian. Gaussian, Inc. maintains a comprehensive, continuously updated bibliography of several hundred references on its website; the entries below are a representative, categorized subset spanning the program's foundational and most consequential theoretical developments.

### 10.1 Program citations (the software itself)
- Frisch, M. J. et al. *Gaussian 16*, Revision C.01 (and successors), Gaussian, Inc., Wallingford, CT.
- Frisch, M. J. et al. *Gaussian 09*, Gaussian, Inc., Wallingford, CT, 2009.
- Frisch, M. J. et al. *Gaussian 03*, Gaussian, Inc., Wallingford, CT, 2003.
- Hehre, W. J.; Lathan, W. A.; Ditchfield, R.; Newton, M. D.; Pople, J. A. *Gaussian 70*, Quantum Chemistry Program Exchange, Program No. 237 (1970).
- Foresman, J. B.; Frisch, Æ. *Exploring Chemistry with Electronic Structure Methods*, 3rd ed., Gaussian, Inc., Wallingford, CT, 2015.

### 10.2 Density functional theory — foundations and exchange-correlation functionals
- Hohenberg, P.; Kohn, W. "Inhomogeneous Electron Gas." *Phys. Rev.* **136**, B864 (1964).
- Kohn, W.; Sham, L. J. "Self-Consistent Equations Including Exchange and Correlation Effects." *Phys. Rev.* **140**, A1133 (1965).
- Becke, A. D. "Density-functional exchange-energy approximation with correct asymptotic behavior." *Phys. Rev. A* **38**, 3098 (1988).
- Becke, A. D. "A new mixing of Hartree–Fock and local density-functional theories." *J. Chem. Phys.* **98**, 1372 (1993).
- Becke, A. D. "Density-functional thermochemistry. III. The role of exact exchange." *J. Chem. Phys.* **98**, 5648 (1993).
- Lee, C.; Yang, W.; Parr, R. G. "Development of the Colle-Salvetti correlation-energy formula into a functional of the electron density." *Phys. Rev. B* **37**, 785 (1988).
- Perdew, J. P.; Burke, K.; Ernzerhof, M. "Generalized Gradient Approximation Made Simple" (PBE functional context); Ernzerhof, M.; Scuseria, G. E. "Assessment of the Perdew–Burke–Ernzerhof exchange-correlation functional." *J. Chem. Phys.* **110**, 5029 (1999).
- Adamo, C.; Barone, V. "Toward reliable density functional methods without adjustable parameters: The PBE0 model." *J. Chem. Phys.* **110**, 6158 (1999).
- Heyd, J.; Scuseria, G. E.; Ernzerhof, M. "Hybrid functionals based on a screened Coulomb potential." *J. Chem. Phys.* **118**, 8207 (2003) (HSE functional).
- Zhao, Y.; Truhlar, D. G. "The M06 suite of density functionals..." *Theor. Chem. Acc.* **120**, 215 (2008) (M06 family).
- Yanai, T.; Tew, D. P.; Handy, N. C. "A new hybrid exchange–correlation functional using the Coulomb-attenuating method (CAM-B3LYP)." *Chem. Phys. Lett.* **393**, 51 (2004).
- Chai, J.-D.; Head-Gordon, M. "Systematic optimization of long-range corrected hybrid density functionals." *J. Chem. Phys.* **128**, 084106 (2008).
- Chai, J.-D.; Head-Gordon, M. "Long-range corrected hybrid density functionals with damped atom–atom dispersion corrections." *Phys. Chem. Chem. Phys.* **10**, 6615 (2008) (ωB97X-D).
- Grimme, S. "Semiempirical GGA-type density functional constructed with a long-range dispersion correction." *J. Comp. Chem.* **27**, 1787 (2006).
- Grimme, S.; Antony, J.; Ehrlich, S.; Krieg, H. "A consistent and accurate ab initio parametrization of density functional dispersion correction (DFT-D) for the 94 elements H–Pu." *J. Chem. Phys.* **132**, 154104 (2010) (DFT-D3).
- Grimme, S.; Ehrlich, S.; Goerigk, L. "Effect of the damping function in dispersion corrected density functional theory." *J. Comp. Chem.* **32**, 1456 (2011) (Becke–Johnson damping).
- Grimme, S. "Semiempirical hybrid density functional with perturbative second-order correlation." *J. Chem. Phys.* **124**, 034108 (2006) (B2PLYP double hybrid).
- Kozuch, S.; Martin, J. M. L. "DSD-PBEP86: In search of the best double-hybrid DFT..." *Phys. Chem. Chem. Phys.* **13**, 20104 (2011).
- Brémond, É.; Adamo, C. "Seeking for parameter-free double-hybrid functionals: The PBE0-DH model." *J. Chem. Phys.* **135**, 024106 (2011).
- Austin, A.; Petersson, G. A.; Frisch, M. J.; Dobek, F. J.; Scalmani, G.; Throssell, K. "A density functional with spherical atom dispersion terms." *J. Chem. Theory Comput.* **8**, 4989 (2012) (APFD).

### 10.3 Hartree–Fock, integrals, and SCF algorithms
- Almlöf, J.; Korsell, K.; Fægri, K. Jr. "Principles for a direct SCF approach to LCAO-MO ab initio calculations." *J. Comp. Chem.* **3**, 385 (1982).
- Dupuis, M.; Rys, J.; King, H. F. "Evaluation of molecular integrals over Gaussian basis functions." *J. Chem. Phys.* **65**, 111 (1976).
- Head-Gordon, M.; Pople, J. A. "A Method for Two-Electron Gaussian Integral and Integral Derivative Evaluation Using Recurrence Relations." *J. Chem. Phys.* **89**, 5777 (1988).
- Bacskay, G. B. "A Quadratically Convergent Hartree-Fock (QC-SCF) Method." *Chem. Phys.* **61**, 385 (1981).
- Kudin, K. N.; Scuseria, G. E.; Cancès, E. "A black-box self-consistent field convergence algorithm." *J. Chem. Phys.* **116**, 8255 (2002).
- Burant, J. C.; Scuseria, G. E.; Frisch, M. J. "Linear scaling method for Hartree-Fock exchange calculations of large molecules." *J. Chem. Phys.* **105**, 8969 (1996) (Gaussian Very Fast Multipole Method).

### 10.4 Post-Hartree-Fock wavefunction methods
- Møller, C.; Plesset, M. S. (foundational MP perturbation theory); Head-Gordon, M.; Pople, J. A.; Frisch, M. J. "MP2 energy evaluation by direct methods." *Chem. Phys. Lett.* **153**, 503 (1988).
- Frisch, M. J.; Head-Gordon, M.; Pople, J. A. "Direct MP2 gradient method." *Chem. Phys. Lett.* **166**, 275 (1990).
- Cížek, J. (foundational coupled-cluster theory), in *Advances in Chemical Physics*, Vol. 14 (1969).
- Purvis, G. D.; Bartlett, R. J. (CCSD formulation); Bartlett, R. J.; Purvis, G. D. "Many-body perturbation theory, coupled-pair many-electron theory, and importance of quadruple excitations for correlation problem." *Int. J. Quantum Chem.* **14**, 561 (1978).
- Lee, T. J.; Rendell, A. P.; Taylor, P. R. "Comparison of the Quadratic Configuration Interaction and Coupled-Cluster Approaches to Electron Correlation Including the Effect of Triple Excitations." *J. Phys. Chem.* **94**, 5463 (1990) (CCSD(T)).
- Handy, N. C.; Pople, J. A.; Head-Gordon, M.; Raghavachari, K.; Trucks, G. W. "Size-consistent Brueckner theory limited to double substitutions." *Chem. Phys. Lett.* **164**, 185 (1989) (BD method).
- Foresman, J. B.; Head-Gordon, M.; Pople, J. A.; Frisch, M. J. "Toward a Systematic Molecular Orbital Theory for Excited States." *J. Phys. Chem.* **96**, 135 (1992) (CIS excited states).
- Cederbaum, L. S. "One-body Green's function for atoms and molecules: Theory and application." *J. Phys. B* **8**, 290 (1975) (electron propagator theory).
- Nakatsuji, H. and coworkers, SAC/SAC-CI method development (multiple papers, e.g., Fukuda, R.; Nakatsuji, H. "Formulation and implementation of direct algorithm for the symmetry adapted cluster and symmetry adapted cluster-configuration interaction method." *J. Chem. Phys.* **128**, 094105 (2008)).

### 10.5 Composite thermochemical methods
- Curtiss, L. A.; Raghavachari, K.; Trucks, G. W.; Pople, J. A. "Gaussian-2 theory for molecular energies of first- and second-row compounds." *J. Chem. Phys.* **94**, 7221 (1991).
- Curtiss, L. A.; Raghavachari, K.; Redfern, P. C.; Rassolov, V.; Pople, J. A. "Gaussian-3 (G3) theory for molecules containing first and second-row atoms." *J. Chem. Phys.* **109**, 7764 (1998).
- Curtiss, L. A.; Redfern, P. C.; Raghavachari, K. "Gaussian-4 theory." *J. Chem. Phys.* **126**, 084108 (2007).
- Curtiss, L. A.; Redfern, P. C.; Raghavachari, K. "Gaussian-4 theory using reduced order perturbation theory." *J. Chem. Phys.* **127**, 124105 (2007) (G4(MP2)).
- Baboul, A. G.; Curtiss, L. A.; Redfern, P. C.; Raghavachari, K. "Gaussian-3 theory using density functional geometries and zero-point energies." *J. Chem. Phys.* **110**, 7650 (1999).
- Montgomery, J. A.; Frisch, M. J.; Ochterski, J. W.; Petersson, G. A. (CBS-QB3 and related complete basis set methods).
- Barnes, E. C.; Petersson, G. A.; Montgomery, J. A.; Frisch, M. J.; Martin, J. M. L. "Unrestricted Coupled Cluster and Brueckner Doubles Variations of W1 Theory." *J. Chem. Theory Comput.* **5**, 2687 (2009).

### 10.6 Geometry optimization and reaction-path algorithms
- Schlegel, H. B. "Berny" geometry optimization algorithm foundational work; Fogarasi, G.; Zhou, X.; Taylor, P.; Pulay, P. "The calculation of ab initio molecular geometries: Efficient optimization by natural internal coordinates." *J. Am. Chem. Soc.* **114**, 8191 (1992).
- Peng, C.; Ayala, P. Y.; Schlegel, H. B.; Frisch, M. J. (redundant internal coordinate optimization).
- Fukui, K. "The path of chemical reactions – The IRC approach." *Acc. Chem. Res.* **14**, 363 (1981).
- Gonzalez, C.; Schlegel, H. B. "An Improved Algorithm for Reaction Path Following." *J. Chem. Phys.* **90**, 2154 (1989); "Reaction Path Following in Mass-Weighted Internal Coordinates." *J. Phys. Chem.* **94**, 5523 (1990).
- Peng, C.; Schlegel, H. B. (Synchronous Transit-Guided Quasi-Newton, QST2/QST3, transition-state search methods).
- Farkas, Ö.; Schlegel, H. B. "Methods for geometry optimization of large molecules." *J. Chem. Phys.* **109**, 7100 (1998) and **111**, 10806 (1999).

### 10.7 ONIOM and QM/MM
- Dapprich, S.; Komáromi, I.; Byun, K. S.; Morokuma, K.; Frisch, M. J. "A New ONIOM Implementation in Gaussian 98. 1. The Calculation of Energies, Gradients and Vibrational Frequencies and Electric Field Derivatives." *J. Mol. Struct. (Theochem)* **462**, 1 (1999).
- Humbel, S.; Sieber, S.; Morokuma, K. "The IMOMO method..." *J. Chem. Phys.* **105**, 1959 (1996).

### 10.8 Solvation models
- Tomasi, J.; Persico, M. (foundational PCM development); Cancès, E.; Mennucci, B.; Tomasi, J. "A new integral equation formalism for the polarizable continuum model." *J. Chem. Phys.* **107**, 3032 (1997) (IEF-PCM).
- Cossi, M.; Rega, N.; Scalmani, G.; Barone, V. "Energies, structures, and electronic properties of molecules in solution with the C-PCM solvation model." *J. Comp. Chem.* **24**, 669 (2003).
- Marenich, A. V.; Cramer, C. J.; Truhlar, D. G. "Universal solvation model based on solute electron density..." (SMD model).
- Improta, R.; Barone, V.; Scalmani, G.; Frisch, M. J. "A state-specific polarizable continuum model time dependent density functional method for excited state calculations in solution." *J. Chem. Phys.* **125**, 054103 (2006).
- Caricato, M. et al. Multiple papers on CCSD-PCM and EOM-CCSD-PCM (2011–2014), extending correlated wavefunction methods into solution.

### 10.9 Vibrational spectroscopy, anharmonicity, and chiroptical properties
- Barone, V. "Anharmonic vibrational properties by a fully automated second-order perturbative approach." *J. Chem. Phys.* **122**, 014108 (2005) (VPT2).
- Bloino, J.; Biczysko, M.; Barone, V. "General perturbative approach for spectroscopy, thermodynamics and kinetics: Methodological background and benchmark studies." *J. Chem. Theory Comput.* **8**, 1015 (2012).
- Baiardi, A.; Bloino, J.; Barone, V. "General Time Dependent Approach to Vibronic Spectroscopy Including Franck-Condon, Herzberg-Teller and Duschinsky Effects." *J. Chem. Theory Comput.* **9**, 4097 (2013).
- Cheeseman, J. R.; Frisch, M. J.; Devlin, F. J.; Stephens, P. J. "Ab Initio Calculation of Atomic Axial Tensors and Vibrational Rotational Strengths Using Density Functional Theory." *Chem. Phys. Lett.* **252**, 211 (1996) (VCD).
- Cheeseman, J. R.; Frisch, M. J. "Basis Set Dependence of Vibrational Raman and Raman Optical Activity Intensities." *J. Chem. Theory Comput.* **7**, 3323 (2011).

### 10.10 NMR properties
- Ditchfield, R. "Self-consistent perturbation theory of diamagnetism. I. Gauge-invariant LCAO method for NMR chemical shifts." *Mol. Phys.* **27**, 789 (1974) (GIAO method).
- Cheeseman, J. R.; Trucks, G. W.; Keith, T. A.; Frisch, M. J. "A Comparison of Models for Calculating Nuclear Magnetic Resonance Shielding Tensors." *J. Chem. Phys.* **104**, 5497 (1996).
- Deng, W.; Cheeseman, J. R.; Frisch, M. J. "Calculation of Nuclear Spin-Spin Coupling Constants of Molecules with First and Second Row Atoms." *J. Chem. Theory Comput.* **2**, 1028 (2006).

### 10.11 Basis sets and pseudopotentials
- Ditchfield, R.; Hehre, W. J.; Pople, J. A. "Self-Consistent Molecular Orbital Methods. 9. Extended Gaussian-type basis for molecular-orbital studies of organic molecules." *J. Chem. Phys.* **54**, 724 (1971).
- Hariharan, P. C.; Pople, J. A. "The influence of polarization functions on molecular orbital hydrogenation energies." *Theor. Chim. Acta* **28**, 213 (1973).
- Krishnan, R.; Binkley, J. S.; Seeger, R.; Pople, J. A. (6-311G basis set family).
- Dunning, T. H. "Gaussian basis sets for use in correlated molecular calculations. I. The atoms boron through neon and hydrogen." *J. Chem. Phys.* **90**, 1007 (1989) (correlation-consistent basis sets).
- Kendall, R. A.; Dunning, T. H.; Harrison, R. J. "Electron affinities of the first-row atoms revisited. Systematic basis sets and wave functions." *J. Chem. Phys.* **96**, 6796 (1992) (augmented correlation-consistent sets).
- Hay, P. J.; Wadt, W. R. "Ab initio effective core potentials for molecular calculations." *J. Chem. Phys.* **82**, 270 and 299 (1985) (LANL2DZ ECPs).
- Dolg, M.; Stoll, H.; Preuss, H. and coworkers (Stuttgart–Dresden energy-consistent pseudopotentials for transition metals, lanthanides, and actinides — multiple papers).

### 10.12 Semi-empirical methods
- Dewar, M. J. S.; Zoebisch, E. G.; Healy, E. F. "AM1: A New General Purpose Quantum Mechanical Molecular Model." *J. Am. Chem. Soc.* **107**, 3902 (1985).
- Dewar, M. J. S.; Thiel, W. "Ground States of Molecules. 38. The MNDO Method: Approximations and Parameters." *J. Am. Chem. Soc.* **99**, 4899 (1977).

### 10.13 Relativistic treatments
- Douglas, M.; Kroll, N. M. "Quantum electrodynamical corrections to the fine structure of helium." *Ann. Phys. (NY)* **82**, 89 (1974) (Douglas–Kroll transformation).
- Hess, B. A. "Relativistic electronic-structure calculations employing a two-component no-pair formalism with external-field projection operators." *Phys. Rev. A* **33**, 3742 (1986) (DKH Hamiltonian).

---

*Note: this document synthesizes information from Gaussian, Inc.'s official technical documentation (capabilities, DFT methods, and references pages at gaussian.com) together with independent literature on Gaussian's methodology, current as of August 2026.*

---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the Gaussian 	Widely used commercial quantum chemistry package supporting a vast range of DFT functionals plus HF and post-HF methods for molecules. Also provide a list of publications related to the package's theory. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
