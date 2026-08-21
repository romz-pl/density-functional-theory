# All-Electron LAPW Codes with Roots in Surface/Thin-Film Electronic Structure: The *flair*/FLAPW Family

## A note on naming, up front

The query name **"FHI96"** conflates two genuinely distinct code lineages, and getting this right matters for anyone trying to actually use or cite these codes:

- **FHI96 / `fhi96md`** is a real, published code (Bockstedte, Kley, Neugebauer & Scheffler, *Comput. Phys. Commun.* 1997, distributed through the CPC Program Library) from the Fritz-Haber-Institut der Max-Planck-Gesellschaft, Berlin. It is a **plane-wave pseudopotential** total-energy/molecular-dynamics code (LDA/GGA, first-principles norm-conserving pseudopotentials), a predecessor in spirit to `fhi98md`/`fhi98PP` and eventually FHI-aims. It is **not** an all-electron LAPW code, and it has no surface/thin-film-specific lineage distinct from the rest of the FHI pseudopotential program — its "surfaces" capability is the generic supercell-slab capability any plane-wave pseudopotential code has.
- The code family that actually matches the description in the request — **all-electron, LAPW-based, with historical roots specifically in surface and thin-film electronic structure** — is the **FLAPW (Full-potential Linearized Augmented Plane Wave) lineage descending from the original Northwestern/Freeman film program**, whose modern direct descendant maintained under a compact, actively developed codebase is ***flair*** (Weinert, Univ. of Wisconsin–Milwaukee, with Oregon State, U. Vienna, and TU Vienna collaborators). Sibling/cousin codes in the same LAPW family tree include **FLEUR** (Jülich), **WIEN2k** (Vienna), **Elk**, and **exciting**.

Given the substance of the request (all-electron LAPW, surface/thin-film roots), this review treats **flair** as the primary subject, places it in the FLAPW family tree, and is explicit about where the FHI96/fhi96md naming diverges so nothing is silently misattributed.

---

## 1. Origins: why LAPW methods have surface/thin-film roots at all

Unlike most all-electron LAPW codes, whose lineage starts from bulk band-structure calculations, the FLAPW family's origin story runs in the *opposite* direction — thin films first, bulk second. This is the detail that makes "roots in surface and thin-film electronic structure" a technically precise, not just marketing, description.

- **1975** — O. K. Andersen formalized the **linear methods in band theory** (the "L" in LAPW), giving APW basis functions and their energy derivatives a rigorous linearized form (*Phys. Rev. B* **12**, 3060).
- **1975** — Koelling and Arbman demonstrated early practical use of the energy-derivative trick in an APW context (*J. Phys. F* **5**, 2041).
- **1979** — **Krakauer, Posternak, and Freeman** published the paper that is the true genesis of this family: *"Linearized augmented plane-wave method for the electronic band structure of thin films,"* *Phys. Rev. B* **19**, 1706 (1979). This is the first LAPW formulation built explicitly around a film geometry — semi-infinite vacuum regions above and below a slab, rather than a 3D-periodic bulk cell — making it directly suited to surface science.
- **1981** — **Wimmer, Freeman, Krakauer, and Weinert**, *"Full-potential self-consistent linearized-augmented-plane-wave method for calculating the electronic structure of molecules and surfaces: O₂ molecule,"* *Phys. Rev. B* **24**, 864 (1981). This introduced the **full-potential** (as opposed to muffin-tin/shape-approximated) treatment into the LAPW-for-films framework — the paper usually cited as the birth of "FLAPW" proper.
- **1982** — **Weinert, Wimmer, and Freeman**, *"Total-energy all-electron density functional method for bulk solids and surfaces,"* *Phys. Rev. B* **26**, 4571 (1982), extended the formalism to total-energy calculations, making structural relaxation and energetics on surfaces tractable.
- **1984** — Jansen and Freeman (Northwestern) then took the film-based FLAPW machinery and reformulated it for **bulk solids**, explicitly noting in the paper that the bulk method follows "as in the thin-film FLAPW approach" (*Phys. Rev. B* **30**, 561) — i.e., bulk FLAPW is historically the derivative case, not the ancestor.

The upshot: the entire FLAPW methodological tradition — and every code descending from it (flair, FLEUR, WIEN2k, Elk, exciting) — inherits a basis-set formalism in which the **film geometry (2D-periodic slab with true semi-infinite vacuum, not a periodically repeated slab-plus-vacuum supercell)** is a first-class citizen of the method, not a bolted-on afterthought. This is what distinguishes the family from bulk-first plane-wave/pseudopotential codes (VASP, Quantum ESPRESSO, CASTEP, fhi96md/fhi98md/FHI-aims) that simulate surfaces via artificial periodic slab supercells.

---

## 2. What "all-electron LAPW" means in this family

- **Basis set**: Space is partitioned into non-overlapping atom-centered **muffin-tin spheres** and an **interstitial region** (plus, for film geometry, **vacuum regions**). Inside spheres, basis functions are numerical radial solutions of the spherical part of the potential times spherical harmonics, augmented by their energy derivatives (the "linearization"); in the interstitial region (and vacuum, for films) they are plane waves. Continuity of value and slope is enforced at the sphere/vacuum boundaries.
- **Full potential**: unlike the older "muffin-tin approximation," which spherically averages the potential inside spheres and flattens it in the interstitial region, **full-potential** LAPW (FLAPW) represents the true, general-shape potential and density everywhere — essential for surfaces, low-symmetry structures, and accurate total energies/forces.
- **All-electron**: no pseudopotential approximation — core and valence electrons are treated on the same footing, typically with core states solved via the (scalar-)relativistic or full Dirac radial equation and valence states via the LAPW variational basis. This gives access to core-level and hyperfine properties (EFGs, hyperfine fields, XPS/XANES-relevant core shifts) that pseudopotential methods cannot directly provide.
- **Local orbitals**: added basis functions (beyond APW + energy derivative) to correctly describe semi-core states and to reduce linearization error — flair is noted as one of the earliest codes to implement explicit core-orthogonalization to handle semi-core states cleanly.

---

## 3. *flair*: the modern direct descendant

### 3.1 Identity and development team

*flair* (stylized lower-case) is described by its maintainers as "an implementation of the Full-potential Linearized Augmented Plane Wave (FLAPW) method for bulk and thin films." Its own documentation is explicit about the lineage: *"flair's roots date back to the original FLAPW code."* Development is a multi-institution collaboration:

| Institution | Lead |
|---|---|
| University of Wisconsin–Milwaukee | Michael Weinert |
| Oregon State University | G. Schneider |
| University of Vienna | R. Podloucky |
| Technical University of Vienna | J. Redinger |

These groups collectively span the people who did the *original* FLAPW development (Weinert co-authored the 1981/1982 foundational papers above) and decades of subsequent surface-, film-, and bulk-electronic-structure work, which is the basis for the site's claim that the collaboration has "a long history of working with the FLAPW method — including its original development."

### 3.2 Scope: bulk *and* films, natively

*flair* explicitly supports **both bulk 3D-periodic systems and thin-film (2D-periodic, semi-infinite-vacuum) geometries** within the same code and basis-set machinery, rather than treating films as a special-cased afterthought. This directly continues the Krakauer–Posternak–Freeman (1979) film-LAPW → Weinert–Wimmer–Freeman (1981/82) full-potential-film lineage, later generalized to bulk by Jansen and Freeman (1984).

### 3.3 Notable/originating capabilities

The *flair* homepage highlights several features it claims to have been the *first* code to offer:
- **k-projected band structures** — resolving bulk-like bands projected onto a 2D surface Brillouin zone, a diagnostic central to surface-state/resonance identification in ARPES-comparison work.
- **Explicit core-orthogonalization** to correctly treat **semi-core states** (states energetically close to valence but spatially core-like), avoiding the "ghost band"/linearization-error problems these states can otherwise cause in LAPW basis sets.

Additional stated output/analysis capabilities:
- Density of states
- Band structures (including the k-projected variant above)
- Charge/spin density plots
- STM (scanning tunneling microscopy) image simulation — a capability whose relevance is itself a marker of the code's surface-science heritage, since STM simulation is of little use for a bulk-only code.

### 3.4 Usability philosophy

Despite having **no GUI**, the developers emphasize ease of use through:
- Flexible, largely structural-only required input (most numerical/technical parameters — e.g., LAPW energy parameters — are computed and automatically updated by the code rather than requiring manual tuning).
- "Intelligent defaults" derived from the collaboration's decades of accumulated FLAPW experience.
- A Python front end for input generation was reported as in development (status as of the cited documentation; not confirmed as released).

### 3.5 Technical/build environment

| Aspect | Detail |
|---|---|
| Parallelism | Serial and parallel (MPI, specifically OpenMPI) |
| Primary development OS | Linux |
| Compilers | gfortran and Intel `ifort` |
| Windows support | Serial-only, gfortran-compiled; no working OpenMPI setup attempted/tested on Windows |
| Distribution | Source code, from the flair homepage; **not open-source in the "download freely" sense** — access requires emailing `weinert@uwm.edu` for information and registration |
| License model | Registration/gatekept access rather than a public GPL-style release |

This registration-gated distribution model is a notable contrast with sibling FLAPW-family codes: **Elk** and **exciting** are both released under the GNU GPL with unrestricted public download, while **FLEUR** (Jülich, MaX EU project) and **WIEN2k** (TU Wien) sit at intermediate points (FLEUR is open-source under Jülich/MaX stewardship; WIEN2k is a paid commercial-academic license). *flair* is the most restrictively distributed of the major FLAPW-lineage codes, consistent with its history as a smaller, tightly-held collaborative research tool rather than a large community-support project.

---

## 4. *flair* in context: the wider FLAPW/LAPW code family tree

| Code | Institution(s) | Bulk | Film (true 2D + vacuum) | License/distribution | Notable distinguishing features |
|---|---|---|---|---|---|
| **flair** | UW–Milwaukee, Oregon State, U. Vienna, TU Vienna | Yes | Yes (native, ancestral capability) | Registration-gated source | k-projected bands, explicit semi-core orthogonalization, STM simulation, minimal-input philosophy |
| **FLEUR** | Forschungszentrum Jülich (Peter Grünberg Institute); MaX EU Centre of Excellence | Yes | Yes | Open source (Jülich/MaX) | Non-collinear magnetism + spin-orbit, external electric fields, spin-dependent transport, 100,000+ lines, 20+ years of Jülich development |
| **WIEN2k** | TU Wien (Blaha, Schwarz et al.) | Yes | Via slab supercell (not native semi-infinite-vacuum film basis in the same sense) | Commercial-academic license | Very large user base; extensively used for hyperfine properties, EFGs, and structural optimization (`fhi95force`-style force/relaxation methods were historically implemented into WIEN by collaborators, illustrating cross-pollination between the FHI pseudopotential-code group and the WIEN LAPW-code group) |
| **Elk** | Originally Karl-Franzens-Universität Graz (EXCITING EU RTN) | Yes | Not a primary design target | GNU GPL, fully open | Non-collinear spins as a free vector field, libxc interface, "as simple as possible" developer philosophy |
| **exciting** | Karl-Franzens-Universität Graz and collaborators | Yes | Not a primary design target | GNU GPL, fully open | Emphasis on excited-state (BSE/GW-adjacent) properties beyond ground-state DFT; developer-friendly codebase goal |

All five codes trace their basis-set formalism to the same 1975–1982 Andersen/Koelling-Arbman/Krakauer-Posternak-Freeman/Wimmer-Freeman-Krakauer-Weinert lineage described in Section 1, and all cite essentially the same core reference set (Andersen 1975; Koelling & Arbman 1975; Krakauer, Posternak & Freeman 1979; Wimmer, Freeman, Krakauer & Weinert 1981; Weinert, Wimmer & Freeman 1982), confirmed independently by both the FLEUR and *flair*-adjacent reference pages. FLEUR in particular is the code that has carried the **film-native** capability furthest in parallel with *flair*, given its Jülich lineage's own deep roots in surface magnetism.

---

## 5. Where FHI-associated codes fit (and don't fit) relative to this family

To close the loop on the "FHI96" naming question from Section 0:

- **`fhi96md`** (1997, Bockstedte/Kley/Neugebauer/Scheffler) → plane-wave, norm-conserving pseudopotential, LDA/GGA, FORTRAN 77, ~14,000 lines, designed to run on then-modest hardware (IBM RS/6000, Pentium PCs). It handles surfaces the way essentially all plane-wave pseudopotential codes do: periodically repeated slab supercells with vacuum padding. It is a precursor to `fhi98md`/`fhi98PP`, which in turn fed methodologically (not code-wise) into the much later, unrelated-in-implementation **FHI-aims** (2003–present), an all-electron code — but one based on **numeric atom-centered orbitals (NAOs)**, not an LAPW basis. FHI-aims is "all-electron" like flair/FLEUR/WIEN2k/Elk/exciting, but it is architecturally a different method entirely (NAO vs. augmented-plane-wave), and its own developer history does not claim thin-film/surface-LAPW ancestry — its lineage is the FHI pseudopotential-code tradition generalized to an all-electron NAO basis.
- **FHI-gap** (Jiang, Gómez-Ábal, Li, Meisenbichler, Ambrosch-Draxl & Scheffler, 2013) is a further, separate FHI-associated code: a **GW quasiparticle code built on top of WIEN2k's augmented-plane-wave output**, illustrating a genuine (if indirect and much later) point of contact between the FHI computational tradition and the LAPW family — but as an add-on package consuming WIEN2k wavefunctions, not as a standalone LAPW ground-state code in its own right.

None of these FHI-lineage codes is the "FHI96" the request seems to be pointing at if the intent was an LAPW/surface-rooted code; that description matches *flair* and its FLAPW siblings, not the FHI pseudopotential/NAO lineage.

---

## 6. Summary assessment

- **If the goal is an all-electron LAPW code whose method genuinely originates in surface/thin-film physics**, *flair* is the correct and most direct match: it is maintained by the same personnel lineage (Weinert et al.) responsible for the original 1981–1982 full-potential film-LAPW papers, natively treats both bulk and true semi-infinite-vacuum film geometries in one code, and offers surface-diagnostic features (k-projected bands, STM simulation) that only make sense for a code built around surface physics from the start.
- **Its closest living relatives** are FLEUR (most similar in film-native scope and magnetism sophistication, but with a much larger, openly-distributed Jülich-maintained codebase) and WIEN2k (much larger user base and hyperfine-property track record, but film treatment via supercell rather than native vacuum basis).
- **Its main practical limitation relative to those peers** is distribution: source access requires direct registration with the UWM group rather than a public GPL download, which constrains its adoption relative to Elk, exciting, and FLEUR.
- **"FHI96"/`fhi96md`** is a real, separately-documented FHI pseudopotential-plane-wave code and should not be conflated with this LAPW family; if the user's actual target was an FHI-branded all-electron code, **FHI-aims** (NAO-based, not LAPW) or **FHI-gap** (a WIEN2k-based GW add-on) are the genuine FHI-associated all-electron options, though neither carries the surface/thin-film-native LAPW ancestry that *flair* and FLEUR do.


---

> [!NOTE]
> 
> Generated by Claude.ai
>
> Model: Sonet 5
>
> Prompt: Create an exhaustive review of the FHI96 / FLAIR 	All-electron LAPW-based codes with historical roots in surface and thin-film electronic structure studies. Show the output in Markdown format. Do not copy the output of the exported files into the chat.
